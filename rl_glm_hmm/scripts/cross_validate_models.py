#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cross_validate_models.py

Session-fold cross-validation wrapper for the fast hybrid RL-HMM.

This uses the same train/test fold logic as the online RL-GLM CV script:
    - train on sessions where fold != held-out fold
    - score on sessions where fold == held-out fold
    - report fold-level held-out LL/trial
    - report animal-level held-out contributions

Current model variant matches the uploaded runner/module:
    RL-HMM with state-specific alpha and split beta_pos/beta_neg,
    plus extra hybrid predictors.

Required CSV columns:
    animalID, iOrig, jOrig, currChoice, currReward, fold
    plus the extra predictor columns listed in EXTRA_COLS.

If your main CSV lacks fold, pass --fold-csv with columns:
    animalID, iOrig, fold
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# Locate the project-specific RL-GLM-HMM module relative to this script.
MODEL_DIR = Path(__file__).resolve().parent.parent

if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))

import rl_glm_hmm as fastmod
from rl_glm_hmm import RLHMMHybrid

# ----------------------------- defaults -----------------------------

EXTRA_COLS = [
    "prevChoice_prevReward",
    "prevChoice_prevUnreward",
    "animalBias",
    "recentChoiceFrac_z",
    "prevChoice_switchFrac_z",
]

META_COLS = ["animalID", "iOrig", "jOrig", "Y", "currChoice", "currReward"] + EXTRA_COLS


# ----------------------------- utilities -----------------------------

def ensure_dir(path: str | Path) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def standard_error(x: Sequence[float]) -> float:
    x = np.asarray(x, dtype=float)
    good = np.isfinite(x)
    if np.sum(good) <= 1:
        return np.nan
    return float(np.nanstd(x[good], ddof=1) / np.sqrt(np.sum(good)))


def logit(p: np.ndarray | float) -> np.ndarray | float:
    p = np.clip(p, 1e-6, 1.0 - 1e-6)
    return np.log(p / (1.0 - p))


def get_final_ll(lls: Sequence[float]) -> float:
    s = pd.Series(lls).dropna()
    return float(s.iloc[-1]) if len(s) else np.nan


def build_session_boundaries(df: pd.DataFrame) -> np.ndarray:
    sess = [0]
    for _, g in df.groupby(["animalID", "iOrig"], sort=False):
        sess.append(sess[-1] + len(g))
    return np.asarray(sess, dtype=np.int64)


def merge_folds_if_needed(df: pd.DataFrame, fold_csv: Optional[str]) -> pd.DataFrame:
    if "fold" in df.columns:
        return df
    if fold_csv is None:
        raise ValueError("Input CSV has no 'fold' column. Supply --fold-csv.")

    foldT = pd.read_csv(fold_csv)
    required = ["animalID", "iOrig", "fold"]
    missing = [c for c in required if c not in foldT.columns]
    if missing:
        raise ValueError(f"Fold CSV missing required columns: {missing}")

    foldT = foldT[required].drop_duplicates()
    out = df.merge(foldT, on=["animalID", "iOrig"], how="left", validate="many_to_one")
    if out["fold"].isna().any():
        bad = out.loc[out["fold"].isna(), ["animalID", "iOrig"]].drop_duplicates().head(10)
        raise ValueError(f"Some sessions lacked folds. First bad sessions:\n{bad}")
    return out


def zscore_column_if_needed(df: pd.DataFrame, raw_col: str, z_col: str) -> pd.DataFrame:
    if z_col in df.columns:
        return df
    if raw_col not in df.columns:
        return df
    x = pd.to_numeric(df[raw_col], errors="coerce").to_numpy(dtype=float)
    mu = np.nanmean(x)
    sd = np.nanstd(x)
    if not np.isfinite(sd) or sd == 0:
        df[z_col] = 0.0
    else:
        df[z_col] = (x - mu) / sd
    print(f"Created {z_col} from {raw_col} using full-data z-scoring.")
    return df


def prepare_dataframe(csv_path: str, fold_csv: Optional[str], extra_cols: Sequence[str]) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df = merge_folds_if_needed(df, fold_csv)

    required_core = ["animalID", "iOrig", "jOrig", "currChoice", "currReward", "fold"]
    missing = [c for c in required_core if c not in df.columns]
    if missing:
        raise ValueError(f"Trial CSV missing required columns: {missing}")

    if "animalBias" not in df.columns and "bias" in df.columns:
        df["animalBias"] = df["bias"]
        print("WARNING: copied input column 'bias' to 'animalBias'.")

    df = zscore_column_if_needed(df, "recentChoiceFrac", "recentChoiceFrac_z")
    df = zscore_column_if_needed(df, "prevChoice_switchFrac", "prevChoice_switchFrac_z")

    missing = [c for c in extra_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required hybrid predictors: {missing}")

    df = df.sort_values(["animalID", "iOrig", "jOrig"]).reset_index(drop=True)
    df = df[df["currChoice"].isin([-1, 1])].copy().reset_index(drop=True)
    df["Y"] = (df["currChoice"] == 1).astype(int)

    numeric_cols = ["currReward", "fold"] + list(extra_cols)
    before = len(df)
    for c in numeric_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=numeric_cols).reset_index(drop=True)
    after = len(df)
    if after < before:
        print(f"Dropped {before-after} rows with missing/nonfinite RL-HMM variables.")

    return df


def get_arrays(df: pd.DataFrame, extra_cols: Sequence[str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    y = df["Y"].to_numpy(dtype=np.int64)
    rewards = df["currReward"].to_numpy(dtype=np.float64)
    X_extra = df[list(extra_cols)].to_numpy(dtype=np.float64)
    sess = build_session_boundaries(df)
    return y, rewards, X_extra, sess


# ----------------------------- model helpers -----------------------------

def warm_start_alpha_and_betas(
    theta_init: np.ndarray,
    k: int,
    seed: int,
    fit_alpha: bool,
    split_beta: bool,
    alpha_single: Sequence[float],
    beta_single: Sequence[float],
    jitter_alpha: float = 0.05,
    jitter_beta: float = 0.03,
) -> np.ndarray:
    """Warm-start alpha and beta(s); leave bias/history random."""
    rng = np.random.default_rng(seed)
    theta = np.asarray(theta_init, dtype=float).copy()

    alpha_single = np.asarray(alpha_single, dtype=float)
    beta_single = np.asarray(beta_single, dtype=float)
    if alpha_single.size != k or beta_single.size != k:
        raise ValueError(f"Warm-start alpha_single and beta_single must have length k={k}.")

    ix = 0
    if fit_alpha:
        alpha_jit = alpha_single + rng.normal(0.0, jitter_alpha, size=k)
        alpha_jit = np.clip(alpha_jit, 0.05, 0.95)
        theta[ix:ix+k] = logit(alpha_jit)
        ix += k

    ix += k  # bias stays random

    beta_jit = beta_single + rng.normal(0.0, jitter_beta, size=k)
    theta[ix:ix+k] = beta_jit
    ix += k
    if split_beta:
        beta_jit2 = beta_single + rng.normal(0.0, jitter_beta, size=k)
        theta[ix:ix+k] = beta_jit2

    return theta


def unpack_theta_for_table(theta: np.ndarray, k: int, extra_cols: Sequence[str], fit_alpha: bool, alpha_fixed: float, split_beta: bool) -> pd.DataFrame:
    theta = np.asarray(theta)
    ix = 0

    if fit_alpha:
        alpha = 1.0 / (1.0 + np.exp(-theta[ix:ix+k]))
        ix += k
    else:
        alpha = np.full(k, alpha_fixed, dtype=float)

    bias = theta[ix:ix+k]
    ix += k
    if split_beta:
        beta_pos = theta[ix:ix+k]
        ix += k
        beta_neg = theta[ix:ix+k]
        ix += k
    else:
        beta_single = theta[ix:ix+k]
        ix += k
        beta_pos = beta_single
        beta_neg = beta_single

    n_extra = len(extra_cols)
    if n_extra:
        gamma_extra = theta[ix:ix+k*n_extra].reshape(k, n_extra)
    else:
        gamma_extra = np.zeros((k, 0))

    rows = []
    for s in range(k):
        row = {
            "state": s + 1,
            "alpha": alpha[s],
            "bias": bias[s],
            "beta_pos": beta_pos[s],
            "beta_neg": beta_neg[s],
        }
        if not split_beta:
            row["beta_single"] = beta_pos[s]
        for j, col in enumerate(extra_cols):
            row[col] = gamma_extra[s, j]
        rows.append(row)
    return pd.DataFrame(rows)


def get_state_order_by_prev_reward(theta: np.ndarray, k: int, extra_cols: Sequence[str], fit_alpha: bool, split_beta: bool) -> np.ndarray:
    if "prevChoice_prevReward" not in extra_cols:
        return np.arange(k)

    ix = 0
    if fit_alpha:
        ix += k
    ix += k  # bias
    ix += k  # beta_pos or beta_single
    if split_beta:
        ix += k  # beta_neg

    n_extra = len(extra_cols)
    gamma_extra = theta[ix:ix+k*n_extra].reshape(k, n_extra)
    score = gamma_extra[:, list(extra_cols).index("prevChoice_prevReward")]
    return np.argsort(-score)


def reorder_theta(theta: np.ndarray, order: np.ndarray, k: int, n_extra: int, fit_alpha: bool, split_beta: bool) -> np.ndarray:
    theta = np.asarray(theta, dtype=float)
    parts = []
    ix = 0

    if fit_alpha:
        parts.append(theta[ix:ix+k][order])
        ix += k

    bias = theta[ix:ix+k][order]
    ix += k
    beta_pos = theta[ix:ix+k][order]
    ix += k
    parts.extend([bias, beta_pos])
    if split_beta:
        beta_neg = theta[ix:ix+k][order]
        ix += k
        parts.append(beta_neg)

    if n_extra:
        gamma_extra = theta[ix:ix+k*n_extra].reshape(k, n_extra)
        parts.append(gamma_extra[order, :].reshape(-1))

    return np.concatenate(parts)


def reorder_fit_result(res: Dict[str, object], order: np.ndarray, extra_cols: Sequence[str], fit_alpha: bool, split_beta: bool) -> Dict[str, object]:
    res = res.copy()
    k = int(res["k"])
    n_extra = len(extra_cols)
    old_to_new = np.empty_like(order)
    old_to_new[order] = np.arange(k)

    res["A_fit"] = res["A_fit"][order][:, order]
    res["pi0_fit"] = res["pi0_fit"][order]
    res["gamma"] = res["gamma"][:, order]
    res["zhat"] = old_to_new[res["zhat"]]
    res["theta_fit"] = reorder_theta(res["theta_fit"], order, k, n_extra, fit_alpha, split_beta)
    res["state_order"] = order.copy()
    return res


def make_model(n: int, k: int, n_extra: int, reward_scale: float, q_init: float,
               fit_alpha: bool, alpha_fixed: float, split_beta: bool, verbose: bool) -> RLHMMHybrid:
    return RLHMMHybrid(
        n=n,
        c=2,
        k=k,
        n_extra=n_extra,
        reward_scale=reward_scale,
        q_init=q_init,
        fit_alpha=fit_alpha,
        alpha_fixed=alpha_fixed,
        split_beta=split_beta,
        verbose=verbose,
    )


def score_hmm_loglik(
    y: np.ndarray,
    rewards: np.ndarray,
    X_extra: np.ndarray,
    sess: np.ndarray,
    A: np.ndarray,
    theta: np.ndarray,
    pi0: np.ndarray,
    k: int,
    reward_scale: float,
    q_init: float,
    fit_alpha: bool,
    alpha_fixed: float,
    split_beta: bool,
) -> float:
    """Held-out HMM LL, resetting forward/backward and Q traces at session boundaries."""
    scorer = make_model(
        n=len(y),
        k=k,
        n_extra=X_extra.shape[1],
        reward_scale=reward_scale,
        q_init=q_init,
        fit_alpha=fit_alpha,
        alpha_fixed=alpha_fixed,
        split_beta=split_beta,
        verbose=False,
    )
    phi = scorer.compute_phi(y, rewards, sess, theta, X_extra=X_extra)

    ll = 0.0
    for s0, s1 in zip(sess[:-1], sess[1:]):
        sl = slice(s0, s1)
        ll_s, _, _, _ = scorer.forwardPass(y[sl], A, phi[sl, :, :], pi0=pi0)
        ll += ll_s
    return float(ll)


def fit_one_model(
    train_df: pd.DataFrame,
    k: int,
    seed: int,
    extra_cols: Sequence[str],
    reward_scale: float,
    q_init: float,
    fit_alpha: bool,
    alpha_fixed: float,
    split_beta: bool,
    maxiter: int,
    tol: float,
    fit_init_states: bool,
    l2: float,
    l2_extra: float,
    obs_maxiter: int,
    warm_alpha: Sequence[float],
    warm_beta: Sequence[float],
    warm_start: bool,
    verbose: bool,
) -> Dict[str, object]:
    np.random.seed(seed)
    y, rewards, X_extra, sess = get_arrays(train_df, extra_cols)

    model = make_model(
        n=len(y),
        k=k,
        n_extra=X_extra.shape[1],
        reward_scale=reward_scale,
        q_init=q_init,
        fit_alpha=fit_alpha,
        alpha_fixed=alpha_fixed,
        split_beta=split_beta,
        verbose=verbose,
    )

    A_init, theta_init, pi0_init = model.generate_params(
        transitions=["dirichlet", 10, 1],
        state_priors="uniform",
        seed_scale=0.2,
    )

    if warm_start:
    	theta_init = warm_start_alpha_and_betas(
        	theta_init=theta_init,
        	k=k,
        	seed=seed,
        	fit_alpha=fit_alpha,
        	split_beta=split_beta,
        	alpha_single=warm_alpha,
        	beta_single=warm_beta,
    	)

    lls, A_fit, theta_fit, pi0_fit = model.fit(
        y=y,
        rewards=rewards,
        A=A_init,
        theta=theta_init,
        pi0=pi0_init,
        fit_init_states=fit_init_states,
        maxiter=maxiter,
        tol=tol,
        sess=sess,
        X_extra=X_extra,
        l2=l2,
        l2_extra=l2_extra,
        obs_maxiter=obs_maxiter,
    )

    return {
        "model": model,
        "A_fit": A_fit.copy(),
        "theta_fit": theta_fit.copy(),
        "pi0_fit": np.ravel(pi0_fit).copy(),
        "lls": np.asarray(lls, dtype=float).copy(),
        "final_ll": get_final_ll(lls),
        "gamma": model.pStates.copy(),
        "zhat": model.states.astype(int).copy(),
        "k": k,
        "seed": seed,
        "n_train_trials": len(y),
    }


def score_by_animal(best_res: Dict[str, object], test_df: pd.DataFrame, fold, model_name: str,
                    extra_cols: Sequence[str], reward_scale: float, q_init: float,
                    fit_alpha: bool, alpha_fixed: float, split_beta: bool) -> List[Dict[str, object]]:
    rows = []
    k = int(best_res["k"])
    for animal_id, animalT in test_df.groupby("animalID", sort=False):
        y, rewards, X_extra, sess = get_arrays(animalT, extra_cols)
        ll = score_hmm_loglik(
            y=y,
            rewards=rewards,
            X_extra=X_extra,
            sess=sess,
            A=best_res["A_fit"],
            theta=best_res["theta_fit"],
            pi0=best_res["pi0_fit"],
            k=k,
            reward_scale=reward_scale,
            q_init=q_init,
            fit_alpha=fit_alpha,
            alpha_fixed=alpha_fixed,
            split_beta=split_beta,
        )
        rows.append({
            "fold": fold,
            "model": model_name,
            "animalID": animal_id,
            "testLL": ll,
            "testNLL": -ll,
            "nScoredTrials": len(y),
            "testLL_perTrial": ll / len(y) if len(y) else np.nan,
        })
    return rows


# ----------------------------- CV runner -----------------------------

def run_cv(
    df: pd.DataFrame,
    out_dir: str,
    ks: Sequence[int],
    n_seeds: int,
    seed_base: int,
    extra_cols: Sequence[str],
    reward_scale: float,
    q_init: float,
    fit_alpha: bool,
    alpha_fixed: float,
    split_beta: bool,
    maxiter: int,
    tol: float,
    fit_init_states: bool,
    l2: float,
    l2_extra: float,
    obs_maxiter: int,
    warm_alpha: Sequence[float],
    warm_beta: Sequence[float],
    warm_start: bool,
    verbose: bool,
) -> None:
    ensure_dir(out_dir)

    folds = sorted(pd.unique(df["fold"]))
    fold_rows: List[Dict[str, object]] = []
    seed_rows: List[Dict[str, object]] = []
    param_rows: List[Dict[str, object]] = []
    animal_rows: List[Dict[str, object]] = []

    for fold in folds:
        trainT = df[df["fold"] != fold].copy().reset_index(drop=True)
        testT = df[df["fold"] == fold].copy().reset_index(drop=True)
        print(f"\nFold {fold} | train n={len(trainT)} | test n={len(testT)}", flush=True)

        for k in ks:
            beta_label = "splitBeta" if split_beta else "singleBeta"
            alpha_label = "stateAlpha" if fit_alpha else f"fixedAlpha{alpha_fixed:g}"
            model_name = f"RLHMM_K{k}_{alpha_label}_{beta_label}_hybrid"
            print(f"  fitting {model_name} with {n_seeds} seeds ...", flush=True)

            fold_seed_results = []
            for seed_i in range(n_seeds):
                actual_seed = seed_base + int(fold) * 1000 + seed_i if str(fold).replace('.', '', 1).isdigit() else seed_base + seed_i
                t0 = time.time()
                res = fit_one_model(
                    train_df=trainT,
                    k=k,
                    seed=actual_seed,
                    extra_cols=extra_cols,
                    reward_scale=reward_scale,
                    q_init=q_init,
                    fit_alpha=fit_alpha,
                    alpha_fixed=alpha_fixed,
                    split_beta=split_beta,
                    maxiter=maxiter,
                    tol=tol,
                    fit_init_states=fit_init_states,
                    l2=l2,
                    l2_extra=l2_extra,
                    obs_maxiter=obs_maxiter,
                    warm_alpha=warm_alpha,
                    warm_beta=warm_beta,
                    warm_start=warm_start,
                    verbose=verbose,
                )

                order = get_state_order_by_prev_reward(res["theta_fit"], k, extra_cols, fit_alpha, split_beta)
                res = reorder_fit_result(res, order, extra_cols, fit_alpha, split_beta)
                elapsed = time.time() - t0
                fold_seed_results.append(res)

                seed_rows.append({
                    "fold": fold,
                    "model": model_name,
                    "k": k,
                    "seed_index": seed_i,
                    "seed": actual_seed,
                    "trainLL": res["final_ll"],
                    "trainLL_perTrial": res["final_ll"] / res["n_train_trials"],
                    "nTrainTrials": res["n_train_trials"],
                    "elapsed_sec": elapsed,
                    "state_order": ",".join(map(str, res["state_order"])),
                })
                print(
                    f"    fold {fold} | k={k} | seed {seed_i+1}/{n_seeds} done "
                    f"(seed={actual_seed}, train LL/trial={res['final_ll']/res['n_train_trials']:.6f}, {elapsed:.1f}s)",
                    flush=True,
                )

            best_res = max(fold_seed_results, key=lambda r: r["final_ll"])
            best_seed = int(best_res["seed"])

            y_train, r_train, X_train, sess_train = get_arrays(trainT, extra_cols)
            y_test, r_test, X_test, sess_test = get_arrays(testT, extra_cols)

            # Re-score unpenalized HMM LL on train and test using the selected best seed.
            train_ll = score_hmm_loglik(
                y_train, r_train, X_train, sess_train,
                best_res["A_fit"], best_res["theta_fit"], best_res["pi0_fit"],
                k, reward_scale, q_init, fit_alpha, alpha_fixed, split_beta,
            )
            test_ll = score_hmm_loglik(
                y_test, r_test, X_test, sess_test,
                best_res["A_fit"], best_res["theta_fit"], best_res["pi0_fit"],
                k, reward_scale, q_init, fit_alpha, alpha_fixed, split_beta,
            )

            fold_rows.append({
                "fold": fold,
                "model": model_name,
                "k": k,
                "best_seed": best_seed,
                "train_n_trials": len(trainT),
                "test_n_trials": len(testT),
                "train_n_sessions": trainT[["animalID", "iOrig"]].drop_duplicates().shape[0],
                "test_n_sessions": testT[["animalID", "iOrig"]].drop_duplicates().shape[0],
                "trainLL": train_ll,
                "trainLL_perTrial": train_ll / len(y_train),
                "testLL": test_ll,
                "testNLL": -test_ll,
                "nScoredTrials": len(y_test),
                "testLL_perTrial": test_ll / len(y_test),
            })

            theta_table = unpack_theta_for_table(best_res["theta_fit"], k, extra_cols, fit_alpha, alpha_fixed, split_beta)
            for _, row in theta_table.iterrows():
                out = {"fold": fold, "model": model_name, "k": k, "best_seed": best_seed}
                out.update(row.to_dict())
                param_rows.append(out)

            animal_rows.extend(
                score_by_animal(best_res, testT, fold, model_name, extra_cols, reward_scale, q_init, fit_alpha, alpha_fixed, split_beta)
            )

            # Save per-fold best-fit artifacts.
            fit_dir = Path(out_dir) / f"fold{fold}" / f"k{k}_bestSeed{best_seed}"
            ensure_dir(fit_dir)
            pd.DataFrame({"iter": np.arange(len(best_res["lls"])), "loglik": best_res["lls"]}).to_csv(fit_dir / "loglik_trace.csv", index=False)
            pd.DataFrame(best_res["A_fit"]).to_csv(fit_dir / "transition_matrix.csv", index=False)
            theta_table.to_csv(fit_dir / "rl_params.csv", index=False)

            # Save train posterior only for the selected best seed.
            train_post = trainT[[c for c in META_COLS if c in trainT.columns]].copy()
            for state_idx in range(best_res["gamma"].shape[1]):
                train_post[f"stateProb_{state_idx+1}"] = best_res["gamma"][:, state_idx]
            train_post["stateMAP"] = best_res["zhat"] + 1
            train_post.to_csv(fit_dir / "train_posteriors.csv", index=False)

            print(
                f"    selected seed={best_seed} | test LL/trial={test_ll/len(y_test):.6f}",
                flush=True,
            )

    foldResults = pd.DataFrame(fold_rows)
    seedResults = pd.DataFrame(seed_rows)
    params = pd.DataFrame(param_rows)
    animalContribs = pd.DataFrame(animal_rows)

    agg_rows = []
    for model_name, g in foldResults.groupby("model", sort=False):
        vals = g["testLL_perTrial"].to_numpy(dtype=float)
        agg_rows.append({
            "model": model_name,
            "mean_testLL_perTrial": float(np.nanmean(vals)),
            "sem_testLL_perTrial": standard_error(vals),
            "n_folds": int(np.sum(np.isfinite(vals))),
        })
    foldAggregate = pd.DataFrame(agg_rows)

    animalModelSummary = (
        animalContribs
        .groupby(["animalID", "model"], as_index=False)
        .agg(
            totalTestLL=("testLL", "sum"),
            totalTestNLL=("testNLL", "sum"),
            totalScoredTrials=("nScoredTrials", "sum"),
            nTestFoldsRepresented=("fold", pd.Series.nunique),
        )
    )
    animalModelSummary["testLL_perTrial"] = animalModelSummary["totalTestLL"] / animalModelSummary["totalScoredTrials"]

    foldResults.to_csv(Path(out_dir) / "cv_foldResults_RLHMM.csv", index=False)
    foldAggregate.to_csv(Path(out_dir) / "cv_foldAggregate_RLHMM.csv", index=False)
    seedResults.to_csv(Path(out_dir) / "cv_seedResults_RLHMM.csv", index=False)
    params.to_csv(Path(out_dir) / "cv_params_RLHMM.csv", index=False)
    animalContribs.to_csv(Path(out_dir) / "cv_animalHeldoutContribs_RLHMM.csv", index=False)
    animalModelSummary.to_csv(Path(out_dir) / "cv_animalModelSummary_RLHMM.csv", index=False)

    info = {
        "ks": list(ks),
        "n_seeds": n_seeds,
        "seed_base": seed_base,
        "reward_scale": reward_scale,
        "q_init": q_init,
        "fit_alpha": fit_alpha,
        "alpha_fixed": alpha_fixed,
        "alpha_mode": "state-specific bounded 0.05-0.95" if fit_alpha else "fixed",
        "beta_mode": "split beta_pos/beta_neg" if split_beta else "single beta_deltaQ",
        "maxiter": maxiter,
        "tol": tol,
        "fit_init_states": fit_init_states,
        "l2": l2,
        "l2_extra": l2_extra,
        "obs_maxiter": obs_maxiter,
        "extra_cols": list(extra_cols),
        "warm_alpha": list(map(float, warm_alpha)),
        "warm_beta": list(map(float, warm_beta)),
        "warm_start": bool(warm_start),
        "numba_available": bool(getattr(fastmod, "NUMBA_AVAILABLE", False)),
        "rlhmm_module_path": getattr(fastmod, "__file__", "unknown"),
        "n_trials_after_filtering": int(len(df)),
        "n_sessions_after_filtering": int(df[["animalID", "iOrig"]].drop_duplicates().shape[0]),
    }
    with open(Path(out_dir) / "fit_info_RLHMM_CV.json", "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2)

    print("\nDone. Fold aggregate:")
    print(foldAggregate.to_string(index=False))
    print(f"\nSaved to: {out_dir}")


# ----------------------------- CLI -----------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fit state-alpha hybrid RL-HMM with session-fold CV; supports split-beta or single-beta.")
    p.add_argument("--csv", required=True, help="Trial/design CSV.")
    p.add_argument("--fold-csv", default=None, help="Optional session fold CSV with animalID,iOrig,fold if --csv lacks fold.")
    p.add_argument("--out-dir", required=True, help="Output directory.")
    p.add_argument("--ks", nargs="+", type=int, default=[2])
    p.add_argument("--n-seeds", type=int, default=5)
    p.add_argument("--seed-base", type=int, default=9000)
    p.add_argument("--maxiter", type=int, default=100)
    p.add_argument("--tol", type=float, default=1e-3)
    p.add_argument("--obs-maxiter", type=int, default=100)
    p.add_argument("--reward-scale", type=float, default=8.0)
    p.add_argument("--q-init", type=float, default=4.0)
    p.add_argument("--l2", type=float, default=1e-3)
    p.add_argument("--l2-extra", type=float, default=1e-3)
    p.add_argument("--fit-alpha", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--split-beta", action=argparse.BooleanOptionalAction, default=True,
                   help="Use beta_pos/beta_neg if true; use one beta_deltaQ per state if false.")
    p.add_argument("--alpha-fixed", type=float, default=0.5)
    p.add_argument("--fit-init-states", action=argparse.BooleanOptionalAction, default=False)
    p.add_argument("--verbose", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--warm-start", action=argparse.BooleanOptionalAction, default=True,
                   help="If true, warm-start alpha/beta from supplied values. If false, use random initialization.")
    p.add_argument("--warm-alpha", nargs="+", type=float, default=[0.4293581, 0.5108002],
                   help="Warm-start alpha values, length must match K.")
    p.add_argument("--warm-beta", nargs="+", type=float, default=[0.1276689, 0.1352828],
                   help="Warm-start beta values for both beta_pos and beta_neg, length must match K.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    df = prepare_dataframe(args.csv, args.fold_csv, EXTRA_COLS)

    print("=" * 70)
    print("HYBRID RL-HMM SPLIT-BETA STATE-ALPHA CV")
    print("=" * 70)
    print(f"csv        : {args.csv}")
    print(f"fold_csv   : {args.fold_csv}")
    print(f"out_dir    : {args.out_dir}")
    print(f"n trials   : {len(df)}")
    print(f"n sessions : {df[['animalID', 'iOrig']].drop_duplicates().shape[0]}")
    print(f"folds      : {sorted(pd.unique(df['fold']))}")
    print(f"ks         : {args.ks}")
    print(f"n seeds    : {args.n_seeds}")
    print(f"split_beta : {args.split_beta}")
    print(f"numba      : {getattr(fastmod, 'NUMBA_AVAILABLE', False)}")
    print(f"module     : {getattr(fastmod, '__file__', 'unknown')}")
    print("=" * 70)

    run_cv(
        df=df,
        out_dir=args.out_dir,
        ks=args.ks,
        n_seeds=args.n_seeds,
        seed_base=args.seed_base,
        extra_cols=EXTRA_COLS,
        reward_scale=args.reward_scale,
        q_init=args.q_init,
        fit_alpha=args.fit_alpha,
        alpha_fixed=args.alpha_fixed,
        split_beta=args.split_beta,
        maxiter=args.maxiter,
        tol=args.tol,
        fit_init_states=args.fit_init_states,
        l2=args.l2,
        l2_extra=args.l2_extra,
        obs_maxiter=args.obs_maxiter,
        warm_alpha=args.warm_alpha,
        warm_beta=args.warm_beta,
        warm_start=args.warm_start,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
