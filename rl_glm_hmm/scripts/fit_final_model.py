#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import argparse

# Locate the project-specific RL-GLM-HMM module relative to this script.
MODEL_DIR = Path(__file__).resolve().parent.parent

if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))

import rl_glm_hmm as fastmod
from rl_glm_hmm import RLHMMHybrid


# =========================================================
# SETTINGS
# =========================================================

Ks = [3]
n_seeds = 10
maxiter = 100
tol = 1e-3
fit_init_states = False

reward_scale = 8.0
q_init = 4.0

# Ridge penalties
l2 = 1e-3
l2_extra = 1e-3

obs_maxiter = 100

fit_alpha = True
alpha_fixed = 0.5

# Initialization mode
#   False = random alpha/beta initialization from generate_params; best for exploratory K=3+ fits
#   True  = warm-start alpha/beta from warm_alpha_single / warm_beta_single below
use_warm_start = False

# Optional warm-start values. These must either be empty or have length k.
# For K=2, the values below are your stable nested single-beta solution.
warm_alpha_single = np.array([0.4293581, 0.5108002], dtype=float)
warm_beta_single  = np.array([0.1276689, 0.1352828], dtype=float)

# Hybrid predictors to add to RL softmax.
# Start with clean WS/LS terms only.
extra_cols = [
    "prevChoice_prevReward",
    "prevChoice_prevUnreward",
    "animalBias",
    "recentChoiceFrac_z",
    "prevChoice_switchFrac_z",
]

meta_cols = ["animalID", "iOrig", "jOrig", "Y", "currChoice", "currReward"] + extra_cols

# =========================================================
# HELPERS
# =========================================================

def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def build_session_boundaries(df):
    sess = [0]
    for _, g in df.groupby(["animalID", "iOrig"], sort=False):
        sess.append(sess[-1] + len(g))
    return np.array(sess, dtype=int)


def get_final_ll(lls):
    s = pd.Series(lls).dropna()
    return float(s.iloc[-1]) if len(s) > 0 else np.nan


def logit(p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


def warm_start_alpha_and_betas(
    theta_init,
    k,
    seed,
    alpha_single=None,
    beta_single=None,
    jitter_alpha=0.05,
    jitter_beta=0.03,
):
    """
    Optionally warm-start ONLY alpha and beta_pos/beta_neg.

    This function is now K-flexible: alpha_single and beta_single must either
    both be supplied with length k, or this will error clearly. If you do not
    want a warm start, set use_warm_start = False and this function will not be
    called.

    Everything else in theta_init is left as the random initialization from
    generate_params:
        - bias
        - prevChoice_prevReward
        - prevChoice_prevUnreward
        - animalBias
        - recentChoiceFrac_z
        - prevChoice_switchFrac_z
    """

    rng = np.random.default_rng(seed)
    theta = np.asarray(theta_init, dtype=float).copy()

    if alpha_single is None or beta_single is None:
        raise ValueError(
            "Warm start requested, but alpha_single/beta_single were not supplied. "
            "Either set use_warm_start=False or provide arrays of length k."
        )

    alpha_single = np.asarray(alpha_single, dtype=float)
    beta_single = np.asarray(beta_single, dtype=float)

    if alpha_single.size != k or beta_single.size != k:
        raise ValueError(
            f"Warm start requested for K={k}, but alpha_single has length "
            f"{alpha_single.size} and beta_single has length {beta_single.size}. "
            "Either set use_warm_start=False for random initialization, or provide "
            "warm-start arrays with one value per state."
        )

    ix = 0

    if fit_alpha:
        alpha_jit = alpha_single + rng.normal(0.0, jitter_alpha, size=k)
        alpha_jit = np.clip(alpha_jit, 0.05, 0.95)
        theta[ix:ix+k] = logit(alpha_jit)
        ix += k

    ix += k  # bias stays random

    beta_pos = beta_single + rng.normal(0.0, jitter_beta, size=k)
    beta_neg = beta_single + rng.normal(0.0, jitter_beta, size=k)

    theta[ix:ix+k] = beta_pos
    ix += k
    theta[ix:ix+k] = beta_neg
    ix += k

    # Remaining gamma_extra terms stay random.
    return theta

def fit_one_model(Y, rewards, X_extra, sess, k, seed):

    np.random.seed(seed)

    n = len(Y)
    c = 2

    model = RLHMMHybrid(
        n=n,
        c=c,
        k=k,
        n_extra=X_extra.shape[1],
        reward_scale=reward_scale,
        q_init=q_init,
        fit_alpha=fit_alpha,
        split_beta=True,
        alpha_fixed=alpha_fixed,
        verbose=True,  # faster; seed-level progress still prints below
    )

    A_init, theta_init, pi0_init = model.generate_params(
        transitions=['dirichlet', 10, 1],
        state_priors='uniform',
        seed_scale=0.2,
    )

    if use_warm_start:
        theta_init = warm_start_alpha_and_betas(
            theta_init=theta_init,
            k=k,
            seed=seed,
            alpha_single=warm_alpha_single,
            beta_single=warm_beta_single,
            jitter_alpha=0.05,
            jitter_beta=0.03,
        )

    lls, A_fit, theta_fit, pi0_fit = model.fit(
        y=Y,
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
        "lls": np.array(lls, dtype=float).copy(),
        "final_ll": get_final_ll(lls),
        "gamma": model.pStates.copy(),
        "zhat": model.states.astype(int).copy(),
        "k": k,
        "seed": seed,
    }


def unpack_theta_for_table(theta, k, extra_cols):

    theta = np.asarray(theta)
    ix = 0

    if fit_alpha:
        alpha = 1.0 / (1.0 + np.exp(-theta[ix:ix+k]))
        ix += k
    else:
        alpha = np.full(k, alpha_fixed, dtype=float)

    bias = theta[ix:ix+k]
    ix += k

    beta_pos = theta[ix:ix+k]
    ix += k

    beta_neg = theta[ix:ix+k]
    ix += k

    n_extra = len(extra_cols)
    if n_extra > 0:
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

        for j, col in enumerate(extra_cols):
            row[col] = gamma_extra[s, j]

        rows.append(row)

    return pd.DataFrame(rows)


def get_state_order_by_prev_reward(theta, k):
    """
    Sort states so State 1 has the highest prevChoice_prevReward coefficient.
    This is only label alignment after fitting; it does not affect the objective.
    """

    ix = 0

    if fit_alpha:
        ix += k

    ix += k  # bias
    ix += k  # beta_pos
    ix += k  # beta_neg

    n_extra = len(extra_cols)
    gamma_extra = theta[ix:ix+k*n_extra].reshape(k, n_extra)

    col = extra_cols.index("prevChoice_prevReward")
    score = gamma_extra[:, col]

    return np.argsort(-score)


def reorder_fit_result(res, order):

    res = res.copy()
    k = res["k"]
    n_extra = len(extra_cols)

    old_to_new = np.empty_like(order)
    old_to_new[order] = np.arange(k)

    res["A_fit"] = res["A_fit"][order][:, order]
    res["pi0_fit"] = res["pi0_fit"][order]
    res["gamma"] = res["gamma"][:, order]
    res["zhat"] = old_to_new[res["zhat"]]

    theta = res["theta_fit"].copy()
    parts = []
    ix = 0

    if fit_alpha:
        alpha = theta[ix:ix+k][order]
        parts.append(alpha)
        ix += k

    bias = theta[ix:ix+k][order]
    ix += k

    beta_pos = theta[ix:ix+k][order]
    ix += k

    beta_neg = theta[ix:ix+k][order]
    ix += k

    parts.extend([bias, beta_pos, beta_neg])

    if n_extra > 0:
        gamma_extra = theta[ix:ix+k*n_extra].reshape(k, n_extra)
        gamma_extra = gamma_extra[order, :]
        parts.append(gamma_extra.reshape(-1))

    res["theta_fit"] = np.concatenate(parts)
    res["state_order"] = order.copy()

    return res


def save_posteriors(df, gamma, zhat, out_csv):

    cols = [c for c in meta_cols if c in df.columns]
    out = df[cols].copy()

    for state_idx in range(gamma.shape[1]):
        out[f"stateProb_{state_idx+1}"] = gamma[:, state_idx]

    out["stateMAP"] = zhat + 1
    out.to_csv(out_csv, index=False)


def save_loglik_trace(lls, out_csv):
    pd.DataFrame({
        "iter": np.arange(len(lls)),
        "loglik": lls,
    }).to_csv(out_csv, index=False)


def save_transition_matrix(A_fit, out_csv):
    pd.DataFrame(A_fit).to_csv(out_csv, index=False)


def save_state_occupancy(gamma, out_csv):

    occ = gamma.mean(axis=0)

    pd.DataFrame({
        "state": np.arange(1, len(occ) + 1),
        "occupancy_mean_posterior": occ,
    }).to_csv(out_csv, index=False)


def summarize_fit_result(res):

    row = {
        "k": res["k"],
        "seed": res["seed"],
        "final_ll": res["final_ll"],
    }

    occ = res["gamma"].mean(axis=0)

    for s in range(res["k"]):
        row[f"occupancy_state{s+1}"] = occ[s]

    return row

def parse_args():
    parser = argparse.ArgumentParser(
        description="Fit the final three-state RL-GLM-HMM to the complete dataset."
    )
    parser.add_argument("--csv", required=True, help="Trial-level input CSV.")
    parser.add_argument("--out-dir", required=True, help="Output directory.")
    return parser.parse_args()

# =========================================================
# MAIN
# =========================================================

def main():

    args = parse_args()
    csv_path = args.csv
    out_root = args.out_dir

    ensure_dir(out_root)

    df = pd.read_csv(csv_path)

    df = df.sort_values(
        ["animalID", "iOrig", "jOrig"]
    ).reset_index(drop=True)

    required = [
        "animalID",
        "iOrig",
        "jOrig",
        "currChoice",
        "currReward",
    ] + extra_cols

    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Match RL-HMM binary valid-choice filtering.
    df = df[df["currChoice"].isin([-1, 1])].copy().reset_index(drop=True)

    # Y: 1 push, 0 pull.
    df["Y"] = (df["currChoice"] == 1).astype(int)

    Y = df["Y"].to_numpy(dtype=int)
    rewards = df["currReward"].to_numpy(dtype=float)
    X_extra = df[extra_cols].to_numpy(dtype=float).copy()

    # Replace inf with nan and catch missing predictors explicitly.
    X_extra[~np.isfinite(X_extra)] = np.nan
    if np.isnan(X_extra).any():
        bad_cols = []
        for j, col in enumerate(extra_cols):
            if np.isnan(X_extra[:, j]).any():
                bad_cols.append(col)
        raise ValueError(
            "NaN/Inf found in hybrid predictors after filtering valid choices: "
            f"{bad_cols}. Handle before fitting."
        )

    sess = build_session_boundaries(df)

    print("=" * 60)
    print("HYBRID RL-HMM MODEL DATA SUMMARY")
    print("=" * 60)
    print(f"n trials      : {len(Y)}")
    print(f"n sessions    : {len(sess)-1}")
    print(f"Ks            : {Ks}")
    print(f"n seeds       : {n_seeds}")
    print(f"reward_scale  : {reward_scale}")
    print(f"q_init        : {q_init}")
    print(f"fit_alpha     : {fit_alpha}")
    print("alpha_mode    : state-specific, bounded 0.05-0.95" if fit_alpha else "alpha_mode    : fixed")
    print(f"alpha_fixed   : {alpha_fixed}")
    if use_warm_start:
        print("init_mode     : warm-start alpha/beta; random bias/history")
        print(f"warm_alpha    : {warm_alpha_single}")
        print(f"warm_beta     : {warm_beta_single}")
    else:
        print("init_mode     : random initialization for alpha/beta/bias/history")
    print("state_sort    : descending prevChoice_prevReward")
    print(f"extra_cols    : {extra_cols}")
    print("=" * 60)

    all_results = []

    for k in Ks:

        for seed in range(n_seeds):

            print(f"\nFitting hybrid RL-HMM: k={k}, seed={seed} ...", flush=True)
            t_seed = time.time()

            res = fit_one_model(
                Y=Y,
                rewards=rewards,
                X_extra=X_extra,
                sess=sess,
                k=k,
                seed=seed,
            )

            order = get_state_order_by_prev_reward(
                res["theta_fit"],
                k,
            )

            res = reorder_fit_result(res, order)

            all_results.append(res)

            theta_table = unpack_theta_for_table(
                res["theta_fit"],
                k,
                extra_cols,
            )

            print(f"  seed done in {time.time() - t_seed:.1f}s")
            print(f"  final LL = {res['final_ll']:.6f}")
            print(f"  state order old->new by prevReward sensitivity = {order}")
            print(f"  occupancies = {np.round(res['gamma'].mean(axis=0),4)}")
            print(theta_table.to_string(index=False))

    summary_rows = []

    for res in all_results:

        k = res["k"]
        seed = res["seed"]

        fit_dir = os.path.join(
            out_root,
            f"k{k}",
            f"seed{seed}",
        )

        ensure_dir(fit_dir)

        save_posteriors(
            df,
            res["gamma"],
            res["zhat"],
            os.path.join(fit_dir, "rlhmm_posteriors.csv"),
        )

        save_transition_matrix(
            res["A_fit"],
            os.path.join(fit_dir, "transition_matrix.csv"),
        )

        save_state_occupancy(
            res["gamma"],
            os.path.join(fit_dir, "state_occupancy.csv"),
        )

        save_loglik_trace(
            res["lls"],
            os.path.join(fit_dir, "loglik_trace.csv"),
        )

        unpack_theta_for_table(
            res["theta_fit"],
            k,
            extra_cols,
        ).to_csv(
            os.path.join(fit_dir, "rl_params.csv"),
            index=False,
        )

        pd.DataFrame({
            "k": [k],
            "seed": [seed],
            "final_ll": [res["final_ll"]],
            "state_order": [",".join(map(str, res["state_order"]))],
            "alpha_mode": ["state-specific bounded 0.05-0.95"],
            "init_mode": ["warm alpha/beta only" if use_warm_start else "random alpha/beta/bias/history"],
            "use_warm_start": [use_warm_start],
            "state_sort": ["descending prevChoice_prevReward"],
            "extra_cols": [";".join(extra_cols)],
        }).to_csv(
            os.path.join(fit_dir, "fit_info.csv"),
            index=False,
        )

        summary_rows.append(
            summarize_fit_result(res)
        )

    summary = pd.DataFrame(summary_rows).sort_values(
        ["k", "final_ll"],
        ascending=[True, False],
    )

    # Identify the highest-likelihood initialization for each K.
    summary["selected_best_fit"] = False
    best_indices = summary.groupby("k")["final_ll"].idxmax()
    summary.loc[best_indices, "selected_best_fit"] = True

    # Save results from every initialization.
    summary.to_csv(
        os.path.join(out_root, "all_fits_summary.csv"),
        index=False,
    )

    # Save a separate record of the initialization selected for analysis.
    summary.loc[summary["selected_best_fit"]].to_csv(
        os.path.join(out_root, "selected_best_fits.csv"),
        index=False,
    )

    print("\nDone.")
    print(summary.to_string(index=False))
    print("\nSelected best fit(s):")
    print(
        summary.loc[summary["selected_best_fit"]].to_string(index=False)
    )
    print(f"\nSaved to: {out_root}")


if __name__ == "__main__":
    main()
