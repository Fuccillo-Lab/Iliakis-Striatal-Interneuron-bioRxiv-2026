#!/usr/bin/env python3
"""Assemble Figure S15C's paired animal-level held-out model scores.

Inputs are the animal-model summaries produced by the Q-backbone, fixed-DeltaQ
GLM, online RL-GLM, and two RL-GLM-HMM cross-validation runners. The script
selects the six submitted model specifications without changing any scores.
The output has the same columns expected by plot_rlhmm_figureS15.m.

The optional reference is the archived full masterLL.csv. It is used only to
check the assembled values and is never read to construct the output.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


MODEL_FILES = (
    ("qf", "QF"),
    ("glm", "M5_value_bias_WSLS_choiceHx_switchiness"),
    ("glm", "M6_splitValue_bias_WSLS_choiceHx_switchiness"),
    ("online", "M6_online_splitValue_bias_WSLS_choiceHx_switchiness"),
    ("k2", "RLHMM_K2_stateAlpha_splitBeta_hybrid"),
    ("k3", "RLHMM_K3_stateAlpha_splitBeta_hybrid"),
)
COLS = (
    "animalID", "model", "totalTestNLL", "totalTestLL",
    "totalScoredTrials", "nTestFoldsRepresented", "testLL_perTrial",
)
NUMERIC = COLS[2:]


def read_summary(path: Path, description: str) -> pd.DataFrame:
    if not path.is_file():
        raise ValueError(f"{description}: file not found: {path}")
    table = pd.read_csv(path)
    missing = set(COLS) - set(table.columns)
    if missing:
        raise ValueError(f"{description}: missing columns: {sorted(missing)}")
    if table.duplicated(["animalID", "model"]).any():
        raise ValueError(f"{description}: duplicate animalID/model pairs")
    for col in NUMERIC:
        table[col] = pd.to_numeric(table[col], errors="raise")
        if not np.isfinite(table[col]).all():
            raise ValueError(f"{description}: non-finite {col}")
    if (table.totalScoredTrials <= 0).any() or (table.nTestFoldsRepresented <= 0).any():
        raise ValueError(f"{description}: non-positive scored-trial or fold count")
    if not np.allclose(table.totalTestNLL, -table.totalTestLL, rtol=0, atol=1e-6):
        raise ValueError(f"{description}: totalTestNLL differs from -totalTestLL")
    if not np.allclose(
        table.testLL_perTrial,
        table.totalTestLL / table.totalScoredTrials,
        rtol=0,
        atol=1e-8,
    ):
        raise ValueError(f"{description}: testLL_perTrial differs from totalTestLL / totalScoredTrials")
    return table


def assemble(paths: dict[str, Path]) -> pd.DataFrame:
    sources = {name: read_summary(path, name) for name, path in paths.items()}
    selected = []
    for name, model in MODEL_FILES:
        rows = sources[name].loc[sources[name].model == model, list(COLS)].copy()
        if rows.empty:
            raise ValueError(f"{name}: missing model {model}")
        selected.append(rows)

    baseline = selected[0].set_index("animalID")
    for rows in selected[1:]:
        by_animal = rows.set_index("animalID")
        if set(by_animal.index) != set(baseline.index):
            raise ValueError(f"{rows.model.iloc[0]}: animal IDs do not match QF")
        counts = by_animal.loc[baseline.index, "nTestFoldsRepresented"]
        if not counts.equals(baseline.nTestFoldsRepresented):
            raise ValueError(f"{rows.model.iloc[0]}: represented-fold counts differ from QF")

    result = pd.concat(selected, ignore_index=True)
    result["model"] = pd.Categorical(
        result.model, categories=[model for _, model in MODEL_FILES], ordered=True
    )
    result = result.sort_values(["animalID", "model"]).reset_index(drop=True)
    result["model"] = result.model.astype(str)
    return result


def compare_reference(result: pd.DataFrame, reference_path: Path) -> None:
    reference = read_summary(reference_path, "reference masterLL")
    reference = reference.loc[
        reference.model.isin([model for _, model in MODEL_FILES]), list(COLS)
    ]
    keys = ["animalID", "model"]
    joined = result.merge(reference, on=keys, how="outer", indicator=True, suffixes=("_new", "_reference"))
    if len(joined) != len(result) or not joined._merge.eq("both").all():
        raise ValueError("reference masterLL has different S15C animal/model pairs")
    for col in NUMERIC:
        tolerance = 1e-8 if col == "testLL_perTrial" else 1e-6
        if not np.allclose(joined[f"{col}_new"], joined[f"{col}_reference"], rtol=0, atol=tolerance):
            raise ValueError(f"reference masterLL differs in {col}")
    print(f"Verified all {len(result)} animal/model rows against {reference_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qf", required=True, type=Path, help="Q-backbone cv_animalModelSummary.csv")
    parser.add_argument("--glm", required=True, type=Path, help="Fixed-DeltaQ cv_animalModelSummary_DQPolicyGLM.csv")
    parser.add_argument("--online", required=True, type=Path, help="Online RL-GLM animal-model summary")
    parser.add_argument("--k2", required=True, type=Path, help="Two-state RL-GLM-HMM animal-model summary")
    parser.add_argument("--k3", required=True, type=Path, help="Three-state RL-GLM-HMM animal-model summary")
    parser.add_argument("--out", required=True, type=Path, help="Output CSV for Figure S15C")
    parser.add_argument("--reference", type=Path, help="Optional archived masterLL.csv for verification")
    args = parser.parse_args()

    result = assemble({key: getattr(args, key) for key in ("qf", "glm", "online", "k2", "k3")})
    if args.reference:
        compare_reference(result, args.reference)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, index=False)
    print(f"Wrote {len(result)} rows for {result.animalID.nunique()} animals to {args.out}")


if __name__ == "__main__":
    main()
