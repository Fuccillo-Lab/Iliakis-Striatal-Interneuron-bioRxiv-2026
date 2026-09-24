# Optogenetic inhibition behavior

This directory contains the MATLAB analyses for the optogenetic behavioral
results in Figure 4C–F and Figure S14 of:

> Iliakis EA et al. *Striatal interneuron microcircuits gate reinforcement to
> stabilize adaptive choice.* bioRxiv (2026).
> [https://doi.org/10.64898/2026.09.21.753300](https://doi.org/10.64898/2026.09.21.753300)

Figure 4A is artwork, Figure 4B was assembled from Figure 2 data, and the
remaining Figure 4 panels belong to other experimental datasets.

## Contents

- `plot_figure4_behavior.m` reproduces Figure 4C–F.
- `plot_figureS14_policy.m` reproduces Figure S14A–D.
- `README.md` documents the shared inputs, analyses, and outputs.

## Data availability

The data are not stored in this GitHub repository. The analysis dataset will
be archived separately on Zenodo.

**Zenodo record:** forthcoming

After downloading and extracting the Zenodo archive, set `dataRoot` below to
the local `OptogeneticsData` directory. Both scripts use the same four inputs:

| File | Role |
| --- | --- |
| `RLGLMHMM/ModelInput/modelTfinal_for_rlhmm.csv` | Model-ready trial table produced by `prepare_rlhmm_table.m` |
| `RLGLMHMM/ModelInput/rlhmm_plotting_metadata.csv` | Trial-keyed phase and delta-Q metadata |
| `RLGLMHMM/FinalFit/rlhmm_posteriors.csv` | Trial keys retained in the selected final RL-GLM-HMM fit |
| `animals.csv` | Animal cell-type and experimental-condition assignments |

The required columns and validation checks are documented in each function's
header. The posterior file is used only to retain the exact historical trial
cohort; state probabilities are not used for these figures.

## Usage

From the repository root in MATLAB:

```matlab
addpath("opto_behavior");

dataRoot = "path/to/OptogeneticsData";
outputFolder = "path/to/output/folder";

trialFile = fullfile(dataRoot, "RLGLMHMM", "ModelInput", ...
    "modelTfinal_for_rlhmm.csv");
metadataFile = fullfile(dataRoot, "RLGLMHMM", "ModelInput", ...
    "rlhmm_plotting_metadata.csv");
posteriorFile = fullfile(dataRoot, "RLGLMHMM", "FinalFit", ...
    "rlhmm_posteriors.csv");
animalFile = fullfile(dataRoot, "animals.csv");

results4 = plot_figure4_behavior( ...
    trialFile, metadataFile, posteriorFile, animalFile, outputFolder);

resultsS14 = plot_figureS14_policy( ...
    trialFile, metadataFile, posteriorFile, animalFile, outputFolder);
```

MATLAB's Statistics and Machine Learning Toolbox is required. Figure 4 uses
`fitglme`; Figure S14 uses `fitglm` and `fitlme`.

## Figure S14 analysis

Panels A–B show P(push) across delta-Q bins (`-9:2:9`) for SST and TH cohorts.
Trials are first averaged within animal and delta-Q bin; a bin contributes to
the group mean only when that animal has at least 10 trials in it. Animals are
then weighted equally in each mean and SEM.

Panel C is a schematic of a conventional single-slope policy and the split-beta
policy used in panel D.

For panel D, positive and negative delta-Q values are entered as separate
predictors in a per-animal, per-phase logistic model:

```text
P(push) ~ beta_pos * max(deltaQ, 0) + beta_neg * min(deltaQ, 0)
```

Both predictors are scaled by the pooled delta-Q standard deviation within
cell type. Fits require at least 40 trials, push probability between 0.05 and
0.95, at least 15 trials on each side of zero, and at least three push and
three pull choices on each side. Fits with either absolute slope above 6 are
excluded as unstable. Each slope is then tested with:

```text
slope ~ condition * phase + (1 | animalID)
```

Significance markers denote inhibition-versus-control coefficients with
*p* < 0.05. At PRE, this is the group contrast; at later phases, it is the
condition-by-phase interaction relative to PRE.

## Outputs

Figure 4 exports its composite figure and animal-, group-, interaction-, and
coefficient-level CSV files as documented in `plot_figure4_behavior.m`.

Figure S14 exports:

- `figureS14_panelsA_B_psychometric.png` and `.pdf`
- `figureS14_panelC_schematic.png` and `.pdf`
- `figureS14_panelD_split_beta.png` and `.pdf`
- `figureS14_psychometric_animal_bin.csv`
- `figureS14_psychometric_group_summary.csv`
- `figureS14_split_beta_animal_phase.csv`
- `figureS14_split_beta_group_summary.csv`
- `figureS14_split_beta_interaction_stats.csv`
- `figureS14_split_beta_coefficient_stats.csv`
- `figureS14_split_beta_qc_summary.csv`

The three S14 panel files are separate because the final manuscript composite
was assembled during figure layout. The returned `resultsS14` structure also
contains all figure handles and underlying tables.
