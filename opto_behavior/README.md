# Optogenetic inhibition behavior

This directory contains the MATLAB analysis used to generate the behavioral results in Figure 4C–F of:

> Iliakis EA et al. *Striatal interneuron microcircuits gate reinforcement to stabilize adaptive choice.* bioRxiv (2026). [https://doi.org/10.64898/2026.09.21.753300](https://doi.org/10.64898/2026.09.21.753300)

Figure 4A is artwork, Figure 4B was assembled from Figure 2 data, and the remaining panels belong to other experimental datasets.

## Contents

- `plot_figure4_behavior.m` reproduces the data-driven behavioral panels in Figure 4C–F.
- `README.md` documents the required inputs and analysis outputs.

## Data availability

The data required to run this analysis are not stored in this GitHub repository. The analysis dataset will be archived separately on Zenodo.

**Zenodo record:** forthcoming

After downloading and extracting the Zenodo archive, set `dataRoot` in the example below to the local `OptogeneticsData` directory. The expected inputs are:

| File | Role |
| --- | --- |
| `RLGLMHMM/ModelInput/modelTfinal_for_rlhmm.csv` | Model-ready trial table produced by `prepare_rlhmm_table.m` |
| `RLGLMHMM/ModelInput/rlhmm_plotting_metadata.csv` | Trial-keyed phase metadata |
| `RLGLMHMM/FinalFit/rlhmm_posteriors.csv` | Trial keys retained in the selected final RL-GLM-HMM fit |
| `animals.csv` | Animal cell-type and experimental-condition assignments |

The required columns and validation checks for each input are documented in the header of `plot_figure4_behavior.m`.

## Usage

From the repository root in MATLAB:

```matlab
addpath("opto_behavior");

dataRoot = "path/to/OptogeneticsData";
outputFolder = "path/to/output/folder";

results4 = plot_figure4_behavior( ...
    fullfile(dataRoot, "RLGLMHMM", "ModelInput", "modelTfinal_for_rlhmm.csv"), ...
    fullfile(dataRoot, "RLGLMHMM", "ModelInput", "rlhmm_plotting_metadata.csv"), ...
    fullfile(dataRoot, "RLGLMHMM", "FinalFit", "rlhmm_posteriors.csv"), ...
    fullfile(dataRoot, "animals.csv"), ...
    outputFolder);
```

MATLAB's Statistics and Machine Learning Toolbox is required for `fitglme`.

## Analysis details

The analysis and the RL-GLM-HMM models use the same optogenetic experiment. The posterior file is used only to retain the exact historical analysis cohort through its trial keys and choice variable (`Y`); state probabilities are not used to generate Figure 4C–F.

The function first averages trials within animal × phase and then gives animals equal weight in each group mean and SEM. For Figure 4C/E, it plots P(pull on the next trial) following a rewarded push, the complement of the source analysis's push win-stay measure. For Figure 4D/F, it plots P(pull on the next trial) following a rewarded pull. The statistical tests preserve the source analysis's separate consecutive-trial filtering and random intercepts for animal and session.

Significance markers indicate coefficients with *p* < 0.05. At PRE, the marked coefficient is the inhibition-versus-control group contrast. At later phases, the marked coefficient is the condition-by-phase interaction relative to PRE.

## Outputs

When `outputFolder` is supplied, the function exports:

- `figure4_panelsC_to_F.png`
- `figure4_panelsC_to_F.pdf`
- `figure4_panelsC_to_F_animal_bin.csv`
- `figure4_panelsC_to_F_group_summary.csv`
- `figure4_panelsC_to_F_interaction_stats.csv`
- `figure4_panelsC_to_F_coefficient_stats.csv`

The returned `results4` structure also contains the figure handle and the underlying animal-level, group-summary, interaction-test, and coefficient-test tables.
