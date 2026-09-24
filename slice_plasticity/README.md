# Slice plasticity

This directory contains the MATLAB analysis used to generate the AMPA/NMDA ratio results in Figure 4H–I of:

> Iliakis EA et al. *Striatal interneuron microcircuits gate reinforcement to stabilize adaptive choice.* bioRxiv (2026). [https://doi.org/10.64898/2026.09.21.753300](https://doi.org/10.64898/2026.09.21.753300)

Figure 4G illustrates the experimental strategy. `analyze_slice_plasticity.m` reproduces the SST and TH quantitative panels.

## Data availability

The data required to run this analysis are not stored in this GitHub repository. They are archived on Zenodo with the other study data.

**Zenodo record:** forthcoming

The analysis expects two cell-level CSV files in the Zenodo archive:

| File | Description |
| --- | --- |
| `Plasticity/SSTANR.csv` | SPN AMPA/NMDA ratios after repeated SST interneuron inhibition or control treatment |
| `Plasticity/THANR.csv` | SPN AMPA/NMDA ratios after repeated TH interneuron inhibition or control treatment |

Each row represents one recorded cell. Both files require the columns `animal_id`, `condition`, and `ampa_nmda`. Conditions are encoded as `ctrl` and `ib`.

## Usage

From the repository root in MATLAB:

```matlab
addpath("slice_plasticity");

dataRoot = "path/to/Plasticity";
outputFolder = "path/to/output/folder";

results = analyze_slice_plasticity( ...
    fullfile(dataRoot, "SSTANR.csv"), ...
    fullfile(dataRoot, "THANR.csv"), ...
    outputFolder);
```

MATLAB's Statistics and Machine Learning Toolbox is required for `ttest2`.

## Analysis

The function averages AMPA/NMDA ratios across cells within each animal, making animal the inferential unit. It then compares control and inhibition animal means separately for the SST and TH experiments using two-sided Welch unequal-variance t-tests, matching the analyses reported in the manuscript.

Individual cells are shown as transparent points for visualization. Filled points represent animal means, and error bars show mean ± SEM across animals. The two panels use a common y-axis that includes every supplied cell-level observation.

## Outputs

When `outputFolder` is supplied, the function exports:

- `figure4_panelsH_to_I.png`
- `figure4_panelsH_to_I.pdf`
- `figure4_panelsH_to_I_cell_data.csv`
- `figure4_panelsH_to_I_animal_summary.csv`
- `figure4_panelsH_to_I_welch_tests.csv`

The returned `results` structure contains the figure handle, combined cell-level data, animal-level summary, and Welch-test table.
