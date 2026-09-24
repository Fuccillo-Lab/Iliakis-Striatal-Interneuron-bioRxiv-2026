# Optogenetic inhibition behavior

This directory collects the behavioral analyses of the optogenetic inhibition experiment. It currently covers Figure 4C-F; related supplementary analyses can be added here under the same experiment. Figure 4A is artwork, Figure 4B was assembled from Figure 2 data, and the slice results in Figure 4H-I belong with the slice experiment.

`scripts/plot_figure4_behavior.m` reproduces the data-driven behavioral panels C-F.

From the repository root in MATLAB:

```matlab
addpath("opto_behavior/scripts");
results4 = plot_figure4_behavior( ...
    "rl_glm_hmm/data/modelTfinal_for_rlhmm_allopto.csv", ...
    "rl_glm_hmm/data/figure5_trial_metadata.csv", ...
    "rl_glm_hmm/data/rlhmm_posteriors.csv", ...
    "rl_glm_hmm/data/animals.csv", ...
    "opto_behavior/outputs");
```

This analysis and the RL-GLM-HMM models use the same optogenetic experiment. The current Figure 4 implementation reads inputs prepared within `rl_glm_hmm/`: the first two CSVs come from `rl_glm_hmm/scripts/prepare_rlhmm_table.m`. The selected final fit generates `rlhmm_posteriors.csv`; only its trial keys and `Y` are used to retain the historical analysis cohort, not its state probabilities. The roster supplies cell type and inhibition/control assignment. The required input columns and exact phase mapping are documented in the MATLAB function header. These trial-level inputs are not currently included in the repository.

The function exports an animal-by-phase table, group means/SEMs, the binomial mixed-model interaction results, and PNG/PDF panels. For C/E, it plots P(pull) after a rewarded push (the complement of the source script's push win-stay). For D/F, it plots P(pull) after a rewarded pull. The tests preserve the source script's separate consecutive-trial filtering and random intercepts for animal and session. 

MATLAB's Statistics and Machine Learning Toolbox is required for `fitglme`.
