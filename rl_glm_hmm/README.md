# RL-GLM-HMM analysis



This directory contains the RL-GLM-HMM implementation used for the behavioral analyses in Iliakis et al. (2026).



## Contents



* `rl_glm_hmm.py`: Project-specific `RLHMMHybrid` implementation.

* `upstream_glmhmm/`: Pinned Git submodule containing the upstream GLM-HMM infrastructure developed by Iris Stone.

* `scripts/`: Data preparation, model fitting, cross-validation, and figure-generation scripts:
  - `prepare_rlhmm_table.m`: Constructs the analysis-ready model input and Figure 5 trial metadata.
  - `cross_validate_models.py`: Performs session-fold model cross-validation.
  - `modelEval_QBackbone_population_sessionCV.m`: Fits Q, QF, and QDF baselines with the saved session folds.
  - `build_eval_DQPolicyGLM_QF_v2.m`: Constructs the fixed-Î”Q policy-GLM design and evaluates M1â€“M6.
  - `run_online_rlglm_M5_M6_cv.py`: Fits the online RL-GLM M5/M6 benchmarks with the saved session folds.
  - `assemble_s15c_masterll.py`: Selects and validates the six animal-level model scores plotted in Figure S15C.
  - `fit_final_model.py`: Fits the selected model specification using multiple random initializations.
  - `plot_rlhmm_figure5.m`: Reproduces the computational panels and statistics for Figure 5.
  - `plot_rlhmm_figureS15.m`: Reproduces the computational panels and statistics for Figure S15.

* `requirements.txt`: Python dependencies required for the project-specific implementation and analysis scripts.

* `data/`: Archived session-fold assignments, animal-level cross-validation summaries for Figure S15C, and the full historical `masterLL.csv` used for the submitted plot. Trial-level input data and posterior probabilities are not included here.



## Code provenance



This implementation builds on the GLM-HMM and HMM inference infrastructure developed by Iris Stone:



https://github.com/irisstone/glmhmm



The upstream repository is included as the Git submodule `upstream\_glmhmm` and is pinned to the following commit:



```text

860c30d68155401381552e48e19c79bb78b5a0ea

```



The files within `upstream_glmhmm/` originate from Iris Stoneâ€™s repository and are not original code from this project. No substantive changes were made to the pinned upstream files.



The project-specific implementation is contained in `rl_glm_hmm.py`. It subclasses the upstream `HMM` class and reuses its forward-backward inference, initial-state update, and parameter-initialization routines. The expectation-maximization procedure was adapted from the structure of the upstream GLM-HMM fitting implementation.



Project-specific extensions include:



* state-specific Q-learning traces;

* state-specific learning rates;

* bounded learning-rate estimation;

* separate policy coefficients for positive and negative action-value differences;

* additional choice- and outcome-history predictors;

* explicit handling of session boundaries;

* sessionwise transition-matrix updates;

* bounded and regularized parameter optimization; and

* optional Numba acceleration for model fitting.

## Model fitting and initialization

Models were evaluated using session-fold cross-validation, with complete sessions assigned to either the training or held-out dataset. To reduce sensitivity to local optima, each model was fitted from multiple parameter initializations. Within each cross-validation fold, the initialization with the highest final log-likelihood was selected for held-out scoring.

The archived `data/session_folds_allopto.csv` contains the actual five-fold assignment for 942 sessions from 38 animals. Supply it as `--fold-csv` to the RL-GLM-HMM cross-validation runner when the trial table has no fold column. The supplied historical `make_session_folds_perAnimal.m` version targets other datasets and does not reproduce this assignment; use the archived CSV for the submitted comparison.

The three baseline runners use the same model-ready trial table (created by `prepare_rlhmm_table.m`) and the saved folds. To regenerate their animal-level summaries from that table, set the trial CSV path and run the MATLAB functions from the repository root:

```matlab
addpath("rl_glm_hmm/scripts");
trialCsv = "rl_glm_hmm/data/modelTfinal_for_rlhmm_allopto.csv";
foldCsv = "rl_glm_hmm/data/session_folds_allopto.csv";
modelEval_QBackbone_population_sessionCV(trialCsv, foldCsv, ...
    "rl_glm_hmm/outputs/QBackbone_population_sessionCV_allopto");
build_eval_DQPolicyGLM_QF_v2(trialCsv, foldCsv, ...
    "rl_glm_hmm/outputs/DQPolicyGLM_QF_allopto");
```

The Q-backbone runner requires MATLAB Optimization Toolbox. The fixed-Î”Q runner also requires Statistics and Machine Learning Toolbox. Their primary outputs for panel C are `cv_animalModelSummary.csv` and `cv_animalModelSummary_DQPolicyGLM.csv`, respectively. The prepared trial CSV is not included in this repository and must be created from the source trial exports and calendar metadata.

The online RL-GLM runner uses the Python environment above:

```bash
python rl_glm_hmm/scripts/run_online_rlglm_M5_M6_cv.py \
  --csv rl_glm_hmm/data/modelTfinal_for_rlhmm_allopto.csv \
  --fold-csv rl_glm_hmm/data/session_folds_allopto.csv \
  --out-dir rl_glm_hmm/outputs/onlineRLGLM
```

Its archived output is `data/cv_animalModelSummary_onlineRLGLM.csv`. The RL-GLM-HMM runner is `scripts/cross_validate_models.py`; run `--help` for available fitting settings. A K3 run with its default two-element warm-start vectors needs `--no-warm-start` or three-element `--warm-alpha` and `--warm-beta` values. The original fold-fit configuration files were not available for this release, so the precise historical initialization options for the archived K2/K3 results have not been verified. The archived animal-level summaries reproduce the submitted panel C scores independently of rerunning optimization.

The selected three-state model specification was subsequently fitted to the complete dataset using ten random initializations. The initialization with the highest final log-likelihood was used for the reported posterior state probabilities and downstream analyses.

### Figure S15C score assembly

The QF, fixed-Î”Q M5/M6 GLMs, online M6 RL-GLM, and two/three-state RL-GLM-HMM results are first summarized as one row per animal and model, pooling held-out log likelihood and trial counts across folds. `data/masterLL.csv` archives the complete historical 16-model comparison. The following command reconstructs its six rows per animal used in panel C directly from the five archived source summaries:

```bash
python rl_glm_hmm/scripts/assemble_s15c_masterll.py \
  --qf rl_glm_hmm/data/cv_animalModelSummary_QBackbone.csv \
  --glm rl_glm_hmm/data/cv_animalModelSummary_DQPolicyGLM.csv \
  --online rl_glm_hmm/data/cv_animalModelSummary_onlineRLGLM.csv \
  --k2 rl_glm_hmm/data/cv_animalModelSummary_RLHMM_K2.csv \
  --k3 rl_glm_hmm/data/cv_animalModelSummary_RLHMM_K3.csv \
  --reference rl_glm_hmm/data/masterLL.csv \
  --out rl_glm_hmm/outputs/masterLL_panelC.csv
```

Run this from the repository root. `--reference` checks the output against the archived submitted scores and can be omitted when assembling newly generated results. The output contains 228 animal-by-model rows; it can replace `masterLL.csv` in the MATLAB S15 plotting call below. It does not refit models. The fixed-Î”Q GLM analysis fit its QF backbone and constructed the Î”Q and animal-bias predictors using the full dataset before cross-validating GLM coefficients; this historical convention is retained in the archived scores. Its models score 121,830 trials, while QF and the online/RL-GLM-HMM models score 123,714 trials.

## Cloning the repository



Because the upstream dependency is stored as a Git submodule, clone the complete paper repository using:



```bash

git clone --recurse-submodules https://github.com/Fuccillo-Lab/Iliakis-Striatal-Interneuron-bioRxiv-2026.git

```



If the repository has already been cloned without its submodules, initialize them from the repository root using:



```bash

git submodule update --init --recursive

```



The project-specific implementation locates `upstream\_glmhmm` relative to its own file location, so no manual modification of the Python path should be required.



## Software environment



The project-specific RL-GLM-HMM implementation was tested using Python 3.12.10. Its direct Python dependencies are recorded in `requirements.txt`.

Install the dependencies using:

```bash
python -m pip install -r rl_glm_hmm/requirements.txt
```
The upstream env.yml describes the historical environment associated with the original GLM-HMM repository. It is retained as part of the pinned upstream submodule but does not define the environment used for the analyses in Iliakis et al.

Virtual-environment directories such as .venv/ are not included in the repository.


The code under `upstream_glmhmm/` remains a separate upstream repository and is not relicensed by this project. Any license subsequently applied to the paper repository applies only to original project-specific materials unless explicitly stated otherwise.


## Reproducing the figures

Run the MATLAB plotting functions from the root of the paper repository. Add the analysis scripts to the MATLAB path:

```matlab
addpath("rl_glm_hmm/scripts");
```

Set `dataDir` to the directory containing the analysis-ready CSV files:

```matlab
dataDir = "rl_glm_hmm/data";
```

### Main Figure 5

```matlab
results5 = plot_rlhmm_figure5( ...
    fullfile(dataDir, "figure5_trial_metadata.csv"), ...
    fullfile(dataDir, "rlhmm_posteriors.csv"), ...
    fullfile(dataDir, "rl_params.csv"), ...
    fullfile(dataDir, "transition_matrix.csv"), ...
    fullfile(dataDir, "animals.csv"), ...
    "rl_glm_hmm/outputs/figure5");
```

This function reproduces the analytical panels in Figure 5 and writes figure source files and supporting statistics tables to the specified output directory.

### Supplemental Figure S15

```matlab
resultsS15 = plot_rlhmm_figureS15( ...
    fullfile(dataDir, "masterLL.csv"), ...
    fullfile(dataDir, "rlhmm_posteriors.csv"), ...
    fullfile(dataDir, "animals.csv"), ...
    "rl_glm_hmm/outputs/figureS15");
```

This function reproduces the data-driven panels Câ€“L of Figure S15. Panels A and B are conceptual schematics and are not generated programmatically. The six example sessions in panels Dâ€“I are identified explicitly in the function to reproduce the submitted figure.

The final output-directory argument is optional. If it is omitted, MATLAB opens the figures and returns the numerical results without writing files to disk.

Detailed descriptions of the required inputs, historical analysis conventions, outputs, and software dependencies are available from MATLAB:

```matlab
help plot_rlhmm_figure5
help plot_rlhmm_figureS15
```

The Statistics and Machine Learning Toolbox is required for the statistical tests and mixed-effects models used by these functions.
