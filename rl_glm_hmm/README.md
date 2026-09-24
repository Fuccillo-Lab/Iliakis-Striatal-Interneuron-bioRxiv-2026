# RL-GLM-HMM analysis



This directory contains the RL-GLM-HMM implementation used for the behavioral analyses in Iliakis et al. (2026).



## Contents



* `rl_glm_hmm.py`: Project-specific `RLHMMHybrid` implementation.

* `upstream_glmhmm/`: Pinned Git submodule containing the upstream GLM-HMM infrastructure developed by Iris Stone.

* `scripts/`: Data preparation, model fitting, cross-validation, and figure-generation scripts:
  - `prepare_rlhmm_table.m`: Constructs the analysis-ready model input and Figure 5 trial metadata.
  - `cross_validate_models.py`: Performs session-fold model cross-validation.
  - `fit_final_model.py`: Fits the selected model specification using multiple random initializations.
  - `plot_rlhmm_figure5.m`: Reproduces the computational panels and statistics for Figure 5.
  - `plot_rlhmm_figureS15.m`: Reproduces the computational panels and statistics for Figure S15.

* `requirements.txt`: Python dependencies required for the project-specific implementation and analysis scripts.



## Code provenance



This implementation builds on the GLM-HMM and HMM inference infrastructure developed by Iris Stone:



https://github.com/irisstone/glmhmm



The upstream repository is included as the Git submodule `upstream\_glmhmm` and is pinned to the following commit:



```text

860c30d68155401381552e48e19c79bb78b5a0ea

```



The files within `upstream_glmhmm/` originate from Iris Stone’s repository and are not original code from this project. No substantive changes were made to the pinned upstream files.



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

The selected three-state model specification was subsequently fitted to the complete dataset using ten random initializations. The initialization with the highest final log-likelihood was used for the reported posterior state probabilities and downstream analyses.

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

This function reproduces the data-driven panels C–L of Figure S15. Panels A and B are conceptual schematics and are not generated programmatically. The six example sessions in panels D–I are identified explicitly in the function to reproduce the submitted figure.

The final output-directory argument is optional. If it is omitted, MATLAB opens the figures and returns the numerical results without writing files to disk.

Detailed descriptions of the required inputs, historical analysis conventions, outputs, and software dependencies are available from MATLAB:

```matlab
help plot_rlhmm_figure5
help plot_rlhmm_figureS15
```

The Statistics and Machine Learning Toolbox is required for the statistical tests and mixed-effects models used by these functions.

