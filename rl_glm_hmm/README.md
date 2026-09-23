# RL-GLM-HMM analysis



This directory contains the RL-GLM-HMM implementation used for the behavioral analyses in Iliakis et al. (2026).



## Contents



* `rl_glm_hmm.py`: Project-specific `RLHMMHybrid` implementation.

* `upstream_glmhmm/`: Pinned Git submodule containing the upstream GLM-HMM infrastructure developed by Iris Stone.

* `scripts/`: Data preparation, model fitting, model evaluation, and figure-generation scripts associated with the paper. These scripts will be documented as they are added.

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



The upstream `env.yml` describes the historical software environment used for the original GLM-HMM repository. It does not define the complete environment used for the analyses in Iliakis et al.



A project-specific `requirements.txt` will specify the dependencies and versions required to reproduce the present analyses. The local virtual-environment directory itself is not included in the repository.





The code under `upstream_glmhmm/` remains a separate upstream repository and is not relicensed by this project. Any license subsequently applied to the paper repository applies only to original project-specific materials unless explicitly stated otherwise.



