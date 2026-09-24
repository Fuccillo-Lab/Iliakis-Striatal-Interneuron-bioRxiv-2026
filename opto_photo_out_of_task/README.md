# Out-of-task optogenetics and photometry

This folder contains MATLAB analyses of photometry recordings collected during optogenetic stimulation outside the behavioral task.

These data contribute to:

- **Figure 3D–G:** effects of TH-interneuron excitation on SST, D1-SPN, and A2A-SPN activity
- **Figure S12B:** validation of eNpHR3.0-mediated inhibition in SST and TH interneurons

## Repository contents

```text
opto_photo_out_of_task/
├── README.md
├── plot_figureS12_panelB.m
└── .gitignore
```

Code for the Figure 3 out-of-task excitation analysis will be added separately.

## Data availability

The processed photometry data are distributed with the accompanying Zenodo data release rather than through GitHub because the HDF5 file is approximately 1 GB.

Required files:

- `stream_key_build_clean_02.csv`
- `streams_build_clean_02.h5`

Zenodo record: **[link to be added]**

Place both files in the same local data directory. The MATLAB functions read the HDF5 file in place and do not modify it.

Each `/streams/<stream_id>/` group in the HDF5 file contains the processed signals and metadata required for event alignment, including:

- `timestamps`
- `ttl_timestamps`
- `dff_debleached_artcorr_z`

The processed dF/F signal was generated using 5-Hz low-pass filtering, correction of the continuous-light artifact, debleaching, robust 415-to-470 regression, and whole-stream Z-scoring.

## Figure S12B

`plot_figureS12_panelB.m` regenerates the stimulation-aligned SST-eNpHR3.0/SST-GCaMP and TH-eNpHR3.0/TH-GCaMP traces shown in Figure S12B.

### Run in MATLAB

```matlab
addpath("path/to/Iliakis-Striatal-Interneuron-bioRxiv-2026/opto_photo_out_of_task")

dataFolder = "path/to/downloaded/OptoPhotoOutOfTask";
outputFolder = fullfile(dataFolder,"figure_outputs");

resultsS12B = plot_figureS12_panelB( ...
    fullfile(dataFolder,"stream_key_build_clean_02.csv"), ...
    fullfile(dataFolder,"streams_build_clean_02.h5"), ...
    outputFolder);
```

If the HDF5 file has a different local filename, pass its exact path as the second argument.

### Analysis

The default analysis:

1. Selects successfully imported, usable, continuous-stimulation eNpHR3.0 recordings in which the optogenetically manipulated and recorded cell types match.
2. Selects the 10 mW stimulation condition.
3. Aligns the Z-scored dF/F signal to stimulation onset from −2.5 to +5 s.
4. Baseline-centers each epoch using the −2 to 0 s interval.
5. Takes the median across stimulation epochs within each stream.
6. Collapses bilateral streams within each animal.
7. Plots the animal-level mean ± SEM.

Each animal therefore contributes equally to the group trace.

The available processed dataset contains:

- SST-eNpHR3.0/SST-GCaMP: **N = 3 animals**
- TH-eNpHR3.0/TH-GCaMP: **N = 2 animals**



### Optional parameters

```matlab
resultsS12B = plot_figureS12_panelB( ...
    keyFile,h5File,outputFolder, ...
    'PowerMw',10, ...
    'SignalName',"dff_debleached_artcorr_z", ...
    'EpochWindow',[-2.5 5], ...
    'BaselineWindow',[-2 0], ...
    'StimWindow',[0 1], ...
    'EpochDt',0.025, ...
    'Aggregation',"median");
```

### Outputs

The function exports:

- `figureS12_panelB_opto_photometry.pdf`
- `figureS12_panelB_opto_photometry.png`
- `figureS12_panelB_opto_photometry_animal_traces.csv`
- `figureS12_panelB_opto_photometry_group_summary.csv`
- `figureS12_panelB_opto_photometry_stream_summary.csv`

The CSV files record the animal-level traces, plotted group summaries, and stream-selection provenance.

## Figure 3D–G

Figure 3D–G uses related out-of-task recordings to measure the effects of TH-interneuron excitation on:

- SST interneurons
- D1 dopamine-receptor-expressing SPNs
- A2A dopamine-receptor-expressing SPNs

The corresponding plotting function and detailed analysis documentation will be added to this folder.

## Software requirements

The Figure S12B plotting function uses base MATLAB HDF5, table, and graphics functions. It does not require the Statistics and Machine Learning Toolbox.
