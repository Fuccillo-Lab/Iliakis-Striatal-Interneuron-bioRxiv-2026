# Out-of-task optogenetics and photometry

This folder contains MATLAB analyses of photometry recordings collected
during optogenetic stimulation outside the behavioral task.

These data contribute to:

- **Figure 3D-G:** effects of TH-interneuron excitation on SST, D1-SPN, and
  A2A-SPN activity
- **Figure S12B:** validation of eNpHR3.0-mediated inhibition in SST and TH
  interneurons

Panel 3D is an experimental schematic. The functions here regenerate the
data-driven panels 3E-G and S12B.

## Repository contents

```text
opto_photo_out_of_task/
├── README.md
├── plot_figure3_out_of_task.m
├── plot_figureS12_panelB.m
└── .gitignore
```

## Data availability

The processed photometry data are distributed with the accompanying Zenodo
data release rather than through GitHub because the HDF5 file is approximately
1 GB.

Required files:

- `stream_key_build_clean_02.csv`
- `streams_build_clean_02.h5`

Zenodo record: **[link to be added]**

Place both files in the same local data directory. The functions read the HDF5
file in place and do not modify it.

Each `/streams/<stream_id>/` group contains processed signals and event times,
including:

- `timestamps`
- `ttl_timestamps`
- `dff_debleached_artcorr`
- `dff_debleached_artcorr_z`

The processed dF/F signal was generated using 5-Hz low-pass filtering,
continuous-light artifact correction, debleaching, and robust 415-to-470
regression. The `_z` dataset was additionally Z-scored across the full stream.

## Run in MATLAB

```matlab
addpath("path/to/Iliakis-Striatal-Interneuron-bioRxiv-2026/opto_photo_out_of_task")

dataFolder = "path/to/downloaded/OptoPhotoOutOfTask";
outputFolder = fullfile(dataFolder,"figure_outputs");

results3 = plot_figure3_out_of_task( ...
    fullfile(dataFolder,"stream_key_build_clean_02.csv"), ...
    fullfile(dataFolder,"streams_build_clean_02.h5"), ...
    outputFolder);

resultsS12B = plot_figureS12_panelB( ...
    fullfile(dataFolder,"stream_key_build_clean_02.csv"), ...
    fullfile(dataFolder,"streams_build_clean_02.h5"), ...
    outputFolder);
```

If the HDF5 file has a different local filename, pass its exact path as the
second argument.

## Figure 3E-G

`plot_figure3_out_of_task.m` compares 5-mW TH-ChrimsonR stimulation with
opsin-negative controls while recording SST, D1, or A2A populations.

The default analysis:

1. Selects successfully imported, usable, 20-Hz pulsed-stimulation recordings
   at 5 mW (1-s train; 20-ms pulse width).
2. Selects TH-ChrimsonR animals and opsin-negative controls matched by recorded
   cell type. Opsin-negative controls are pooled across recombinase backgrounds,
   matching the original analysis and manuscript sample sizes.
3. Aligns `dff_debleached_artcorr` to stimulation onset from -2.5 to +5 s.
4. Baseline-centers each epoch using -2 to 0 s.
5. Takes the mean across epochs within each stream and then across bilateral
   streams within each animal.
6. Plots the animal-level mean +/- SEM from -1 to +3 s.
7. Tests excitation against control from 0 to 1.5 s using a two-sided,
   animal-level cluster-mass permutation test with 5,000 label permutations.

The stream key produces the manuscript sample sizes:

- SST recording: TH opsin `N = 8`; control `N = 4`
- D1 recording: TH opsin `N = 7`; control `N = 8`
- A2A recording: TH opsin `N = 6`; control `N = 8`

The sample-size labels are calculated from the included HDF5 data and are not
hard-coded in the figure.

Figure 3 outputs:

- `figure3_out_of_task_photometry.pdf`
- `figure3_out_of_task_photometry.png`
- `figure3_out_of_task_photometry_animal_traces.csv`
- `figure3_out_of_task_photometry_animal_responses.csv`
- `figure3_out_of_task_photometry_group_summary.csv`
- `figure3_out_of_task_photometry_stream_summary.csv`
- `figure3_out_of_task_photometry_cluster_statistics.csv`

The permutation test is implemented directly in the plotting function and
uses a fixed random seed by default. It does not require the Statistics and
Machine Learning Toolbox.

## Figure S12B

`plot_figureS12_panelB.m` regenerates the stimulation-aligned
SST-eNpHR3.0/SST-GCaMP and TH-eNpHR3.0/TH-GCaMP traces.

The default analysis:

1. Selects usable, continuous-stimulation eNpHR3.0 recordings at 10 mW in
   which the manipulated and recorded cell types match.
2. Aligns `dff_debleached_artcorr_z` to stimulation onset from -2.5 to +5 s.
3. Baseline-centers each epoch using -2 to 0 s.
4. Takes the median across epochs within each stream and then across bilateral
   streams within each animal.
5. Plots the animal-level mean +/- SEM.

The available processed dataset contains:

- SST-eNpHR3.0/SST-GCaMP: `N = 3` animals
- TH-eNpHR3.0/TH-GCaMP: `N = 2` animals

Figure S12B outputs:

- `figureS12_panelB_opto_photometry.pdf`
- `figureS12_panelB_opto_photometry.png`
- `figureS12_panelB_opto_photometry_animal_traces.csv`
- `figureS12_panelB_opto_photometry_group_summary.csv`
- `figureS12_panelB_opto_photometry_stream_summary.csv`

## Software requirements

The plotting functions use base MATLAB HDF5, table, mathematical, and graphics
functions. No additional MATLAB toolbox is required.

## Large files

The `.gitignore` prevents processed HDF5 files and generated outputs from being
committed to GitHub. These files should remain in the external data directory
distributed through Zenodo.
