# Anatomical targeting

This directory contains the MATLAB code used to generate the Allen-atlas
targeting panels in Figures S1, S3, S8, S10, and S12 of:

> Iliakis EA et al. *Striatal interneuron microcircuits gate reinforcement to
> stabilize adaptive choice.* bioRxiv (2026).
> [https://doi.org/10.64898/2026.09.21.753300](https://doi.org/10.64898/2026.09.21.753300)

The scripts reproduce the coordinate-based coronal panels. Histology images,
experimental schematics, panel lettering, and the final multi-panel layouts
were assembled separately for the manuscript.

## Contents

- `plot_anatomy_supplements.m` is the public entry point and generates all
  anatomy panels from the five animal metadata tables.
- `plot_allen_coronal_sections.m` is the shared renderer that maps bilateral
  histological coordinates onto Allen CCFv3 coronal sections.
- `.gitignore` prevents local NRRD atlas files from being committed.

## Data availability

The study metadata are not stored in this GitHub repository. They will be
archived separately on Zenodo with the other study data.

**Zenodo record:** forthcoming

The anatomy analysis expects these files from the Zenodo archive:

| File | Figures | Description |
| --- | --- | --- |
| `animalsMuscimol.csv` | S1 | Muscimol cannula placements |
| `animalsPhotometry.csv` | S3, S8 | PV, SST, TH, D1-SPN, and A2A-SPN photometry placements |
| `animalsOPInTask.csv` | S10 | In-task TH excitation with SST photometry |
| `animalsOPOutOfTask.csv` | S10 | Out-of-task TH excitation with SST or SPN photometry |
| `animalsOptoBehavior.csv` | S12 | SST and TH inhibition cohorts |

All tables contain histological fiber coordinates in millimeters relative to
bregma. The default schema uses `animalID`, `good_targeting`, and the bilateral
coordinate columns:

```text
left_histo_fiber_ap   left_histo_fiber_dv   left_histo_fiber_ml
right_histo_fiber_ap  right_histo_fiber_dv  right_histo_fiber_ml
```

The photometry tables additionally use `usable_data_left` and
`usable_data_right`. `animalsOptoBehavior.csv` identifies animals with
`animal_id`; the entry-point function maps this field to `animalID` internally.

## Allen CCFv3 template

The 10-Âµm Allen CCFv3 average template is an external dependency and is not
redistributed in this repository or in the Zenodo archive. Download
[`average_template_10.nrrd`](https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/average_template_10.nrrd)
from the
[Allen Institute download directory](https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/)
and provide its local path when running the analysis.

Please cite the atlas as:

> Wang Q, Ding S-L, Li Y, et al. The Allen Mouse Brain Common Coordinate
> Framework: A 3D Reference Atlas. *Cell*. 2020;181(4):936â€“953.e20.
> [https://doi.org/10.1016/j.cell.2020.04.007](https://doi.org/10.1016/j.cell.2020.04.007)

An annotation volume is not required for the manuscript panels. The lower-level
renderer accepts an optional `AnnotationFile` if anatomical boundaries are
desired.

## Usage

From the repository root in MATLAB:

```matlab
addpath("anatomy");

dataRoot = "path/to/Zenodo/anatomy";
atlasFile = "path/to/average_template_10.nrrd";
outputFolder = "path/to/output/folder";

results = plot_anatomy_supplements( ...
    fullfile(dataRoot, "animalsMuscimol.csv"), ...
    fullfile(dataRoot, "animalsPhotometry.csv"), ...
    fullfile(dataRoot, "animalsOPInTask.csv"), ...
    fullfile(dataRoot, "animalsOPOutOfTask.csv"), ...
    fullfile(dataRoot, "animalsOptoBehavior.csv"), ...
    atlasFile, ...
    outputFolder);
```

A function named `nrrdread` must be available on the MATLAB path. Recent
MATLAB releases provide it through Medical Imaging Toolbox; a compatible NRRD
reader may also be used. No Statistics and Machine Learning Toolbox functions
are used. Figures are exported with `exportgraphics` when available and with a
`print` fallback on older MATLAB releases.

## Panel selection

| Figure | Selection |
| --- | --- |
| S1A | Muscimol animals with `good_targeting` marked as included |
| S3Câ€“E | mDMS photometry animals with `photo_cell` equal to `pv`, `ltsi`, or `thin` |
| S8B, D | mDMS photometry animals with `photo_cell` equal to `d1` or `a2a` |
| S10B, E, G | Out-of-task animals with TH-Chrimson stimulation and SST, D1, or A2A photometry |
| S10C | In-task animals with TH-Chrimson stimulation and SST photometry |
| S12Câ€“D | SST or TH control/inhibition animals included in behavior or slice experiments |

For photometry panels, left- and right-hemisphere inclusion is determined by
`usable_data_left` and `usable_data_right`. Included placements are shown as
white circles connected by blue lines; excluded hemispheres are shown as gray
crosses. Bilateral placements are assigned together to the nearest AP panel
using the animal's mean AP coordinate. Panels are spaced by 0.30 mm around the
+0.75-mm anchor used in the manuscript.

The Allen volume is interpreted in AP Ã— DV Ã— ML orientation. Coordinates are
converted using the SHARP-Track/SHARCQ bregma convention of `[5.4, 0, 5.7]` mm
from the volume origin.

## Outputs

Each anatomy subpanel is exported separately as a vector PDF and a 600-dpi
PNG. For every subpanel, the function also writes:

- `*_sites.csv`, containing the plotted hemisphere-level coordinates, panel
  assignments, inclusion flags, and Allen voxel locations.
- `*_animals.csv`, containing one row per animal with its mean AP coordinate,
  assigned panel, number of plotted sites, and inclusion summary.

The returned `results` structure is organized by figure and panel and contains
the figure handles and both underlying tables.
