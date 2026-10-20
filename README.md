# Remote Photoplethysmography (rPPG) — Signal Processing

Signal-processing work on **remote photoplethysmography**: recovering the blood-volume pulse from facial video, using a 4-D signal (height × width × RGB channel × time) recorded from each person.

## What was done

1. **Pre-processing** — downsampling the 4-D data from 320×510 to 32×51 spatial locations and filtering all three colour channels.
2. **Pulse averaging** — for each location in the green channel, splitting the signal into individual pulses at the minima, aligning them and averaging to get a template pulse ("mask") per location, then a total mask per person.
3. **Validation** — comparing every location with its local mask and the total mask using correlation and convolution, producing maps of signal quality across the face.
4. **Time–frequency analysis** — frequency-domain plots and continuous wavelet transforms of the signals.

## Contents

| File | Description |
|---|---|
| `RPPG_project_week1.m` | 2-D signal plots, frequency domain, wavelet transform, first pulse averaging |
| `RPPG_project_week2.m` | Pulse detection with peak finding and average-pulse construction |
| `RPPG_project_week4.m` | Full 4-D pipeline: per-location filtering, averaging and validation |
| `week1.pptx`, `Week 2.pptm`, `Week 3.pptx` | Weekly progress slides |
| `Results.pptx` | Summary of the method and results |

Related: my bachelor's thesis on full-body skin segmentation for robust rPPG — [SkinMap](https://github.com/rosebaekfany/SkinMap-A-Novel-Weighted-Full-Body-Skin-Segmentation).

## Tools

MATLAB
