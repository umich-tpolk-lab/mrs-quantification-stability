# MRS quantification, longitudinal stability, and change

Reproducibility package for **“How quantification choices shape MRS measures of GABA and Glx: Longitudinal stability and change in older adults.”**

## What this release contains

This repository contains a manuscript-specific, de-identified extract of the derived MRS measurements and spectrum-level QC metrics needed to reproduce the reported analyses. The release is restricted to the 56 participants and two sessions analyzed in the paper. Public participant IDs were newly generated for this release and are not the study's internal IDs.

The original master datasets contained hundreds of additional records and many unrelated variables. Those data are intentionally not included.

## Quick start

From the repository root:

```r
install.packages("metafor")   # only external package required for the meta-analysis
source("run_all.R")
```

The main numerical results are written to `outputs/` as CSV files. The two supplied figure scripts also write the variance-decomposition and longitudinal-change figures there.

The analysis uses 10,000 bootstrap resamples and a fixed random seed.

## Reproduced analyses

`R/reproduce_manuscript.R` implements:

1. Spectrum-level QC using creatine linewidth >15 Hz plus the >8-MAD non-physical-fit rule.
2. Bilateral measures from QC-passing hemispheres.
3. ICC(C,1) longitudinal stability with 10,000-resample bootstrap confidence intervals.
4. Paired bootstrap comparisons among quantification strategies.
5. Between- and within-subject variance decomposition.
6. Variance explained by voxel tissue composition and the CSF/tGM axes.
7. Tissue-residualized ICCs.
8. Longitudinal metabolite change and tissue-fraction change.
9. Bonferroni-adjusted p-values for the longitudinal-change families.
10. Linewidth-threshold sensitivity (12–20 Hz).
11. Shorter-interval sensitivity using the ≤4.02-year half of the sample.
12. The reported random-effects/meta-regression and multivariate meta-analysis, if `metafor` is installed.

## Expected structural checks

The release should contain:

- 56 participants.
- 112 participant-session rows.
- 668 available spectra.
- At the 15-Hz linewidth threshold, 51 GABA+ and 53 Glx spectra excluded.
- Cross-ROI mean ICCs at 15 Hz approximately:
  - GABA+: UNC 0.535, CSF 0.314, TC 0.293, ATC 0.301, CR 0.398.
  - Glx: UNC 0.482, CSF 0.227, TC 0.207, ATC 0.193, CR 0.378.

## Files

- `data/mrs_measurements.csv` — de-identified, manuscript-only hemisphere-level derived MRS/tissue data.
- `data/gannet_qc_metrics.csv` — de-identified spectrum-level QC metrics.
- `data/participants.csv` — public participant ID and rounded inter-session interval only.
- `R/reproduce_manuscript.R` — primary numerical analysis.
- `R/make_figure4_variance.R` — supplied variance figure script.
- `R/make_figure5_longitudinal.R` — supplied longitudinal figure script.
- `run_all.R` — one-command entry point.
