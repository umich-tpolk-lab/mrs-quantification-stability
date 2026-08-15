# Manuscript target values for release verification

These are **validation targets only**. They are not read by the analysis code.

## QC
- Available spectra: 668
- Excluded GABA+: 51
- Excluded Glx: 53

## Cross-ROI mean ICC(C,1)
| Metabolite | UNC | CSF | TC | ATC | CR |
|---|---:|---:|---:|---:|---:|
| GABA+ | 0.535 | 0.314 | 0.293 | 0.301 | 0.398 |
| Glx | 0.482 | 0.227 | 0.207 | 0.193 | 0.378 |

## Variance decomposition (cross-ROI mean CV%)
- GABA+ UNC: between 10.6, within 9.9
- GABA+ ATC: between 6.3, within 9.8
- Glx UNC: between 9.8, within 10.1
- Glx ATC: between 4.7, within 9.6

## Tissue-residualized ICC
- GABA+: UNC 0.535 → tissue-residualized 0.264; ATC 0.301
- Glx: UNC 0.482 → tissue-residualized 0.160; ATC 0.193

## Longitudinal percent change
GABA+ (AUD, SM, VV):
- UNC: -13.5%, -16.9%, -16.5%
- CR: -6.0%, -6.2%, -9.4%

Glx (AUD, SM, VV):
- UNC: -7.1%, -1.8%, -6.9%
- CR: +0.8%, +11.5%, +1.0%

## Meta-analysis values reported in the manuscript
These should be checked after running the new `metafor` section locally:
- UNC–CSF pooled ΔICC: 0.221 [0.151, 0.290]
- UNC–ATC pooled ΔICC: 0.253 [0.165, 0.341]
- Multivariate UNC–CSF: 0.217 [0.154, 0.280]
- Multivariate UNC–ATC: 0.261 [0.188, 0.333]

## Participant demographics
- Wave 1 age: 70.1 ± 4.6 years, range 65–82
- Wave 2 age: 74.1 ± 4.6 years, range 68–86
- Inter-session interval: 4.0 ± 0.6 years, range 3.1–5.1
