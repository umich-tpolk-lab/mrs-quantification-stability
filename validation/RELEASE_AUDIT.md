# Release audit

This draft package was constructed from the original files supplied for the manuscript.

## Data minimization performed

- Restricted the master data from 420 source rows to the **56 manuscript participants and waves 1–2 only (112 rows)**.
- Removed all original `Subject` and `SubNum` identifiers.
- Generated new public identifiers (`sub-001` … `sub-056`) using a randomized mapping. The mapping is **not included** in the release.
- Removed exact ages, sex, names, RA fields, visit dates/times, recruitment variables, clinical/drug variables, internal paths, raw filenames, and all unrelated master-sheet variables.
- Retained only the five quantification strategies in the manuscript (UNC, CSF, TC, ATC, CR), the six hemisphere-region voxels, and the required tissue fractions.
- Reduced the QC file to scientific QC metrics and public participant/session/location keys.

## Structural validation

- Participants: **56**
- Participant-session rows: **112**
- Available spectrum-level QC rows: **668**
- GABA+ exclusions at 15 Hz + 8-MAD rule: **51**
- Glx exclusions at 15 Hz + 8-MAD rule: **53**
- Short-interval subset (≤4.02 years): **28 participants**
- Mean interval in that subset: **3.51 years**

## Primary ICC check

A Python implementation of the same ICC(C,1) and QC logic was used as an independent structural check. It reproduces the manuscript's point estimates, including the cross-ROI means:

- GABA+: UNC 0.535, CSF 0.314, TC 0.293, ATC 0.301, CR 0.398.
- Glx: UNC 0.482, CSF 0.227, TC 0.207, ATC 0.193, CR 0.378.

The detailed check is in `primary_icc_point_estimates_python_check.csv`.

## Items that still require author-side verification before public release

1. **Run the R package locally.** R is not installed in the environment that built this draft, so the R files could not be executed here. The core formulas were independently checked in Python, but the final public release should be run in R 4.6.0.
2. **Meta-analysis.** The supplied development script did not contain the `metafor` meta-analysis reported in the current manuscript. A transparent implementation has been added to `R/reproduce_manuscript.R`; its exact output should be checked against the manuscript values.
3. **Figure significance stars.** The supplied longitudinal figure script contains hard-coded raw p-values. The manuscript methods state that p-values were Bonferroni-corrected within families. Confirm whether the figure should display raw or adjusted significance.
4. **Demographic fields.** Exact age and sex remain deliberately excluded from the public data because they are unnecessary for the primary analyses. Using the exact 56-participant analysis roster, the source master sheet gives Wave 1 age 70.1 ± 4.6 years (range 65–82), Wave 2 age 74.1 ± 4.6 years (range 68–86), and an inter-session interval of 4.0 ± 0.6 years (range 3.1–5.1; rounded to one decimal). These aggregate values have been used to correct the manuscript.
5. **IRB/consent approval.** Confirm that unrestricted public sharing of this de-identified individual-level derived dataset is permitted.
6. **Licenses and repository metadata.** Choose code/data licenses and replace the GitHub placeholder in `CITATION.cff`.
