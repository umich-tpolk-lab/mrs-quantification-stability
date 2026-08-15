# Changes made after checking the current manuscript

- Figure 5 significance markers now use the Bonferroni-adjusted p-values specified in Methods.
- `make_longitudinal_figure.R` now reads the numerical output files instead of hard-coding figure values.
- `make_variance_figure.R` now reads the numerical output files instead of hard-coding Figure 4 values.
- The primary analysis remains configured for 10,000 bootstrap iterations.
- The meta-analysis implementation remains flagged for local R verification against the manuscript targets.
- Participant demographics were recalculated from the exact 56-participant analysis roster: Wave 1 age 70.1 ± 4.6 years (65–82), Wave 2 age 74.1 ± 4.6 years (68–86), and inter-session interval 4.0 ± 0.6 years (3.1–5.1). Exact ages remain excluded from the public data because they are not required for reproduction.
