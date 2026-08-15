#!/usr/bin/env Rscript
# Run from repository root.
source(file.path("R","reproduce_manuscript.R"))
source(file.path("R","make_figure4_variance.R"))
source(file.path("R","make_figure5_longitudinal.R"))
cat("All analyses and figure scripts completed.\n")
