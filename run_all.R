#!/usr/bin/env Rscript
# Run from repository root.
source(file.path("R","reproduce_manuscript.R"))
source(file.path("R","make_variance_figure.R"))
source(file.path("R","make_longitudinal_figure.R"))
cat("All analyses and figure scripts completed.\n")
