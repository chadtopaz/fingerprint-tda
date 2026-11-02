# Install the R package dependencies of the fingerprint-TDA replication code (release 1.1.0).
# Run once from the repository root:  source("DEPENDENCIES.R")
#
# R >= 4.3. Every reported run used R 4.5.3 on macOS (arm64); the pipeline fixes all seeds and
# iteration counts, and parallel feature extraction uses future.apply with future.seed = TRUE
# (L'Ecuyer-CMRG streams), so results are independent of the number of workers.

cran <- c("tidyverse", "readr", "jsonlite", "png", "TDA", "TDAvec", "digest",
          "future", "future.apply", "matrixStats", "pROC", "progressr")
to_install <- setdiff(cran, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

# EBImage is on Bioconductor (image I/O and Euclidean distance transforms).
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("EBImage", quietly = TRUE)) BiocManager::install("EBImage")

# Optional: the cubical-construction comparison (Code_cubical/) reaches Python GUDHI through reticulate.
if (!requireNamespace("reticulate", quietly = TRUE)) install.packages("reticulate")

# Persistent homology of the reported (triangulated) construction uses the GUDHI library bundled
# inside the TDA package: no separate installation. Python environments: see README.md, "Software".

cat("Dependencies installed.\n")
cat("Record your exact environment with sessionInfo() for the replication record.\n")
