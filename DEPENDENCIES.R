# Install R package dependencies for the Fingerprint TDA replication.
# Run once from the repository root:  source("DEPENDENCIES.R")

cran <- c("tidyverse", "TDA", "TDAvec", "digest",
          "future", "future.apply", "matrixStats", "pROC", "progressr")
to_install <- setdiff(cran, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

# EBImage is on Bioconductor (image I/O and Euclidean distance transforms).
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("EBImage", quietly = TRUE)) BiocManager::install("EBImage")

# Persistent homology uses the GUDHI library, which ships inside the TDA
# package -- no separate installation is required.

cat("Dependencies installed.\n")
cat("Tip: record your exact environment with sessionInfo() for the replication record.\n")
