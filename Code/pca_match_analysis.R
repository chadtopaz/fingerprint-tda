############################################################
# PCA-MATCHED COMPARISON OF GLOBAL BETTI vs GLOBAL PI
#
# The two global topological baselines in the main pipeline achieve nearly
# identical AUC (~0.847) but differ at low FAR (e.g., TAR@1e-3: Betti 6.6%,
# PI 3.0%). This gap could plausibly reflect any of: representation
# (Betti curve vs persistence image), dimensionality (100 vs 800),
# distance metric (L^1 vs L^2), or kernel bandwidth.
#
# This script disentangles dimensionality + metric by projecting each of the
# Betti and PI vectors onto its own PCA basis estimated from training-fold
# data, then scoring pairs with L^2 distance on the matched-dim projections.
# Comparing AUC and TAR@FAR for the matched-dim methods partially reduces the
# dimensionality and metric confounds, though representation-specific choices
# such as smoothing and scaling remain.
#
# Usage (from project root):
#   source("Code/pca_match_analysis.R")
#
# The script sources Code/fingerprintTDA.R with the main entry point suppressed,
# loads the cached feature bank, and uses the same fold partition and pair
# enumeration as the main pipeline. Cross-validation results are written to
#   outputs/pca_match_fold_metrics.csv  (per-fold AUC, EER, TAR@FAR)
#   outputs/pca_match_summary.csv       (mean +/- SD across folds)
#
# Config (top of file) sets the PCA target dimension k.
############################################################

# Suppress fingerprintTDA.R's main entry point only while we source it; restore
# whatever the caller had before so a subsequent source("Code/fingerprintTDA.R")
# in the same R session runs the full pipeline as expected.
.old_fp_suppress <- getOption("fp.suppress_main")
options(fp.suppress_main = TRUE)
source("Code/fingerprintTDA.R")  # provides helpers, CONFIG, and the cached bank
options(fp.suppress_main = .old_fp_suppress)
rm(.old_fp_suppress)

# ----------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------

PCA_CFG <- list(
  k                 = 50,    # PCA target dimensionality (matched for both)
  output_fold_csv   = file.path(CONFIG$output_dir, "pca_match_fold_metrics.csv"),
  output_summary_csv = file.path(CONFIG$output_dir, "pca_match_summary.csv")
)

dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------
# Load cached feature bank
# ----------------------------------------------------------------

paths <- get_paths(CONFIG)
stopifnot(length(paths) > 0)
ids <- map_chr(paths, get_identity)
folds <- make_identity_folds(observed_ids = ids,
                             n_folds = CONFIG$n_folds,
                             seed = CONFIG$cv_seed)

cat("============================================\n")
cat("PCA-MATCHED Betti vs PI comparison\n")
cat("============================================\n\n")
cat("Images:", length(paths), "\n")
cat("Identities:", length(unique(ids)), "\n")
cat("Folds:", CONFIG$n_folds, "\n")
cat("PCA target dimension k:", PCA_CFG$k, "\n\n")

bank <- load_or_compute_bank(paths, CONFIG)

betti_mat <- bank_to_matrix(bank, "betti_global_100")   # n_images x 100
pi_mat    <- bank_to_matrix(bank, "pi_global")          # n_images x 800

stopifnot(nrow(betti_mat) == length(paths),
          nrow(pi_mat)    == length(paths))

# ----------------------------------------------------------------
# Per-fold scoring helpers
# ----------------------------------------------------------------

# Fit PCA on training rows of `mat`, project all rows to k components.
# Returns the n_images x k projected matrix.
pca_project_all <- function(mat, train_idx, k) {
  fit <- prcomp(mat[train_idx, , drop = FALSE], center = TRUE, scale. = FALSE)
  k_use <- min(k, ncol(fit$rotation))
  # Project ALL rows (train + test) using the training-derived center and rotation
  centered <- sweep(mat, 2, fit$center, FUN = "-")
  centered %*% fit$rotation[, seq_len(k_use), drop = FALSE]
}

# ----------------------------------------------------------------
# Cross-validation loop
# ----------------------------------------------------------------

per_fold_rows <- vector("list", length(folds))

for (fold_id in seq_along(folds)) {
  test_ids  <- folds[[fold_id]]
  train_ids <- setdiff(unique(ids), test_ids)
  train_idx <- which(ids %in% train_ids)
  test_idx  <- which(ids %in% test_ids)

  pair_train <- build_pair_mat_global(train_idx)
  pair_test  <- build_pair_mat_global(test_idx)
  y_train <- compute_y_from_pair_mat(pair_train, ids)
  y_test  <- compute_y_from_pair_mat(pair_test,  ids)

  cat(sprintf("Fold %d/%d:  train ids=%d, test ids=%d, train pairs=%d, test pairs=%d\n",
              fold_id, length(folds), length(train_ids), length(test_ids),
              ncol(pair_train), ncol(pair_test)))

  # PCA projection (training-fit, applied to all)
  betti_pca <- pca_project_all(betti_mat, train_idx, PCA_CFG$k)
  pi_pca    <- pca_project_all(pi_mat,    train_idx, PCA_CFG$k)

  # Score train + test pairs with L^2 on matched-dim PCA features
  s_betti_train <- score_pairs_l2(pair_train, betti_pca)
  s_betti_test  <- score_pairs_l2(pair_test,  betti_pca)
  s_pi_train    <- score_pairs_l2(pair_train, pi_pca)
  s_pi_test     <- score_pairs_l2(pair_test,  pi_pca)

  # Metrics
  res_betti <- eval_fold_method(sprintf("Betti_PCA%d_L2", PCA_CFG$k),
                                y_train, s_betti_train, y_test, s_betti_test,
                                CONFIG$fars_report)
  res_pi    <- eval_fold_method(sprintf("PI_PCA%d_L2", PCA_CFG$k),
                                y_train, s_pi_train, y_test, s_pi_test,
                                CONFIG$fars_report)

  for (res in list(res_betti, res_pi)) {
    tar_row <- res$tar_tbl %>%
      pivot_wider(names_from = FAR, values_from = TAR, names_prefix = "TAR_FAR_")
    per_fold_rows[[length(per_fold_rows) + 1]] <- tibble(
      fold   = fold_id,
      method = res$method,
      AUC    = res$auc,
      EER    = res$eer
    ) %>%
      bind_cols(tar_row)
  }

  # EER is stored as a decimal fraction; multiply by 100 for the percent label.
  cat(sprintf("    Betti_PCA%d_L2:  AUC=%.3f  EER=%.1f%%\n",
              PCA_CFG$k, res_betti$auc, 100 * res_betti$eer))
  cat(sprintf("    PI_PCA%d_L2:     AUC=%.3f  EER=%.1f%%\n",
              PCA_CFG$k, res_pi$auc, 100 * res_pi$eer))
}

# ----------------------------------------------------------------
# Aggregate and write
# ----------------------------------------------------------------

per_fold <- bind_rows(per_fold_rows)
write_csv(per_fold, PCA_CFG$output_fold_csv)

summary_tbl <- per_fold %>%
  group_by(method) %>%
  summarize(
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_sd   = sd(AUC, na.rm = TRUE),
    EER_mean = mean(EER, na.rm = TRUE),
    EER_sd   = sd(EER, na.rm = TRUE),
    across(starts_with("TAR_FAR_"),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd   = ~sd(.x,   na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups  = "drop"
  )
write_csv(summary_tbl, PCA_CFG$output_summary_csv)

cat("\n============================================\n")
cat("Summary across folds (mean +/- SD):\n")
cat("============================================\n")
print(summary_tbl)
cat("\nWrote:\n")
cat("  ", PCA_CFG$output_fold_csv, "\n", sep = "")
cat("  ", PCA_CFG$output_summary_csv, "\n", sep = "")
