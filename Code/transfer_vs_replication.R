############################################################
# STEP 1.1 (C1, transfer) — TRANSFER minus REPLICATION paired differences on DB3, from the
# three bootstrap runs (same seed -> identical finger multiplicities, asserted):
#   <fp.transfer_tgt>/analysis/bootstrap_replicates.rds           replication (trained within DB3; default Outputs_db3_v2)
#   Outputs_db3_transfer/analysis/bootstrap_replicates.rds        transfer, frozen DB1 thresholds (primary)
#   Outputs_db3_transfer_recal/analysis/bootstrap_replicates.rds  transfer, target-calibrated thresholds (secondary)
#
# Run from the PROJECT ROOT after the three bootstraps:
#   source("Code/transfer_vs_replication.R")
#
# Labelling (STEP_1.1_DB3_SPEC.md §5, QA amendment 13 Oct 2026): a cell is "identical by
# construction" ONLY where the identity is structural and was asserted on the score vectors by
# transfer_db1_to_db3.R (training-free methods: AUC/EER in both variants; everything in the
# target-calibrated variant). Any other cell with all-zero replicate differences is reported
# as DEGENERATE, not as an identity. Intervals are pointwise 95 % percentile, conditional on
# the fitted models, thresholds and fold split.
# Writes only Outputs_db3_transfer/analysis/ci_transfer_minus_replication.csv.
############################################################

REP_DIR <- getOption("fp.transfer_tgt", "Outputs_db3_v2"); FRZ_DIR <- "Outputs_db3_transfer"; RCL_DIR <- "Outputs_db3_transfer_recal"   # the replication run (DB3), as in transfer_db1_to_db3.R
LEVEL <- 0.95; alpha <- (1 - LEVEL) / 2
section <- function(s) cat("\n", s, "\n", strrep("-", nchar(s)), "\n", sep = "")
load_boot <- function(d) { f <- file.path(d, "analysis", "bootstrap_replicates.rds"); if (!file.exists(f)) stop("missing ", f, " - run the bootstrap on ", d, " first."); readRDS(f) }
B_rep <- load_boot(REP_DIR); B_frz <- load_boot(FRZ_DIR); B_rcl <- load_boot(RCL_DIR)
report <- readRDS(file.path(FRZ_DIR, "transfer_report.rds"))
training_free <- report$training_free
stopifnot(all(report$identical_scores[training_free]))            # asserted by transfer_db1_to_db3.R on every fold
source("Code/transfer_scope.R")
# the transfer bootstraps carry the 16 methods of the transfer scope; the replication bootstrap may also carry the clamped
# control (a 17th method in the _v2 runs), which is outside the experiment and is dropped here, nothing else may differ
stopifnot(setequal(B_frz$methods, report$methods), setequal(B_rcl$methods, report$methods), length(report$methods) == 16L)
extra_rep <- setdiff(B_rep$methods, report$methods)
if (!all(extra_rep %in% TRANSFER_SCOPE_OUTSIDE))
  stop("replication bootstrap carries method(s) that are neither in the transfer scope nor the declared clamped control: ", paste(setdiff(extra_rep, TRANSFER_SCOPE_OUTSIDE), collapse = ", "), call. = FALSE)
if (!is.null(report$outside_scope) && !setequal(extra_rep, report$outside_scope))
  stop("methods outside the scope differ between the transfer report (", paste(report$outside_scope, collapse = ", "), ") and the replication bootstrap (", paste(extra_rep, collapse = ", "), ")", call. = FALSE)
stopifnot(all(report$methods %in% B_rep$methods))
cat("  transfer scope: 16 methods;", if (length(extra_rep)) paste("replication also carries", paste(extra_rep, collapse = ", "), "(outside the experiment, dropped)") else "no method outside the scope", "\n")

# ------------------------------------------------------------------ same resamples, same fingers
section("Shared resamples")
fingers_of <- function(d) { m <- readRDS(file.path(d, "run_manifest.rds")); lapply(m$folds, function(mf) sort(unique(m$ids[mf$test_idx]))) }
for (Bt in list(B_frz, B_rcl)) {
  stopifnot(identical(Bt$seed, B_rep$seed), identical(Bt$B, B_rep$B), identical(Bt$level, B_rep$level),
            identical(Bt$fars, B_rep$fars), all(Bt$methods %in% B_rep$methods),
            identical(Bt$multiplicities, B_rep$multiplicities))      # same finger multiplicities, replicate by replicate
}
stopifnot(identical(fingers_of(FRZ_DIR), fingers_of(REP_DIR)), identical(fingers_of(RCL_DIR), fingers_of(REP_DIR)))   # same finger ordering per fold
methods <- report$methods; fars <- B_rep$fars; n_folds <- ncol(B_rep$AUC[[1]])
cat(sprintf("  B = %d, seed = %d; multiplicities identical across the three bootstraps; finger ordering identical in all %d folds.\n", B_rep$B, B_rep$seed, n_folds))

# ------------------------------------------------------------------ differences
ci <- function(x) stats::quantile(x, c(alpha, 1 - alpha), names = FALSE, na.rm = TRUE)
reps_of <- function(Bx, m, metric) {                               # replicate values of the 5-fold mean
  if (metric == "AUC") return(rowMeans(Bx$AUC[[m]]))
  if (metric == "EER") return(rowMeans(Bx$EER[[m]]))
  k <- match(sub("^(TAR|realisedFAR)@", "", metric), as.character(fars))
  if (startsWith(metric, "TAR@")) rowMeans(Bx$TAR[[m]][, , k, drop = TRUE]) else rowMeans(Bx$realisedFAR[[m]][, , k, drop = TRUE])
}
point_of <- function(Bx, m, metric) {
  pf <- Bx$point[[m]]
  if (metric == "AUC") return(mean(vapply(pf, `[[`, numeric(1), "auc")))
  if (metric == "EER") return(mean(vapply(pf, `[[`, numeric(1), "eer")))
  k <- match(sub("^(TAR|realisedFAR)@", "", metric), as.character(fars))
  if (startsWith(metric, "TAR@")) mean(vapply(pf, function(p) p$tar[k], numeric(1))) else mean(vapply(pf, function(p) p$far[k], numeric(1)))
}
metrics <- c("AUC", "EER", as.vector(rbind(paste0("TAR@", fars), paste0("realisedFAR@", fars))))
rows <- list()
for (cal in c("frozen", "recal")) {
  Bt <- if (cal == "frozen") B_frz else B_rcl
  for (m in methods) for (metric in metrics) {
    d  <- reps_of(Bt, m, metric) - reps_of(B_rep, m, metric)
    q  <- ci(d)
    pr <- point_of(B_rep, m, metric); pt <- point_of(Bt, m, metric)
    structural <- (m %in% training_free) && (metric %in% c("AUC", "EER") || cal == "recal")
    if (structural && !(all(d == 0) && pt == pr))
      stop("structural identity violated: ", cal, " ", m, " ", metric, " - transfer and replication differ although the scores were asserted identical.", call. = FALSE)
    label <- if (structural) "identical by construction" else if (all(d == 0)) "degenerate: all replicate differences zero (not an identity)" else ""
    rows[[length(rows) + 1]] <- tibble::tibble(
      calibration = cal, method = m, metric = metric,
      point_replication = pr, point_transfer = pt, difference = pt - pr,
      ci_low = q[1], ci_high = q[2], excludes_zero = (q[1] > 0 | q[2] < 0) & !structural,
      label = label, interval_type = "pointwise 95% percentile, transfer minus replication on identical resamples")
  }
}
diff_df <- dplyr::bind_rows(rows)
out_file <- file.path(FRZ_DIR, "analysis", "ci_transfer_minus_replication.csv")
readr::write_csv(diff_df, out_file)

# ------------------------------------------------------------------ headline (pre-specified, STEP_1.1_DB3_SPEC.md section 5)
fmt <- function(df) as.data.frame(dplyr::mutate(df, dplyr::across(c(point_replication, point_transfer, difference, ci_low, ci_high), ~ round(.x, 4))))
key_metrics <- c("AUC", "EER", "TAR@0.001", "realisedFAR@0.001", "TAR@0.01", "realisedFAR@0.01")
section("1. Transfer (frozen DB1 thresholds) minus replication - Fusion_eq, V23_full")
print(fmt(diff_df[diff_df$calibration == "frozen" & diff_df$method %in% c("Fusion_eq", "V23_full") & diff_df$metric %in% key_metrics,
                  c("method", "metric", "point_replication", "point_transfer", "difference", "ci_low", "ci_high", "excludes_zero")]), row.names = FALSE)
cat("  (AUC/EER: representation + model transfer; TAR is read with the realised FAR on the same row.)\n")
section("2. Transfer (target-calibrated thresholds) minus replication - Fusion_eq, V23_full, TAR only")
print(fmt(diff_df[diff_df$calibration == "recal" & diff_df$method %in% c("Fusion_eq", "V23_full") & diff_df$metric %in% c("TAR@0.001", "realisedFAR@0.001", "TAR@0.01", "realisedFAR@0.01"),
                  c("method", "metric", "point_replication", "point_transfer", "difference", "ci_low", "ci_high", "excludes_zero")]), row.names = FALSE)
section("3. Within transfer (frozen): does the fusion advantage survive transfer? (from bootstrap_ci.R on the transfer directory)")
pd <- readr::read_csv(file.path(FRZ_DIR, "analysis", "ci_paired_differences.csv"), show_col_types = FALSE)
w <- pd[pd$method_a == "Fusion_eq" & pd$method_b %in% c("V23_full", "OT") & pd$metric %in% key_metrics, c("method_a", "method_b", "metric", "difference", "ci_low", "ci_high", "excludes_zero")]
print(as.data.frame(dplyr::mutate(w, dplyr::across(c(difference, ci_low, ci_high), ~ round(.x, 4)))), row.names = FALSE)
section("Other trained methods (frozen), AUC and TAR@0.001 / 0.01")
print(fmt(diff_df[diff_df$calibration == "frozen" & diff_df$method %in% c("Fusion", "V23_noTopK", "V23_globalOnly", "V23_localOnly") & diff_df$metric %in% c("AUC", "TAR@0.001", "TAR@0.01"),
                  c("method", "metric", "point_replication", "point_transfer", "difference", "ci_low", "ci_high", "excludes_zero")]), row.names = FALSE)
section("Labels")
print(table(diff_df$calibration, ifelse(diff_df$label == "", "interval reported", diff_df$label)))
cat("\nWritten:", out_file, "\n")
