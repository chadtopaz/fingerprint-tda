############################################################
# STEP 1.4 (C6) — bootstrap every sensitivity condition with the unchanged Code/bootstrap_ci.R
# (same seed -> the same finger resamples as the clean run) and compute PERTURBED - CLEAN paired
# differences of the 5-fold mean on identical pairs and resamples, pointwise 95 %.
#
# Run from the PROJECT ROOT after Code/sensitivity/score_perturbed.R (~1.3 min per condition):
#   source("Code/sensitivity/run_sensitivity_bootstraps.R")
# Reads <fp.sens_src>/analysis/bootstrap_replicates.rds (clean; default Outputs_regen_v2) and Outputs_sensitivity/<cond>/;
# writes Outputs_sensitivity/<cond>/analysis/ (via bootstrap_ci.R) and
# Outputs_sensitivity/summary/ci_perturbed_minus_clean.csv.
############################################################

OUT_ROOT <- getOption("fp.sens_out", "Outputs_sensitivity"); CLEAN_DIR <- getOption("fp.sens_src", "Outputs_regen_v2")
LEVEL <- 0.95; alpha <- (1 - LEVEL) / 2
section <- function(s) cat("\n", s, "\n", strrep("-", nchar(s)), "\n", sep = "")
REQUESTED <- getOption("fp.sens_requested", c("filter_m2", "filter_m4", "scale_055", "scale_075", "segment_020", "segment_040"))   # Tier B; Tier R via run_realistic.R
# every requested condition must be COMPLETE (marker + index + caches), QA 13 Oct: an interrupted scoring run cannot silently drop a condition
complete <- vapply(REQUESTED, function(cond) {
  d <- file.path(OUT_ROOT, cond)
  if (!file.exists(file.path(d, "SENSITIVITY_COMPLETE.txt")) || !file.exists(file.path(d, "fold_cache_index.csv"))) return(FALSE)
  idx <- readr::read_csv(file.path(d, "fold_cache_index.csv"), show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  nrow(idx) == 75 && all(file.exists(file.path(d, "cache", "folds", idx$file))) && file.exists(file.path(d, "run_manifest.rds"))
}, logical(1))
if (!all(complete)) stop("requested conditions not complete: ", paste(REQUESTED[!complete], collapse = ", "), " - run Code/sensitivity/score_perturbed.R first.", call. = FALSE)
conds <- REQUESTED
stopifnot(file.exists(file.path(CLEAN_DIR, "analysis", "bootstrap_replicates.rds")))
# all requested conditions must share one pairing protocol (pairing + probe-role assignment), QA 14 Oct
`%||%` <- function(x, y) if (!is.null(x)) x else y
proto <- lapply(conds, function(cond) { sv <- readRDS(file.path(OUT_ROOT, cond, "run_manifest.rds"))$sensitivity
  list(pairing = sv$pairing %||% "both", sha = sv$probe_assignment_sha %||% NA_character_) })
if (length(unique(vapply(proto, function(p) paste(p$pairing, p$sha), character(1)))) != 1) stop("the requested conditions were scored under different pairing protocols: ",
  paste(sprintf("%s=%s/%s", conds, vapply(proto, `[[`, character(1), "pairing"), substr(vapply(proto, `[[`, character(1), "sha"), 1, 8)), collapse = ", "), call. = FALSE)
cat(sprintf("  pairing protocol of all %d conditions: '%s'%s\n", length(conds), proto[[1]]$pairing, if (!is.na(proto[[1]]$sha)) sprintf(" (probe-role assignment sha256 %s...)", substr(proto[[1]]$sha, 1, 12)) else ""))

# ------------------------------------------------------------------ 1. bootstraps (skip only those whose signature matches)
# Reuse signature (audit patch, 31 Oct 2026): an existing bootstrap is reused only if it was computed from
# the SAME fold caches (md5s of every file in the condition's index) with the same settings as the clean
# replicates it is compared with; previously the mere existence of bootstrap_replicates.rds was enough.
boot_signature <- function(d) {
  idx <- readr::read_csv(file.path(d, "fold_cache_index.csv"), show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  files <- sort(idx$file)
  list(version = 1L, cache_md5 = setNames(unname(tools::md5sum(file.path(d, "cache", "folds", files))), files),
       clean_dir = CLEAN_DIR, clean_replicates_md5 = unname(tools::md5sum(file.path(CLEAN_DIR, "analysis", "bootstrap_replicates.rds"))),
       bootstrap_script_sha256 = digest::digest(file = "Code/bootstrap_ci.R", algo = "sha256"))
}
for (cond in conds) {
  d <- file.path(OUT_ROOT, cond); sig_file <- file.path(d, "analysis", "bootstrap_signature.rds")
  if (file.exists(file.path(d, "analysis", "bootstrap_replicates.rds"))) {
    if (file.exists(sig_file) && identical(readRDS(sig_file), boot_signature(d))) { cat("bootstrap present with matching signature, skipped:", cond, "\n"); next }
    if (!file.exists(sig_file) && isTRUE(getOption("fp.sens_allow_unsigned_reuse", FALSE))) { cat("bootstrap present WITHOUT signature, reused (fp.sens_allow_unsigned_reuse):", cond, "\n"); next }
    stop("bootstrap of ", cond, " exists but ", if (file.exists(sig_file)) "its signature does not match the current caches/settings" else "has no reuse signature",
         ". Remove <cond>/analysis to recompute it, or set options(fp.sens_allow_unsigned_reuse = TRUE) for pre-patch directories.", call. = FALSE)
  }
  section(paste("Bootstrap:", cond))
  options(fp.bootstrap_dir = d)
  source("Code/bootstrap_ci.R", local = new.env())      # the script's own globals stay in a fresh environment
  saveRDS(boot_signature(d), sig_file)
}
options(fp.bootstrap_dir = NULL)

# ------------------------------------------------------------------ 2. perturbed - clean on identical resamples
section("Perturbed - clean paired differences (5-fold mean, identical resamples, pointwise 95 %)")
B0 <- readRDS(file.path(CLEAN_DIR, "analysis", "bootstrap_replicates.rds"))
fingers_of <- function(d) { m <- readRDS(file.path(d, "run_manifest.rds")); lapply(m$folds, function(mf) sort(unique(m$ids[mf$test_idx]))) }
f0 <- fingers_of(CLEAN_DIR); fars <- B0$fars
ci <- function(x) stats::quantile(x, c(alpha, 1 - alpha), names = FALSE, na.rm = TRUE)
reps_of <- function(Bx, m, metric) {
  if (metric == "AUC") return(rowMeans(Bx$AUC[[m]])); if (metric == "EER") return(rowMeans(Bx$EER[[m]]))
  k <- match(sub("^(TAR|realisedFAR)@", "", metric), as.character(fars))
  if (startsWith(metric, "TAR@")) rowMeans(Bx$TAR[[m]][, , k, drop = TRUE]) else rowMeans(Bx$realisedFAR[[m]][, , k, drop = TRUE])
}
point_of <- function(Bx, m, metric) {
  pf <- Bx$point[[m]]
  if (metric == "AUC") return(mean(vapply(pf, `[[`, numeric(1), "auc"))); if (metric == "EER") return(mean(vapply(pf, `[[`, numeric(1), "eer")))
  k <- match(sub("^(TAR|realisedFAR)@", "", metric), as.character(fars))
  if (startsWith(metric, "TAR@")) mean(vapply(pf, function(p) p$tar[k], numeric(1))) else mean(vapply(pf, function(p) p$far[k], numeric(1)))
}
metrics <- c("AUC", "EER", as.vector(rbind(paste0("TAR@", fars), paste0("realisedFAR@", fars))))
rows <- list()
for (cond in conds) {
  d <- file.path(OUT_ROOT, cond); Bc <- readRDS(file.path(d, "analysis", "bootstrap_replicates.rds"))
  stopifnot(identical(Bc$seed, B0$seed), identical(Bc$B, B0$B), identical(Bc$fars, B0$fars), identical(Bc$multiplicities, B0$multiplicities),
            identical(fingers_of(d), f0), all(Bc$methods %in% B0$methods))
  for (m in Bc$methods) for (metric in metrics) {
    dd <- reps_of(Bc, m, metric) - reps_of(B0, m, metric); q <- ci(dd); p0 <- point_of(B0, m, metric); p1 <- point_of(Bc, m, metric)
    rows[[length(rows) + 1]] <- tibble::tibble(condition = cond, method = m, metric = metric, point_clean = p0, point_perturbed = p1, difference = p1 - p0,
                                               ci_low = q[1], ci_high = q[2], excludes_zero = (q[1] > 0 | q[2] < 0),
                                               label = if (all(dd == 0)) "degenerate: all replicate differences zero (not an identity)" else "",
                                               interval_type = "pointwise 95% percentile, perturbed minus clean on identical resamples")
  }
}
diff_df <- dplyr::bind_rows(rows)
pe <- file.path(OUT_ROOT, "summary", "point_estimates_by_condition.csv")
if (file.exists(pe)) {
  ag <- readr::read_csv(pe, show_col_types = FALSE) |> dplyr::distinct(condition, tier, agreement_median, agreement_q10, agreement_min)
  diff_df <- dplyr::left_join(diff_df, ag, by = "condition")
}
dir.create(file.path(OUT_ROOT, "summary"), showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(OUT_ROOT, "summary", "ci_perturbed_minus_clean.csv"); readr::write_csv(diff_df, out_file)
cat("  resamples identical to the clean run in every condition; finger ordering identical.\n")

# ------------------------------------------------------------------ 3. headline
order_c <- c("erode1", "dilate1", "breaks50", "breaks150", "bridges50", "bridges150", "speckle0.005", "speckle0.02",
             "filter_m2", "filter_m4", "scale_055", "scale_075", "segment_020", "segment_040",
             "dry_typ", "dry_up", "wet_typ", "wet_up", "noise_typ", "noise_up", "blur_typ", "blur_up")   # Tier R (spec v3)
fmt <- function(df) as.data.frame(dplyr::mutate(df, dplyr::across(c(point_clean, point_perturbed, difference, ci_low, ci_high), ~ round(.x, 4))))
for (m in c("Fusion_eq", "V23_full", "OT", "GlobalBetti_L1")) {
  section(paste("Perturbed - clean:", m, "(AUC, EER, TAR@1e-3 with realised FAR, TAR@1e-2 with realised FAR)"))
  h <- diff_df[diff_df$method == m & diff_df$metric %in% c("AUC", "EER", "TAR@0.001", "realisedFAR@0.001", "TAR@0.01", "realisedFAR@0.01"), ]
  h$condition <- factor(h$condition, levels = c(intersect(order_c, unique(h$condition)), setdiff(unique(h$condition), order_c)))
  h <- h[order(h$condition, match(h$metric, metrics)), c("condition", "metric", "point_clean", "point_perturbed", "difference", "ci_low", "ci_high", "excludes_zero")]
  print(fmt(h), row.names = FALSE)
}
cat("\nWritten:", out_file, "\n")
