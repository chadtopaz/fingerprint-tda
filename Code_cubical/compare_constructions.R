############################################################
# STEP 1.5 (D2/D3) — TRIANGULATED (reported; Outputs_regen) vs CUBICAL (Outputs_cubical)
# on the canonical DB1 inputs: identical images, folds, pairs, seeds and bootstrap resamples.
#
# Run from the PROJECT ROOT after the cubical run AND its bootstrap
#   (options(fp.bootstrap_dir = "Outputs_cubical"); source("Code/bootstrap_ci.R")):
#   source("Code_cubical/compare_constructions.R")
#
# A. Pipeline level: per (fold, method) the two runs' scores; the persistence-independent
#    methods (DT_L2, DT_L2_150, DT_L2_150_aligned) must be bit-identical (structural identity,
#    asserted); for every method, cubical - triangulated paired differences of the 5-fold mean
#    on IDENTICAL finger resamples (multiplicities and finger ordering asserted), pointwise 95 %.
# B. Diagram level (D3): for all 800 images, the global ridge and valley H1 diagrams under both
#    constructions: counts, bottleneck distance delta, and the number / persistence share of
#    features with persistence <= 2*delta in each diagram. In an optimal bottleneck matching a
#    point matched to the diagonal has persistence <= 2*delta, so these counts and shares are
#    UPPER BOUNDS on the diagonal-matched contributions of each diagram - not measurements of
#    "construction-specific" contributions: matched points, including long-lived ones, may
#    also move by up to delta between the two diagrams (QA, 13 Oct). Shares are reported as the
#    median per-image share. Filtration range 0..maxscale px.
# Writes only Outputs_cubical/analysis/{ci_cubical_minus_triangulated.csv, diagram_comparison.csv,
# diagram_comparison_summary.txt, comparison_report.rds}.
############################################################

options(fp.suppress_main = TRUE)
source("Code/fingerprintTDA.R")
source("Code_cubical/cubical_ph.R"); PY <- configure_cubical_python()

# A = triangulated run, B = its cubical counterpart. DB1 by default; DB3 follow-up:
#   options(fp.compare_a = "Outputs_db3", fp.compare_b = "Outputs_cubical_db3"); source("Code_cubical/compare_constructions.R")
A_DIR <- getOption("fp.compare_a", "Outputs_regen_v2"); B_DIR <- getOption("fp.compare_b", "Outputs_cubical"); OUT <- file.path(B_DIR, "analysis")
LEVEL <- 0.95; alpha <- (1 - LEVEL) / 2
section <- function(s) cat("\n", s, "\n", strrep("-", nchar(s)), "\n", sep = "")
RESULTS <- list()
report <- function(label, status, note = "") { RESULTS[[length(RESULTS) + 1]] <<- list(label = label, status = status, note = note)
  cat(sprintf("  [%-9s] %-70s %s\n", status, label, note)); invisible(status) }
chk <- function(cond, label, note = "") report(label, if (isTRUE(cond)) "EXACT" else "FAIL", note)
stopifnot(dir.exists(OUT), file.exists(file.path(OUT, "bootstrap_replicates.rds")), file.exists(file.path(A_DIR, "analysis", "bootstrap_replicates.rds")))

# ------------------------------------------------------------------ A. manifests and caches
section("A1. Same run, different complex")
mA <- readRDS(file.path(A_DIR, "run_manifest.rds")); mB <- readRDS(file.path(B_DIR, "run_manifest.rds"))
cA <- mA$config; cB <- mB$config
chk(identical(cA$polarity %||% "auto", cB$polarity %||% "auto"), "same effective binary polarity rule (absent = auto)", cA$polarity %||% "auto")
for (f in c("output_dir", "export_to_manuscript", "complex_type", "n_workers", "polarity")) { cA[[f]] <- NULL; cB[[f]] <- NULL }
chk(identical(cA, cB), "configurations identical apart from output_dir / export flag / complex_type / n_workers")
chk(identical(mB$config$complex_type, "cubical") && is.null(mA$config$complex_type), "A = triangulated (default), B = cubical")
chk(identical(mA$files, mB$files) && identical(mA$md5, mB$md5) && identical(mA$ids, mB$ids) && identical(mA$folds, mB$folds), "same files, md5s, identities, folds, pair orderings")
read_index <- function(d) readr::read_csv(file.path(d, "fold_cache_index.csv"), show_col_types = FALSE,
                                          col_types = readr::cols(.default = readr::col_character(), fold = readr::col_integer(), loaded_from_cache = readr::col_logical()))
iA <- read_index(A_DIR); iB <- read_index(B_DIR)
chk(setequal(paste(iA$fold, iA$method), paste(iB$fold, iB$method)) && !any(iB$loaded_from_cache), "same (fold, method) set; cubical run computed everything")
chk(unique(iA$fold_key) != unique(iB$fold_key), "different fold keys (complex in the key)", paste(unique(iA$fold_key), unique(iB$fold_key)))
load_cache <- function(d, idx, f, m) { hit <- idx[idx$fold == f & idx$method == m, ]; stopifnot(nrow(hit) == 1); readRDS(file.path(d, "cache", "folds", hit$file)) }
methods <- unique(iA$method); n_folds <- mA$config$n_folds
structural <- c("DT_L2", "DT_L2_150", "DT_L2_150_aligned")           # persistence-independent -> identical by construction
same_vec <- function(a, b) length(a) == length(b) && all(a == b)
section("A2. Persistence-independent methods bit-identical; score correlation for the others")
corr_rows <- list()
for (f in seq_len(n_folds)) for (m in methods) {
  a <- load_cache(A_DIR, iA, f, m); b <- load_cache(B_DIR, iB, f, m)
  stopifnot(identical(a$y_test, b$y_test), identical(a$y_train, b$y_train))
  if (m %in% structural) chk(same_vec(a$s_test, b$s_test) && same_vec(a$s_train, b$s_train) && same_vec(unname(a$thresholds), unname(b$thresholds)),
                             sprintf("fold %d %-18s scores and thresholds identical (structural)", f, m))
  corr_rows[[length(corr_rows) + 1]] <- tibble(fold = f, method = m, test_score_corr = cor(a$s_test, b$s_test), train_score_corr = cor(a$s_train, b$s_train))
}
corr_df <- bind_rows(corr_rows)
if (any(vapply(RESULTS, function(r) r$status == "FAIL", logical(1)))) stop("structural identity violated - see FAIL lines above.", call. = FALSE)
cs <- corr_df %>% group_by(method) %>% summarize(test_score_corr = mean(test_score_corr), .groups = "drop") %>% arrange(match(method, methods))
cat("  mean test-score correlation cubical vs triangulated, per method:\n"); print(as.data.frame(mutate(cs, test_score_corr = round(test_score_corr, 4))), row.names = FALSE)

# ------------------------------------------------------------------ A3. paired differences on identical resamples
section("A3. Cubical - triangulated paired differences (5-fold mean, identical resamples, pointwise 95 %)")
BA <- readRDS(file.path(A_DIR, "analysis", "bootstrap_replicates.rds")); BB <- readRDS(file.path(OUT, "bootstrap_replicates.rds"))
fingers_of <- function(m) lapply(m$folds, function(mf) sort(unique(m$ids[mf$test_idx])))
stopifnot(identical(BA$seed, BB$seed), identical(BA$B, BB$B), identical(BA$fars, BB$fars), setequal(BA$methods, BB$methods),
          identical(BA$multiplicities, BB$multiplicities), identical(fingers_of(mA), fingers_of(mB)))
cat("  multiplicities identical across the two bootstraps; finger ordering identical in all folds.\n")
fars <- BA$fars
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
for (m in BA$methods) for (metric in metrics) {
  d <- reps_of(BB, m, metric) - reps_of(BA, m, metric); q <- ci(d)
  pa <- point_of(BA, m, metric); pb <- point_of(BB, m, metric)
  is_struct <- m %in% structural
  if (is_struct && !(all(d == 0) && pa == pb)) stop("structural identity violated in the bootstrap: ", m, " ", metric, call. = FALSE)
  rows[[length(rows) + 1]] <- tibble(method = m, metric = metric, point_triangulated = pa, point_cubical = pb, difference = pb - pa,
                                     ci_low = q[1], ci_high = q[2], excludes_zero = (q[1] > 0 | q[2] < 0) & !is_struct,
                                     label = if (is_struct) "identical by construction" else if (all(d == 0)) "degenerate: all replicate differences zero (not an identity)" else "",
                                     interval_type = "pointwise 95% percentile, cubical minus triangulated on identical resamples")
}
diff_df <- bind_rows(rows); write_csv(diff_df, file.path(OUT, "ci_cubical_minus_triangulated.csv"))
show <- c("GlobalBetti_L1", "GlobalPI_L2", "V23_full", "OT", "Fusion", "Fusion_eq")
hd <- diff_df %>% filter(method %in% show, metric %in% c("AUC", "EER", "TAR@0.001", "realisedFAR@0.001", "TAR@0.01", "realisedFAR@0.01")) %>%
  mutate(method = factor(method, levels = show)) %>% arrange(method) %>% mutate(across(where(is.numeric), ~ round(.x, 4))) %>% select(-label, -interval_type)
print(as.data.frame(hd), row.names = FALSE)
cat("\n  Within-cubical pre-specified comparisons (from bootstrap_ci.R on Outputs_cubical): do the conclusions hold under the cubical complex?\n")
pdB <- read_csv(file.path(OUT, "ci_paired_differences.csv"), show_col_types = FALSE); pdA <- read_csv(file.path(A_DIR, "analysis", "ci_paired_differences.csv"), show_col_types = FALSE)
keyp <- pdB %>% filter(method_a %in% c("Fusion_eq", "GlobalBetti_L1", "V23_full") & method_b %in% c("V23_full", "OT", "DT_L2_150", "DT_L2_150_aligned"),
                       metric %in% c("AUC", "TAR@0.001", "TAR@0.01")) %>% select(method_a, method_b, metric, difference, ci_low, ci_high, excludes_zero) %>%
  inner_join(pdA %>% select(method_a, method_b, metric, difference_tri = difference, excl_tri = excludes_zero), by = c("method_a", "method_b", "metric"))
print(as.data.frame(mutate(keyp, across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)

# ------------------------------------------------------------------ B. diagram level (D3)
section("B. Global H1 diagrams of all 800 images under both constructions (ridge and valley filtrations)")
cfg <- mA$config; paths <- mA$paths; stopifnot(all(file.exists(paths)))
future::plan(future::multisession, workers = CONFIG$n_workers)
t0 <- Sys.time()
diag_rows <- progressr::with_progress({
  p <- progressr::progressor(steps = length(paths))
  future.apply::future_lapply(seq_along(paths), function(i) {
    bin <- read_binary_fingerprint(paths[i])
    out <- list()
    for (ch in c("ridge", "valley")) {
      fld <- downsample_matrix(if (ch == "ridge") compute_dt_to_ridge(bin) else compute_dt_to_valley(bin), cfg$target_size)
      dT <- as.matrix(compute_persistence(fld, 1, cfg$maxscale, "triangulated")); dC <- compute_persistence(fld, 1, cfg$maxscale, "cubical")
      hT <- dT[dT[, 1] == 1, , drop = FALSE]; hC <- dC[dC[, 1] == 1, , drop = FALSE]
      pT <- hT[, 3] - hT[, 2]; pC <- hC[, 3] - hC[, 2]
      delta <- TDA::bottleneck(dT, dC, dimension = 1)
      out[[ch]] <- tibble(image = basename(paths[i]), channel = ch,
                          n_h1_tri = nrow(hT), n_h1_cub = nrow(hC), bottleneck = delta,
                          max_pers_tri = if (length(pT)) max(pT) else 0, max_pers_cub = if (length(pC)) max(pC) else 0,
                          n_le_2delta_tri = sum(pT <= 2 * delta + 1e-9), n_le_2delta_cub = sum(pC <= 2 * delta + 1e-9),
                          pers_share_le_2delta_tri = if (sum(pT) > 0) sum(pT[pT <= 2 * delta + 1e-9]) / sum(pT) else NA_real_,
                          pers_share_le_2delta_cub = if (sum(pC) > 0) sum(pC[pC <= 2 * delta + 1e-9]) / sum(pC) else NA_real_,
                          total_pers_tri = sum(pT), total_pers_cub = sum(pC))
    }
    p(); bind_rows(out)
  }, future.seed = TRUE)
})
future::plan(future::sequential)
dd <- bind_rows(diag_rows); write_csv(dd, file.path(OUT, "diagram_comparison.csv"))
cat(sprintf("  %d images x 2 channels in %.1f min\n", length(paths), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
summ <- dd %>% group_by(channel) %>% summarize(
  n = n(), mean_h1_tri = mean(n_h1_tri), mean_h1_cub = mean(n_h1_cub), ratio_cub_tri = mean(n_h1_cub / n_h1_tri),
  identical_counts = sum(n_h1_tri == n_h1_cub),
  bottleneck_median = median(bottleneck), bottleneck_q90 = quantile(bottleneck, 0.9), bottleneck_max = max(bottleneck),
  two_delta_max = 2 * max(bottleneck), max_pers_tri_median = median(max_pers_tri),
  excess_cub_mean = mean(n_h1_cub - n_h1_tri), n_le_2delta_cub_mean = mean(n_le_2delta_cub), n_le_2delta_tri_mean = mean(n_le_2delta_tri),
  pers_share_le_2delta_cub_median = median(pers_share_le_2delta_cub, na.rm = TRUE), pers_share_le_2delta_tri_median = median(pers_share_le_2delta_tri, na.rm = TRUE),
  .groups = "drop")
print(as.data.frame(mutate(summ, across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
lines <- c(
  sprintf("Global H1 diagrams, %d images, filtration range 0..%g px (distance transform, %dx%d grid).", length(paths), cfg$maxscale, cfg$target_size, cfg$target_size),
  "Bottleneck distance delta between the triangulated and cubical H1 diagram of the SAME image and channel.",
  "A point matched to the diagonal in an optimal bottleneck matching has persistence <= 2*delta, so n_le_2delta_* and",
  "pers_share_le_2delta_* are UPPER BOUNDS on the diagonal-matched contributions of each diagram; they do not measure",
  "construction-specific contributions (matched points, long-lived ones included, may move by up to delta).",
  "pers_share_le_2delta_*_median = median over images of the per-image share of total persistence in features with persistence <= 2*delta.",
  capture.output(print(as.data.frame(mutate(summ, across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)))
writeLines(lines, file.path(OUT, "diagram_comparison_summary.txt"))
saveRDS(list(created = Sys.time(), A = A_DIR, B = B_DIR, results = RESULTS, score_corr = corr_df, differences = diff_df, diagram_summary = summ),
        file.path(OUT, "comparison_report.rds"))
st <- vapply(RESULTS, `[[`, character(1), "status"); print(table(factor(st, levels = c("EXACT", "FAIL"))))
cat("\nOutputs written to", OUT, "\n")
