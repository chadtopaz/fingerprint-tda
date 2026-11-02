############################################################
# STEP 2.1 + 2.2 — Finger-level bootstrap intervals; fusion weights;
#                  equal-weight fusion
#
# Run from the PROJECT ROOT after step 0.7:
#   source("Code/bootstrap_ci.R")
#
# Reads only Outputs_repro/ (caches + manifest + archived fold metrics).
# Writes only Outputs_repro/analysis/.
#
# Resampling unit: the FINGER. For each replicate and fold, 20 fingers are
# drawn with replacement from the fold's 20 and tabulated into multiplicities
# m_i. Every OBSERVED pair keeps its score and gets a weight:
#     genuine pair of finger i         -> m_i
#     impostor pair, fingers i != j    -> m_i * m_j
# Fingers with m_i = 0 drop out. Metrics are weighted means / weighted
# U-statistics normalised by total weight. No pair is invented: no image is
# compared with itself and copies of one finger never form an impostor pair.
# Thresholds are held FIXED at the training-calibrated values.
# The same multiplicities are applied to every method, so paired differences
# are computed on identical resamples.
#
# Intervals are CONDITIONAL on the fitted models, thresholds and fold split:
# they quantify test-finger sampling variability of the reported number and
# do not capture refitting or re-splitting uncertainty. They are POINTWISE.
#
# Revision step 1.2 integration (2026-10-11): methods are whatever the pipeline
# wrote to fold_metrics_auc_eer.csv (now including DT_L2_150_aligned and the
# pipeline method Fusion_eq); per-fold caches are resolved through the pipeline's
# fold_cache_index.csv (method-specific file suffixes), not by rebuilding names.
# Fusion_eq is no longer constructed here; it is checked against its components
# and, when present, against the frozen pre-integration reference in
# Outputs_repro/snapshot_before_new_methods/analysis/fusion_eq_scores.rds.
# This script no longer writes fusion_eq_scores.rds.
############################################################

options(fp.suppress_main = TRUE)
source("Code/fingerprintTDA.R")

# Run directory. Default = the archived-input run; the regenerated-input run (step 0.4) is
# analysed with  options(fp.bootstrap_dir = "Outputs_regen"); source("Code/bootstrap_ci.R")
REPRO_DIR <- getOption("fp.bootstrap_dir", "Outputs_repro")
OUT       <- file.path(REPRO_DIR, "analysis")
B         <- 2000L
SEED      <- 20261010L
LEVEL     <- 0.95
TOL       <- 1e-9
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------ reporting
RESULTS <- list()
report <- function(label, status, note = "") {
  RESULTS[[length(RESULTS) + 1]] <<- list(label = label, status = status, note = note)
  cat(sprintf("  [%-9s] %-60s %s\n", status, label, note)); invisible(status)
}
cmp <- function(a, b, label, tol = TOL) {
  a <- as.numeric(unlist(a)); b <- as.numeric(unlist(b))
  if (length(a) != length(b)) return(report(label, "FAIL", sprintf("length %d vs %d", length(a), length(b))))
  if (any(!is.finite(a)) || any(!is.finite(b))) return(report(label, "FAIL", "non-finite"))
  if (all(a == b)) return(report(label, "EXACT", ""))
  m <- max(abs(a - b))
  if (m <= tol) report(label, "TOLERANCE", sprintf("max|diff| = %.3g", m)) else report(label, "FAIL", sprintf("max|diff| = %.3g", m))
}
section <- function(s) cat("\n", s, "\n", strrep("-", nchar(s)), "\n", sep = "")
stop_if_fail <- function(what) {
  st <- vapply(RESULTS, `[[`, character(1), "status")
  if (any(st == "FAIL")) stop(what, ": a check FAILED (see above). Stopping before any interval is computed.", call. = FALSE)
}

# ============================================================ weighted metrics
# Structures that do not depend on the weights, computed once per (fold, method).
prep_method <- function(s, gen, thr) {
  iG <- which(gen); iI <- which(!gen)
  sG <- s[iG]; sI <- s[iI]
  oG <- order(sG); oI <- order(sI)
  sGs <- sG[oG]; sIs <- sI[oI]
  u <- sort(unique(s))
  list(
    iG = iG, iI = iI, oG = oG, oI = oI,
    # AUC: per genuine score, number of impostor scores strictly below / at most equal
    n_lt = findInterval(sG, sIs, left.open = TRUE),
    n_le = findInterval(sG, sIs),
    # EER on pROC's threshold set: per unique value u_k, counts of scores <= u_k by class
    nI_le_u = findInterval(u, sIs),
    nG_le_u = findInterval(u, sGs),
    K = length(u),
    # fixed thresholds: acceptance (score >= threshold) by class
    accG = lapply(thr, function(t) sG >= t),
    accI = lapply(thr, function(t) sI >= t)
  )
}

# Weighted AUC (Mann-Whitney, ties 1/2), weighted EER mirroring pROC's discrete
# rule (thresholds -Inf, midpoints of consecutive unique scores, +Inf; first
# index minimising |FPR - FNR|; (FPR+FNR)/2), weighted TAR and realised FAR at
# the fixed thresholds.
weighted_metrics <- function(P, w) {
  wG <- w[P$iG]; wI <- w[P$iI]
  WG <- sum(wG); WI <- sum(wI)
  if (WG <= 0 || WI <= 0)   # only if every draw was the same finger; metrics undefined
    return(list(auc = NA_real_, eer = NA_real_, tar = rep(NA_real_, length(P$accG)), far = rep(NA_real_, length(P$accI))))
  cwI <- c(0, cumsum(wI[P$oI]))         # cwI[k+1] = total weight of the k smallest impostors
  cwG <- c(0, cumsum(wG[P$oG]))
  w_lt <- cwI[P$n_lt + 1L]; w_le <- cwI[P$n_le + 1L]
  auc  <- sum(wG * (w_lt + 0.5 * (w_le - w_lt))) / (WG * WI)
  K <- P$K
  # Written in pROC's arithmetic order (fpr = 1 - specificity, fnr = 1 - sensitivity).
  # Algebraically equivalent forms differ in the last bit, which can change the
  # which.min() tie-break between ROC points and hence the EER. Do not "simplify".
  fpr_mid <- 1 - cwI[P$nI_le_u + 1L] / WI               # 1 - spec,  spec = weight(impostors <= u_k) / WI
  fnr_mid <- 1 - (WG - cwG[P$nG_le_u + 1L]) / WG        # 1 - sens,  sens = weight(genuine   >  u_k) / WG
  fpr <- c(1, fpr_mid[-K], 0); fnr <- c(0, fnr_mid[-K], 1)
  idx <- which.min(abs(fpr - fnr)); eer <- (fpr[idx] + fnr[idx]) / 2
  tar <- vapply(P$accG, function(a) sum(wG[a]) / WG, numeric(1))
  far <- vapply(P$accI, function(a) sum(wI[a]) / WI, numeric(1))
  list(auc = auc, eer = eer, tar = tar, far = far)
}

# Reference implementation on explicitly repeated records (for the weighting self-test).
materialised_metrics <- function(s, gen, w, thr) {
  rs <- rep(s, w); rg <- rep(gen, w)
  roc <- pROC::roc(response = as.integer(rg), predictor = as.numeric(rs), quiet = TRUE, direction = "<")
  sens <- roc$sensitivities; fpr <- 1 - roc$specificities; fnr <- 1 - sens
  i <- which.min(abs(fpr - fnr))
  list(auc = as.numeric(pROC::auc(roc)), eer = (fpr[i] + fnr[i]) / 2,
       tar = vapply(thr, function(t) mean(rs[rg] >= t),  numeric(1)),
       far = vapply(thr, function(t) mean(rs[!rg] >= t), numeric(1)))
}

# ============================================================ self-test 1: weighting
section("Self-test 1: weighting scheme vs explicitly repeated records (ties included)")
{
  ids6 <- c("A","A","B","B","C","C"); pm6 <- combn(6, 2)
  f1 <- ids6[pm6[1, ]]; f2 <- ids6[pm6[2, ]]; gen6 <- f1 == f2
  # 15 pairs in combn order; genuine at positions 1 (A-A), 10 (B-B), 15 (C-C).
  # Ties: genuine 0.6 == impostor 0.6; impostor 0.1 x2, 0.5 x2; genuine 0.6 x2.
  s6 <- c(0.9, 0.1, 0.6, 0.3, 0.4, 0.2, 0.5, 0.7, 0.1, 0.6, 0.3, 0.8, 0.2, 0.5, 0.6)
  thr6 <- c(far_a = 0.5, far_b = 0.65, far_c = 0.95)
  # Regression case (QA): with the algebraically-equivalent-but-differently-rounded EER
  # formulas, unit multiplicities gave EER 0.625 against pROC's 0.7083333.
  s6_reg <- c(.7, .2, .2, .8, .5, .3, .7, 0, .3, 0, .5, .8, .1, .1, .1)
  mults <- list(c(A = 2, B = 1, C = 0), c(A = 1, B = 0, C = 3), c(A = 1, B = 1, C = 1), c(A = 3, B = 2, C = 1))
  for (case in list(list(name = "ties", s = s6), list(name = "QA-regression", s = s6_reg))) {
    P6 <- prep_method(case$s, gen6, thr6)
    for (mv in mults) {
      w6 <- unname(mv[f1] * ifelse(gen6, 1, mv[f2]))
      wm <- weighted_metrics(P6, w6); mm <- materialised_metrics(case$s, gen6, w6, thr6)
      lab <- sprintf("%-13s m=(%s)", case$name, paste(mv, collapse = ","))
      cmp(wm$auc, mm$auc, paste(lab, "AUC")); cmp(wm$eer, mm$eer, paste(lab, "EER"))
      cmp(wm$tar, mm$tar, paste(lab, "TAR")); cmp(wm$far, mm$far, paste(lab, "realised FAR"))
    }
  }
}
stop_if_fail("Self-test 1")

# ============================================================ load
section(paste("Loading", REPRO_DIR))
manifest <- readRDS(file.path(REPRO_DIR, "run_manifest.rds"))
CFG   <- manifest$config
fars  <- CFG$fars_report
far_names <- paste0("far_", fars)
fkey  <- substr(fold_cache_key(CFG), 1, 8)

# Per-fold caches are resolved through the pipeline's fold_cache_index.csv, which
# records the file each (fold, method) was read from or written to in the run that
# produced the metrics CSVs (method-specific suffixes, e.g. DT_L2_150_aligned).
index_file  <- file.path(REPRO_DIR, "fold_cache_index.csv")
if (!file.exists(index_file)) stop("fold_cache_index.csv not found in ", REPRO_DIR, ": run the pipeline (Code/run_revision_methods.R) first.")
cache_index <- readr::read_csv(index_file, show_col_types = FALSE, col_types = readr::cols(.default = readr::col_character(), fold = readr::col_integer(), loaded_from_cache = readr::col_logical()))
if (!all(cache_index$fold_key == fkey)) stop("fold_cache_index.csv was written under a different configuration key than the manifest's config implies.")
if (any(duplicated(cache_index[, c("fold", "method")]))) stop("fold_cache_index.csv has duplicate (fold, method) entries.")
load_fold <- function(f, m) {
  hit <- cache_index[cache_index$fold == f & cache_index$method == m, ]
  if (nrow(hit) != 1) stop("fold_cache_index.csv has no entry for fold ", f, " method ", m)
  p <- file.path(REPRO_DIR, "cache", "folds", hit$file)
  if (!file.exists(p)) stop("missing fold cache: ", p)
  cc <- readRDS(p)
  if (!is.na(hit$extra_key) && !identical(cc$aligned_key, hit$extra_key))
    stop("fold cache ", p, " does not carry the extra key recorded in the index")
  cc
}
arch_auc <- readr::read_csv(file.path(REPRO_DIR, "fold_metrics_auc_eer.csv"), show_col_types = FALSE)
arch_tar <- readr::read_csv(file.path(REPRO_DIR, "fold_metrics_tar.csv"),     show_col_types = FALSE)
methods_pipeline <- unique(arch_auc$method)
if (!setequal(methods_pipeline, unique(cache_index$method))) stop("methods in fold_metrics_auc_eer.csv and fold_cache_index.csv differ")
if (!"Fusion_eq" %in% methods_pipeline) stop("Fusion_eq is expected as a PIPELINE method (revision step 2.2b); run Code/run_revision_methods.R first.")
n_folds <- CFG$n_folds
cat("  methods:", length(methods_pipeline), "| folds:", n_folds, "| FARs:", paste(fars, collapse = ", "), "\n")

# Frozen pre-integration Fusion_eq reference (made once by run_revision_methods.R).
FROZEN_FUSION_EQ <- file.path(REPRO_DIR, "snapshot_before_new_methods", "analysis", "fusion_eq_scores.rds")
frozen_eq <- if (file.exists(FROZEN_FUSION_EQ)) readRDS(FROZEN_FUSION_EQ) else NULL

# ------------------------------------------------------------ per-fold data
FOLD <- vector("list", n_folds)
for (f in seq_len(n_folds)) {
  mf  <- manifest$folds[[f]]
  pm  <- mf$pair_test
  fid1 <- manifest$ids[pm[1, ]]; fid2 <- manifest$ids[pm[2, ]]
  fingers <- sort(unique(manifest$ids[mf$test_idx]))
  F1 <- match(fid1, fingers); F2 <- match(fid2, fingers)
  gen <- F1 == F2
  D <- list(F1 = F1, F2 = F2, gen = gen, n_fingers = length(fingers), fingers = fingers, M = list())
  for (m in methods_pipeline) {
    cc <- load_fold(f, m)
    D$M[[m]] <- list(s_test = cc$s_test, s_train = cc$s_train, y_train = cc$y_train, y_test = cc$y_test,
                     thr = cc$thresholds[far_names], best_w = cc$best_w, zfit_ot = cc$zfit_ot, zfit_v = cc$zfit_v)
  }
  FOLD[[f]] <- D
}
methods <- unique(methods_pipeline)   # Fusion_eq is a pipeline method; nothing is appended

# ============================================================ self-test 2: labels, thresholds, point estimates
section("Self-test 2: labels from manifest, stored thresholds, and unit-multiplicity point estimates vs archived")
for (f in seq_len(n_folds)) {
  D <- FOLD[[f]]
  for (m in methods_pipeline) {
    cmp(as.integer(D$gen), D$M[[m]]$y_test, sprintf("fold %d %-17s y_test == genuine mask from manifest", f, m))
    imp <- D$M[[m]]$s_train[D$M[[m]]$y_train == 0]
    cmp(vapply(fars, function(fa) as.numeric(quantile(imp, probs = 1 - fa, na.rm = TRUE)), numeric(1)),
        D$M[[m]]$thr, sprintf("fold %d %-17s stored thresholds", f, m))
  }
  # Fusion_eq (pipeline method) == w = 0.5 fusion of the cached components under the
  # tuned Fusion's stored normalisation (same deterministic fit of the same inputs) ...
  fu <- D$M[["Fusion"]]; ot <- D$M[["OT"]]; v <- D$M[["V23_full"]]; eq <- D$M[["Fusion_eq"]]
  cmp(eq$s_train, apply_fusion(ot$s_train, v$s_train, 0.5, fu$zfit_ot, fu$zfit_v), sprintf("fold %d Fusion_eq s_train == 0.5-fusion of components", f))
  cmp(eq$s_test,  apply_fusion(ot$s_test,  v$s_test,  0.5, fu$zfit_ot, fu$zfit_v), sprintf("fold %d Fusion_eq s_test  == 0.5-fusion of components", f))
  cmp(unlist(eq$zfit_ot), unlist(fu$zfit_ot), sprintf("fold %d Fusion_eq zfit_ot == Fusion zfit_ot", f))
  cmp(unlist(eq$zfit_v),  unlist(fu$zfit_v),  sprintf("fold %d Fusion_eq zfit_v  == Fusion zfit_v", f))
  # ... and == the frozen pre-integration reference, when it is present.
  if (!is.null(frozen_eq)) {
    cmp(eq$s_train, frozen_eq$s_train[[f]], sprintf("fold %d Fusion_eq s_train == frozen reference", f))
    cmp(eq$s_test,  frozen_eq$s_test[[f]],  sprintf("fold %d Fusion_eq s_test  == frozen reference", f))
    cmp(eq$thr,     frozen_eq$thresholds[[f]][far_names], sprintf("fold %d Fusion_eq thresholds == frozen reference", f))
  } else {
    report(sprintf("fold %d Fusion_eq vs frozen reference", f), "SKIPPED", "no snapshot_before_new_methods/analysis/fusion_eq_scores.rds")
  }
}
PREP <- lapply(seq_len(n_folds), function(f) lapply(FOLD[[f]]$M, function(M) prep_method(M$s_test, FOLD[[f]]$gen, M$thr)))
point <- list()   # point[[m]][[f]] = unit-weight metrics
for (f in seq_len(n_folds)) {
  w1 <- rep(1, length(FOLD[[f]]$gen))
  for (m in methods) {
    pt <- weighted_metrics(PREP[[f]][[m]], w1); point[[m]][[f]] <- pt
    if (m %in% methods_pipeline) {
      a <- arch_auc[arch_auc$fold == f & arch_auc$method == m, ]
      cmp(pt$auc, a$AUC, sprintf("fold %d %-15s AUC vs archived", f, m))
      cmp(pt$eer, a$EER, sprintf("fold %d %-15s EER vs archived", f, m))
      tt <- arch_tar[arch_tar$fold == f & arch_tar$method == m, ]
      tt <- tt[match(round(fars, 10), round(tt$FAR, 10)), ]   # robust to CSV float round-trip
      cmp(pt$tar, tt$TAR, sprintf("fold %d %-15s TAR vs archived", f, m))
    }
  }
}
stop_if_fail("Self-test 2")
readr::write_csv(
  dplyr::bind_rows(lapply(RESULTS, as.data.frame)),
  file.path(OUT, "selfcheck_point_estimates.csv"))

# ============================================================ sparse-tail diagnostics
section("Acceptance counts at the fixed thresholds (sparse-tail diagnostics)")
acc_rows <- list()
for (f in seq_len(n_folds)) for (m in methods) {
  D <- FOLD[[f]]; M <- D$M[[m]]
  for (k in seq_along(fars)) {
    acc <- M$s_test >= M$thr[k]
    g <- acc & D$gen; i <- acc & !D$gen
    fp <- paste(pmin(D$F1[i], D$F2[i]), pmax(D$F1[i], D$F2[i]))
    acc_rows[[length(acc_rows) + 1]] <- tibble::tibble(
      fold = f, method = m, FAR = fars[k], threshold = unname(M$thr[k]),
      n_genuine = sum(D$gen), n_impostor = sum(!D$gen),
      accepted_genuine = sum(g), distinct_fingers_genuine = length(unique(D$F1[g])),
      accepted_impostor = sum(i), distinct_finger_pairs_impostor = length(unique(fp)))
  }
}
acc_df <- dplyr::bind_rows(acc_rows)
readr::write_csv(acc_df, file.path(OUT, "acceptance_counts.csv"))
sparse <- acc_df[acc_df$accepted_genuine <= 5, ]
cat(sprintf("  %d (fold, method, FAR) cells have <= 5 accepted genuine pairs; e.g.:\n", nrow(sparse)))
print(head(sparse[, c("fold","method","FAR","accepted_genuine","distinct_fingers_genuine","accepted_impostor")], 8))

# ============================================================ bootstrap
section(sprintf("Bootstrap: B = %d, seed = %d, finger-level multiplicities shared across methods", B, SEED))
set.seed(SEED)
nF <- vapply(FOLD, `[[`, integer(1), "n_fingers")
MULT <- lapply(seq_len(n_folds), function(f) {
  M <- matrix(0L, nrow = B, ncol = nF[f])
  for (b in seq_len(B)) M[b, ] <- tabulate(sample.int(nF[f], nF[f], replace = TRUE), nbins = nF[f])
  M
})
# result arrays: per method, [B, fold] for auc/eer; [B, fold, far] for tar/far
AUC <- EER <- lapply(methods, function(m) matrix(NA_real_, B, n_folds)); names(AUC) <- names(EER) <- methods
TAR <- FARr <- lapply(methods, function(m) array(NA_real_, c(B, n_folds, length(fars)))); names(TAR) <- names(FARr) <- methods
t0 <- Sys.time()
for (f in seq_len(n_folds)) {
  D <- FOLD[[f]]
  for (b in seq_len(B)) {
    mv <- MULT[[f]][b, ]
    w  <- mv[D$F1] * ifelse(D$gen, 1, mv[D$F2])
    for (m in methods) {
      r <- weighted_metrics(PREP[[f]][[m]], w)
      AUC[[m]][b, f] <- r$auc; EER[[m]][b, f] <- r$eer
      TAR[[m]][b, f, ] <- r$tar; FARr[[m]][b, f, ] <- r$far
    }
  }
  cat(sprintf("  fold %d done (%.1f min elapsed)\n", f, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# ============================================================ summaries
alpha <- (1 - LEVEL) / 2
ci <- function(x) stats::quantile(x, c(alpha, 1 - alpha), names = FALSE, na.rm = TRUE)

# replicate value of the reported quantity = mean over folds
rep_mean <- function(arr3_or_mat, k = NULL) if (is.null(k)) rowMeans(arr3_or_mat) else rowMeans(arr3_or_mat[, , k, drop = TRUE])

metric_table <- function(idx_B = seq_len(B)) {
  rows <- list()
  for (m in methods) {
    pf <- point[[m]]
    add <- function(metric, reps, pts) {
      q <- ci(reps[idx_B])
      rows[[length(rows) + 1]] <<- tibble::tibble(
        method = m, metric = metric,
        point = mean(pts), fold_sd = stats::sd(pts),
        ci_low = q[1], ci_high = q[2],
        degenerate = isTRUE(q[1] == q[2]))
    }
    add("AUC", rep_mean(AUC[[m]]), vapply(pf, `[[`, numeric(1), "auc"))
    add("EER", rep_mean(EER[[m]]), vapply(pf, `[[`, numeric(1), "eer"))
    for (k in seq_along(fars)) {
      add(sprintf("TAR@%g", fars[k]),          rep_mean(TAR[[m]],  k), vapply(pf, function(p) p$tar[k], numeric(1)))
      add(sprintf("realisedFAR@%g", fars[k]),  rep_mean(FARr[[m]], k), vapply(pf, function(p) p$far[k], numeric(1)))
    }
  }
  dplyr::bind_rows(rows)
}
ci_df <- metric_table()
readr::write_csv(ci_df, file.path(OUT, "ci_metrics.csv"))

# per-fold intervals (supplementary)
pf_rows <- list()
for (m in methods) for (f in seq_len(n_folds)) {
  addpf <- function(metric, reps, pt) { q <- ci(reps); pf_rows[[length(pf_rows) + 1]] <<- tibble::tibble(method = m, fold = f, metric = metric, point = pt, ci_low = q[1], ci_high = q[2]) }
  addpf("AUC", AUC[[m]][, f], point[[m]][[f]]$auc); addpf("EER", EER[[m]][, f], point[[m]][[f]]$eer)
  for (k in seq_along(fars)) {
    addpf(sprintf("TAR@%g", fars[k]),         TAR[[m]][, f, k],  point[[m]][[f]]$tar[k])
    addpf(sprintf("realisedFAR@%g", fars[k]), FARr[[m]][, f, k], point[[m]][[f]]$far[k])
  }
}
readr::write_csv(dplyr::bind_rows(pf_rows), file.path(OUT, "ci_metrics_per_fold.csv"))

# endpoint stability: first half of replicates vs all
half <- metric_table(seq_len(B %/% 2))
stab <- dplyr::inner_join(half, ci_df, by = c("method", "metric"), suffix = c(".half", ".full"))
stab <- dplyr::mutate(stab,
  d_low = abs(ci_low.full - ci_low.half), d_high = abs(ci_high.full - ci_high.half),
  width = ci_high.full - ci_low.full,
  rel_change = pmax(d_low, d_high) / ifelse(width > 0, width, NA_real_))
readr::write_csv(dplyr::select(stab, method, metric, ci_low.half, ci_high.half, ci_low.full, ci_high.full, width, rel_change),
                 file.path(OUT, "endpoint_stability.csv"))

# paired differences — pre-specified; POINTWISE intervals
pairs <- list(
  c("Fusion_eq", "V23_full"), c("Fusion_eq", "OT"),           # decisive for R4-9
  c("Fusion", "Fusion_eq"),                                   # reported; zero-in-interval is inconclusive
  c("Fusion", "V23_full"), c("Fusion", "OT"), c("V23_full", "OT"),
  c("GlobalBetti_L1", "DT_L2_150"), c("GlobalPI_L2", "DT_L2_150"),
  c("V23_full", "GlobalBetti_L1"), c("V23_noTopK", "V23_full"),
  # aligned geometry baseline (R4-6), pre-specified before the integration run:
  c("DT_L2_150_aligned", "DT_L2_150"),                        # does rigid alignment strengthen the geometric baseline?
  c("GlobalBetti_L1", "DT_L2_150_aligned"), c("GlobalPI_L2", "DT_L2_150_aligned"),   # both topological summaries vs the aligned geometry
  c("V23_full", "DT_L2_150_aligned"), c("OT", "DT_L2_150_aligned"), c("Fusion_eq", "DT_L2_150_aligned"),
  # matched-transformation geometry baseline (audit patch, 31 Oct 2026): the 150x150 fields clamped to [0, 20]
  # and rounded to 0.1 before RMS normalisation, as the topological methods receive them; present only in runs
  # made with the patched pipeline, skipped elsewhere:
  c("DT_L2_150_clamped", "DT_L2_150"),                                             # does the matched transformation change the geometry baseline?
  c("GlobalBetti_L1", "DT_L2_150_clamped"), c("GlobalPI_L2", "DT_L2_150_clamped"),   # the untrained summaries vs the matched baseline
  c("V23_full", "DT_L2_150_clamped"), c("OT", "DT_L2_150_clamped"), c("Fusion_eq", "DT_L2_150_clamped"),
  # reference matcher (R4-7, revision step 1.3), pre-specified in STEP_1.3_REFERENCE_MATCHER_SPEC.md;
  # present only in the run-shaped reference directories, skipped elsewhere:
  c("Fusion_eq", "NBIS_gray"), c("V23_full", "NBIS_gray"),                        # the gap, calibrated on our protocol
  c("Fus_NBISgray_FusionEq", "NBIS_gray"), c("Fus_NBISgray_TopoLR", "NBIS_gray"),  # does topology add to the minutiae matcher?
  c("Fus_NBISgray_FusionEq", "Fusion_eq"),                                         # does the matcher add to topology? (expected: yes)
  c("NBIS_bin", "NBIS_gray"))                                                      # the matcher on our binary maps vs on grayscale
pairs_skipped <- Filter(function(pr) !all(pr %in% methods), pairs)
if (length(pairs_skipped)) cat("  skipped (method absent):", paste(vapply(pairs_skipped, paste, character(1), collapse = " vs "), collapse = "; "), "\n")
pairs <- Filter(function(pr) all(pr %in% methods), pairs)
diff_rows <- list()
for (pr in pairs) {
  a <- pr[1]; b <- pr[2]
  addd <- function(metric, ra, rb, pa, pb) {
    d <- ra - rb; q <- ci(d)
    diff_rows[[length(diff_rows) + 1]] <<- tibble::tibble(
      method_a = a, method_b = b, metric = metric,
      difference = mean(pa) - mean(pb), ci_low = q[1], ci_high = q[2],
      excludes_zero = (q[1] > 0 | q[2] < 0), interval_type = "pointwise 95% percentile")
  }
  addd("AUC", rep_mean(AUC[[a]]), rep_mean(AUC[[b]]), vapply(point[[a]], `[[`, numeric(1), "auc"), vapply(point[[b]], `[[`, numeric(1), "auc"))
  addd("EER", rep_mean(EER[[a]]), rep_mean(EER[[b]]), vapply(point[[a]], `[[`, numeric(1), "eer"), vapply(point[[b]], `[[`, numeric(1), "eer"))
  for (k in seq_along(fars)) {
    addd(sprintf("TAR@%g", fars[k]), rep_mean(TAR[[a]], k), rep_mean(TAR[[b]], k),
         vapply(point[[a]], function(p) p$tar[k], numeric(1)), vapply(point[[b]], function(p) p$tar[k], numeric(1)))
    addd(sprintf("realisedFAR@%g", fars[k]), rep_mean(FARr[[a]], k), rep_mean(FARr[[b]], k),
         vapply(point[[a]], function(p) p$far[k], numeric(1)), vapply(point[[b]], function(p) p$far[k], numeric(1)))
  }
}
diff_df <- dplyr::bind_rows(diff_rows)
readr::write_csv(diff_df, file.path(OUT, "ci_paired_differences.csv"))

# fusion weights (tuned Fusion's per-fold w; Fusion_eq is w = 0.5 by definition)
fw <- tibble::tibble(fold = seq_len(n_folds), best_w = vapply(FOLD, function(D) as.numeric(D$M[["Fusion"]]$best_w), numeric(1)))
readr::write_csv(fw, file.path(OUT, "fusion_weights.csv"))
fe_rows <- lapply(seq_len(n_folds), function(f) {
  p <- point[["Fusion_eq"]][[f]]
  tibble::tibble(fold = f, AUC = p$auc, EER = p$eer,
                 !!!setNames(as.list(p$tar), paste0("TAR@", fars)),
                 !!!setNames(as.list(p$far), paste0("realisedFAR@", fars)))
})
readr::write_csv(dplyr::bind_rows(fe_rows), file.path(OUT, "fusion_eq_fold_metrics.csv"))
# NOTE: fusion_eq_scores.rds is no longer written here. The pre-integration reference is
# frozen in snapshot_before_new_methods/analysis/ and Fusion_eq now lives in the pipeline caches.

saveRDS(list(B = B, seed = SEED, level = LEVEL, multiplicities = MULT, AUC = AUC, EER = EER, TAR = TAR, realisedFAR = FARr,
             methods = methods, fars = fars, point = point),
        file.path(OUT, "bootstrap_replicates.rds"))
writeLines(c(
  sprintf("B = %d replicates; seed = %d; %.0f%% percentile intervals, POINTWISE (not simultaneous)", B, SEED, 100 * LEVEL),
  "Resampling unit: finger (20 per test fold), drawn with replacement; multiplicities m_i.",
  "Weights: genuine pair of finger i -> m_i; impostor pair (i != j) -> m_i * m_j; m_i = 0 -> dropped.",
  "No pair invented; copies of one finger never form an impostor pair; metrics normalised by total weight.",
  "Thresholds FIXED at training-calibrated values: quantile(train impostor scores, 1 - FAR), type 7; accept if score >= threshold.",
  "Same multiplicities for all methods (paired differences on identical resamples).",
  "Intervals are conditional on fitted models, thresholds and the fold split; headline interval is on the 5-fold mean.",
  "All methods, including Fusion_eq and DT_L2_150_aligned, are read from the pipeline's per-fold caches via fold_cache_index.csv.",
  "Fusion_eq: pipeline method, w = 0.5 with the training-impostor robust-z normalisation; checked against its components and the frozen pre-integration reference.",
  "DT_L2_150_aligned: rigid-search geometric baseline (coarse angular search with local refinement and exhaustive translation search); score = -d_sym.",
  "Tuned fusion weights optimise training TAR at nominal FAR 1e-3; held-out effects assessed empirically."),
  file.path(OUT, "bootstrap_settings.txt"))

# ============================================================ headline print
section("Headline (5-fold mean, 95% pointwise interval)")
show <- intersect(c("DT_L2_150", "DT_L2_150_aligned", "GlobalBetti_L1", "GlobalPI_L2", "V23_full", "OT", "Fusion", "Fusion_eq",
                    "NBIS_gray", "NBIS_bin", "Fus_NBISgray_TopoLR", "Fus_NBISgray_FusionEq"), methods)
hd <- ci_df[ci_df$method %in% show & ci_df$metric %in% c("AUC", "EER", "TAR@0.001", "realisedFAR@0.001", "TAR@0.01"), ]
hd$method <- factor(hd$method, levels = show); hd <- hd[order(hd$method, hd$metric), ]
print(as.data.frame(dplyr::mutate(hd, dplyr::across(c(point, fold_sd, ci_low, ci_high), ~ round(.x, 4)))), row.names = FALSE)

section("Fusion weights and the R4-9 comparisons")
cat(sprintf("  tuned w per fold: %s  (mean %.3f, SD %.3f, range %.2f-%.2f)\n",
            paste(sprintf("%.2f", fw$best_w), collapse = ", "), mean(fw$best_w), sd(fw$best_w), min(fw$best_w), max(fw$best_w)))
key <- diff_df[diff_df$method_a %in% c("Fusion_eq", "Fusion") & diff_df$method_b %in% c("V23_full", "OT", "Fusion_eq"), ]
cat("  (TAR differences are read together with the realisedFAR differences at the same FAR:\n")
cat("   a TAR gain alongside a higher realised FAR is partly a less stringent operating point.)\n")
print(as.data.frame(dplyr::mutate(key, dplyr::across(c(difference, ci_low, ci_high), ~ round(.x, 4)))), row.names = FALSE)

if ("DT_L2_150_aligned" %in% methods) {
  section("Aligned geometry baseline (R4-6): pre-specified paired differences")
  al <- diff_df[(diff_df$method_a == "DT_L2_150_aligned" | diff_df$method_b == "DT_L2_150_aligned") &
                diff_df$metric %in% c("AUC", "EER", "TAR@0.001", "realisedFAR@0.001", "TAR@0.01"), ]
  print(as.data.frame(dplyr::mutate(al, dplyr::across(c(difference, ci_low, ci_high), ~ round(.x, 4)))), row.names = FALSE)
}

if ("NBIS_gray" %in% methods) {
  section("Reference matcher (R4-7): pre-specified paired differences")
  rf <- diff_df[(diff_df$method_a %in% c("NBIS_gray", "NBIS_bin", "Fus_NBISgray_TopoLR", "Fus_NBISgray_FusionEq") |
                 diff_df$method_b %in% c("NBIS_gray", "NBIS_bin")) &
                diff_df$metric %in% c("AUC", "EER", "TAR@0.001", "realisedFAR@0.001", "TAR@0.01", "realisedFAR@0.01"), ]
  print(as.data.frame(dplyr::mutate(rf, dplyr::across(c(difference, ci_low, ci_high), ~ round(.x, 4)))), row.names = FALSE)
}

section("Degenerate intervals (lo == hi) — a [0,0] interval is NOT evidence of a zero rate")
dg <- ci_df[ci_df$degenerate, c("method", "metric", "point", "ci_low", "ci_high")]
if (nrow(dg)) print(as.data.frame(dg), row.names = FALSE) else cat("  none\n")

section("Endpoint stability (B/2 vs B): cells with relative endpoint change > 10% of interval width")
us <- stab[!is.na(stab$rel_change) & stab$rel_change > 0.10, c("method", "metric", "ci_low.half", "ci_low.full", "ci_high.half", "ci_high.full", "rel_change")]
if (nrow(us)) print(as.data.frame(dplyr::mutate(us, dplyr::across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE) else cat("  none\n")

cat("\nOutputs written to", OUT, "\n")
