############################################################
# VERIFY A NEW RUN AGAINST A CANONICAL RUN (audit patch, 31 Oct 2026; patch 1b: strict mode).
# Generalises Code/verify_regen_reproduction.R: the new run (B) may carry ADDITIONAL methods and
# additional bank fields (here DT_L2_150_clamped / dt_vec_150_clamped); everything the two runs
# share must be identical, object by object.
#
#   options(fp.verify_a = "Outputs_regen", fp.verify_b = "Outputs_regen_v2"); source("Code/verify_run_vs_canonical.R")
#   options(fp.verify_expect_field = "dt_vec_150_clamped", fp.verify_expect_method = "DT_L2_150_clamped")   # defaults; "" = no expectation
#
# Checks:
#   A  manifests: same files, md5s, identities, folds and pair orderings; same EFFECTIVE configuration (output_dir and
#      the operational fields n_workers / future_globals_gb excluded; complex_type / polarity read as defaults when absent)
#   B  feature banks: every field of A's bank entries identical in B's (EXACT); the expected new field present in
#      every entry of B and sanity-checked (length 2*target_size^2, finite, unit RMS); input sidecars identical
#   C  per (fold, method) shared by both: s_train, s_test, labels, thresholds, TopoLR model, fusion fits EXACT;
#      B computed everything itself (nothing loaded from cache); the expected new method present in every fold
#      with finite scores of the right length (its fold AUCs are printed)
#   D  aligned pair caches: same key, identical scores and diagnostics (EXACT)
#   E  output CSVs on the shared methods: text columns identical, numeric columns identical up to 1e-9 after a
#      CSV round trip (TOLERANCE is accepted here and ONLY here); scores_raw pooled scores EXACT
# Core objects (B, C, D, scores_raw) are compared bit-for-bit: any difference is a FAIL, whatever its size.
# The report is saved to <B>/verification_vs_<A>.rds and then the script STOPS with an error if any check
# failed, so a sourcing driver cannot continue past a failed verification.
############################################################

# Dependencies: base R and readr only (the verifier deliberately does not source the pipeline: it compares
# stored objects and must not depend on the code under test).
if (!requireNamespace("readr", quietly = TRUE)) stop("package readr is required", call. = FALSE)

A_DIR <- getOption("fp.verify_a", "Outputs_regen")
B_DIR <- getOption("fp.verify_b", "Outputs_regen_v2")
EXPECT_FIELD  <- getOption("fp.verify_expect_field",  "dt_vec_150_clamped")
EXPECT_METHOD <- getOption("fp.verify_expect_method", "DT_L2_150_clamped")
CSV_TOL <- 1e-9
stopifnot(dir.exists(A_DIR), dir.exists(B_DIR), file.exists(file.path(A_DIR, "run_manifest.rds")), file.exists(file.path(B_DIR, "run_manifest.rds")))

RESULTS <- list()
report <- function(label, status, note = "") {
  RESULTS[[length(RESULTS) + 1]] <<- list(label = label, status = status, note = note)
  cat(sprintf("  [%-9s] %-70s %s\n", status, label, note)); invisible(status)
}
# exact comparison of numeric content (core objects): EXACT or FAIL, never TOLERANCE
cmp_exact <- function(a, b, label) {
  a <- as.numeric(unlist(a)); b <- as.numeric(unlist(b))
  if (length(a) != length(b)) return(report(label, "FAIL", sprintf("length %d vs %d", length(a), length(b))))
  if (length(a) == 0) return(report(label, "FAIL", "empty"))
  if (any(!is.finite(a)) || any(!is.finite(b))) return(report(label, "FAIL", "non-finite"))
  if (all(a == b)) report(label, "EXACT", "") else report(label, "FAIL", sprintf("%d of %d values differ, max|diff| = %.3g", sum(a != b), length(a), max(abs(a - b))))
}
# comparison after a CSV round trip (section E only): EXACT, TOLERANCE (<= CSV_TOL) or FAIL
cmp_csv <- function(a, b, label, tol = CSV_TOL) {
  a <- as.numeric(unlist(a)); b <- as.numeric(unlist(b))
  if (length(a) != length(b)) return(report(label, "FAIL", sprintf("length %d vs %d", length(a), length(b))))
  if (length(a) == 0) return(report(label, "FAIL", "empty"))
  na_a <- !is.finite(a); na_b <- !is.finite(b)
  if (!identical(na_a, na_b)) return(report(label, "FAIL", "non-finite pattern differs"))
  a <- a[!na_a]; b <- b[!na_b]
  if (all(a == b)) return(report(label, "EXACT", ""))
  m <- max(abs(a - b)); if (m <= tol) report(label, "TOLERANCE", sprintf("max|diff| = %.3g", m)) else report(label, "FAIL", sprintf("max|diff| = %.3g", m))
}
chk <- function(cond, label, note = "") report(label, if (isTRUE(cond)) "EXACT" else "FAIL", note)
section <- function(s) cat("\n", s, "\n", strrep("-", nchar(s)), "\n", sep = "")
rd <- function(d, f) readr::read_csv(file.path(d, f), show_col_types = FALSE, progress = FALSE)
# sort a data frame by whichever of the candidate key columns it has; the key list is UNNAMED so that a column
# called "method" cannot bind order()'s own `method` argument
sort_by_keys <- function(d, candidates = c("method", "fold", "FAR")) {
  keys <- intersect(candidates, names(d)); if (!length(keys)) return(d)
  d[do.call(order, unname(as.list(d[keys]))), , drop = FALSE]
}

section("A. Manifests and configuration")
mA <- readRDS(file.path(A_DIR, "run_manifest.rds")); mB <- readRDS(file.path(B_DIR, "run_manifest.rds"))
chk(identical(mA$files, mB$files) && identical(mA$md5, mB$md5), "same input files and md5s")
chk(identical(mA$ids, mB$ids) && identical(mA$folds, mB$folds), "same identities, folds and pair orderings")
# configuration: compared as EFFECTIVE scientific settings — output_dir and the operational fields (worker count,
# future globals limit) are dropped; complex_type / polarity, absent from pre-1.5 manifests, are read as their
# defaults ("triangulated", "auto"); fields are compared by name, not position (same rule as Code/run_canonical.R)
`%||%` <- function(a, b) if (is.null(a)) b else a
OPERATIONAL_FIELDS <- c("n_workers", "future_globals_gb")
norm_cfg <- function(cfg) { cfg$output_dir <- NULL; for (f in OPERATIONAL_FIELDS) cfg[[f]] <- NULL
  cfg$complex_type <- cfg$complex_type %||% "triangulated"; cfg$polarity <- cfg$polarity %||% "auto"; cfg[order(names(cfg))] }
cA <- norm_cfg(mA$config); cB <- norm_cfg(mB$config)
cfg_diff <- c(names(cB)[!vapply(names(cB), function(n) identical(cA[[n]], cB[[n]]), logical(1))], setdiff(names(cA), names(cB)))
chk(identical(cA, cB), "same effective configuration (apart from output_dir and operational fields)", if (length(cfg_diff)) paste("differs in:", paste(cfg_diff, collapse = ", ")) else "")
for (f in OPERATIONAL_FIELDS) if (!identical(mA$config[[f]], mB$config[[f]])) cat(sprintf("  (operational field %s differs: %s vs %s; not a scientific setting)\n", f, format(mA$config[[f]]), format(mB$config[[f]])))

section("B. Feature banks")
bank_of <- function(d) { f <- list.files(file.path(d, "cache"), pattern = "^bank_[0-9a-f]{12}\\.rds$", full.names = TRUE); stopifnot(length(f) == 1); f }
bA <- bank_of(A_DIR); bB <- bank_of(B_DIR)
chk(basename(bA) == basename(bB), "same bank key", basename(bA))
bankA <- readRDS(bA); bankB <- readRDS(bB)
chk(length(bankA) == length(bankB) && length(bankB) == length(mB$files), "same number of bank entries (= number of images)", sprintf("%d", length(bankB)))
fA <- names(bankA[[1]]); fB <- names(bankB[[1]])
chk(all(fA %in% fB), "every field of A present in B", if (!all(fA %in% fB)) paste("missing:", paste(setdiff(fA, fB), collapse = ", ")) else "")
chk(all(vapply(bankB, function(e) identical(names(e), fB), logical(1))), "every B entry has the same field names")
extra <- setdiff(fB, fA)
cat("  B-only bank fields:", if (length(extra)) paste(extra, collapse = ", ") else "(none)", "\n")
n_pairs <- min(length(bankA), length(bankB))
same_shared <- all(vapply(seq_len(n_pairs), function(i) all(vapply(fA, function(fn) identical(bankA[[i]][[fn]], bankB[[i]][[fn]]), logical(1))), logical(1)))
chk(same_shared, "all shared bank fields identical in every entry (bit-for-bit)")
if (!same_shared) for (fn in fA) {
  d <- max(vapply(seq_len(n_pairs), function(i) { a <- bankA[[i]][[fn]]; b <- bankB[[i]][[fn]]
    if (is.null(a) || is.null(b) || length(a) != length(b)) Inf else max(abs(as.numeric(a) - as.numeric(b))) }, numeric(1)))
  cat(sprintf("     field %-22s max |A - B| = %s\n", fn, format(d)))
}
if (nzchar(EXPECT_FIELD)) {
  present <- all(vapply(bankB, function(e) !is.null(e[[EXPECT_FIELD]]), logical(1)))
  chk(present, sprintf("expected new bank field %s present in every B entry", EXPECT_FIELD))
  chk(!(EXPECT_FIELD %in% fA), sprintf("expected new field %s is NEW (absent from A)", EXPECT_FIELD))
  if (present && EXPECT_FIELD == "dt_vec_150_clamped") {
    ts <- mB$config$target_size
    ok_len <- all(vapply(bankB, function(e) length(e$dt_vec_150_clamped) == 2 * ts^2, logical(1)))
    ok_fin <- all(vapply(bankB, function(e) all(is.finite(e$dt_vec_150_clamped)), logical(1)))
    rms <- vapply(bankB, function(e) sqrt(mean(e$dt_vec_150_clamped^2)), numeric(1))
    chk(ok_len && ok_fin, "dt_vec_150_clamped: length 2*target_size^2 and finite in every entry")
    chk(max(abs(rms - 1)) < 1e-9, "dt_vec_150_clamped: unit RMS in every entry", sprintf("max |rms - 1| = %.2g", max(abs(rms - 1))))
    # (the arithmetic of the field against the recomputed 150x150 distance fields is checked by Code/test_hardening.R T7)
    chk(!all(vapply(bankB, function(e) identical(e$dt_vec_150_clamped, e$dt_vec_150), logical(1))), "dt_vec_150_clamped differs from dt_vec_150 (clipping occurred somewhere)")
  }
}
chk(identical(readRDS(sub("bank_", "inputs_bank_", bA)), readRDS(sub("bank_", "inputs_bank_", bB))), "bank input sidecars identical")
rm(bankA, bankB); invisible(gc())                     # the banks are large; the remaining sections read one object pair at a time

section("C. Per-fold caches")
idx_of <- function(d) readr::read_csv(file.path(d, "fold_cache_index.csv"), show_col_types = FALSE,
                                      col_types = readr::cols(.default = readr::col_character(), fold = readr::col_integer(), loaded_from_cache = readr::col_logical()))
iA <- idx_of(A_DIR); iB <- idx_of(B_DIR)
keyA <- paste(iA$fold, iA$method); keyB <- paste(iB$fold, iB$method)
chk(all(keyA %in% keyB), "every (fold, method) of A present in B", sprintf("%d of %d rows", sum(keyA %in% keyB), nrow(iA)))
chk(!anyDuplicated(keyB), "no duplicated (fold, method) rows in B's cache index")
chk(!any(iB$loaded_from_cache), "B computed every method itself (nothing loaded from cache)")
b_only <- setdiff(unique(iB$method), unique(iA$method))
cat("  B-only methods:", if (length(b_only)) paste(b_only, collapse = ", ") else "(none)", "\n")
for (k in seq_len(nrow(iA))) {
  f <- iA$fold[k]; m <- iA$method[k]
  rowB <- iB[iB$fold == f & iB$method == m, ]
  if (nrow(rowB) != 1) { report(sprintf("fold %d %-18s present in B", f, m), "FAIL"); next }
  ca <- readRDS(file.path(A_DIR, "cache", "folds", iA$file[k])); cb <- readRDS(file.path(B_DIR, "cache", "folds", rowB$file))
  cmp_exact(ca$s_train, cb$s_train, sprintf("fold %d %-18s s_train", f, m)); cmp_exact(ca$s_test, cb$s_test, sprintf("fold %d %-18s s_test", f, m))
  chk(identical(ca$y_train, cb$y_train) && identical(ca$y_test, cb$y_test), sprintf("fold %d %-18s labels", f, m))
  cmp_exact(ca$thresholds, cb$thresholds, sprintf("fold %d %-18s thresholds", f, m))
  if (!is.null(ca$model)) { cmp_exact(ca$model$w, cb$model$w, sprintf("fold %d %-18s model w", f, m)); cmp_exact(ca$model$b, cb$model$b, sprintf("fold %d %-18s model b", f, m))
                            cmp_exact(unlist(ca$model$standardizer), unlist(cb$model$standardizer), sprintf("fold %d %-18s standardizer", f, m)) }
  if (!is.null(ca$best_w)) { cmp_exact(ca$best_w, cb$best_w, sprintf("fold %d %-18s fusion weight", f, m)); cmp_exact(unlist(ca$zfit_ot), unlist(cb$zfit_ot), sprintf("fold %d %-18s zfit_ot", f, m)); cmp_exact(unlist(ca$zfit_v), unlist(cb$zfit_v), sprintf("fold %d %-18s zfit_v", f, m)) }
  if (!is.null(ca$aligned_key)) chk(identical(ca$aligned_key, cb$aligned_key), sprintf("fold %d %-18s aligned key", f, m))
}
if (nzchar(EXPECT_METHOD)) {
  chk(EXPECT_METHOD %in% b_only, sprintf("expected new method %s present in B and absent from A", EXPECT_METHOD))
  n_folds_B <- length(mB$folds)
  chk(sum(iB$method == EXPECT_METHOD) == n_folds_B, sprintf("%s present in all %d folds of B", EXPECT_METHOD, n_folds_B))
}
if (length(b_only)) {
  aucB <- rd(B_DIR, "fold_metrics_auc_eer.csv")
  for (m in b_only) {
    r <- aucB[aucB$method == m, ]; r <- r[order(r$fold), ]
    cat(sprintf("  B-only method %-18s fold AUCs: %s  (mean %.4f)\n", m, paste(sprintf("%.7f", r$AUC), collapse = " "), mean(r$AUC)))
    for (f in sort(unique(iB$fold[iB$method == m]))) {
      cb <- readRDS(file.path(B_DIR, "cache", "folds", iB$file[iB$fold == f & iB$method == m]))
      chk(all(is.finite(cb$s_train)) && all(is.finite(cb$s_test)) && length(cb$s_test) == ncol(mB$folds[[f]]$pair_test) && length(cb$s_train) == ncol(mB$folds[[f]]$pair_train),
          sprintf("%-18s fold %d scores finite, right lengths", m, f))
    }
  }
}

section("D. Aligned pair caches")
alA <- list.files(file.path(A_DIR, "cache"), pattern = "^aligned_[0-9a-f]{12}\\.rds$", full.names = TRUE)
alB <- list.files(file.path(B_DIR, "cache"), pattern = "^aligned_[0-9a-f]{12}\\.rds$", full.names = TRUE)
chk(length(alA) == 1 && length(alB) == 1 && basename(alA) == basename(alB), "one aligned pair cache each, same key")
if (length(alA) == 1 && length(alB) == 1) {
  a <- readRDS(alA); b <- readRDS(alB)
  chk(identical(a$key, b$key) && identical(a$files, b$files) && identical(a$pairs, b$pairs), "aligned cache identity (key, files, pairs)")
  cmp_exact(a$score, b$score, "aligned scores")
  num <- names(a$diag)[vapply(a$diag, is.numeric, logical(1))]
  for (col in num) cmp_exact(a$diag[[col]], b$diag[[col]], paste("aligned diag", col))
}

section("E. Output CSVs (shared methods)")
shared <- intersect(unique(iA$method), unique(iB$method))
for (f in c("fold_metrics_auc_eer.csv", "fold_metrics_tar.csv", "summary_auc_eer_mean_sd.csv", "summary_tar_mean_sd.csv", "aligned_boundary_rates.csv")) {
  if (!file.exists(file.path(A_DIR, f)) || !file.exists(file.path(B_DIR, f))) { report(paste(f, "present in both runs"), "FAIL"); next }
  a <- rd(A_DIR, f); b <- rd(B_DIR, f)
  if ("method" %in% names(a)) { a <- a[a$method %in% shared, ]; b <- b[b$method %in% shared, ] }
  a <- sort_by_keys(a); b <- sort_by_keys(b)
  if (!identical(names(a), names(b)) || nrow(a) != nrow(b)) { report(paste(f, "same columns and row count (shared methods)"), "FAIL", sprintf("%d x %d vs %d x %d", nrow(a), ncol(a), nrow(b), ncol(b))); next }
  numc <- names(a)[vapply(a, is.numeric, logical(1))]; chrc <- setdiff(names(a), c(numc, "aligned_key"))
  ok_chr <- all(vapply(chrc, function(cn) identical(as.character(a[[cn]]), as.character(b[[cn]])), logical(1)))
  if (length(numc)) cmp_csv(unlist(a[numc]), unlist(b[numc]), paste(f, "numeric columns (shared methods)"))
  chk(ok_chr, paste(f, "rows and text columns (shared methods)"))
}
cat("  reading scores_raw.csv (both runs) ...\n")
a <- rd(A_DIR, "scores_raw.csv"); b <- rd(B_DIR, "scores_raw.csv"); b <- b[b$method %in% shared, ]
a <- a[a$method %in% shared, ]
chk(nrow(a) == nrow(b) && identical(a$method, b$method) && identical(a$fold, b$fold) && identical(a$y, b$y), "scores_raw.csv structure identical on the shared methods")
if (nrow(a) == nrow(b)) cmp_exact(a$score, b$score, "scores_raw.csv pooled test scores (shared methods)")

section("Summary")
st <- vapply(RESULTS, `[[`, character(1), "status"); print(table(factor(st, levels = c("EXACT", "TOLERANCE", "FAIL"))))
fails <- RESULTS[st == "FAIL"]; if (length(fails)) { cat("\nFAILED:\n"); for (x in fails) cat("  -", x$label, x$note, "\n") }
out <- file.path(B_DIR, paste0("verification_vs_", basename(A_DIR), ".rds"))
saveRDS(list(created = Sys.time(), canonical = A_DIR, new_run = B_DIR, b_only_methods = b_only, b_only_bank_fields = extra,
             expected_field = EXPECT_FIELD, expected_method = EXPECT_METHOD, n_fail = length(fails), results = RESULTS), out)
cat("\nReport saved to", out, "\n")
if (length(fails)) stop(length(fails), " verification check(s) FAILED (see report); the new run must not be used.", call. = FALSE)
cat("VERIFICATION PASSED:", sum(st == "EXACT"), "exact,", sum(st == "TOLERANCE"), "within CSV tolerance, 0 failures.\n")
