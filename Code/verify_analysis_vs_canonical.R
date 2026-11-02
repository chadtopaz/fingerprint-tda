############################################################
# VERIFY A NEW BOOTSTRAP ANALYSIS AGAINST A CANONICAL ONE (audit patch, 31 Oct 2026).
# The new analysis directory (B) may contain additional methods and additional paired comparisons;
# every row the two analyses share (method x metric; method_a x method_b x metric) must be identical,
# because the finger resamples depend only on the seed and the fold structure and the shared scores
# are bit-identical (verified by Code/verify_run_vs_canonical.R).
#
#   options(fp.verify_a = "Outputs_reference_db1", fp.verify_b = "Outputs_reference_db1_v2"); source("Code/verify_analysis_vs_canonical.R")
#
# Writes only <B>/analysis/verification_vs_<A>.txt.
############################################################

A_DIR <- getOption("fp.verify_a", "Outputs_reference_db1")
B_DIR <- getOption("fp.verify_b", "Outputs_reference_db1_v2")
rd <- function(d, f) readr::read_csv(file.path(d, "analysis", f), show_col_types = FALSE, progress = FALSE)
out_lines <- character(0)
say <- function(...) { s <- sprintf(...); cat(s, "\n"); out_lines <<- c(out_lines, s) }
n_fail <- 0L
compare <- function(file, keys) {
  a <- rd(A_DIR, file); b <- rd(B_DIR, file)
  ka <- do.call(paste, c(a[keys], sep = "|")); kb <- do.call(paste, c(b[keys], sep = "|"))
  shared <- intersect(ka, kb); only_b <- setdiff(kb, ka); only_a <- setdiff(ka, kb)
  a <- a[match(shared, ka), ]; b <- b[match(shared, kb), ]
  numc <- names(a)[vapply(a, is.numeric, logical(1))]
  bad <- 0L
  for (cn in numc) {
    x <- a[[cn]]; y <- b[[cn]]
    same <- (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
    if (!all(same)) { bad <- bad + sum(!same); say("  MISMATCH %s column %s: %d of %d shared rows differ (max |diff| %.3g)", file, cn, sum(!same), length(same), max(abs(x - y), na.rm = TRUE)) }
  }
  chr <- setdiff(names(a), numc)
  for (cn in chr) { x <- as.character(a[[cn]]); y <- as.character(b[[cn]]); same <- (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
    if (!all(same)) { bad <- bad + sum(!same); say("  MISMATCH %s text column %s: %d rows differ", file, cn, sum(!same)) } }
  say("%-32s shared rows %5d  identical: %s   B-only rows: %d   A-only rows: %d", file, length(shared), if (bad == 0) "YES" else "NO", length(only_b), length(only_a))
  if (length(only_a)) say("    A-only keys (unexpected): %s", paste(head(only_a, 5), collapse = "; "))
  if (bad > 0 || length(only_a)) n_fail <<- n_fail + 1L
  invisible(list(only_b = only_b))
}
say("Analysis comparison: A = %s, B = %s", A_DIR, B_DIR)
r1 <- compare("ci_metrics.csv", c("method", "metric"))
r2 <- compare("ci_paired_differences.csv", c("method_a", "method_b", "metric"))
r3 <- compare("ci_metrics_per_fold.csv", c("method", "fold", "metric"))
r4 <- compare("acceptance_counts.csv", c("fold", "method", "FAR"))
r5 <- compare("endpoint_stability.csv", c("method", "metric"))
# settings
sa <- readLines(file.path(A_DIR, "analysis", "bootstrap_settings.txt")); sb <- readLines(file.path(B_DIR, "analysis", "bootstrap_settings.txt"))
say("bootstrap_settings.txt identical: %s", if (identical(sa, sb)) "YES" else "NO (see files)")
if (!identical(sa, sb)) n_fail <- n_fail + 1L
# new rows, summarised
newm <- unique(sub("\\|.*$", "", r1$only_b)); say("B-only methods in ci_metrics: %s", if (length(newm)) paste(newm, collapse = ", ") else "(none)")
if (length(r2$only_b)) say("B-only paired comparisons: %s", paste(unique(sub("\\|[^|]*$", "", r2$only_b)), collapse = "; "))
if (length(newm)) {
  b <- rd(B_DIR, "ci_metrics.csv"); b <- b[b$method %in% newm & b$metric %in% c("AUC", "EER", "TAR@0.001", "realisedFAR@0.001", "TAR@0.01"), ]
  for (i in seq_len(nrow(b))) say("    %-18s %-18s %.4f [%.4f, %.4f]", b$method[i], b$metric[i], b$point[i], b$ci_low[i], b$ci_high[i])
  d <- rd(B_DIR, "ci_paired_differences.csv"); d <- d[(d$method_a %in% newm | d$method_b %in% newm) & d$metric %in% c("AUC", "EER", "TAR@0.001", "TAR@0.01"), ]
  for (i in seq_len(nrow(d))) say("    %-16s - %-18s %-10s %+.4f [%+.4f, %+.4f] excludes zero: %s", d$method_a[i], d$method_b[i], d$metric[i], d$difference[i], d$ci_low[i], d$ci_high[i], d$excludes_zero[i])
}
say("RESULT: %s", if (n_fail == 0L) "all shared analysis rows identical" else sprintf("%d file(s) with differences", n_fail))
writeLines(out_lines, file.path(B_DIR, "analysis", paste0("verification_vs_", basename(A_DIR), ".txt")))
if (n_fail > 0L) stop("analysis comparison found differences on shared rows", call. = FALSE)
