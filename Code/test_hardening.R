############################################################
# REGRESSION TESTS FOR THE AUDIT PATCH (31 Oct 2026). Run from the project root:
#   source("Code/test_hardening.R")
# Writes nothing outside a temporary directory. Every check prints PASS/FAIL and the script stops
# at the end if any check failed.
#
#   T1  a NULL diagram is rejected by the vectorisers (never turned into zero features)
#   T2  a forced TDA::gridDiag failure propagates out of compute_persistence() with context
#   T3  ... and out of extract_feature_bank_one() / precompute_feature_bank(): the failing image is
#       named, load_or_compute_bank() stops, and NO bank file is written to the (temporary) output dir
#   T4  a valid diagram with no H1 features is accepted (zero Betti curve, zero persistence image)
#   T5  sample_anchor_points() gives the same anchors whatever RNG kind, normal kind and SAMPLE KIND the
#       session is in (Mersenne-Twister, L'Ecuyer-CMRG, sample.kind = "Rounding", no .Random.seed), and
#       leaves the session's RNG kind and state unchanged; control: an unpinned sample.kind would change the draw
#   T6  the anchors recovered from the canonical bank (fp.test_run, default Outputs_regen_v2, if present) equal those drawn by
#       sample_anchor_points() in this serial session, for several images, as a plain integer matrix
#       without attributes (identical()): figures/validators now reproduce the analysis anchors
#   T7  dt_vec_150_clamped of a freshly extracted entry equals the clamp-round-normalise of the
#       150x150 fields recomputed here, and the entry's other fields equal the canonical bank's
#   T8  (patch 1b) verifier failure injection: a symlink mirror of the reference run passes
#       Code/verify_run_vs_canonical.R with 0 failures; the same mirror with ONE test score of ONE fold cache
#       perturbed by 1e-12 is rejected (the script stops after saving its report, which names the object)
#   T9  (patch 1b) figure overlay coordinates: Code/make_figures.R's pixel_xy() equals the image() cell
#       centres for every cell of 300x300, 40x60 and 7x3 matrices, and a deliberately wrong mapping
#       ((i-1)/(n-1), the image() default) is caught by check_pixel_coordinates()
# T8 and T9 use no TDA computation; T8 needs the reference run and symlink support (macOS/Linux).
############################################################

options(fp.suppress_main = TRUE)
source("Code/fingerprintTDA.R")

FAILS <- 0L
check <- function(cond, label) { ok <- isTRUE(cond); cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", label)); if (!ok) FAILS <<- FAILS + 1L; invisible(ok) }
err_of <- function(expr) tryCatch({ expr; NULL }, error = function(e) conditionMessage(e))

cat("\nT1 NULL diagrams are rejected\n")
check(!is.null(err_of(compute_betti_curve(NULL, 1, 50, 20))), "compute_betti_curve(NULL) stops")
check(!is.null(err_of(diag_to_pi_vec(NULL, 1, 20, 20, 0.5))),  "diag_to_pi_vec(NULL) stops")

cat("\nT2 forced backend failure propagates from compute_persistence()\n")
orig_gridDiag <- TDA::gridDiag
force_fail <- function(...) stop("forced failure (test)")
assignInNamespace("gridDiag", force_fail, ns = "TDA")
m <- err_of(compute_persistence(matrix(runif(100), 10, 10), maxdim = 1, maxscale = 20))
check(!is.null(m) && grepl("forced failure", m, fixed = TRUE) && grepl("triangulated persistence", m, fixed = TRUE), paste("error propagated with context:", substr(m %||% "", 1, 80)))

cat("\nT3 ... and from the bank precompute, with the image named and no bank written\n")
# The forced gridDiag failure is visible only in this session (the precompute's workers load TDA afresh),
# so the per-image propagation is tested here directly and the precompute's aggregation/stop path with a
# failure that every worker sees (an unreadable input): the recorded error, the named image, no bank file.
CANON <- "Data/FVC2000_binary_regen"
if (dir.exists(CANON)) {
  cfg <- CONFIG; cfg$binary_dir <- CANON; cfg$output_dir <- tempfile("hardening_test_"); cfg$n_workers <- 2L; cfg$export_to_manuscript <- FALSE
  paths <- head(get_paths(cfg), 2)
  m1 <- err_of(extract_feature_bank_one(paths[1], cfg, make_disk_mask(cfg$ot$radius)))
  check(!is.null(m1) && grepl("forced failure", m1, fixed = TRUE), "extract_feature_bank_one() stops on the forced persistence failure")
  assignInNamespace("gridDiag", orig_gridDiag, ns = "TDA")
  bad <- c(paths[1], file.path(CANON, "does_not_exist_bin.png"))
  m2 <- err_of(load_or_compute_bank(bad, cfg))
  check(!is.null(m2) && grepl("failed precompute", m2, fixed = TRUE) && grepl("does_not_exist_bin.png", m2, fixed = TRUE), paste("load_or_compute_bank() stops and names the failing image:", substr(m2 %||% "", 1, 100)))
  check(length(list.files(file.path(cfg$output_dir, "cache"), pattern = "^bank_")) == 0, "no bank file written after the failure")
  future::plan(future::sequential)
} else cat("  (skipped: ", CANON, " not present)\n", sep = "")
assignInNamespace("gridDiag", orig_gridDiag, ns = "TDA")
check(identical(TDA::gridDiag, orig_gridDiag), "TDA::gridDiag restored")

cat("\nT4 a valid diagram without H1 features is accepted\n")
d0 <- compute_persistence(matrix(0, 10, 10), maxdim = 1, maxscale = 20)      # constant field: H0 only
check(!is.null(d0) && is.matrix(as.matrix(d0)), "constant field gives a diagram")
check(all(compute_betti_curve(d0, 1, 50, 20) == 0), "its H1 Betti curve is zero")
check(all(diag_to_pi_vec(d0, 1, 20, 20, 0.5) == 0), "its H1 persistence image is zero")

cat("\nT5 anchors independent of the session RNG kind, normal kind and sample kind; session RNG untouched\n")
set.seed(1); bin_ds <- matrix(as.integer(runif(150 * 150) < 0.3), 150, 150)
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(99); state_before <- .Random.seed; kind_before <- RNGkind()
a_mt <- sample_anchor_points(bin_ds, 64, 12345L)
check(identical(RNGkind(), kind_before) && identical(.Random.seed, state_before), "session RNG kind and state restored after the draw (Mersenne-Twister session)")
RNGkind("L'Ecuyer-CMRG"); set.seed(7); a_le <- sample_anchor_points(bin_ds, 64, 12345L)
check(identical(a_mt, a_le), "same anchors from a Mersenne-Twister session and an L'Ecuyer-CMRG session")
suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding")); set.seed(3); kind_r <- RNGkind(); state_r <- .Random.seed
a_rd <- sample_anchor_points(bin_ds, 64, 12345L)
check(identical(a_mt, a_rd), "same anchors from a session with sample.kind = 'Rounding' (the sample kind is pinned inside the helper)")
check(identical(RNGkind(), kind_r) && identical(.Random.seed, state_r), "session RNG kind (incl. sample.kind = 'Rounding') and state restored after the draw")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
# control: without the pin, 'Rounding' changes the draw, so the pin is doing real work
RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(12345L); idx0 <- which(bin_ds == 1L, arr.ind = TRUE); p1 <- idx0[sample.int(nrow(idx0), 64), ]
suppressWarnings(RNGkind("L'Ecuyer-CMRG", "Inversion", "Rounding")); set.seed(12345L); p2 <- idx0[sample.int(nrow(idx0), 64), ]
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
check(!identical(p1, p2), "control: an unpinned sample.kind = 'Rounding' would have changed the anchors")
if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
a_ns <- sample_anchor_points(bin_ds, 64, 12345L)
check(identical(a_mt, a_ns) && !exists(".Random.seed", envir = globalenv(), inherits = FALSE), "same anchors from a session without .Random.seed, and none is left behind")

cat("\nT6 figure/validator anchors == analysis anchors (canonical bank)\n")
REF <- getOption("fp.test_run", "Outputs_regen_v2")      # the reference DB1 run (any complete run of the pipeline)
bank_file <- list.files(file.path(REF, "cache"), pattern = "^bank_[0-9a-f]{12}\\.rds$", full.names = TRUE)
if (length(bank_file) == 1 && file.exists(file.path(REF, "run_manifest.rds"))) {
  man <- readRDS(file.path(REF, "run_manifest.rds")); bank <- readRDS(bank_file); cfgR <- man$config
  for (i in c(1L, 9L, 400L, 800L)) {
    bin <- read_binary_fingerprint(file.path(cfgR$binary_dir, man$files[i]), polarity = cfgR$polarity %||% "auto")
    bin_dsT <- downsample_matrix(bin, cfgR$target_size)
    a_now <- sample_anchor_points(bin_dsT, cfgR$v23$n_anchors, seed_from_string(man$files[i]))
    a_bank <- anchors_from_bank_entry(bank[[i]], cfgR)
    check(identical(unname(a_now), a_bank), sprintf("%s: anchors drawn now == anchors recovered from the canonical bank entry (identical())", man$files[i]))
    check(is.matrix(a_bank) && is.integer(a_bank) && identical(names(attributes(a_bank)), "dim"), sprintf("%s: recovered anchors are a plain integer matrix (no scaled:center attribute, no dimnames)", man$files[i]))
  }
} else cat("  (skipped: canonical run ", REF, " not present)\n", sep = "")

cat("\nT7 fresh extraction == canonical entry on shared fields; clamped field as specified\n")
if (length(bank_file) == 1 && dir.exists(cfgR$binary_dir)) {
  i <- 2L; p <- file.path(cfgR$binary_dir, man$files[i])
  e_new <- extract_feature_bank_one(p, cfgR, make_disk_mask(cfgR$ot$radius)); e_old <- bank[[i]]
  shared <- intersect(names(e_old), names(e_new))
  check(all(vapply(shared, function(fn) identical(e_old[[fn]], e_new[[fn]]), logical(1))), sprintf("%s: all %d shared fields identical to the canonical entry", man$files[i], length(shared)))
  bin <- read_binary_fingerprint(p, polarity = cfgR$polarity %||% "auto")
  r150 <- pad_or_crop_to_size(downsample_matrix(compute_dt_to_ridge(bin), cfgR$target_size), cfgR$target_size)
  v150 <- pad_or_crop_to_size(downsample_matrix(compute_dt_to_valley(bin), cfgR$target_size), cfgR$target_size)
  cr <- function(m) { x <- pmin(pmax(m, 0), cfgR$maxscale); x[!is.finite(x)] <- cfgR$maxscale; round(x, 1) }
  v <- c(as.numeric(cr(r150)), as.numeric(cr(v150))); v <- v / sqrt(mean(v^2))
  check(identical(v, e_new$dt_vec_150_clamped), "dt_vec_150_clamped == clamp(0,20) -> round(0.1) -> RMS-normalise of the 150x150 fields")
  check(!identical(e_new$dt_vec_150_clamped, e_new$dt_vec_150), "clamped vector differs from the unclamped one (fields above 20 exist)")
} else cat("  (skipped)\n")

cat("\nT8 verifier failure injection (symlink mirror of ", REF, ")\n", sep = "")
mirror_run <- function(src, dst, mutate_file = NULL, mutate = identity) {
  # mirror a run directory with symlinks; copy and mutate one fold cache if asked
  dir.create(file.path(dst, "cache", "folds"), recursive = TRUE, showWarnings = FALSE)
  for (f in list.files(src, recursive = FALSE, all.files = FALSE)) {
    if (f == "cache") next
    if (!dir.exists(file.path(src, f))) file.symlink(normalizePath(file.path(src, f)), file.path(dst, f))
  }
  for (f in list.files(file.path(src, "cache"), recursive = FALSE)) {
    if (f == "folds") next
    file.symlink(normalizePath(file.path(src, "cache", f)), file.path(dst, "cache", f))
  }
  for (f in list.files(file.path(src, "cache", "folds"))) {
    if (!is.null(mutate_file) && f == mutate_file) { obj <- readRDS(file.path(src, "cache", "folds", f)); saveRDS(mutate(obj), file.path(dst, "cache", "folds", f)) }
    else file.symlink(normalizePath(file.path(src, "cache", "folds", f)), file.path(dst, "cache", "folds", f))
  }
  invisible(dst)
}
run_verifier <- function(a, b) {
  old <- options(fp.verify_a = a, fp.verify_b = b, fp.verify_expect_field = "", fp.verify_expect_method = ""); on.exit(options(old))
  sink(tempfile()); on.exit(sink(), add = TRUE)                                      # the verifier prints ~500 lines; keep the test log short
  res <- tryCatch({ source("Code/verify_run_vs_canonical.R", local = new.env(parent = globalenv())); list(ok = TRUE, msg = "") },
                  error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  rep_file <- file.path(b, paste0("verification_vs_", basename(a), ".rds"))
  res$report <- if (file.exists(rep_file)) readRDS(rep_file) else NULL
  res
}
if (dir.exists(REF) && file.exists(file.path(REF, "fold_cache_index.csv")) && .Platform$OS.type == "unix") {
  idxR <- read.csv(file.path(REF, "fold_cache_index.csv"), stringsAsFactors = FALSE)
  target <- idxR$file[idxR$method == "GlobalBetti_L1" & idxR$fold == 2][1]
  m_clean <- mirror_run(REF, tempfile("mirror_clean_"))
  r_clean <- run_verifier(REF, m_clean)
  st_clean <- if (!is.null(r_clean$report)) table(factor(vapply(r_clean$report$results, `[[`, character(1), "status"), levels = c("EXACT", "TOLERANCE", "FAIL"))) else NULL
  check(r_clean$ok && !is.null(st_clean) && st_clean[["FAIL"]] == 0, sprintf("clean mirror passes: %s (EXACT %s, TOLERANCE %s, FAIL %s)", if (r_clean$ok) "no error" else r_clean$msg, st_clean[["EXACT"]] %||% NA, st_clean[["TOLERANCE"]] %||% NA, st_clean[["FAIL"]] %||% NA))
  m_bad <- mirror_run(REF, tempfile("mirror_perturbed_"), mutate_file = target, mutate = function(o) { o$s_test[1] <- o$s_test[1] + 1e-12; o })
  r_bad <- run_verifier(REF, m_bad)
  fails_bad <- if (!is.null(r_bad$report)) Filter(function(x) x$status == "FAIL", r_bad$report$results) else list()
  check(!r_bad$ok && grepl("FAILED", r_bad$msg, fixed = TRUE), paste("perturbed mirror (one score + 1e-12) is rejected with an error:", substr(r_bad$msg, 1, 90)))
  check(!is.null(r_bad$report) && length(fails_bad) >= 1 && any(grepl("fold 2 GlobalBetti_L1", vapply(fails_bad, `[[`, character(1), "label"), fixed = TRUE)) && any(grepl("s_test", vapply(fails_bad, `[[`, character(1), "label"), fixed = TRUE)),
        sprintf("report saved before the stop and names the perturbed object (%d FAIL row(s): %s)", length(fails_bad), paste(head(vapply(fails_bad, `[[`, character(1), "label"), 3), collapse = " | ")))
  check(length(fails_bad) >= 1 && !any(vapply(r_bad$report$results, function(x) x$status == "TOLERANCE" && grepl("s_test", x$label, fixed = TRUE), logical(1))), "the 1e-12 difference is a FAIL, not a TOLERANCE (core objects compared exactly)")
  unlink(c(m_clean, m_bad), recursive = TRUE)
} else cat("  (skipped: ", REF, " not present or no symlink support)\n", sep = "")

cat("\nT9 figure overlay coordinates == image() cell centres\n")
fig_src <- parse("Code/make_figures.R"); fig_env <- new.env()
for (e in fig_src) if (is.call(e) && identical(as.character(e[[1]]), "<-") && as.character(e[[2]]) %in% c("image_coords", "show_image", "pixel_xy", "image_cell_bounds", "check_pixel_coordinates")) eval(e, envir = fig_env)
check(all(vapply(c("image_coords", "show_image", "pixel_xy", "image_cell_bounds", "check_pixel_coordinates"), exists, logical(1), envir = fig_env, inherits = FALSE)), "helper definitions found in Code/make_figures.R")
ok9 <- is.null(err_of({ fig_env$check_pixel_coordinates(matrix(0, 300, 300)); fig_env$check_pixel_coordinates(matrix(0, 40, 60)); fig_env$check_pixel_coordinates(matrix(0, 7, 3)) }))
check(ok9, "pixel_xy() == image() cell centres for every cell (300x300, 40x60, 7x3)")
wrong_env <- new.env(parent = fig_env); wrong_env$pixel_xy <- function(i, j, M) list(x = (i - 1) / (nrow(M) - 1), y = 1 - (j - 1) / (ncol(M) - 1))
f_wrong <- fig_env$check_pixel_coordinates; environment(f_wrong) <- wrong_env      # same checker, wrong mapping
m9 <- err_of(f_wrong(matrix(0, 300, 300)))
check(!is.null(m9) && grepl("disagrees", m9, fixed = TRUE), paste("the image() default mapping (i-1)/(n-1) is caught:", substr(m9 %||% "", 1, 70)))
pdf(NULL); cc9 <- fig_env$show_image(matrix(0, 300, 300), col = c("black", "white"), axes = FALSE, asp = 1); u9 <- par("usr"); dev.off()
check(identical(cc9$x, (1:300 - 0.5) / 300) && u9[1] <= 0 && u9[2] >= 1 && u9[3] <= 0 && u9[4] >= 1, "show_image() passes centres (i-0.5)/n and a window containing [0,1]^2")

cat("\n", if (FAILS == 0L) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", FAILS), "\n", sep = "")
if (FAILS > 0L) stop(FAILS, " hardening check(s) failed", call. = FALSE)
