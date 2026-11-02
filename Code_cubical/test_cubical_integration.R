############################################################
# STEP 1.5 (cubical side-by-side) — integration test of the complex_type switch added to
# Code/fingerprintTDA.R, run BEFORE the cubical run. Nothing is written anywhere.
#
# Run from the PROJECT ROOT in RStudio (~2-4 min):
#   source("Code_cubical/test_cubical_integration.R")
#
#   K1  cache keys of the existing runs are unchanged by the patch (bank key file names and
#       fold keys of the DB1 and DB3 runs (fp.test_run_db1 / fp.test_run_db3, defaults Outputs_regen_v2 / Outputs_db3_v2) reproduced from their manifests; explicit
#       complex_type = "triangulated" gives the same keys; "cubical" gives different ones)
#   T1  triangulated per-image features recomputed with the patched pipeline for 8 DB1 images
#       are bit-identical to the cached canonical bank entries (default path untouched)
#   S1  cubical backend self-test (hand-checkable cases incl. the T-construction check)
#   C1  cubical per-image features for the same 8 images: computed without error; the
#       persistence-independent fields (dt_vec, dt_vec_150, ot_Z anchors) are bit-identical to
#       the triangulated ones; the persistence-derived fields differ; H1 counts and bottleneck
#       distances of the global ridge diagrams are reported
#   P1  the same 8 cubical entries computed through precompute_feature_bank() with 4
#       multisession workers are bit-identical to the sequential ones (reticulate in workers)
############################################################

options(fp.suppress_main = TRUE)
source("Code/fingerprintTDA.R")
source("Code_cubical/cubical_ph.R"); PY <- configure_cubical_python()

PASS <- logical(0)
check <- function(ok, label, detail = "") {
  PASS <<- c(PASS, isTRUE(ok)); cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label))
  if (nzchar(detail)) cat("        -> ", substr(gsub("\\s+", " ", detail), 1, 240), "\n", sep = "")
}
section <- function(s) cat("\n", s, "\n", strrep("-", nchar(s)), "\n", sep = "")

# ---- K1: keys unchanged --------------------------------------------------------------------
section("K1: cache keys of the existing runs are unchanged")
RUN_DB1 <- getOption("fp.test_run_db1", "Outputs_regen_v2"); RUN_DB3 <- getOption("fp.test_run_db3", "Outputs_db3_v2")
for (d in c(RUN_DB1, RUN_DB3)) {
  m <- readRDS(file.path(d, "run_manifest.rds")); cfg <- m$config
  bank_files <- list.files(file.path(d, "cache"), pattern = "^bank_[0-9a-f]{12}\\.rds$")
  idx <- readr::read_csv(file.path(d, "fold_cache_index.csv"), show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  k_absent <- paste0("bank_", substr(bank_cache_key(cfg), 1, 12), ".rds")
  cfg_t <- cfg; cfg_t$complex_type <- "triangulated"
  cfg_c <- cfg; cfg_c$complex_type <- "cubical"
  check(identical(bank_files, k_absent), paste(d, ": bank key from the manifest config (complex_type absent) == existing bank file"), k_absent)
  check(identical(bank_cache_key(cfg_t), bank_cache_key(cfg)) && identical(fold_cache_key(cfg_t), fold_cache_key(cfg)),
        paste(d, ": explicit complex_type = 'triangulated' gives the same bank and fold keys"))
  check(all(idx$fold_key == substr(fold_cache_key(cfg), 1, 8)), paste(d, ": fold key == fold_cache_index.csv"), unique(idx$fold_key))
  check(!identical(bank_cache_key(cfg_c), bank_cache_key(cfg)) && !identical(fold_cache_key(cfg_c), fold_cache_key(cfg)),
        paste(d, ": complex_type = 'cubical' gives different bank and fold keys"), paste0("bank_", substr(bank_cache_key(cfg_c), 1, 12)))
}

# ---- T1: triangulated features unchanged -------------------------------------------------
section("T1: triangulated per-image features (8 DB1 images) == cached canonical bank entries")
m1   <- readRDS(file.path(RUN_DB1, "run_manifest.rds")); cfg1 <- m1$config
bank1 <- readRDS(file.path(RUN_DB1, "cache", paste0("bank_", substr(bank_cache_key(cfg1), 1, 12), ".rds")))
ii <- c(1, 2, 3, 100, 101, 400, 799, 800)                   # fixed set: first/last files and some in between
paths8 <- m1$paths[ii]; stopifnot(all(file.exists(paths8)))
disk <- make_disk_mask(cfg1$ot$radius)
# The production bank is computed in future::multisession workers, where the RNG kind is
# L'Ecuyer-CMRG (future.seed = TRUE); the per-image set.seed() inside the extractor therefore
# draws the anchors from that generator. A sequential session defaults to Mersenne-Twister and
# the same seed would select different anchors (QA, 13 Oct). Use the production generator for the
# sequential comparisons below and restore the session's RNG afterwards; production is untouched.
old_rng <- RNGkind(); RNGkind("L'Ecuyer-CMRG")
restore_rng <- function() do.call(RNGkind, as.list(old_rng))   # called at the end (on.exit has no effect at top level)
tri_seq <- lapply(paths8, function(fp) extract_feature_bank_one(fp, cfg1, disk))
check(all(vapply(seq_along(ii), function(j) identical(tri_seq[[j]], bank1[[ii[j]]]), logical(1))),
      "8 recomputed triangulated entries bit-identical to the cached bank (complex_type absent)")
cfg1t <- cfg1; cfg1t$complex_type <- "triangulated"
tri_seq2 <- lapply(paths8, function(fp) extract_feature_bank_one(fp, cfg1t, disk))
check(identical(tri_seq2, tri_seq), "same with explicit complex_type = 'triangulated'")

# ---- S1: backend self-test -----------------------------------------------------------------
section("S1: cubical backend self-test")
check(isTRUE(cubical_ph_selftest(verbose = TRUE)), "cubical_ph_selftest() passed")

# ---- C1: cubical features, sequential ------------------------------------------------------
section("C1: cubical per-image features for the same 8 images (sequential)")
cfg1c <- cfg1; cfg1c$complex_type <- "cubical"
t0 <- Sys.time()
cub_seq <- lapply(paths8, function(fp) extract_feature_bank_one(fp, cfg1c, disk))
cat(sprintf("  8 images in %.1f s (%.1f s / image, single process)\n", as.numeric(difftime(Sys.time(), t0, units = "secs")), as.numeric(difftime(Sys.time(), t0, units = "secs")) / 8))
check(all(!vapply(cub_seq, is.null, logical(1))), "all 8 cubical entries computed (no NULL)")
same_field <- function(f) all(vapply(1:8, function(j) identical(cub_seq[[j]][[f]], tri_seq[[j]][[f]]), logical(1)))
check(same_field("dt_vec") && same_field("dt_vec_150") && same_field("ot_Z"), "persistence-independent fields identical (dt_vec, dt_vec_150, ot_Z anchors)")
check(!same_field("betti_global_100") && !same_field("pi_global") && !same_field("v23_local_full_300") && !same_field("ot_X_ridge"),
      "persistence-derived fields differ (betti_global_100, pi_global, v23_local_full_300, ot_X_ridge)")
all_finite <- all(vapply(cub_seq, function(e) all(vapply(e[c("betti_global_100", "pi_global", "v23_local_full_300", "ot_X_ridge", "ot_X_valley")], function(v) all(is.finite(v)), logical(1))), logical(1)))
check(all_finite, "cubical persistence-derived fields all finite")
# diagram-level view of the global ridge filtration for the 8 images
cat("  global ridge diagram, H1 count triangulated / cubical, bottleneck distance:\n")
for (j in 1:8) {
  bin <- read_binary_fingerprint(paths8[j]); dtr <- downsample_matrix(compute_dt_to_ridge(bin), cfg1$target_size)
  dT <- compute_persistence(dtr, 1, cfg1$maxscale, "triangulated"); dC <- compute_persistence(dtr, 1, cfg1$maxscale, "cubical")
  hT <- as.matrix(dT)[as.matrix(dT)[, 1] == 1, 2:3, drop = FALSE]; hC <- dC[dC[, 1] == 1, 2:3, drop = FALSE]
  bn <- TDA::bottleneck(as.matrix(dT), dC, dimension = 1)
  cat(sprintf("    %-14s H1 %4d / %4d   bottleneck %.2f\n", basename(paths8[j]), nrow(hT), nrow(hC), bn))
}

# ---- E1: injected backend failures must STOP extraction, never become zero features ---------
section("E1: injected backend failures propagate (global and local persistence calls); a valid diagram without H1 is accepted")
cub_ok <- cubical_diagram
err_of <- function(expr) tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
cubical_diagram <- function(x, ...) if (nrow(x) >= 100) stop("injected GLOBAL failure") else cub_ok(x, ...)   # the 150 x 150 global call
e1 <- err_of(extract_feature_bank_one(paths8[1], cfg1c, disk))
check(!is.null(e1) && grepl("injected GLOBAL failure", e1, fixed = TRUE), "failure in the GLOBAL persistence call propagates out of the extractor", e1)
cubical_diagram <- function(x, ...) if (nrow(x) < 100) stop("injected LOCAL failure") else cub_ok(x, ...)     # the patch-level calls
e2 <- err_of(extract_feature_bank_one(paths8[1], cfg1c, disk))
check(!is.null(e2) && grepl("injected LOCAL failure", e2, fixed = TRUE), "failure in a LOCAL (patch) persistence call propagates out of the extractor", e2)
cfg_e <- cfg1c; cfg_e$n_workers <- 2L
e3 <- err_of(precompute_feature_bank(paths8[1:2], cfg_e)); future::plan(future::sequential)
check(!is.null(e3) && grepl("failed precompute", e3, fixed = TRUE) && grepl("injected LOCAL failure", e3, fixed = TRUE),
      "precompute_feature_bank() stops and names the image and the error", e3)
cubical_diagram <- function(x, ...) { d <- cub_ok(x, ...); d[d[, 1] == 0, , drop = FALSE] }                    # valid diagram, H0 only
r4 <- tryCatch(extract_feature_bank_one(paths8[1], cfg1c, disk), error = function(e) e)
check(!inherits(r4, "error") && all(r4$betti_global_100 == 0) && all(r4$pi_global == 0), "a valid diagram with no H1 features is accepted (zero persistence features, no error)")
cubical_diagram <- cub_ok
check(identical(extract_feature_bank_one(paths8[1], cfg1c, disk), cub_seq[[1]]), "backend restored: entry identical to C1 again")

# ---- P1: cubical in multisession workers ---------------------------------------------------
section("P1: the same 8 cubical entries through precompute_feature_bank() with 4 workers")
cfg1c$n_workers <- 4L
cub_par <- precompute_feature_bank(paths8, cfg1c)
future::plan(future::sequential)
check(identical(cub_par, cub_seq), "parallel cubical entries bit-identical to the sequential ones (reticulate works in workers)")

restore_rng()
check(identical(RNGkind(), old_rng), "session RNG kind restored")
cat(sprintf("\n%d / %d checks passed.\n", sum(PASS), length(PASS)))
if (!all(PASS)) stop("Cubical integration test FAILED - do not run the cubical pipeline until this is resolved.", call. = FALSE)
cat("Integration OK. Next (after QA clearance): source(\"Code_cubical/run_cubical_db1.R\")\n")
