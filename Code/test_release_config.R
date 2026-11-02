############################################################
# RELEASE CONFIGURATION TEST (release 1.1.0, 1 Nov 2026). Run from the project root, production R, ≈ 10 s:
#   source("Code/test_release_config.R")
# Nothing is computed and nothing is written. The release sets CONFIG's defaults to the reported configuration
# (regenerated binaries, frozen alignment grid, no manuscript export). binary_dir enters bank_cache_key() and,
# through it, fold_cache_key(); the aligned settings enter the alignment key — so this test asserts, for both
# databases, that the EFFECTIVE configuration (after Code/run_canonical.R's overrides, exactly as the driver
# applies them) and every cache key equal those of the reported runs, and that the old pre-release defaults
# plus the same overrides give the same effective configuration. If either assertion failed, a release run
# would not reuse or reproduce the reported caches.
#   C1  DB1: effective configuration == Outputs_regen_v2's manifest configuration (scientific fields; output_dir
#       and the operational fields excluded; complex_type / polarity as effective values)
#   C2  DB1: bank key, fold key and aligned-cache key == the keys recorded in Outputs_regen_v2
#   C3  DB1: the release defaults ALONE (no driver override except output_dir) give the same keys
#   C4  DB3: C1 and C2 for Outputs_db3_v2
#   C5  old defaults (Revision/submitted_baseline/fingerprintTDA.R.pre_step5_1Nov, if present) + driver overrides
#       == new defaults + driver overrides, for both databases
############################################################

options(fp.suppress_main = TRUE)
source("Code/fingerprintTDA.R")
source("Code/frozen_grid.R")
FAILS <- 0L
check <- function(cond, label) { ok <- isTRUE(cond); cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", label)); if (!ok) FAILS <<- FAILS + 1L; invisible(ok) }
`%||%` <- function(a, b) if (is.null(a)) b else a
OPERATIONAL <- c("n_workers", "future_globals_gb")
norm_cfg <- function(cfg) { cfg$output_dir <- NULL; for (f in OPERATIONAL) cfg[[f]] <- NULL
  cfg$complex_type <- cfg$complex_type %||% "triangulated"; cfg$polarity <- cfg$polarity %||% "auto"; cfg[order(names(cfg))] }
driver_overrides <- function(cfg, binary_dir, output_dir) {        # exactly what Code/run_canonical.R does
  cfg$binary_dir <- binary_dir; cfg$output_dir <- output_dir; cfg$export_to_manuscript <- FALSE; apply_frozen_grid(cfg)
}
keys_of <- function(cfg) {
  paths <- get_paths(cfg)
  list(bank = substr(bank_cache_key(cfg), 1, 12), fold = substr(fold_cache_key(cfg), 1, 8), aligned = substr(aligned_cache_key(cfg$aligned, paths), 1, 12), n = length(paths))
}
recorded_keys <- function(run) {
  bank <- sub("^bank_([0-9a-f]{12})\\.rds$", "\\1", list.files(file.path(run, "cache"), pattern = "^bank_[0-9a-f]{12}\\.rds$"))
  al <- sub("^aligned_([0-9a-f]{12})\\.rds$", "\\1", list.files(file.path(run, "cache"), pattern = "^aligned_[0-9a-f]{12}\\.rds$"))
  idx <- read.csv(file.path(run, "fold_cache_index.csv"), stringsAsFactors = FALSE, colClasses = c(fold_key = "character"))   # an all-digit key must stay a string
  list(bank = bank, fold = unique(idx$fold_key), aligned = al)
}
SETS <- list(db1 = list(run = "Outputs_regen_v2", binary_dir = "Data/FVC2000_binary_regen"),
             db3 = list(run = "Outputs_db3_v2",   binary_dir = "Data/FVC2002_DB3_binary"))
NEW_DEFAULT <- CONFIG
for (db in names(SETS)) {
  S <- SETS[[db]]
  if (!file.exists(file.path(S$run, "run_manifest.rds"))) { cat("  (", db, ": ", S$run, " not present; skipped)\n", sep = ""); next }
  cat("\n", toupper(db), " against ", S$run, "\n", sep = "")
  man <- readRDS(file.path(S$run, "run_manifest.rds"))
  cfg <- driver_overrides(NEW_DEFAULT, S$binary_dir, S$run)
  a <- norm_cfg(man$config); b <- norm_cfg(cfg)
  diff <- c(names(b)[!vapply(names(b), function(n) identical(a[[n]], b[[n]]), logical(1))], setdiff(names(a), names(b)))
  check(identical(a, b), sprintf("effective configuration == the run's manifest configuration%s", if (length(diff)) paste(" (differs in:", paste(diff, collapse = ", "), ")") else ""))
  k <- keys_of(cfg); r <- recorded_keys(S$run)
  check(k$n == 800, "800 input images found")
  check(identical(k$bank, r$bank), sprintf("bank key %s == recorded %s", k$bank, paste(r$bank, collapse = ",")))
  check(length(r$fold) == 1 && identical(k$fold, r$fold), sprintf("fold key %s == recorded %s", k$fold, paste(r$fold, collapse = ",")))
  check(identical(k$aligned, r$aligned), sprintf("aligned-cache key %s == recorded %s", k$aligned, paste(r$aligned, collapse = ",")))
  if (db == "db1") {
    cfg0 <- NEW_DEFAULT; cfg0$output_dir <- S$run                 # the release defaults alone
    k0 <- keys_of(cfg0)
    check(identical(k0$bank, r$bank) && identical(k0$fold, r$fold) && identical(k0$aligned, r$aligned), "release defaults alone (no driver override) give the DB1 run's bank, fold and aligned keys")
    check(identical(norm_cfg(cfg0), a), "release defaults alone == the DB1 run's effective configuration")
  }
}

cat("\nC5 old defaults + driver overrides == new defaults + driver overrides\n")
OLD <- "Pattern Recognition/Revision/submitted_baseline/fingerprintTDA.R.pre_step5_1Nov"
if (file.exists(OLD)) {
  env_old <- new.env(parent = globalenv())
  # evaluate only the old CONFIG definition (same helpers; the old file's CONFIG expression needs aligned_default_config())
  exprs <- parse(OLD); old_cfg_expr <- NULL
  for (e in exprs) if (is.call(e) && identical(as.character(e[[1]]), "<-") && identical(as.character(e[[2]]), "CONFIG")) old_cfg_expr <- e
  stopifnot(!is.null(old_cfg_expr)); eval(old_cfg_expr, envir = env_old); OLD_DEFAULT <- env_old$CONFIG
  for (db in names(SETS)) {
    S <- SETS[[db]]
    o <- norm_cfg(driver_overrides(OLD_DEFAULT, S$binary_dir, S$run)); n <- norm_cfg(driver_overrides(NEW_DEFAULT, S$binary_dir, S$run))
    diff <- c(names(n)[!vapply(names(n), function(x) identical(o[[x]], n[[x]]), logical(1))], setdiff(names(o), names(n)))
    check(identical(o, n), sprintf("%s: identical effective configuration%s", db, if (length(diff)) paste(" (differs in:", paste(diff, collapse = ", "), ")") else ""))
    if (dir.exists(S$binary_dir)) {
      ko <- keys_of(driver_overrides(OLD_DEFAULT, S$binary_dir, S$run)); kn <- keys_of(driver_overrides(NEW_DEFAULT, S$binary_dir, S$run))
      check(identical(ko, kn), sprintf("%s: identical bank / fold / aligned keys", db))
    }
  }
  check(!identical(norm_cfg(OLD_DEFAULT), norm_cfg(NEW_DEFAULT)), "(the defaults themselves differ, as intended: binary_dir, aligned grid)")
} else cat("  (skipped: pre-release copy not present)\n")

cat("\n", if (FAILS == 0L) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", FAILS), "\n", sep = "")
if (FAILS > 0L) stop(FAILS, " release configuration check(s) failed", call. = FALSE)
