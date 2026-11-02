############################################################
# STEP 1.1 (pipeline guards) — exercises the two stops added to fingerprintTDA.R
# (QA follow-up to step 1.2) WITHOUT computing anything:
#   G1  per-fold caches present, no run manifest           -> stop in run_prl_suite_v2()
#   G2  per-fold caches present, manifest matches, no bank -> stop in load_or_compute_bank()
#   G0  control: caches present, manifest matches, bank present -> bank LOADS (no stop)
#
# Run from the PROJECT ROOT in RStudio (about a minute):
#   source("Code/test_pipeline_guards.R")
#
# Creates the scratch directory Outputs_guard_test/, plants a dummy fold cache in it,
# copies the DB1 manifest / bank / sidecar from Outputs_regen (read only), and removes
# Outputs_guard_test/ at the end. Touches nothing else.
############################################################

options(fp.suppress_main = TRUE)
source("Code/fingerprintTDA.R")

ROOT <- "Outputs_guard_test"
SRC  <- getOption("fp.test_run", "Outputs_regen_v2")      # any complete DB1 run (read only)
if (dir.exists(ROOT)) stop(ROOT, " already exists; it is scratch from an earlier attempt - remove it and rerun.", call. = FALSE)
stopifnot(file.exists(file.path(SRC, "run_manifest.rds")), dir.exists("Data/FVC2000_binary_regen"))

plant <- function(d) {
  dir.create(file.path(d, "cache", "folds"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(note = "dummy fold cache planted by test_pipeline_guards.R"),
          file.path(d, "cache", "folds", "fold1_DT_L2_00000000.rds"))
}
PASS <- logical(0)
check <- function(ok, label, detail = "") {
  PASS <<- c(PASS, isTRUE(ok))
  cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label))
  if (nzchar(detail)) cat("        -> ", substr(gsub("\\s+", " ", detail), 1, 220), "\n", sep = "")
}
expect_stop <- function(expr, pattern, label) {
  msg <- tryCatch({ force(expr); "NO ERROR" }, error = function(e) conditionMessage(e))
  check(msg != "NO ERROR" && grepl(pattern, msg, fixed = TRUE), label, msg)
}
cfg <- CONFIG
cfg$binary_dir <- "Data/FVC2000_binary_regen"   # same images as the Outputs_regen manifest (needed for G2/G0)
cfg$export_to_manuscript <- FALSE

# ---- G1: fold caches without a run manifest ---------------------------------------------
d1 <- file.path(ROOT, "G1_no_manifest"); plant(d1)
c1 <- cfg; c1$output_dir <- d1
cat("\nG1: per-fold caches present, no run manifest\n")
expect_stop(run_prl_suite_v2(c1), "has no run manifest", "run_prl_suite_v2 stops")
check(!file.exists(file.path(d1, "run_manifest.rds")) && length(list.files(file.path(d1, "cache"), pattern = "^bank_")) == 0,
      "nothing written (no manifest, no bank)")

# ---- G2: fold caches with a matching manifest but no bank file ---------------------------
d2 <- file.path(ROOT, "G2_no_bank"); plant(d2)
file.copy(file.path(SRC, "run_manifest.rds"), file.path(d2, "run_manifest.rds"))
c2 <- cfg; c2$output_dir <- d2
cat("\nG2: per-fold caches present, manifest matches the images, no bank file\n")
expect_stop(run_prl_suite_v2(c2), "No cached feature bank for this configuration", "load_or_compute_bank stops")
check(length(list.files(file.path(d2, "cache"), pattern = "^bank_")) == 0, "no bank written")
future::plan(future::sequential)

# ---- G0 (control): caches + manifest + bank present -> the bank loads, no stop -------------
d0 <- file.path(ROOT, "G0_control"); plant(d0)
bank_src <- list.files(file.path(SRC, "cache"), pattern = "^bank_[0-9a-f]{12}\\.rds$", full.names = TRUE)
stopifnot(length(bank_src) == 1)
dir.create(file.path(d0, "cache"), showWarnings = FALSE, recursive = TRUE)
file.copy(bank_src, file.path(d0, "cache", basename(bank_src)))
file.copy(sub("bank_", "inputs_bank_", bank_src), file.path(d0, "cache", sub("bank_", "inputs_bank_", basename(bank_src))))
c0 <- cfg; c0$output_dir <- d0
cat("\nG0 (control): caches + matching manifest + bank present -> load_or_compute_bank must LOAD\n")
paths0 <- get_paths(c0)
key0 <- paste0("bank_", substr(bank_cache_key(c0), 1, 12), ".rds")
check(key0 == basename(bank_src), "bank key for this configuration equals the DB1 run's bank key", key0)
res <- tryCatch(load_or_compute_bank(paths0, c0), error = function(e) e)
check(!inherits(res, "error") && is.list(res) && length(res) == length(paths0), "bank loaded from cache without a stop",
      if (inherits(res, "error")) conditionMessage(res) else sprintf("%d images", length(res)))

# ---- clean up (only what this script created) ---------------------------------------------
unlink(ROOT, recursive = TRUE)
check(!dir.exists(ROOT), paste("scratch directory", ROOT, "removed"))

cat(sprintf("\n%d / %d checks passed.\n", sum(PASS), length(PASS)))
if (!all(PASS)) stop("Guard test FAILED - do not run DB3 until this is resolved.", call. = FALSE)
cat("Guards behave as specified. Next: source(\"Code/run_db3.R\")  (after QA clearance of STEP_1.1_DB3_SPEC.md)\n")
