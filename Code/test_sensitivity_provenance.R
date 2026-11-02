############################################################
# SENSITIVITY PROVENANCE MUTATION TESTS (patch 1b, 31 Oct 2026). Run from the project root:
#   source("Code/test_sensitivity_provenance.R")
# Drives Code/sensitivity/score_perturbed.R in VALIDATE-ONLY mode (options(fp.sens_validate_only = TRUE)):
# provenance is validated and completion signatures are computed, nothing is scored, nothing is written
# outside tempdir(). Mutated inputs are COPIES under tempdir(); the project's calibration.json, manifests
# and generator scripts are never touched. Every check prints PASS/FAIL; the script stops at the end if any failed.
#
#   P1  the real inputs validate: one Tier R and one Tier B condition; the v2 signature carries the
#       provenance identity (calibration hash, chosen family/strength/level, degrade-manifest hash,
#       degraded-raw md5s, degrader hash for Tier R; overrides and enhancer hash for Tier B)
#   P2  a mutated calibration.json (the Tier R condition's strength changed) is refused
#   P3  a mutated degrade manifest (strength changed, calibration intact) is refused
#   P4  a mutated degraded-raw md5 in the degrade manifest (broken raw -> degraded -> binary chain) is refused
#   P5  a modified degrader script (Code/sensitivity/degrade_raw.py + one comment line) is refused
#   P7  a completed output directory is reused only on a matching signature: a fixture built from an existing
#       completed Tier R output (small files copied, fold caches as empty placeholders) with the current v2
#       signature injected is skipped as "completed_reused"; with one byte of the stored provenance changed,
#       or one stored condition-image md5 changed, it is refused ("DIFFERENT inputs"); without any signature
#       it is refused unless fp.sens_allow_unsigned_reuse = TRUE
#   P8  the path overrides are refused outside validate-only mode
############################################################

FAILS <- 0L
check <- function(cond, label) { ok <- isTRUE(cond); cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", label)); if (!ok) FAILS <<- FAILS + 1L; invisible(ok) }
`%||%` <- function(a, b) if (is.null(a)) b else a

COND_R <- "Data/FVC2000_binary_degraded/dry_typ"
COND_B <- "Data/FVC2000_binary_variant/filter_m2"
RAW_R  <- "Data/FVC2000_raw_degraded"
COMPLETED_R <- "Outputs_sensitivity_realistic/dry_typ"      # an existing completed Tier R output (fixture source for P7)
stopifnot(dir.exists(COND_R), dir.exists(RAW_R), file.exists(file.path(RAW_R, "calibration.json")))
have_B <- dir.exists(COND_B)

SENS_OPTS <- c("fp.sens_conditions", "fp.sens_out", "fp.sens_validate_only", "fp.sens_pairing", "fp.sens_tier_r_raw", "fp.sens_degrader_path", "fp.sens_allow_unsigned_reuse", "fp.sens_src")
run_scorer <- function(conds, out_root = tempfile("sens_out_"), validate_only = TRUE, ...) {
  clear_opts <- function() for (o in SENS_OPTS) { opt <- list(NULL); names(opt) <- o; options(opt) }
  old <- options()
  on.exit({ clear_opts(); options(old[SENS_OPTS[SENS_OPTS %in% names(old)]]) }, add = TRUE)   # leave the session's fp.sens_* options exactly as found
  clear_opts()                                                                              # start from a clean slate
  options(fp.sens_conditions = conds, fp.sens_out = out_root, fp.sens_validate_only = validate_only, fp.sens_pairing = "probe")
  extra <- list(...); if (length(extra)) options(extra)
  env <- new.env(parent = globalenv())
  sink(tempfile()); on.exit(sink(), add = TRUE)                                            # the scorer's own log is not needed here
  res <- tryCatch({ sys.source("Code/sensitivity/score_perturbed.R", envir = env); list(ok = TRUE, msg = "", result = env$VALIDATE_ONLY_RESULT) },
                  error = function(e) list(ok = FALSE, msg = conditionMessage(e), result = NULL))
  res
}
msg_of <- function(r) substr(gsub("\\s+", " ", r$msg), 1, 110)

cat("\nP1 real inputs validate; v2 signatures carry the provenance identity\n")
conds1 <- c(COND_R, if (have_B) COND_B)
r1 <- run_scorer(conds1)
check(r1$ok, paste("validate-only run on the real inputs succeeds:", if (r1$ok) "no error" else msg_of(r1)))
sigR <- r1$result[["dry_typ"]]$signature
check(!is.null(sigR) && identical(sigR$version, 2L) && identical(r1$result[["dry_typ"]]$status, "validated_not_scored"), "Tier R condition: status validated_not_scored, signature version 2")
cal_sha <- digest::digest(file = file.path(RAW_R, "calibration.json"), algo = "sha256")
dm <- jsonlite::fromJSON(file.path(RAW_R, "dry_typ", "degrade_manifest.json"))
check(!is.null(sigR$provenance) && identical(sigR$provenance$tier, "R") && identical(sigR$provenance$calibration_sha256, cal_sha) &&
        identical(sigR$provenance$degrader_sha256, digest::digest(file = "Code/sensitivity/degrade_raw.py", algo = "sha256")) &&
        identical(sigR$provenance$chosen$family, dm$family) && isTRUE(all.equal(sigR$provenance$chosen$strength, as.numeric(dm$strength))) &&
        length(sigR$provenance$degraded_raw_md5) == 800,
      "Tier R signature: calibration sha, degrader sha, chosen family/strength and 800 degraded-raw md5s present and correct")
if (have_B) check(identical(r1$result[["filter_m2"]]$signature$provenance$tier, "B") && !is.null(r1$result[["filter_m2"]]$signature$provenance$enhancer_source_sha256), "Tier B signature carries the enhancer hash and the overrides")

# ---- helpers building MUTATED COPIES of the Tier R raw-side manifests (only the two JSON files are needed for validation)
mutated_raw_dir <- function(mutate_cal = NULL, mutate_dm = NULL) {
  d <- tempfile("raw_mut_"); dir.create(file.path(d, "dry_typ"), recursive = TRUE)
  if (is.null(mutate_cal)) file.copy(file.path(RAW_R, "calibration.json"), file.path(d, "calibration.json")) else {
    cal <- jsonlite::fromJSON(file.path(RAW_R, "calibration.json"), simplifyVector = FALSE); cal <- mutate_cal(cal)
    jsonlite::write_json(cal, file.path(d, "calibration.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE) }
  if (is.null(mutate_dm)) file.copy(file.path(RAW_R, "dry_typ", "degrade_manifest.json"), file.path(d, "dry_typ", "degrade_manifest.json")) else {
    dmm <- jsonlite::fromJSON(file.path(RAW_R, "dry_typ", "degrade_manifest.json"), simplifyVector = FALSE); dmm <- mutate_dm(dmm)
    jsonlite::write_json(dmm, file.path(d, "dry_typ", "degrade_manifest.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE) }
  d
}
cat("\nP0 control: byte-identical copies of the two Tier R manifests under another root still validate\n")
r0 <- run_scorer(COND_R, fp.sens_tier_r_raw = mutated_raw_dir())
check(r0$ok, paste("copied manifests validate (so the refusals below are due to the mutations):", if (r0$ok) "no error" else msg_of(r0)))

cat("\nP2 mutated calibration.json is refused\n")
r2 <- run_scorer(COND_R, fp.sens_tier_r_raw = mutated_raw_dir(mutate_cal = function(cal) { cal$chosen$dry_typ$strength <- as.numeric(cal$chosen$dry_typ$strength) + 0.05; cal }))
check(!r2$ok, paste("refused:", msg_of(r2)))

cat("\nP3 mutated degrade manifest (strength) is refused\n")
r3 <- run_scorer(COND_R, fp.sens_tier_r_raw = mutated_raw_dir(mutate_dm = function(dmm) { dmm$strength <- as.numeric(dmm$strength) + 0.25; dmm }))
check(!r3$ok, paste("refused:", msg_of(r3)))

cat("\nP4 mutated degraded-raw md5 (broken chain) is refused\n")
r4 <- run_scorer(COND_R, fp.sens_tier_r_raw = mutated_raw_dir(mutate_dm = function(dmm) { dmm$images[[1]]$output_md5 <- "00000000000000000000000000000000"; dmm }))
check(!r4$ok, paste("refused:", msg_of(r4)))

cat("\nP5 modified degrader script is refused\n")
deg_copy <- tempfile("degrade_raw_", fileext = ".py"); writeLines(c(readLines("Code/sensitivity/degrade_raw.py"), "# mutated copy (test)"), deg_copy)
r5 <- run_scorer(COND_R, fp.sens_degrader_path = deg_copy)
check(!r5$ok, paste("refused:", msg_of(r5)))

cat("\nP7 completed output reused only on a matching signature\n")
if (dir.exists(COMPLETED_R) && file.exists(file.path(COMPLETED_R, "SENSITIVITY_COMPLETE.txt")) && !is.null(sigR)) {
  make_fixture <- function(edit_manifest = identity) {
    root <- tempfile("sens_fixture_"); out <- file.path(root, "dry_typ"); dir.create(file.path(out, "cache", "folds"), recursive = TRUE)
    for (f in c("fold_cache_index.csv", "SENSITIVITY_COMPLETE.txt", "fold_metrics_auc_eer.csv", "fold_metrics_tar.csv", "sensitivity_fold_metrics.csv")) file.copy(file.path(COMPLETED_R, f), file.path(out, f))
    idx <- read.csv(file.path(COMPLETED_R, "fold_cache_index.csv"), stringsAsFactors = FALSE)
    file.create(file.path(out, "cache", "folds", idx$file))                                     # placeholders: is_complete() checks existence only
    man <- readRDS(file.path(COMPLETED_R, "run_manifest.rds")); man$sensitivity$signature <- sigR; man <- edit_manifest(man)
    saveRDS(man, file.path(out, "run_manifest.rds")); root
  }
  r7a <- run_scorer(COND_R, out_root = make_fixture())
  check(r7a$ok && identical(r7a$result[["dry_typ"]]$status, "completed_reused"), paste("fixture with the current signature is reused:", if (r7a$ok) r7a$result[["dry_typ"]]$status %||% "(no status)" else msg_of(r7a)))
  r7b <- run_scorer(COND_R, out_root = make_fixture(function(man) { man$sensitivity$signature$provenance$calibration_sha256 <- "0000"; man }))
  check(!r7b$ok && grepl("DIFFERENT inputs", r7b$msg, fixed = TRUE), paste("stored provenance (calibration sha) altered -> refused:", msg_of(r7b)))
  r7c <- run_scorer(COND_R, out_root = make_fixture(function(man) { man$sensitivity$signature$condition_images$md5[1] <- "deadbeef"; man }))
  check(!r7c$ok && grepl("DIFFERENT inputs", r7c$msg, fixed = TRUE), paste("stored condition-image md5 altered -> refused:", msg_of(r7c)))
  r7d <- run_scorer(COND_R, out_root = make_fixture(function(man) { man$sensitivity$signature <- NULL; man }))
  check(!r7d$ok && grepl("before completion signatures existed", r7d$msg, fixed = TRUE), paste("unsigned completed output -> refused:", msg_of(r7d)))
  r7e <- run_scorer(COND_R, out_root = make_fixture(function(man) { man$sensitivity$signature <- NULL; man }), fp.sens_allow_unsigned_reuse = TRUE)
  check(r7e$ok && identical(r7e$result[["dry_typ"]]$status, "completed_reused"), "unsigned completed output reused only under the explicit escape option")
} else cat("  (skipped: ", COMPLETED_R, " not present or P1 failed)\n", sep = "")

cat("\nP8 path overrides refused outside validate-only mode\n")
r8 <- run_scorer(COND_R, validate_only = FALSE, fp.sens_degrader_path = deg_copy)
check(!r8$ok && grepl("fp.sens_validate_only", r8$msg, fixed = TRUE), paste("refused before anything is scored:", msg_of(r8)))

cat("\n", if (FAILS == 0L) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", FAILS), "\n", sep = "")
if (FAILS > 0L) stop(FAILS, " sensitivity provenance check(s) failed", call. = FALSE)
