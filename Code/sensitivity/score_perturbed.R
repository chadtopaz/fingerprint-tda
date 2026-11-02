############################################################
# FIXED-MODEL SENSITIVITY (R4-10): the clean-trained DB1 verifier applied to perturbed held-out images.
#
# For every condition directory (Tier B: Data/FVC2000_binary_variant/<cond>, the canonical binaries
# regenerated with one altered enhancement parameter; Tier R: Data/FVC2000_binary_degraded/<cond>, the raw
# images degraded at the calibrated strength and enhanced with the canonical procedure) the feature bank is
# computed for the perturbed images, and the canonical run's fitted objects (fp.sens_src, default
# Outputs_regen_v2, fold k) are applied, with NOTHING refitted or recalibrated, to the fold-k TEST
# pairs: TopoLR models, fusion weights and normalisations, and the training-calibrated thresholds
# are all the clean ones. The DB1 training scores stay the clean ones. The 15 variants of the pipeline are
# scored: DT_L2_150_aligned is not included (its rigid search would have to be redone per condition; it is
# geometry, not the representation under test), nor is the clamped geometry control DT_L2_150_clamped.
#
# Each condition is written as a run-shaped directory Outputs_sensitivity/<cond>/ so that
# Code/bootstrap_ci.R runs on it unchanged (same seed -> identical finger resamples as the clean
# run), and Code/sensitivity/run_sensitivity_bootstraps.R then computes perturbed - clean paired
# differences on identical pairs and resamples.
#
# Run from the PROJECT ROOT (~6 min per condition; the six Tier B conditions ~40 min):
#   source("Code/sensitivity/score_perturbed.R")
# The four Tier R conditions are driven by Code/sensitivity/run_realistic.R (probe pairing).
# Options: fp.sens_conditions = character vector of condition directories (default: the six Tier B conditions),
#          fp.sens_out = "Outputs_sensitivity", fp.sens_src = "Outputs_regen_v2", fp.sens_pairing = "both" | "probe".
# Reads the source run (never writes to it). Every requested condition must exist; a condition whose output
# directory exists but is incomplete STOPS the run (remove that directory to redo it); a completed one
# (completion marker + full index + matching completion signature) is skipped. Provenance of every condition is
# validated before anything else: Tier B raw-image md5s == the canonical DB1 raw images, the parameter
# overrides == the ones prescribed for that condition, and the enhancer source hash == the canonical one;
# Tier R: the binaries' preprocess manifest is the CANONICAL procedure (no overrides, zero errors, canonical
# enhancer hash) run on the degraded raw directory; the degraded raw directory's degrade manifest names the
# condition, family and strength recorded in calibration.json (file hash recorded), the degrader hash equals
# the project's degrade_raw.py, and every image chains canonical raw md5 -> degraded raw md5 -> binary md5.
# The binary reader runs with polarity = "ridge255" (white = ridge, never inverted), which coincides with the
# main pipeline's "auto" rule on every reported condition; the scorer asserts that "ridge255" reproduces the
# canonical bank entries on clean images before any condition is scored.
#
# Pairing (fp.sens_pairing):
#   "both"  (default; every image of every test pair is from the condition — Tier B)
#   "probe" (clean enrollment / degraded probe — Tier R: in every test pair one member's features come from the
#            canonical clean bank and the other's from the condition bank; WHICH member is the degraded probe is
#            a fixed seeded draw per pair, independent of impression number, finger number and file order,
#            shared by every condition and method and recorded in the manifests)
############################################################

options(fp.suppress_main = TRUE)
source("Code/fingerprintTDA.R")

SRC_DIR <- getOption("fp.sens_src", "Outputs_regen_v2")          # the clean-trained verifier (the reported DB1 run)
OUT_ROOT <- getOption("fp.sens_out", "Outputs_sensitivity")
PAIRING <- getOption("fp.sens_pairing", "both")
if (!PAIRING %in% c("both", "probe")) stop("fp.sens_pairing must be 'both' or 'probe'", call. = FALSE)
TIER_B <- file.path("Data/FVC2000_binary_variant", c("filter_m2", "filter_m4", "scale_055", "scale_075", "segment_020", "segment_040"))
TIER_R <- file.path("Data/FVC2000_binary_degraded", c("dry_typ", "wet_typ", "noise_typ", "blur_typ"))   # the calibrated ("typical") level, the only one
# Validate-only mode and path overrides (patch 1b) exist for Code/test_sensitivity_provenance.R: in validate-only
# mode every requested condition's provenance is validated and its completion signature computed (or, for a
# completed output, re-checked), nothing is scored and nothing is written; the overrides may point the
# degrader / degraded-raw manifests at mutated copies and are REFUSED outside validate-only mode.
VALIDATE_ONLY <- isTRUE(getOption("fp.sens_validate_only", FALSE))
TIER_R_RAW <- getOption("fp.sens_tier_r_raw", "Data/FVC2000_raw_degraded"); CALIBRATION_FILE <- file.path(TIER_R_RAW, "calibration.json")
DEGRADER_PATH <- getOption("fp.sens_degrader_path", "Code/sensitivity/degrade_raw.py")
if (!VALIDATE_ONLY && (TIER_R_RAW != "Data/FVC2000_raw_degraded" || DEGRADER_PATH != "Code/sensitivity/degrade_raw.py"))
  stop("fp.sens_tier_r_raw / fp.sens_degrader_path may only be set together with fp.sens_validate_only = TRUE (test mode); production scoring uses the project paths.", call. = FALSE)
TIER_R_NAMES <- basename(TIER_R)
DEGRADER_SHA <- digest::digest(file = DEGRADER_PATH, algo = "sha256")
COND_DIRS <- getOption("fp.sens_conditions", TIER_B)
TIER_B_OVERRIDES <- list(filter_m2 = list(ridge_filter_thresh = -2), filter_m4 = list(ridge_filter_thresh = -4),
                         scale_055 = list(relative_scale_factor_x = 0.55, relative_scale_factor_y = 0.55),
                         scale_075 = list(relative_scale_factor_x = 0.75, relative_scale_factor_y = 0.75),
                         segment_020 = list(ridge_segment_thresh = 0.2), segment_040 = list(ridge_segment_thresh = 0.4))
SCORER_SHA <- digest::digest(file = "Code/sensitivity/score_perturbed.R", algo = "sha256")      # this script's own identity, RECORDED in each manifest (not part of the reuse signature)
CANON_BIN_MANIFEST <- jsonlite::fromJSON("Data/FVC2000_binary_regen/preprocess_manifest.json")   # canonical binaries: raw md5s, enhancer hash
section <- function(s) cat("\n", s, "\n", strrep("-", nchar(s)), "\n", sep = "")
t_start <- Sys.time()

# ------------------------------------------------------------------ source run (frozen objects)
section("Source run (clean-trained verifier)")
src <- readRDS(file.path(SRC_DIR, "run_manifest.rds")); CFG <- src$config
fars <- CFG$fars_report; far_names <- paste0("far_", fars); n_folds <- CFG$n_folds
fkey_src <- substr(fold_cache_key(CFG), 1, 8)
idx_src <- readr::read_csv(file.path(SRC_DIR, "fold_cache_index.csv"), show_col_types = FALSE,
                           col_types = readr::cols(.default = readr::col_character(), fold = readr::col_integer(), loaded_from_cache = readr::col_logical()))
stopifnot(all(idx_src$fold_key == fkey_src))
load_src <- function(f, m) { hit <- idx_src[idx_src$fold == f & idx_src$method == m, ]; stopifnot(nrow(hit) == 1); readRDS(file.path(SRC_DIR, "cache", "folds", hit$file)) }
methods <- setdiff(unique(idx_src$method), c("DT_L2_150_aligned", "DT_L2_150_clamped"))
v23_variants <- CFG$v23_variants; ot_variants <- CFG$ot_variants
cat("  ", SRC_DIR, ": fold key ", fkey_src, "; ", length(methods), " methods scored per condition (DT_L2_150_aligned and DT_L2_150_clamped excluded)\n", sep = "")
cat("  output root: ", OUT_ROOT, "; pairing: '", PAIRING, "' (", if (PAIRING == "probe") "clean enrollment / degraded probe" else "both images of every test pair from the condition", ")\n", sep = "")
pipeline_thresholds <- function(s_train, y_train) { imp <- s_train[y_train == 0]
  setNames(vapply(fars, function(f) as.numeric(quantile(imp, probs = 1 - f, na.rm = TRUE)), numeric(1)), far_names) }
metrics_of <- function(y, s, thr) { mm <- roc_metrics(y, s)
  tibble(AUC = mm$auc, EER = mm$eer,
         !!!setNames(as.list(vapply(far_names, function(fn) mean(s[y == 1] >= thr[[fn]]), numeric(1))), paste0("TAR@", fars)),
         !!!setNames(as.list(vapply(far_names, function(fn) mean(s[y == 0] >= thr[[fn]]), numeric(1))), paste0("realisedFAR@", fars))) }
frozen_fields <- function(cc) cc[intersect(names(cc), c("model", "best_w", "zfit_ot", "zfit_v"))]
same_vec <- function(a, b) length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)) && all(a == b)

# the clean caches, loaded once
SRC <- lapply(seq_len(n_folds), function(f) setNames(lapply(methods, function(m) load_src(f, m)), methods))
for (f in seq_len(n_folds)) for (m in methods) {
  cc <- SRC[[f]][[m]]
  stopifnot(same_vec(unname(cc$thresholds[far_names]), unname(pipeline_thresholds(cc$s_train, cc$y_train))))   # stored thresholds == rule
}
if (!VALIDATE_ONLY) { future::plan(future::multisession, workers = CONFIG$n_workers); options(future.globals.maxSize = CFG$future_globals_gb * 1024^3) }

summary_rows <- list()
VALIDATE_ONLY_RESULT <- list()                                    # validate-only mode: per condition, status + signature (read by the test script)

# ---- probe-role assignment for pairing "probe" (QA, 14 Oct): which member of each ORIGINAL test pair is the
# degraded probe is a fixed, seeded draw per pair — not "the later-listed image" — so that the role is not
# tied to impression number, finger number or file order. The draw depends only on the fold and the seed,
# so it is shared by every condition and every method, and it is recorded (vectors + hash) in each
# condition's manifest; a completed condition is reused only if its recorded protocol matches.
PROBE_SEED_BASE <- 20261014L
probe_assignment <- function(f, n_pairs) {                       # isolated RNG state: Mersenne-Twister, seed = base + fold; session RNG restored
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  old_kind <- RNGkind(); on.exit({ do.call(RNGkind, as.list(old_kind)); if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv) }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(PROBE_SEED_BASE + f)
  sample(c(1L, 2L), n_pairs, replace = TRUE)                      # 1 = first member of the stored pair is the probe, 2 = second
}
PROBE_ASSIGN <- if (PAIRING == "probe") lapply(seq_len(n_folds), function(f) probe_assignment(f, ncol(src$folds[[f]]$pair_test))) else NULL
PROBE_SHA <- if (PAIRING == "probe") digest::digest(PROBE_ASSIGN, algo = "sha256") else NA_character_
if (PAIRING == "probe") {
  a1 <- probe_assignment(1L, ncol(src$folds[[1]]$pair_test)); stopifnot(identical(a1, PROBE_ASSIGN[[1]]))    # reproducible within the session
  cat(sprintf("  probe-role assignment: seeded per pair (Mersenne-Twister, seed %d + fold), shared by all conditions and methods; member 2 drawn in %.1f %% of test pairs; sha256 %s...\n",
              PROBE_SEED_BASE, 100 * mean(unlist(PROBE_ASSIGN) == 2L), substr(PROBE_SHA, 1, 12)))
}
protocol_matches <- function(out) {                                 # a completed condition is reused only under the same pairing protocol
  man <- readRDS(file.path(out, "run_manifest.rds")); sv <- man$sensitivity
  identical(sv$pairing %||% "both", PAIRING) && identical(sv$probe_assignment_sha %||% NA_character_, PROBE_SHA)
}
# Completion signature (audit patch, 31 Oct 2026; v2 in patch 1b): a completed condition is reused only if
# the inputs it was scored from are the CURRENT ones — the condition's manifest and its 800 binaries, the
# provenance identity validated by validate_provenance() (degrader / calibration.json hashes,
# chosen family-strength-level, degraded-raw md5s, enhancer hash), the source run's fold caches (the frozen
# models and thresholds) and bank key, and the pairing protocol (the scorer's own hash is recorded). Provenance is
# validated BEFORE the reuse shortcut, so a completed directory is never skipped on stale provenance.
# Before the patch the shortcut accepted marker + index + pairing only, so a regenerated condition set or
# a re-trained source under the same directory names would have been taken as current.
condition_signature <- function(cond_dir, tier, prov) {
  pm_file <- file.path(cond_dir, "preprocess_manifest.json")
  pngs <- sort(list.files(cond_dir, pattern = "_bin\\.png$", full.names = TRUE))
  src_files <- idx_src$file[idx_src$method %in% methods]
  list(version = 2L,
       condition_manifest_sha256 = digest::digest(file = pm_file, algo = "sha256"),
       condition_images = list(files = basename(pngs), md5 = unname(tools::md5sum(pngs))),
       provenance = prov,                                   # v2: degrader/calibration identity, chosen strengths, degraded-raw md5s, enhancer hash
       source_run = SRC_DIR, source_fold_key = fkey_src,
       source_fold_cache_md5 = setNames(unname(tools::md5sum(file.path(SRC_DIR, "cache", "folds", src_files))), src_files),
       source_bank_key = bank_cache_key(CFG),
       pairing = PAIRING, probe_assignment_sha = PROBE_SHA)
}
signature_matches <- function(out, cond_dir, tier, prov) {
  stored <- readRDS(file.path(out, "run_manifest.rds"))$sensitivity$signature
  if (is.null(stored)) {
    if (isTRUE(getOption("fp.sens_allow_unsigned_reuse", FALSE))) { cat("  WARNING: completed before completion signatures existed; reused on the pairing protocol only (fp.sens_allow_unsigned_reuse).\n"); return(TRUE) }
    stop("output directory ", out, " was completed before completion signatures existed, so its inputs cannot be verified against the current ones. ",
         "Score it into a fresh output root, or set options(fp.sens_allow_unsigned_reuse = TRUE) to reuse it on the pairing protocol alone.", call. = FALSE)
  }
  identical(stored, condition_signature(cond_dir, tier, prov))
}
is_complete <- function(out) {
  if (!file.exists(file.path(out, "SENSITIVITY_COMPLETE.txt")) || !file.exists(file.path(out, "fold_cache_index.csv"))) return(FALSE)
  idx <- readr::read_csv(file.path(out, "fold_cache_index.csv"), show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  nrow(idx) == length(methods) * n_folds && all(file.exists(file.path(out, "cache", "folds", idx$file))) &&
    all(file.exists(file.path(out, c("run_manifest.rds", "fold_metrics_auc_eer.csv", "fold_metrics_tar.csv", "sensitivity_fold_metrics.csv"))))
}
# every requested condition directory must exist before anything is scored
missing <- COND_DIRS[!dir.exists(COND_DIRS)]
if (length(missing)) stop("requested condition directories missing: ", paste(missing, collapse = ", "), call. = FALSE)

# ---- polarity = "ridge255" must reproduce the canonical bank on clean images (the canonical run used "auto"; no canonical image exceeds 50 % white)
if (!VALIDATE_ONLY) {
section("Polarity check: ridge255 == canonical on clean images")
cfg_clean <- CFG; cfg_clean$polarity <- "ridge255"; cfg_clean$n_workers <- CONFIG$n_workers; cfg_clean$output_dir <- tempfile("polarity_check_")
bank_canon <- readRDS(file.path(SRC_DIR, "cache", paste0("bank_", substr(bank_cache_key(CFG), 1, 12), ".rds")))
ii <- c(1, 2, 3, 100, 101, 400, 799, 800)
bank_chk <- precompute_feature_bank(src$paths[ii], cfg_clean)
same_on_shared <- function(a, b) { fn <- intersect(names(a), names(b)); length(fn) > 0 && all(vapply(fn, function(k) identical(a[[k]], b[[k]]), logical(1))) }
stopifnot(all(vapply(seq_along(ii), function(j) same_on_shared(bank_chk[[j]], bank_canon[[ii[j]]]), logical(1))))   # fields the canonical bank has must agree (newer fields may be absent there)
stopifnot(bank_cache_key(cfg_clean) != bank_cache_key(CFG))                                      # the rule is part of the cache identity
cat("  8 clean entries under polarity = 'ridge255' bit-identical to the canonical bank; key differs from the canonical key as intended.\n")
} else cat("\nVALIDATE-ONLY MODE: provenance and completion signatures only; nothing is scored or written.\n")

# ---- provenance (QA, 13/14 Oct; patch 1b: run BEFORE the reuse shortcut and folded into the completion signature):
# the condition's manifest must trace to the canonical inputs and to this experiment's prescription. Returns the
# tier, the manifest, the image key and a provenance identity (degrader / calibration hashes, chosen
# strengths, degraded-raw md5s, enhancer hash) that becomes part of the completion signature.
validate_provenance <- function(cond_dir, cond) {
  tier <- if (cond %in% names(TIER_B_OVERRIDES)) "B" else if (cond %in% TIER_R_NAMES) "R" else stop("unknown condition name: ", cond, call. = FALSE)
  pm <- jsonlite::fromJSON(file.path(cond_dir, "preprocess_manifest.json")); pm_rows <- pm$images
  if (tier == "R") {
    # (1) the binaries come from the CANONICAL procedure (library defaults, no overrides) applied to the degraded raw directory
    stopifnot(length(pm$parameter_overrides) == 0, grepl("defaults", pm$procedure, fixed = TRUE), !grepl("sensitivity variant", pm$procedure, fixed = TRUE),
              pm$n_errors == 0, identical(pm$environment$enhancer_source_sha256, CANON_BIN_MANIFEST$environment$enhancer_source_sha256),
              basename(pm$input_dir) == cond)
    # (2) the degraded raw directory: condition, family and strength as recorded in calibration.json; generator hash; canonical raw inputs
    cal <- jsonlite::fromJSON(CALIBRATION_FILE, simplifyVector = FALSE)
    cal_sha <- digest::digest(file = CALIBRATION_FILE, algo = "sha256")
    stopifnot(setequal(names(cal$chosen), TIER_R_NAMES), basename(cal$benchmark$canonical_dir) == "FVC2000_binary_regen")
    dm <- jsonlite::fromJSON(file.path(TIER_R_RAW, cond, "degrade_manifest.json")); dm_rows <- dm$images
    ch <- cal$chosen[[cond]]
    stopifnot(identical(dm$condition, cond), identical(dm$family, ch$family), isTRUE(all.equal(as.numeric(dm$strength), as.numeric(ch$strength))),
              identical(dm$level, ch$level), identical(dm$calibration_sha256, cal_sha), identical(dm$script_sha256, DEGRADER_SHA),
              basename(dm$input_dir) == "FVC2000_raw", nrow(dm_rows) == 800)
    canon_raw <- setNames(CANON_BIN_MANIFEST$images$input_md5, sub("_bin\\.png$", ".tif", CANON_BIN_MANIFEST$images$output))
    stopifnot(setequal(dm_rows$file, names(canon_raw)), identical(dm_rows$input_md5, unname(canon_raw[dm_rows$file])))   # degraded FROM the canonical raw images
    # (3) chain: the binaries' recorded inputs are the degraded raw files
    deg_md5 <- setNames(dm_rows$output_md5, dm_rows$file)
    stopifnot(setequal(pm_rows$input, dm_rows$file), identical(pm_rows$input_md5, unname(deg_md5[pm_rows$input])))
    key <- pm_rows$output
    stopifnot(setequal(key, src$files))
    cat(sprintf("  Tier R: %s, strength %s (%s level); calibration.json sha256 %s...; chain canonical raw -> degraded raw -> binary verified for 800 images\n",
                ch$family, format(ch$strength), ch$level, substr(cal_sha, 1, 12)))
    prov <- list(tier = "R", calibration_sha256 = cal_sha, chosen = list(family = ch$family, strength = as.numeric(ch$strength), level = ch$level),
                 degrade_manifest_sha256 = digest::digest(file = file.path(TIER_R_RAW, cond, "degrade_manifest.json"), algo = "sha256"),
                 degraded_raw_md5 = setNames(dm_rows$output_md5, dm_rows$file)[order(dm_rows$file)],
                 degrader_sha256 = DEGRADER_SHA, enhancer_source_sha256 = pm$environment$enhancer_source_sha256)
  } else {
    exp <- TIER_B_OVERRIDES[[cond]]
    ov <- pm$parameter_overrides
    stopifnot(setequal(names(ov), names(exp)), all(vapply(names(exp), function(k) isTRUE(all.equal(as.numeric(ov[[k]]), as.numeric(exp[[k]]))), logical(1))),
              grepl("sensitivity variant", pm$procedure, fixed = TRUE), pm$n_errors == 0,
              identical(pm$environment$enhancer_source_sha256, CANON_BIN_MANIFEST$environment$enhancer_source_sha256),
              basename(pm$input_dir) == "FVC2000_raw")
    key <- pm_rows$output
    canon_raw <- setNames(CANON_BIN_MANIFEST$images$input_md5, CANON_BIN_MANIFEST$images$output)
    stopifnot(setequal(key, src$files), identical(pm_rows$input_md5, unname(canon_raw[key])))     # regenerated FROM the canonical raw images
    prov <- list(tier = "B", parameter_overrides = ov, enhancer_source_sha256 = pm$environment$enhancer_source_sha256)
  }
  list(tier = tier, pm = pm, pm_rows = pm_rows, key = key, prov = prov)
}

for (cond_dir in COND_DIRS) {
  cond <- basename(cond_dir); out <- file.path(OUT_ROOT, cond)
  section(paste("Condition", cond))
  pv <- validate_provenance(cond_dir, cond)                      # provenance FIRST (patch 1b): also when the condition is already completed
  tier <- pv$tier; pm <- pv$pm; pm_rows <- pv$pm_rows; key <- pv$key
  if (dir.exists(out) && length(list.files(out, recursive = TRUE, all.files = TRUE, no.. = TRUE)) > 0) {
    if (is_complete(out)) {
      if (!protocol_matches(out)) stop("output directory ", out, " was completed under a DIFFERENT pairing protocol (pairing / probe-role assignment). Remove it to redo this condition under the current protocol.", call. = FALSE)
      if (!signature_matches(out, cond_dir, tier, pv$prov)) stop("output directory ", out, " was completed from DIFFERENT inputs (condition images/manifest, provenance identity, source fold caches or bank key differ from the current ones). ",
                                                                "Use a fresh output root for the current inputs.", call. = FALSE)
      cat("  completed earlier (marker + full index, same pairing protocol, provenance re-validated, same input signature), skipped:", out, "\n")
      if (VALIDATE_ONLY) VALIDATE_ONLY_RESULT[[cond]] <- list(status = "completed_reused", tier = tier, signature = condition_signature(cond_dir, tier, pv$prov))
      next
    }
    stop("output directory ", out, " exists but is INCOMPLETE (interrupted run?). Remove it to redo this condition.", call. = FALSE)
  }
  if (VALIDATE_ONLY) {
    VALIDATE_ONLY_RESULT[[cond]] <- list(status = "validated_not_scored", tier = tier, signature = condition_signature(cond_dir, tier, pv$prov))
    cat("  validate-only: provenance OK (tier ", tier, "); completion signature v2 computed; not scored.\n", sep = ""); next
  }
  cfg_t <- CFG; cfg_t$binary_dir <- cond_dir; cfg_t$output_dir <- out; cfg_t$export_to_manuscript <- FALSE; cfg_t$n_workers <- CONFIG$n_workers
  cfg_t$polarity <- "ridge255"
  paths_t <- get_paths(cfg_t)
  stopifnot(length(paths_t) == 800, identical(basename(paths_t), src$files))                    # same names and order as the clean manifest
  md5_t <- unname(tools::md5sum(paths_t)); stopifnot(!identical(md5_t, src$md5))                # ... but different images
  ids_t <- map_chr(paths_t, get_identity); stopifnot(identical(ids_t, src$ids))
  stopifnot(identical(pm_rows$output_md5[match(basename(paths_t), key)], md5_t))                 # files == the condition's own manifest
  # pixel agreement with the canonical binaries, per image (reported with the results)
  agree <- vapply(seq_along(paths_t), function(i) {
    a <- read_binary_fingerprint(src$paths[i], polarity = "ridge255"); b <- read_binary_fingerprint(paths_t[i], polarity = "ridge255"); mean(a == b) }, numeric(1))
  n_white <- sum(vapply(paths_t, function(p) mean(read_binary_fingerprint(p, polarity = "ridge255")) > 0.5, logical(1)))
  cat(sprintf("  provenance OK (tier %s); pixel agreement with canonical: median %.4f, 10th pct %.4f, min %.4f; images > 50 %% white: %d (polarity fixed, no inversion)\n",
              tier, median(agree), quantile(agree, 0.1), min(agree), n_white))
  dir.create(file.path(out, "cache", "folds"), recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  bank_t <- load_or_compute_bank(paths_t, cfg_t)
  if (PAIRING == "probe") {
    # clean enrollment / degraded probe: entries 1..800 = canonical clean bank, 801..1600 = condition bank;
    # in each test pair the member drawn as the probe (PROBE_ASSIGN, seeded per pair, shared by all conditions
    # and methods) is offset by 800 and so read from the condition bank; the other member stays clean. The pair
    # order is preserved; the labels come from the original indices.
    stopifnot(length(bank_canon) == 800, length(bank_t) == 800)
    bank_s <- c(bank_canon, bank_t); offset <- 800L
  } else { bank_s <- bank_t; offset <- 0L }
  dt_t <- bank_to_matrix(bank_s, "dt_vec"); dt150_t <- bank_to_matrix(bank_s, "dt_vec_150")
  betti_t <- bank_to_matrix(bank_s, "betti_global_100"); gpi_t <- bank_to_matrix(bank_s, "pi_global")
  v23_t <- setNames(lapply(v23_variants, function(vv) get_v23_feature_matrix(bank_s, vv)), v23_variants)
  index_rows <- list(); auc_rows <- list(); tar_rows <- list(); fm_rows <- list()
  for (f in seq_len(n_folds)) {
    mf <- src$folds[[f]]; pair_test <- mf$pair_test
    y_test <- compute_y_from_pair_mat(pair_test, src$ids)                                        # labels from the ORIGINAL pair indices
    pair_s <- pair_test
    if (PAIRING == "probe") {
      a <- PROBE_ASSIGN[[f]]; stopifnot(length(a) == ncol(pair_test), all(a %in% c(1L, 2L)))
      pair_s[1, a == 1L] <- pair_s[1, a == 1L] + offset                                           # drawn member 1 is the degraded probe
      pair_s[2, a == 2L] <- pair_s[2, a == 2L] + offset                                           # drawn member 2 is the degraded probe
      stopifnot(all((pair_s[1, ] > 800L) + (pair_s[2, ] > 800L) == 1L))                           # exactly one degraded image per pair
    }
    s <- list()
    s[["DT_L2"]] <- score_pairs_l2(pair_s, dt_t); s[["DT_L2_150"]] <- score_pairs_l2(pair_s, dt150_t)
    s[["GlobalBetti_L1"]] <- score_pairs_l1(pair_s, betti_t); s[["GlobalPI_L2"]] <- score_pairs_l2(pair_s, gpi_t)
    for (vv in v23_variants) s[[vv]] <- predict_logistic_pairs(pair_s, v23_t[[vv]], SRC[[f]][[vv]]$model, chunk_size = 8000)
    for (ii in seq_len(nrow(ot_variants))) { vrow <- ot_variants[ii, ]; s[[vrow$method]] <- score_pairs_ot(pair_s, bank_s, cfg_t, vrow) }
    fu <- SRC[[f]][["Fusion"]]; fe <- SRC[[f]][["Fusion_eq"]]; stopifnot(identical(as.numeric(fe$best_w), 0.5))
    s[["Fusion"]]    <- apply_fusion(s[["OT"]], s[["V23_full"]], fu$best_w, fu$zfit_ot, fu$zfit_v)
    s[["Fusion_eq"]] <- apply_fusion(s[["OT"]], s[["V23_full"]], 0.5,       fe$zfit_ot, fe$zfit_v)
    if (length(setdiff(methods, names(s)))) stop("no sensitivity scorer for method(s): ", paste(setdiff(methods, names(s)), collapse = ", "), call. = FALSE)
    for (m in methods) {
      cc <- SRC[[f]][[m]]; stopifnot(identical(as.integer(cc$y_test), y_test), length(cc$s_test) == length(s[[m]]))
      thr <- cc$thresholds[far_names]
      obj <- c(list(s_train = cc$s_train, s_test = s[[m]], y_train = cc$y_train, y_test = y_test, elapsed = NA_real_, thresholds = thr), frozen_fields(cc),
               list(sensitivity = list(condition = cond, condition_dir = cond_dir, source_run = SRC_DIR, source_fold = f,
                                       note = "clean DB1 training scores and thresholds; test scores on perturbed images under the frozen objects")))
      fn <- sprintf("sens_fold%d_%s_%s.rds", f, m, fkey_src); saveRDS(obj, file.path(out, "cache", "folds", fn))
      index_rows[[length(index_rows) + 1]] <- tibble(fold = f, method = m, file = fn, fold_key = fkey_src, extra_key = NA_character_, loaded_from_cache = FALSE)
      met <- metrics_of(y_test, s[[m]], thr); clean <- metrics_of(y_test, cc$s_test, thr)
      auc_rows[[length(auc_rows) + 1]] <- tibble(fold = f, method = m, AUC = met$AUC, EER = met$EER)
      tar_rows[[length(tar_rows) + 1]] <- tibble(fold = f, method = m, FAR = fars, TAR = as.numeric(met[paste0("TAR@", fars)]))
      fm_rows[[length(fm_rows) + 1]] <- bind_rows(mutate(clean, input = "clean"), mutate(met, input = "perturbed")) %>%
        mutate(condition = cond, fold = f, method = m, score_corr = cor(cc$s_test, s[[m]]), .before = 1)
    }
    cat(sprintf("  fold %d scored (%.1f min)\n", f, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  readr::write_csv(bind_rows(index_rows), file.path(out, "fold_cache_index.csv"))
  readr::write_csv(bind_rows(auc_rows), file.path(out, "fold_metrics_auc_eer.csv"))
  readr::write_csv(bind_rows(tar_rows), file.path(out, "fold_metrics_tar.csv"))
  fm <- bind_rows(fm_rows); readr::write_csv(fm, file.path(out, "sensitivity_fold_metrics.csv"))
  man <- src; man$created <- Sys.time()
  # (the manifest keeps the clean config so that bootstrap_ci.R's fold-key check matches the copied index;
  #  the polarity rule used for the perturbed bank is recorded below and in the condition's own bank key)
  man$sensitivity <- list(condition = cond, tier = tier, condition_dir = cond_dir, condition_manifest = pm[setdiff(names(pm), "images")], source_run = SRC_DIR,
                          methods = methods, polarity = "ridge255", n_images_over_half_white = n_white,
                          pixel_agreement_with_canonical = list(per_image = agree, median = median(agree), q10 = unname(quantile(agree, 0.1)), min = min(agree)),
                          pairing = PAIRING, both_images_of_every_pair_perturbed = (PAIRING == "both"), single_realisation = TRUE,
                          probe_assignment = PROBE_ASSIGN, probe_assignment_sha = PROBE_SHA,
                          signature = condition_signature(cond_dir, tier, pv$prov), scorer_sha256 = SCORER_SHA,
                          probe_protocol = if (PAIRING == "probe") sprintf("per original test pair, the degraded member is drawn by sample(c(1, 2)) under Mersenne-Twister with seed %d + fold (1 = first stored member, 2 = second); the same draw is used by every condition and method; pair order preserved; labels from the original indices", PROBE_SEED_BASE) else NA_character_,
                          spec = "Pattern Recognition/Revision/STEP_1.4_SENSITIVITY_SPEC.md", created = Sys.time())
  saveRDS(man, file.path(out, "run_manifest.rds"))
  writeLines(c("SENSITIVITY DIRECTORY - NOT A PIPELINE RUN", sprintf("Condition %s (%s, tier %s, pairing '%s'); written by Code/sensitivity/score_perturbed.R on %s.", cond, cond_dir, tier, PAIRING, format(Sys.time())),
               "Caches: clean DB1 fold-k training scores and thresholds; test scores on the perturbed images under the frozen fold-k objects.",
               if (PAIRING == "probe") sprintf("Pairing 'probe': in every test pair one image is the canonical clean one and the other the condition's; the degraded member is a seeded per-pair draw (seed %d + fold), identical across conditions and methods, recorded in run_manifest.rds$sensitivity$probe_assignment (sha256 %s).", PROBE_SEED_BASE, PROBE_SHA) else "Pairing 'both': both images of every test pair are the condition's.",
               "Layout exists so that Code/bootstrap_ci.R can analyse it unchanged; DT_L2_150_aligned not included."), file.path(out, "SENSITIVITY_README.txt"))
  writeLines(c(sprintf("condition %s completed %s", cond, format(Sys.time())), sprintf("pixel agreement with canonical: median %.4f q10 %.4f min %.4f", median(agree), quantile(agree, 0.1), min(agree))),
             file.path(out, "SENSITIVITY_COMPLETE.txt"))
  sm <- fm %>% filter(input == "perturbed") %>% group_by(method) %>% summarize(across(c(AUC, EER, `TAR@0.001`, `realisedFAR@0.001`, `TAR@0.01`, `realisedFAR@0.01`, score_corr), mean), .groups = "drop") %>%
    mutate(condition = cond, tier = tier, agreement_median = median(agree), agreement_q10 = unname(quantile(agree, 0.1)), agreement_min = min(agree), .before = 1)
  summary_rows[[length(summary_rows) + 1]] <- sm
  show <- sm %>% filter(method %in% c("GlobalBetti_L1", "V23_full", "OT", "Fusion_eq")) %>% mutate(across(where(is.numeric), ~ round(.x, 4)))
  print(as.data.frame(show), row.names = FALSE)
  cat(sprintf("  condition %s done in %.1f min\n", cond, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
future::plan(future::sequential)
if (length(summary_rows)) {
  dir.create(file.path(OUT_ROOT, "summary"), showWarnings = FALSE, recursive = TRUE)
  prev <- file.path(OUT_ROOT, "summary", "point_estimates_by_condition.csv")
  all_sm <- bind_rows(summary_rows); if (file.exists(prev)) all_sm <- bind_rows(readr::read_csv(prev, show_col_types = FALSE) %>% filter(!condition %in% all_sm$condition), all_sm)
  readr::write_csv(all_sm, prev)
}
if (VALIDATE_ONLY) {
  st <- vapply(VALIDATE_ONLY_RESULT, `[[`, character(1), "status")
  cat(sprintf("\nVALIDATE-ONLY: %d condition(s) validated (%d completed outputs re-checked and reusable, %d not scored). Nothing written.\n",
              length(st), sum(st == "completed_reused"), sum(st == "validated_not_scored")))
} else {
  done <- vapply(file.path(OUT_ROOT, basename(COND_DIRS)), is_complete, logical(1))
  if (!all(done)) stop("conditions NOT complete: ", paste(basename(COND_DIRS)[!done], collapse = ", "), call. = FALSE)
  cat(sprintf("\nAll %d requested conditions complete. Done in %.1f minutes. Next:  source(\"Code/sensitivity/run_sensitivity_bootstraps.R\")\n", length(COND_DIRS), as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}
