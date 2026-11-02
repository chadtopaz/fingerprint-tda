############################################################
# STEP 1.1 (C1, transfer) — DB1 -> DB3 TRANSFER TEST, pre-specified in
# Pattern Recognition/Revision/STEP_1.1_DB3_SPEC.md §5 (fixed 13 Oct 2026 before any DB3 result).
#
# The five fitted DB1 models of the canonical run (Outputs_regen, fold k = 1..5) are applied,
# with NOTHING refitted, retuned or re-standardised, to the DB3 fold-k TEST pairs of the
# replication run (Outputs_db3). Fold-matched: DB1 fold k -> DB3 fold k (the two runs share the
# fold assignment of finger NUMBERS by construction of the shared seed; different people).
#
# Two calibrations, each written as a run-shaped directory so that Code/bootstrap_ci.R runs on it
# UNCHANGED via options(fp.bootstrap_dir = ...):
#   Outputs_db3_transfer/        PRIMARY   : DB1 fold-k thresholds (TAR is read with the realised FAR)
#   Outputs_db3_transfer_recal/  SECONDARY : thresholds re-estimated on the DB3 fold-k TRAINING impostor
#                                            scores under the frozen DB1 model (pipeline rule, type 7)
#
# Per method (all from the DB1 fold-k cache):
#   training-free (DT_L2, DT_L2_150, DT_L2_150_aligned, GlobalBetti_L1, GlobalPI_L2, OT + 4 ablations):
#       scores recomputed from the DB3 bank / DB3 aligned cache and ASSERTED bit-exact equal to the
#       replication's test scores (structural identity); only the thresholds transfer
#   V23_* (TopoLR): predict_logistic_pairs() with the stored model (w, b, standardizer, ops)
#   Fusion: apply_fusion() with stored best_w, zfit_ot, zfit_v;  Fusion_eq: w = 0.5, its own stored zfits
#
# Run from the PROJECT ROOT in RStudio (~10-15 min), after the DB3 replication + its bootstrap:
#   source("Code/transfer_db1_to_db3.R")
# Then:
#   options(fp.bootstrap_dir = "Outputs_db3_transfer");       source("Code/bootstrap_ci.R")
#   options(fp.bootstrap_dir = "Outputs_db3_transfer_recal"); source("Code/bootstrap_ci.R")
#   source("Code/transfer_vs_replication.R")
#
# Reads Outputs_regen/ and Outputs_db3/ and never writes to them (the DB3 bank is loaded through
# its sidecar, the aligned cache through its key; both are asserted untouched afterwards).
# Refuses to run if either output directory already contains files.
############################################################

options(fp.suppress_main = TRUE)
source("Code/fingerprintTDA.R")
source("Code/frozen_grid.R")

SRC_DIR    <- getOption("fp.transfer_src", "Outputs_regen_v2")   # DB1 canonical run: fitted objects, thresholds, training scores
TGT_DIR    <- getOption("fp.transfer_tgt", "Outputs_db3_v2")     # DB3 replication run: manifest (pairs), bank, aligned cache, replication caches
OUT_FROZEN <- "Outputs_db3_transfer"
OUT_RECAL  <- "Outputs_db3_transfer_recal"
section <- function(s) cat("\n", s, "\n", strrep("-", nchar(s)), "\n", sep = "")
t_start <- Sys.time()

for (d in c(OUT_FROZEN, OUT_RECAL))
  if (dir.exists(d) && length(list.files(d, recursive = TRUE, all.files = TRUE, no.. = TRUE)) > 0)
    stop(d, " already contains files; the transfer must be written into new or empty directories.", call. = FALSE)

# ------------------------------------------------------------------ manifests and protocol identity
section("Source (DB1) and target (DB3) runs")
src <- readRDS(file.path(SRC_DIR, "run_manifest.rds"))
tgt <- readRDS(file.path(TGT_DIR, "run_manifest.rds"))
cfg_a <- src$config; cfg_b <- tgt$config
for (fld in c("binary_dir", "output_dir", "export_to_manuscript")) { cfg_a[[fld]] <- NULL; cfg_b[[fld]] <- NULL }
stopifnot(identical(cfg_a, cfg_b))                                   # same protocol and parameters
stopifnot(identical(src$files, tgt$files), identical(src$ids, tgt$ids), !identical(src$md5, tgt$md5))
n_folds <- tgt$config$n_folds
for (k in seq_len(n_folds)) {
  a <- src$folds[[k]]; b <- tgt$folds[[k]]
  stopifnot(identical(a$test_ids, b$test_ids), identical(a$train_ids, b$train_ids),
            identical(a$train_idx, b$train_idx), identical(a$test_idx, b$test_idx),
            identical(a$pair_train, b$pair_train), identical(a$pair_test, b$pair_test))
}
fars <- tgt$config$fars_report; far_names <- paste0("far_", fars)
fkey_tgt <- substr(fold_cache_key(tgt$config), 1, 8)
cat("  protocol identical; file names, identities, fold memberships and pair orderings identical; images differ (md5).\n")
cat("  DB3 fold key:", fkey_tgt, "| FARs:", paste(fars, collapse = ", "), "\n")

read_index <- function(d) {
  idx <- readr::read_csv(file.path(d, "fold_cache_index.csv"), show_col_types = FALSE,
                         col_types = readr::cols(.default = readr::col_character(), fold = readr::col_integer(), loaded_from_cache = readr::col_logical()))
  if (any(duplicated(idx[, c("fold", "method")]))) stop("duplicate (fold, method) in ", d)
  idx
}
idx_src <- read_index(SRC_DIR); idx_tgt <- read_index(TGT_DIR)
stopifnot(all(idx_tgt$fold_key == fkey_tgt))
load_cache <- function(d, idx, f, m) {
  hit <- idx[idx$fold == f & idx$method == m, ]
  if (nrow(hit) != 1) stop("no unique cache for fold ", f, " method ", m, " in ", d)
  cc <- readRDS(file.path(d, "cache", "folds", hit$file))
  if (!is.na(hit$extra_key) && !identical(cc$aligned_key, hit$extra_key)) stop("extra key mismatch: ", hit$file)
  cc
}
v23_variants <- tgt$config$v23_variants
ot_variants  <- tgt$config$ot_variants
source("Code/transfer_scope.R")                                     # the declared 16-method scope; the clamped control is outside it
scope <- transfer_method_scope(unique(idx_src$method), unique(idx_tgt$method), v23_variants, ot_variants$method)
methods <- scope$methods; training_free <- scope$training_free; trained <- scope$trained
cat("  methods in the transfer scope:", length(methods), "(", length(training_free), "training-free,", length(trained), "trained );",
    if (length(scope$outside)) paste("present in the runs but outside the scope (not transferred):", paste(scope$outside, collapse = ", ")) else "no method outside the scope", "\n")

# ------------------------------------------------------------------ DB3 bank (load only) and aligned cache
section("Target feature bank and aligned cache (DB3 cache identities)")
paths_t <- get_paths(tgt$config)
stopifnot(identical(basename(paths_t), tgt$files), identical(unname(tools::md5sum(paths_t)), tgt$md5))
bank_file <- file.path(TGT_DIR, "cache", paste0("bank_", substr(bank_cache_key(tgt$config), 1, 12), ".rds"))
side_file <- sub("bank_", "inputs_bank_", bank_file)
stopifnot(file.exists(bank_file), file.exists(side_file))
stopifnot(identical(readRDS(side_file), list(files = basename(paths_t), md5 = unname(tools::md5sum(paths_t)))))
cache_state <- function() file.info(list.files(file.path(TGT_DIR, "cache"), full.names = TRUE, recursive = TRUE))[, c("size", "mtime")]
state0 <- cache_state()
cfg_run <- tgt$config                                 # the replication's own configuration (OT chunking included, so the
cfg_run$n_workers <- CONFIG$n_workers                 # recomputation is as faithful as possible); n_workers is runtime-only
future::plan(future::multisession, workers = cfg_run$n_workers)
options(future.globals.maxSize = cfg_run$future_globals_gb * 1024^3)
bank_t <- load_or_compute_bank(paths_t, cfg_run)      # must LOAD (sidecar validated above); asserted below
stopifnot(identical(cache_state(), state0))           # nothing in Outputs_db3/cache was written or touched
dt_t    <- bank_to_matrix(bank_t, "dt_vec")
dt150_t <- bank_to_matrix(bank_t, "dt_vec_150")
betti_t <- bank_to_matrix(bank_t, "betti_global_100")
gpi_t   <- bank_to_matrix(bank_t, "pi_global")
v23_t   <- setNames(lapply(v23_variants, function(vv) get_v23_feature_matrix(bank_t, vv)), v23_variants)

stopifnot(identical(apply_frozen_grid(tgt$config)$aligned, tgt$config$aligned))   # replication ran on the frozen grid
al_key  <- aligned_cache_key(tgt$config$aligned, paths_t)
al_file <- file.path(TGT_DIR, "cache", paste0("aligned_", substr(al_key, 1, 12), ".rds"))
stopifnot(file.exists(al_file))
aligned_t <- readRDS(al_file)
stopifnot(identical(aligned_t$key, al_key), identical(aligned_t$files, basename(paths_t)))
al_pidx <- pair_index_matrix(length(paths_t), aligned_t$pairs)
cat("  bank", basename(bank_file), "loaded through its sidecar; aligned cache", basename(al_file), "key verified.\n")

# ------------------------------------------------------------------ helpers
pipeline_thresholds <- function(s_train, y_train) {                 # identical to run_or_load_method()
  imp <- s_train[y_train == 0]
  setNames(vapply(fars, function(f) as.numeric(quantile(imp, probs = 1 - f, na.rm = TRUE)), numeric(1)), far_names)
}
metrics_of <- function(y, s, thr) {
  mm <- roc_metrics(y, s)
  tibble(AUC = mm$auc, EER = mm$eer,
         !!!setNames(as.list(vapply(far_names, function(fn) mean(s[y == 1] >= thr[[fn]]), numeric(1))), paste0("TAR@", fars)),
         !!!setNames(as.list(vapply(far_names, function(fn) mean(s[y == 0] >= thr[[fn]]), numeric(1))), paste0("realisedFAR@", fars)))
}
same_vec <- function(a, b) length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)) && all(a == b)
frozen_fields <- function(cc) cc[intersect(names(cc), c("model", "best_w", "zfit_ot", "zfit_v"))]   # DB1 fitted objects, copied verbatim

OUTS <- list(frozen = OUT_FROZEN, recal = OUT_RECAL)
for (d in OUTS) dir.create(file.path(d, "cache", "folds"), recursive = TRUE, showWarnings = FALSE)
index_rows <- list(frozen = list(), recal = list())
metric_rows <- list()                                   # long table: fold, method, protocol, metrics
identical_scores <- setNames(rep(NA, length(methods)), methods)   # training-free: transfer test scores == replication's (bit-exact)

# ------------------------------------------------------------------ per fold
section("Fold-matched transfer: DB1 fold-k objects -> DB3 fold-k test pairs")
for (k in seq_len(n_folds)) {
  mf <- tgt$folds[[k]]
  pair_test <- mf$pair_test; pair_train <- mf$pair_train
  y_test  <- compute_y_from_pair_mat(pair_test,  tgt$ids)
  y_train <- compute_y_from_pair_mat(pair_train, tgt$ids)
  cat(sprintf("  fold %d: %d test pairs (%d genuine), %d training pairs for the secondary calibration\n", k, ncol(pair_test), sum(y_test == 1), ncol(pair_train)))
  SRC <- setNames(lapply(methods, function(m) load_cache(SRC_DIR, idx_src, k, m)), methods)
  REP <- setNames(lapply(methods, function(m) load_cache(TGT_DIR, idx_tgt, k, m)), methods)
  for (m in methods) stopifnot(identical(as.integer(REP[[m]]$y_test), as.integer(y_test)), identical(as.integer(REP[[m]]$y_train), as.integer(y_train)),
                               length(SRC[[m]]$s_train) == length(y_train), identical(as.integer(SRC[[m]]$y_train), as.integer(compute_y_from_pair_mat(src$folds[[k]]$pair_train, src$ids))))

  # --- test scores on DB3 under the DB1 objects
  s_test <- list()
  s_test[["DT_L2"]]             <- score_pairs_l2(pair_test, dt_t)
  s_test[["DT_L2_150"]]         <- score_pairs_l2(pair_test, dt150_t)
  s_test[["GlobalBetti_L1"]]    <- score_pairs_l1(pair_test, betti_t)
  s_test[["GlobalPI_L2"]]       <- score_pairs_l2(pair_test, gpi_t)
  s_test[["DT_L2_150_aligned"]] <- score_pairs_aligned(pair_test, aligned_t, al_pidx)
  for (ii in seq_len(nrow(ot_variants))) { vrow <- ot_variants[ii, ]; s_test[[vrow$method]] <- score_pairs_ot(pair_test, bank_t, cfg_run, vrow) }
  for (m in training_free) {                            # structural identity, asserted on the score vectors
    ok <- same_vec(s_test[[m]], REP[[m]]$s_test)
    identical_scores[m] <- if (is.na(identical_scores[m])) ok else identical_scores[m] && ok
    if (!ok) stop("training-free method ", m, " fold ", k, ": recomputed DB3 test scores are not bit-identical to the replication's.", call. = FALSE)
  }
  for (vv in v23_variants) { stopifnot(!is.null(SRC[[vv]]$model)); s_test[[vv]] <- predict_logistic_pairs(pair_test, v23_t[[vv]], SRC[[vv]]$model, chunk_size = 8000) }
  fu <- SRC[["Fusion"]]; fe <- SRC[["Fusion_eq"]]
  stopifnot(!is.null(fu$best_w), !is.null(fu$zfit_ot), !is.null(fu$zfit_v), identical(as.numeric(fe$best_w), 0.5), !is.null(fe$zfit_ot), !is.null(fe$zfit_v))
  s_test[["Fusion"]]    <- apply_fusion(s_test[["OT"]], s_test[["V23_full"]], fu$best_w, fu$zfit_ot, fu$zfit_v)
  s_test[["Fusion_eq"]] <- apply_fusion(s_test[["OT"]], s_test[["V23_full"]], 0.5,       fe$zfit_ot, fe$zfit_v)

  # --- secondary calibration: DB3 fold-k TRAINING pairs under the same frozen objects
  s_train_recal <- list()
  for (m in training_free) s_train_recal[[m]] <- REP[[m]]$s_train          # training-free: the replication's own training scores
  for (vv in v23_variants) s_train_recal[[vv]] <- predict_logistic_pairs(pair_train, v23_t[[vv]], SRC[[vv]]$model, chunk_size = 8000)
  s_train_recal[["Fusion"]]    <- apply_fusion(s_train_recal[["OT"]], s_train_recal[["V23_full"]], fu$best_w, fu$zfit_ot, fu$zfit_v)
  s_train_recal[["Fusion_eq"]] <- apply_fusion(s_train_recal[["OT"]], s_train_recal[["V23_full"]], 0.5,       fe$zfit_ot, fe$zfit_v)

  for (m in methods) {
    thr_frozen <- SRC[[m]]$thresholds[far_names]
    stopifnot(same_vec(unname(thr_frozen), unname(pipeline_thresholds(SRC[[m]]$s_train, SRC[[m]]$y_train))))   # stored DB1 thresholds == rule
    thr_recal  <- pipeline_thresholds(s_train_recal[[m]], y_train)
    thr_rep    <- REP[[m]]$thresholds[far_names]
    if (m %in% training_free) stopifnot(same_vec(unname(thr_recal), unname(thr_rep)))                         # recal == replication for training-free
    extra_key <- if (m == "DT_L2_150_aligned") al_key else NULL
    suffix    <- if (is.null(extra_key)) "" else paste0("_", substr(extra_key, 1, 8))
    base <- c(list(s_test = s_test[[m]], y_test = y_test, elapsed = NA_real_), frozen_fields(SRC[[m]]),
              if (!is.null(extra_key)) list(aligned_key = al_key, aligned_cache_file = basename(al_file)))
    objs <- list(
      frozen = c(base, list(s_train = SRC[[m]]$s_train, y_train = SRC[[m]]$y_train, thresholds = thr_frozen,
                            transfer = list(source_run = SRC_DIR, source_fold = k, target_run = TGT_DIR, target_fold = k,
                                            calibration = "frozen: DB1 fold-k thresholds; s_train/y_train are the DB1 fold-k training scores"))),
      recal  = c(base, list(s_train = s_train_recal[[m]], y_train = y_train, thresholds = thr_recal,
                            transfer = list(source_run = SRC_DIR, source_fold = k, target_run = TGT_DIR, target_fold = k,
                                            calibration = "target-calibrated: thresholds from DB3 fold-k training impostor scores under the frozen DB1 model"))))
    for (v in names(OUTS)) {
      fn <- sprintf("transfer_fold%d_%s_%s%s.rds", k, m, fkey_tgt, suffix)
      saveRDS(objs[[v]], file.path(OUTS[[v]], "cache", "folds", fn))
      index_rows[[v]][[length(index_rows[[v]]) + 1]] <- tibble(fold = k, method = m, file = fn, fold_key = fkey_tgt,
                                                              extra_key = extra_key %||% NA_character_, loaded_from_cache = FALSE)
    }
    metric_rows[[length(metric_rows) + 1]] <- bind_rows(
      mutate(metrics_of(y_test, REP[[m]]$s_test, thr_rep), protocol = "replication"),
      mutate(metrics_of(y_test, s_test[[m]],     thr_frozen), protocol = "transfer_frozen"),
      mutate(metrics_of(y_test, s_test[[m]],     thr_recal),  protocol = "transfer_recal")) %>%
      mutate(fold = k, method = m, .before = 1)
  }
}
future::plan(future::sequential)
stopifnot(identical(cache_state(), state0), all(identical_scores[training_free]))
cat("  training-free methods: transfer test scores bit-identical to the replication's in every fold (structural identity).\n")

# ------------------------------------------------------------------ run-shaped directories
section("Writing run-shaped transfer directories")
fm <- bind_rows(metric_rows)
for (v in names(OUTS)) {
  d <- OUTS[[v]]
  write_csv(bind_rows(index_rows[[v]]), file.path(d, "fold_cache_index.csv"))
  proto <- if (v == "frozen") "transfer_frozen" else "transfer_recal"
  sub_ <- fm %>% filter(protocol == proto)
  write_csv(sub_ %>% select(fold, method, AUC, EER), file.path(d, "fold_metrics_auc_eer.csv"))
  write_csv(sub_ %>% select(fold, method, starts_with("TAR@")) %>%
              pivot_longer(starts_with("TAR@"), names_to = "FAR", values_to = "TAR") %>%
              mutate(FAR = as.numeric(sub("^TAR@", "", FAR))) %>% arrange(fold, match(method, methods), FAR),
            file.path(d, "fold_metrics_tar.csv"))
  man <- tgt
  man$created <- Sys.time()
  man$transfer <- list(kind = "DB1 -> DB3 transfer (fold-matched)", calibration = v, source_run = SRC_DIR, target_run = TGT_DIR,
                       source_manifest_created = src$created, target_manifest_created = tgt$created,
                       spec = "Pattern Recognition/Revision/STEP_1.1_DB3_SPEC.md section 5", created = Sys.time())
  saveRDS(man, file.path(d, "run_manifest.rds"))
  writeLines(c(
    "TRANSFER DIRECTORY - NOT A PIPELINE RUN",
    sprintf("Written by Code/transfer_db1_to_db3.R on %s.", format(Sys.time())),
    paste0("DB1 -> DB3 transfer test, fold-matched: the fitted objects of DB1 fold k (", SRC_DIR, ") applied to the DB3 fold-k test pairs (", TGT_DIR, ")."),
    sprintf("Calibration: %s.", if (v == "frozen") "FROZEN - DB1 fold-k thresholds; cache s_train/y_train are the DB1 fold-k training scores (so the stored thresholds equal their training quantiles)"
                                 else "TARGET-CALIBRATED - thresholds re-estimated on the DB3 fold-k training impostor scores under the frozen DB1 model; cache s_train/y_train are those DB3 training scores"),
    "Nothing was refitted, retuned or re-standardised; model/best_w/zfit_* are the DB1 objects, copied verbatim.",
    "The manifest is the DB3 replication manifest (ids, folds, pairs) plus a `transfer` element; the layout exists so that Code/bootstrap_ci.R can analyse it unchanged.",
    "Training-free methods: transfer test scores are bit-identical to the replication's (asserted); only the thresholds differ in the frozen variant, nothing in the recal variant.",
    "Design: Pattern Recognition/Revision/STEP_1.1_DB3_SPEC.md section 5."), file.path(d, "TRANSFER_README.txt"))
  cat("  ", d, ": ", nrow(bind_rows(index_rows[[v]])), " caches, index, manifest, metric CSVs, README\n", sep = "")
}
write_csv(fm, file.path(OUT_FROZEN, "transfer_fold_metrics.csv"))
saveRDS(list(created = Sys.time(), source_run = SRC_DIR, target_run = TGT_DIR, assignment = "fold-matched (DB1 fold k -> DB3 fold k)",
             methods = methods, training_free = training_free, trained = trained, outside_scope = scope$outside, identical_scores = identical_scores,
             aligned_key = al_key, bank_file = basename(bank_file), fold_key = fkey_tgt, fars = fars,
             elapsed_min = as.numeric(difftime(Sys.time(), t_start, units = "mins"))),
        file.path(OUT_FROZEN, "transfer_report.rds"))

# ------------------------------------------------------------------ point estimates (intervals come from the bootstrap)
section("5-fold means: replication vs transfer (frozen thresholds) vs transfer (target-calibrated) - point estimates only")
summ <- fm %>% group_by(method, protocol) %>% summarize(across(c(AUC, EER, `TAR@0.001`, `realisedFAR@0.001`, `TAR@0.01`, `realisedFAR@0.01`), mean), .groups = "drop")
show <- c("DT_L2_150", "DT_L2_150_aligned", "GlobalBetti_L1", "GlobalPI_L2", "V23_full", "OT", "Fusion", "Fusion_eq")
out <- summ %>% filter(method %in% show) %>% mutate(method = factor(method, levels = show), protocol = factor(protocol, levels = c("replication", "transfer_frozen", "transfer_recal"))) %>%
  arrange(method, protocol) %>% mutate(across(where(is.numeric), ~ round(.x, 4)))
print(as.data.frame(out), row.names = FALSE)
cat("\n  (training-free methods: AUC/EER identical across protocols by construction; recal == replication.)\n")
cat(sprintf("\nDone in %.1f minutes. Next:\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("  options(fp.bootstrap_dir = \"Outputs_db3_transfer\");       source(\"Code/bootstrap_ci.R\")\n")
cat("  options(fp.bootstrap_dir = \"Outputs_db3_transfer_recal\"); source(\"Code/bootstrap_ci.R\")\n")
cat("  source(\"Code/transfer_vs_replication.R\")\n")
