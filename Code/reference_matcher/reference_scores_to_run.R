############################################################
# STEP 1.3 (C3, R4-7) — fold the NBIS reference-matcher scores into a RUN-SHAPED directory
# next to a pipeline run, so that Code/bootstrap_ci.R gives intervals and paired differences
# for the reference rows and the pre-named fusions on the SAME folds, pairs and resamples.
#
# Inputs (options; defaults = DB1):
#   fp.ref_src          pipeline run (all its methods are carried over) (default "Outputs_regen_v2"; DB3: "Outputs_db3_v2")
#   fp.ref_scores_gray  run_nbis_matcher.py output, grayscale kind    (default "Outputs_reference/nbis_db1_gray")
#   fp.ref_scores_bin   run_nbis_matcher.py output, binary kind       (default "Outputs_reference/nbis_db1_bin"; "" = none)
#   fp.ref_out          run-shaped output directory                   (default "Outputs_reference_db1_v2")
# DB3:  options(fp.ref_src = "Outputs_db3_v2", fp.ref_scores_gray = "Outputs_reference/nbis_db3_gray",
#               fp.ref_scores_bin = "Outputs_reference/nbis_db3_bin", fp.ref_out = "Outputs_reference_db3_v2", fp.ref_raw_dir = "Data/FVC2002_DB3_raw")
#
# Methods added (names pre-specified in STEP_1.3_REFERENCE_MATCHER_SPEC.md):
#   NBIS_gray             mindtct + bozorth3 on the raw grayscale images; pair score = mean of the two orderings  (PRIMARY reference)
#   NBIS_bin              the same matcher given OUR binary ridge maps (inverted)                                   (secondary; caveat)
#   Fus_NBISgray_TopoLR   equal-weight (w = 0.5) robust-z fusion of NBIS_gray with TopoLR (V23_full), pipeline rule
#   Fus_NBISgray_FusionEq equal-weight robust-z fusion of NBIS_gray with Fusion_eq
# Thresholds for every added method: quantile(training impostor scores, 1 - FAR), type 7 (the pipeline rule);
# bozorth3 scores are integers, so realised FAR at these thresholds is the operating-point check.
# Every pipeline cache of the source run is COPIED in unchanged (same fold key), with its index rows and metric rows
# (16 methods in the submitted pipeline; 17 with DT_L2_150_clamped after the audit patch of 31 Oct 2026).
#
# Run from the PROJECT ROOT (~1 min):   source("Code/reference_matcher/reference_scores_to_run.R")
# Then:  options(fp.bootstrap_dir = "Outputs_reference_db1_v2"); source("Code/bootstrap_ci.R")
############################################################

options(fp.suppress_main = TRUE)
source("Code/fingerprintTDA.R")

SRC_DIR  <- getOption("fp.ref_src", "Outputs_regen_v2")
GRAY_DIR <- getOption("fp.ref_scores_gray", "Outputs_reference/nbis_db1_gray")
BIN_DIR  <- getOption("fp.ref_scores_bin",  "Outputs_reference/nbis_db1_bin")
OUT_DIR  <- getOption("fp.ref_out", "Outputs_reference_db1_v2")
RAW_DIR  <- getOption("fp.ref_raw_dir", "Data/FVC2000_raw")      # the grayscale images that were scored (md5 cross-check); DB3: "Data/FVC2002_DB3_raw"
section <- function(s) cat("\n", s, "\n", strrep("-", nchar(s)), "\n", sep = "")
if (dir.exists(OUT_DIR) && length(list.files(OUT_DIR, recursive = TRUE, all.files = TRUE, no.. = TRUE)) > 0)
  stop(OUT_DIR, " already contains files; use a new or empty directory.", call. = FALSE)

# ------------------------------------------------------------------ source run
section("Source run")
man <- readRDS(file.path(SRC_DIR, "run_manifest.rds")); CFG <- man$config
fars <- CFG$fars_report; far_names <- paste0("far_", fars); n_folds <- CFG$n_folds
fkey <- substr(fold_cache_key(CFG), 1, 8)
idx <- readr::read_csv(file.path(SRC_DIR, "fold_cache_index.csv"), show_col_types = FALSE,
                       col_types = readr::cols(.default = readr::col_character(), fold = readr::col_integer(), loaded_from_cache = readr::col_logical()))
stopifnot(all(idx$fold_key == fkey), !any(duplicated(idx[, c("fold", "method")])))
load_cache <- function(f, m) { hit <- idx[idx$fold == f & idx$method == m, ]; stopifnot(nrow(hit) == 1); readRDS(file.path(SRC_DIR, "cache", "folds", hit$file)) }
stems <- sub("_bin\\.png$", "", man$files); stopifnot(!any(duplicated(stems)), length(stems) == 800)
cat("  ", SRC_DIR, ": ", length(man$files), " images, ", length(unique(idx$method)), " methods, fold key ", fkey, "\n", sep = "")

# ------------------------------------------------------------------ reference score matrices
read_scores <- function(d, kind) {
  mj <- jsonlite::fromJSON(file.path(d, "manifest.json")); stopifnot(mj$kind == kind, mj$n_images == 800, mj$fingers == 100)
  im <- readr::read_csv(file.path(d, "images.csv"), show_col_types = FALSE, col_types = readr::cols(stem = "c", file = "c", md5 = "c", nfiq = "c"))
  sc <- readr::read_csv(file.path(d, "scores_pairs.csv"), show_col_types = FALSE)
  # Verify the files against the manifest BEFORE anything is used (QA, 13 Oct): recorded md5s of all
  # three CSVs; every pair index, file name and finger label; finite scores; stored mean == average
  # of the two orderings; genuine flag consistent. Protects against files mixed up in transfer.
  for (fn in c("images.csv", "scores_pairs.csv", "self_scores.csv"))
    if (!identical(unname(tools::md5sum(file.path(d, fn))), mj$outputs[[fn]])) stop(fn, " in ", d, " does not match the md5 recorded in its manifest", call. = FALSE)
  stopifnot(nrow(im) == 800, nrow(sc) == 800 * 799 / 2, setequal(im$stem, stems), identical(as.integer(im$index), 1:800))
  ex <- t(utils::combn(800L, 2L))
  stopifnot(identical(as.integer(sc$i), ex[, 1]), identical(as.integer(sc$j), ex[, 2]),
            identical(sc$file_i, im$file[sc$i]), identical(sc$file_j, im$file[sc$j]),
            identical(as.integer(sc$finger_i), as.integer(im$finger[sc$i])), identical(as.integer(sc$finger_j), as.integer(im$finger[sc$j])),
            identical(as.integer(sc$genuine), as.integer(sc$finger_i == sc$finger_j)),
            all(is.finite(sc$score_ij)), all(is.finite(sc$score_ji)), all(is.finite(sc$score_mean)),
            all(sc$score_mean == (sc$score_ij + sc$score_ji) / 2))
  # The matcher's table is in numeric (finger, impression) order; the pipeline manifest is in R's
  # lexicographic file order. Map every pipeline position to its row in the matcher's table.
  pos <- match(stems, im$stem); stopifnot(!anyNA(pos))
  if (kind == "grayscale") {                                                      # md5s of the raw images actually scored
    stopifnot(identical(unname(tools::md5sum(file.path(RAW_DIR, im$file))), im$md5))
  } else stopifnot(identical(im$md5[pos], man$md5))                               # binary kind: exactly the pipeline's input files
  Mref <- matrix(NA_real_, 800, 800); Mref[cbind(sc$i, sc$j)] <- sc$score_mean; Mref[cbind(sc$j, sc$i)] <- sc$score_mean
  stopifnot(!anyNA(Mref[upper.tri(Mref)]))
  M <- Mref[pos, pos]                                                             # now indexed like the pipeline manifest
  list(M = M, manifest = mj, images = im, pos = pos)
}
section("Reference scores")
gray <- read_scores(GRAY_DIR, "grayscale"); cat("  grayscale: ", GRAY_DIR, " (", gray$manifest$nbis$version, "; minutiae median ", gray$manifest$minutiae$median,
                                                 ", ", gray$manifest$minutiae$n_below_10, " image(s) below 10; orderings identical in ", round(100 * gray$manifest$pairs$frac_orderings_identical, 1), " % of pairs)\n", sep = "")
bin <- if (nzchar(BIN_DIR)) read_scores(BIN_DIR, "binary") else NULL
if (!is.null(bin)) cat("  binary   : ", BIN_DIR, " (minutiae median ", bin$manifest$minutiae$median, ", ", bin$manifest$minutiae$n_below_10, " image(s) below 10)\n", sep = "")

# ------------------------------------------------------------------ build the run-shaped directory
section("Writing the run-shaped directory")
dir.create(file.path(OUT_DIR, "cache", "folds"), recursive = TRUE, showWarnings = FALSE)
pipeline_thresholds <- function(s_train, y_train) { imp <- s_train[y_train == 0]
  setNames(vapply(fars, function(f) as.numeric(quantile(imp, probs = 1 - f, na.rm = TRUE)), numeric(1)), far_names) }
metrics_of <- function(y, s, thr) { mm <- roc_metrics(y, s)
  tibble(AUC = mm$auc, EER = mm$eer,
         !!!setNames(as.list(vapply(far_names, function(fn) mean(s[y == 1] >= thr[[fn]]), numeric(1))), paste0("TAR@", fars)),
         !!!setNames(as.list(vapply(far_names, function(fn) mean(s[y == 0] >= thr[[fn]]), numeric(1))), paste0("realisedFAR@", fars))) }
# copy the pipeline caches and index rows unchanged
file.copy(file.path(SRC_DIR, "cache", "folds", idx$file), file.path(OUT_DIR, "cache", "folds", idx$file))
index_rows <- list(idx)
auc_rows <- list(readr::read_csv(file.path(SRC_DIR, "fold_metrics_auc_eer.csv"), show_col_types = FALSE))
tar_rows <- list(readr::read_csv(file.path(SRC_DIR, "fold_metrics_tar.csv"), show_col_types = FALSE))
ref_rows <- list()
add_method <- function(f, m, s_train, s_test, y_train, y_test, extras = list()) {
  thr <- pipeline_thresholds(s_train, y_train)
  obj <- c(list(s_train = s_train, s_test = s_test, y_train = y_train, y_test = y_test, elapsed = NA_real_, thresholds = thr), extras)
  fn <- sprintf("fold%d_%s_%s.rds", f, m, fkey); saveRDS(obj, file.path(OUT_DIR, "cache", "folds", fn))
  index_rows[[length(index_rows) + 1]] <<- tibble(fold = f, method = m, file = fn, fold_key = fkey, extra_key = NA_character_, loaded_from_cache = FALSE)
  met <- metrics_of(y_test, s_test, thr)
  auc_rows[[length(auc_rows) + 1]] <<- tibble(fold = f, method = m, AUC = met$AUC, EER = met$EER)
  tar_rows[[length(tar_rows) + 1]] <<- tibble(fold = f, method = m, FAR = fars, TAR = as.numeric(met[paste0("TAR@", fars)]))
  ref_rows[[length(ref_rows) + 1]] <<- mutate(met, fold = f, method = m, .before = 1)
}
for (f in seq_len(n_folds)) {
  mf <- man$folds[[f]]; ptr <- mf$pair_train; pte <- mf$pair_test
  y_train <- compute_y_from_pair_mat(ptr, man$ids); y_test <- compute_y_from_pair_mat(pte, man$ids)
  v23 <- load_cache(f, "V23_full"); feq <- load_cache(f, "Fusion_eq")
  stopifnot(identical(as.integer(v23$y_train), y_train), identical(as.integer(v23$y_test), y_test))
  g_tr <- gray$M[cbind(ptr[1, ], ptr[2, ])]; g_te <- gray$M[cbind(pte[1, ], pte[2, ])]
  add_method(f, "NBIS_gray", g_tr, g_te, y_train, y_test, list(reference = "NBIS 5.0.0 mindtct+bozorth3, raw grayscale, mean of orderings"))
  if (!is.null(bin)) add_method(f, "NBIS_bin", bin$M[cbind(ptr[1, ], ptr[2, ])], bin$M[cbind(pte[1, ], pte[2, ])], y_train, y_test,
                                list(reference = "NBIS 5.0.0 mindtct+bozorth3 on the pipeline's binary ridge maps (inverted), mean of orderings"))
  z_g <- robust_z_fit_train_impostors(g_tr, y_train)
  for (comp in c("TopoLR", "FusionEq")) {
    cc <- if (comp == "TopoLR") v23 else feq
    z_c <- robust_z_fit_train_impostors(cc$s_train, y_train)
    add_method(f, paste0("Fus_NBISgray_", comp),
               apply_fusion(g_tr, cc$s_train, 0.5, z_g, z_c), apply_fusion(g_te, cc$s_test, 0.5, z_g, z_c), y_train, y_test,
               list(best_w = 0.5, zfit_ot = z_g, zfit_v = z_c, components = c("NBIS_gray", if (comp == "TopoLR") "V23_full" else "Fusion_eq"),
                    note = "apply_fusion(s_a = NBIS_gray, s_b = component, w = 0.5): w weights the component's robust z, 1 - w the matcher's"))
  }
  cat(sprintf("  fold %d: NBIS_gray train/test pairs %d/%d, thresholds %s\n", f, length(g_tr), length(g_te),
              paste(sprintf("%.1f", pipeline_thresholds(g_tr, y_train)), collapse = " / ")))
}
readr::write_csv(bind_rows(index_rows), file.path(OUT_DIR, "fold_cache_index.csv"))
readr::write_csv(bind_rows(auc_rows), file.path(OUT_DIR, "fold_metrics_auc_eer.csv"))
readr::write_csv(bind_rows(tar_rows), file.path(OUT_DIR, "fold_metrics_tar.csv"))
ref_df <- bind_rows(ref_rows); readr::write_csv(ref_df, file.path(OUT_DIR, "reference_fold_metrics.csv"))
man_out <- man; man_out$created <- Sys.time()
man_out$reference <- list(kind = "NBIS reference matcher rows + pre-named equal-weight fusions added to the pipeline run", source_run = SRC_DIR,
                          grayscale_scores = GRAY_DIR, binary_scores = if (!is.null(bin)) BIN_DIR else NA, gray_manifest = gray$manifest,
                          bin_manifest = if (!is.null(bin)) bin$manifest else NULL, spec = "Pattern Recognition/Revision/STEP_1.3_REFERENCE_MATCHER_SPEC.md", created = Sys.time())
saveRDS(man_out, file.path(OUT_DIR, "run_manifest.rds"))
writeLines(c("REFERENCE-MATCHER DIRECTORY - NOT A PIPELINE RUN",
             sprintf("Written by Code/reference_matcher/reference_scores_to_run.R on %s from %s.", format(Sys.time()), SRC_DIR),
             sprintf("The %d pipeline caches per fold are copied unchanged; NBIS_gray, NBIS_bin and the two equal-weight fusions are added with pipeline-rule thresholds.", length(unique(idx$method))),
             "Layout exists so that Code/bootstrap_ci.R can analyse all methods on identical folds, pairs and resamples.",
             "Design and pre-specified comparisons: Pattern Recognition/Revision/STEP_1.3_REFERENCE_MATCHER_SPEC.md"), file.path(OUT_DIR, "REFERENCE_README.txt"))

section("5-fold means of the added methods (point estimates; intervals from the bootstrap)")
print(as.data.frame(ref_df %>% group_by(method) %>% summarize(across(c(AUC, EER, `TAR@0.001`, `realisedFAR@0.001`, `TAR@0.01`, `realisedFAR@0.01`), mean), .groups = "drop") %>%
                      mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)
cat("\nNext:  options(fp.bootstrap_dir = \"", OUT_DIR, "\"); source(\"Code/bootstrap_ci.R\")\n", sep = "")
