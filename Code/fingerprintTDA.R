library(tidyverse)
library(EBImage)
library(TDA)
library(TDAvec)
library(pROC)
library(future)
library(future.apply)
library(progressr)
library(parallel)
library(matrixStats)
library(digest)

# Aligned geometry baseline engine (revision step 1.2). Like every other path in
# this file, resolved relative to the PROJECT ROOT (the working directory); the
# lowercase fallback mirrors pca_match_analysis.R / validate_downsampling.R.
.aligned_src <- if (file.exists(file.path("Code", "aligned_baseline.R"))) file.path("Code", "aligned_baseline.R") else file.path("code", "aligned_baseline.R")
if (!file.exists(.aligned_src)) stop("aligned_baseline.R not found under Code/ (or code/): source this file from the project root.", call. = FALSE)
source(.aligned_src)
rm(.aligned_src)

############################################################
# TDA FINGERPRINT PIPELINE (Pattern Recognition revision, 2026)
#
# Entry points: Code/run_canonical.R (canonical DB1 / DB3 runs on the regenerated binaries, frozen
# alignment grid), then Code/bootstrap_ci.R. The main entry at the end of this file is the legacy
# single-call driver and is suppressed with options(fp.suppress_main = TRUE) by every driver.
#
# Audit patch (31 Oct 2026): persistence failures are propagated (never vectorised as zero features);
# anchors are drawn under an explicitly set RNG kind; ot$geo_scale enters the bank key when changed;
# matched-transformation geometry baseline DT_L2_150_clamped added (bank field dt_vec_150_clamped).
#
# Historical note (first version): "MERGED V2 (PRL-READY, NO LEAKAGE); speedups 1-6":
# Changes from original:
#  1. Vectorized baseline pair scoring (score_pairs_l2, score_pairs_l1)
#  2. Increased OT chunk size (300 -> 5000)
#  3. Avoid allocating P matrix in Sinkhorn (compute cost directly)
#  4. Sinkhorn early stopping (OPT-IN: set sinkhorn_tol > 0 to enable)
#  5. Replaced outer() calls with matrix() construction
#  6. Reduced V23 standardizer sample size (50k -> 20k)
#
# NOTE: Early stopping (speedup #4) is DISABLED by default (tol=0)
#       to preserve exact reproducibility with fixed-iteration Sinkhorn.
#       Set sinkhorn_tol > 0 (e.g., 1e-4) to enable adaptive stopping.
############################################################

# ============================================================
# CONFIG
# ============================================================

# Release 1.1.0 (1 Nov 2026): the defaults below are the reported configuration — the regenerated
# binaries (Code/preprocess.py), the frozen alignment grid (Code/frozen_grid.R) and no manuscript export —
# so that Code/run_canonical.R reproduces the paper without hidden overrides. The drivers still set the
# directories explicitly; Code/test_release_config.R asserts that the effective configuration and every
# cache key of the reported runs are unchanged by these defaults.
.frozen_src <- if (file.exists(file.path("Code", "frozen_grid.R"))) file.path("Code", "frozen_grid.R") else file.path("code", "frozen_grid.R")
if (!file.exists(.frozen_src)) stop("frozen_grid.R not found under Code/ (or code/): source this file from the project root.", call. = FALSE)
source(.frozen_src)
rm(.frozen_src)

CONFIG <- list(
  binary_dir = "Data/FVC2000_binary_regen",
  n_id = 100,
  
  # CV
  n_folds = 5,
  cv_seed = 20260110,
  
  # Common PH / geometry
  maxscale = 20,
  target_size = 150,
  
  # Baselines
  baselines = list(
    dt_downsample_size = 50,
    global_pi_res = 20,
    global_pi_sigma = 0.5
  ),
  
  # OT (frozen method definition)
  ot = list(
    n_anchors = 64,
    radius = 12,
    pi_res = 30,
    pi_sigma = 0.25,
    use_dim = 1,
    sinkhorn_reg = 0.12,
    sinkhorn_iter = 60,
    sinkhorn_tol = 0,              # 0 = disabled; set > 0 (e.g. 1e-4) to enable early stopping
    sinkhorn_check_every = 10,     # Only used when sinkhorn_tol > 0
    dustbin = TRUE,
    dustbin_cost = 0.4,
    eps_dustbin = 0.20,
    beta_geo = 6,
    geo_scale = "target",
    chunk_size_pairs = 5000
  ),
  
  # OT ablations (scoring-time changes only)
  ot_variants = tribble(
    ~method,          ~use_ridge, ~use_valley, ~beta_geo, ~dustbin,
    "OT",             TRUE,       TRUE,        6,         TRUE,
    "OT_noGeo",       TRUE,       TRUE,        0,         TRUE,
    "OT_noDustbin",   TRUE,       TRUE,        6,         FALSE,
    "OT_ridgeOnly",   TRUE,       FALSE,       6,         TRUE,
    "OT_valleyOnly",  FALSE,      TRUE,        6,         TRUE
  ),
  
  # V23
  v23 = list(
    n_points = 50,
    n_anchors = 64,
    radius = 12,
    top_k = 10,
    symmetric_ops = c("diff", "min", "max", "product"),
    lambda = 0.01,
    max_iter = 200,
    tol = 0,                       # 0 = disabled; set > 0 to enable early stopping
    batch_size = 10000,
    lr0 = 0.1,
    seed = 123,
    n_pairs_for_standardize = 20000L
  ),
  
  v23_variants = c("V23_full", "V23_noTopK", "V23_globalOnly", "V23_localOnly"),
  
  # Fusion
  fusion_w_grid = seq(0, 1, by = 0.02),
  
  # FAR operating points
  fars_report = c(1e-3, 2e-3, 5e-3, 1e-2),
  
  # Parallel
  n_workers = { nc <- suppressWarnings(parallel::detectCores()); if (is.na(nc)) 4L else max(1L, as.integer(nc) - 2L) },   # operational, not scientific; NA-safe
  future_globals_gb = 32,

  # Aligned geometry baseline (revision step 1.2): rigid-search grid for
  # DT_L2_150_aligned. Developed on DB1 under the pre-specified expansion rule,
  # then FROZEN before any other database is run (Code/frozen_grid.R: theta_max 50,
  # t_max 40; the engine's own default stays +-30 degrees). Not part of fold_cache_key():
  # the aligned method's fold caches carry their own identity (the pair-cache key).
  aligned = apply_frozen_grid(list(aligned = aligned_default_config()))$aligned,

  # Output
  output_dir = "Outputs",
  # Copy roc_pooled.pdf into Manuscript/ at the end of a run. FALSE: reproduction and
  # experimental runs never touch a manuscript figure (the figures come from make_figures.R).
  export_to_manuscript = FALSE,
  # Bump to invalidate per-fold score caches. v2: caches carry fitted objects
  # (TopoLR model, fusion normalisation, thresholds) for fixed-model reuse.
  fold_cache_version = 2L,

  # Filtration complex for sublevel-set persistence (revision step 1.5, D2/D3):
  #   "triangulated" (default; every reported result) -> TDA::gridDiag, which builds a
  #                   lower-star filtration on the triangulated pixel grid;
  #   "cubical"      -> GUDHI Bitmap_cubical_complex (pixels as top cells, lower cells take
  #                   the minimum), via Code_cubical/cubical_ph.R + reticulate; supplementary
  #                   side-by-side only. Enters bank_cache_key() only when not "triangulated",
  #                   so every existing cache key is unchanged.
  complex_type = "triangulated",

  # Polarity rule of the binary reader (revision step 1.4): "auto" (default; every reported run)
  # inverts an image with > 50 % white pixels; "ridge255" trusts white = ridge and never inverts
  # (used only by the sensitivity experiment). Enters bank_cache_key() only when not "auto".
  polarity = "auto",

  # Debug
  debug_fast = FALSE,
  debug_n_id = 20,
  debug_max_pairs_train = 20000,
  debug_max_pairs_test  = 20000
)

# ============================================================
# BASIC HELPERS
# ============================================================

get_identity <- function(filename) {
  basename(filename) %>%
    str_extract("^\\d+") %>%
    { if_else(is.na(.) | . == "", NA_character_, .) }
}

seed_from_string <- function(s) {
  # A SHA-256 hash of the input string, with the leading 7 hex digits reduced
  # to a positive 32-bit integer seed. The previous "sum of UTF-8 bytes"
  # scheme collided extensively on FVC2000 filenames (e.g., "1_2_bin.png"
  # and "2_1_bin.png" hashed to the same value), violating the
  # one-seed-per-image design.
  h <- digest::digest(s, algo = "sha256", serialize = FALSE)
  as.integer((strtoi(substr(h, 1, 7), base = 16L) %% 2147483646L) + 1L)
}

read_binary_fingerprint <- function(path, polarity = "auto") {
  img <- EBImage::readImage(path)
  if (length(dim(img)) == 3) img <- EBImage::channel(img, "gray")
  mat <- as.matrix(img)
  mat[!is.finite(mat)] <- 0
  bin <- as.integer(mat >= 0.5)
  dim(bin) <- dim(mat)
  # polarity (revision step 1.4): "auto" (default; every reported run) inverts an image whose
  # white fraction exceeds 0.5, a guard against inverted inputs; "ridge255" trusts the fixed
  # convention white = ridge and never inverts. The sensitivity experiment uses "ridge255"
  # because dilation pushes 113 of 800 canonical images past 0.5 white, which under "auto"
  # would silently swap ridge and valley (QA, 13 Oct). No canonical image exceeds 0.5, so the
  # two settings coincide on clean inputs (asserted by the sensitivity scorer).
  if (identical(polarity, "auto")) { if (mean(bin) > 0.5) bin <- 1L - bin }
  else if (!identical(polarity, "ridge255")) stop("unknown polarity: ", polarity, call. = FALSE)
  bin
}

compute_dt_to_ridge <- function(bin) {
  inv <- 1L - bin
  storage.mode(inv) <- "integer"
  as.matrix(EBImage::distmap(inv))
}

compute_dt_to_valley <- function(bin) {
  storage.mode(bin) <- "integer"
  as.matrix(EBImage::distmap(bin))
}

downsample_matrix <- function(mat, target_size) {
  nr <- nrow(mat); nc <- ncol(mat)
  if (nr <= target_size && nc <= target_size) return(mat)
  step <- max(ceiling(nr / target_size), ceiling(nc / target_size))
  mat[seq(1, nr, by = step), seq(1, nc, by = step), drop = FALSE]
}

pad_or_crop_to_size <- function(mat, target) {
  out <- matrix(0, target, target)
  nr <- min(nrow(mat), target)
  nc <- min(ncol(mat), target)
  out[1:nr, 1:nc] <- mat[1:nr, 1:nc]
  out
}

make_disk_mask <- function(r) {
  xs <- -r:r
  ys <- -r:r
  grid <- expand.grid(x = xs, y = ys)
  matrix((grid$x^2 + grid$y^2) <= r^2, nrow = length(xs), byrow = FALSE)
}

extract_patch_with_mask <- function(mat, cx, cy, r, mask_full) {
  nr <- nrow(mat); nc <- ncol(mat)
  
  r1 <- max(1, cx - r); r2 <- min(nr, cx + r)
  c1 <- max(1, cy - r); c2 <- min(nc, cy + r)
  
  patch <- mat[r1:r2, c1:c2, drop = FALSE]
  
  mr1 <- 1 + (r1 - (cx - r))
  mr2 <- mr1 + (r2 - r1)
  mc1 <- 1 + (c1 - (cy - r))
  mc2 <- mc1 + (c2 - c1)
  
  mask <- mask_full[mr1:mr2, mc1:mc2, drop = FALSE]
  list(patch = patch, mask = mask)
}

# Anchor RNG kind (audit patch, 31 Oct 2026). The feature bank is computed inside
# future_lapply(..., future.seed = TRUE), whose workers run under L'Ecuyer-CMRG, so every reported
# bank drew its anchors with set.seed(seed) under THAT kind. A serial call (figure script, validators,
# tests) used to run under the session default (Mersenne-Twister) and drew DIFFERENT anchors from the
# same seed. The kind is now set explicitly (and restored afterwards), so the draw no longer depends
# on the caller's RNG state; production banks are unchanged because the kind is the one they used.
# The normal and sample kinds are pinned too ("Inversion", "Rejection": R's defaults, in force in the workers).
ANCHOR_RNG_KIND <- "L'Ecuyer-CMRG"

sample_anchor_points <- function(bin_ds, n_anchors, seed) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  suppressWarnings(RNGkind(ANCHOR_RNG_KIND, "Inversion", "Rejection"))   # kind, normal.kind AND sample.kind pinned (a caller's sample.kind = "Rounding" would change the draw)
  set.seed(seed)
  idx <- which(bin_ds == 1L, arr.ind = TRUE)
  if (nrow(idx) == 0) return(matrix(numeric(0), ncol = 2))
  pick <- sample.int(nrow(idx), n_anchors, replace = (nrow(idx) < n_anchors))
  idx[pick, , drop = FALSE]
}

# The anchors of a cached bank entry, recovered from its stored OT coordinates (centred on the
# anchor mean and divided by target_size under the default geo_scale): rows = (row, col) on the
# downsampled grid. Used by the figure script so that illustrations show the analysis anchors.
anchors_from_bank_entry <- function(entry, config) {
  Z <- entry$ot_Z; ctr <- attr(Z, "scaled:center")
  if (is.null(ctr)) stop("bank entry has no anchor centre attribute", call. = FALSE)
  scale_by <- if (identical(config$ot$geo_scale %||% "target", "target")) config$target_size else config$ot$radius
  A <- sweep(Z * scale_by, 2, ctr, "+")
  if (max(abs(A - round(A))) > 1e-6) stop("recovered anchor coordinates are not integers (geo_scale or centre mismatch?)", call. = FALSE)
  matrix(as.integer(round(A)), ncol = 2L)      # a plain matrix: no dimnames, no scaled:center attribute
}

compute_persistence <- function(filt_mat, maxdim = 1, maxscale = 20, complex_type = "triangulated") {
  x <- pmin(pmax(filt_mat, 0), maxscale)
  x[!is.finite(x)] <- maxscale
  x <- round(x, 1)

  # A backend failure is PROPAGATED on BOTH constructions (cubical since QA 13 Oct; triangulated
  # since the audit patch of 31 Oct 2026 — before that the triangulated path returned NULL, which
  # the vectorisers turned into all-zero features, indistinguishable from an image without
  # features). precompute_feature_bank() reports the failing image and stops; no bank is written.
  # A valid diagram with no H1 features is not an error and is accepted as such.
  # Revision step 1.5: optional cubical construction (same clipped, rounded input; same P x 3
  # diagram layout). cubical_diagram() comes from Code_cubical/cubical_ph.R, which the cubical
  # driver sources.
  if (identical(complex_type, "cubical")) {
    if (!exists("cubical_diagram", mode = "function"))
      stop("complex_type = 'cubical' requires Code_cubical/cubical_ph.R to be sourced first.", call. = FALSE)
    return(tryCatch(cubical_diagram(x, maxdim = maxdim, inf_replace = maxscale),
                    error = function(e) stop("cubical persistence failed: ", conditionMessage(e), call. = FALSE)))
  }
  if (!identical(complex_type, "triangulated")) stop("unknown complex_type: ", complex_type, call. = FALSE)

  d <- tryCatch({
    TDA::gridDiag(
      FUNvalues = x,
      maxdimension = maxdim,
      sublevel = TRUE,
      library = "GUDHI",
      printProgress = FALSE
    )$diagram
  }, error = function(e) stop("triangulated persistence (TDA::gridDiag) failed: ", conditionMessage(e), call. = FALSE))
  if (is.null(d)) stop("triangulated persistence (TDA::gridDiag) returned no diagram", call. = FALSE)
  d
}

compute_betti_curve <- function(diag, hom_dim, n_points, maxscale) {
  curve <- rep(0, n_points)
  # A NULL diagram is a propagated-failure marker, never a valid empty diagram (gridDiag always
  # returns at least the H0 class); it must not be vectorised as a zero curve (audit, 31 Oct 2026).
  if (is.null(diag)) stop("compute_betti_curve: NULL diagram (a persistence failure was not propagated)", call. = FALSE)
  
  mat <- as.matrix(diag)
  idx <- mat[, 1] == hom_dim & is.finite(mat[, 2]) & is.finite(mat[, 3])
  if (!any(idx)) return(curve)
  
  birth <- mat[idx, 2]
  death <- mat[idx, 3]
  
  t_seq <- seq(0, maxscale, length.out = n_points)
  alive <- outer(birth, t_seq, "<=") & outer(death, t_seq, ">")
  colSums(alive)
}

diag_to_pi_vec <- function(diag, hom_dim, maxscale, pi_res, pi_sigma) {
  n_pix <- pi_res * pi_res
  if (is.null(diag)) stop("diag_to_pi_vec: NULL diagram (a persistence failure was not propagated)", call. = FALSE)
  
  M <- as.matrix(diag)
  if (ncol(M) < 3) return(rep(0, n_pix))
  
  keep <- M[, 1] == hom_dim & is.finite(M[, 2]) & is.finite(M[, 3]) & (M[, 3] > M[, 2])
  if (!any(keep)) return(rep(0, n_pix))
  
  birth <- M[keep, 2]
  pers  <- M[keep, 3] - M[keep, 2]
  
  birth <- pmin(pmax(birth, 0), maxscale)
  pers  <- pmin(pmax(pers,  0), maxscale)
  
  D <- cbind(rep(hom_dim, length(birth)), birth, pers)
  xSeq <- seq(0, maxscale, length.out = pi_res + 1)
  ySeq <- seq(0, maxscale, length.out = pi_res + 1)
  
  v <- TDAvec::computePersistenceImage(D = D, homDim = hom_dim, xSeq = xSeq, ySeq = ySeq, sigma = pi_sigma)
  v[!is.finite(v)] <- 0
  as.numeric(v)
}

`%||%` <- function(x, y) if (!is.null(x)) x else y

# Build symmetric pair features for a chunk of pairs, vectorized.
make_pair_features_block <- function(feats, i_idx, j_idx, ops) {
  A <- feats[i_idx, , drop = FALSE]
  B <- feats[j_idx, , drop = FALSE]
  
  blocks <- list()
  if ("diff" %in% ops)    blocks <- c(blocks, list(abs(A - B)))
  if ("min" %in% ops)     blocks <- c(blocks, list(pmin(A, B)))
  if ("max" %in% ops)     blocks <- c(blocks, list(pmax(A, B)))
  if ("product" %in% ops) blocks <- c(blocks, list(A * B))
  
  do.call(cbind, blocks)
}

# ============================================================
# METRICS (NO LEAKAGE FOR TAR@FAR)
# ============================================================

roc_metrics <- function(y_true, y_score) {
  if (length(unique(y_true)) < 2) return(list(auc = NA, eer = NA, roc = NULL))
  roc_obj <- pROC::roc(response = y_true, predictor = as.numeric(y_score), quiet = TRUE, direction = "<")
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  sens <- roc_obj$sensitivities
  fpr  <- 1 - roc_obj$specificities
  fnr  <- 1 - sens
  eer_idx <- which.min(abs(fpr - fnr))
  eer <- (fpr[eer_idx] + fnr[eer_idx]) / 2
  
  list(auc = auc_val, eer = eer, roc = roc_obj)
}

tar_at_far_train_threshold <- function(y_train, s_train, y_test, s_test, far) {
  imp_train <- s_train[y_train == 0]
  gen_test  <- s_test[y_test == 1]
  thr <- as.numeric(quantile(imp_train, probs = 1 - far, na.rm = TRUE))
  tar <- mean(gen_test >= thr, na.rm = TRUE)
  list(TAR = tar, threshold = thr)
}

robust_z_fit_train_impostors <- function(scores_train, y_train) {
  imp <- scores_train[y_train == 0]
  med <- median(imp, na.rm = TRUE)
  sc  <- mad(imp, constant = 1, na.rm = TRUE)
  if (!is.finite(sc) || sc < 1e-12) sc <- sd(imp, na.rm = TRUE)
  if (!is.finite(sc) || sc < 1e-12) sc <- 1
  list(med = med, sc = sc)
}

robust_z_apply <- function(scores, fit) {
  z <- (scores - fit$med) / fit$sc
  z[!is.finite(z)] <- 0
  z
}

# ============================================================
# DATA LOADING + CV SPLITS
# ============================================================

get_paths <- function(config) {
  n_id_use <- if (isTRUE(config$debug_fast)) config$debug_n_id else config$n_id
  list.files(config$binary_dir, pattern = "\\.png$", full.names = TRUE, ignore.case = TRUE) %>%
    keep(~{
      id <- get_identity(.x)
      !is.na(id) && as.integer(id) >= 1 && as.integer(id) <= n_id_use
    }) %>%
    sort()
}

make_identity_folds <- function(observed_ids, n_folds, seed) {
  set.seed(seed)
  uniq <- sort(unique(observed_ids))
  uniq <- sample(uniq, length(uniq), replace = FALSE)
  split(uniq, rep(1:n_folds, length.out = length(uniq)))
}

build_pair_mat_global <- function(idx_global) {
  mm <- combn(length(idx_global), 2)
  rbind(idx_global[mm[1, ]], idx_global[mm[2, ]])
}

compute_y_from_pair_mat <- function(pair_mat, ids) {
  as.integer(ids[pair_mat[1, ]] == ids[pair_mat[2, ]])
}

maybe_subsample_pairs <- function(pair_mat, max_pairs, seed) {
  if (is.null(max_pairs) || ncol(pair_mat) <= max_pairs) return(pair_mat)
  set.seed(seed)
  keep <- sample.int(ncol(pair_mat), size = max_pairs, replace = FALSE)
  pair_mat[, keep, drop = FALSE]
}

# ============================================================
# FEATURE BANK (PRECOMPUTE PER IMAGE ONCE)
# ============================================================

pool_curves_topk <- function(curves_mat, top_k) {
  if (is.null(curves_mat) || nrow(curves_mat) == 0) {
    return(list(mean_all = numeric(0), mean_topk = numeric(0), sd_all = numeric(0)))
  }
  nA <- nrow(curves_mat)
  k <- min(top_k, nA)
  
  mean_all <- colMeans(curves_mat)
  anchor_scores <- rowSums(curves_mat)
  top_idx <- order(anchor_scores, decreasing = TRUE)[1:k]
  mean_topk <- colMeans(curves_mat[top_idx, , drop = FALSE])
  
  sd_all <- matrixStats::colSds(curves_mat)
  sd_all[!is.finite(sd_all)] <- 0
  
  list(mean_all = mean_all, mean_topk = mean_topk, sd_all = sd_all)
}

pool_curves_no_topk <- function(curves_mat) {
  if (is.null(curves_mat) || nrow(curves_mat) == 0) {
    return(list(mean_all = numeric(0), sd_all = numeric(0)))
  }
  mean_all <- colMeans(curves_mat)
  sd_all <- matrixStats::colSds(curves_mat)
  sd_all[!is.finite(sd_all)] <- 0
  list(mean_all = mean_all, sd_all = sd_all)
}

extract_feature_bank_one <- function(path, config, disk_mask) {
  bin <- read_binary_fingerprint(path, polarity = config$polarity %||% "auto")
  ridge_area <- mean(bin)
  norm_factor <- max(ridge_area, 0.01)
  
  dt_ridge <- compute_dt_to_ridge(bin)
  dt_valley <- compute_dt_to_valley(bin)
  
  # DT baseline vector: fixed size
  sz <- config$baselines$dt_downsample_size
  ridge_ds <- downsample_matrix(dt_ridge, sz)
  valley_ds <- downsample_matrix(dt_valley, sz)
  ridge_fixed <- pad_or_crop_to_size(ridge_ds, sz)
  valley_fixed <- pad_or_crop_to_size(valley_ds, sz)
  dt_vec <- c(as.numeric(ridge_fixed), as.numeric(valley_fixed))
  dt_vec[!is.finite(dt_vec)] <- 0
  rms <- sqrt(mean(dt_vec^2))
  dt_vec <- dt_vec / if (is.finite(rms) && rms > 1e-12) rms else 1
  
  # DT baseline vector at matched PH resolution (150x150)
  ridge_ds150  <- downsample_matrix(dt_ridge,  config$target_size)
  valley_ds150 <- downsample_matrix(dt_valley, config$target_size)
  ridge_fixed_150  <- pad_or_crop_to_size(ridge_ds150,  config$target_size)
  valley_fixed_150 <- pad_or_crop_to_size(valley_ds150, config$target_size)
  dt_vec_150 <- c(as.numeric(ridge_fixed_150), as.numeric(valley_fixed_150))
  dt_vec_150[!is.finite(dt_vec_150)] <- 0
  rms150 <- sqrt(mean(dt_vec_150^2))
  dt_vec_150 <- dt_vec_150 / if (is.finite(rms150) && rms150 > 1e-12) rms150 else 1

  # Matched-transformation geometry baseline (audit, 31 Oct 2026): the SAME 150x150 fields, but
  # clamped to [0, maxscale], non-finite -> maxscale, rounded to 0.1 BEFORE the RMS normalisation —
  # exactly what compute_persistence() does to the fields the topological methods see — so that
  # DT_L2_150_clamped and the topological methods start from identical values. dt_vec_150 (the
  # original baseline, unclamped) is kept unchanged.
  clamp_round <- function(m, T) { x <- pmin(pmax(m, 0), T); x[!is.finite(x)] <- T; round(x, 1) }
  dt_vec_150_clamped <- c(as.numeric(clamp_round(ridge_fixed_150, config$maxscale)), as.numeric(clamp_round(valley_fixed_150, config$maxscale)))
  rms150c <- sqrt(mean(dt_vec_150_clamped^2))
  dt_vec_150_clamped <- dt_vec_150_clamped / if (is.finite(rms150c) && rms150c > 1e-12) rms150c else 1

  # Common downsample for PH computations (reuse the 150x150 DTs)
  dt_ridge_dsT  <- ridge_ds150
  dt_valley_dsT <- valley_ds150
  bin_dsT       <- downsample_matrix(bin,       config$target_size)
  
  # Global persistence for Betti & PI
  diag_ridge_global  <- compute_persistence(dt_ridge_dsT,  maxdim = 1, maxscale = config$maxscale, complex_type = config$complex_type %||% "triangulated")
  diag_valley_global <- compute_persistence(dt_valley_dsT, maxdim = 1, maxscale = config$maxscale, complex_type = config$complex_type %||% "triangulated")
  
  # Global Betti (100d)
  betti_ridge  <- compute_betti_curve(diag_ridge_global,  1, config$v23$n_points, config$maxscale)
  betti_valley <- compute_betti_curve(diag_valley_global, 1, config$v23$n_points, config$maxscale)
  betti_global_100 <- c(betti_ridge, betti_valley) / norm_factor
  betti_global_100[!is.finite(betti_global_100)] <- 0
  
  # Global PI baseline (separate params)
  gpi_res <- config$baselines$global_pi_res
  gpi_sig <- config$baselines$global_pi_sigma
  pi_ridge_global  <- diag_to_pi_vec(diag_ridge_global,  1, config$maxscale, gpi_res, gpi_sig)
  pi_valley_global <- diag_to_pi_vec(diag_valley_global, 1, config$maxscale, gpi_res, gpi_sig)
  pi_global <- c(pi_ridge_global, pi_valley_global)
  nn <- sqrt(sum(pi_global^2))
  if (!is.finite(nn) || nn < 1e-12) nn <- 1
  pi_global <- pi_global / nn
  
  # Anchors (shared for local V23 + OT)
  seed <- seed_from_string(basename(path))
  anchors <- sample_anchor_points(bin_dsT, config$v23$n_anchors, seed)
  nA <- nrow(anchors)
  
  # Local V23 curves
  ridge_curves  <- matrix(0, nrow = nA, ncol = config$v23$n_points)
  valley_curves <- matrix(0, nrow = nA, ncol = config$v23$n_points)
  
  # Local OT PIs
  d_pi <- config$ot$pi_res * config$ot$pi_res
  X_ridge  <- matrix(0, nrow = nA, ncol = d_pi)
  X_valley <- matrix(0, nrow = nA, ncol = d_pi)
  
  for (k in seq_len(nA)) {
    cx <- anchors[k, 1]; cy <- anchors[k, 2]
    
    pr <- extract_patch_with_mask(dt_ridge_dsT, cx, cy, config$ot$radius, disk_mask)
    rp <- pr$patch
    if (nrow(rp) >= 5 && ncol(rp) >= 5) {
      rp[!pr$mask] <- config$maxscale
      diag_r <- compute_persistence(rp, maxdim = 1, maxscale = config$maxscale, complex_type = config$complex_type %||% "triangulated")
      ridge_curves[k, ] <- compute_betti_curve(diag_r, 1, config$v23$n_points, config$maxscale)
      X_ridge[k, ] <- diag_to_pi_vec(diag_r, 1, config$maxscale, config$ot$pi_res, config$ot$pi_sigma)
    }
    
    pv <- extract_patch_with_mask(dt_valley_dsT, cx, cy, config$ot$radius, disk_mask)
    vp <- pv$patch
    if (nrow(vp) >= 5 && ncol(vp) >= 5) {
      vp[!pv$mask] <- config$maxscale
      diag_v <- compute_persistence(vp, maxdim = 1, maxscale = config$maxscale, complex_type = config$complex_type %||% "triangulated")
      valley_curves[k, ] <- compute_betti_curve(diag_v, 1, config$v23$n_points, config$maxscale)
      X_valley[k, ] <- diag_to_pi_vec(diag_v, 1, config$maxscale, config$ot$pi_res, config$ot$pi_sigma)
    }
  }
  
  # normalize OT rows (L2)
  if (nA > 0) {
    nr <- sqrt(rowSums(X_ridge^2));  nr[nr < 1e-12 | !is.finite(nr)] <- 1
    nv <- sqrt(rowSums(X_valley^2)); nv[nv < 1e-12 | !is.finite(nv)] <- 1
    X_ridge  <- X_ridge  / nr
    X_valley <- X_valley / nv
  }
  
  # OT coords (centered + scaled)
  Z <- if (nA > 0) cbind(anchors[, 1], anchors[, 2]) else matrix(0, nrow = 0, ncol = 2)
  Z <- as.matrix(scale(Z, center = TRUE, scale = FALSE))
  if (config$ot$geo_scale == "target") Z <- Z / config$target_size else Z <- Z / config$ot$radius
  
  # OT precomputes
  X_both <- if (nA > 0) cbind(X_ridge, X_valley) else matrix(0, nrow = 0, ncol = 0)
  
  n2_ridge <- if (nA > 0) rowSums(X_ridge^2)  else numeric(0)
  n2_valley <- if (nA > 0) rowSums(X_valley^2) else numeric(0)
  n2_both  <- if (nA > 0) rowSums(X_both^2)   else numeric(0)
  z2       <- if (nA > 0) rowSums(Z^2)        else numeric(0)
  
  # V23 local pooled features
  if (nA == 0) {
    local_full_300 <- rep(0, 6 * config$v23$n_points)
    local_noTopK_200 <- rep(0, 4 * config$v23$n_points)
  } else {
    ridge_pool  <- pool_curves_topk(ridge_curves,  config$v23$top_k)
    valley_pool <- pool_curves_topk(valley_curves, config$v23$top_k)
    local_full_300 <- c(
      ridge_pool$mean_all, ridge_pool$mean_topk, ridge_pool$sd_all,
      valley_pool$mean_all, valley_pool$mean_topk, valley_pool$sd_all
    )
    
    ridge_nt  <- pool_curves_no_topk(ridge_curves)
    valley_nt <- pool_curves_no_topk(valley_curves)
    local_noTopK_200 <- c(ridge_nt$mean_all, ridge_nt$sd_all, valley_nt$mean_all, valley_nt$sd_all)
  }
  
  local_full_300 <- local_full_300 / norm_factor
  local_full_300[!is.finite(local_full_300)] <- 0
  local_noTopK_200 <- local_noTopK_200 / norm_factor
  local_noTopK_200[!is.finite(local_noTopK_200)] <- 0
  
  list(
    dt_vec = dt_vec,
    dt_vec_150 = dt_vec_150,
    dt_vec_150_clamped = dt_vec_150_clamped,
    betti_global_100 = betti_global_100,
    pi_global = pi_global,
    v23_global_100 = betti_global_100,
    v23_local_full_300 = local_full_300,
    v23_local_noTopK_200 = local_noTopK_200,
    
    # OT fields
    ot_X_ridge = X_ridge,
    ot_X_valley = X_valley,
    ot_Z = Z,
    ot_X_both = X_both,
    ot_n2_ridge = n2_ridge,
    ot_n2_valley = n2_valley,
    ot_n2_both = n2_both,
    ot_z2 = z2
  )
}

precompute_feature_bank <- function(paths, config) {
  future::plan(future::multisession, workers = config$n_workers)
  progressr::handlers(progressr::handler_txtprogressbar(style = 3))
  options(future.globals.maxSize = config$future_globals_gb * 1024^3)
  
  disk_mask <- make_disk_mask(config$ot$radius)
  
  cat("============================================\n")
  cat("PRECOMPUTE FEATURE BANK\n")
  cat("============================================\n\n")
  cat("Images:", length(paths), "\n")
  cat("Workers:", config$n_workers, "\n\n")
  
  t0 <- Sys.time()
  bank <- progressr::with_progress({
    p <- progressr::progressor(steps = length(paths))
    future.apply::future_lapply(paths, function(fp) {
      # A failed image is recorded WITH its error message (revision step 1.5; previously NULL),
      # so the stop below can say why. Successful entries are unchanged.
      out <- tryCatch(extract_feature_bank_one(fp, config, disk_mask),
                      error = function(e) structure(list(error = conditionMessage(e)), class = "fp_precompute_error"))
      p()
      out
    }, future.seed = TRUE)
  })
  cat("Done. Time:", round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n\n")

  failed <- map_lgl(bank, function(b) is.null(b) || inherits(b, "fp_precompute_error"))
  if (any(failed)) {
    bad  <- which(failed)
    msgs <- vapply(bank[head(bad, 3)], function(b) if (inherits(b, "fp_precompute_error")) b$error else "NULL", character(1))
    stop("Some images failed precompute. First bad indices: ", paste(head(bad, 5), collapse = ", "),
         " (", paste(basename(paths[head(bad, 3)]), collapse = ", "), "); first error(s): ", paste(unique(msgs), collapse = " | "))
  }
  
  bank
}

# Cache key: hash config params that affect feature extraction
bank_cache_key <- function(config) {
  key_parts <- list(
    target_size      = config$target_size,
    maxscale         = config$maxscale,
    dt_ds_size       = config$baselines$dt_downsample_size,
    gpi_res          = config$baselines$global_pi_res,
    gpi_sigma        = config$baselines$global_pi_sigma,
    ot_n_anchors     = config$ot$n_anchors,
    ot_radius        = config$ot$radius,
    ot_pi_res        = config$ot$pi_res,
    ot_pi_sigma      = config$ot$pi_sigma,
    ot_use_dim       = config$ot$use_dim,
    v23_n_points     = config$v23$n_points,
    v23_n_anchors    = config$v23$n_anchors,
    v23_radius       = config$v23$radius,
    v23_top_k        = config$v23$top_k,
    binary_dir       = normalizePath(config$binary_dir, mustWork = FALSE),
    n_id             = config$n_id,
    debug_fast       = config$debug_fast %||% FALSE,
    debug_n_id       = config$debug_n_id %||% NA,
    # Hash the per-image seed function so any change to the seed scheme
    # invalidates the cached feature bank automatically (anchors depend on it).
    seed_fn_hash     = digest::digest(deparse(seed_from_string), algo = "sha256")
  )
  # Revision step 1.5: the complex enters the key ONLY for the non-default construction, so
  # the key (and hence every existing bank and fold cache) of a triangulated run is unchanged.
  ct <- config$complex_type %||% "triangulated"
  if (!identical(ct, "triangulated")) key_parts$complex_type <- ct
  # Revision step 1.4: likewise the polarity rule enters the key only when it is not the default.
  pol <- config$polarity %||% "auto"
  if (!identical(pol, "auto")) key_parts$polarity <- pol
  # Audit patch (31 Oct 2026): ot$geo_scale changes the stored anchor coordinates (ot_Z) and was
  # not in the key, so a changed scale could reuse a stale bank; it enters the key when not default.
  gs <- config$ot$geo_scale %||% "target"
  if (!identical(gs, "target")) key_parts$ot_geo_scale <- gs
  digest::digest(key_parts, algo = "sha256")
}

# Every field a bank entry must carry for the current method set (audit patch, 31 Oct 2026). A bank
# computed before a field existed cannot be reused silently: run_prl_suite_v2() stops if one is missing.
BANK_REQUIRED_FIELDS <- c("dt_vec", "dt_vec_150", "dt_vec_150_clamped", "betti_global_100", "pi_global",
                          "v23_global_100", "v23_local_full_300", "v23_local_noTopK_200",
                          "ot_X_ridge", "ot_X_valley", "ot_Z", "ot_X_both", "ot_n2_ridge", "ot_n2_valley", "ot_n2_both", "ot_z2")
check_bank_fields <- function(bank, where = "") {
  miss <- setdiff(BANK_REQUIRED_FIELDS, names(bank[[1]]))
  if (length(miss))
    stop("feature bank", if (nzchar(where)) paste0(" (", where, ")") else "", " lacks field(s) ", paste(miss, collapse = ", "),
         ": it predates the current method set. Use a fresh output directory so the bank is recomputed.", call. = FALSE)
  invisible(TRUE)
}

# Fold-score cache key: hashes the feature-bank key AND the training/scoring/
# evaluation parameters that determine per-fold method results. Without this,
# changing a training hyperparameter (e.g., v23$batch_size) while reusing an
# old feature bank would silently load stale fold scores from cache.
fold_cache_key <- function(config) {
  digest::digest(list(
    fold_cache_version = config$fold_cache_version %||% 1L,
    feature_bank   = bank_cache_key(config),
    cv_seed        = config$cv_seed,
    n_folds        = config$n_folds,
    fars_report    = config$fars_report,
    v23            = config$v23,
    ot             = config$ot,
    ot_variants    = config$ot_variants,
    v23_variants   = config$v23_variants,
    fusion_w_grid  = config$fusion_w_grid
  ), algo = "sha256")
}

load_or_compute_bank <- function(paths, config) {
  cache_dir <- file.path(config$output_dir, "cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  key <- bank_cache_key(config)
  cache_file <- file.path(cache_dir, paste0("bank_", substr(key, 1, 12), ".rds"))

  # Input sidecar (revision step 1.2 integration): the ordered basenames and md5
  # sums of the images the cached bank was computed from. A cached bank is reused
  # only if this sidecar exists and matches the current images exactly (the old
  # check compared image COUNTS only). The sidecar is named so that it does NOT
  # match the bank discovery pattern ^bank_.*\.rds$ used by verify_reproduction.R
  # and validate_aligned_baseline.R.
  inputs_file <- file.path(cache_dir, paste0("inputs_bank_", substr(key, 1, 12), ".rds"))
  inputs_now  <- list(files = basename(paths), md5 = unname(tools::md5sum(paths)))

  # Per-fold score caches in cache/folds were computed from whatever bank was in force
  # then; any recompute of the bank that is not shown to be identical must not proceed
  # while such caches exist (they would be reused by configuration key).
  legacy_fold_caches <- list.files(file.path(cache_dir, "folds"), pattern = "\\.rds$")
  stop_if_fold_caches <- function(why) {
    if (length(legacy_fold_caches) > 0)
      stop(why, " ", length(legacy_fold_caches), " per-fold caches exist in ", file.path(cache_dir, "folds"),
           " and would be reused with a bank they were not computed from. Use a fresh output directory ",
           "for new inputs, or investigate before deleting anything.", call. = FALSE)
  }

  # Guard (revision step 1.1, QA follow-up to 1.2): per-fold caches whose bank file is
  # absent. The bank would be recomputed and the fold caches reused by configuration key,
  # without any check that they were computed from this bank. An output directory holds
  # one run; a missing bank next to fold caches means that run is no longer intact.
  if (!file.exists(cache_file))
    stop_if_fold_caches(paste0("No cached feature bank for this configuration (", basename(cache_file), ") exists, yet"))

  old_bank <- NULL   # a sidecar-less (pre-1.2) cached bank, kept for comparison after recompute
  if (file.exists(cache_file)) {
    cat("============================================\n")
    cat("LOADING CACHED FEATURE BANK\n")
    cat("  ", cache_file, "\n")
    cat("============================================\n\n")
    if (file.exists(inputs_file)) {
      if (identical(readRDS(inputs_file), inputs_now)) {
        bank <- readRDS(cache_file)
        if (length(bank) == length(paths)) {
          cat("Cache valid: input manifest matches (", length(bank), "images ).\n\n")
          return(bank)
        }
        stop_if_fold_caches("Cached bank size mismatch despite a matching input manifest (corrupt cache?).")
        cat("Cache size mismatch (", length(bank), "vs", length(paths), ") despite a matching input manifest. Recomputing.\n\n")
      } else {
        stop_if_fold_caches("Cached bank was computed from DIFFERENT images (input manifest mismatch).")
        cat("Cached bank was computed from DIFFERENT images (input manifest mismatch). Recomputing.\n\n")
      }
    } else {
      cat("Cached bank has no input manifest (pre-1.2 cache). Recomputing and comparing with the cached bank.\n\n")
      old_bank <- tryCatch(readRDS(cache_file), error = function(e) NULL)
      if (is.null(old_bank)) stop_if_fold_caches("Sidecar-less cached bank could not be read, so it cannot be compared with a recompute.")
    }
  }

  bank <- precompute_feature_bank(paths, config)

  if (!is.null(old_bank)) {
    same <- identical(bank, old_bank)
    cat("Recomputed bank identical to the sidecar-less cached bank:", same, "\n")
    if (!same) {
      # Per-field report, then stop if anything downstream could reuse the old bank.
      fields <- union(names(bank[[1]]), names(old_bank[[1]]))
      for (fn in fields) {
        d <- tryCatch(max(vapply(seq_along(bank), function(i) {
          a <- bank[[i]][[fn]]; b <- old_bank[[i]][[fn]]
          if (is.null(a) || is.null(b) || length(a) != length(b)) Inf else max(abs(as.numeric(a) - as.numeric(b)))
        }, numeric(1))), error = function(e) NA_real_)
        cat(sprintf("   field %-22s max |new - cached| = %s\n", fn, format(d)))
      }
      stop_if_fold_caches("Recomputed feature bank differs from the sidecar-less cached bank (differences printed above).")
      cat("No per-fold caches present; proceeding with the recomputed bank.\n\n")
    }
  }

  cat("Saving feature bank cache to", cache_file, "\n")
  saveRDS(bank, cache_file)
  saveRDS(inputs_now, inputs_file)
  cat("Cache and input manifest saved.\n\n")

  bank
}

bank_to_matrix <- function(bank, field_name) {
  do.call(rbind, lapply(bank, `[[`, field_name))
}

get_v23_feature_matrix <- function(bank, variant) {
  global_100 <- bank_to_matrix(bank, "v23_global_100")
  local_300  <- bank_to_matrix(bank, "v23_local_full_300")
  local_200  <- bank_to_matrix(bank, "v23_local_noTopK_200")
  
  if (variant == "V23_full") return(cbind(global_100, local_300))
  if (variant == "V23_noTopK") return(cbind(global_100, local_200))
  if (variant == "V23_globalOnly") return(global_100)
  if (variant == "V23_localOnly") return(local_300)
  stop("Unknown V23 variant: ", variant)
}

# ============================================================
# SCORING: BASELINES (VECTORIZED - Speedup #1)
# ============================================================

score_pairs_l2 <- function(pair_mat, feats_mat, chunk_size = 10000L) {
  n_pairs <- ncol(pair_mat)
  scores <- numeric(n_pairs)
  
  starts <- seq(1L, n_pairs, by = chunk_size)
  ends <- pmin(starts + chunk_size - 1L, n_pairs)
  
  
  for (ci in seq_along(starts)) {
    s <- starts[ci]; e <- ends[ci]
    i_idx <- pair_mat[1, s:e]
    j_idx <- pair_mat[2, s:e]
    
    A <- feats_mat[i_idx, , drop = FALSE]
    B <- feats_mat[j_idx, , drop = FALSE]
    
    # Vectorized L2 distance
    diff <- A - B
    scores[s:e] <- -sqrt(rowSums(diff * diff))
  }
  
  scores
}

score_pairs_l1 <- function(pair_mat, feats_mat, chunk_size = 10000L) {
  n_pairs <- ncol(pair_mat)
  scores <- numeric(n_pairs)
  
  starts <- seq(1L, n_pairs, by = chunk_size)
  ends <- pmin(starts + chunk_size - 1L, n_pairs)
  
  for (ci in seq_along(starts)) {
    s <- starts[ci]; e <- ends[ci]
    i_idx <- pair_mat[1, s:e]
    j_idx <- pair_mat[2, s:e]
    
    A <- feats_mat[i_idx, , drop = FALSE]
    B <- feats_mat[j_idx, , drop = FALSE]
    
    # Vectorized L1 distance
    scores[s:e] <- -rowSums(abs(A - B))
  }
  
  scores
}

# ============================================================
# OT SCORING (Speedups #3, #4, #5)
# ============================================================

sinkhorn_cost_variant <- function(A, B,
                                  reg, niter, tol, check_every,
                                  dustbin, dustbin_cost, eps_dustbin,
                                  beta_geo,
                                  use_ridge, use_valley) {
  AZ <- A$ot_Z; BZ <- B$ot_Z
  n0 <- nrow(AZ); m0 <- nrow(BZ)
  if (n0 == 0 || m0 == 0) return(Inf)
  
  
  if (use_ridge && use_valley) {
    AX <- A$ot_X_both;  BX <- B$ot_X_both
    aa <- A$ot_n2_both; bb <- B$ot_n2_both
  } else if (use_ridge) {
    AX <- A$ot_X_ridge; BX <- B$ot_X_ridge
    aa <- A$ot_n2_ridge; bb <- B$ot_n2_ridge
  } else if (use_valley) {
    AX <- A$ot_X_valley; BX <- B$ot_X_valley
    aa <- A$ot_n2_valley; bb <- B$ot_n2_valley
  } else {
    return(Inf)
  }
  
  # Descriptor cost matrix: ||x_i - y_j||^2 = aa + bb - 2 x_i·y_j
  # Speedup #5: Avoid outer() - use rep + matrix instead
  G <- AX %*% t(BX)
  aa_mat <- matrix(aa, nrow = n0, ncol = m0, byrow = FALSE)
  bb_mat <- matrix(bb, nrow = n0, ncol = m0, byrow = TRUE)
  Cx <- aa_mat + bb_mat - 2 * G
  Cx <- pmax(Cx, 0)
  
  # Geometry cost matrix
  za <- A$ot_z2
  zb <- B$ot_z2
  Hz <- AZ %*% t(BZ)
  za_mat <- matrix(za, nrow = n0, ncol = m0, byrow = FALSE)
  zb_mat <- matrix(zb, nrow = n0, ncol = m0, byrow = TRUE)
  Cz <- za_mat + zb_mat - 2 * Hz
  Cz <- pmax(Cz, 0)
  
  C <- Cx + beta_geo * Cz
  
  if (dustbin) {
    C_aug <- matrix(dustbin_cost, nrow = n0 + 1, ncol = m0 + 1)
    C_aug[1:n0, 1:m0] <- C
    C_aug[n0 + 1, m0 + 1] <- 0
    C <- C_aug
    
    a <- c(rep((1 - eps_dustbin) / n0, n0), eps_dustbin)
    b <- c(rep((1 - eps_dustbin) / m0, m0), eps_dustbin)
    n <- n0 + 1
    m <- m0 + 1
  } else {
    a <- rep(1 / n0, n0)
    b <- rep(1 / m0, m0)
    n <- n0
    m <- m0
  }
  
  # Sinkhorn kernel
  K <- exp(-C / reg)
  K[!is.finite(K)] <- 0
  if (sum(K) == 0) return(Inf)
  
  # Precompute C * K for final cost (Speedup #3 preparation)
  CK <- C * K
  
  u <- rep(1, n)
  v <- rep(1, m)
  
  # Sinkhorn iterations (with optional early stopping when tol > 0)
  use_early_stop <- tol > 0
  
  for (it in seq_len(niter)) {
    if (use_early_stop) u_old <- u
    
    Kv <- K %*% v
    Kv[Kv < 1e-300] <- 1e-300
    u <- a / as.numeric(Kv)
    
    KTu <- crossprod(K, u)
    KTu[KTu < 1e-300] <- 1e-300
    v <- b / as.numeric(KTu)
    
    # Check convergence every check_every iterations (only if tol > 0)
    if (use_early_stop && (it %% check_every == 0)) {
      rel_change <- max(abs(u - u_old)) / (max(abs(u_old)) + 1e-10)
      if (rel_change < tol) break
    }
  }
  
  # Speedup #3: Compute cost without allocating P
  # sum(P * C) = sum((u %o% v) * K * C) = sum(u * (CK %*% v))
  cost <- sum(u * (CK %*% v))
  cost
}

score_pairs_ot <- function(pair_mat, bank, config, variant_row) {
  n_pairs <- ncol(pair_mat)
  
  # Speedup #2: Use larger chunk size
  chunk_size <- config$ot$chunk_size_pairs %||% 5000L
  chunk_size <- max(500L, chunk_size)
  
  starts <- seq(1L, n_pairs, by = chunk_size)
  ends <- pmin(starts + chunk_size - 1L, n_pairs)
  
  # Extract config values once
  reg <- config$ot$sinkhorn_reg
  niter <- config$ot$sinkhorn_iter
  tol <- config$ot$sinkhorn_tol %||% 0
  check_every <- config$ot$sinkhorn_check_every %||% 10L
  dustbin_cost <- config$ot$dustbin_cost
  eps_dustbin <- config$ot$eps_dustbin
  
  use_ridge <- variant_row$use_ridge
  use_valley <- variant_row$use_valley
  beta_geo <- variant_row$beta_geo
  dustbin <- variant_row$dustbin
  
  res_list <- progressr::with_progress({
    p <- progressr::progressor(steps = length(starts))
    future.apply::future_lapply(seq_along(starts), function(ci) {
      s <- starts[ci]; e <- ends[ci]
      out <- numeric(e - s + 1L)
      
      for (k in seq_along(out)) {
        idx <- s + k - 1L
        i <- pair_mat[1, idx]; j <- pair_mat[2, idx]
        cost <- sinkhorn_cost_variant(
          bank[[i]], bank[[j]],
          reg = reg,
          niter = niter,
          tol = tol,
          check_every = check_every,
          dustbin = dustbin,
          dustbin_cost = dustbin_cost,
          eps_dustbin = eps_dustbin,
          beta_geo = beta_geo,
          use_ridge = use_ridge,
          use_valley = use_valley
        )
        out[k] <- -cost
      }
      
      p()
      list(s = s, e = e, out = out)
    }, future.seed = TRUE)
  })
  
  scores <- numeric(n_pairs)
  for (r in res_list) scores[r$s:r$e] <- r$out
  scores
}

# ============================================================
# V23 PAIR LOGISTIC (TRAIN ON TRAIN PAIRS ONLY)
# ============================================================

compute_symmetric_features <- function(v1, v2, ops) {
  out <- list()
  if ("diff" %in% ops)    out$diff    <- abs(v1 - v2)
  if ("min" %in% ops)     out$min     <- pmin(v1, v2)
  if ("max" %in% ops)     out$max     <- pmax(v1, v2)
  if ("product" %in% ops) out$product <- v1 * v2
  do.call(c, out)
}

# Speedup #6: Reduced default sample size (50k -> 20k in CONFIG)
estimate_pair_standardizer <- function(pair_mat, feats, ops, n_sample, seed, chunk_pairs = 2000L) {
  set.seed(seed)
  n_pairs <- ncol(pair_mat)
  idx_all <- sample.int(n_pairs, size = min(n_sample, n_pairs), replace = FALSE)
  
  d <- ncol(feats)
  p <- d * length(ops)
  
  sum_x  <- rep(0, p)
  sum_x2 <- rep(0, p)
  n_seen <- 0L
  
  starts <- seq(1L, length(idx_all), by = chunk_pairs)
  ends <- pmin(starts + chunk_pairs - 1L, length(idx_all))
  
  for (ci in seq_along(starts)) {
    s <- starts[ci]; e <- ends[ci]
    idx <- idx_all[s:e]
    i_idx <- pair_mat[1, idx]
    j_idx <- pair_mat[2, idx]
    
    X <- make_pair_features_block(feats, i_idx, j_idx, ops)
    
    sum_x  <- sum_x  + colSums(X)
    sum_x2 <- sum_x2 + colSums(X * X)
    n_seen <- n_seen + nrow(X)
  }
  
  mu <- sum_x / n_seen
  var <- (sum_x2 - n_seen * mu * mu) / max(1, (n_seen - 1))
  sdv <- sqrt(pmax(var, 0))
  sdv[!is.finite(sdv) | sdv < 1e-10] <- 1
  
  list(mu = mu, sd = sdv)
}

fit_logistic_balanced_pairs <- function(pair_mat, y, feats, ops, standardizer,
                                        lambda, max_iter, tol, batch_size, lr0,
                                        seed, verbose = FALSE) {
  set.seed(seed)
  
  mu <- standardizer$mu
  sdv <- standardizer$sd
  p <- length(mu)
  
  pos_idx <- which(y == 1)
  neg_idx <- which(y == 0)
  
  n_pos <- length(pos_idx)
  n_neg <- length(neg_idx)
  
  # Balanced mini-batches with replacement, per the manuscript's TopoLR spec:
  # 5,000 genuine + 5,000 impostor pairs drawn independently with replacement.
  # The previous implementation used min(half_batch, n_pos), which silently
  # collapsed the genuine half to whatever was available (often all 2,240
  # training-fold genuine pairs sampled without replacement), violating the
  # stated training procedure.
  half_batch  <- batch_size %/% 2L
  n_pos_batch <- half_batch
  n_neg_batch <- batch_size - half_batch

  w <- rep(0, p)
  b <- 0

  lr <- lr0
  prev_loss <- Inf
  use_early_stop <- is.finite(tol) && tol > 0

  build_X_std <- function(pair_indices) {
    i_idx <- pair_mat[1, pair_indices]
    j_idx <- pair_mat[2, pair_indices]
    X <- make_pair_features_block(feats, i_idx, j_idx, ops)
    sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
  }

  for (iter in 1:max_iter) {
    pos_sample <- pos_idx[sample.int(n_pos, n_pos_batch, replace = TRUE)]
    neg_sample <- neg_idx[sample.int(n_neg, n_neg_batch, replace = TRUE)]
    batch_idx <- c(pos_sample, neg_sample)
    
    Xb <- build_X_std(batch_idx)
    yb <- y[batch_idx]
    bsz <- length(batch_idx)
    
    z <- Xb %*% w + b
    z <- pmin(pmax(z, -500), 500)
    pred <- 1 / (1 + exp(-z))
    pred <- pmax(pmin(pred, 1 - 1e-10), 1e-10)
    
    error <- as.numeric(pred) - yb
    grad_w <- (t(Xb) %*% error) / bsz + lambda * w
    grad_b <- mean(error)
    
    w <- w - lr * as.numeric(grad_w)
    b <- b - lr * grad_b
    
    if (iter %% 20 == 0) {
      eval_pos <- pos_idx[sample.int(n_pos, n_pos_batch, replace = TRUE)]
      eval_neg <- neg_idx[sample.int(n_neg, n_neg_batch, replace = TRUE)]
      eval_idx <- c(eval_pos, eval_neg)

      Xe <- build_X_std(eval_idx)
      ye <- y[eval_idx]

      ze <- Xe %*% w + b
      ze <- pmin(pmax(ze, -500), 500)
      pre <- 1 / (1 + exp(-ze))
      pre <- pmax(pmin(pre, 1 - 1e-10), 1e-10)

      loss <- -mean(ye * log(pre) + (1 - ye) * log(1 - pre)) + 0.5 * lambda * sum(w^2)

      if (verbose) cat(sprintf("      iter=%d  loss=%.4f  lr=%.4f\n", iter, loss, lr))
      if (use_early_stop && abs(prev_loss - loss) < tol) break
      prev_loss <- loss
      lr <- lr * 0.9
    }
  }
  
  list(w = as.vector(w), b = b, standardizer = standardizer, ops = ops)
}

predict_logistic_pairs <- function(pair_mat, feats, model, chunk_size = 8000L) {
  n_pairs <- ncol(pair_mat)
  
  w <- model$w
  b <- model$b
  mu <- model$standardizer$mu
  sdv <- model$standardizer$sd
  ops <- model$ops
  
  scores <- numeric(n_pairs)
  
  starts <- seq(1L, n_pairs, by = chunk_size)
  ends <- pmin(starts + chunk_size - 1L, n_pairs)
  
  for (ci in seq_along(starts)) {
    s <- starts[ci]; e <- ends[ci]
    idx <- s:e
    
    i_idx <- pair_mat[1, idx]
    j_idx <- pair_mat[2, idx]
    
    X <- make_pair_features_block(feats, i_idx, j_idx, ops)
    X <- sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
    
    z <- X %*% w + b
    z <- pmin(pmax(z, -500), 500)
    scores[s:e] <- as.numeric(1 / (1 + exp(-z)))
  }
  
  scores
}

# ============================================================
# FUSION (TUNE w ON TRAIN ONLY; APPLY TO TEST)
# ============================================================

tune_fusion_w_on_train <- function(y_train, s_train_ot, s_train_v, w_grid) {
  zfit_ot <- robust_z_fit_train_impostors(s_train_ot, y_train)
  zfit_v  <- robust_z_fit_train_impostors(s_train_v,  y_train)
  
  z_ot <- robust_z_apply(s_train_ot, zfit_ot)
  z_v  <- robust_z_apply(s_train_v,  zfit_v)
  
  best_w <- NA_real_
  best_tar <- -Inf
  
  for (w in w_grid) {
    s <- w * z_v + (1 - w) * z_ot
    imp <- s[y_train == 0]
    gen <- s[y_train == 1]
    thr <- as.numeric(quantile(imp, probs = 1 - 1e-3, na.rm = TRUE))
    tar <- mean(gen >= thr, na.rm = TRUE)
    if (tar > best_tar) {
      best_tar <- tar
      best_w <- w
    }
  }
  
  list(best_w = best_w, zfit_ot = zfit_ot, zfit_v = zfit_v)
}

apply_fusion <- function(s_ot, s_v, w, zfit_ot, zfit_v) {
  z_ot <- robust_z_apply(s_ot, zfit_ot)
  z_v  <- robust_z_apply(s_v,  zfit_v)
  w * z_v + (1 - w) * z_ot
}

# ============================================================
# PER-FOLD EVALUATION CORE
# ============================================================

eval_fold_method <- function(method_name, y_train, s_train, y_test, s_test, fars) {
  mm <- roc_metrics(y_test, s_test)
  
  tar_tbl <- tibble(FAR = fars, TAR = NA_real_)
  for (ii in seq_along(fars)) {
    out <- tar_at_far_train_threshold(y_train, s_train, y_test, s_test, fars[ii])
    tar_tbl$TAR[ii] <- out$TAR
  }
  
  list(method = method_name, auc = mm$auc, eer = mm$eer, tar_tbl = tar_tbl)
}

# ============================================================
# MAIN RUNNER
# ============================================================

run_prl_suite_v2 <- function(config = CONFIG) {
  dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create("Manuscript", showWarnings = FALSE, recursive = TRUE)
  
  paths <- get_paths(config)
  stopifnot(length(paths) > 0)
  ids <- map_chr(paths, get_identity)
  stopifnot(all(!is.na(ids)))
  
  folds <- make_identity_folds(observed_ids = ids, n_folds = config$n_folds, seed = config$cv_seed)

  # ---- Run manifest: records exactly which images, identities, fold memberships
  # and pair orderings this run used, so outputs can be compared positionally.
  # Pair matrices are recorded BEFORE any debug_fast subsampling.
  dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)
  manifest_folds <- lapply(seq_along(folds), function(f) {
    test_ids  <- folds[[f]]
    train_ids <- setdiff(unique(ids), test_ids)
    train_idx <- which(ids %in% train_ids)
    test_idx  <- which(ids %in% test_ids)
    list(fold = f, test_ids = test_ids, train_ids = train_ids,
         train_idx = train_idx, test_idx = test_idx,
         pair_train = build_pair_mat_global(train_idx),
         pair_test  = build_pair_mat_global(test_idx))
  })
  run_manifest <- list(
    created      = Sys.time(),
    paths        = paths,
    files        = basename(paths),
    md5          = unname(tools::md5sum(paths)),
    ids          = ids,
    folds        = manifest_folds,
    pairs_are_pre_subsampling = TRUE,
    acceptance_rule = "score >= threshold; threshold = quantile(train impostor scores, 1 - FAR), type 7",
    config       = config,
    locale       = Sys.getlocale(),
    session_info = utils::sessionInfo()
  )

  # ---- Input-change guard (revision step 1.2 integration). Per-fold score caches are
  # keyed by CONFIGURATION (fold_cache_key) but not by image identity, so if the images,
  # their order, the identities or the fold/pair assignments changed since the previous
  # run in this output directory, reusing the existing caches would silently mix inputs.
  # Compare with the previous manifest BEFORE overwriting it; stop on mismatch whenever
  # legacy per-fold caches are present.
  manifest_file      <- file.path(config$output_dir, "run_manifest.rds")
  legacy_fold_caches <- list.files(file.path(config$output_dir, "cache", "folds"), pattern = "\\.rds$")
  if (file.exists(manifest_file)) {
    prev <- readRDS(manifest_file)
    same_inputs <- identical(prev$files, run_manifest$files) &&
                   identical(prev$md5,   run_manifest$md5)   &&
                   identical(prev$ids,   run_manifest$ids)   &&
                   identical(prev$folds, run_manifest$folds)
    if (same_inputs) {
      cat("Previous run manifest found (", format(prev$created), "): identical images, md5s, identities, folds and pair orderings.\n",
          "  ", length(legacy_fold_caches), " existing per-fold caches may be reused (by configuration key).\n", sep = "")
    } else if (length(legacy_fold_caches) > 0) {
      stop("Inputs differ from the previous run manifest in '", config$output_dir, "' (files, md5, identities, ",
           "folds or pair orderings), but ", length(legacy_fold_caches), " per-fold caches from that run are present ",
           "and would be reused. Use a fresh output directory for the new inputs (or move the old caches away).",
           call. = FALSE)
    } else {
      cat("Previous run manifest differs from the current inputs; no per-fold caches present, proceeding fresh.\n")
    }
    run_manifest$previous_manifest_created <- prev$created
    run_manifest$inputs_identical_to_previous <- same_inputs
  } else if (length(legacy_fold_caches) > 0) {
    # Guard (revision step 1.1, QA follow-up to 1.2): without a manifest the images these
    # caches were computed from cannot be verified, so they must not be reused.
    stop(length(legacy_fold_caches), " per-fold caches are present in '", file.path(config$output_dir, "cache", "folds"),
         "' but '", config$output_dir, "' has no run manifest, so the images they were computed from cannot be ",
         "verified (they would be reused by configuration key alone). Use a fresh output directory for this run, ",
         "or restore the manifest of the run that produced them.", call. = FALSE)
  }

  saveRDS(run_manifest, manifest_file)
  cat("Run manifest written to", manifest_file, "\n\n")
  
  cat("============================================\n")
  cat("DATA SUMMARY\n")
  cat("============================================\n\n")
  cat("Images:", length(paths), "\n")
  cat("Identities:", length(unique(ids)), "\n")
  cat("Folds:", config$n_folds, "\n")
  cat("Workers:", config$n_workers, "\n")
  cat("OT chunk size:", config$ot$chunk_size_pairs, "\n")
  
  sink_tol <- config$ot$sinkhorn_tol %||% 0
  if (sink_tol > 0) {
    cat("Sinkhorn early stopping: ENABLED (tol =", sink_tol, 
        ", check every", config$ot$sinkhorn_check_every %||% 10, "iters)\n\n")
  } else {
    cat("Sinkhorn early stopping: DISABLED (fixed", config$ot$sinkhorn_iter, "iters)\n\n")
  }
  
  # Ensure the configured parallel plan is active even when the feature bank is
  # loaded from cache, so later OT scoring uses config$n_workers (runtime reproducibility;
  # results are unaffected). precompute_feature_bank() also sets this, harmlessly, on a cache miss.
  future::plan(future::multisession, workers = config$n_workers)
  options(future.globals.maxSize = config$future_globals_gb * 1024^3)

  bank <- load_or_compute_bank(paths, config)

  # Feature matrices
  check_bank_fields(bank, config$output_dir)
  dt_mat <- bank_to_matrix(bank, "dt_vec")
  dt_mat_150 <- bank_to_matrix(bank, "dt_vec_150")
  dt_mat_150_clamped <- bank_to_matrix(bank, "dt_vec_150_clamped")
  betti_mat <- bank_to_matrix(bank, "betti_global_100")
  gpi_mat <- bank_to_matrix(bank, "pi_global")
  
  v23_feat_cache <- setNames(
    lapply(config$v23_variants, function(vv) get_v23_feature_matrix(bank, vv)),
    config$v23_variants
  )

  # ---- Aligned geometry baseline (revision step 1.2): the all-pairs rigid search is a
  # SHARED, training-free computation (~2 h on DB1), done once per (grid, images) identity
  # and cached in cache/aligned_<key12>.rds; the fold loop below only LOOKS UP scores.
  # Its elapsed time is therefore recorded here (aligned_runtime.csv), not in the per-fold
  # runtime rows, which for this method measure lookup only.
  aligned_key_now    <- aligned_cache_key(config$aligned, paths)
  aligned_cache_file <- file.path(config$output_dir, "cache", paste0("aligned_", substr(aligned_key_now, 1, 12), ".rds"))
  aligned_run_start  <- Sys.time()
  aligned <- load_or_compute_aligned(bank, paths, config)
  stopifnot(identical(aligned$key, aligned_key_now), identical(aligned$files, basename(paths)))
  # The engine stamps `created` only when it computes, so this is TRUE exactly when the search
  # ran in this call (including the recompute-on-identity-mismatch path inside the engine).
  aligned_computed_now <- isTRUE(aligned$created >= aligned_run_start)
  aligned_pidx <- pair_index_matrix(length(paths), aligned$pairs)

  write_csv(tibble(
    method = "DT_L2_150_aligned", aligned_key = aligned$key, cache_file = basename(aligned_cache_file),
    computed_this_run = aligned_computed_now, shared_search_minutes = aligned$elapsed_min,
    n_workers_if_computed_this_run = if (aligned_computed_now) as.integer(config$n_workers) else NA_integer_,
    n_images = length(aligned$files), n_pairs = nrow(aligned$pairs), canvas = aligned$canvas,
    theta_max = aligned$config$theta_max, theta_step = aligned$config$theta_step, refine_half = aligned$config$refine_half,
    t_max = aligned$config$t_max, min_overlap_frac = aligned$config$min_overlap_frac, scorer_version = aligned$config$scorer_version,
    created = aligned$created,
    note = "fold_runtimes_seconds.csv rows for this method are score LOOKUP time only"),
    file.path(config$output_dir, "aligned_runtime.csv"))

  # Per-pair diagnostics and the boundary-hit rates used by the pre-specified expansion
  # rule (which is evaluated on GENUINE pairs only; impostor rows are context).
  aligned_pair_type <- ifelse(ids[aligned$pairs$i] == ids[aligned$pairs$j], "genuine", "impostor")
  aligned_diag <- cbind(aligned$diag, pair_type = aligned_pair_type, stringsAsFactors = FALSE)
  write_csv(aligned_diag, file.path(config$output_dir, "aligned_pair_diagnostics.csv"))
  aligned_boundary <- aligned_diag %>%
    group_by(pair_type) %>%
    summarize(n = n(),
              theta_boundary_rate = mean(theta_at_bound),
              t_boundary_rate     = mean(t_at_bound),
              content_min_lt_0.5  = mean(content_min < 0.5),
              median_omega_frac   = median(omega) / config$aligned$frame^2,
              median_abs_theta    = median(abs(theta)),
              median_max_abs_t    = median(pmax(abs(tx), abs(ty))),
              median_d_sym        = median(d_sym),
              median_d_unaligned  = median(d_unaligned),
              frac_d_improved     = mean(d_sym < d_unaligned),
              .groups = "drop") %>%
    mutate(theta_max = config$aligned$theta_max, t_max = config$aligned$t_max, aligned_key = aligned$key)
  write_csv(aligned_boundary, file.path(config$output_dir, "aligned_boundary_rates.csv"))
  cat("Aligned baseline boundary-hit rates (expansion rule applies to genuine pairs only):\n")
  print(as.data.frame(aligned_boundary[, c("pair_type", "n", "theta_boundary_rate", "t_boundary_rate", "content_min_lt_0.5", "median_omega_frac")]), row.names = FALSE)
  cat("\n")

  runtime_rows <- list()
  fold_records <- list()
  plot_scores <- list()
  fold_cache_index <- list()   # (fold, method, cache file, extra key) for every method touched in this run

  # Fold-level score caching helpers
  fold_cache_dir <- file.path(config$output_dir, "cache", "folds")
  dir.create(fold_cache_dir, showWarnings = FALSE, recursive = TRUE)

  # Fold-score cache key includes BOTH feature-bank parameters AND training/
  # evaluation parameters, so changing v23/ot/fusion settings invalidates
  # stale per-fold score caches.
  fold_key <- substr(fold_cache_key(config), 1, 8)

  # Optional method-specific identity (revision step 1.2): a method may pass an
  # `extra_key` whose first 8 characters are appended to its cache file name, so that
  # a change in ITS OWN inputs (e.g. the aligned search grid) invalidates its caches
  # and nobody else's. Only DT_L2_150_aligned uses this (its full pair-cache key).
  fold_cache_path <- function(fold_id, method, extra_key = NULL) {
    suffix <- if (is.null(extra_key)) "" else paste0("_", substr(extra_key, 1, 8))
    file.path(fold_cache_dir, sprintf("fold%d_%s_%s%s.rds", fold_id, method, fold_key, suffix))
  }

  save_fold_cache <- function(fold_id, method, s_train, s_test, y_train, y_test, elapsed, extras = list(), extra_key = NULL) {
    obj <- c(list(s_train = s_train, s_test = s_test,
                  y_train = y_train, y_test = y_test, elapsed = elapsed), extras)
    saveRDS(obj, fold_cache_path(fold_id, method, extra_key))
  }

  load_fold_cache <- function(fold_id, method, extra_key = NULL) {
    fp <- fold_cache_path(fold_id, method, extra_key)
    if (file.exists(fp)) readRDS(fp) else NULL
  }

  record_fold_cache <- function(fold_id, method, extra_key, from_cache) {
    fold_cache_index[[length(fold_cache_index) + 1]] <<- tibble(
      fold = fold_id, method = method,
      file = basename(fold_cache_path(fold_id, method, extra_key)),
      fold_key = fold_key, extra_key = extra_key %||% NA_character_, loaded_from_cache = from_cache)
  }

  # Helper: run or load a method for a fold
  run_or_load_method <- function(fold_id, method, compute_fn, y_train, y_test, extra_key = NULL) {
    cached <- load_fold_cache(fold_id, method, extra_key)
    if (!is.null(cached)) {
      cat(sprintf("  [cached] %s\n", method))
      res <- eval_fold_method(method, cached$y_train, cached$s_train, cached$y_test, cached$s_test, config$fars_report)
      runtime_rows[[length(runtime_rows) + 1]] <<- tibble(fold = fold_id, method = method, seconds = cached$elapsed)
      fold_records[[length(fold_records) + 1]] <<- list(fold = fold_id, res = res)
      plot_scores[[length(plot_scores) + 1]] <<- tibble(fold = fold_id, method = method, y = cached$y_test, score = cached$s_test)
      record_fold_cache(fold_id, method, extra_key, TRUE)
      return(list(s_train = cached$s_train, s_test = cached$s_test))
    }

    t0 <- Sys.time()
    scores <- compute_fn()
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    extras <- scores[setdiff(names(scores), c("s_train", "s_test"))]
    # Training-calibrated thresholds at each reported FAR, for fixed-model reuse.
    # Identical to tar_at_far_train_threshold(); acceptance rule is score >= threshold,
    # threshold = quantile(training impostor scores, 1 - FAR), R default type 7.
    imp_train <- scores$s_train[y_train == 0]
    extras$thresholds <- setNames(
      vapply(config$fars_report,
             function(f) as.numeric(quantile(imp_train, probs = 1 - f, na.rm = TRUE)),
             numeric(1)),
      paste0("far_", config$fars_report))
    save_fold_cache(fold_id, method, scores$s_train, scores$s_test, y_train, y_test, elapsed, extras, extra_key)
    record_fold_cache(fold_id, method, extra_key, FALSE)

    runtime_rows[[length(runtime_rows) + 1]] <<- tibble(fold = fold_id, method = method, seconds = elapsed)
    res <- eval_fold_method(method, y_train, scores$s_train, y_test, scores$s_test, config$fars_report)
    fold_records[[length(fold_records) + 1]] <<- list(fold = fold_id, res = res)
    plot_scores[[length(plot_scores) + 1]] <<- tibble(fold = fold_id, method = method, y = y_test, score = scores$s_test)
    scores
  }

  for (fold_id in seq_along(folds)) {
    test_ids <- folds[[fold_id]]
    train_ids <- setdiff(unique(ids), test_ids)

    train_idx <- which(ids %in% train_ids)
    test_idx  <- which(ids %in% test_ids)

    pair_train <- build_pair_mat_global(train_idx)
    pair_test  <- build_pair_mat_global(test_idx)

    if (isTRUE(config$debug_fast)) {
      pair_train <- maybe_subsample_pairs(pair_train, config$debug_max_pairs_train, seed = 1000 + fold_id)
      pair_test  <- maybe_subsample_pairs(pair_test,  config$debug_max_pairs_test,  seed = 2000 + fold_id)
    }

    y_train <- compute_y_from_pair_mat(pair_train, ids)
    y_test  <- compute_y_from_pair_mat(pair_test,  ids)

    cat("============================================\n")
    cat(sprintf("FOLD %d/%d\n", fold_id, config$n_folds))
    cat("============================================\n\n")
    cat("Train IDs:", length(train_ids), " Train images:", length(train_idx), " Train pairs:", ncol(pair_train), "\n")
    cat("Test  IDs:", length(test_ids),  " Test images:",  length(test_idx),  " Test pairs:",  ncol(pair_test),  "\n")
    cat("Test impostors:", sum(y_test == 0), " (FAR=1e-3 => ~", round(sum(y_test == 0) * 1e-3), ")\n\n", sep = "")

    # ---------------- Baseline: DT_L2 ----------------
    run_or_load_method(fold_id, "DT_L2", function() {
      list(s_train = score_pairs_l2(pair_train, dt_mat),
           s_test  = score_pairs_l2(pair_test,  dt_mat))
    }, y_train, y_test)

    # ---------------- Baseline: DT_L2_150 (matched resolution) ----------------
    run_or_load_method(fold_id, "DT_L2_150", function() {
      list(s_train = score_pairs_l2(pair_train, dt_mat_150),
           s_test  = score_pairs_l2(pair_test,  dt_mat_150))
    }, y_train, y_test)

    # ---------------- Baseline: DT_L2_150_clamped (matched transformation; audit, 31 Oct 2026) ----------------
    # The 150x150 fields clamped to [0, 20] and rounded to 0.1 before RMS normalisation, as the
    # topological methods receive them; otherwise identical to DT_L2_150. Training-free.
    run_or_load_method(fold_id, "DT_L2_150_clamped", function() {
      list(s_train = score_pairs_l2(pair_train, dt_mat_150_clamped),
           s_test  = score_pairs_l2(pair_test,  dt_mat_150_clamped))
    }, y_train, y_test)

    # ---------------- Baseline: DT_L2_150_aligned (rigid search; revision step 1.2) ----------------
    # Same inputs as DT_L2_150; score = -d_sym from the shared pair cache computed above
    # (coarse angular search with local refinement and exhaustive translation search).
    # Training-free; thresholds come from training-impostor quantiles like every method.
    # The cache file name carries the pair-cache key, and the cache itself records the
    # full key and configuration (large diagnostics stay in the pair cache).
    run_or_load_method(fold_id, "DT_L2_150_aligned", function() {
      list(s_train = score_pairs_aligned(pair_train, aligned, aligned_pidx),
           s_test  = score_pairs_aligned(pair_test,  aligned, aligned_pidx),
           aligned_key        = aligned$key,
           aligned_config     = aligned$config,
           aligned_cache_file = basename(aligned_cache_file),
           canvas             = aligned$canvas,
           coarse_angles      = aligned$coarse_angles,
           shared_search_minutes = aligned$elapsed_min,
           timing_note = "elapsed = score lookup only; the shared alignment search time is in aligned_runtime.csv")
    }, y_train, y_test, extra_key = aligned$key)

    # ---------------- Baseline: GlobalBetti_L1 ----------------
    run_or_load_method(fold_id, "GlobalBetti_L1", function() {
      list(s_train = score_pairs_l1(pair_train, betti_mat),
           s_test  = score_pairs_l1(pair_test,  betti_mat))
    }, y_train, y_test)

    # ---------------- Baseline: GlobalPI_L2 ----------------
    run_or_load_method(fold_id, "GlobalPI_L2", function() {
      list(s_train = score_pairs_l2(pair_train, gpi_mat),
           s_test  = score_pairs_l2(pair_test,  gpi_mat))
    }, y_train, y_test)

    # ---------------- V23 variants (trained per fold) ----------------
    v23_full_train <- NULL
    v23_full_test  <- NULL

    for (vvar in config$v23_variants) {
      sc <- run_or_load_method(fold_id, vvar, function() {
        feats_v23 <- v23_feat_cache[[vvar]]
        standardizer <- estimate_pair_standardizer(
          pair_mat = pair_train, feats = feats_v23,
          ops = config$v23$symmetric_ops,
          n_sample = config$v23$n_pairs_for_standardize,
          seed = config$v23$seed + fold_id
        )
        model <- fit_logistic_balanced_pairs(
          pair_mat = pair_train, y = y_train, feats = feats_v23,
          ops = config$v23$symmetric_ops, standardizer = standardizer,
          lambda = config$v23$lambda, max_iter = config$v23$max_iter,
          tol = config$v23$tol, batch_size = config$v23$batch_size,
          lr0 = config$v23$lr0, seed = config$v23$seed + fold_id,
          verbose = FALSE
        )
        list(s_train = predict_logistic_pairs(pair_train, feats_v23, model, chunk_size = 8000),
             s_test  = predict_logistic_pairs(pair_test,  feats_v23, model, chunk_size = 8000),
             model   = model)
      }, y_train, y_test)

      if (vvar == "V23_full") {
        v23_full_train <- sc$s_train
        v23_full_test  <- sc$s_test
      }
    }

    # ---------------- OT variants (no training) ----------------
    ot_main_train <- NULL
    ot_main_test  <- NULL

    for (ii in seq_len(nrow(config$ot_variants))) {
      vrow <- config$ot_variants[ii, ]

      sc <- run_or_load_method(fold_id, vrow$method, function() {
        list(s_train = score_pairs_ot(pair_train, bank, config, vrow),
             s_test  = score_pairs_ot(pair_test,  bank, config, vrow))
      }, y_train, y_test)

      if (vrow$method == "OT") {
        ot_main_train <- sc$s_train
        ot_main_test  <- sc$s_test
      }
    }

    # ---------------- Fusion (V23_full + OT) ----------------
    if (is.null(v23_full_train) || is.null(ot_main_train)) stop("Missing V23_full or OT main scores for fusion in fold ", fold_id)

    run_or_load_method(fold_id, "Fusion", function() {
      tune <- tune_fusion_w_on_train(y_train, ot_main_train, v23_full_train, config$fusion_w_grid)
      cat(sprintf("    Fusion weight (fold %d): w = %.2f\n", fold_id, tune$best_w))
      list(s_train = apply_fusion(ot_main_train, v23_full_train, tune$best_w, tune$zfit_ot, tune$zfit_v),
           s_test  = apply_fusion(ot_main_test,  v23_full_test,  tune$best_w, tune$zfit_ot, tune$zfit_v),
           best_w  = tune$best_w,
           zfit_ot = tune$zfit_ot,
           zfit_v  = tune$zfit_v)
    }, y_train, y_test)

    # ---------------- Fusion_eq: equal weighting, w = 0.5 (revision step 2.2b) ----------------
    # Same robust-z normalisation (training-impostor median/MAD) as the tuned Fusion;
    # the only difference is the fixed weight. No tuning, hence no selection on training TAR.
    run_or_load_method(fold_id, "Fusion_eq", function() {
      zfit_ot <- robust_z_fit_train_impostors(ot_main_train, y_train)
      zfit_v  <- robust_z_fit_train_impostors(v23_full_train, y_train)
      list(s_train = apply_fusion(ot_main_train, v23_full_train, 0.5, zfit_ot, zfit_v),
           s_test  = apply_fusion(ot_main_test,  v23_full_test,  0.5, zfit_ot, zfit_v),
           best_w  = 0.5,
           zfit_ot = zfit_ot,
           zfit_v  = zfit_v)
    }, y_train, y_test)

    cat(sprintf("  Fold %d done.\n\n", fold_id))
  }

  # Which cache file each (fold, method) was read from or written to in THIS run.
  # Downstream scripts (bootstrap_ci.R, verify_integration.R) resolve caches through
  # this index rather than rebuilding file names (method-specific suffixes).
  fold_cache_index_df <- bind_rows(fold_cache_index)
  write_csv(fold_cache_index_df, file.path(config$output_dir, "fold_cache_index.csv"))

  # ============================================================
  # AGGREGATE RESULTS
  # ============================================================

  perf_rows <- list()
  tar_rows <- list()
  
  for (x in fold_records) {
    r <- x$res
    perf_rows[[length(perf_rows) + 1]] <- tibble(
      fold = x$fold,
      method = r$method,
      AUC = r$auc,
      EER = r$eer
    )
    tar_rows[[length(tar_rows) + 1]] <- r$tar_tbl %>%
      mutate(fold = x$fold, method = r$method) %>%
      select(fold, method, FAR, TAR)
  }
  
  perf_df <- bind_rows(perf_rows)
  tar_df  <- bind_rows(tar_rows)
  runtime_df <- bind_rows(runtime_rows)
  scores_df <- bind_rows(plot_scores)
  write_csv(scores_df, file.path(config$output_dir, "scores_raw.csv"))

  # Extract and report fusion weights per fold
  fusion_weights <- sapply(seq_len(config$n_folds), function(f) {
    cached <- load_fold_cache(f, "Fusion")
    if (!is.null(cached) && !is.null(cached$best_w)) cached$best_w else NA_real_
  })
  if (any(!is.na(fusion_weights))) {
    cat("\nFusion weights per fold:", paste(sprintf("%.2f", fusion_weights), collapse = ", "), "\n")
    cat(sprintf("  Mean = %.2f, SD = %.3f\n", mean(fusion_weights, na.rm = TRUE), sd(fusion_weights, na.rm = TRUE)))
    write_csv(tibble(fold = seq_along(fusion_weights), best_w = fusion_weights),
              file.path(config$output_dir, "fusion_weights_per_fold.csv"))
  }

  perf_summary <- perf_df %>%
    group_by(method) %>%
    summarize(
      AUC_mean = mean(AUC, na.rm = TRUE),
      AUC_sd   = sd(AUC, na.rm = TRUE),
      EER_mean = mean(EER, na.rm = TRUE),
      EER_sd   = sd(EER, na.rm = TRUE),
      .groups = "drop"
    )
  
  tar_summary <- tar_df %>%
    group_by(method, FAR) %>%
    summarize(
      TAR_mean = mean(TAR, na.rm = TRUE),
      TAR_sd   = sd(TAR, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Per-fold runtime rows measure what happened INSIDE run_or_load_method(). For
  # DT_L2_150_aligned that is a score lookup from the shared pair cache, so its rows
  # are labelled as such here and in the summary; the alignment search time itself is
  # in aligned_runtime.csv (shared_search_minutes) and in each aligned fold cache.
  runtime_df <- runtime_df %>%
    mutate(measures = ifelse(method == "DT_L2_150_aligned",
                             "score lookup only (shared alignment search: see aligned_runtime.csv)",
                             "per-fold scoring (incl. training where applicable)"))
  runtime_summary <- runtime_df %>%
    group_by(method) %>%
    summarize(
      seconds_mean = mean(seconds, na.rm = TRUE),
      seconds_sd   = sd(seconds, na.rm = TRUE),
      measures     = first(measures),
      .groups = "drop"
    )
  if ("DT_L2_150_aligned" %in% runtime_summary$method) {
    runtime_summary$shared_precompute_minutes <- ifelse(runtime_summary$method == "DT_L2_150_aligned",
                                                        aligned$elapsed_min, NA_real_)
  }

  # Write outputs
  write_csv(perf_df, file.path(config$output_dir, "fold_metrics_auc_eer.csv"))
  write_csv(tar_df,  file.path(config$output_dir, "fold_metrics_tar.csv"))
  write_csv(runtime_df, file.path(config$output_dir, "fold_runtimes_seconds.csv"))

  write_csv(perf_summary, file.path(config$output_dir, "summary_auc_eer_mean_sd.csv"))
  write_csv(tar_summary,  file.path(config$output_dir, "summary_tar_mean_sd.csv"))
  write_csv(runtime_summary, file.path(config$output_dir, "summary_runtime_mean_sd.csv"))

  # Publication-friendly display table (order; all other methods follow alphabetically)
  main_methods <- c("DT_L2", "DT_L2_150", "DT_L2_150_clamped", "DT_L2_150_aligned", "GlobalBetti_L1", "GlobalPI_L2",
                    "V23_full", "OT", "Fusion", "Fusion_eq")
  
  tar_wide <- tar_summary %>%
    mutate(FAR_label = paste0("TAR@", format(FAR, scientific = TRUE))) %>%
    select(method, FAR_label, TAR_mean, TAR_sd) %>%
    pivot_wider(names_from = FAR_label, values_from = c(TAR_mean, TAR_sd))
  
  display_table <- perf_summary %>%
    left_join(tar_wide, by = "method") %>%
    mutate(
      AUC = sprintf("%.3f ± %.3f", AUC_mean, AUC_sd),
      EER = sprintf("%.1f%% ± %.1f%%", 100 * EER_mean, 100 * EER_sd)
    )
  
  for (far in config$fars_report) {
    lbl <- paste0("TAR@", format(far, scientific = TRUE))
    mcol <- paste0("TAR_mean_", lbl)
    scol <- paste0("TAR_sd_", lbl)
    if (mcol %in% names(display_table) && scol %in% names(display_table)) {
      display_table[[lbl]] <- sprintf("%.1f%% ± %.1f%%", 100 * display_table[[mcol]], 100 * display_table[[scol]])
    }
  }
  
  display_table <- display_table %>%
    select(method, AUC, EER, any_of(paste0("TAR@", format(config$fars_report, scientific = TRUE)))) %>%
    arrange(match(method, c(main_methods, sort(setdiff(unique(perf_summary$method), main_methods)))))
  
  write_csv(display_table, file.path(config$output_dir, "display_table_main.csv"))
  
  # ============================================================
  # FIGURES
  # ============================================================

  # Publication-ready method labels and colorblind-friendly palette
  method_labels <- c(
    "DT_L2"             = "DT L2 (50\u00b2)",
    "DT_L2_150"         = "DT L2 (150\u00b2)",
    "DT_L2_150_clamped" = "DT L2 (150\u00b2, clamped)",
    "DT_L2_150_aligned" = "DT L2 (150\u00b2, aligned)",
    "GlobalBetti_L1"    = "Betti L1",
    "GlobalPI_L2"       = "PI L2",
    "V23_full"          = "TopoLR",
    "OT"                = "LPOT",
    "Fusion"            = "Fusion",
    "Fusion_eq"         = "Fusion (equal)"
  )

  # Compute pooled AUCs for legend. (Figure COMPOSITION for the revised paper is a later,
  # separate decision; the new methods are included here so they can be judged.)
  roc_plot_methods <- c("DT_L2", "DT_L2_150", "DT_L2_150_aligned", "GlobalBetti_L1", "GlobalPI_L2",
                        "V23_full", "OT", "Fusion", "Fusion_eq")
  pooled_aucs <- scores_df %>%
    filter(method %in% roc_plot_methods) %>%
    group_by(method) %>%
    summarize(auc = as.numeric(pROC::auc(pROC::roc(response = y, predictor = as.numeric(score), quiet = TRUE, direction = "<"))), .groups = "drop")

  # Build legend labels with AUC
  auc_legend <- setNames(
    paste0(method_labels[pooled_aucs$method], " (", sprintf("%.3f", pooled_aucs$auc), ")"),
    pooled_aucs$method
  )

  # Wong (2011) / Okabe-Ito colorblind-safe palette (the 8 colours, plus black)
  cb_palette <- c(
    "DT_L2"             = "#999999",
    "DT_L2_150"         = "#E69F00",
    "DT_L2_150_aligned" = "#F0E442",
    "GlobalBetti_L1"    = "#56B4E9",
    "GlobalPI_L2"       = "#009E73",
    "V23_full"          = "#CC79A7",
    "OT"                = "#0072B2",
    "Fusion"            = "#D55E00",
    "Fusion_eq"         = "#000000"
  )

  roc_long <- scores_df %>%
    filter(method %in% roc_plot_methods) %>%
    group_by(method) %>%
    group_modify(~{
      yy <- .x$y
      ss <- .x$score
      rr <- pROC::roc(response = yy, predictor = as.numeric(ss), quiet = TRUE, direction = "<")
      tibble(FPR = 1 - rr$specificities, TPR = rr$sensitivities)
    }) %>%
    ungroup()

  chance_df <- data.frame(method = "Chance",
                          FPR = 10^seq(-4, 0, length.out = 200))
  chance_df$TPR <- chance_df$FPR
  roc_long <- dplyr::bind_rows(roc_long, chance_df)

  roc_levels    <- c(roc_plot_methods, "Chance")
  roc_colors    <- c(cb_palette, "Chance" = "grey55")
  roc_labels    <- c(auc_legend, "Chance" = "Chance")
  roc_linetypes <- c(setNames(rep("solid", length(roc_plot_methods)), roc_plot_methods),
                     "Chance" = "22")  # short on/off dash so it reads as dashed in the short legend key
  legend_rows   <- if (length(roc_levels) > 8) 3 else 2

  p_roc <- ggplot(roc_long, aes(x = FPR, y = TPR, color = method, linetype = method)) +
    geom_line(linewidth = 0.8) +
    scale_x_log10(
      limits = c(1e-4, 1),
      breaks = 10^(-4:0),
      labels = c(expression(10^{-4}), expression(10^{-3}), expression(10^{-2}), expression(10^{-1}), expression(10^{0}))
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    annotation_logticks(sides = "b") +
    scale_color_manual(values = roc_colors, labels = roc_labels, breaks = roc_levels) +
    scale_linetype_manual(values = roc_linetypes, labels = roc_labels, breaks = roc_levels) +
    theme_minimal(base_size = 11) +
    labs(x = "False Accept Rate (FAR)", y = "True Accept Rate (TAR)", color = NULL, linetype = NULL) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8),
          legend.key.width = grid::unit(0.9, "cm"),
          panel.grid.minor.y = element_blank()) +
    guides(color = guide_legend(nrow = legend_rows), linetype = guide_legend(nrow = legend_rows))

  ggsave(file.path(config$output_dir, "roc_pooled.png"), p_roc, width = 6.151, height = 4.5, dpi = 300)
  ggsave(file.path(config$output_dir, "roc_pooled.pdf"), p_roc, width = 6.151, height = 4.5)
  if (isTRUE(config$export_to_manuscript %||% TRUE)) {
    file.copy(file.path(config$output_dir, "roc_pooled.pdf"), "Manuscript/roc_pooled.pdf", overwrite = TRUE)
    cat("  ROC figure (PDF) copied to Manuscript/\n")
  } else {
    cat("  (export_to_manuscript = FALSE: Manuscript/roc_pooled.pdf left untouched)\n")
  }

  tar_plot_methods <- c("GlobalBetti_L1", "V23_full", "OT", "Fusion", "Fusion_eq")
  tar_plot_df <- tar_summary %>% filter(method %in% tar_plot_methods)

  p_tar <- ggplot(tar_plot_df, aes(x = FAR, y = TAR_mean, color = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = pmax(0, TAR_mean - TAR_sd),
                      ymax = pmin(1, TAR_mean + TAR_sd)),
                  width = 0.0) +
    scale_x_log10() +
    scale_color_manual(values = cb_palette[tar_plot_methods], labels = method_labels[tar_plot_methods]) +
    theme_minimal(base_size = 12) +
    labs(x = "FAR (log scale)", y = "TAR", color = NULL) +
    theme(legend.position = "bottom")

  ggsave(file.path(config$output_dir, "tar_vs_far.png"), p_tar, width = 8, height = 6, dpi = 300)
  ggsave(file.path(config$output_dir, "tar_vs_far.pdf"), p_tar, width = 8, height = 6)
  
  cat("============================================\n")
  cat("DONE.\n")
  cat("Outputs written to: ", config$output_dir, "\n", sep = "")
  cat("Main display table saved as display_table_main.csv\n\n")
  
  print(display_table, n = 50)
  
  invisible(list(
    config = config,
    paths = paths,
    ids = ids,
    folds = folds,
    bank = bank,
    perf_df = perf_df,
    tar_df = tar_df,
    runtime_df = runtime_df,
    perf_summary = perf_summary,
    tar_summary = tar_summary,
    runtime_summary = runtime_summary,
    display_table = display_table
  ))
}

# ============================================================
# PASS 1: 50-ID SANITY / REGRESSION CHECK
# ============================================================

# CONFIG$debug_fast <- FALSE
# CONFIG$n_id <- 50
# 
# # Only keep the main OT variant (saves time, no ablation noise)
# CONFIG$ot_variants <- CONFIG$ot_variants %>%
#   filter(method == "OT")
# 
# # Smaller batch for faster V23 training in sanity run
# CONFIG$v23$batch_size <- 2000
# 
# # Keep Sinkhorn early stopping OFF for reproducibility
# CONFIG$ot$sinkhorn_tol <- 0
# 
# # Run
# results_v2 <- run_prl_suite_v2(CONFIG)

# ============================================================
# FULL PUBLICATION-GRADE RUN (100 IDs)
# - Main text: DT_L2, GlobalBetti_L1, GlobalPI_L2, V23_full, OT, Fusion
# - Appendices: V23 ablations + OT ablations (all variants)
# - Reproducible OT: fixed Sinkhorn iterations (early stopping OFF)
# ============================================================

CONFIG$debug_fast <- FALSE
CONFIG$n_id <- 100

# Keep ALL variants for appendix/ablation tables
# (Leave as defined in CONFIG: 4 V23 variants + 5 OT variants)
# CONFIG$v23_variants unchanged
# CONFIG$ot_variants unchanged

# Reproducibility: fixed-iteration solvers (no adaptive stopping)
CONFIG$ot$sinkhorn_tol <- 0
CONFIG$v23$tol         <- 0

# "Paper-grade" V23 training stability
CONFIG$v23$batch_size <- 10000

# (Optional but recommended) keep seeds exactly as set in CONFIG
# CONFIG$cv_seed and CONFIG$v23$seed already fixed above

# Entry point. Skipped when the file is source()d with the suppress option set
# (e.g., by validate_downsampling.R or pca_match_analysis.R, which only need
# the helper functions and cached feature bank).
if (!isTRUE(getOption("fp.suppress_main", FALSE))) {
  results_pub <- run_prl_suite_v2(CONFIG)
}