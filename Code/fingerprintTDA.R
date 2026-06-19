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

############################################################
# TDA FINGERPRINT — MERGED V2 (PRL-READY, NO LEAKAGE)
# OPTIMIZED VERSION — Speedups 1-6 implemented
#
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

CONFIG <- list(
  binary_dir = "Data/FVC2000_binary",
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
  n_workers = max(1, parallel::detectCores() - 2),
  future_globals_gb = 32,
  
  # Output
  output_dir = "Outputs",
  
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

read_binary_fingerprint <- function(path) {
  img <- EBImage::readImage(path)
  if (length(dim(img)) == 3) img <- EBImage::channel(img, "gray")
  mat <- as.matrix(img)
  mat[!is.finite(mat)] <- 0
  bin <- as.integer(mat >= 0.5)
  dim(bin) <- dim(mat)
  if (mean(bin) > 0.5) bin <- 1L - bin
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

sample_anchor_points <- function(bin_ds, n_anchors, seed) {
  set.seed(seed)
  idx <- which(bin_ds == 1L, arr.ind = TRUE)
  if (nrow(idx) == 0) return(matrix(numeric(0), ncol = 2))
  pick <- sample.int(nrow(idx), n_anchors, replace = (nrow(idx) < n_anchors))
  idx[pick, , drop = FALSE]
}

compute_persistence <- function(filt_mat, maxdim = 1, maxscale = 20) {
  x <- pmin(pmax(filt_mat, 0), maxscale)
  x[!is.finite(x)] <- maxscale
  x <- round(x, 1)
  
  tryCatch({
    TDA::gridDiag(
      FUNvalues = x,
      maxdimension = maxdim,
      sublevel = TRUE,
      library = "GUDHI",
      printProgress = FALSE
    )$diagram
  }, error = function(e) NULL)
}

compute_betti_curve <- function(diag, hom_dim, n_points, maxscale) {
  curve <- rep(0, n_points)
  if (is.null(diag)) return(curve)
  
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
  if (is.null(diag)) return(rep(0, n_pix))
  
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
  bin <- read_binary_fingerprint(path)
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

  # Common downsample for PH computations (reuse the 150x150 DTs)
  dt_ridge_dsT  <- ridge_ds150
  dt_valley_dsT <- valley_ds150
  bin_dsT       <- downsample_matrix(bin,       config$target_size)
  
  # Global persistence for Betti & PI
  diag_ridge_global  <- compute_persistence(dt_ridge_dsT,  maxdim = 1, maxscale = config$maxscale)
  diag_valley_global <- compute_persistence(dt_valley_dsT, maxdim = 1, maxscale = config$maxscale)
  
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
      diag_r <- compute_persistence(rp, maxdim = 1, maxscale = config$maxscale)
      ridge_curves[k, ] <- compute_betti_curve(diag_r, 1, config$v23$n_points, config$maxscale)
      X_ridge[k, ] <- diag_to_pi_vec(diag_r, 1, config$maxscale, config$ot$pi_res, config$ot$pi_sigma)
    }
    
    pv <- extract_patch_with_mask(dt_valley_dsT, cx, cy, config$ot$radius, disk_mask)
    vp <- pv$patch
    if (nrow(vp) >= 5 && ncol(vp) >= 5) {
      vp[!pv$mask] <- config$maxscale
      diag_v <- compute_persistence(vp, maxdim = 1, maxscale = config$maxscale)
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
      out <- tryCatch(extract_feature_bank_one(fp, config, disk_mask), error = function(e) NULL)
      p()
      out
    }, future.seed = TRUE)
  })
  cat("Done. Time:", round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n\n")
  
  if (any(map_lgl(bank, is.null))) {
    bad <- which(map_lgl(bank, is.null))
    stop("Some images failed precompute. First bad indices: ", paste(head(bad, 5), collapse = ", "))
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
  digest::digest(key_parts, algo = "sha256")
}

# Fold-score cache key: hashes the feature-bank key AND the training/scoring/
# evaluation parameters that determine per-fold method results. Without this,
# changing a training hyperparameter (e.g., v23$batch_size) while reusing an
# old feature bank would silently load stale fold scores from cache.
fold_cache_key <- function(config) {
  digest::digest(list(
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

  if (file.exists(cache_file)) {
    cat("============================================\n")
    cat("LOADING CACHED FEATURE BANK\n")
    cat("  ", cache_file, "\n")
    cat("============================================\n\n")
    bank <- readRDS(cache_file)
    # Verify cache matches current image set
    if (length(bank) == length(paths)) {
      cat("Cache valid:", length(bank), "images.\n\n")
      return(bank)
    } else {
      cat("Cache size mismatch (", length(bank), "vs", length(paths), "). Recomputing.\n\n")
    }
  }

  bank <- precompute_feature_bank(paths, config)

  cat("Saving feature bank cache to", cache_file, "\n")
  saveRDS(bank, cache_file)
  cat("Cache saved.\n\n")

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
  dt_mat <- bank_to_matrix(bank, "dt_vec")
  dt_mat_150 <- bank_to_matrix(bank, "dt_vec_150")
  betti_mat <- bank_to_matrix(bank, "betti_global_100")
  gpi_mat <- bank_to_matrix(bank, "pi_global")
  
  v23_feat_cache <- setNames(
    lapply(config$v23_variants, function(vv) get_v23_feature_matrix(bank, vv)),
    config$v23_variants
  )
  
  runtime_rows <- list()
  fold_records <- list()
  plot_scores <- list()

  # Fold-level score caching helpers
  fold_cache_dir <- file.path(config$output_dir, "cache", "folds")
  dir.create(fold_cache_dir, showWarnings = FALSE, recursive = TRUE)

  # Fold-score cache key includes BOTH feature-bank parameters AND training/
  # evaluation parameters, so changing v23/ot/fusion settings invalidates
  # stale per-fold score caches.
  fold_key <- substr(fold_cache_key(config), 1, 8)

  fold_cache_path <- function(fold_id, method) {
    file.path(fold_cache_dir, sprintf("fold%d_%s_%s.rds", fold_id, method, fold_key))
  }

  save_fold_cache <- function(fold_id, method, s_train, s_test, y_train, y_test, elapsed, extras = list()) {
    obj <- c(list(s_train = s_train, s_test = s_test,
                  y_train = y_train, y_test = y_test, elapsed = elapsed), extras)
    saveRDS(obj, fold_cache_path(fold_id, method))
  }

  load_fold_cache <- function(fold_id, method) {
    fp <- fold_cache_path(fold_id, method)
    if (file.exists(fp)) readRDS(fp) else NULL
  }

  # Helper: run or load a method for a fold
  run_or_load_method <- function(fold_id, method, compute_fn, y_train, y_test) {
    cached <- load_fold_cache(fold_id, method)
    if (!is.null(cached)) {
      cat(sprintf("  [cached] %s\n", method))
      res <- eval_fold_method(method, cached$y_train, cached$s_train, cached$y_test, cached$s_test, config$fars_report)
      runtime_rows[[length(runtime_rows) + 1]] <<- tibble(fold = fold_id, method = method, seconds = cached$elapsed)
      fold_records[[length(fold_records) + 1]] <<- list(fold = fold_id, res = res)
      plot_scores[[length(plot_scores) + 1]] <<- tibble(fold = fold_id, method = method, y = cached$y_test, score = cached$s_test)
      return(list(s_train = cached$s_train, s_test = cached$s_test))
    }

    t0 <- Sys.time()
    scores <- compute_fn()
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    extras <- scores[setdiff(names(scores), c("s_train", "s_test"))]
    save_fold_cache(fold_id, method, scores$s_train, scores$s_test, y_train, y_test, elapsed, extras)

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
             s_test  = predict_logistic_pairs(pair_test,  feats_v23, model, chunk_size = 8000))
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
           best_w  = tune$best_w)
    }, y_train, y_test)

    cat(sprintf("  Fold %d done.\n\n", fold_id))
  }
  
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
  
  runtime_summary <- runtime_df %>%
    group_by(method) %>%
    summarize(
      seconds_mean = mean(seconds, na.rm = TRUE),
      seconds_sd   = sd(seconds, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Write outputs
  write_csv(perf_df, file.path(config$output_dir, "fold_metrics_auc_eer.csv"))
  write_csv(tar_df,  file.path(config$output_dir, "fold_metrics_tar.csv"))
  write_csv(runtime_df, file.path(config$output_dir, "fold_runtimes_seconds.csv"))
  
  write_csv(perf_summary, file.path(config$output_dir, "summary_auc_eer_mean_sd.csv"))
  write_csv(tar_summary,  file.path(config$output_dir, "summary_tar_mean_sd.csv"))
  write_csv(runtime_summary, file.path(config$output_dir, "summary_runtime_mean_sd.csv"))
  
  # Publication-friendly display table
  main_methods <- c("DT_L2", "DT_L2_150", "GlobalBetti_L1", "GlobalPI_L2", "V23_full", "OT", "Fusion")
  
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
    "DT_L2"          = "DT L2 (50\u00b2)",
    "DT_L2_150"      = "DT L2 (150\u00b2)",
    "GlobalBetti_L1" = "Betti L1",
    "GlobalPI_L2"    = "PI L2",
    "V23_full"       = "TopoLR",
    "OT"             = "LPOT",
    "Fusion"         = "Fusion"
  )

  # Compute pooled AUCs for legend
  roc_plot_methods <- c("DT_L2", "DT_L2_150", "GlobalBetti_L1", "GlobalPI_L2", "V23_full", "OT", "Fusion")
  pooled_aucs <- scores_df %>%
    filter(method %in% roc_plot_methods) %>%
    group_by(method) %>%
    summarize(auc = as.numeric(pROC::auc(pROC::roc(response = y, predictor = as.numeric(score), quiet = TRUE, direction = "<"))), .groups = "drop")

  # Build legend labels with AUC
  auc_legend <- setNames(
    paste0(method_labels[pooled_aucs$method], " (", sprintf("%.3f", pooled_aucs$auc), ")"),
    pooled_aucs$method
  )

  # Wong (2011) colorblind-safe palette (7 colors)
  cb_palette <- c(
    "DT_L2"          = "#999999",
    "DT_L2_150"      = "#E69F00",
    "GlobalBetti_L1" = "#56B4E9",
    "GlobalPI_L2"    = "#009E73",
    "V23_full"       = "#CC79A7",
    "OT"             = "#0072B2",
    "Fusion"         = "#D55E00"
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
    guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))

  ggsave(file.path(config$output_dir, "roc_pooled.png"), p_roc, width = 6.151, height = 4.5, dpi = 300)
  ggsave(file.path(config$output_dir, "roc_pooled.pdf"), p_roc, width = 6.151, height = 4.5)
  file.copy(file.path(config$output_dir, "roc_pooled.pdf"), "Manuscript/roc_pooled.pdf", overwrite = TRUE)
  cat("  ROC figure (PDF) copied to Manuscript/\n")

  tar_plot_methods <- c("GlobalBetti_L1", "V23_full", "OT", "Fusion")
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