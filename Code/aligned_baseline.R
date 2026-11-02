############################################################
# ALIGNED GEOMETRY BASELINE — engine  (step 1.2, spec v3)
#
# Strengthens DT_L2_150 by allowing a rigid transform of one image before the
# distance is taken over the region where both images have valid data.
#
#   inputs  : the pipeline's normalised two-channel 150x150 DT fields
#             (bank field "dt_vec_150", column-major, ridge then valley)
#   search  : rotation theta (coarse grid = multiples of theta_step anchored at
#             0 within the bound, PLUS both endpoints; then 1-degree refinement
#             over +-refine_half around EACH direction's own coarse optimum),
#             integer translations |tx|,|ty| <= t_max
#   distance: d_T = sqrt( SSD / (2 |Omega|) ) over the overlap Omega, 2 channels
#   score   : -d_sym,   d_sym = min{ d_directed(a->b), d_directed(b->a) }
#
# Conventions (stated once, used everywhere):
#   * canvas L x L, pixel centres at integer coordinates 1..L, centre c=(L+1)/2;
#     the 150x150 frame is embedded centred; matrices are [row = y, col = x].
#   * rotation of the MOVING image by theta degrees about the canvas centre,
#     rendered on the full canvas (no intermediate crop); bilinear interpolation
#     of field and frame mask; mask thresholded at mask_thresh; the field is set
#     to exactly 0 wherever the thresholded mask is 0.
#   * then integer translation t = (tx, ty):  moved(x) = b_theta(x - t).
#   * SSD(t) for ALL t by one inverse FFT per (pair, theta, direction):
#       SSD = sum_c [ a_c^2 (x) M_b  +  M_a (x) b_c^2  -  2 a_c (x) b_c ](t)
#     with (f (x) g)(t) = sum_x f(x) g(x - t) = IFFT( F(f) * Conj(F(g)) ),
#     lag t at array index (t mod L) + 1 on each axis.
#   * L is chosen so the circular correlation equals the linear one for EVERY
#     lag in the window; direct-sum checks therefore hold at any grid point.
#   * transform direction: moving image -> fixed image.
############################################################

suppressPackageStartupMessages({ library(future); library(future.apply); library(digest) })

ALIGNED_SCORER_VERSION <- "1.0"

aligned_default_config <- function() list(
  theta_max        = 30,     # degrees; expansion rule may raise to 40 / 50
  theta_step       = 3,
  refine_half      = 2,      # 1-degree refinement over +-refine_half around each direction's coarse optimum
  t_max            = 40,     # px; expansion rule may raise to 60 / 80
  min_overlap_frac = 0.60,
  frame            = 150,
  mask_thresh      = 0.999,
  interp           = "bilinear",
  block_size       = 60,     # images per block; affects speed/memory only, never results
  scorer_version   = ALIGNED_SCORER_VERSION
)

# ------------------------------------------------------------------ grids
coarse_angles <- function(theta_max, step) {
  k <- floor(theta_max / step)
  sort(unique(c(-theta_max, seq(-k * step, k * step, by = step), theta_max)))
}
refinement_angles <- function(theta_star, cfg, coarse) {
  cand <- seq(theta_star - cfg$refine_half, theta_star + cfg$refine_half, by = 1)
  cand <- cand[abs(cand) <= cfg$theta_max]            # clipped to the current bound
  setdiff(cand, coarse)                                # not already evaluated
}
canvas_size <- function(cfg) {
  th <- min(cfg$theta_max, 45) * pi / 180
  need <- max(cfg$frame + 2 * cfg$t_max, ceiling(cfg$frame * (abs(cos(th)) + abs(sin(th)))) + 2)
  cands <- c(256L, 320L, 384L, 512L, 640L, 768L, 1024L)
  L <- cands[which(cands >= need)[1]]
  if (is.na(L)) stop("grid too large for supported canvas sizes")
  L
}
lag_window <- function(L, t_max) {
  t_max <- as.integer(t_max)
  if (t_max == 0L) return(list(idx = 1L, lags = 0L))            # identity-only search
  list(idx = c(seq_len(t_max + 1L), (L - t_max + 1L):L), lags = c(0L:t_max, (-t_max):(-1L)))
}

# ------------------------------------------------------------------ geometry
embed_center <- function(m, L) {
  f <- nrow(m); o <- (L - f) %/% 2
  out <- matrix(0, L, L); out[(o + 1):(o + f), (o + 1):(o + f)] <- m; out
}
frame_mask <- function(frame, L) embed_center(matrix(1, frame, frame), L)

# forward rotation matrix acting on centred (x, y) column vectors (the map source -> destination)
rot_mat <- function(theta_deg) { th <- theta_deg * pi / 180; matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2, 2) }

bilinear_sample <- function(img, sx, sy) {
  L <- nrow(img)
  x0 <- floor(sx); y0 <- floor(sy); fx <- sx - x0; fy <- sy - y0
  get <- function(xx, yy) {
    ok <- xx >= 1 & xx <= L & yy >= 1 & yy <= L
    v <- numeric(length(xx)); v[ok] <- img[cbind(yy[ok], xx[ok])]; v
  }
  (1 - fx) * (1 - fy) * get(x0, y0) + fx * (1 - fy) * get(x0 + 1, y0) +
    (1 - fx) * fy * get(x0, y0 + 1) + fx * fy * get(x0 + 1, y0 + 1)
}
nearest_sample <- function(img, sx, sy) {
  L <- nrow(img); xx <- round(sx); yy <- round(sy)
  ok <- xx >= 1 & xx <= L & yy >= 1 & yy <= L
  v <- numeric(length(xx)); v[ok] <- img[cbind(yy[ok], xx[ok])]; v
}
# destination -> source coordinates for rotation by theta about the canvas centre
rotation_source_coords <- function(L, theta_deg) {
  c0 <- (L + 1) / 2; th <- theta_deg * pi / 180
  xs <- rep(seq_len(L), each = L); ys <- rep(seq_len(L), times = L)   # column-major traversal: col = x, row = y
  u <- xs - c0; v <- ys - c0
  list(sx = cos(th) * u + sin(th) * v + c0, sy = -sin(th) * u + cos(th) * v + c0)
}
rotated_mask <- function(L, theta_deg, cfg, sc = NULL) {
  if (is.null(sc)) sc <- rotation_source_coords(L, theta_deg)
  m <- bilinear_sample(frame_mask(cfg$frame, L), sc$sx, sc$sy)
  matrix(as.numeric(m > cfg$mask_thresh), L, L)
}
rotate_field <- function(canvas_field, L, mask_rot, sc) matrix(bilinear_sample(canvas_field, sc$sx, sc$sy), L, L) * mask_rot
# integer shift with zero fill:  out(x) = f(x - t)
shift_field <- function(f, tx, ty) {
  L <- nrow(f); out <- matrix(0, L, L)
  rs <- max(1, 1 + ty):min(L, L + ty); cs <- max(1, 1 + tx):min(L, L + tx)
  out[rs, cs] <- f[rs - ty, cs - tx]; out
}

# ------------------------------------------------------------------ fields
fields_from_bank_entry <- function(entry, frame) {
  v <- entry$dt_vec_150; n <- frame * frame
  list(r = matrix(v[1:n], frame, frame), v = matrix(v[(n + 1):(2 * n)], frame, frame))
}
ridge_mask_from_fields <- function(fld) (fld$r == 0) * 1   # ridge pixels: ridge-DT exactly 0

# ------------------------------------------------------------------ spectra
fixed_spectra <- function(fld, L) {
  ar <- embed_center(fld$r, L); av <- embed_center(fld$v, L)
  list(Fr = fft(ar), Fv = fft(av), Fr2 = fft(ar^2), Fv2 = fft(av^2))
}
moving_spectra <- function(fld, L, mask_rot, sc) {
  br <- rotate_field(embed_center(fld$r, L), L, mask_rot, sc)
  bv <- rotate_field(embed_center(fld$v, L), L, mask_rot, sc)
  list(Fr = fft(br), Fv = fft(bv), Fr2 = fft(br^2), Fv2 = fft(bv^2))
}
ixcorr <- function(S) Re(fft(S, inverse = TRUE)) / length(S)

overlap_window <- function(mask_rot, FMa, win) {
  Om <- ixcorr(FMa * Conj(fft(mask_rot)))
  round(Om[win$idx, win$idx, drop = FALSE])
}
ssd_window <- function(Fa, Fb, FMa, FMb, win) {
  S <- Fa$Fr2 * Conj(FMb) + FMa * Conj(Fb$Fr2) - 2 * Fa$Fr * Conj(Fb$Fr) +
       Fa$Fv2 * Conj(FMb) + FMa * Conj(Fb$Fv2) - 2 * Fa$Fv * Conj(Fb$Fv)
  ixcorr(S)[win$idx, win$idx, drop = FALSE]
}
best_translation <- function(SSDw, Omw, win, min_overlap) {
  d2 <- SSDw / (2 * Omw)
  d2[Omw < min_overlap] <- Inf
  d2[d2 < 0] <- 0
  k <- which.min(d2)
  r <- (k - 1) %% nrow(d2) + 1; cc <- (k - 1) %/% nrow(d2) + 1     # row <-> ty, col <-> tx
  c(d = sqrt(d2[k]), tx = win$lags[cc], ty = win$lags[r], omega = Omw[k])
}

# ------------------------------------------------------------------ one task: requests at one theta
# req: data.frame(i, j, direction); fields: named list (by image index) of 150x150 field pairs.
# direction 1: i fixed, j moving;  direction 2: j fixed, i moving.
# Memory: fixed-side spectra cached for the task's fixed images; moving-side spectra
# built once per moving image and discarded.
# Optimised inner loop. By linearity of the FFT and of conjugation,
#   S = sum_c [F(a_c^2) conj(F(M_b)) + F(M_a) conj(F(b_c^2)) - 2 F(a_c) conj(F(b_c))]
#     = A_i + B_j - 2 [F(a_r) conj(F(b_r)) + F(a_v) conj(F(b_v))]
# with A_i = (F(a_r^2) + F(a_v^2)) conj(F(M_b))  (per fixed image, per theta)
# and  B_j = F(M_a) conj(F(b_r^2) + F(b_v^2))    (per moving image, per theta).
# Identical numbers to ssd_window() (the reference formulation, kept for validation),
# with ~3x fewer array operations per evaluation. A periodic gc() bounds the heap
# growth from the 1 MB temporaries that each evaluation creates.
evaluate_requests_at_theta <- function(req, fields, theta_deg, cfg, L, gc_every = 100L) {
  sc <- rotation_source_coords(L, theta_deg)
  mask_rot <- rotated_mask(L, theta_deg, cfg, sc)
  FMa <- fft(frame_mask(cfg$frame, L)); cFMb <- Conj(fft(mask_rot))
  win <- lag_window(L, cfg$t_max)
  Omw <- overlap_window(mask_rot, FMa, win)
  min_ov <- cfg$min_overlap_frac * cfg$frame^2
  nL2 <- as.numeric(L)^2
  fixed_id  <- ifelse(req$direction == 1L, req$i, req$j)
  moving_id <- ifelse(req$direction == 1L, req$j, req$i)
  FX <- new.env()
  get_fixed <- function(k) {
    kk <- as.character(k)
    if (is.null(FX[[kk]])) {
      f <- fields[[kk]]; ar <- embed_center(f$r, L); av <- embed_center(f$v, L)
      FX[[kk]] <- list(Fr = fft(ar), Fv = fft(av), A = (fft(ar^2) + fft(av^2)) * cFMb)
    }
    FX[[kk]]
  }
  out <- matrix(NA_real_, nrow(req), 4); n_done <- 0L
  for (mj in unique(moving_id)) {
    f <- fields[[as.character(mj)]]
    br <- rotate_field(embed_center(f$r, L), L, mask_rot, sc)
    bv <- rotate_field(embed_center(f$v, L), L, mask_rot, sc)
    cFbr <- Conj(fft(br)); cFbv <- Conj(fft(bv))
    B <- FMa * Conj(fft(br^2) + fft(bv^2))
    for (r in which(moving_id == mj)) {
      Fa <- get_fixed(fixed_id[r])
      S  <- Fa$A + B - 2 * (Fa$Fr * cFbr + Fa$Fv * cFbv)
      SSDw <- Re(fft(S, inverse = TRUE)[win$idx, win$idx, drop = FALSE]) / nL2
      out[r, ] <- best_translation(SSDw, Omw, win, min_ov)
      n_done <- n_done + 1L
      if (n_done %% gc_every == 0L) invisible(gc(full = FALSE))
    }
  }
  data.frame(i = req$i, j = req$j, direction = req$direction, theta = theta_deg,
             d = out[, 1], tx = out[, 2], ty = out[, 3], omega = out[, 4])
}

# ------------------------------------------------------------------ parallel driver
pair_table <- function(n) { mm <- combn(n, 2); data.frame(i = mm[1, ], j = mm[2, ]) }
block_of <- function(idx, block_size) (idx - 1L) %/% block_size + 1L

# TOP-LEVEL worker: everything it needs arrives as arguments, so `future` exports only
# the small helper functions, never the caller's frame (all_fields, tasks, ...).
worker_score_task <- function(tk, theta_deg, cfg, L, measure_mem = FALSE) {
  if (measure_mem) invisible(gc(reset = TRUE))
  res <- evaluate_requests_at_theta(tk$req, tk$fields, theta_deg, cfg, L)
  if (measure_mem) {
    g <- gc()                                                      # column layout varies by R version (7 columns on 4.5.x)
    peak_mb_col <- match("max used", colnames(g)) + 1L             # the "(Mb)" column immediately after "max used"
    attr(res, "peak_heap_mb") <- sum(g[, peak_mb_col])             # R heap, megabytes
  }
  res
}
build_score_tasks <- function(req, all_fields, block_size) {
  key <- paste(block_of(req$i, block_size), block_of(req$j, block_size))
  groups <- split(seq_len(nrow(req)), key)
  lapply(groups, function(rows) {
    sub <- req[rows, c("i", "j", "direction"), drop = FALSE]
    list(req = sub, fields = all_fields[as.character(unique(c(sub$i, sub$j)))])
  })
}
# Partition requests into block-pair tasks; each task receives ONLY its images' fields.
run_requests_parallel <- function(req, all_fields, theta_deg, cfg, L, measure_mem = FALSE) {
  tasks <- build_score_tasks(req, all_fields, cfg$block_size)
  # future >= 1.58 attaches 'where' metadata to exported globals that retains the caller's
  # environment (and with it the whole task collection). Drop it for the duration of the call.
  old <- options(future.globals.keepWhere = FALSE); on.exit(options(old), add = TRUE)
  res <- future.apply::future_lapply(tasks, worker_score_task, theta_deg = theta_deg, cfg = cfg, L = L,
                                     measure_mem = measure_mem, future.seed = NULL)
  out <- do.call(rbind, res)
  if (measure_mem) attr(out, "peak_heap_mb") <- max(vapply(res, function(r) attr(r, "peak_heap_mb") %||% NA_real_, numeric(1)), na.rm = TRUE)
  out
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# Measure what `future` ACTUALLY exports for each worker: create a lazy future for a
# representative task and serialize its resolved globals (environments and metadata
# included). Reports the export with and without the intended task payload. Run this
# BEFORE plan(multisession); it needs no workers.
measure_export <- function(expr_quote, task_obj, envir = parent.frame()) {
  old <- options(future.globals.keepWhere = FALSE); on.exit(options(old), add = TRUE)
  # evaluate in the CALLER's environment so the task object and settings named in the
  # expression are found and resolved as globals, exactly as they would be in production
  f <- eval(bquote(future::future(.(expr_quote), lazy = TRUE)), envir = envir)
  g <- f[["globals"]]
  if (is.null(g)) stop("could not access future globals; cannot verify worker export")
  total <- as.numeric(object.size(serialize(g, NULL))) / 1e6
  task  <- as.numeric(object.size(serialize(task_obj, NULL))) / 1e6
  list(total_mb = total, task_mb = task, overhead_mb = total - task, names = names(g))
}
payload_report <- function(req, all_fields, cfg, L) {
  tasks <- build_score_tasks(req, all_fields, cfg$block_size)
  tk <- tasks[[which.max(vapply(tasks, function(x) length(x$fields), integer(1)))]]
  sc <- measure_export(quote(worker_score_task(tk, theta_deg = 9, cfg = cfg, L = L)), tk)
  dg <- data.frame(i = req$i[seq_len(min(50, nrow(req)))], j = req$j[seq_len(min(50, nrow(req)))], direction = 1L, theta = 9, tx = 5, ty = -7)
  ct <- build_content_tasks(dg, all_fields, cfg$block_size)$tasks
  tkc <- ct[[which.max(vapply(ct, function(x) length(x$fields), integer(1)))]]
  cc <- measure_export(quote(worker_content_task(tkc, cfg = cfg, L = L)), tkc)
  list(score = sc, content = cc, n_score_tasks = length(tasks))
}

update_best <- function(best, res, pid_of_row, theta) {
  for (dir in 1:2) {
    sel <- res$direction == dir
    rr <- res[sel, ]; pid <- pid_of_row[sel]
    upd <- rr$d < best[[dir]]$d[pid]; k <- pid[upd]
    best[[dir]]$d[k] <- rr$d[upd]; best[[dir]]$theta[k] <- theta
    best[[dir]]$tx[k] <- rr$tx[upd]; best[[dir]]$ty[k] <- rr$ty[upd]; best[[dir]]$omega[k] <- rr$omega[upd]
  }
  best
}

# Full computation for all unordered pairs. Returns the cache object.
compute_aligned_scores <- function(bank, paths, cfg, verbose = TRUE, pairs = NULL) {
  n <- length(bank); L <- canvas_size(cfg)
  all_fields <- setNames(lapply(bank, fields_from_bank_entry, frame = cfg$frame), as.character(seq_len(n)))
  pt <- if (is.null(pairs)) pair_table(n) else pairs
  P <- nrow(pt)
  coarse <- coarse_angles(cfg$theta_max, cfg$theta_step)
  say <- function(...) if (verbose) cat(...)
  say(sprintf("Aligned baseline: %d images, %d pairs, canvas %d, coarse angles %d [%s..%s], t_max %d, block %d\n",
              n, P, L, length(coarse), min(coarse), max(coarse), cfg$t_max, cfg$block_size))

  best <- lapply(1:2, function(d) data.frame(d = rep(Inf, P), theta = NA_real_, tx = NA_integer_, ty = NA_integer_, omega = NA_real_))
  req_all <- rbind(data.frame(pt, direction = 1L), data.frame(pt, direction = 2L))
  pid_all <- rep(seq_len(P), 2)
  row_key <- function(df) paste(df$i, df$j, df$direction)

  t0 <- Sys.time()
  for (th in coarse) {
    res <- run_requests_parallel(req_all, all_fields, th, cfg, L)
    pid <- pid_all[match(row_key(res), row_key(req_all))]
    best <- update_best(best, res, pid, th)
    say(sprintf("  coarse theta %+3d done  (%.1f min)\n", th, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }

  # refinement, each direction's own neighbourhood, grouped by requested angle (vectorised bookkeeping)
  FR <- do.call(rbind, lapply(1:2, function(dir) {
    ts <- best[[dir]]$theta
    do.call(rbind, lapply(setdiff(seq(-cfg$refine_half, cfg$refine_half), 0), function(off) {
      phi <- ts + off
      keep <- abs(phi) <= cfg$theta_max & !(phi %in% coarse)
      data.frame(i = pt$i[keep], j = pt$j[keep], direction = rep.int(dir, sum(keep)), phi = phi[keep], pid = which(keep))
    }))
  }))
  for (phi in if (!is.null(FR) && nrow(FR) > 0L) sort(unique(FR$phi)) else numeric(0)) {
    sub <- FR[FR$phi == phi, ]
    res <- run_requests_parallel(sub[, c("i", "j", "direction")], all_fields, phi, cfg, L)
    pid <- sub$pid[match(row_key(res), row_key(sub))]
    best <- update_best(best, res, pid, phi)
    say(sprintf("  fine theta %+3d (%d requests) done  (%.1f min)\n", phi, nrow(sub), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }

  # symmetric reduction and diagnostics
  win_dir <- ifelse(best[[2]]$d < best[[1]]$d, 2L, 1L)
  pick <- function(col) ifelse(win_dir == 1L, best[[1]][[col]], best[[2]][[col]])
  diag <- data.frame(i = pt$i, j = pt$j,
                     d_sym = pmin(best[[1]]$d, best[[2]]$d), direction = win_dir,
                     theta = pick("theta"), tx = pick("tx"), ty = pick("ty"), omega = pick("omega"),
                     d_dir1 = best[[1]]$d, theta1 = best[[1]]$theta, tx1 = best[[1]]$tx, ty1 = best[[1]]$ty,
                     d_dir2 = best[[2]]$d, theta2 = best[[2]]$theta, tx2 = best[[2]]$tx, ty2 = best[[2]]$ty)
  diag$theta_at_bound <- abs(diag$theta) >= cfg$theta_max
  diag$t_at_bound     <- pmax(abs(diag$tx), abs(diag$ty)) >= cfg$t_max
  diag$d_unaligned <- vapply(seq_len(P), function(p) {
    a <- all_fields[[as.character(pt$i[p])]]; b <- all_fields[[as.character(pt$j[p])]]
    sqrt(mean(c((a$r - b$r)^2, (a$v - b$v)^2)))
  }, numeric(1))
  co <- content_overlap_all(diag, all_fields, cfg, L)
  diag$content_fixed <- co$fixed; diag$content_moving <- co$moving; diag$content_min <- pmin(co$fixed, co$moving)

  files <- basename(paths); md5 <- unname(tools::md5sum(paths))
  key <- aligned_cache_key(cfg, paths)
  list(key = key, config = cfg, canvas = L, coarse_angles = coarse, files = files, md5 = md5,
       pairs = pt, score = -diag$d_sym, diag = diag,
       elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")), created = Sys.time())
}

# content overlap under the winning transform, each image in its own coordinates,
# denominator = that image's full original ridge count.
#   fixed image : its ridge pixels (fixed coordinates) that lie in Omega
#   moving image: its ORIGINAL ridge pixels whose forward-transformed location
#                 (rotation about the canvas centre, then +t, rounded) lies in Omega.
# Counting source pixels, not resampled destination pixels, keeps the fraction in [0, 1]
# (nearest-neighbour resampling can duplicate a source pixel).
content_overlap_pair <- function(fld_fixed, fld_moving, theta, tx, ty, cfg, L) {
  FM <- frame_mask(cfg$frame, L)
  sc <- rotation_source_coords(L, theta); mrot <- rotated_mask(L, theta, cfg, sc)
  Om <- FM * shift_field(mrot, tx, ty)
  rf <- embed_center(ridge_mask_from_fields(fld_fixed), L)
  f_fixed <- sum(rf * Om) / max(1, sum(rf))
  rm <- embed_center(ridge_mask_from_fields(fld_moving), L)
  src <- which(rm == 1, arr.ind = TRUE)                     # [row = y, col = x] in canvas coords
  if (nrow(src) == 0) return(c(fixed = f_fixed, moving = 0))
  c0 <- (L + 1) / 2; R <- rot_mat(theta)
  xy <- R %*% rbind(src[, 2] - c0, src[, 1] - c0)             # forward rotation of centred (x, y)
  xd <- round(xy[1, ] + c0 + tx); yd <- round(xy[2, ] + c0 + ty)
  ok <- xd >= 1 & xd <= L & yd >= 1 & yd <= L
  inside <- logical(nrow(src)); inside[ok] <- Om[cbind(yd[ok], xd[ok])] == 1
  c(fixed = f_fixed, moving = sum(inside) / nrow(src))
}
worker_content_task <- function(tk, cfg, L) {
  out <- matrix(NA_real_, nrow(tk$dg), 2)
  for (q in seq_len(nrow(tk$dg))) {
    dg <- tk$dg[q, ]
    fi <- if (dg$direction == 1L) dg$i else dg$j
    mj <- if (dg$direction == 1L) dg$j else dg$i
    out[q, ] <- content_overlap_pair(tk$fields[[as.character(fi)]], tk$fields[[as.character(mj)]], dg$theta, dg$tx, dg$ty, cfg, L)
  }
  out
}
build_content_tasks <- function(diag, all_fields, block_size) {
  # group by BLOCK PAIR (as for scoring) so a task never carries more than two blocks of images
  key <- paste(block_of(diag$i, block_size), block_of(diag$j, block_size))
  groups <- split(seq_len(nrow(diag)), key)
  list(groups = groups,
       tasks = lapply(groups, function(rows) list(dg = diag[rows, ], fields = all_fields[as.character(unique(c(diag$i[rows], diag$j[rows])))])))
}
content_overlap_all <- function(diag, all_fields, cfg, L) {
  P <- nrow(diag)
  bt <- build_content_tasks(diag, all_fields, cfg$block_size); groups <- bt$groups; tasks <- bt$tasks
  old <- options(future.globals.keepWhere = FALSE); on.exit(options(old), add = TRUE)
  res <- future.apply::future_lapply(tasks, worker_content_task, cfg = cfg, L = L, future.seed = NULL)
  M <- do.call(rbind, res); ord <- unlist(groups)
  fixed <- numeric(P); moving <- numeric(P); fixed[ord] <- M[, 1]; moving[ord] <- M[, 2]
  stopifnot(all(fixed >= 0 & fixed <= 1), all(moving >= 0 & moving <= 1))
  list(fixed = fixed, moving = moving)
}

# ------------------------------------------------------------------ cache identity + lookup
# Full identity: scorer version + configuration (all result-affecting fields) + canvas +
# ordered image filenames + content hashes.
aligned_cache_key <- function(cfg, paths) {
  digest::digest(list(scorer_version = cfg$scorer_version,
                      config = cfg[setdiff(names(cfg), "block_size")],
                      canvas = canvas_size(cfg),
                      files = basename(paths), md5 = unname(tools::md5sum(paths))), algo = "sha256")
}
load_or_compute_aligned <- function(bank, paths, config) {
  cfg <- config$aligned; key <- aligned_cache_key(cfg, paths)
  dir.create(file.path(config$output_dir, "cache"), showWarnings = FALSE, recursive = TRUE)
  f <- file.path(config$output_dir, "cache", paste0("aligned_", substr(key, 1, 12), ".rds"))
  if (file.exists(f)) {
    obj <- readRDS(f)
    if (identical(obj$key, key) && identical(obj$files, basename(paths))) { cat("Loaded aligned pair scores:", f, "\n"); return(obj) }
    cat("Aligned cache identity mismatch; recomputing.\n")
  }
  obj <- compute_aligned_scores(bank, paths, cfg)
  saveRDS(obj, f); cat("Saved aligned pair scores:", f, "\n")
  obj
}
pair_index_matrix <- function(n, pt) {
  M <- matrix(NA_integer_, n, n)
  M[cbind(pt$i, pt$j)] <- seq_len(nrow(pt)); M[cbind(pt$j, pt$i)] <- seq_len(nrow(pt)); M
}
score_pairs_aligned <- function(pair_mat, aligned, pidx = NULL) {
  if (is.null(pidx)) pidx <- pair_index_matrix(length(aligned$files), aligned$pairs)
  idx <- pidx[cbind(pair_mat[1, ], pair_mat[2, ])]
  if (anyNA(idx)) stop("pair not found in aligned cache")
  aligned$score[idx]
}

# ------------------------------------------------------------------ explicit reference implementation (validation only)
# One (theta, tx, ty): explicit transform + direct masked sums. Returns ssd, omega, d.
direct_distance <- function(fld_fixed, fld_moving, theta_deg, tx, ty, cfg, L, pre = NULL) {
  if (is.null(pre)) {
    sc <- rotation_source_coords(L, theta_deg); mrot <- rotated_mask(L, theta_deg, cfg, sc)
    pre <- list(mrot = mrot,
                br = rotate_field(embed_center(fld_moving$r, L), L, mrot, sc),
                bv = rotate_field(embed_center(fld_moving$v, L), L, mrot, sc))
  }
  FM <- frame_mask(cfg$frame, L)
  ar <- embed_center(fld_fixed$r, L); av <- embed_center(fld_fixed$v, L)
  br <- shift_field(pre$br, tx, ty); bv <- shift_field(pre$bv, tx, ty)
  Om <- FM * shift_field(pre$mrot, tx, ty)
  ssd <- sum(Om * ((ar - br)^2 + (av - bv)^2)); om <- sum(Om)
  c(ssd = ssd, omega = om, d = if (om > 0) sqrt(ssd / (2 * om)) else NA_real_)
}
# Explicit enumeration of the prescribed procedure for ONE unordered pair under cfg
# (rotate once per angle; loop lags). Intended for a reduced grid.
direct_procedure <- function(fa, fb, cfg, L) {
  coarse <- coarse_angles(cfg$theta_max, cfg$theta_step); min_ov <- cfg$min_overlap_frac * cfg$frame^2
  lags <- (-cfg$t_max):cfg$t_max
  one_dir <- function(ff, fm) {
    eval_at <- function(th) {
      sc <- rotation_source_coords(L, th); mrot <- rotated_mask(L, th, cfg, sc)
      pre <- list(mrot = mrot, br = rotate_field(embed_center(fm$r, L), L, mrot, sc), bv = rotate_field(embed_center(fm$v, L), L, mrot, sc))
      best <- c(d = Inf, tx = NA_real_, ty = NA_real_, omega = NA_real_)
      for (tx in lags) for (ty in lags) {
        r <- direct_distance(ff, fm, th, tx, ty, cfg, L, pre)
        if (r[["omega"]] >= min_ov && is.finite(r[["d"]]) && r[["d"]] < best[["d"]])
          best <- c(d = r[["d"]], tx = tx, ty = ty, omega = r[["omega"]])
      }
      best
    }
    cb <- lapply(coarse, eval_at); ds <- vapply(cb, `[[`, numeric(1), "d"); k <- which.min(ds)
    best <- c(cb[[k]], theta = coarse[k])
    for (phi in refinement_angles(coarse[k], cfg, coarse)) { r <- eval_at(phi); if (r[["d"]] < best[["d"]]) best <- c(r, theta = phi) }
    best
  }
  d1 <- one_dir(fa, fb); d2 <- one_dir(fb, fa)
  list(dir1 = d1, dir2 = d2, d_sym = min(d1[["d"]], d2[["d"]]))
}
