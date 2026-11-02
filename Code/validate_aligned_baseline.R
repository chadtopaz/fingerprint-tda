############################################################
# STEP 1.2 — Validate the aligned-baseline engine BEFORE any full run
#
# Run from the PROJECT ROOT (needs a run directory with a feature bank; options(fp.validate_run), default Outputs_regen_v2):
#   source("Code/validate_aligned_baseline.R")
#
# Checks (spec v3):
#   1  direct masked sums vs FFT at theta != 0, boundary lags, and returned optima
#   2  synthetic pose recovery, correctly posed (inverse transform; fractional
#      translation tolerated; asymmetric patterns; no re-normalisation)
#   3  identity check: search restricted to identity == DT_L2_150 / sqrt(45000)
#   4  monotonicity on the sample: d_sym <= unaligned RMS (necessary, not sufficient)
#   5  swap check: d_sym(a,b) == d_sym(b,a) exactly
#   6  mask count (subsumed by 1; kept as smoke test)
#   7  SELECTION: independent explicit enumeration of the coarse+fine procedure
#      (reduced grid) vs the engine under the same reduced configuration
#   8  benchmark of the actual parallel implementation; go/no-go at 2 h
############################################################

options(fp.suppress_main = TRUE)
source("Code/fingerprintTDA.R")
source("Code/aligned_baseline.R")

REPRO_DIR <- getOption("fp.validate_run", "Outputs_regen_v2")   # any run directory with a feature bank (the engine checks use its images)
TOL_REL   <- 1e-8
RUNTIME_BUDGET_H <- 3      # operational budget (QA: keep the 3-degree grid; ~3 h acceptable if memory is manageable)
MEM_FRACTION     <- 0.70   # projected aggregate worker heap must stay below this fraction of physical RAM
physical_ram_gb <- function() {
  gb <- tryCatch({
    if (Sys.info()[["sysname"]] == "Darwin") as.numeric(system("sysctl -n hw.memsize", intern = TRUE)) / 1e9
    else { l <- readLines("/proc/meminfo", n = 1); as.numeric(gsub("[^0-9]", "", l)) * 1024 / 1e9 }
  }, error = function(e) NA_real_)
  gb
}
set.seed(20261012)

RESULTS <- list()
report <- function(label, status, note = "") {
  RESULTS[[length(RESULTS) + 1]] <<- list(label = label, status = status, note = note)
  cat(sprintf("  [%-9s] %-66s %s\n", status, label, note)); invisible(status)
}
cmp_rel <- function(a, b, label, tol = TOL_REL) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b)) return(report(label, "FAIL", "length"))
  if (any(!is.finite(a)) || any(!is.finite(b))) return(report(label, "FAIL", "non-finite"))
  if (all(a == b)) return(report(label, "EXACT", ""))
  rel <- max(abs(a - b) / pmax(abs(b), 1e-12))
  if (rel <= tol) report(label, "TOLERANCE", sprintf("max rel diff %.2e", rel)) else report(label, "FAIL", sprintf("max rel diff %.2e", rel))
}
section <- function(s) cat("\n", s, "\n", strrep("-", nchar(s)), "\n", sep = "")

# ------------------------------------------------------------------ setup
cfg <- aligned_default_config()
CFG <- CONFIG; CFG$aligned <- cfg
options(future.globals.maxSize = CFG$future_globals_gb * 1024^3)
# NOTE: plan(multisession) is set AFTER the payload preflight (section 0) below.

paths <- get_paths(CFG)
bank_file <- list.files(file.path(REPRO_DIR, "cache"), pattern = "^bank_[0-9a-f]{12}\\.rds$", full.names = TRUE)
stopifnot(length(bank_file) == 1)
bank <- readRDS(bank_file)
stopifnot(length(bank) == length(paths))
n <- length(bank); L <- canvas_size(cfg)
fields <- setNames(lapply(bank, fields_from_bank_entry, frame = cfg$frame), as.character(seq_len(n)))
cat(sprintf("images %d | canvas %d | coarse angles %s | t_max %d | block %d | workers %d\n",
            n, L, paste(range(coarse_angles(cfg$theta_max, cfg$theta_step)), collapse = ".."), cfg$t_max, cfg$block_size, CFG$n_workers))

# ================================================================== 0
section("0. PREFLIGHT — actual serialized export per worker (lazy futures; before any workers start)")
pr <- payload_report(data.frame(pair_table(min(n, 2 * cfg$block_size)), direction = 1L), fields, cfg, L)
cat(sprintf("  scoring worker : export %.2f MB = task payload %.2f MB + overhead %.3f MB   [%s]\n",
            pr$score$total_mb, pr$score$task_mb, pr$score$overhead_mb, paste(head(pr$score$names, 10), collapse = ", ")))
cat(sprintf("  content worker : export %.2f MB = task payload %.2f MB + overhead %.3f MB\n",
            pr$content$total_mb, pr$content$task_mb, pr$content$overhead_mb))
report("scoring export overhead beyond the task < 5 MB (no captured environments)", if (pr$score$overhead_mb < 5) "EXACT" else "FAIL", sprintf("%.3f MB", pr$score$overhead_mb))
report("content export overhead beyond the task < 5 MB", if (pr$content$overhead_mb < 5) "EXACT" else "FAIL", sprintf("%.3f MB", pr$content$overhead_mb))
report("largest scoring task payload < 200 MB", if (pr$score$task_mb < 200) "EXACT" else "FAIL", sprintf("%.1f MB", pr$score$task_mb))
report("largest content task payload < 200 MB", if (pr$content$task_mb < 200) "EXACT" else "FAIL", sprintf("%.1f MB", pr$content$task_mb))
if (any(vapply(RESULTS, `[[`, character(1), "status") == "FAIL"))
  stop("PREFLIGHT FAILED: worker export is not bounded. Not starting ", CFG$n_workers, " workers.", call. = FALSE)
future::plan(future::multisession, workers = CFG$n_workers)

# helper: engine quantities at one (fixed, moving, theta) without the search
engine_at <- function(fi, mj, theta) {
  sc <- rotation_source_coords(L, theta); mrot <- rotated_mask(L, theta, cfg, sc)
  FMa <- fft(frame_mask(cfg$frame, L)); FMb <- fft(mrot); win <- lag_window(L, cfg$t_max)
  list(SSDw = ssd_window(fixed_spectra(fields[[fi]], L), moving_spectra(fields[[mj]], L, mrot, sc), FMa, FMb, win),
       Omw = overlap_window(mrot, FMa, win), win = win)
}
lag_at <- function(win, tx, ty) c(match(ty, win$lags), match(tx, win$lags))   # row = ty, col = tx

# ================================================================== 1
section("1. Direct masked sums vs FFT (theta != 0, boundary lags, optima)")
sample_pairs <- cbind(sample(n, 10), sample(n, 10)); sample_pairs <- sample_pairs[sample_pairs[, 1] != sample_pairs[, 2], , drop = FALSE]
angles <- c(0, 9, -21, cfg$theta_max, -cfg$theta_max)
lags   <- list(c(0, 0), c(cfg$t_max, 0), c(0, -cfg$t_max), c(-cfg$t_max, cfg$t_max), c(17, -23), c(-5, 31))
for (r in seq_len(min(6, nrow(sample_pairs)))) {
  fi <- as.character(sample_pairs[r, 1]); mj <- as.character(sample_pairs[r, 2])
  for (th in angles) {
    E <- engine_at(fi, mj, th)
    for (lg in lags) {
      D <- direct_distance(fields[[fi]], fields[[mj]], th, lg[1], lg[2], cfg, L)
      k <- lag_at(E$win, lg[1], lg[2])
      cmp_rel(E$Omw[k[1], k[2]], D[["omega"]], sprintf("pair %s-%s th %+3d t(%+3d,%+3d) |Omega|", fi, mj, th, lg[1], lg[2]), tol = 0)
      cmp_rel(E$SSDw[k[1], k[2]], D[["ssd"]],   sprintf("pair %s-%s th %+3d t(%+3d,%+3d) SSD",     fi, mj, th, lg[1], lg[2]))
    }
    # at the returned optimum for this theta (reference formulation vs direct sums)
    bt <- best_translation(E$SSDw, E$Omw, E$win, cfg$min_overlap_frac * cfg$frame^2)
    D <- direct_distance(fields[[fi]], fields[[mj]], th, bt[["tx"]], bt[["ty"]], cfg, L)
    cmp_rel(bt[["d"]], D[["d"]], sprintf("pair %s-%s th %+3d optimum t(%+3d,%+3d) d", fi, mj, th, bt[["tx"]], bt[["ty"]]))
    # the OPTIMISED production path must return the same optimum as the reference formulation
    op <- evaluate_requests_at_theta(data.frame(i = as.integer(fi), j = as.integer(mj), direction = 1L), fields[c(fi, mj)], th, cfg, L)
    same_t <- op$tx == bt[["tx"]] && op$ty == bt[["ty"]]
    cmp_rel(op$d, bt[["d"]], sprintf("pair %s-%s th %+3d optimised path d == reference d%s", fi, mj, th, if (same_t) "" else "  (argmin differs!)"))
    if (!same_t) report(sprintf("pair %s-%s th %+3d optimised path argmin == reference argmin", fi, mj, th), "FAIL", sprintf("(%d,%d) vs (%d,%d)", op$tx, op$ty, bt[["tx"]], bt[["ty"]]))
  }
}

# ================================================================== 2
section("2. Synthetic pose recovery (asymmetric pattern; target = inverse transform)")
make_synthetic_fields <- function(seed) {
  set.seed(seed); f <- cfg$frame; bin <- matrix(0L, f, f)
  for (k in 1:6) {                      # asymmetric union of random ellipses
    cx <- runif(1, 25, f - 25); cy <- runif(1, 25, f - 25); a <- runif(1, 6, 18); b <- runif(1, 3, 10); ph <- runif(1, 0, pi)
    xs <- matrix(rep(seq_len(f), each = f), f, f); ys <- matrix(rep(seq_len(f), times = f), f, f)
    u <- (xs - cx) * cos(ph) + (ys - cy) * sin(ph); v <- -(xs - cx) * sin(ph) + (ys - cy) * cos(ph)
    bin[(u / a)^2 + (v / b)^2 <= 1] <- 1L
  }
  dtr <- compute_dt_to_ridge(bin); dtv <- compute_dt_to_valley(bin)
  vec <- c(as.numeric(dtr), as.numeric(dtv)); vec[!is.finite(vec)] <- 0
  vec <- vec / sqrt(mean(vec^2))         # normalised once, as the pipeline does
  list(r = matrix(vec[1:(f * f)], f, f), v = matrix(vec[(f * f + 1):(2 * f * f)], f, f))
}
apply_forward <- function(fld, theta, tx, ty) {
  # Forward transform on the canvas, then crop back to the frame: content leaving the
  # frame is lost, as for a real impression. Amplitudes are preserved and there is NO
  # re-normalisation. Pixels of the frame that received no content (vacated strip) are
  # filled with the channel's mean value rather than 0, because 0 in a DT field means
  # "ridge pixel" and would plant spurious ridges along the frame edge.
  sc <- rotation_source_coords(L, theta); mrot <- rotated_mask(L, theta, cfg, sc)
  o <- (L - cfg$frame) %/% 2; crop <- function(m) m[(o + 1):(o + cfg$frame), (o + 1):(o + cfg$frame)]
  got <- crop(shift_field(mrot, tx, ty))                      # 1 where moved content landed
  fill <- function(m, mu) { m[got == 0] <- mu; m }
  list(r = fill(crop(shift_field(rotate_field(embed_center(fld$r, L), L, mrot, sc), tx, ty)), mean(fld$r)),
       v = fill(crop(shift_field(rotate_field(embed_center(fld$v, L), L, mrot, sc), tx, ty)), mean(fld$v)))
}
cases <- list(c(0, 12, -7), c(9, 0, 0), c(-15, 10, 5), c(7.4, -13, 9), c(21, 18, -20), c(-26, -6, 14),
              c(-30, 5, -3), c(30, -4, 6))        # last two: recovered angle sits ON the angular bound (regression for empty refinement)
for (cs in cases) {
  th <- cs[1]; tx <- cs[2]; ty <- cs[3]
  base <- make_synthetic_fields(seed = round(1000 + th * 10 + tx))
  moved <- apply_forward(base, th, tx, ty)
  tiny_bank <- list(list(dt_vec_150 = c(as.numeric(base$r), as.numeric(base$v))), list(dt_vec_150 = c(as.numeric(moved$r), as.numeric(moved$v))))
  res <- compute_aligned_scores(tiny_bank, paths[1:2], cfg, verbose = FALSE, pairs = data.frame(i = 1L, j = 2L))
  dg <- res$diag
  # direction 1: image 1 (base) fixed, image 2 (moved) moving  -> should recover the INVERSE of (th, t)
  exp_theta <- -th
  exp_t <- as.numeric(-rot_mat(-th) %*% c(tx, ty))
  ok_th <- abs(dg$theta1 - exp_theta) <= 1
  ok_t  <- max(abs(c(dg$tx1, dg$ty1) - round(exp_t))) <= 1          # spec: within 1 px of the NEAREST GRID POINT of the fractional inverse
  note <- sprintf("got th %+.0f t(%+d,%+d); expect th %+.1f t(%+.2f,%+.2f); d/d_unaligned %.3f",
                  dg$theta1, dg$tx1, dg$ty1, exp_theta, exp_t[1], exp_t[2], dg$d_dir1 / dg$d_unaligned)
  report(sprintf("synthetic T*=(%+.1f, %+d, %+d)", th, tx, ty), if (ok_th && ok_t) "EXACT" else "FAIL", note)
  if (abs(th) >= cfg$theta_max)
    report(sprintf("  boundary case flags theta_at_bound"), if (isTRUE(dg$theta_at_bound)) "EXACT" else "FAIL", sprintf("theta = %+.0f", dg$theta))
  rng_ok <- all(c(dg$content_fixed, dg$content_moving) >= 0 & c(dg$content_fixed, dg$content_moving) <= 1)
  report(sprintf("  content fractions in [0,1]"), if (rng_ok) "EXACT" else "FAIL", sprintf("fixed %.3f moving %.3f", dg$content_fixed, dg$content_moving))
}
# real fields: report only (ridge periodicity can make optima ambiguous)
cat("  real-field recovery (reported, not enforced):\n")
for (k in sample(n, 4)) {
  th <- 12; tx <- -9; ty <- 15
  moved <- apply_forward(fields[[as.character(k)]], th, tx, ty)
  tiny_bank <- list(list(dt_vec_150 = c(as.numeric(fields[[as.character(k)]]$r), as.numeric(fields[[as.character(k)]]$v))),
                    list(dt_vec_150 = c(as.numeric(moved$r), as.numeric(moved$v))))
  dg <- compute_aligned_scores(tiny_bank, paths[1:2], cfg, verbose = FALSE, pairs = data.frame(i = 1L, j = 2L))$diag
  exp_t <- as.numeric(-rot_mat(-th) %*% c(tx, ty))
  cat(sprintf("    image %-4d got th %+.0f t(%+d,%+d); expect th %+d t(%+.2f,%+.2f); d/d_unaligned %.3f\n",
              k, dg$theta1, dg$tx1, dg$ty1, -th, exp_t[1], exp_t[2], dg$d_dir1 / dg$d_unaligned))
}

# ================================================================== 2b
section("2b. Content retention: known-count test (source-pixel counting; must be in [0,1])")
{
  f <- cfg$frame; mk <- function(px) { r <- matrix(1, f, f); r[px] <- 0; list(r = r, v = matrix(1, f, f)) }   # ridge pixels where r == 0
  # moving image: 10 isolated ridge pixels on row 75 at columns 10, 20, ..., 100
  px_mov <- cbind(75, seq(10, 100, by = 10)); mov <- mk(px_mov)
  fix <- mk(cbind(c(30, 60, 90, 120), c(20, 70, 100, 130)))      # 4 ridge pixels at frame columns 20, 70, 100, 130
  # theta = 0, shift +60 in x. Omega = frame columns 61..150 (where the shifted frame still overlaps).
  #   moving pixels at x = 10..90 land at 70..150 (inside); x = 100 lands at 160 (outside)  -> 9/10
  #   fixed pixels at x = 70, 100, 130 lie in Omega; x = 20 does not                        -> 3/4
  co <- content_overlap_pair(fix, mov, 0, 60, 0, cfg, L)
  report("known-count: moving retention 9/10 at (0, +60, 0)", if (isTRUE(all.equal(co[["moving"]], 0.9))) "EXACT" else "FAIL", sprintf("got %.3f", co[["moving"]]))
  report("known-count: fixed retention 3/4 at (0, +60, 0)",   if (isTRUE(all.equal(co[["fixed"]], 0.75))) "EXACT" else "FAIL", sprintf("got %.3f", co[["fixed"]]))
  # 30-degree rotation of a single pixel: must be exactly 0 or 1, never 2
  one_px <- mk(cbind(41, 41)); co2 <- content_overlap_pair(fix, one_px, 30, 0, 0, cfg, L)   # off-centre: the defective code reported 2.0 here
  report("single pixel at 30 deg: retention in {0,1} (not duplicated)", if (co2[["moving"]] %in% c(0, 1)) "EXACT" else "FAIL", sprintf("got %.3f", co2[["moving"]]))
}

# ================================================================== 3
section("3. Identity check: search restricted to the identity == DT_L2_150 / sqrt(45000)")
cfg_id <- cfg; cfg_id$theta_max <- 0; cfg_id$refine_half <- 0; cfg_id$t_max <- 0
pp <- pair_table(n)[sample(choose(n, 2), 50), ]
res_id <- compute_aligned_scores(bank, paths, cfg_id, verbose = FALSE, pairs = pp)
dt_mat_150 <- bank_to_matrix(bank, "dt_vec_150")
ref <- -score_pairs_l2(rbind(pp$i, pp$j), dt_mat_150) / sqrt(ncol(dt_mat_150))
cmp_rel(res_id$diag$d_sym, ref, "identity-restricted d_sym vs unaligned RMS (50 pairs)")

# ================================================================== 4, 5
section("4/5. Monotonicity and swap on a sample under the full configuration")
pp2 <- pair_table(n)[sample(choose(n, 2), 40), ]
res_a <- compute_aligned_scores(bank, paths, cfg, verbose = FALSE, pairs = pp2)
res_b <- compute_aligned_scores(bank, paths, cfg, verbose = FALSE, pairs = data.frame(i = pp2$j, j = pp2$i))
mono <- all(res_a$diag$d_sym <= res_a$diag$d_unaligned + 1e-12)
report("monotonicity d_sym <= d_unaligned (40 pairs)", if (mono) "EXACT" else "FAIL",
       sprintf("max(d_sym - d_unaligned) = %.3e", max(res_a$diag$d_sym - res_a$diag$d_unaligned)))
report("swap: d_sym(a,b) == d_sym(b,a)", if (identical(res_a$diag$d_sym, res_b$diag$d_sym)) "EXACT" else "FAIL",
       sprintf("max |diff| = %.3e", max(abs(res_a$diag$d_sym - res_b$diag$d_sym))))
report("swap: directed distances exchange", if (isTRUE(all.equal(res_a$diag$d_dir1, res_b$diag$d_dir2)) && isTRUE(all.equal(res_a$diag$d_dir2, res_b$diag$d_dir1))) "EXACT" else "FAIL", "")
cat(sprintf("  directed disagreement |d_dir1 - d_dir2|: median %.4f, max %.4f; direction-2 wins %d/40\n",
            median(abs(res_a$diag$d_dir1 - res_a$diag$d_dir2)), max(abs(res_a$diag$d_dir1 - res_a$diag$d_dir2)), sum(res_a$diag$direction == 2L)))

# ================================================================== 7
section("7. SELECTION: explicit enumeration of coarse+fine vs engine, reduced grid")
cfg_red <- cfg; cfg_red$theta_max <- 9; cfg_red$t_max <- 8
L_red <- canvas_size(cfg_red)
pp3 <- pair_table(n)[sample(choose(n, 2), 5), ]
res_red <- compute_aligned_scores(bank, paths, cfg_red, verbose = FALSE, pairs = pp3)
for (q in seq_len(nrow(pp3))) {
  fa <- fields[[as.character(pp3$i[q])]]; fb <- fields[[as.character(pp3$j[q])]]
  en <- direct_procedure(fa, fb, cfg_red, L_red)
  dg <- res_red$diag[q, ]
  cmp_rel(dg$d_dir1, en$dir1[["d"]], sprintf("pair %d-%d directed d (dir 1)", pp3$i[q], pp3$j[q]))
  cmp_rel(dg$d_dir2, en$dir2[["d"]], sprintf("pair %d-%d directed d (dir 2)", pp3$i[q], pp3$j[q]))
  cmp_rel(dg$d_sym,  en$d_sym,       sprintf("pair %d-%d symmetric d",        pp3$i[q], pp3$j[q]))
  same_arg <- dg$theta1 == en$dir1[["theta"]] && dg$tx1 == en$dir1[["tx"]] && dg$ty1 == en$dir1[["ty"]]
  cat(sprintf("    pair %d-%d dir1 argmin engine (%+.0f,%+d,%+d) vs enumeration (%+.0f,%+.0f,%+.0f)%s\n",
              pp3$i[q], pp3$j[q], dg$theta1, dg$tx1, dg$ty1, en$dir1[["theta"]], en$dir1[["tx"]], en$dir1[["ty"]],
              if (same_arg) "" else "   <- differs (near-tie?)"))
}

# ================================================================== 8
section("8. Benchmark of the actual parallel implementation: coarse pass, grouped refinement, diagnostics")
P_full   <- choose(n, 2)
coarse   <- coarse_angles(cfg$theta_max, cfg$theta_step); n_coarse <- length(coarse)
fine_set <- setdiff(seq(-cfg$theta_max, cfg$theta_max), coarse); n_fine <- length(fine_set)
# (a) dense coarse pass on 8 blocks (saturates the workers), with measured peak worker memory
nb  <- min(n, 8 * cfg$block_size); ptb <- pair_table(nb)
reqb <- rbind(data.frame(ptb, direction = 1L), data.frame(ptb, direction = 2L))
tb  <- system.time(resb <- run_requests_parallel(reqb, fields[as.character(seq_len(nb))], 9, cfg, L, measure_mem = TRUE))[["elapsed"]]
mem_coarse <- attr(resb, "peak_heap_mb")
T_coarse_full <- tb * (2 * P_full) / nrow(reqb)                      # one coarse angle, all pairs, both directions
# (b) grouped refinement: a fine angle receives ~ (2*refine_half / n_fine) of all (pair,direction) requests
dens <- (2 * cfg$refine_half) / n_fine
reqf <- reqb[sample(nrow(reqb), ceiling(dens * nrow(reqb))), ]
tf  <- system.time(resf <- run_requests_parallel(reqf, fields[as.character(seq_len(nb))], 10, cfg, L, measure_mem = TRUE))[["elapsed"]]
mem_fine <- attr(resf, "peak_heap_mb")
T_fine_full <- tf * (dens * 2 * P_full) / nrow(reqf)                 # one fine angle, all its requests
# (c) content diagnostics on a sample
sel <- sample(nrow(ptb), min(2000, nrow(ptb)))                       # random pairs spanning the block pairs -> many parallel tasks
dgb <- data.frame(i = ptb$i[sel], j = ptb$j[sel], direction = 1L, theta = 9, tx = 5, ty = -7)
n_diag_tasks <- length(build_content_tasks(dgb, fields[as.character(seq_len(nb))], cfg$block_size)$groups)
td  <- system.time(co <- content_overlap_all(dgb, fields[as.character(seq_len(nb))], cfg, L))[["elapsed"]]
T_diag_full <- td * P_full / nrow(dgb)                                 # wall time measured with %d parallel tasks
proj_h <- (n_coarse * T_coarse_full + n_fine * T_fine_full + T_diag_full) / 3600
cat(sprintf("  coarse : %d evals in %.1f s wall  -> %.1f s per angle at full size   | peak worker mem %.0f MB\n", nrow(reqb), tb, T_coarse_full, mem_coarse))
cat(sprintf("  fine   : %d evals in %.1f s wall  -> %.1f s per fine angle at full size (%d fine angles) | peak %.0f MB\n", nrow(reqf), tf, T_fine_full, n_fine, mem_fine))
cat(sprintf("  diag   : %d pairs over %d parallel tasks in %.1f s wall  -> %.1f s at full size\n", nrow(dgb), n_diag_tasks, td, T_diag_full))
cat(sprintf("  PROJECTED full run: %d coarse x %.0f s + %d fine x %.0f s + diag %.0f s = %.2f h  on %d workers\n",
            n_coarse, T_coarse_full, n_fine, T_fine_full, T_diag_full, proj_h, CFG$n_workers))
cat(sprintf("  memory: peak R heap per worker ~%.0f MB x %d workers = %.1f GB (R heap only; excludes the process baseline); scoring task payload %.1f MB (section 0)\n",
            max(mem_coarse, mem_fine), CFG$n_workers, max(mem_coarse, mem_fine) * CFG$n_workers / 1000, pr$score$task_mb))
ram <- physical_ram_gb(); peak_gb <- max(mem_coarse, mem_fine) / 1000; agg_gb <- peak_gb * CFG$n_workers
rec_workers <- if (is.finite(ram)) max(1L, floor(MEM_FRACTION * ram / peak_gb)) else NA_integer_
cat(sprintf("  physical RAM %.0f GB | projected aggregate worker heap %.1f GB (%.0f%% of RAM) | max workers at %.0f%%: %s\n",
            ram, agg_gb, 100 * agg_gb / ram, 100 * MEM_FRACTION, rec_workers))
report(sprintf("runtime go/no-go (<= %.1f h projected, all stages)", RUNTIME_BUDGET_H), if (proj_h <= RUNTIME_BUDGET_H) "EXACT" else "FAIL",
       if (proj_h <= RUNTIME_BUDGET_H) "GO" else "NO-GO on runtime: decide explicitly (reduce block size / workers first; coarsening theta_step is a scientific change)")
report(sprintf("memory go/no-go (aggregate heap <= %.0f%% of RAM)", 100 * MEM_FRACTION),
       if (is.finite(ram) && agg_gb <= MEM_FRACTION * ram) "EXACT" else "FAIL",
       if (is.finite(ram) && agg_gb <= MEM_FRACTION * ram) "GO" else sprintf("NO-GO on memory: set CONFIG$n_workers <- %s (or reduce block_size) and re-benchmark", rec_workers))

# ================================================================== verdict
section("VERDICT")
st <- vapply(RESULTS, `[[`, character(1), "status")
cat(sprintf("  %d checks: %d EXACT, %d TOLERANCE, %d FAIL\n", length(st), sum(st == "EXACT"), sum(st == "TOLERANCE"), sum(st == "FAIL")))
if (any(st == "FAIL")) { cat("  Failures:\n"); for (r in RESULTS) if (r$status == "FAIL") cat("   -", r$label, ":", r$note, "\n") }
cat("  ==>", if (any(st == "FAIL")) "DO NOT RUN THE FULL COMPUTATION" else "ENGINE VALIDATED — proceed to the integration patch and full run", "\n")
VAL_OUT <- getOption("fp.validate_out", "Outputs_validation"); dir.create(VAL_OUT, showWarnings = FALSE, recursive = TRUE)   # the run directory is never written to
saveRDS(list(results = RESULTS, config = cfg, canvas = L, when = Sys.time()), file.path(VAL_OUT, "aligned_validation_report.rds"))
