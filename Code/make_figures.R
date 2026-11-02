############################################################
# FIGURE GENERATION FOR FINGERPRINT TDA MANUSCRIPT
#
# Generates:
#   fig_pipeline.pdf  Pipeline illustration (binary → DT → PD → Betti)
#                    for both ridge and valley channels
#   fig_patches.pdf  Local patch visualization: 64 anchor locations
#                    on a fingerprint, plus three example patches with
#                    their persistence diagrams
#   fig_score_dists.pdf  Score distributions (TopoLR vs LPOT vs equal-weight Fusion)
#   fig_roc.pdf          Two-panel ROC (FVC2000 DB1 | FVC2002 DB3), pooled test
#                        scores of the principal methods incl. the NBIS reference
#
# Revision (Oct 2026): reads the canonical regenerated inputs and runs
# (Data/FVC2000_binary_regen, Outputs_regen, Outputs_db3, Outputs_reference_*)
# and writes to Outputs_figures/, copying each PDF to "Pattern Recognition/".
# The frozen run directories are never written to.
#
# Usage (from the project root): Rscript Code/make_figures.R
############################################################

# Audit patch (31 Oct 2026): the figure script no longer carries its own copies of the pipeline's helper
# functions; it sources Code/fingerprintTDA.R (main entry suppressed) so that persistence, downsampling,
# patches and anchors are computed by exactly the code that produced the results. The anchors of the
# example image are READ FROM THE CANONICAL FEATURE BANK (anchors_from_bank_entry) and cross-checked
# against sample_anchor_points(), so the patch figure shows the analysis anchors. Images are displayed
# in their natural orientation (EBImage stores (x, y); see show_image()). Patch 1b: show_image() passes
# explicit cell centres to image() so that overlays (pixel_xy) sit on the right pixels, verified by
# check_pixel_coordinates(); the orientation self-test is isolated (options(fp.fig_orientation_test = TRUE)
# runs the test only and regenerates no figure).

options(fp.suppress_main = TRUE)
source("Code/fingerprintTDA.R")        # tidyverse, EBImage, TDA, TDAvec, pROC, digest and the pipeline helpers
library(pROC)

# ============================================================
# CONFIG
# ============================================================

output_dir     <- "Outputs_figures"            # figure files only (never a run directory)
manuscript_dir <- "Pattern Recognition"         # the revised manuscript's directory
run_db1_dir    <- getOption("fp.fig_run_db1", "Outputs_regen_v2")   # canonical DB1 run (scores_raw.csv, feature bank)
run_db3_dir    <- getOption("fp.fig_run_db3", "Outputs_db3_v2")     # DB3 replication (scores_raw.csv)
ref_db1_dir    <- getOption("fp.fig_ref_db1", "Outputs_reference_db1_v2")   # reference matcher fold caches (NBIS_gray), DB1
ref_db3_dir    <- getOption("fp.fig_ref_db3", "Outputs_reference_db3_v2")   # same, DB3
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(manuscript_dir), dir.exists(run_db1_dir), dir.exists(run_db3_dir))
binary_dir <- "Data/FVC2000_binary_regen"       # canonical regenerated binaries (Code/preprocess.py)
example_image <- file.path(binary_dir, "2_1_bin.png")
stopifnot(file.exists(example_image))

export_figure <- function(name) {
  file.copy(file.path(output_dir, name), file.path(manuscript_dir, name), overwrite = TRUE)
  cat(sprintf("  -> %s/%s (copied to %s/)\n", output_dir, name, manuscript_dir))
}

# The canonical run's configuration and feature bank: every parameter below comes from the manifest
# (maxscale, target_size, Betti grid, anchors, radius), not from constants typed here.
run_manifest <- readRDS(file.path(run_db1_dir, "run_manifest.rds"))
CFG_RUN      <- run_manifest$config
stopifnot(identical(normalizePath(CFG_RUN$binary_dir, mustWork = FALSE), normalizePath(binary_dir, mustWork = FALSE)))
bank_file <- list.files(file.path(run_db1_dir, "cache"), pattern = "^bank_[0-9a-f]{12}\\.rds$", full.names = TRUE)
stopifnot(length(bank_file) == 1)
BANK <- readRDS(bank_file)
example_idx <- match(basename(example_image), run_manifest$files); stopifnot(!is.na(example_idx))

maxscale       <- CFG_RUN$maxscale
target_size    <- CFG_RUN$target_size
n_betti_points <- CFG_RUN$v23$n_points

# ------------------------------------------------------------
# Display helpers. EBImage::readImage() returns the image as M[x, y] (x = column of the PNG, left to
# right; y = row of the PNG, top to bottom). base R's image(z) draws z[i, j] at x = i, y = j with y
# increasing UPWARDS, so the natural orientation is image(M[, ncol(M):1]); the previous
# image(t(M[nrow(M):1, ])) reflected the print across the diagonal.
#
# Cell geometry (patch 1b). show_image() passes the cell CENTRES to image() explicitly, together with
# the window [0, 1] x [0, 1]: the n cells along an axis are centred at (i - 0.5)/n and cover
# [(i - 1)/n, i/n]. (image(z) without coordinates would centre them at (i - 1)/(n - 1) and clip the
# half-cells at the margins, up to half a pixel away from (i - 0.5)/n.) pixel_xy() gives the plot
# coordinates of matrix cell (i, j) of M under show_image(M) for overlays; check_pixel_coordinates()
# verifies numerically, against the geometry image() derives from the coordinates it is given, that
# the two agree for every cell.
# ------------------------------------------------------------
image_coords <- function(M) list(x = (seq_len(nrow(M)) - 0.5) / nrow(M), y = (seq_len(ncol(M)) - 0.5) / ncol(M))
show_image <- function(M, ...) {
  cc <- image_coords(M); args <- list(...)
  if (is.null(args$xlim)) args$xlim <- c(0, 1)
  if (is.null(args$ylim)) args$ylim <- c(0, 1)
  do.call(image, c(list(x = cc$x, y = cc$y, z = M[, ncol(M):1, drop = FALSE]), args))
  invisible(cc)
}
pixel_xy <- function(i, j, M) list(x = (i - 0.5) / nrow(M), y = 1 - (j - 0.5) / ncol(M))
# image.default()'s cell boundaries from a vector of centres: midpoints, extended by half a step at both ends
image_cell_bounds <- function(x) { dx <- 0.5 * diff(x); c(x[1L] - dx[1L], x[-length(x)] + dx, x[length(x)] + dx[length(x) - 1L]) }
check_pixel_coordinates <- function(M, tol = 1e-12) {
  cc <- image_coords(M); bx <- image_cell_bounds(cc$x); by <- image_cell_bounds(cc$y)
  n <- nrow(M); m <- ncol(M)
  stopifnot(abs(bx[1] - 0) < tol, abs(bx[n + 1] - 1) < tol, abs(by[1] - 0) < tol, abs(by[m + 1] - 1) < tol)   # the window [0, 1]^2 is exactly the outer cell boundaries
  ii <- rep(seq_len(n), times = m); jj <- rep(seq_len(m), each = n)                  # every cell
  pxy <- pixel_xy(ii, jj, M)
  k <- m + 1L - jj                                                                   # displayed column of original column j (z = M[, m:1])
  ok_x <- abs(pxy$x - cc$x[ii]) < tol & abs(pxy$x - 0.5 * (bx[ii] + bx[ii + 1L])) < tol
  ok_y <- abs(pxy$y - cc$y[k])  < tol & abs(pxy$y - 0.5 * (by[k]  + by[k + 1L]))  < tol
  if (!all(ok_x & ok_y)) stop(sprintf("pixel_xy() disagrees with the image() cell centres for %d of %d cells", sum(!(ok_x & ok_y)), n * m), call. = FALSE)
  invisible(TRUE)
}

# Orientation self-test (options(fp.fig_orientation_test = TRUE); patch 1b: ISOLATED — the test runs, writes
# its two PNGs into Outputs_figures/orientation_test/ and the script ends; no figure is regenerated and nothing
# is copied to the manuscript directory). An asymmetric synthetic pattern (an "F" whose stem is on the LEFT
# and whose bars point RIGHT and sit at the TOP) is written as a PNG, read back through the pipeline's
# reader, displayed with show_image(), and three of its pixels are marked through pixel_xy(). In the saved
# PNG the F must read correctly and every marker must sit on the pattern (the 40 x 60 canvas is drawn
# square, so the F appears stretched; its reading is what matters). A second panel shows the example
# fingerprint with 20 ridge pixels marked. The numerical coordinate check runs first, on both matrices.
ORIENTATION_TEST_ONLY <- isTRUE(getOption("fp.fig_orientation_test", FALSE))
if (ORIENTATION_TEST_ONLY) {
  test_dir <- file.path(output_dir, "orientation_test"); dir.create(test_dir, showWarnings = FALSE, recursive = TRUE)
  tst <- matrix(0L, 60, 40)                       # as a PNG: 60 rows (y), 40 columns (x)
  tst[8:52, 8:11] <- 1L; tst[8:11, 8:30] <- 1L; tst[28:31, 8:24] <- 1L
  png_path <- file.path(test_dir, "orientation_test_input.png")
  png::writePNG(matrix(as.numeric(tst), nrow(tst)), png_path)
  Mt <- read_binary_fingerprint(png_path, polarity = "ridge255")   # M[x, y]
  stopifnot(dim(Mt) == c(40, 60))
  Mb <- read_binary_fingerprint(example_image, polarity = CFG_RUN$polarity %||% "auto")
  check_pixel_coordinates(Mt); check_pixel_coordinates(Mb); check_pixel_coordinates(matrix(0, 7, 3))
  cat("  numerical coordinate check: pixel_xy() == image() cell centres for every cell (40 x 60, ", nrow(Mb), " x ", ncol(Mb), ", 7 x 3)\n", sep = "")
  marks <- list(c(x = 9L, y = 9L), c(x = 29L, y = 9L), c(x = 9L, y = 51L))   # top-left corner, end of the top bar, bottom of the stem (PNG col, row)
  png(file.path(test_dir, "orientation_test.png"), width = 900, height = 450)
  par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
  show_image(Mt, col = c("black", "white"), axes = FALSE, asp = 1, main = "F: stem left, bars right, top; markers on it")
  for (mk in marks) { stopifnot(Mt[mk["x"], mk["y"]] == 1L); pxy <- pixel_xy(mk["x"], mk["y"], Mt); points(pxy$x, pxy$y, pch = 21, bg = "red", col = "white", cex = 2) }
  show_image(Mb, col = c("black", "white"), axes = FALSE, asp = 1, main = "example print; red = 20 ridge pixels")
  set.seed(1); rp <- which(Mb == 1L, arr.ind = TRUE)[sample.int(sum(Mb == 1L), 20), ]
  for (k in seq_len(nrow(rp))) { pxy <- pixel_xy(rp[k, 1], rp[k, 2], Mb); points(pxy$x, pxy$y, pch = 21, bg = "red", col = "white", cex = 1.2) }
  dev.off()
  cat("Orientation test written:", file.path(test_dir, "orientation_test.png"), "\n")
  cat("Orientation test only: no figures regenerated, nothing copied to", manuscript_dir, "\n")
}

if (!ORIENTATION_TEST_ONLY) {   # ---- everything below regenerates the manuscript figures (closed at the end of the file)

# ============================================================
# FIGURE 1: PIPELINE ILLUSTRATION (RIDGE + VALLEY CHANNELS)
# ============================================================

cat("Generating Figure 1: Pipeline illustration...\n")

bin <- read_binary_fingerprint(example_image, polarity = CFG_RUN$polarity %||% "auto")
dt_ridge  <- compute_dt_to_ridge(bin)
dt_valley <- compute_dt_to_valley(bin)

# Downsample for PH (same as main pipeline)
dt_ridge_ds  <- downsample_matrix(dt_ridge,  target_size)
dt_valley_ds <- downsample_matrix(dt_valley, target_size)

# Compute persistence on both channels
diag_ridge  <- compute_persistence(dt_ridge_ds,  maxdim = 1, maxscale = maxscale)
diag_valley <- compute_persistence(dt_valley_ds, maxdim = 1, maxscale = maxscale)

# Extract H1 birth/death pairs from a diagram
extract_h1 <- function(diag) {
  m <- as.matrix(diag)
  i <- m[, 1] == 1 & is.finite(m[, 2]) & is.finite(m[, 3])
  list(birth = m[i, 2], death = m[i, 3])
}
h1_r <- extract_h1(diag_ridge)
h1_v <- extract_h1(diag_valley)

# Compute Betti curves for both channels
betti_r <- compute_betti_curve(diag_ridge,  1, n_betti_points, maxscale)
betti_v <- compute_betti_curve(diag_valley, 1, n_betti_points, maxscale)
t_seq   <- seq(0, maxscale, length.out = n_betti_points)

# Common y-limit so the two Betti curves are visually comparable
betti_ymax <- max(c(betti_r, betti_v))

# Viridis-like palette for heatmaps
heat_cols <- colorRampPalette(c("#440154", "#31688E", "#35B779", "#FDE725"))(256)

# --- Create multi-panel figure: 2 rows x 4 cols, top-left = binary ---
# Row 1: (a) binary    (b) DT_ridge    (c) H1_ridge    (d) Betti_ridge
# Row 2: [blank]       (e) DT_valley   (f) H1_valley   (g) Betti_valley

pdf(file.path(output_dir, "fig_pipeline.pdf"), width = 10, height = 6.0)

layout(matrix(c(1, 2, 3, 4,
                8, 5, 6, 7), nrow = 2, byrow = TRUE),
       widths  = c(1, 1, 1.15, 1.15),
       heights = c(1, 1))

# (a) Binary fingerprint
par(mar = c(0.5, 0.3, 2.2, 0.3), mgp = c(2, 0.7, 0),
    cex.main = 1.4, cex.axis = 1.2, cex.lab = 1.3)
show_image(bin, col = c("black", "white"),
      axes = FALSE, asp = 1, main = expression("(A) Binary image"))

# (b) Ridge distance transform
par(mar = c(0.5, 0.3, 2.2, 0.3))
dt_ridge_clip <- pmin(dt_ridge, maxscale)
show_image(dt_ridge_clip,
      col = heat_cols, axes = FALSE, asp = 1,
      main = expression("(B) " * DT[ridge]))

# (c) Ridge H1 persistence diagram
par(mar = c(3.5, 3.5, 2.2, 0.5))
plot(h1_r$birth, h1_r$death, pch = 16, cex = 1.0,
     col = adjustcolor("#0072B2", alpha.f = 0.6),
     xlim = c(0, maxscale), ylim = c(0, maxscale), asp = 1,
     xlab = "Birth", ylab = "Death",
     main = expression("(C) " * H[1] * " of " * DT[ridge]))
abline(0, 1, lty = 2, col = "grey50")

# (d) Ridge Betti curve
par(mar = c(3.5, 3.5, 2.2, 0.5))
plot(t_seq, betti_r, type = "l", lwd = 2, col = "#D55E00",
     xlim = c(0, maxscale), ylim = c(0, betti_ymax),
     xlab = "Filtration value t", ylab = expression(beta[1](t)),
     main = expression("(D) Betti curve, ridge"))

# (e) Valley distance transform
par(mar = c(0.5, 0.3, 2.2, 0.3))
dt_valley_clip <- pmin(dt_valley, maxscale)
show_image(dt_valley_clip,
      col = heat_cols, axes = FALSE, asp = 1,
      main = expression("(E) " * DT[valley]))

# (f) Valley H1 persistence diagram
par(mar = c(3.5, 3.5, 2.2, 0.5))
plot(h1_v$birth, h1_v$death, pch = 16, cex = 1.0,
     col = adjustcolor("#0072B2", alpha.f = 0.6),
     xlim = c(0, maxscale), ylim = c(0, maxscale), asp = 1,
     xlab = "Birth", ylab = "Death",
     main = expression("(F) " * H[1] * " of " * DT[valley]))
abline(0, 1, lty = 2, col = "grey50")

# (g) Valley Betti curve
par(mar = c(3.5, 3.5, 2.2, 0.5))
plot(t_seq, betti_v, type = "l", lwd = 2, col = "#D55E00",
     xlim = c(0, maxscale), ylim = c(0, betti_ymax),
     xlab = "Filtration value t", ylab = expression(beta[1](t)),
     main = expression("(G) Betti curve, valley"))

# Colour scales for the two distance transform panels (B and E). Each channel
# keeps its own range, since valley distances are much smaller than ridge
# distances; a single shared scale would wash out panel (E). Drawn last so that
# this lands in the empty bottom-left cell of the layout.
rng_ridge  <- range(dt_ridge_clip)
rng_valley <- range(dt_valley_clip)

par(mar = c(0.5, 1.4, 0.5, 1.4))
plot(NA, xlim = c(0, 1), ylim = c(0, 1),
     axes = FALSE, xlab = "", ylab = "", main = "")

draw_dt_scale <- function(y0, y1, rng, lab) {
  rasterImage(as.raster(matrix(heat_cols, nrow = 1)),
              0.04, y0, 0.96, y1, interpolate = TRUE)
  rect(0.04, y0, 0.96, y1, border = "black", lwd = 0.7)
  at <- pretty(rng, n = 4)
  at <- at[at >= rng[1] & at <= rng[2]]
  xs <- if (diff(rng) > 0) 0.04 + (at - rng[1]) / diff(rng) * 0.92 else 0.5
  axis(1, at = xs, labels = at, pos = y0, lwd = 0, lwd.ticks = 0.7,
       tcl = -0.25, cex.axis = 1.05, mgp = c(2, 0.30, 0))
  text(0.5, y1 + 0.055, lab, cex = 1.15, xpd = NA)
}

draw_dt_scale(0.66, 0.76, rng_ridge,  expression("(B) " * DT[ridge] * ", pixels"))
draw_dt_scale(0.28, 0.38, rng_valley, expression("(E) " * DT[valley] * ", pixels"))

dev.off()
export_figure("fig_pipeline.pdf")


# ============================================================
# FIGURE 2: LOCAL PATCH VISUALIZATION
# ============================================================

cat("Generating Figure 2: Local patch visualization...\n")

# Configuration of the main pipeline (from the run manifest)
patch_radius     <- CFG_RUN$ot$radius
n_anchors_patch  <- CFG_RUN$v23$n_anchors

# Downsampled binary (the anchor sampling domain)
bin_ds <- downsample_matrix(bin, target_size)
storage.mode(bin_ds) <- "integer"

# The ANALYSIS anchors of this image: recovered from the canonical feature bank entry (centre attribute
# of the stored OT coordinates) and cross-checked against a fresh draw by the pipeline's helper, which
# now fixes the RNG kind (audit patch, 31 Oct 2026). Before the patch the figure re-drew anchors under
# the session's default RNG and showed a different set from the analysis.
anchors <- anchors_from_bank_entry(BANK[[example_idx]], CFG_RUN)
anchors_check <- sample_anchor_points(bin_ds, n_anchors_patch, seed_from_string(basename(example_image)))
stopifnot(identical(unname(anchors_check), anchors), nrow(anchors) == n_anchors_patch, all(bin_ds[anchors] == 1L))
cat("  anchors: recovered from the canonical bank and reproduced by sample_anchor_points() (", nrow(anchors), " anchors)\n", sep = "")

disk_mask <- make_disk_mask(patch_radius)

# Extract all 64 ridge patches (with disk mask applied) and compute their PDs
patches_ridge <- vector("list", n_anchors_patch)
patches_diag  <- vector("list", n_anchors_patch)
patch_persistence <- numeric(n_anchors_patch)  # total H1 persistence per patch
for (k in seq_len(n_anchors_patch)) {
  cx <- anchors[k, 1]; cy <- anchors[k, 2]
  pr <- extract_patch_with_mask(dt_ridge_ds, cx, cy,
                                patch_radius, disk_mask)
  pm <- pr$patch
  pm[!pr$mask] <- maxscale  # disk-mask: outside-disk pixels → clamp value
  patches_ridge[[k]] <- pm
  d <- compute_persistence(pm, maxdim = 1, maxscale = maxscale)
  patches_diag[[k]] <- d
  if (!is.null(d)) {
    m <- as.matrix(d)
    i <- m[, 1] == 1 & is.finite(m[, 2]) & is.finite(m[, 3])
    if (any(i)) patch_persistence[k] <- sum(pmax(0, m[i, 3] - m[i, 2]))
  }
}

# Select three representative patches: rich (max persistence), medium (median),
# sparse (min, typically near-empty diagram)
ord <- order(patch_persistence, decreasing = TRUE)
sel_idx <- as.integer(c(
  rich   = ord[1],
  medium = ord[ceiling(length(ord) / 2)],
  sparse = ord[length(ord)]
))
sel_letters <- c("B", "C", "D")
sel_descs   <- c("rich", "medium", "sparse")

# --- Create figure ---

pdf(file.path(output_dir, "fig_patches.pdf"), width = 9.5, height = 5.2)

layout(matrix(c(1, 2, 3, 4,
                1, 5, 6, 7), nrow = 2, byrow = TRUE),
       widths = c(2, 1, 1, 1))

# (a) Fingerprint with 64 anchors overlaid; three selected anchors highlighted
par(mar = c(0.5, 0.3, 2.2, 0.3), mgp = c(2, 0.7, 0),
    cex.main = 1.2, cex.axis = 0.9, cex.lab = 1.0)
show_image(bin_ds,
      col = c("black", "white"), axes = FALSE, asp = 1,
      main = "(A) 64 anchor locations")
# Overlay: anchors are matrix cells (i, j) of bin_ds; pixel_xy() maps them to the plot coordinates of show_image()
check_pixel_coordinates(bin_ds)        # numerical check: overlay coordinates == image() cell centres (patch 1b)
axy <- pixel_xy(anchors[, 1], anchors[, 2], bin_ds)
ax <- axy$x; ay <- axy$y
points(ax, ay, pch = 21, bg = "#D55E00", col = "white",
       cex = 0.9, lwd = 0.6)
# Highlight the three selected anchors with large white-filled circles
# (dark blue border, dark blue bold letter) for high-contrast labels
points(ax[sel_idx], ay[sel_idx],
       pch = 21, bg = "white", col = "#0072B2", cex = 3.5, lwd = 2.2)
text(ax[sel_idx], ay[sel_idx], labels = sel_letters,
     col = "#0072B2", cex = 1.4, font = 2)

# Top row, cols 2-4: the three patches with disk masks applied.
# For display, reconstruct each patch with the out-of-disk pixels set to NA so
# they render as background (white) rather than the bright clamp color, and
# scale the colormap to each patch's actual interior data range so the full
# viridis gradient is used.
for (j in seq_along(sel_idx)) {
  k <- sel_idx[j]
  cx <- anchors[k, 1]; cy <- anchors[k, 2]
  pr_disp <- extract_patch_with_mask(dt_ridge_ds, cx, cy,
                                     patch_radius, disk_mask)
  patch_display <- pr_disp$patch
  patch_display[!pr_disp$mask] <- NA
  zmax <- suppressWarnings(max(patch_display, na.rm = TRUE))
  if (!is.finite(zmax) || zmax <= 0) zmax <- 1
  par(mar = c(0.5, 0.3, 2.2, 0.3))
  show_image(patch_display,
        col = heat_cols, axes = FALSE, asp = 1,
        zlim = c(0, zmax),
        main = sprintf("(%s) %s patch", sel_letters[j], sel_descs[j]))
  # Thin grey ring marking the disk boundary
  symbols(0.5, 0.5, circles = 0.5, inches = FALSE, add = TRUE,
          fg = "grey50", lty = 1, lwd = 0.8)
}

# Bottom row, cols 2-4: persistence diagrams of the three patches
pd_letters <- c("E", "F", "G")
for (j in seq_along(sel_idx)) {
  k <- sel_idx[j]
  d <- patches_diag[[k]]
  par(mar = c(3.5, 3.5, 2.2, 0.5))
  plot(NA, NA, xlim = c(0, maxscale), ylim = c(0, maxscale), asp = 1,
       xlab = "Birth", ylab = "Death",
       main = bquote("(" * .(pd_letters[j]) * ") " * H[1] *
                     " of (" * .(sel_letters[j]) * ")"))
  abline(0, 1, lty = 2, col = "grey50")
  if (!is.null(d)) {
    m <- as.matrix(d)
    i <- m[, 1] == 1 & is.finite(m[, 2]) & is.finite(m[, 3])
    if (any(i)) {
      points(m[i, 2], m[i, 3], pch = 16, cex = 1.4,
             col = adjustcolor("#0072B2", alpha.f = 0.7))
    }
  }
}

dev.off()
export_figure("fig_patches.pdf")


# ============================================================
# FIGURE 3: SCORE DISTRIBUTIONS (COMPLEMENTARITY)
# ============================================================

cat("Generating Figure 3: Score distributions...\n")

scores_file <- file.path(run_db1_dir, "scores_raw.csv")
if (!file.exists(scores_file)) {
  cat("  WARNING: scores_raw.csv not found. Run the canonical DB1 pipeline first.\n")
  cat("  Skipping Figure 3.\n")
} else {
  scores <- read_csv(scores_file, show_col_types = FALSE)

  # The primary fusion of the revised paper is the equal-weight rule (Fusion_eq).
  fig_methods <- c("V23_full", "OT", "Fusion_eq")
  fig_labels <- c("V23_full" = "TopoLR", "OT" = "LPOT", "Fusion_eq" = "Fusion")

  scores_fig <- scores %>%
    filter(method %in% fig_methods) %>%
    mutate(
      label = fig_labels[method],
      label = factor(label, levels = c("TopoLR", "LPOT", "Fusion")),
      pair_type = ifelse(y == 1, "Genuine", "Impostor")
    )

  # Compute impostor quantile thresholds (pooled across folds) at FAR=1e-3
  thresholds <- scores_fig %>%
    filter(pair_type == "Impostor") %>%
    group_by(label) %>%
    summarize(thr = quantile(score, probs = 0.999), .groups = "drop")

  pdf(file.path(output_dir, "fig_score_dists.pdf"), width = 6.151, height = 2.8)

  par(mfrow = c(1, 3), mar = c(4, 3.8, 2.2, 1.2), mgp = c(2.2, 0.6, 0), cex.main = 1.1, cex.axis = 0.9, cex.lab = 1.0)

  for (meth in c("TopoLR", "LPOT", "Fusion")) {
    sub <- scores_fig %>% filter(label == meth)
    gen <- sub$score[sub$pair_type == "Genuine"]
    imp <- sub$score[sub$pair_type == "Impostor"]
    thr <- thresholds$thr[thresholds$label == meth]

    # Compute density
    d_gen <- density(gen, n = 512, adjust = 1.2)
    d_imp <- density(imp, n = 512, adjust = 1.2)

    # Set x-limits: use data range but trim extreme outliers
    all_scores <- c(gen, imp)
    xlims <- quantile(all_scores, c(0.005, 0.995))
    # Ensure threshold is visible
    xlims[2] <- max(xlims[2], thr * 1.05)
    xlims <- xlims + c(-1, 1) * diff(xlims) * 0.05  # 5% padding
    # Trim left edge where density is negligible (3% of peak)
    ymax <- max(c(d_gen$y, d_imp$y))
    left_active <- min(
      d_gen$x[which(d_gen$y > ymax * 0.03)[1]],
      d_imp$x[which(d_imp$y > ymax * 0.03)[1]]
    )
    xlims[1] <- max(xlims[1], left_active - diff(c(left_active, xlims[2])) * 0.08)
    ylims <- c(0, max(c(d_gen$y, d_imp$y)) * 1.15)

    plot(d_imp$x, d_imp$y, type = "l", lwd = 2, col = "#E69F00",
         xlim = xlims, ylim = ylims,
         xlab = "Similarity score", ylab = "Density", main = meth)
    polygon(d_imp$x, d_imp$y, col = adjustcolor("#E69F00", 0.2), border = NA)
    lines(d_gen$x, d_gen$y, lwd = 2, col = "#0072B2")
    polygon(d_gen$x, d_gen$y, col = adjustcolor("#0072B2", 0.2), border = NA)

    # Threshold line
    abline(v = thr, lty = 2, lwd = 1.5, col = "grey30")
    text(thr, ylims[2] * 0.93, expression("FAR=" * 10^{-3}), cex = 0.80, pos = 2, col = "grey30")

    # Legend (only on first panel, placed low-left where curves are sparse)
    if (meth == "TopoLR") {
      legend(xlims[1] + diff(xlims) * 0.22, ylims[2] * 0.75,
             legend = c("Genuine", "Impostor"),
             col = c("#0072B2", "#E69F00"), lwd = 2, cex = 0.80, bty = "n")
    }
  }

  dev.off()
  export_figure("fig_score_dists.pdf")
}


# ============================================================
# FIGURE 4: TWO-PANEL ROC (FVC2000 DB1 | FVC2002 DB3)
# ============================================================
# Pooled test-set scores across the five folds of each database (the pipeline's
# scores_raw.csv) plus the NBIS reference matcher read from the reference-run fold
# caches (s_test / y_test, identical pairs and folds). The legend AUCs are computed
# on the pooled scores and serve this visualization only; the paper's estimates are
# the fold means with finger-level intervals (Code/bootstrap_ci.R).

cat("Generating Figure 4: two-panel ROC...\n")

roc_methods <- c("DT_L2_150", "DT_L2_150_aligned", "GlobalBetti_L1", "GlobalPI_L2",
                 "V23_full", "OT", "Fusion_eq", "NBIS_gray")
roc_col <- c("DT_L2_150" = "#E69F00", "DT_L2_150_aligned" = "#E69F00",
             "GlobalBetti_L1" = "#56B4E9", "GlobalPI_L2" = "#009E73",
             "V23_full" = "#CC79A7", "OT" = "#0072B2", "Fusion_eq" = "#000000",
             "NBIS_gray" = "#D55E00")
roc_lty <- c("DT_L2_150" = 1, "DT_L2_150_aligned" = 2, "GlobalBetti_L1" = 1, "GlobalPI_L2" = 1,
             "V23_full" = 1, "OT" = 1, "Fusion_eq" = 1, "NBIS_gray" = 1)
roc_lwd <- c("DT_L2_150" = 1.3, "DT_L2_150_aligned" = 1.3, "GlobalBetti_L1" = 1.3, "GlobalPI_L2" = 1.3,
             "V23_full" = 1.3, "OT" = 1.3, "Fusion_eq" = 1.9, "NBIS_gray" = 1.3)
# Legend labels as plotmath expressions; the pooled AUCs of both databases are
# filled in as "DB1 | DB3" so that one legend serves both panels.
roc_label <- function(m, auc1, auc3) {
  a <- sprintf("%.3f | %.3f", auc1, auc3)
  switch(m,
    "DT_L2_150"         = bquote(DT~L^2~(.(a))),
    "DT_L2_150_aligned" = bquote(DT~L^2*", aligned"~(.(a))),
    "GlobalBetti_L1"    = bquote(Betti~L^1~(.(a))),
    "GlobalPI_L2"       = bquote(PI~L^2~(.(a))),
    "V23_full"          = bquote("TopoLR"~(.(a))),
    "OT"                = bquote("LPOT"~(.(a))),
    "Fusion_eq"         = bquote("Fusion (equal weight)"~(.(a))),
    "NBIS_gray"         = bquote("NBIS (grayscale)"~(.(a))))
}

# Pooled scores of one database: pipeline methods from scores_raw.csv, NBIS from the
# reference-run fold caches listed in its fold_cache_index.csv.
pooled_scores <- function(run_dir, ref_dir) {
  sc <- read_csv(file.path(run_dir, "scores_raw.csv"), show_col_types = FALSE) %>%
    filter(method %in% roc_methods) %>%
    transmute(method, y = as.integer(y), score = as.numeric(score))
  idx <- read_csv(file.path(ref_dir, "fold_cache_index.csv"), show_col_types = FALSE) %>%
    filter(method == "NBIS_gray") %>% arrange(fold)
  stopifnot(nrow(idx) == 5L)
  nb <- bind_rows(lapply(seq_len(nrow(idx)), function(i) {
    cf <- readRDS(file.path(ref_dir, "cache", "folds", idx$file[i]))
    stopifnot(length(cf$s_test) == length(cf$y_test))
    tibble(method = "NBIS_gray", y = as.integer(cf$y_test), score = as.numeric(cf$s_test))
  }))
  bind_rows(sc, nb)
}

# ROC curve (FAR, TAR) and pooled AUC via pROC, as in the pipeline's pooled figure.
# Only points with FAR > 0 can be drawn on the log axis; the curve therefore starts at
# the smallest positive FAR (1 / number of pooled impostors, below the displayed range)
# and the plotting window clips it. For the PDF the curve is thinned on the plotted
# scale: a point is kept when it moves at least `tol` in log10(FAR) or in TAR from the
# last kept point (first and last points always kept), so low-FAR detail, where
# consecutive points are far apart in log10(FAR), is preserved and only the dense
# high-FAR part is thinned.
roc_curve <- function(y, score, tol = 0.002) {
  rr <- pROC::roc(response = y, predictor = score, quiet = TRUE, direction = "<")
  far <- 1 - rr$specificities; tar <- rr$sensitivities
  o <- order(far, tar); far <- far[o]; tar <- tar[o]
  pos <- far > 0; far <- far[pos]; tar <- tar[pos]
  n <- length(far); keep <- logical(n); lf <- log10(far)
  last_lf <- -Inf; last_t <- -Inf
  for (i in seq_len(n)) {
    if (lf[i] - last_lf >= tol || tar[i] - last_t >= tol) {
      keep[i] <- TRUE; last_lf <- lf[i]; last_t <- tar[i]
    }
  }
  keep[c(1L, n)] <- TRUE
  list(far = far[keep], tar = tar[keep], auc = as.numeric(pROC::auc(rr)))
}

# Draws one panel and returns the pooled AUC of every method (named vector).
draw_roc_panel <- function(pooled, title) {
  far_min <- 1e-4
  plot(NA, NA, log = "x", xlim = c(far_min, 1), ylim = c(0, 1), xaxt = "n",
       xlab = "False accept rate (FAR)", ylab = "True accept rate (TAR)", main = title)
  axis(1, at = 10^(-4:0),
       labels = c(expression(10^{-4}), expression(10^{-3}), expression(10^{-2}),
                  expression(10^{-1}), expression(10^{0})))
  abline(h = seq(0, 1, 0.2), v = 10^(-4:0), col = "grey90", lwd = 0.5)
  cf <- 10^seq(log10(far_min), 0, length.out = 100)
  lines(cf, cf, lty = 3, col = "grey55", lwd = 1)        # chance
  aucs <- setNames(numeric(length(roc_methods)), roc_methods)
  for (m in roc_methods) {
    d <- pooled %>% filter(method == m)
    stopifnot(nrow(d) > 0)
    rc <- roc_curve(d$y, d$score)
    lines(rc$far, rc$tar, col = roc_col[m], lty = roc_lty[m], lwd = roc_lwd[m])  # clipped at xlim
    aucs[m] <- rc$auc
  }
  aucs
}

pooled_db1 <- pooled_scores(run_db1_dir, ref_db1_dir)
pooled_db3 <- pooled_scores(run_db3_dir, ref_db3_dir)

pdf(file.path(output_dir, "fig_roc.pdf"), width = 7.2, height = 4.2)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(3.1, 1.1))
par(mar = c(3.4, 3.6, 2.0, 0.8), mgp = c(2.1, 0.6, 0),
    cex.main = 1.05, cex.axis = 0.95, cex.lab = 1.0)
auc_db1 <- draw_roc_panel(pooled_db1, "(A) FVC2000 DB1")
auc_db3 <- draw_roc_panel(pooled_db3, "(B) FVC2002 DB3")
# Shared legend below the panels (one entry per method; AUCs as "DB1 | DB3"), in two
# columns so that the longest entry fits with room to spare (three columns clipped it).
par(mar = c(0, 0.5, 0, 0.5))
plot.new()
labs <- lapply(roc_methods, function(m) roc_label(m, auc_db1[m], auc_db3[m]))
lg <- legend("center", legend = c(as.expression(labs), "Chance"),
             col = c(roc_col[roc_methods], "grey55"), lty = c(roc_lty[roc_methods], 3),
             lwd = c(roc_lwd[roc_methods], 1), ncol = 2, cex = 0.88, bty = "n",
             seg.len = 2.4, x.intersp = 0.8)
# Guard against clipping: the legend box must lie inside the panel's [0, 1] window.
stopifnot(lg$rect$left >= 0, lg$rect$left + lg$rect$w <= 1,
          lg$rect$top <= 1, lg$rect$top - lg$rect$h >= 0)
dev.off()
export_figure("fig_roc.pdf")

# Record the pooled AUCs next to the figure (they are for the legend only; the
# paper's estimates are fold means with intervals).
write_csv(tibble(method = roc_methods, auc_pooled_db1 = as.numeric(auc_db1[roc_methods]),
                 auc_pooled_db3 = as.numeric(auc_db3[roc_methods])),
          file.path(output_dir, "fig_roc_pooled_auc.csv"))

cat("Done.\n")
}   # ---- end of figure generation (if (!ORIENTATION_TEST_ONLY))
