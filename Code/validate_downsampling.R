############################################################
# DOWNSAMPLE VALIDATION
#
# The main pipeline computes the Euclidean distance transform on the full
# 300x300 binary image and then downsamples to a 150x150 grid by taking
# every second pixel. This script verifies that the downsampled persistent
# homology is close to the full-resolution PH on a small random sample of
# images, justifying the downsampling step on tractability grounds without
# materially changing the topological summary.
#
# Method: for N random images, compute H1 persistence diagrams of BOTH the
# ridge and valley distance transforms at both resolutions, and report the
# bottleneck distance between full-resolution and downsampled for each channel.
# Filtration values are in the same units (pixel distance) in both cases,
# so the distances are directly comparable. Both channels are validated because
# the pipeline uses persistent homology of both DT_ridge and DT_valley.
#
# Usage (from project root):
#   source("Code/validate_downsampling.R")
#
# Outputs:
#   outputs/downsample_validation.csv  one row per image, with columns
#       image, and for ridge and valley each: n_h1_*_full, n_h1_*_ds, bottleneck_*
#   Console summary: mean and max bottleneck distance.
############################################################

# Suppress fingerprintTDA.R's main entry point only while we source it; restore
# whatever the caller had before so a subsequent source("Code/fingerprintTDA.R")
# in the same R session runs the full pipeline as expected.
.old_fp_suppress <- getOption("fp.suppress_main")
options(fp.suppress_main = TRUE)
source("Code/fingerprintTDA.R")  # for helpers; main run is suppressed
options(fp.suppress_main = .old_fp_suppress)
rm(.old_fp_suppress)

# ----------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------

VALIDATE <- list(
  n_images       = 10,
  seed           = 20260517,
  output_csv     = file.path(CONFIG$output_dir, "downsample_validation.csv"),
  maxscale       = CONFIG$maxscale,
  target_size_ds = CONFIG$target_size   # 150
)

dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------

# Extract H1 birth/death pairs as a 2-column matrix; returns 0-row matrix
# when there are no finite H1 features.
h1_pairs <- function(diag) {
  if (is.null(diag)) return(matrix(numeric(0), ncol = 2,
                                   dimnames = list(NULL, c("Birth", "Death"))))
  m <- as.matrix(diag)
  i <- m[, 1] == 1 & is.finite(m[, 2]) & is.finite(m[, 3])
  if (!any(i)) return(matrix(numeric(0), ncol = 2,
                             dimnames = list(NULL, c("Birth", "Death"))))
  out <- cbind(Birth = m[i, 2], Death = m[i, 3])
  out
}

bottleneck_h1 <- function(diag_a, diag_b) {
  pa <- h1_pairs(diag_a)
  pb <- h1_pairs(diag_b)
  # TDA::bottleneck wraps Diag1, Diag2 with the dimension argument.
  # Construct minimal "diagram" objects with a dimension column expected by TDA.
  to_diag <- function(p) {
    if (nrow(p) == 0) return(matrix(numeric(0), ncol = 3,
                                    dimnames = list(NULL, c("dimension", "Birth", "Death"))))
    cbind(dimension = 1, Birth = p[, "Birth"], Death = p[, "Death"])
  }
  da <- to_diag(pa); db <- to_diag(pb)
  # If either diagram is empty, bottleneck distance equals the max half-persistence
  # of the non-empty side, since every off-diagonal point matches to the diagonal.
  if (nrow(da) == 0 && nrow(db) == 0) return(0)
  if (nrow(da) == 0) return(max((db[, "Death"] - db[, "Birth"]) / 2))
  if (nrow(db) == 0) return(max((da[, "Death"] - da[, "Birth"]) / 2))
  TDA::bottleneck(da, db, dimension = 1)
}

# ----------------------------------------------------------------
# Sample N images
# ----------------------------------------------------------------

paths <- get_paths(CONFIG)
stopifnot(length(paths) >= VALIDATE$n_images)
set.seed(VALIDATE$seed)
sample_paths <- sample(paths, VALIDATE$n_images)

cat(sprintf("Validating downsampling on %d random images (seed = %d)\n\n",
            VALIDATE$n_images, VALIDATE$seed))

# ----------------------------------------------------------------
# Compute and compare
# ----------------------------------------------------------------

# For one channel's full-resolution distance transform, compute the bottleneck
# distance between full-res and downsampled H1 diagrams, plus the H1 counts.
channel_bottleneck <- function(dt_full) {
  dt_ds     <- downsample_matrix(dt_full, VALIDATE$target_size_ds)  # 150 x 150
  diag_full <- compute_persistence(dt_full, maxdim = 1, maxscale = VALIDATE$maxscale)
  diag_ds   <- compute_persistence(dt_ds,   maxdim = 1, maxscale = VALIDATE$maxscale)
  list(n_full = nrow(h1_pairs(diag_full)),
       n_ds   = nrow(h1_pairs(diag_ds)),
       bn     = bottleneck_h1(diag_full, diag_ds))
}

rows <- vector("list", VALIDATE$n_images)
for (i in seq_along(sample_paths)) {
  p <- sample_paths[i]
  bin <- read_binary_fingerprint(p)
  r <- channel_bottleneck(compute_dt_to_ridge(bin))    # ridge channel
  v <- channel_bottleneck(compute_dt_to_valley(bin))   # valley channel

  rows[[i]] <- tibble(
    image             = basename(p),
    n_h1_ridge_full   = r$n_full,
    n_h1_ridge_ds     = r$n_ds,
    bottleneck_ridge  = r$bn,
    n_h1_valley_full  = v$n_full,
    n_h1_valley_ds    = v$n_ds,
    bottleneck_valley = v$bn
  )
  cat(sprintf("  %s  ridge: H1 %4d->%4d bn=%.3f   valley: H1 %4d->%4d bn=%.3f\n",
              basename(p), r$n_full, r$n_ds, r$bn, v$n_full, v$n_ds, v$bn))
}

results <- bind_rows(rows)
write_csv(results, VALIDATE$output_csv)

cat("\n")
cat(sprintf("Summary over %d images:\n", VALIDATE$n_images))
cat(sprintf("  ridge  bottleneck:  mean = %.3f,  max = %.3f,  median = %.3f\n",
            mean(results$bottleneck_ridge),
            max(results$bottleneck_ridge),
            median(results$bottleneck_ridge)))
cat(sprintf("  valley bottleneck:  mean = %.3f,  max = %.3f,  median = %.3f\n",
            mean(results$bottleneck_valley),
            max(results$bottleneck_valley),
            median(results$bottleneck_valley)))
cat(sprintf("\nWrote %s\n", VALIDATE$output_csv))
