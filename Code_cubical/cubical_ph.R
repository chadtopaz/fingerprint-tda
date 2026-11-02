############################################################
# CUBICAL PERSISTENT HOMOLOGY (revision step 1.5, D2/D3)
#
# Sublevel-set persistence on a TRUE cubical complex, the construction the
# submitted manuscript described:
#   "each pixel takes its function value, and each edge and vertex takes the
#    minimum value over the cells it bounds"
# = GUDHI's Bitmap_cubical_complex with top_dimensional_cells = the pixel values
# (T-construction), reached through reticulate.
#
# TDA::gridDiag does NOT do this: its documentation says it "constructs a
# filtration by triangulating the grid"; library = "GUDHI" only selects the
# backend that computes persistence of that simplicial filtration. Every
# reported result uses gridDiag (the text is being corrected); this file exists
# for the supplementary side-by-side.
#
# cubical_diagram() returns the same shape as TDA::gridDiag(...)$diagram:
# a P x 3 matrix (dimension, Birth, Death), infinite deaths replaced.
#
# Parallel use: the pipeline calls this inside future::multisession workers.
# No Python object is kept in a global that could be exported to a worker as a
# dead pointer: the module handle is cached per process and re-imported if the
# cached handle is a null external pointer (what a serialised handle becomes).
# The Python to use is fixed with RETICULATE_PYTHON before reticulate starts
# (configure_cubical_python(); workers inherit the environment variable).
############################################################

# Python with the gudhi module (release 1.1.0): configured, never hard-coded — in order of precedence
# options(fp.python = "<path>"), the environment variable FP_PYTHON, or an already set RETICULATE_PYTHON.
configure_cubical_python <- function() {
  py <- getOption("fp.python", "")
  if (!nzchar(py)) py <- Sys.getenv("FP_PYTHON", "")
  if (!nzchar(py)) py <- Sys.getenv("RETICULATE_PYTHON", "")
  if (!nzchar(py)) stop("no Python configured for the cubical construction: set options(fp.python = <python with gudhi>) or the environment variable FP_PYTHON (see README, 'Cubical construction').", call. = FALSE)
  py <- path.expand(py)
  if (!file.exists(py)) stop("configured Python not found: ", py, call. = FALSE)
  Sys.setenv(RETICULATE_PYTHON = py)          # workers inherit the environment variable
  invisible(py)
}

.cub_env <- new.env(parent = emptyenv())

gudhi_handle <- function() {
  h <- .cub_env$gudhi
  if (!is.null(h) && !reticulate::py_is_null_xptr(h)) return(h)
  if (!reticulate::py_module_available("gudhi"))
    stop("Python module 'gudhi' not found in ", reticulate::py_config()$python,
         "\n  install it with:  pip install gudhi   (in that Python), or set RETICULATE_PYTHON to a Python that has it.",
         call. = FALSE)
  .cub_env$gudhi <- reticulate::import("gudhi", delay_load = FALSE)
  .cub_env$gudhi
}

#' Sublevel-set persistence diagram on the cubical complex of a pixel grid.
#' @param x           numeric matrix of function values on the pixels (top-dimensional cells)
#' @param maxdim      maximum homology dimension returned
#' @param inf_replace value substituted for an infinite death (default: max(x))
#' @return P x 3 numeric matrix with columns dimension, Birth, Death
cubical_diagram <- function(x, maxdim = 1, inf_replace = NULL) {
  stopifnot(is.matrix(x), is.numeric(x))
  if (is.null(inf_replace)) inf_replace <- max(x, na.rm = TRUE)
  gd <- gudhi_handle()
  cc <- gd$CubicalComplex(top_dimensional_cells = x)     # lower cells take the MINIMUM over incident top cells
  cc$compute_persistence()
  out <- vector("list", maxdim + 1L)
  for (d in 0:maxdim) {
    iv <- cc$persistence_intervals_in_dimension(as.integer(d))
    if (length(iv) == 0) next
    iv <- matrix(as.numeric(iv), ncol = 2)
    iv[!is.finite(iv[, 2]), 2] <- inf_replace
    out[[d + 1L]] <- cbind(d, iv[, 1], iv[, 2])
  }
  res <- do.call(rbind, out)
  if (is.null(res)) res <- matrix(numeric(0), ncol = 3)
  colnames(res) <- c("dimension", "Birth", "Death")
  res
}

############################################################
# SELF-TESTS — hand-checkable cases; run once per machine before any result is used.
############################################################
cubical_ph_selftest <- function(verbose = TRUE) {
  ok <- TRUE
  say <- function(...) if (verbose) cat(...)
  say("Cubical PH self-test (", reticulate::py_config()$python, ", gudhi ", as.character(gudhi_handle()$`__version__`), ")\n", sep = "")

  # 1. A ring of eight 0-pixels around a single 9: the sublevel set at t = 0 is an annulus ->
  #    exactly one H1 class born at 0, which dies when the centre enters at t = 9.
  ring <- matrix(0, 3, 3); ring[2, 2] <- 9
  d1 <- cubical_diagram(ring, maxdim = 1); h1 <- d1[d1[, 1] == 1, , drop = FALSE]
  t1 <- nrow(h1) == 1 && isTRUE(all.equal(unname(h1[1, 2]), 0)) && isTRUE(all.equal(unname(h1[1, 3]), 9))
  say("  [", if (t1) "PASS" else "FAIL", "] ring: one H1 class born 0, dies 9\n"); ok <- ok && t1

  # 2. A constant field: one component, no loops.
  d2 <- cubical_diagram(matrix(3, 5, 5), maxdim = 1)
  t2 <- sum(d2[, 1] == 1) == 0 && sum(d2[, 1] == 0) == 1
  say("  [", if (t2) "PASS" else "FAIL", "] constant field: one component, no loops\n"); ok <- ok && t2

  # 3. Transpose invariance (relabelling axes cannot change homology).
  set.seed(1); a <- matrix(runif(36), 6, 6)
  srt <- function(m) m[order(m[, 1], m[, 2], m[, 3]), , drop = FALSE]
  t3 <- isTRUE(all.equal(srt(cubical_diagram(a, 1)), srt(cubical_diagram(t(a), 1)), tolerance = 1e-10))
  say("  [", if (t3) "PASS" else "FAIL", "] transpose invariance\n"); ok <- ok && t3

  # 4. T-construction check. Two DIAGONALLY adjacent 0-pixels in a field of 9s share one
  #    vertex, and in the T-construction that vertex takes the minimum over its four pixels,
  #    i.e. 0, so the two pixels are connected at t = 0: ONE H0 class born at 0 (the
  #    zero-length merge is dropped by GUDHI's default min_persistence = 0). A V-construction
  #    with 4-connectivity would give TWO classes born at 0, one dying at 9.
  diag2 <- matrix(9, 3, 3); diag2[1, 1] <- 0; diag2[2, 2] <- 0
  d4 <- cubical_diagram(diag2, maxdim = 1); h0 <- d4[d4[, 1] == 0, , drop = FALSE]
  n_born0 <- sum(abs(h0[, 2]) < 1e-12)
  t4 <- n_born0 == 1
  say("  [", if (t4) "PASS" else "FAIL", "] T-construction: diagonal 0-pixels connected through their shared vertex (one H0 class born at 0; got ", n_born0, ")\n"); ok <- ok && t4

  # 5. (informational) the two constructions differ on a real-sized random field.
  if (requireNamespace("TDA", quietly = TRUE)) {
    set.seed(2); f <- round(matrix(runif(400, 0, 20), 20, 20), 1)
    tri <- as.matrix(TDA::gridDiag(FUNvalues = f, maxdimension = 1, sublevel = TRUE, library = "GUDHI", printProgress = FALSE)$diagram)
    say("  [INFO] 20 x 20 random field: H1 classes triangulated = ", sum(tri[, 1] == 1), ", cubical = ", sum(cubical_diagram(f, 1)[, 1] == 1), "\n")
  }
  say(if (ok) "ALL TESTS PASSED\n" else "*** SELF-TEST FAILED - do not use results ***\n")
  invisible(ok)
}
