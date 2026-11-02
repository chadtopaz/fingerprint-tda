############################################################
# CANONICAL RUN DRIVER (audit patch, 31 Oct 2026) — the full pipeline on one database with the
# frozen configuration of the paper, from the regenerated binaries, into a FRESH output directory.
# This is the clean entry point for reproduction: it needs no historical output directory.
#
#   options(fp.db = "db1")                                 # or "db3"
#   options(fp.out_dir = "Outputs_regen_v2")               # default: db1 -> Outputs_regen_v2, db3 -> Outputs_db3_v2
#   options(fp.reference_run = "Outputs_regen")            # optional, default as below; "" = none
#   options(fp.reuse_aligned = TRUE)                       # copy the reference run's aligned pair cache (saves the ~2.4 h search)
#   options(fp.n_workers = 20)                             # optional: worker count (operational only; default detectCores() - 2, NA-safe)
#   source("Code/run_canonical.R")
#
# What it does, in order:
#   1. input integrity: 800 binaries whose md5s equal the preprocessing manifest (0 errors);
#   2. configuration: CONFIG of Code/fingerprintTDA.R with binary_dir / output_dir set, export_to_manuscript
#      FALSE and the FROZEN alignment grid (Code/frozen_grid.R);
#   3. pre-flight against the reference run, when one is given and exists: same configuration apart from the
#      directories, same file names, identities and identity folds — so the new run is like-for-like;
#   4. optional reuse of the reference run's aligned pair cache (identity re-validated by the engine from its
#      key and file list; the search is repeated if it does not match);
#   5. a COLD run of every method (the output directory must be new or empty), then run_record.rds.
# Afterwards:  options(fp.verify_a = <reference run>, fp.verify_b = <this run>); source("Code/verify_run_vs_canonical.R")
#              options(fp.bootstrap_dir = <this run>); source("Code/bootstrap_ci.R")
# Writes only into the output directory. Never touches Data/, the reference run, or the manuscript.
############################################################

options(fp.suppress_main = TRUE)
source("Code/fingerprintTDA.R")
source("Code/frozen_grid.R")

DB <- match.arg(getOption("fp.db", "db1"), c("db1", "db3"))
SET <- list(
  db1 = list(label = "FVC2000 DB1", binary_dir = "Data/FVC2000_binary_regen", default_out = "Outputs_regen_v2", default_ref = "Outputs_regen"),
  db3 = list(label = "FVC2002 DB3 A", binary_dir = "Data/FVC2002_DB3_binary", default_out = "Outputs_db3_v2",   default_ref = "Outputs_db3")
)[[DB]]
OUT_DIR <- getOption("fp.out_dir", SET$default_out)
REF_DIR <- getOption("fp.reference_run", SET$default_ref)
REUSE_ALIGNED <- isTRUE(getOption("fp.reuse_aligned", TRUE))
section <- function(s) cat("\n", s, "\n", strrep("-", nchar(s)), "\n", sep = "")
# worker count: operational, never scientific (future_lapply(future.seed = TRUE) seeds every element, so results do not
# depend on it); validated override, default detectCores() - 2 with NA handled in CONFIG
OPERATIONAL_FIELDS <- c("n_workers", "future_globals_gb")
if (!is.null(getOption("fp.n_workers"))) {
  nw <- suppressWarnings(as.integer(getOption("fp.n_workers")))
  if (length(nw) != 1L || is.na(nw) || nw < 1L) stop("fp.n_workers must be a single integer >= 1", call. = FALSE)
  nc <- suppressWarnings(parallel::detectCores())
  if (!is.na(nc) && nw > nc) stop("fp.n_workers (", nw, ") exceeds the detected core count (", nc, ")", call. = FALSE)
  CONFIG$n_workers <- nw
}
stopifnot(is.numeric(CONFIG$n_workers), length(CONFIG$n_workers) == 1L, !is.na(CONFIG$n_workers), CONFIG$n_workers >= 1)

# ---- cold-run guard
if (dir.exists(OUT_DIR) && length(list.files(OUT_DIR, recursive = TRUE, all.files = TRUE, no.. = TRUE)) > 0)
  stop(OUT_DIR, " already contains files; this driver makes a COLD run into a new or empty directory.", call. = FALSE)
stopifnot(isTRUE(CONFIG$fold_cache_version >= 2L))

# ---- 1. input integrity
section(paste("Inputs:", SET$label))
pngs <- list.files(SET$binary_dir, pattern = "_bin\\.png$")
stopifnot(length(pngs) == 800)
pm <- jsonlite::fromJSON(file.path(SET$binary_dir, "preprocess_manifest.json"), simplifyVector = TRUE)
stopifnot(pm$n_images == 800, pm$n_errors == 0, all(pm$images$status == "ok"), setequal(pm$images$output, pngs))
md5_now <- unname(tools::md5sum(file.path(SET$binary_dir, pm$images$output)))
if (!identical(md5_now, pm$images$output_md5))
  stop("binaries differ from their preprocessing manifest (output_md5): ", sum(md5_now != pm$images$output_md5), " files.", call. = FALSE)
cat("  ", SET$binary_dir, ": 800 binaries, md5s equal to preprocess_manifest.json (0 errors); enhancer ",
    pm$environment$fingerprint_enhancer, ", source sha ", substr(pm$environment$enhancer_source_sha256, 1, 12), "...\n", sep = "")

# ---- 2. configuration
CONFIG$binary_dir           <- SET$binary_dir
CONFIG$output_dir           <- OUT_DIR
CONFIG$export_to_manuscript <- FALSE
CONFIG                      <- apply_frozen_grid(CONFIG)

# ---- 3. pre-flight against the reference run (optional)
paths <- get_paths(CONFIG)
stopifnot(length(paths) == 800)
ids   <- map_chr(paths, get_identity)
folds <- make_identity_folds(observed_ids = ids, n_folds = CONFIG$n_folds, seed = CONFIG$cv_seed)
if (nzchar(REF_DIR) && file.exists(file.path(REF_DIR, "run_manifest.rds"))) {
  section(paste("Pre-flight against the reference run", REF_DIR))
  ref <- readRDS(file.path(REF_DIR, "run_manifest.rds"))
  cfg_a <- ref$config; cfg_b <- CONFIG
  for (fld in c("binary_dir", "output_dir", "export_to_manuscript", OPERATIONAL_FIELDS)) { cfg_a[[fld]] <- NULL; cfg_b[[fld]] <- NULL }
  if (!identical(as.integer(ref$config$n_workers), as.integer(CONFIG$n_workers))) cat("  (worker count differs from the reference run: ", ref$config$n_workers, " vs ", CONFIG$n_workers, "; operational only)\n", sep = "")
  # fields compared as EFFECTIVE values (absent = default), then removed
  for (fld in c("complex_type", "polarity")) {
    dflt <- c(complex_type = "triangulated", polarity = "auto")[[fld]]
    va <- cfg_a[[fld]] %||% dflt; vb <- cfg_b[[fld]] %||% dflt
    if (!identical(va, vb)) stop(fld, " differs from the reference run: ", va, " vs ", vb, call. = FALSE)
    cfg_a[[fld]] <- NULL; cfg_b[[fld]] <- NULL
  }
  if (!identical(cfg_a, cfg_b)) {
    diff_fields <- names(cfg_b)[!vapply(names(cfg_b), function(n) identical(cfg_a[[n]], cfg_b[[n]]), logical(1))]
    stop("CONFIG differs from the reference run's configuration in: ",
         paste(c(diff_fields, setdiff(names(cfg_a), names(cfg_b))), collapse = ", "), call. = FALSE)
  }
  stopifnot(identical(basename(paths), ref$files), identical(ids, ref$ids))
  stopifnot(identical(unname(lapply(folds, as.character)), unname(lapply(ref$folds, function(f) as.character(f$test_ids)))))
  same_imgs <- identical(unname(tools::md5sum(paths)), ref$md5)
  cat("  configuration identical apart from the directories and the operational fields (", paste(OPERATIONAL_FIELDS, collapse = ", "), "); 800 file names, identities and 5 identity folds identical;",
      if (same_imgs) "images identical (same md5s)." else "images DIFFER from the reference run (expected only for a different database).", "\n")
  # ---- 4. aligned pair cache reuse
  if (REUSE_ALIGNED) {
    key_now <- aligned_cache_key(CONFIG$aligned, paths)
    src_f <- file.path(REF_DIR, "cache", paste0("aligned_", substr(key_now, 1, 12), ".rds"))
    if (file.exists(src_f)) {
      obj <- readRDS(src_f)
      if (identical(obj$key, key_now) && identical(obj$files, basename(paths))) {
        dir.create(file.path(OUT_DIR, "cache"), recursive = TRUE, showWarnings = FALSE)
        file.copy(src_f, file.path(OUT_DIR, "cache", basename(src_f)))
        cat("  aligned pair cache reused from the reference run (key ", substr(key_now, 1, 12), "...; identity re-validated by the engine at load).\n", sep = "")
      } else cat("  reference aligned cache does not match this run's key/files; the search will be recomputed (~2.4 h).\n")
    } else cat("  no matching aligned cache in the reference run; the search will be recomputed (~2.4 h).\n")
  }
} else if (nzchar(REF_DIR)) {
  cat("  (reference run ", REF_DIR, " not found: no pre-flight comparison, aligned search will run)\n", sep = "")
}

section(paste("Canonical run:", SET$label))
cat("  input  :", normalizePath(SET$binary_dir), "\n")
cat("  output :", OUT_DIR, "\n")
cat("  aligned grid (frozen): theta_max =", CONFIG$aligned$theta_max, " t_max =", CONFIG$aligned$t_max, "\n")
cat("  fold_cache_version:", CONFIG$fold_cache_version, " export_to_manuscript:", CONFIG$export_to_manuscript, " workers:", CONFIG$n_workers, "\n")
cat("  methods: the 16 of the submitted pipeline + DT_L2_150_clamped (matched-transformation geometry baseline)\n\n")

t0 <- Sys.time()
results <- run_prl_suite_v2(CONFIG)
cat(sprintf("\nDone in %.1f minutes.\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

br <- read_csv(file.path(OUT_DIR, "aligned_boundary_rates.csv"), show_col_types = FALSE)
rt <- read_csv(file.path(OUT_DIR, "aligned_runtime.csv"), show_col_types = FALSE)
saveRDS(list(dataset = SET$label, binary_dir = SET$binary_dir, preprocess_manifest_created = pm$created,
             grid = FROZEN_GRID, aligned_key = rt$aligned_key, aligned_cache_reused = isTRUE(!rt$computed_this_run),
             shared_search_minutes = rt$shared_search_minutes, boundary_rates = br,
             reference_run = if (nzchar(REF_DIR)) REF_DIR else NA_character_, driver = "Code/run_canonical.R", created = Sys.time()),
        file.path(OUT_DIR, "run_record.rds"))
cat("Next:  options(fp.verify_a = \"", REF_DIR, "\", fp.verify_b = \"", OUT_DIR, "\"); source(\"Code/verify_run_vs_canonical.R\")\n",
    "       options(fp.bootstrap_dir = \"", OUT_DIR, "\"); source(\"Code/bootstrap_ci.R\")\n", sep = "")
