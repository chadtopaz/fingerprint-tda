############################################################
# STEP 1.5 (D2/D3) — FULL PIPELINE ON THE CUBICAL COMPLEX, canonical DB1 inputs, as a
# supplementary side-by-side with the reported (triangulated) run in Outputs_regen.
#
# Everything is the canonical DB1 run's configuration (its manifest config) with ONE change:
# complex_type = "cubical" (GUDHI Bitmap_cubical_complex via Code_cubical/cubical_ph.R).
# Same images, identities, folds, seeds, parameters, frozen aligned grid, Fusion_eq default.
#
# Pre-registered rule (STEP_1.5_CUBICAL_SPEC.md): the REPORTED method stays triangulated
# whatever this run shows; the comparison is descriptive and goes to the supplement.
#
# The aligned geometry baseline does not involve persistence (it works on the binary images'
# distance fields, bank field dt_vec_150, identical under both constructions), and its pair
# cache is keyed by (grid, image md5s) only; the canonical cache is therefore copied in and
# re-validated by the engine's own key/files check instead of recomputing the 2.4-h search.
#
# Run from the PROJECT ROOT in RStudio (~1.5 h: cubical bank a few minutes, fold loop ~80 min),
# after Code_cubical/test_cubical_integration.R has passed:
#   source("Code_cubical/run_cubical_db1.R")
# Writes only into Outputs_cubical/ (must be new or empty). Never touches Outputs_regen/.
############################################################

options(fp.suppress_main = TRUE)
source("Code/fingerprintTDA.R")
source("Code_cubical/cubical_ph.R"); PY <- configure_cubical_python()
source("Code/frozen_grid.R")

# Source run (triangulated) and output directory. DB1 by default; for the DB3 follow-up (spec section 10):
#   options(fp.cubical_src = "Outputs_db3", fp.cubical_out = "Outputs_cubical_db3"); source("Code_cubical/run_cubical_db1.R")
SRC_DIR <- getOption("fp.cubical_src", "Outputs_regen_v2")
CUB_DIR <- getOption("fp.cubical_out", "Outputs_cubical")
if (dir.exists(CUB_DIR) && length(list.files(CUB_DIR, recursive = TRUE, all.files = TRUE, no.. = TRUE)) > 0)
  stop(CUB_DIR, " already contains files; this must be a cold run into a new or empty directory.", call. = FALSE)
stopifnot(isTRUE(cubical_ph_selftest(verbose = FALSE)))
cat("Cubical backend self-test passed (", reticulate::py_config()$python, ").\n", sep = "")

# ---- configuration: the canonical run's own config with the complex switched
src <- readRDS(file.path(SRC_DIR, "run_manifest.rds"))
stopifnot(identical(src$config$complex_type %||% "triangulated", "triangulated"))   # the canonical run's effective construction
cfg <- src$config
cfg$complex_type         <- "cubical"
cfg$output_dir           <- CUB_DIR
cfg$export_to_manuscript <- FALSE
cfg$n_workers            <- CONFIG$n_workers
stopifnot(cfg$binary_dir %in% c("Data/FVC2000_binary_regen", "Data/FVC2002_DB3_binary"), identical(apply_frozen_grid(cfg)$aligned, cfg$aligned))
# the file's CONFIG must not have drifted from the canonical run either
fa <- apply_frozen_grid(CONFIG); fb <- cfg
stopifnot(identical(src$config$polarity %||% "auto", "auto"))                                       # the source run used the default polarity rule
for (fld in c("binary_dir", "output_dir", "export_to_manuscript", "complex_type", "polarity")) { fa[[fld]] <- NULL; fb[[fld]] <- NULL }
stopifnot(identical(fa, fb))
paths <- get_paths(cfg)
stopifnot(identical(basename(paths), src$files), identical(unname(tools::md5sum(paths)), src$md5))
cat("Configuration = the source run (", SRC_DIR, ", ", cfg$binary_dir, ") with complex_type = 'cubical'; 800 files and md5s identical to its manifest.\n", sep = "")
cat("  bank key (cubical):", substr(bank_cache_key(cfg), 1, 12), " vs canonical", substr(bank_cache_key(src$config), 1, 12),
    "| fold key:", substr(fold_cache_key(cfg), 1, 8), "vs", substr(fold_cache_key(src$config), 1, 8), "\n")

# ---- aligned pair cache: same key as the canonical run (grid + images), copied in, re-validated by the engine
al_key <- aligned_cache_key(cfg$aligned, paths)
al_src <- file.path(SRC_DIR, "cache", paste0("aligned_", substr(al_key, 1, 12), ".rds"))
stopifnot(file.exists(al_src))
al_obj <- readRDS(al_src); stopifnot(identical(al_obj$key, al_key), identical(al_obj$files, basename(paths)))
dir.create(file.path(CUB_DIR, "cache"), recursive = TRUE, showWarnings = FALSE)
file.copy(al_src, file.path(CUB_DIR, "cache", basename(al_src)))
cat("Aligned pair cache", basename(al_src), "copied from", SRC_DIR, "(persistence-independent; key verified).\n\n")

t0 <- Sys.time()
results_cub <- run_prl_suite_v2(cfg)
cat(sprintf("\nDone in %.1f minutes.\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
rt <- read_csv(file.path(CUB_DIR, "aligned_runtime.csv"), show_col_types = FALSE)
stopifnot(identical(rt$aligned_key, al_key), !isTRUE(rt$computed_this_run))
saveRDS(list(complex_type = "cubical", source_run = SRC_DIR, grid = FROZEN_GRID, aligned_key = al_key, aligned_cache_reused = TRUE,
             python = reticulate::py_config()$python, gudhi = as.character(gudhi_handle()$`__version__`), created = Sys.time()),
        file.path(CUB_DIR, "run_record.rds"))
cat("Next:  options(fp.bootstrap_dir = \"", CUB_DIR, "\"); source(\"Code/bootstrap_ci.R\")\n", sep = "")
cat("Then:  options(fp.compare_a = \"", SRC_DIR, "\", fp.compare_b = \"", CUB_DIR, "\"); source(\"Code_cubical/compare_constructions.R\")\n", sep = "")
