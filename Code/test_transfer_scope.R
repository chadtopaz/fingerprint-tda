############################################################
# TRANSFER SCOPE REGRESSION (release 1.1.0, 2 Nov 2026). Run from the project root, ≈ 5 s, nothing written:
#   source("Code/test_transfer_scope.R")
# The transfer experiment covers the 16 methods of the submitted pipeline; the _v2 runs carry a 17th, the clamped
# geometry control, which is outside it. Before this test the transfer script stopped on the _v2 runs.
#   S1  the real DB1/DB3 runs (fp.transfer_src / fp.transfer_tgt, defaults _v2): scope = 16 methods (10 training-free,
#       6 trained), outside = exactly the clamped control
#   S2  a 16-method pair (the clamped control absent): scope = 16, nothing outside
#   S3  a method of the scope missing from the runs -> stop naming it
#   S4  an unknown extra method -> stop naming it
#   S5  source and target method sets differ -> stop
#   S6  the archived transfer bootstrap (Outputs_db3_transfer, if present) carries exactly the 16 scope methods
############################################################
`%||%` <- function(x, y) if (!is.null(x)) x else y
source("Code/transfer_scope.R")
FAILS <- 0L
check <- function(cond, label) { ok <- isTRUE(cond); cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", label)); if (!ok) FAILS <<- FAILS + 1L; invisible(ok) }
err_of <- function(expr) tryCatch({ expr; NULL }, error = function(e) conditionMessage(e))
SRC <- getOption("fp.transfer_src", "Outputs_regen_v2"); TGT <- getOption("fp.transfer_tgt", "Outputs_db3_v2")
idx <- function(d) read.csv(file.path(d, "fold_cache_index.csv"), stringsAsFactors = FALSE)
v23 <- c("V23_full", "V23_noTopK", "V23_globalOnly", "V23_localOnly"); ot <- c("OT", "OT_noGeo", "OT_noDustbin", "OT_ridgeOnly", "OT_valleyOnly")

cat("S1 real runs\n")
if (file.exists(file.path(SRC, "fold_cache_index.csv")) && file.exists(file.path(TGT, "fold_cache_index.csv"))) {
  ms <- unique(idx(SRC)$method); mt <- unique(idx(TGT)$method)
  tcfg <- readRDS(file.path(TGT, "run_manifest.rds"))$config
  sc <- transfer_method_scope(ms, mt, tcfg$v23_variants, tcfg$ot_variants$method)
  check(length(sc$methods) == 16 && length(sc$training_free) == 10 && length(sc$trained) == 6, sprintf("scope of 16 methods (%d training-free, %d trained) from %s / %s", length(sc$training_free), length(sc$trained), SRC, TGT))
  check(identical(sc$outside, "DT_L2_150_clamped"), sprintf("outside the scope: %s", paste(sc$outside, collapse = ", ")))
} else cat("  (skipped: runs not present)\n")

base16 <- c("DT_L2", "DT_L2_150", "DT_L2_150_aligned", "GlobalBetti_L1", "GlobalPI_L2", ot, v23, "Fusion", "Fusion_eq")
cat("S2 sixteen-method pair\n")
sc <- transfer_method_scope(base16, base16, v23, ot); check(length(sc$methods) == 16 && length(sc$outside) == 0, "scope 16, nothing outside")
cat("S3 missing scope method\n")
m <- err_of(transfer_method_scope(setdiff(base16, "OT_noGeo"), setdiff(base16, "OT_noGeo"), v23, ot)); check(!is.null(m) && grepl("OT_noGeo", m), paste("stops:", substr(m %||% "", 1, 80)))
cat("S4 unknown extra method\n")
m <- err_of(transfer_method_scope(c(base16, "Mystery"), c(base16, "Mystery"), v23, ot)); check(!is.null(m) && grepl("Mystery", m), paste("stops:", substr(m %||% "", 1, 80)))
cat("S5 source/target differ\n")
m <- err_of(transfer_method_scope(base16, c(base16, "DT_L2_150_clamped"), v23, ot)); check(!is.null(m) && grepl("different method sets", m), paste("stops:", substr(m %||% "", 1, 80)))
cat("S6 archived transfer bootstrap\n")
if (file.exists("Outputs_db3_transfer/analysis/bootstrap_replicates.rds")) {
  B <- readRDS("Outputs_db3_transfer/analysis/bootstrap_replicates.rds"); check(setequal(B$methods, base16), sprintf("Outputs_db3_transfer bootstrap carries exactly the 16 scope methods (%d)", length(B$methods)))
} else cat("  (skipped)\n")
cat("\n", if (FAILS == 0L) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", FAILS), "\n", sep = "")
if (FAILS > 0L) stop(FAILS, " transfer scope check(s) failed", call. = FALSE)
