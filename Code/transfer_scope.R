############################################################
# TRANSFER METHOD SCOPE (release 1.1.0, 2 Nov 2026). The DB1 -> DB3 transfer experiment of the paper covers the 16
# methods of the submitted pipeline: 10 training-free (the three DT L2 baselines incl. the aligned one, the two global
# summaries, the five LPOT variants) and 6 trained (the four TopoLR variants, the tuned and the equal-weight fusion).
# The clamped geometry control DT_L2_150_clamped, added in the 31 Oct audit, is training-free and OUTSIDE this
# experiment: its DB3 value is its replication value, reported in the DB3 tables; nothing of it is transferred.
# Both transfer scripts and Code/test_transfer_scope.R use this one definition.
TRANSFER_SCOPE_OUTSIDE <- "DT_L2_150_clamped"
transfer_method_scope <- function(methods_src, methods_tgt, v23_variants, ot_methods) {
  training_free <- c("DT_L2", "DT_L2_150", "DT_L2_150_aligned", "GlobalBetti_L1", "GlobalPI_L2", ot_methods)
  trained       <- c(v23_variants, "Fusion", "Fusion_eq")
  scope <- c(training_free, trained)
  if (!setequal(methods_src, methods_tgt))
    stop("source and target runs carry different method sets: only in source ", paste(setdiff(methods_src, methods_tgt), collapse = ", "),
         "; only in target ", paste(setdiff(methods_tgt, methods_src), collapse = ", "), call. = FALSE)
  missing <- setdiff(scope, methods_tgt)
  if (length(missing)) stop("transfer scope: method(s) missing from the runs: ", paste(missing, collapse = ", "), call. = FALSE)
  outside <- setdiff(methods_tgt, scope)
  unknown <- setdiff(outside, TRANSFER_SCOPE_OUTSIDE)
  if (length(unknown)) stop("method(s) in the runs that are neither in the transfer scope nor the declared clamped control: ", paste(unknown, collapse = ", "), call. = FALSE)
  if (length(scope) != 16L) stop("transfer scope must have 16 methods, has ", length(scope), call. = FALSE)
  list(methods = scope, training_free = training_free, trained = trained, outside = outside)
}
