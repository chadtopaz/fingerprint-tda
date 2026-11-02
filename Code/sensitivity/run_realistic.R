############################################################
# Tier R driver (R4-10): the clean-trained DB1 verifier applied to the ACQUISITION DEGRADATIONS
# (Data/FVC2000_binary_degraded/<cond>: raw-domain degradations at the calibrated strength, passed through the
# canonical enhancer), with the CLEAN-ENROLLMENT / DEGRADED-PROBE pairing as primary.
#
# Run from the PROJECT ROOT (~40 min scoring + ~6 min bootstraps for the four conditions):
#   source("Code/sensitivity/run_realistic.R")
# Secondary (both images of every pair degraded), only if QA asks for it:
#   options(fp.realistic_pairing = "both"); source("Code/sensitivity/run_realistic.R")
# Writes Outputs_sensitivity_realistic/ (pairing "probe") or Outputs_sensitivity_realistic_both/ (pairing "both"):
# one run-shaped directory per condition + analysis/ (bootstrap_ci.R unchanged) + summary/ with the
# degraded - clean paired differences on the clean run's finger resamples. Never writes to the source run
# or to Outputs_sensitivity (the Tier B record).
############################################################
pairing <- getOption("fp.realistic_pairing", "probe")
stopifnot(pairing %in% c("probe", "both"))
conds <- c("dry_typ", "wet_typ", "noise_typ", "blur_typ")      # the calibrated level of each family
out_root <- if (pairing == "probe") "Outputs_sensitivity_realistic" else "Outputs_sensitivity_realistic_both"
old <- options(fp.sens_conditions = file.path("Data/FVC2000_binary_degraded", conds), fp.sens_out = out_root,
               fp.sens_pairing = pairing, fp.sens_requested = conds)
cat("\nTier R — realistic acquisition degradations; pairing '", pairing, "'; output ", out_root, "\n", sep = "")
tryCatch({
  source("Code/sensitivity/score_perturbed.R")                 # provenance checks, feature banks, frozen-model scoring (stops on any failure)
  source("Code/sensitivity/run_sensitivity_bootstraps.R")      # bootstraps + degraded - clean paired differences + headline
}, finally = options(old))                                     # the options never leak into a later plain source() of either script
cat("\nTier R complete: ", out_root, "/summary/ci_perturbed_minus_clean.csv\n", sep = "")
