TRANSFER DIRECTORY - NOT A PIPELINE RUN
Written by Code/transfer_db1_to_db3.R on 2026-10-13 10:38:49.
DB1 -> DB3 transfer test, fold-matched: the fitted objects of DB1 fold k (Outputs_regen) applied to the DB3 fold-k test pairs (Outputs_db3).
Calibration: FROZEN - DB1 fold-k thresholds; cache s_train/y_train are the DB1 fold-k training scores (so the stored thresholds equal their training quantiles).
Nothing was refitted, retuned or re-standardised; model/best_w/zfit_* are the DB1 objects, copied verbatim.
The manifest is the DB3 replication manifest (ids, folds, pairs) plus a `transfer` element; the layout exists so that Code/bootstrap_ci.R can analyse it unchanged.
Training-free methods: transfer test scores are bit-identical to the replication's (asserted); only the thresholds differ in the frozen variant, nothing in the recal variant.
Design: Pattern Recognition/Revision/STEP_1.1_DB3_SPEC.md section 5.
