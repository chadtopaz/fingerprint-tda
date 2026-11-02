#!/usr/bin/env python3
"""
build_results_bundle.py — the released Results/ bundle: the aggregate, pixel-free files from which every table of
the article and supplement is formatted (Code/make_tables.py --source-root Results), plus the per-run summary
files, the calibration record of the degradation experiment and the preprocessing manifests. Each file is copied
to Results/<its project-relative path>, so every script's path contract holds unchanged under --source-root.

    python3 Code/build_results_bundle.py build  [--target Results] [--final]   # copy the allowlist from the project; write INVENTORY.json
    python3 Code/build_results_bundle.py check  [--target Results] [--final]   # verify a bundle against its INVENTORY.json and the allowlist
    python3 Code/build_results_bundle.py list                                  # print the allowlist

--final adds the CONTENT SCOPE GATE of the public release (without it a build is labelled "preparatory" in its inventory
and is never release acceptance): the Tier B summary holds exactly the six enhancement variants x 15 methods x 10 metrics
(900 rows) and the Tier R summary exactly the four calibrated degradations (600 rows), with unique keys; the point-estimate
tables exactly 6 x 15 and 4 x 15; calibration.json's chosen set is exactly the four conditions at level "typical"; the
NFIQ record holds exactly the clean directory and the four conditions; the verify and distribution records exactly the
four; the four degrade manifests record the hash of the shipped calibration.json. A final build fails if any of this
does not hold (the hashed calibration metadata is never hand-filtered to make it hold: the chain is regenerated).

Policy (release step 5, 1 Nov 2026): nothing derived from the licensed images at the pixel or pair level is
shipped — no image, binary map, degraded image, feature bank, fold cache, pair-score file, bootstrap replicate
file or NBIS template/score file. Only files of the kinds below, each of which carries estimates, counts, settings,
hashes or per-image scalar descriptors. The inventory records each file's sha256, size and source run, and
`check` fails on any missing, altered or unlisted file and on any forbidden file kind anywhere under the target.
"""
import argparse, hashlib, json, os, shutil, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ANALYSIS = ["analysis/ci_metrics.csv", "analysis/ci_paired_differences.csv", "analysis/ci_metrics_per_fold.csv",
            "analysis/acceptance_counts.csv", "analysis/endpoint_stability.csv", "analysis/bootstrap_settings.txt",
            "analysis/fusion_eq_fold_metrics.csv", "analysis/fusion_weights.csv", "analysis/selfcheck_point_estimates.csv"]
RUN = ["fold_metrics_auc_eer.csv", "fold_metrics_tar.csv", "summary_auc_eer_mean_sd.csv", "summary_tar_mean_sd.csv",
       "fusion_weights_per_fold.csv", "aligned_boundary_rates.csv", "aligned_runtime.csv", "fold_runtimes_seconds.csv",
       "summary_runtime_mean_sd.csv", "fold_cache_index.csv"]
REF = ["fold_metrics_auc_eer.csv", "fold_metrics_tar.csv", "reference_fold_metrics.csv", "fold_cache_index.csv", "REFERENCE_README.txt"]
TIER_B = ["filter_m2", "filter_m4", "scale_055", "scale_075", "segment_020", "segment_040"]
TIER_R = ["dry_typ", "wet_typ", "noise_typ", "blur_typ"]

# (project-relative path, source description); every entry is REQUIRED unless it ends with "?" (optional)
ALLOWLIST = []
def add(paths, source):
    for p in paths:
        ALLOWLIST.append((p, source))

for d, src in (("Outputs_regen_v2", "canonical DB1 run (17 methods)"), ("Outputs_db3_v2", "canonical DB3 run (17 methods)")):
    add([d + "/" + f for f in RUN + ANALYSIS], src)
for d, src in (("Outputs_reference_db1_v2", "DB1 run + NBIS reference matcher + two NBIS fusions"), ("Outputs_reference_db3_v2", "DB3 run + NBIS reference matcher + two NBIS fusions")):
    add([d + "/" + f for f in REF + ANALYSIS], src)
add(["Outputs_regen/summary_runtime_mean_sd.csv", "Outputs_regen/fold_runtimes_seconds.csv", "Outputs_regen/aligned_runtime.csv"], "canonical DB1 run that measured the runtimes and performed the alignment search (bit-identical scores to Outputs_regen_v2)")
add(["Outputs_sensitivity/summary/ci_perturbed_minus_clean.csv", "Outputs_sensitivity/summary/point_estimates_by_condition.csv"], "fixed-model sensitivity, six enhancement-parameter variants (both images re-processed)")
add(["Outputs_sensitivity/%s/sensitivity_fold_metrics.csv" % c for c in TIER_B] + ["Outputs_sensitivity/%s/fold_metrics_auc_eer.csv" % c for c in TIER_B] + ["Outputs_sensitivity/%s/fold_metrics_tar.csv" % c for c in TIER_B], "fixed-model sensitivity, per condition")
add(["Outputs_sensitivity/%s/" % c + f for c in TIER_B for f in ANALYSIS], "fixed-model sensitivity, per-condition bootstrap")
add(["Outputs_sensitivity_realistic/summary/ci_perturbed_minus_clean.csv", "Outputs_sensitivity_realistic/summary/point_estimates_by_condition.csv"], "fixed-model sensitivity, four calibrated degradations of the probe")
add(["Outputs_sensitivity_realistic/%s/sensitivity_fold_metrics.csv" % c for c in TIER_R] + ["Outputs_sensitivity_realistic/%s/fold_metrics_auc_eer.csv" % c for c in TIER_R] + ["Outputs_sensitivity_realistic/%s/fold_metrics_tar.csv" % c for c in TIER_R], "fixed-model sensitivity, per condition")
add(["Outputs_sensitivity_realistic/%s/" % c + f for c in TIER_R for f in ANALYSIS], "fixed-model sensitivity, per-condition bootstrap")
for d, src in (("Outputs_cubical", "cubical construction, DB1"), ("Outputs_cubical_db3", "cubical construction, DB3")):
    add([d + "/analysis/ci_cubical_minus_triangulated.csv", d + "/analysis/diagram_comparison.csv", d + "/analysis/diagram_comparison_summary.txt",
         d + "/fold_metrics_auc_eer.csv", d + "/fold_metrics_tar.csv"] + [d + "/" + f for f in ANALYSIS], src)
add(["Outputs_db3_transfer/analysis/ci_transfer_minus_replication.csv", "Outputs_db3_transfer/transfer_fold_metrics.csv", "Outputs_db3_transfer/fold_metrics_auc_eer.csv", "Outputs_db3_transfer/fold_metrics_tar.csv", "Outputs_db3_transfer/TRANSFER_README.txt"] + ["Outputs_db3_transfer/" + f for f in ANALYSIS], "DB1-fitted models applied to DB3, frozen thresholds")
add(["Outputs_db3_transfer_recal/fold_metrics_auc_eer.csv", "Outputs_db3_transfer_recal/fold_metrics_tar.csv"] + ["Outputs_db3_transfer_recal/" + f for f in ANALYSIS], "DB1-fitted models applied to DB3, thresholds re-estimated on DB3 training impostors")
add(["Data/FVC2000_raw_degraded/" + f for f in ("calibration.json", "calibration_benchmark.json", "calibration_scan.json", "calibration_distribution.json", "calibration_nfiq.json", "calibration_verify.json", "calibration_stats_clean.csv")], "calibration of the degradation strengths (descriptor statistics; no pixels)")
add(["Data/FVC2000_raw_degraded/%s/degrade_manifest.json" % c for c in TIER_R], "degradation manifests (per-image hashes, seeds and scalar summaries; required by make_tables.py)")
add(["Data/FVC2000_binary_degraded/%s/preprocess_manifest.json" % c for c in TIER_R] + ["Data/FVC2000_binary_variant/%s/preprocess_manifest.json" % c for c in TIER_B], "preprocessing manifests of the sensitivity conditions")
add(["Data/FVC2000_binary_regen/preprocess_manifest.json", "Data/FVC2002_DB3_binary/preprocess_manifest.json"], "preprocessing manifests of the canonical binaries (enhancer version, per-image hashes): the check for a reproduced preprocessing")
add(["Outputs_figures/fig_roc_pooled_auc.csv"], "pooled AUCs of the ROC figure legend")

METHODS = ["DT_L2", "DT_L2_150", "GlobalBetti_L1", "GlobalPI_L2", "V23_full", "V23_globalOnly", "V23_localOnly", "V23_noTopK",
           "OT", "OT_noGeo", "OT_noDustbin", "OT_ridgeOnly", "OT_valleyOnly", "Fusion_eq", "Fusion"]
METRICS = ["AUC", "EER", "TAR@0.001", "TAR@0.002", "TAR@0.005", "TAR@0.01", "realisedFAR@0.001", "realisedFAR@0.002", "realisedFAR@0.005", "realisedFAR@0.01"]


def scope_problems(root):
    """The content scope gate of the public release; returns a list of problems (empty = pass). root = where the files sit."""
    import csv
    problems = []
    def keys(path, cols):
        rows = list(csv.DictReader(open(os.path.join(root, path))))
        ks = [tuple(r[c] for c in cols) for r in rows]
        if len(set(ks)) != len(ks):
            problems.append("%s: duplicate keys" % path)
        return set(ks), rows
    for summ, conds, tag in (("Outputs_sensitivity/summary", TIER_B, "Tier B"), ("Outputs_sensitivity_realistic/summary", TIER_R, "Tier R")):
        exp3 = {(c, m, t) for c in conds for m in METHODS for t in METRICS}; exp2 = {(c, m) for c in conds for m in METHODS}
        k3, _ = keys(summ + "/ci_perturbed_minus_clean.csv", ["condition", "method", "metric"])
        if k3 != exp3:
            problems.append("%s ci_perturbed_minus_clean.csv: key set != %d conditions x 15 methods x 10 metrics (present %d, expected %d; extra %s)" % (tag, len(conds), len(k3), len(exp3), sorted(k3 - exp3)[:3]))
        k2, _ = keys(summ + "/point_estimates_by_condition.csv", ["condition", "method"])
        if k2 != exp2:
            problems.append("%s point_estimates_by_condition.csv: key set != %d x 15 (present %d; extra %s)" % (tag, len(conds), len(k2), sorted(k2 - exp2)[:3]))
    cal_path = os.path.join(root, "Data/FVC2000_raw_degraded/calibration.json")
    cal = json.load(open(cal_path))
    if sorted(cal["chosen"]) != sorted(TIER_R) or any(cal["chosen"][c].get("level") != "typical" for c in TIER_R if c in cal["chosen"]):
        problems.append("calibration.json: chosen set is %s (expected exactly the four _typ conditions at level typical)" % sorted(cal["chosen"]))
    cal_sha = sha256(cal_path)
    for c in TIER_R:
        dm = json.load(open(os.path.join(root, "Data/FVC2000_raw_degraded/%s/degrade_manifest.json" % c)))
        if dm.get("calibration_sha256") != cal_sha:
            problems.append("%s degrade manifest records a calibration hash that is not the shipped calibration.json's" % c)
    nfiq = json.load(open(os.path.join(root, "Data/FVC2000_raw_degraded/calibration_nfiq.json")))
    if sorted(nfiq) != sorted(["FVC2000_raw"] + TIER_R):
        problems.append("calibration_nfiq.json: entries %s (expected the clean directory and the four conditions)" % sorted(nfiq))
    ver = json.load(open(os.path.join(root, "Data/FVC2000_raw_degraded/calibration_verify.json")))
    if sorted(ver) != sorted(TIER_R):
        problems.append("calibration_verify.json: entries %s" % sorted(ver))
    dist = json.load(open(os.path.join(root, "Data/FVC2000_raw_degraded/calibration_distribution.json")))
    if sorted(dist["conditions"]) != sorted(TIER_R):
        problems.append("calibration_distribution.json: conditions %s" % sorted(dist["conditions"]))
    return problems


FORBIDDEN_EXT = (".png", ".tif", ".tiff", ".jpg", ".jpeg", ".rds", ".rda", ".RData", ".npy", ".npz", ".xyt", ".min", ".qm", ".dm", ".hcm")
FORBIDDEN_NAMES = ("scores_raw.csv", "aligned_pair_diagnostics.csv", "bootstrap_replicates.rds")
FORBIDDEN_DIRS = ("cache",)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def forbidden(rel):
    base = os.path.basename(rel)
    return base in FORBIDDEN_NAMES or base.lower().endswith(FORBIDDEN_EXT) or any(part in FORBIDDEN_DIRS for part in rel.split("/")[:-1])


def build(target, final=False):
    missing = [p for p, _ in ALLOWLIST if not os.path.exists(os.path.join(ROOT, p))]
    if missing:
        sys.exit("cannot build: %d allowlisted file(s) missing from the project:\n  " % len(missing) + "\n  ".join(missing))
    bad = [p for p, _ in ALLOWLIST if forbidden(p)]
    if bad:
        sys.exit("allowlist contains forbidden kinds: %s" % bad)
    if os.path.exists(target) and os.listdir(target):
        sys.exit("target %s exists and is not empty; remove it first (the bundle is rebuilt from scratch)" % target)
    if final:
        sp = scope_problems(ROOT)
        if sp:
            sys.exit("final build refused, the project's files are outside the release scope:\n  " + "\n  ".join(sp))
    inv = []
    for rel, src in ALLOWLIST:
        s = os.path.join(ROOT, rel); t = os.path.join(target, rel)
        os.makedirs(os.path.dirname(t), exist_ok=True); shutil.copyfile(s, t)
        inv.append(dict(path=rel, sha256=sha256(t), bytes=os.path.getsize(t), source=src))
    meta = dict(created=time.strftime("%Y-%m-%dT%H:%M:%S%z"), scope="final" if final else "preparatory (not release acceptance)", n_files=len(inv),
                policy=__doc__.split("Policy")[1].strip(), files=inv)
    json.dump(meta, open(os.path.join(target, "INVENTORY.json"), "w"), indent=1)
    print("built %s (%s): %d files, %.1f MB" % (target, meta["scope"], len(inv), sum(f["bytes"] for f in inv) / 1e6))


def check(target, final=False):
    inv_path = os.path.join(target, "INVENTORY.json")
    if not os.path.exists(inv_path):
        sys.exit("no INVENTORY.json under %s" % target)
    inv = json.load(open(inv_path))
    listed = {f["path"]: f for f in inv["files"]}
    allow = {p for p, _ in ALLOWLIST}
    problems = []
    if set(listed) != allow:
        problems.append("inventory paths differ from the allowlist: missing %s, extra %s" % (sorted(allow - set(listed))[:5], sorted(set(listed) - allow)[:5]))
    present = set()
    for dp, dn, fn in os.walk(target):
        for f in fn:
            rel = os.path.relpath(os.path.join(dp, f), target).replace(os.sep, "/")
            if rel == "INVENTORY.json":
                continue
            present.add(rel)
            if forbidden(rel):
                problems.append("forbidden file kind in the bundle: " + rel)
            if rel not in listed:
                problems.append("file not in the inventory: " + rel)
            elif sha256(os.path.join(dp, f)) != listed[rel]["sha256"]:
                problems.append("hash differs from the inventory: " + rel)
    for rel in listed:
        if rel not in present:
            problems.append("inventoried file missing: " + rel)
    if final:
        if inv.get("scope") != "final":
            problems.append("inventory scope is %r, not a final build" % inv.get("scope"))
        problems += scope_problems(target)
    if problems:
        print("\n".join(problems[:50])); sys.exit("bundle check FAILED (%d problem(s))" % len(problems))
    print("bundle check passed%s: %d files, all hashes match, no forbidden kind, inventory == allowlist" % (" (FINAL scope gate included)" if final else " (preparatory; no scope gate)", len(listed)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("mode", choices=["build", "check", "list"]); ap.add_argument("--target", default=os.path.join(ROOT, "Results"))
    ap.add_argument("--final", action="store_true", help="apply the public-release content scope gate")
    A = ap.parse_args()
    if A.mode == "list":
        for p, s in ALLOWLIST:
            print("%-90s  %s" % (p, s))
    elif A.mode == "build":
        build(A.target, A.final)
    else:
        check(A.target, A.final)
