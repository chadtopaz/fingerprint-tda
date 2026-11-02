#!/usr/bin/env python3
"""
compare_sensitivity_summaries.py — exact scope gate for a regenerated sensitivity run (release step 5, 1 Nov 2026).

    python3 Code/sensitivity/compare_sensitivity_summaries.py --new <new root> --old <archived root> --conditions a,b,...

For summary/ci_perturbed_minus_clean.csv and summary/point_estimates_by_condition.csv of the NEW root:
  * the key set (condition, method, metric) [or (condition, method)] must be EXACTLY the Cartesian product of the
    listed conditions, the 15 pipeline variants and the 10 metrics (no extra, no missing, no duplicate row);
  * every row must equal, column by column, the archived row with the same key (numbers compared as the decimal
    strings written by the analysis; labels and flags as strings); the archived file may hold other conditions.
Fails with a non-zero exit on any difference. Nothing is written.
"""
import argparse, csv, sys

METHODS = ["DT_L2", "DT_L2_150", "GlobalBetti_L1", "GlobalPI_L2", "V23_full", "V23_globalOnly", "V23_localOnly", "V23_noTopK",
           "OT", "OT_noGeo", "OT_noDustbin", "OT_ridgeOnly", "OT_valleyOnly", "Fusion_eq", "Fusion"]
METRICS = ["AUC", "EER", "TAR@0.001", "TAR@0.002", "TAR@0.005", "TAR@0.01", "realisedFAR@0.001", "realisedFAR@0.002", "realisedFAR@0.005", "realisedFAR@0.01"]


def load(path, key):
    rows = list(csv.DictReader(open(path)))
    d = {}
    for r in rows:
        k = tuple(r[c] for c in key)
        if k in d:
            sys.exit("duplicate key %s in %s" % (k, path))
        d[k] = r
    return d, (rows[0].keys() if rows else [])


def compare(name, key, expected_keys, new_root, old_root):
    new, cols_new = load("%s/summary/%s" % (new_root, name), key)
    old, cols_old = load("%s/summary/%s" % (old_root, name), key)
    problems = []
    if set(new) != expected_keys:
        problems.append("%s: key set differs from the expected product: %d expected, %d present; missing %s; extra %s" % (
            name, len(expected_keys), len(new), sorted(expected_keys - set(new))[:3], sorted(set(new) - expected_keys)[:3]))
    if list(cols_new) != list(cols_old):
        problems.append("%s: column set differs: %s vs %s" % (name, list(cols_new), list(cols_old)))
    n_cells = 0
    for k in sorted(expected_keys):
        if k not in new or k not in old:
            if k not in old:
                problems.append("%s: archived file lacks %s" % (name, k))
            continue
        for c in cols_new:
            n_cells += 1
            if new[k][c] != old[k].get(c):
                problems.append("%s: %s column %s: new %r vs archived %r" % (name, k, c, new[k][c], old[k].get(c)))
                if len(problems) > 60:
                    break
    return problems, n_cells


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--new", required=True); ap.add_argument("--old", required=True); ap.add_argument("--conditions", required=True)
    A = ap.parse_args(); conds = [c.strip() for c in A.conditions.split(",")]
    exp3 = {(c, m, t) for c in conds for m in METHODS for t in METRICS}
    exp2 = {(c, m) for c in conds for m in METHODS}
    p1, n1 = compare("ci_perturbed_minus_clean.csv", ["condition", "method", "metric"], exp3, A.new, A.old)
    p2, n2 = compare("point_estimates_by_condition.csv", ["condition", "method"], exp2, A.new, A.old)
    problems = p1 + p2
    if problems:
        print("\n".join(problems[:60])); sys.exit("SENSITIVITY SUMMARY COMPARISON FAILED (%d problem(s))" % len(problems))
    print("sensitivity summaries identical to the archived ones on %s: %d + %d keys (exactly the expected %d + %d), %d cells compared" % (
        ", ".join(conds), len(exp3), len(exp2), len(conds) * len(METHODS) * len(METRICS), len(conds) * len(METHODS), n1 + n2))


if __name__ == "__main__":
    main()
