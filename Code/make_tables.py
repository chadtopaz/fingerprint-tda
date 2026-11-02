#!/usr/bin/env python3
"""
Generate every results table of the revised manuscript and supplement as LaTeX
fragments, from the frozen analysis CSVs (bootstrap_ci.R outputs).

Run from the project root:
    python3 Code/make_tables.py                          # sources under the project root (the run directories)
    python3 Code/make_tables.py --source-root Results    # sources under Results/ (the released aggregate bundle; same relative paths)
    python3 Code/make_tables.py --out <dir>              # write the fragments elsewhere (default Pattern Recognition/tables)
Writes  Pattern Recognition/tables/*.tex  (included by PR_main.tex / PR_supplement.tex
via \\input).  Nothing is computed here; the script only formats numbers that the
analysis scripts wrote. The released Results/ bundle (Code/build_results_bundle.py) holds every
input file below at its project-relative path, so a clean checkout reproduces every table
without the licensed images or any run directory.

Inputs (all produced by Code/bootstrap_ci.R on the corresponding run directories; see SRC below):
  Outputs_reference_db1_v2/analysis/ci_metrics.csv   DB1: 17 pipeline methods + NBIS + two NBIS fusions
  Outputs_reference_db3_v2/analysis/ci_metrics.csv   DB3: same
  Outputs_regen_v2 / Outputs_db3_v2: analysis/ci_paired_differences.csv, fusion_weights_per_fold.csv (tuned fusion)
  Outputs_regen/summary_runtime_mean_sd.csv          runtime per fold (the DB1 run that performed the alignment search)
  Outputs_sensitivity/summary/ci_perturbed_minus_clean.csv            Tier B (parameter variants)
  Outputs_sensitivity_realistic/summary/ci_perturbed_minus_clean.csv  calibrated degradations
  Outputs_cubical*, Outputs_db3_transfer                              construction comparison, transfer
  Data/FVC2000_raw_degraded/calibration*.json                         degradation calibration
The *_v2 run directories (31 Oct 2026) reproduce the canonical runs bit for bit on every shared object
(Code/verify_run_vs_canonical.R, Code/verify_analysis_vs_canonical.R) and add the matched-transformation
geometry baseline DT_L2_150_clamped. Runtime stays with the canonical DB1 run, which measured the 16 original
methods and the alignment search under the documented conditions.
Intervals: finger-level multiplicity-weighted bootstrap, B = 2000, pointwise 95 %,
thresholds fixed at the training-calibrated values (see Code/bootstrap_ci.R).

Formatting policy. Every number is formatted from the stored floating-point estimate with Python's %-formatting
(correctly rounded from the binary value, so an estimate that sits on a decimal .x5 boundary, e.g. a multiple of
1/2800, resolves according to its binary representation, not by a decimal half-even rule). Differences and
interval endpoints are set in math mode so that minus signs are typographic minuses. A difference whose
magnitude rounds to zero is printed unsigned ("0.000"); an endpoint of an interval that EXCLUDES zero but would
print as +-0.000 is printed with the fewest additional decimals that show its sign (e.g. "-0.00003").
"""
import argparse, csv, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_ap = argparse.ArgumentParser(description="format the manuscript tables from the stored analysis files")
_ap.add_argument("--source-root", default=ROOT, help="directory under which the source files sit at their project-relative paths (default: the project root; the released bundle: Results)")
_ap.add_argument("--out", default=os.path.join(ROOT, "Pattern Recognition", "tables"), help="output directory for the .tex fragments")
_args = _ap.parse_args()
SRC_ROOT = os.path.abspath(_args.source_root)
OUT = os.path.abspath(_args.out)
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- sources (run directories)
SRC = {
    "db1": "Outputs_reference_db1_v2",      # pipeline methods + NBIS + NBIS fusions, DB1 (bootstrap analysis)
    "db3": "Outputs_reference_db3_v2",      # same, DB3
    "run_db1": "Outputs_regen_v2",          # canonical DB1 pipeline run: paired differences, tuned fusion weights
    "run_db3": "Outputs_db3_v2",            # canonical DB3 pipeline run
    "runtime": "Outputs_regen",             # per-fold timings (the run that performed the alignment search)
    "sens_B": "Outputs_sensitivity",        # preprocessing-parameter variants
    "sens_R": "Outputs_sensitivity_realistic",   # calibrated degradations of the probe
    "cubical_db1": "Outputs_cubical", "cubical_db3": "Outputs_cubical_db3",
    "transfer": "Outputs_db3_transfer",
    "calibration": os.path.join("Data", "FVC2000_raw_degraded"),
}
# A source directory absent under the source root is taken from the released aggregate bundle, Results/<same path>,
# when that exists there: a reproduction that repeated the runs but not, say, the cubical comparison, or that could not
# repeat the frozen runtime measurement of the canonical DB1 run (Outputs_regen/summary_runtime_mean_sd.csv, kept as the
# reported timing), still formats every table; each substitution is printed.
RESULTS_FALLBACK = os.path.join(SRC_ROOT, "Results")
for k, v in list(SRC.items()):
    if not os.path.isdir(os.path.join(SRC_ROOT, v)):
        if os.path.isdir(os.path.join(RESULTS_FALLBACK, v)):
            SRC[k] = os.path.join("Results", v); print("source %s: %s not present, using the released Results/%s" % (k, v, v))
        else:
            sys.exit("missing source directory for %s: %s (under %s; also not under %s)" % (k, v, SRC_ROOT, RESULTS_FALLBACK))

# ---------------------------------------------------------------- display names
NAME = {
    "DT_L2": r"DT $L^2$ ($50^2$)",
    "DT_L2_150": r"DT $L^2$ ($150^2$)",
    "DT_L2_150_clamped": r"DT $L^2$ ($150^2$, clamped)",
    "DT_L2_150_aligned": r"aligned DT $L^2$ ($150^2$)",
    "GlobalBetti_L1": r"Betti $L^1$",
    "GlobalPI_L2": r"PI $L^2$",
    "V23_full": r"TopoLR",
    "V23_globalOnly": r"TopoLR, global only",
    "V23_localOnly": r"TopoLR, local only",
    "V23_noTopK": r"TopoLR, no top-$K$",
    "OT": r"LPOT",
    "OT_noGeo": r"LPOT, no geometry",
    "OT_noDustbin": r"LPOT, no dustbin",
    "OT_ridgeOnly": r"LPOT, ridge only",
    "OT_valleyOnly": r"LPOT, valley only",
    "Fusion_eq": r"Fusion (equal weight)",
    "Fusion": r"Fusion (tuned weight)",
    "NBIS_gray": r"NBIS (grayscale)",
    "NBIS_bin": r"NBIS (binary)",
    "Fus_NBISgray_TopoLR": r"NBIS $+$ TopoLR",
    "Fus_NBISgray_FusionEq": r"NBIS $+$ Fusion",
}
MAIN_NAME = dict(NAME, **{"DT_L2_150": r"DT $L^2$", "DT_L2_150_clamped": r"DT $L^2$, clamped", "DT_L2_150_aligned": r"DT $L^2$, aligned"})
MAIN_DISCRIM = ["DT_L2_150", "DT_L2_150_clamped", "DT_L2_150_aligned", "GlobalBetti_L1", "GlobalPI_L2", "V23_full", "OT",
                "Fusion_eq", "NBIS_gray", "Fus_NBISgray_TopoLR"]
MAIN_OPER = ["GlobalBetti_L1", "V23_full", "OT", "Fusion_eq", "NBIS_gray", "Fus_NBISgray_TopoLR"]
SUPP_ALL = ["DT_L2", "DT_L2_150", "DT_L2_150_clamped", "DT_L2_150_aligned", "GlobalBetti_L1", "GlobalPI_L2",
            "V23_full", "V23_globalOnly", "V23_localOnly", "V23_noTopK",
            "OT", "OT_noGeo", "OT_noDustbin", "OT_ridgeOnly", "OT_valleyOnly",
            "Fusion_eq", "Fusion", "NBIS_gray", "NBIS_bin", "Fus_NBISgray_TopoLR", "Fus_NBISgray_FusionEq"]
FARS = ["0.001", "0.002", "0.005", "0.01"]
FARLAB = {"0.001": r"10^{-3}", "0.002": r"2{\times}10^{-3}", "0.005": r"5{\times}10^{-3}", "0.01": r"10^{-2}"}


def read_ci(path):
    d = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            d[(r["method"], r["metric"])] = r
    return d


def f3(x):
    return "%.3f" % float(x)


def p1(x):
    return "%.1f" % (100 * float(x))


def p2(x):
    return "%.2f" % (100 * float(x))


def cell(d, m, metric, fmt):
    r = d.get((m, metric))
    if r is None:
        return "---"
    return "%s [%s, %s]" % (fmt(r["point"]), fmt(r["ci_low"]), fmt(r["ci_high"]))


def fmt_signed(v, nd, show_sign_if_tiny=False):
    """Signed value with nd decimals, in math mode. A value whose magnitude rounds to zero prints unsigned
    ("0.000"), unless show_sign_if_tiny (an endpoint of an interval that excludes zero): then the number of
    decimals grows until the printed value is nonzero, so the sign of the endpoint is visible."""
    t = "%+.*f" % (nd, v)
    if float(t) == 0.0:
        if show_sign_if_tiny and v != 0.0:
            k = nd
            while float(t) == 0.0 and k < nd + 8:
                k += 1
                t = "%+.*f" % (k, v)
        else:
            t = "%.*f" % (nd, 0.0)
    return "$%s$" % t


def diff_cell(r, s=1, nd=3):
    """difference [ci_low, ci_high] of a paired-difference row, scaled by s (100 for percentage points)."""
    d, lo, hi = s * float(r["difference"]), s * float(r["ci_low"]), s * float(r["ci_high"])
    excludes_zero = lo > 0 or hi < 0
    return "%s [%s, %s]" % (fmt_signed(d, nd), fmt_signed(lo, nd, excludes_zero), fmt_signed(hi, nd, excludes_zero))


def write(name, text, colsep="4pt"):
    with open(os.path.join(OUT, name), "w") as f:
        f.write("%% generated by Code/make_tables.py -- do not edit by hand\n\\setlength{\\tabcolsep}{%s}\n" % colsep + text)
    print("wrote " + os.path.relpath(os.path.join(OUT, name), ROOT))


db1 = read_ci(os.path.join(SRC_ROOT, SRC["db1"], "analysis", "ci_metrics.csv"))
db3 = read_ci(os.path.join(SRC_ROOT, SRC["db3"], "analysis", "ci_metrics.csv"))
for m in MAIN_DISCRIM:
    if (m, "AUC") not in db1 or (m, "AUC") not in db3:
        sys.exit("method %s missing from the DB1/DB3 analyses" % m)

# ---------------------------------------------------------------- Table 1: discrimination, DB1 and DB3
rows = []
for m in MAIN_DISCRIM:
    rows.append("%s & %s & %s & %s & %s \\\\" % (MAIN_NAME[m], cell(db1, m, "AUC", f3), cell(db1, m, "EER", p1),
                                                 cell(db3, m, "AUC", f3), cell(db3, m, "EER", p1)))
write("tab_main_discrimination.tex", colsep="3pt", text=r"""\begin{tabular}{lcccc}
\toprule
 & \multicolumn{2}{c}{FVC2000 DB1} & \multicolumn{2}{c}{FVC2002 DB3} \\
\cmidrule(lr){2-3}\cmidrule(lr){4-5}
Method & AUC & EER (\%) & AUC & EER (\%) \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
""")

# ---------------------------------------------------------------- Table 2: operating points, DB1
rows = []
for m in MAIN_OPER:
    rows.append("%s & %s & %s & %s & %s \\\\" % (NAME[m], cell(db1, m, "TAR@0.001", p1), cell(db1, m, "realisedFAR@0.001", p2),
                                                 cell(db1, m, "TAR@0.01", p1), cell(db1, m, "realisedFAR@0.01", p2)))
write("tab_main_operating.tex", r"""\begin{tabular}{lcccc}
\toprule
 & \multicolumn{2}{c}{nominal FAR $=10^{-3}$} & \multicolumn{2}{c}{nominal FAR $=10^{-2}$} \\
\cmidrule(lr){2-3}\cmidrule(lr){4-5}
Method & TAR (\%) & realized FAR (\%) & TAR (\%) & realized FAR (\%) \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
""")

# ---------------------------------------------------------------- supplement: all methods, both datasets (two panels each)
for tag, d in (("db1", db1), ("db3", db3)):
    present = [m for m in SUPP_ALL if (m, "AUC") in d]
    rows = ["%s & %s & %s \\\\" % (NAME[m], cell(d, m, "AUC", f3), cell(d, m, "EER", p1)) for m in present]
    write("tab_supp_discrim_%s.tex" % tag, r"""\begin{tabular}{lcc}
\toprule
Method & AUC & EER (\%) \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
""")
    rows = ["%s & %s & %s & %s & %s \\\\" % (NAME[m], cell(d, m, "TAR@0.001", p1), cell(d, m, "realisedFAR@0.001", p2),
                                            cell(d, m, "TAR@0.01", p1), cell(d, m, "realisedFAR@0.01", p2)) for m in present]
    write("tab_supp_oper_%s.tex" % tag, r"""\begin{tabular}{lcccc}
\toprule
 & \multicolumn{2}{c}{nominal FAR $=10^{-3}$} & \multicolumn{2}{c}{nominal FAR $=10^{-2}$} \\
\cmidrule(lr){2-3}\cmidrule(lr){4-5}
Method & TAR (\%) & realized FAR (\%) & TAR (\%) & realized FAR (\%) \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
""")
    # all four operating points, TAR only (realized FAR in the CSV)
    rows = []
    for m in SUPP_ALL:
        if (m, "AUC") not in d:
            continue
        rows.append("%s & " % NAME[m] + " & ".join(cell(d, m, "TAR@" + far, p1) for far in FARS) + r" \\")
    write("tab_supp_tar_%s.tex" % tag, r"""\begin{tabular}{lcccc}
\toprule
Method & """ + " & ".join("TAR@$%s$" % FARLAB[far] for far in FARS) + r""" \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
""")

# ---------------------------------------------------------------- supplement: runtime
rt = list(csv.DictReader(open(os.path.join(SRC_ROOT, SRC["runtime"], "summary_runtime_mean_sd.csv"))))
order = ["GlobalBetti_L1", "GlobalPI_L2", "DT_L2", "DT_L2_150", "DT_L2_150_aligned", "V23_full", "OT", "Fusion_eq", "Fusion"]
rtd = {r["method"]: r for r in rt}
rows = []
for m in order:
    r = rtd[m]
    note = ""
    if r.get("shared_precompute_minutes") not in (None, "", "NA"):
        note = r" (alignment search %.0f min, shared)" % float(r["shared_precompute_minutes"])
    rows.append("%s & %.1f $\\pm$ %.1f%s \\\\" % (NAME[m], float(r["seconds_mean"]), float(r["seconds_sd"]), note))
write("tab_supp_runtime.tex", r"""\begin{tabular}{lc}
\toprule
Method & Seconds per fold \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
""")

# ---------------------------------------------------------------- supplement: tuned fusion weights
fw = list(csv.DictReader(open(os.path.join(SRC_ROOT, SRC["run_db1"], "fusion_weights_per_fold.csv"))))
ws = [float(r["best_w"]) for r in fw]
# macro file, \input in both preambles: it must set no length (a \setlength{\tabcolsep} here would silently change
# the document default); the unused mean/SD macros are gone
with open(os.path.join(OUT, "fusion_weights.tex"), "w") as f:
    f.write("%% generated by Code/make_tables.py -- do not edit by hand\n%% per-fold tuned fusion weights (DB1)\n\\newcommand{\\fusionweights}{%s}\n" % ", ".join("%.2f" % w for w in ws))
print("wrote " + os.path.relpath(os.path.join(OUT, "fusion_weights.tex"), ROOT))

# ---------------------------------------------------------------- Table 3 + supplement: sensitivity (ten conditions)
def read_sens(path):
    d = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            d[(r["condition"], r["method"], r["metric"])] = r
    return d


B = read_sens(os.path.join(SRC_ROOT, SRC["sens_B"], "summary", "ci_perturbed_minus_clean.csv"))
R = read_sens(os.path.join(SRC_ROOT, SRC["sens_R"], "summary", "ci_perturbed_minus_clean.csv"))
# short labels for the tables; the conditions table (supplement) spells them out
COND = [("filter_m2", "filter threshold $-2$", B), ("filter_m4", "filter threshold $-4$", B),
        ("scale_055", "Gabor scale 0.55", B), ("scale_075", "Gabor scale 0.75", B),
        ("segment_020", "segmentation 0.20", B), ("segment_040", "segmentation 0.40", B),
        ("dry_typ", "contrast loss", R), ("wet_typ", "ridge spread", R),
        ("noise_typ", "noise", R), ("blur_typ", "blur", R)]
# the strengths of the calibrated degradations come from calibration.json (never typed here), and each
# condition's degrade manifest must record the same strength (consistency guard for a regenerated calibration)
import json
CAL = os.path.join(SRC_ROOT, SRC["calibration"])
cal = json.load(open(os.path.join(CAL, "calibration.json")))


def strength_of(cond):
    """The condition's strength from calibration.json; the condition's degrade manifest is REQUIRED and must
    record the same family and strength (a table is never generated for a condition whose generated images
    cannot be tied to the calibration)."""
    ch = cal["chosen"][cond]
    dm_path = os.path.join(CAL, cond, "degrade_manifest.json")
    if not os.path.exists(dm_path):
        sys.exit("%s: degrade manifest missing (%s); cannot tie the condition's images to calibration.json" % (cond, dm_path))
    dm = json.load(open(dm_path))
    if abs(float(dm["strength"]) - float(ch["strength"])) > 1e-12 or dm["family"] != ch["family"]:
        sys.exit("%s: degrade manifest strength %s / family %s differ from calibration.json (%s / %s)" % (cond, dm["strength"], dm["family"], ch["strength"], ch["family"]))
    return float(ch["strength"])


def fmt_strength(x, decimals_at_least=0):
    t = "%g" % x
    if decimals_at_least and "." not in t:
        t += "." + "0" * decimals_at_least
    return t


COND_LONG = {"filter_m2": "enhancement ridge-filter threshold $-2$ (default $-3$)",
             "filter_m4": "enhancement ridge-filter threshold $-4$ (default $-3$)",
             "scale_055": "enhancement Gabor scale factor 0.55 (default 0.65)",
             "scale_075": "enhancement Gabor scale factor 0.75 (default 0.65)",
             "segment_020": "enhancement segmentation threshold 0.20 (default 0.10)",
             "segment_040": "enhancement segmentation threshold 0.40 (default 0.10)",
             "dry_typ": "patchy contrast loss of the raw probe image, strength %s" % fmt_strength(strength_of("dry_typ")),
             "wet_typ": "ridge spread with blur of the raw probe image, strength %s" % fmt_strength(strength_of("wet_typ")),
             "noise_typ": "additive Gaussian noise on the raw probe image, $\\sigma = %s$ gray levels" % fmt_strength(strength_of("noise_typ")),
             "blur_typ": "Gaussian blur of the raw probe image, $\\sigma = %s$ px" % fmt_strength(strength_of("blur_typ"), 1)}


def sens_rows(method):
    """Condition | AUC | dAUC [CI] | TAR@1e-3 (rFAR) | TAR@1e-2 (rFAR); agreement and EER live in the conditions table / CSVs."""
    rows = []
    a0 = R[("dry_typ", method, "AUC")]
    t30 = R[("dry_typ", method, "TAR@0.001")]; f30 = R[("dry_typ", method, "realisedFAR@0.001")]
    t20 = R[("dry_typ", method, "TAR@0.01")]; f20 = R[("dry_typ", method, "realisedFAR@0.01")]
    rows.append("clean (reference) & %s & --- & %s (%s) & %s (%s) \\\\" % (
        f3(a0["point_clean"]), p1(t30["point_clean"]), p2(f30["point_clean"]), p1(t20["point_clean"]), p2(f20["point_clean"])))
    rows.append(r"\midrule")
    for i, (c, lab, D) in enumerate(COND):
        if i == 6:
            rows.append(r"\midrule")
        a = D[(c, method, "AUC")]
        t3 = D[(c, method, "TAR@0.001")]; f3_ = D[(c, method, "realisedFAR@0.001")]
        t2 = D[(c, method, "TAR@0.01")]; f2 = D[(c, method, "realisedFAR@0.01")]
        rows.append("%s & %s & %s & %s (%s) & %s (%s) \\\\" % (
            lab, f3(a["point_perturbed"]), diff_cell(a),
            p1(t3["point_perturbed"]), p2(f3_["point_perturbed"]), p1(t2["point_perturbed"]), p2(f2["point_perturbed"])))
    return rows


HEAD = r"""\begin{tabular}{lcccc}
\toprule
Condition & AUC & $\Delta$AUC [95\% CI] & TAR@$10^{-3}$ (rFAR) & TAR@$10^{-2}$ (rFAR) \\
\midrule
"""
write("tab_main_sensitivity.tex", HEAD + "\n".join(sens_rows("Fusion_eq")) + "\n\\bottomrule\n\\end{tabular}\n")
for m in ["V23_full", "OT", "GlobalBetti_L1", "GlobalPI_L2", "DT_L2_150", "Fusion"]:   # aligned baseline not in the sensitivity runs
    write("tab_supp_sensitivity_%s.tex" % m, HEAD + "\n".join(sens_rows(m)) + "\n\\bottomrule\n\\end{tabular}\n")

# conditions table: definition, pixel agreement with the canonical binaries (median / 10th percentile / min), Fusion_eq EER
rows = []
for i, (c, lab, D) in enumerate(COND):
    if i == 6:
        rows.append(r"\midrule")
    a = D[(c, "Fusion_eq", "AUC")]; e = D[(c, "Fusion_eq", "EER")]
    rows.append("%s & %s & %.3f & %.3f & %s \\\\" % (lab, COND_LONG[c], float(a["agreement_median"]),
                float(a["agreement_min"]), p1(e["point_perturbed"])))
write("tab_supp_sensitivity_conditions.tex", r"""\begin{tabular}{lp{118pt}ccc}
\toprule
 & & \multicolumn{2}{c}{pixel agreement} & \\
\cmidrule(lr){3-4}
Condition & Definition & median & min & Fusion EER (\%) \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
""")

# ---------------------------------------------------------------- supplement: ablations (DB1, intervals)
ABL = ["V23_full", "V23_globalOnly", "V23_localOnly", "V23_noTopK", "OT", "OT_noGeo", "OT_noDustbin", "OT_ridgeOnly", "OT_valleyOnly"]
rows = [r"\multicolumn{4}{l}{\emph{TopoLR}} \\"]
for m in ABL:
    if m == "OT":
        rows.append(r"\midrule")
        rows.append(r"\multicolumn{4}{l}{\emph{LPOT}} \\")
    rows.append("\\quad %s & %s & %s & %s \\\\" % (NAME[m], cell(db1, m, "AUC", f3), cell(db1, m, "EER", p1), cell(db1, m, "TAR@0.001", p1)))
write("tab_supp_ablations.tex", r"""\begin{tabular}{lccc}
\toprule
Variant & AUC & EER (\%) & TAR@$10^{-3}$ (\%) \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
""")

# ================================================================ supplement: construction comparison (cubical - triangulated)
SHORT = dict(NAME, **{"Fusion_eq": "Fusion (equal)", "Fusion": "Fusion (tuned)"})
PRINCIPAL = ["GlobalBetti_L1", "GlobalPI_L2", "V23_full", "OT", "Fusion_eq"]
blocks = []
for tag, d in (("FVC2000 DB1", SRC["cubical_db1"]), ("FVC2002 DB3", SRC["cubical_db3"])):
    rows_in = list(csv.DictReader(open(os.path.join(SRC_ROOT, d, "analysis", "ci_cubical_minus_triangulated.csv"))))
    look = {(r["method"], r["metric"]): r for r in rows_in}
    rows = ["\\multicolumn{5}{l}{\\emph{%s}} \\\\" % tag]
    for m in PRINCIPAL:
        a = look[(m, "AUC")]; e = look[(m, "EER")]; t3 = look[(m, "TAR@0.001")]; t2 = look[(m, "TAR@0.01")]
        rows.append("%s & %s & %s & %s & %s \\\\" % (SHORT[m], diff_cell(a), diff_cell(e, 100, 1), diff_cell(t3, 100, 1), diff_cell(t2, 100, 1)))
    blocks.append("\n".join(rows))
write("tab_supp_cubical.tex", colsep="3pt", text=r"""\begin{tabular}{lcccc}
\toprule
Method & $\Delta$AUC & $\Delta$EER (pp) & $\Delta$TAR@$10^{-3}$ (pp) & $\Delta$TAR@$10^{-2}$ (pp) \\
\midrule
""" + "\n\\addlinespace\n".join(blocks) + r"""
\bottomrule
\end{tabular}
""")

# ================================================================ supplement: transfer of the DB1 models to DB3
tr = list(csv.DictReader(open(os.path.join(SRC_ROOT, SRC["transfer"], "analysis", "ci_transfer_minus_replication.csv"))))
look = {(r["calibration"], r["method"], r["metric"]): r for r in tr}
rows_a, rows_b = [], []
for m in ["V23_full", "Fusion_eq", "Fusion"]:
    a = look[("frozen", m, "AUC")]; e = look[("frozen", m, "EER")]
    f3_ = look[("frozen", m, "realisedFAR@0.001")]; t3f = look[("frozen", m, "TAR@0.001")]
    t3r = look[("recal", m, "TAR@0.001")]; t2r = look[("recal", m, "TAR@0.01")]
    f3r = look[("recal", m, "realisedFAR@0.001")]; f2r = look[("recal", m, "realisedFAR@0.01")]
    rows_a.append("%s & %.3f $\\to$ %.3f & %s & %.1f $\\to$ %.1f & %s \\\\" % (
        SHORT[m], float(a["point_replication"]), float(a["point_transfer"]), diff_cell(a), 100 * float(e["point_replication"]), 100 * float(e["point_transfer"]), diff_cell(e, 100, 1)))
    rows_b.append("%s & %.1f (%.2f) & %.1f (%.2f) & %s & %.1f (%.2f) & %s \\\\" % (
        SHORT[m], 100 * float(t3f["point_transfer"]), 100 * float(f3_["point_transfer"]),
        100 * float(t3r["point_transfer"]), 100 * float(f3r["point_transfer"]), diff_cell(t3r, 100, 1),
        100 * float(t2r["point_transfer"]), 100 * float(f2r["point_transfer"]), diff_cell(t2r, 100, 1)))
write("tab_supp_transfer_a.tex", r"""\begin{tabular}{lcccc}
\toprule
 & \multicolumn{2}{c}{AUC} & \multicolumn{2}{c}{EER (\%)} \\
\cmidrule(lr){2-3}\cmidrule(lr){4-5}
Method & repl.\ $\to$ transfer & $\Delta$ [95\% CI] & repl.\ $\to$ transfer & $\Delta$ (pp) \\
\midrule
""" + "\n".join(rows_a) + r"""
\bottomrule
\end{tabular}
""")
write("tab_supp_transfer_b.tex", colsep="3pt", text=r"""\begin{tabular}{lccccc}
\toprule
 & frozen & \multicolumn{4}{c}{target-calibrated} \\
\cmidrule(lr){2-2}\cmidrule(lr){3-6}
 & $10^{-3}$ & \multicolumn{2}{c}{$10^{-3}$} & \multicolumn{2}{c}{$10^{-2}$} \\
\cmidrule(lr){2-2}\cmidrule(lr){3-4}\cmidrule(lr){5-6}
Method & TAR (rFAR) & TAR (rFAR) & $\Delta$ & TAR (rFAR) & $\Delta$ \\
\midrule
""" + "\n".join(rows_b) + r"""
\bottomrule
\end{tabular}
""")

# ================================================================ supplement: tuned versus equal-weight fusion
cols = {}
for tag, d in (("DB1", SRC["run_db1"]), ("DB3", SRC["run_db3"])):
    pdiff = {(r["method_a"], r["method_b"], r["metric"]): r for r in csv.DictReader(open(os.path.join(SRC_ROOT, d, "analysis", "ci_paired_differences.csv")))}
    wts = [float(r["best_w"]) for r in csv.DictReader(open(os.path.join(SRC_ROOT, d, "fusion_weights_per_fold.csv")))]
    g = lambda met, s, nd: diff_cell(pdiff[("Fusion", "Fusion_eq", met)], s, nd)
    cols[tag] = [", ".join("%.2f" % w for w in wts), g("AUC", 1, 3), g("EER", 100, 1), g("TAR@0.001", 100, 1), g("TAR@0.01", 100, 1)]
LAB = ["tuned weights $w$ (folds 1--5)", "$\\Delta$AUC", "$\\Delta$EER (pp)", "$\\Delta$TAR@$10^{-3}$ (pp)", "$\\Delta$TAR@$10^{-2}$ (pp)"]
rows = ["%s & %s & %s \\\\" % (LAB[i], cols["DB1"][i], cols["DB3"][i]) for i in range(5)]
write("tab_supp_fusion_tuning.tex", r"""\begin{tabular}{lcc}
\toprule
 & FVC2000 DB1 & FVC2002 DB3 \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
""")

# ================================================================ supplement: the matcher on our binary maps
cols = {}
for tag, d, cm in (("DB1", SRC["db1"], db1), ("DB3", SRC["db3"], db3)):
    pdiff = {(r["method_a"], r["method_b"], r["metric"]): r for r in csv.DictReader(open(os.path.join(SRC_ROOT, d, "analysis", "ci_paired_differences.csv")))}
    g = lambda met, s, nd: diff_cell(pdiff[("NBIS_bin", "NBIS_gray", met)], s, nd)
    cols[tag] = [cell(cm, "NBIS_bin", "AUC", f3), cell(cm, "NBIS_bin", "EER", p1), g("AUC", 1, 3), g("EER", 100, 1), g("TAR@0.001", 100, 1), g("TAR@0.01", 100, 1)]
LAB = ["AUC, binary input", "EER (\\%), binary input", "$\\Delta$AUC", "$\\Delta$EER (pp)", "$\\Delta$TAR@$10^{-3}$ (pp)", "$\\Delta$TAR@$10^{-2}$ (pp)"]
rows = ["%s & %s & %s \\\\" % (LAB[i], cols["DB1"][i], cols["DB3"][i]) for i in range(6)]
rows.insert(2, "\\addlinespace\n\\multicolumn{3}{l}{\\emph{binary $-$ grayscale input (paired)}} \\\\")
write("tab_supp_nbis_binary.tex", r"""\begin{tabular}{lcc}
\toprule
 & FVC2000 DB1 & FVC2002 DB3 \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
""")

# ================================================================ supplement: calibration of the degradation strengths
dist = json.load(open(os.path.join(CAL, "calibration_distribution.json")))
nfiq = json.load(open(os.path.join(CAL, "calibration_nfiq.json")))
ver = json.load(open(os.path.join(CAL, "calibration_verify.json")))
DESC = [("r", "ridge fraction"), ("fg", "print area (fraction)"), ("c0", "ridge pieces (count)"), ("c1", "enclosed valleys (count)"), ("iso", "specks (count)"), ("w", "ridge width (px)")]
SYM = {"r": "$r$", "fg": "$fg$", "c0": "$c_0$", "c1": "$c_1$", "iso": "$iso$", "w": "$w$"}
b = cal["benchmark"]
rows = []
for k, lab in DESC:
    rows.append("%s & %s & %.3g & [%.3g, %.3g] & %.3g & %.3g \\\\" % (lab, SYM[k], b["clean_median"][k], b["clean_band_5_95"][k][0], b["clean_band_5_95"][k][1], b["T50_within_finger"][k], b["T90_within_finger"][k]))
write("tab_supp_calibration_descriptors.tex", r"""\begin{tabular}{llcccc}
\toprule
 & & \multicolumn{2}{c}{clean images} & \multicolumn{2}{c}{within-finger $|\Delta|$} \\
\cmidrule(lr){3-4}\cmidrule(lr){5-6}
Descriptor & symbol & median & 5--95\% band & median & 90th pct \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
""")
FAM = [("dry_typ", "contrast loss", "dry"), ("wet_typ", "ridge spread", "wet"), ("noise_typ", "noise", "noise"), ("blur_typ", "blur", "blur")]
rows = []
for c, lab, fam in FAM:
    ch = cal["chosen"][c]; dc = dist["conditions"][c]
    if ch["level"] != "typical" or ch["family"] != fam:
        sys.exit("%s: calibration.json reports level %s / family %s" % (c, ch["level"], ch["family"]))
    # binding descriptor(s): every statistic whose median-rule ratio equals 1.00 at the chosen strength (ties are
    # reported together; calibration.json's binding_statistic field names only the first in dictionary order)
    ratios = cal["scan"]["%s@%g" % (fam, ch["strength"])]["ratio_T50"]
    binding = [k for k, _ in DESC if abs(ratios[k] - 1.0) < 1e-9]
    if ch["binding_statistic"] not in binding:
        sys.exit("%s: binding_statistic %s is not at ratio 1.00 (%s)" % (c, ch["binding_statistic"], ratios))
    worst_k = max(DESC, key=lambda kv: dc[kv[0]]["frac_over_T90"])
    rows.append("%s & %s & %s & %.2f & %.0f\\%% (%s) & %.2f \\\\" % (
        lab, fmt_strength(ch["strength"]), ", ".join(SYM[k] for k in binding), ver[c]["R50"], 100 * dc[worst_k[0]]["frac_over_T90"], SYM[worst_k[0]], nfiq[c]["mean"]))
write("tab_supp_calibration_chosen.tex", r"""\begin{tabular}{lccccc}
\toprule
Family & strength & binding & $R_{50}$ & beyond 90th pct & NFIQ \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
""")
# scan table: median-rule ratio R50 for every grid strength (admissible flag)
rows = []
for c, lab, fam in FAM:
    cells = []
    for g in cal["grid"][fam]:
        key = "%s@%g" % (fam, g)
        e = cal["scan"].get(key)
        if e is None:
            cells.append("%g: --" % g); continue
        adm = e["n_errors"] == 0 and e["n_warnings"] == 0 and e["in_band"]
        cells.append("%g: %.2f%s" % (g, e["R50"], "" if adm else "$^{\\dagger}$"))
    rows.append("%s & %s \\\\" % (lab, "; ".join(cells)))
write("tab_supp_calibration_scan.tex", r"""\begin{tabular}{lp{250pt}}
\toprule
Family & strength: $R_{50}$ on the calibration subset ($^{\dagger}$ inadmissible) \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
""")
print("done")
