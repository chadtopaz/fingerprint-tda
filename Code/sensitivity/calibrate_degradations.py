#!/usr/bin/env python3
"""
Calibration of the acquisition degradations (sensitivity experiment, R4-10). Chooses the strength of each
degradation family (Code/sensitivity/degrade_raw.py) by a rule that uses only the inputs (raw images and
canonical binaries) — never a verification result.

RULE. Monitored statistics of an enhanced binary B (ridge = white):
    r    ridge fraction                       fg   print-area proxy: fraction covered by the morphological
    c0   number of 8-connected ridge pieces        closing of the ridge set with a radius-10 disk
    c1   number of enclosed valley regions    iso  number of ridge pieces of <= 2 pixels (specks)
         (4-connected valley components not     w    mean ridge width = ridge area / skeleton length (px)
         touching the image border)
Natural yardstick: for every finger, all 28 pairs of its 8 canonical impressions; |delta stat| pooled over
the 100 fingers (2,800 pairs) -> T50_k (median) and T90_k (90th percentile) per statistic.
For a family at strength s, on the calibration subset (one impression per finger, fixed seed), the
degraded raw image is enhanced with the canonical procedure (library defaults) and
    D_k(s) = median over the subset of |stat_k(degraded) - stat_k(clean)|.
    chosen strength (condition <family>_typ, level "typical") = the largest grid strength with
        max_k D_k(s) / T50_k <= 1
    (the degraded image differs from its own clean version by no more than two impressions of the same
    finger typically do), subject to: no enhancement failure on the subset at s, no subset image whose ridge
    fraction falls outside 0.05-0.5 (the enhancer's sanity flag), and the subset's median ridge fraction inside
    the clean database's [5th, 95th] percentile band. The grid for each family is fixed below; if even the
    smallest grid strength fails the rule, the family is reported as inadmissible (not retuned).
The ratio to the 90th percentile, D_k(s) / T90_k, is recorded in every scan record as a diagnostic and
selects nothing. The realism ratios of the full 800-image sets are recomputed after generation and reported;
no retuning. NFIQ (NBIS nfiq, 1 best .. 5 worst) of the raw images is reported as an external quality descriptor.

Stages (project root; the pinned preprocessing environment):
  benchmark            -> <work>/calibration_benchmark.json (+ calibration_stats_clean.csv)
  scan --family F --strengths a,b,c   (appends to <work>/calibration_scan.json; chunk to stay under time caps)
  select               -> <work>/calibration.json  (chosen strengths; consumed by degrade_raw.py)
  verify --conditions  -> realism ratios of the generated full sets (needs the enhanced binaries)
  distribution --conditions -> <work>/calibration_distribution.json: per-image |delta stat| of the full sets
                          against the clean binaries (median, 90th percentile, fraction beyond T50 / T90)
  nfiq --dirs ...      -> NFIQ distributions of raw directories -> <work>/calibration_nfiq.json
  subset --file F --keep a,b,... -> keep only the listed top-level entries of a JSON record file (e.g. the NFIQ file
                          restricted to the clean directory and the calibrated conditions); writes F in place
  montage              -> <work>/calibration_montage.png (one finger: raw and binary, clean vs each condition)
"""
import argparse, hashlib, json, os, re, subprocess, sys, tempfile, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from degrade_raw import degrade, seed_for, FAMILIES   # noqa: E402

STATS = ("r", "fg", "c0", "c1", "iso", "w")
GRID = {"dry": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.05, 1.1, 1.15, 1.2, 1.4, 1.6, 1.8, 2.0],
        "wet": [0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0],
        "noise": [4, 6, 8, 10, 12, 16, 20, 25, 30, 40, 50],
        "blur": [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0]}
SUBSET_SEED = 20261014
EIGHT = np.ones((3, 3), dtype=bool)
FOUR = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=bool)
DISK10 = (lambda r: (np.add.outer(np.arange(-r, r + 1) ** 2, np.arange(-r, r + 1) ** 2) <= r * r))(10)


def binary_stats(B):
    """B: boolean ridge mask (True = ridge)."""
    from skimage.morphology import skeletonize
    lab, n0 = ndi.label(B, structure=EIGHT)
    sizes = np.bincount(lab.ravel())[1:] if n0 else np.array([], dtype=int)
    vlab, nv = ndi.label(~B, structure=FOUR)
    border = np.unique(np.concatenate([vlab[0, :], vlab[-1, :], vlab[:, 0], vlab[:, -1]]))
    c1 = int(nv - len(set(border.tolist()) - {0}))
    skel = skeletonize(B)
    w = float(B.sum() / max(int(skel.sum()), 1))
    fg = float(ndi.binary_closing(B, structure=DISK10).mean())
    return dict(r=float(B.mean()), fg=fg, c0=int(n0), c1=c1, iso=int((sizes <= 2).sum()), w=w)


def read_bin(p):
    a = np.array(Image.open(p))
    return a == 255


def enhance_arr(img_u8):
    import fingerprint_enhancer, warnings
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter("always")
        out = fingerprint_enhancer.enhance_fingerprint(img_u8)
    out = np.asarray(out)
    if out.dtype != bool or out.shape != img_u8.shape: raise ValueError("unexpected enhancer output")
    frac = float(out.mean())
    warn = "ridge fraction > 0.5" if frac > 0.5 else ("ridge fraction < 0.05" if frac < 0.05 else None)   # preprocess.py's warnings
    return out, warn


def _scan_one(job):
    """Degrade + enhance + stats for one image. job = (raw path, family, strength, cond tag)."""
    p, fam, s, tag = job
    a = np.array(Image.open(p))
    deg = degrade(a, fam, s, np.random.default_rng(seed_for(os.path.basename(p), tag)))
    try:
        B, warn = enhance_arr(deg)
    except Exception as e:                                  # noqa: BLE001
        return dict(file=os.path.basename(p), error=f"{type(e).__name__}: {e}")
    d = binary_stats(B); d["file"] = os.path.basename(p); d["warning"] = warn
    d["mean_abs_change_gray"] = float(np.abs(deg.astype(float) - a.astype(float)).mean())
    return d


def _stats_one(p):
    d = binary_stats(read_bin(p)); d["file"] = os.path.basename(p); return d


def finger_imp(name):
    m = re.match(r"(\d+)_(\d+)", name); return int(m.group(1)), int(m.group(2))


def load_json(p, default):
    return json.load(open(p)) if os.path.exists(p) else default


def stage_benchmark(A):
    bins = sorted(p for p in Path(A.canonical).iterdir() if re.fullmatch(r"\d+_\d+_bin\.png", p.name))
    if len(bins) != 800: sys.exit("expected 800 canonical binaries")
    with ProcessPoolExecutor(A.workers) as ex: rows = list(ex.map(_stats_one, map(str, bins), chunksize=8))
    import csv
    with open(Path(A.work) / "calibration_stats_clean.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["file"] + list(STATS)); w.writeheader()
        for r in rows: w.writerow({k: r[k] for k in ["file"] + list(STATS)})
    by_f = {}
    for r in rows: by_f.setdefault(finger_imp(r["file"])[0], []).append(r)
    diffs = {k: [] for k in STATS}
    for f, lst in by_f.items():
        for i in range(len(lst)):
            for j in range(i + 1, len(lst)):
                for k in STATS: diffs[k].append(abs(lst[i][k] - lst[j][k]))
    T50 = {k: float(np.median(diffs[k])) for k in STATS}; T90 = {k: float(np.percentile(diffs[k], 90)) for k in STATS}
    band = {k: [float(np.percentile([r[k] for r in rows], 5)), float(np.percentile([r[k] for r in rows], 95))] for k in STATS}
    med = {k: float(np.median([r[k] for r in rows])) for k in STATS}
    rng = np.random.default_rng(SUBSET_SEED)
    subset = [f"{f}_{int(rng.integers(1, 9))}.tif" for f in sorted(by_f)]
    bench = dict(created=time.strftime("%Y-%m-%dT%H:%M:%S%z"), canonical_dir=str(Path(A.canonical).resolve()), n_images=len(rows), n_within_finger_pairs=sum(len(v) for v in diffs.values()) // len(STATS),
                 statistics=list(STATS), clean_median=med, clean_band_5_95=band, T50_within_finger=T50, T90_within_finger=T90,
                 subset_seed=SUBSET_SEED, subset=subset, rule=__doc__.split("RULE.")[1].split("Stages")[0].strip())
    json.dump(bench, open(Path(A.work) / "calibration_benchmark.json", "w"), indent=1)
    print("clean medians:", {k: round(v, 4) for k, v in med.items()})
    print("T50 within-finger |delta|:", {k: round(v, 4) for k, v in T50.items()})
    print("T90 within-finger |delta|:", {k: round(v, 4) for k, v in T90.items()})
    print("clean 5-95 band:", {k: [round(x, 4) for x in v] for k, v in band.items()})
    print(f"subset: {len(subset)} images (one impression per finger, seed {SUBSET_SEED})")


def stage_scan(A):
    bench = json.load(open(Path(A.work) / "calibration_benchmark.json"))
    clean = {r["file"]: r for r in _read_clean(A)}
    scan_path = Path(A.work) / "calibration_scan.json"; scan = load_json(scan_path, {})
    fam = A.family
    if fam not in FAMILIES: sys.exit(f"family must be one of {FAMILIES}")
    strengths = [float(x) for x in A.strengths.split(",")]
    for s in strengths:
        tag = f"{fam}@{s:g}"
        jobs = [(str(Path(A.raw) / name), fam, s, tag) for name in bench["subset"]]
        t0 = time.time()
        with ProcessPoolExecutor(A.workers) as ex: res = list(ex.map(_scan_one, jobs, chunksize=4))
        errs = [r for r in res if "error" in r]; ok = [r for r in res if "error" not in r]
        warns = sum(1 for r in ok if r["warning"])
        D = {}; ratio50 = {}; ratio90 = {}
        for k in STATS:
            d = [abs(r[k] - clean[r["file"].replace(".tif", "_bin.png")][k]) for r in ok]
            D[k] = float(np.median(d)) if d else float("nan")
            ratio50[k] = D[k] / bench["T50_within_finger"][k] if bench["T50_within_finger"][k] > 0 else float("inf")
            ratio90[k] = D[k] / bench["T90_within_finger"][k] if bench["T90_within_finger"][k] > 0 else float("inf")
        r_med = float(np.median([r["r"] for r in ok])) if ok else float("nan")
        rec = dict(family=fam, strength=s, n=len(res), n_errors=len(errs), n_warnings=warns, D_median_abs_change=D,
                   ratio_T50=ratio50, ratio_T90=ratio90, R50=max(ratio50.values()), R90=max(ratio90.values()),
                   ridge_fraction_median=r_med, in_band=bool(bench["clean_band_5_95"]["r"][0] <= r_med <= bench["clean_band_5_95"]["r"][1]),
                   mean_abs_change_gray_median=float(np.median([r["mean_abs_change_gray"] for r in ok])) if ok else float("nan"),
                   seconds=round(time.time() - t0, 1))
        scan[tag] = rec
        json.dump(scan, open(scan_path, "w"), indent=1)
        print(f"{fam:>5} s={s:<5g} R50={rec['R50']:.2f} R90={rec['R90']:.2f} err={len(errs)} warn={warns} r_med={r_med:.3f} in_band={rec['in_band']} "
              f"|dgray|={rec['mean_abs_change_gray_median']:.1f}  binding: " + ", ".join(f"{k}:{ratio50[k]:.2f}" for k in STATS) + f"  ({rec['seconds']} s)")


def _read_clean(A):
    import csv
    rows = list(csv.DictReader(open(Path(A.work) / "calibration_stats_clean.csv")))
    for r in rows:
        for k in STATS: r[k] = float(r[k])
    return rows


def stage_select(A):
    bench = json.load(open(Path(A.work) / "calibration_benchmark.json")); scan = load_json(Path(A.work) / "calibration_scan.json", {})
    chosen = {}; table = []
    for fam in FAMILIES:
        recs = sorted([v for v in scan.values() if v["family"] == fam], key=lambda v: v["strength"])
        grid_done = [v["strength"] for v in recs]
        if sorted(grid_done) != sorted(GRID[fam]): sys.exit(f"{fam}: scanned strengths {grid_done} != prescribed grid {GRID[fam]}")
        # admissible: no enhancement failure (n_errors) and no image flagged by the enhancer's ridge-fraction sanity
        # check (n_warnings counts ridge fractions outside 0.05-0.5), subset median ridge fraction in the clean band
        adm = [v for v in recs if v["n_errors"] == 0 and v["n_warnings"] == 0 and v["in_band"]]
        cand = [v for v in adm if v["R50"] <= 1.0]
        if cand:
            v = max(cand, key=lambda v: v["strength"])
            chosen[f"{fam}_typ"] = dict(family=fam, strength=v["strength"], level="typical",
                                        R50=v["R50"], R90=v["R90"], binding_statistic=max(v["ratio_T50"], key=v["ratio_T50"].get))
            table.append((f"{fam}_typ", v["strength"], v["R50"], v["R90"], chosen[f"{fam}_typ"]["binding_statistic"]))
        else:
            table.append((f"{fam}_typ", None, None, None, "INADMISSIBLE"))
    # the rule recorded in calibration.json is THIS script's rule (the benchmark file keeps the text it was written with)
    bench_rec = dict(bench); bench_rec["rule_as_recorded_at_benchmark_time"] = bench_rec.pop("rule", None)
    cal = dict(created=time.strftime("%Y-%m-%dT%H:%M:%S%z"), rule=__doc__.split("RULE.")[1].split("Stages")[0].strip(), grid=GRID, benchmark=bench_rec, scan=scan, chosen=chosen,
               scan_sha256=hashlib.sha256(json.dumps(scan, sort_keys=True).encode()).hexdigest())
    json.dump(cal, open(Path(A.work) / "calibration.json", "w"), indent=1)
    print(f"{'condition':>10} {'strength':>9} {'R50':>6} {'R90':>6}  binding")
    for row in table: print(f"{row[0]:>10} {str(row[1]):>9} {('%.2f' % row[2]) if row[2] is not None else '   -':>6} {('%.2f' % row[3]) if row[3] is not None else '   -':>6}  {row[4]}")
    print("written:", Path(A.work) / "calibration.json")


def stage_verify(A):
    """Realism ratios of the generated FULL sets (800 images per condition) against the benchmark; also the
    per-image pixel agreement of the degraded-input binaries with the canonical binaries."""
    bench = json.load(open(Path(A.work) / "calibration_benchmark.json")); clean = {r["file"]: r for r in _read_clean(A)}
    out = load_json(Path(A.work) / "calibration_verify.json", {})
    for cond in A.conditions.split(","):
        bd = Path(A.binaries) / cond
        bins = sorted(p for p in bd.iterdir() if re.fullmatch(r"\d+_\d+_bin\.png", p.name))
        if len(bins) != 800: sys.exit(f"{cond}: expected 800 binaries, found {len(bins)}")
        with ProcessPoolExecutor(A.workers) as ex: rows = list(ex.map(_stats_one, map(str, bins), chunksize=8))
        D = {k: float(np.median([abs(r[k] - clean[r["file"]][k]) for r in rows])) for k in STATS}
        R50 = max(D[k] / bench["T50_within_finger"][k] for k in STATS); R90 = max(D[k] / bench["T90_within_finger"][k] for k in STATS)
        agree = [float((read_bin(p) == read_bin(Path(A.canonical) / p.name)).mean()) for p in bins]
        iso_med = float(np.median([r["iso"] for r in rows])); r_med = float(np.median([r["r"] for r in rows]))
        out[cond] = dict(n=len(rows), D_median_abs_change=D, R50=R50, R90=R90, ridge_fraction_median=r_med, iso_median=iso_med,
                         clean_iso_median=float(np.median([clean[r["file"]]["iso"] for r in rows])),
                         agreement=dict(median=float(np.median(agree)), q10=float(np.percentile(agree, 10)), min=float(min(agree))),
                         n_over_half_white=int(sum(r["r"] > 0.5 for r in rows)))
        print(f"{cond:>10}: R50 {R50:.2f}  R90 {R90:.2f}  ridge frac {r_med:.3f}  specks median {iso_med:.0f} (clean {out[cond]['clean_iso_median']:.0f})  "
              f"agreement median {out[cond]['agreement']['median']:.3f} min {out[cond]['agreement']['min']:.3f}  >50% white: {out[cond]['n_over_half_white']}")
        json.dump(out, open(Path(A.work) / "calibration_verify.json", "w"), indent=1)


def stage_distribution(A):
    """Distribution of the per-image |delta| of each statistic (degraded-input binary vs the clean binary of the same
    image) over the FULL 800-image set of each condition: median, 90th percentile, and the fraction of images whose
    |delta| exceeds the natural within-finger T50 / T90 (strictly greater). Added 31 Oct 2026 (audit item B5): this
    stage produces calibration_distribution.json, which the supplement's calibration table reads; the file had been
    written by an ad-hoc snippet with exactly this computation."""
    bench = json.load(open(Path(A.work) / "calibration_benchmark.json")); clean = {r["file"]: r for r in _read_clean(A)}
    T50 = bench["T50_within_finger"]; T90 = bench["T90_within_finger"]
    out_path = Path(A.work) / "calibration_distribution.json"
    out = load_json(out_path, dict(note="per-image |delta| of each statistic vs clean, full 800-image sets; natural within-finger T50/T90 from calibration_benchmark.json",
                                   T50=T50, T90=T90, conditions={}))
    out["T50"] = T50; out["T90"] = T90
    for cond in A.conditions.split(","):
        bd = Path(A.binaries) / cond
        bins = sorted(p for p in bd.iterdir() if re.fullmatch(r"\d+_\d+_bin\.png", p.name))
        if len(bins) != 800: sys.exit(f"{cond}: expected 800 binaries, found {len(bins)}")
        with ProcessPoolExecutor(A.workers) as ex: rows = list(ex.map(_stats_one, map(str, bins), chunksize=8))
        rec = {}
        for k in STATS:
            d = np.array([abs(r[k] - clean[r["file"]][k]) for r in rows], dtype=float)
            rec[k] = dict(median=float(np.median(d)), q90=float(np.percentile(d, 90)),
                          frac_over_T50=float(np.mean(d > T50[k])), frac_over_T90=float(np.mean(d > T90[k])))
        out["conditions"][cond] = rec
        print(f"{cond:>10}: " + "  ".join(f"{k}: med {rec[k]['median']:.4g} q90 {rec[k]['q90']:.4g} >T50 {100 * rec[k]['frac_over_T50']:.1f}% >T90 {100 * rec[k]['frac_over_T90']:.2f}%" for k in STATS))
        json.dump(out, open(out_path, "w"), indent=1)
    print("written:", out_path)


def _nfiq_one(job):
    nfiq_bin, p = job
    a = np.array(Image.open(p))
    with tempfile.NamedTemporaryFile(suffix=".raw", delete=False) as f:
        f.write(np.ascontiguousarray(a).tobytes()); tmp = f.name
    try:
        r = subprocess.run([nfiq_bin, "-raw", "300,300,8,500", tmp], capture_output=True, text=True, check=True)
        return int(r.stdout.strip().split()[-1])
    finally:
        os.unlink(tmp)


def stage_nfiq(A):
    out = load_json(Path(A.work) / "calibration_nfiq.json", {})
    for d in A.dirs.split(","):
        files = sorted(p for p in Path(d).iterdir() if re.fullmatch(r"\d+_\d+\.tif", p.name))
        with ProcessPoolExecutor(A.workers) as ex: q = list(ex.map(_nfiq_one, [(A.nfiq, str(p)) for p in files], chunksize=8))
        dist = {str(k): int(sum(1 for x in q if x == k)) for k in range(1, 6)}
        out[os.path.basename(d.rstrip("/"))] = dict(n=len(q), distribution=dist, mean=float(np.mean(q)), median=float(np.median(q)))
        print(f"{os.path.basename(d.rstrip('/')):>20}: NFIQ 1..5 = {list(dist.values())}  mean {np.mean(q):.2f}")
        json.dump(out, open(Path(A.work) / "calibration_nfiq.json", "w"), indent=1)


def stage_subset(A):
    """Keep only the listed top-level entries of a JSON record file (written in place; used to restrict the NFIQ
    record to the clean directory and the calibrated conditions). Refuses to touch calibration.json, whose hash is
    recorded in every degrade manifest."""
    if os.path.basename(A.file) == "calibration.json": sys.exit("calibration.json is never edited; regenerate it with the select stage")
    rec = json.load(open(A.file)); keep = [k.strip() for k in A.keep.split(",")]
    missing = [k for k in keep if k not in rec]
    if missing: sys.exit(f"entries not in {A.file}: {missing}")
    out = {k: rec[k] for k in rec if k in keep}
    json.dump(out, open(A.file, "w"), indent=1); print(f"{A.file}: kept {list(out)} (dropped {[k for k in rec if k not in keep]})")


def stage_montage(A):
    cal = json.load(open(Path(A.work) / "calibration.json")); conds = list(cal["chosen"])
    name = A.image
    raw = np.array(Image.open(Path(A.raw) / f"{name}.tif")); cb = read_bin(Path(A.canonical) / f"{name}_bin.png")
    panels = [("clean raw", raw), ("clean binary", (cb * 255).astype(np.uint8))]
    for c in conds:
        dr = Path(A.work) / c / f"{name}.tif"; db = Path(A.binaries) / c / f"{name}_bin.png"
        if dr.exists(): panels.append((f"{c} raw", np.array(Image.open(dr))))
        if db.exists(): panels.append((f"{c} binary", (read_bin(db) * 255).astype(np.uint8)))
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    n = len(panels); cols = 6; rows = int(np.ceil(n / cols))
    fig, axs = plt.subplots(rows, cols, figsize=(2.6 * cols, 2.8 * rows))
    for ax, (t, im) in zip(axs.ravel(), panels): ax.imshow(im, cmap="gray", vmin=0, vmax=255); ax.set_title(t, fontsize=8); ax.axis("off")
    for ax in axs.ravel()[n:]: ax.axis("off")
    plt.tight_layout(); plt.savefig(Path(A.work) / f"calibration_montage_{name}.png", dpi=110); print("written:", Path(A.work) / f"calibration_montage_{name}.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["benchmark", "scan", "select", "verify", "distribution", "nfiq", "subset", "montage"])
    ap.add_argument("--canonical", default="Data/FVC2000_binary_regen"); ap.add_argument("--raw", default="Data/FVC2000_raw")
    ap.add_argument("--work", default="Data/FVC2000_raw_degraded"); ap.add_argument("--binaries", default="Data/FVC2000_binary_degraded")
    ap.add_argument("--family"); ap.add_argument("--strengths"); ap.add_argument("--conditions"); ap.add_argument("--dirs")
    ap.add_argument("--file"); ap.add_argument("--keep")
    ap.add_argument("--nfiq", default="/tmp/nbis_install/bin/nfiq"); ap.add_argument("--image", default="1_1")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2)))
    A = ap.parse_args(); Path(A.work).mkdir(parents=True, exist_ok=True)
    {"benchmark": stage_benchmark, "scan": stage_scan, "select": stage_select, "verify": stage_verify, "distribution": stage_distribution, "nfiq": stage_nfiq, "subset": stage_subset, "montage": stage_montage}[A.stage](A)


if __name__ == "__main__":
    main()
