#!/usr/bin/env python3
"""
Revision step 1.3 (C3, R4-7) - classical reference matcher: NIST NBIS 5.0.0 mindtct + bozorth3
on every pair of images of a database, written as a score table with a manifest.

Procedure (pre-specified in STEP_1.3_REFERENCE_MATCHER_SPEC.md):
  grayscale kind : the raw FVC TIFF (8-bit, 300x300, 500 ppi) -> raw bytes -> cjpegl (LOSSLESS JPEG, the
                   container mindtct reads) -> mindtct (standard NIST conventions: no -b, no -m1) -> .xyt
                   (x, y, theta, quality) ; nfiq quality (1 best .. 5 worst) recorded per image
  binary kind    : our own binary ridge map <f>_<i>_bin.png (ridge = 255) INVERTED so that ridges are
                   dark as in a grayscale print, then the same chain (secondary row; not the matcher's
                   design condition)
  scoring        : bozorth3 with its defaults (-n 150 minutiae max; minminutiae = 10, below which the
                   score is 0), every probe against the full gallery (-p probe -G gallery.lis), so BOTH
                   orderings of every unordered pair are scored; the pre-specified pair score is the
                   MEAN of the two orderings (bozorth3 is not exactly symmetric); both are kept.

Usage (from the project root):
  python3 Code/reference_matcher/run_nbis_matcher.py --images Data/FVC2000_raw --kind grayscale \
      --out Outputs_reference/nbis_db1_gray --nbis-bin /path/to/nbis/bin [--workers 4] [--fingers 10]
Outputs: images.csv, scores_pairs.csv (i<j), self_scores.csv, xyt/ (templates), scores_chunks/ (raw per-probe
scores), manifest.json. Refuses to run into a non-empty output directory unless --overwrite.
Stages: --stage all (default) does everything; on a machine with a per-command time limit use
  --stage extract, then --stage score --probes 1-200 (etc., any partition of 1..N), then --stage merge.
The merge checks that every probe was scored exactly once; the written tables are identical either way.
"""
import argparse, csv, hashlib, json, os, platform, re, subprocess, sys, tempfile, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from PIL import Image

NBIS_VERSION = "5.0.0"
NBIS_ZIP_SHA256 = "0adf8ab0f6b0e4208de50ca00ba21d3d77112ecd66288757ddfed21f6bee92c3"
BUILD_NOTE = "NIST nbis_v5_0_0.zip; setup.sh --without-X11 --STDLIBS --64; rules.mak ARCH_FLAG := -fPIC -fcommon (no source changes)"

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def parse_name(path, kind):
    m = re.match(r"^(\d+)_(\d+)" + (r"_bin\.png$" if kind == "binary" else r"\.tif$"), path.name)
    return (int(m.group(1)), int(m.group(2))) if m else None

def run(cmd, cwd=None):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed ({r.returncode}): {r.stderr.strip()[:300]}")
    return r.stdout

def extract_one(args):
    """TIFF/PNG -> raw -> JPEGL -> mindtct (+ nfiq for grayscale). Returns a dict for images.csv."""
    idx, path, kind, xyt_dir, nbis = args
    path = Path(path); stem = path.name.replace("_bin.png", "").replace(".tif", "")
    im = Image.open(path); a = np.array(im)
    if a.ndim != 2 or a.dtype != np.uint8 or a.shape != (300, 300):
        raise RuntimeError(f"{path.name}: expected 300x300 8-bit grayscale, got {a.shape} {a.dtype}")
    if kind == "binary":
        vals = np.unique(a)
        if not set(vals.tolist()) <= {0, 255}:
            raise RuntimeError(f"{path.name}: binary image must contain only 0/255, got {vals[:5]}")
        a = (255 - a).astype(np.uint8)            # ridge (255 in our maps) -> dark, as in a grayscale print
    with tempfile.TemporaryDirectory() as td:
        raw = Path(td) / f"{stem}.raw"; a.tofile(raw)
        run([f"{nbis}/cjpegl", "jpl", str(raw), "-raw_in", "300,300,8,500"], cwd=td)
        jpl = Path(td) / f"{stem}.jpl"
        oroot = Path(xyt_dir) / stem
        run([f"{nbis}/mindtct", str(jpl), str(oroot)])
        nfiq = run([f"{nbis}/nfiq", str(jpl)]).strip() if kind == "grayscale" else ""
    xyt = Path(xyt_dir) / f"{stem}.xyt"
    n_min = sum(1 for line in open(xyt) if line.strip())
    for ext in (".min", ".dm", ".hcm", ".lcm", ".lfm", ".qm", ".brw"):   # keep only .xyt and .min
        p = Path(xyt_dir) / f"{stem}{ext}"
        if ext != ".min" and p.exists(): p.unlink()
    return dict(index=idx, file=path.name, stem=stem, md5=md5(path), n_minutiae=n_min, nfiq=nfiq)

def score_probe(args):
    """bozorth3: one probe against the whole gallery; returns {gallery_stem: score}."""
    probe_xyt, gallery_lis, nbis = args
    out = run([f"{nbis}/bozorth3", "-A", "outfmt=pgs", "-p", probe_xyt, "-G", gallery_lis])
    res = {}
    for line in out.splitlines():
        parts = line.split()
        if len(parts) != 3: continue
        res[Path(parts[1]).stem] = int(parts[2])
    return Path(probe_xyt).stem, res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--images", required=True); ap.add_argument("--kind", choices=["grayscale", "binary"], required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--nbis-bin", required=True)
    ap.add_argument("--workers", type=int, default=4); ap.add_argument("--fingers", type=int, default=100, help="use fingers 1..N (dry runs)")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--stage", choices=["all", "extract", "score", "merge"], default="all",
                    help="all = extract + score + merge in one go; or run the stages separately (score in probe chunks) on a "
                         "machine with a per-command time limit")
    ap.add_argument("--probes", default=None, help="for --stage score: 1-based inclusive probe range a-b (default: all)")
    A = ap.parse_args()
    t0 = time.time()
    out = Path(A.out); xyt_dir = out / "xyt"
    if A.stage in ("all", "extract") and out.exists() and any(out.iterdir()) and not A.overwrite:
        sys.exit(f"{out} is not empty; use --overwrite")
    nbis = str(Path(A.nbis_bin).resolve())
    for t in ("mindtct", "bozorth3", "nfiq", "cjpegl"):
        if not os.access(f"{nbis}/{t}", os.X_OK): sys.exit(f"missing NBIS tool {t} in {nbis}")
    out.mkdir(parents=True, exist_ok=True); xyt_dir.mkdir(exist_ok=True)      # only once the tools are known to exist

    files = []
    for p in sorted(Path(A.images).iterdir()):
        k = parse_name(p, A.kind)
        if k and 1 <= k[0] <= A.fingers and 1 <= k[1] <= 8: files.append((k, p))
    files.sort()
    expected = 8 * A.fingers
    if len(files) != expected: sys.exit(f"expected {expected} images for fingers 1..{A.fingers}, found {len(files)}")
    print(f"{len(files)} {A.kind} images from {A.images}; NBIS {NBIS_VERSION} at {nbis}; {A.workers} workers")

    chunk_dir = out / "scores_chunks"
    tfile = out / "stage_timings.json"
    def note_timing(stage, seconds):
        t = json.load(open(tfile)) if tfile.exists() else []
        t.append({"stage": stage, "seconds": round(seconds, 1), "at": time.strftime("%Y-%m-%dT%H:%M:%S%z")})
        json.dump(t, open(tfile, "w"))
    # ---------------------------------------------------------------- stage 1: minutiae extraction
    if A.stage in ("all", "extract"):
        # (QA, 13 Oct) a restarted extraction invalidates every earlier score chunk and template: remove
        # them, so a later merge can never combine new templates with old scores.
        import shutil
        for stale in (chunk_dir, xyt_dir):
            if stale.exists(): shutil.rmtree(stale)
        xyt_dir.mkdir()
        for stale in ("scores_pairs.csv", "self_scores.csv", "manifest.json", "images.csv", "gallery.lis", "stage_timings.json"):
            if (out / stale).exists(): (out / stale).unlink()
        t1 = time.time()
        with ProcessPoolExecutor(A.workers) as ex:
            rows = list(ex.map(extract_one, [(i + 1, str(p), A.kind, str(xyt_dir), nbis) for i, (k, p) in enumerate(files)]))
        for r, (k, p) in zip(rows, files): r["finger"], r["impression"] = k
        with open(out / "images.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["index", "file", "stem", "finger", "impression", "md5", "n_minutiae", "nfiq"]); w.writeheader(); w.writerows(rows)
        (out / "gallery.lis").write_text("".join(f"{xyt_dir / r['stem']}.xyt\n" for r in rows))
        n_min = np.array([r["n_minutiae"] for r in rows])
        print(f"  extraction done in {time.time() - t1:.0f} s: minutiae per image min {n_min.min()} median {int(np.median(n_min))} max {n_min.max()}; "
              f"{int((n_min < 10).sum())} image(s) below bozorth3's minminutiae = 10")
        note_timing("extract", time.time() - t1)
        if A.stage == "extract":
            print(f"extract stage done in {time.time() - t0:.0f} s -> {out}"); return
    rows = list(csv.DictReader(open(out / "images.csv")))
    for r in rows:
        r["index"], r["finger"], r["impression"], r["n_minutiae"] = int(r["index"]), int(r["finger"]), int(r["impression"]), int(r["n_minutiae"])
    stems = [r["stem"] for r in rows]; n = len(stems)
    if n != len(files) or [r["file"] for r in rows] != [p.name for k, p in files]: sys.exit("images.csv does not match the image set")
    gallery_lis = out / "gallery.lis"

    # ---------------------------------------------------------------- stage 2: all-pairs scoring (both orderings), optionally in probe chunks
    if A.stage in ("all", "score"):
        a, b = (1, n) if A.probes is None else (int(A.probes.split("-")[0]), int(A.probes.split("-")[1]))
        if not (1 <= a <= b <= n): sys.exit(f"bad --probes range {A.probes}")
        t2 = time.time()
        with ProcessPoolExecutor(A.workers) as ex:
            results = dict(ex.map(score_probe, [(f"{xyt_dir / s_}.xyt", str(gallery_lis), nbis) for s_ in stems[a - 1:b]]))
        for s_ in stems[a - 1:b]:
            if set(results[s_].keys()) != set(stems): sys.exit(f"probe {s_}: incomplete gallery scores")
        chunk_dir.mkdir(exist_ok=True)
        json.dump({"probes": f"{a}-{b}", "scores": {s_: [results[s_][g] for g in stems] for s_ in stems[a - 1:b]}}, open(chunk_dir / f"probes_{a:04d}-{b:04d}.json", "w"))
        print(f"  scored probes {a}-{b} against all {n} gallery templates in {time.time() - t2:.0f} s")
        note_timing(f"score {a}-{b}", time.time() - t2)
        if A.stage == "score":
            print(f"score stage done in {time.time() - t0:.0f} s -> {chunk_dir}"); return

    # ---------------------------------------------------------------- stage 3: merge chunks, write tables and manifest
    timings = json.load(open(tfile)) if tfile.exists() else []
    S = np.full((n, n), -1, dtype=np.int64); covered = np.zeros(n, dtype=bool); chunks = sorted(chunk_dir.glob("probes_*.json"))
    for cf in chunks:
        d = json.load(open(cf))
        for s_, sc_row in d["scores"].items():
            i = stems.index(s_)
            if covered[i]: sys.exit(f"probe {s_} appears in more than one chunk")
            S[i, :] = sc_row; covered[i] = True
    if not covered.all(): sys.exit(f"{int((~covered).sum())} probe(s) not scored; chunks present: {[c.name for c in chunks]}")
    if (S < 0).any(): sys.exit("missing scores in the bozorth3 output")
    finger = np.array([r["finger"] for r in rows]); n_min = np.array([r["n_minutiae"] for r in rows])
    with open(out / "scores_pairs.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["i", "j", "file_i", "file_j", "finger_i", "finger_j", "genuine", "score_ij", "score_ji", "score_mean"])
        for i in range(n):
            for j in range(i + 1, n):
                w.writerow([i + 1, j + 1, rows[i]["file"], rows[j]["file"], finger[i], finger[j], int(finger[i] == finger[j]),
                            S[i, j], S[j, i], (S[i, j] + S[j, i]) / 2])
    with open(out / "self_scores.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["i", "file", "self_score"]); [w.writerow([i + 1, rows[i]["file"], S[i, i]]) for i in range(n)]
    iu = np.triu_indices(n, 1)
    asym = np.abs(S[iu] - S.T[iu])
    gen = finger[iu[0]] == finger[iu[1]]
    Sm = (S[iu] + S.T[iu]) / 2                      # the pair score: mean of the two orderings (as written to scores_pairs.csv)
    print(f"  merged {len(chunks)} chunk file(s): {len(iu[0])} unordered pairs ({int(gen.sum())} genuine); "
          f"orderings identical in {np.mean(asym == 0) * 100:.1f} % of pairs, max |difference| {asym.max()}; "
          f"mean pair score genuine {Sm[gen].mean():.1f} vs impostor {Sm[~gen].mean():.1f}")

    manifest = dict(
        created=time.strftime("%Y-%m-%dT%H:%M:%S%z"), kind=A.kind, images_dir=str(Path(A.images).resolve()), n_images=n, fingers=A.fingers,
        procedure=("raw 8-bit -> cjpegl (lossless) -> mindtct (standard conventions, no -b, no -m1) -> bozorth3 defaults "
                   "(-n 150, minminutiae 10), every probe vs full gallery (-p/-G); pair score = mean of the two orderings"
                   + ("; binary maps inverted (ridge dark) before extraction" if A.kind == "binary" else "; nfiq per image")),
        scoring_chunks=[c.name for c in chunks],
        nbis=dict(version=NBIS_VERSION, source_zip_sha256=NBIS_ZIP_SHA256, build=BUILD_NOTE, bin=nbis,
                  tool_md5={t: md5(f"{nbis}/{t}") for t in ("mindtct", "bozorth3", "nfiq", "cjpegl")}),
        environment=dict(python=platform.python_version(), platform=platform.platform(), machine=platform.machine(),
                         numpy=np.__version__, pillow=Image.__version__),
        script_sha256=hashlib.sha256(open(__file__, "rb").read()).hexdigest(),
        minutiae=dict(min=int(n_min.min()), median=float(np.median(n_min)), max=int(n_min.max()), n_below_10=int((n_min < 10).sum())),
        pairs=dict(n_unordered=int(len(iu[0])), n_genuine=int(gen.sum()), frac_orderings_identical=float(np.mean(asym == 0)), max_ordering_difference=int(asym.max()),
                   mean_genuine=float(Sm[gen].mean()), mean_impostor=float(Sm[~gen].mean())),
        outputs={f: md5(out / f) for f in ("images.csv", "scores_pairs.csv", "self_scores.csv")},
        elapsed_seconds_this_invocation=round(time.time() - t0, 1),
        stage_timings=timings)
    json.dump(manifest, open(out / "manifest.json", "w"), indent=2)
    print(f"done in {time.time() - t0:.0f} s -> {out}")

if __name__ == "__main__":
    main()
