#!/usr/bin/env python3
"""
check_binary_set.py — completeness and plausibility checks on a directory of enhanced binary
images produced by Code/preprocess.py (revision step 0.3; written for FVC2002 DB3, works for
any FVC "A" set: 100 fingers x 8 impressions).

Hard requirements (any failure -> exit 1):
  * exactly the expected filenames <finger>_<impression>_bin.png for fingers 1..n_fingers and
    impressions 1..n_impressions, nothing missing, nothing extra;
  * every file decodes, has the expected dimensions, and contains only the values {0, 255};
  * preprocess_manifest.json present, n_errors == 0, one 'ok' row per image.
Review flags (never discard; listed for human inspection):
  * ridge fraction outside [--flag-low, --flag-high] (default 0.10-0.45);
  * manifest warnings.
Writes review_list.csv and two contact sheets: a fixed random sample (seed 20261011, raw next
to binary) and all flagged images, so segmentation quality can be judged by eye across the set.

Usage (project root):
    python Code/check_binary_set.py --binary Data/FVC2002_DB3_binary --raw Data/FVC2002_DB3_raw --out Outputs_preprocess_check/db3
"""
import argparse
import csv
import json
import os
import re
import sys

import numpy as np
import cv2


def contact_sheet(pairs, path, tile=300, cols=4):
    """pairs: list of (label, raw_gray or None, binary) -> grid image with raw|binary per cell."""
    if not pairs:
        return
    rows = int(np.ceil(len(pairs) / cols))
    sheet = np.full((rows * (tile + 24), cols * (2 * tile + 8), 3), 40, np.uint8)
    for k, (label, raw, b) in enumerate(pairs):
        r, c = divmod(k, cols)
        y0, x0 = r * (tile + 24), c * (2 * tile + 8)
        cell = np.zeros((tile, 2 * tile, 3), np.uint8)
        if raw is not None:
            cell[:, :tile] = np.stack([cv2.resize(raw, (tile, tile), interpolation=cv2.INTER_AREA)] * 3, -1)
        cell[:, tile:] = np.stack([cv2.resize(b, (tile, tile), interpolation=cv2.INTER_NEAREST)] * 3, -1)
        sheet[y0 + 24:y0 + 24 + tile, x0:x0 + 2 * tile] = cell
        cv2.putText(sheet, label, (x0 + 4, y0 + 18), cv2.FONT_HERSHEY_SIMPLEX, 0.55, (255, 255, 255), 1, cv2.LINE_AA)
    cv2.imwrite(path, sheet)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", required=True)
    ap.add_argument("--raw", default=None, help="raw image directory (for the contact sheets)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-fingers", type=int, default=100)
    ap.add_argument("--n-impressions", type=int, default=8)
    ap.add_argument("--height", type=int, default=300)
    ap.add_argument("--width", type=int, default=300)
    ap.add_argument("--flag-low", type=float, default=0.10)
    ap.add_argument("--flag-high", type=float, default=0.45)
    ap.add_argument("--sample", type=int, default=16)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    expected = {f"{f}_{i}_bin.png" for f in range(1, a.n_fingers + 1) for i in range(1, a.n_impressions + 1)}
    present = {n for n in os.listdir(a.binary) if n.endswith("_bin.png")}
    failures = []
    missing, extra = sorted(expected - present), sorted(present - expected)
    if missing:
        failures.append(f"{len(missing)} expected files missing, e.g. {missing[:5]}")
    if extra:
        failures.append(f"{len(extra)} unexpected files, e.g. {extra[:5]}")

    man_path = os.path.join(a.binary, "preprocess_manifest.json")
    manifest = None
    if not os.path.exists(man_path):
        failures.append("preprocess_manifest.json missing")
    else:
        manifest = json.load(open(man_path))
        if manifest.get("n_errors", 1) != 0:
            failures.append(f"manifest reports {manifest.get('n_errors')} processing errors")
        ok_rows = {r["output"] for r in manifest["images"] if r.get("status") == "ok"}
        if ok_rows != expected:
            failures.append(f"manifest 'ok' rows ({len(ok_rows)}) do not match the expected file set")

    rows, flagged = [], []
    for name in sorted(expected & present, key=lambda s: [int(x) for x in re.findall(r"\d+", s)]):
        img = cv2.imread(os.path.join(a.binary, name), cv2.IMREAD_UNCHANGED)
        row = {"image": name}
        if img is None:
            failures.append(f"{name}: does not decode"); row["status"] = "undecodable"; rows.append(row); continue
        if img.shape != (a.height, a.width):
            failures.append(f"{name}: shape {img.shape}, expected {(a.height, a.width)}")
        vals = np.unique(img)
        if not set(vals.tolist()) <= {0, 255}:
            failures.append(f"{name}: values {vals[:6]} not in {{0, 255}}")
        frac = float((img == 255).mean())
        row.update({"status": "ok", "ridge_fraction": round(frac, 4), "flag": ""})
        if frac < a.flag_low or frac > a.flag_high:
            row["flag"] = f"ridge fraction {frac:.3f} outside [{a.flag_low}, {a.flag_high}] — review"
        if manifest is not None:
            mrow = next((r for r in manifest["images"] if r["output"] == name), None)
            if mrow and mrow.get("warning"):
                row["flag"] = (row["flag"] + "; " if row["flag"] else "") + "manifest: " + mrow["warning"]
        if row["flag"]:
            flagged.append(name)
        rows.append(row)

    with open(os.path.join(a.out, "review_list.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["image", "status", "ridge_fraction", "flag"])
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in w.fieldnames})

    def load_pair(name):
        b = cv2.imread(os.path.join(a.binary, name), 0)
        raw = None
        if a.raw:
            rp = os.path.join(a.raw, name.replace("_bin.png", ".tif"))
            raw = cv2.imread(rp, 0) if os.path.exists(rp) else None
        return raw, b

    rng = np.random.default_rng(20261011)
    names_ok = [r["image"] for r in rows if r.get("status") == "ok"]
    sample = sorted(rng.choice(names_ok, size=min(a.sample, len(names_ok)), replace=False).tolist(),
                    key=lambda s: [int(x) for x in re.findall(r"\d+", s)])
    contact_sheet([(n.replace("_bin.png", ""),) + load_pair(n) for n in sample], os.path.join(a.out, "contact_sheet_random_sample.png"))
    fl_rows = {r["image"]: r for r in rows}
    contact_sheet([(f"{n.replace('_bin.png','')}  rf={fl_rows[n]['ridge_fraction']}",) + load_pair(n) for n in flagged],
                  os.path.join(a.out, "contact_sheet_flagged.png"))

    fracs = np.array([r["ridge_fraction"] for r in rows if r.get("status") == "ok"])
    lines = [
        f"files: expected {len(expected)}, present {len(present)}, missing {len(missing)}, extra {len(extra)}",
        f"manifest: {'present' if manifest else 'MISSING'}; errors: {manifest.get('n_errors') if manifest else 'n/a'}; warnings: {manifest.get('n_warnings') if manifest else 'n/a'}",
        "ridge fraction: min %.3f  p5 %.3f  median %.3f  p95 %.3f  max %.3f" % tuple(np.percentile(fracs, [0, 5, 50, 95, 100])) if fracs.size else "ridge fraction: n/a",
        f"flagged for review (not excluded): {len(flagged)}  -> review_list.csv, contact_sheet_flagged.png",
        f"random sample contact sheet: {len(sample)} images -> contact_sheet_random_sample.png",
        "HARD FAILURES: " + (str(len(failures)) if failures else "none"),
    ] + ["  - " + x for x in failures]
    with open(os.path.join(a.out, "check_summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
