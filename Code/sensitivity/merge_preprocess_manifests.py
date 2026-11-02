#!/usr/bin/env python3
"""
Merge the preprocess_manifest.json files of several preprocess.py chunk runs (same input directory,
same parameters, disjoint file subsets) into one manifest, and move the chunk outputs into one
directory. Needed only where a single preprocess.py run cannot complete within a per-command time
limit (the sandbox); the merged manifest is identical in content to a one-shot run's, plus a
`chunks` record. Revision step 1.4 (Tier B variants).

  python3 Code/sensitivity/merge_preprocess_manifests.py --chunks <dir1> <dir2> ... --out <final dir> [--expect 800]
"""
import argparse, json, os, shutil, sys
from pathlib import Path

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--chunks", nargs="+", required=True); ap.add_argument("--out", required=True); ap.add_argument("--expect", type=int, default=800)
A = ap.parse_args()
mans = [json.load(open(Path(c) / "preprocess_manifest.json")) for c in A.chunks]
keys_same = ["procedure", "parameters", "parameter_overrides", "input_dir"]   # pattern differs by chunk, by design
for k in keys_same + ["environment"]:
    vals = [json.dumps(m[k], sort_keys=True) for m in mans]
    if len(set(vals)) != 1: sys.exit(f"chunks disagree on {k}")
rows = [r for m in mans for r in m["images"]]
names = [r["output"] for r in rows]
if len(names) != len(set(names)): sys.exit("duplicate outputs across chunks")
if len(rows) != A.expect: sys.exit(f"expected {A.expect} images across chunks, found {len(rows)}")
out = Path(A.out); out.mkdir(parents=True, exist_ok=True)
if any(p.suffix == ".png" for p in out.iterdir()): sys.exit(f"{out} already contains PNGs")
for c in A.chunks:
    for p in Path(c).glob("*_bin.png"): shutil.move(str(p), str(out / p.name))
rows.sort(key=lambda r: (int(r["output"].split("_")[0]), int(r["output"].split("_")[1])))
fr = [r["ridge_fraction"] for r in rows if r["status"] == "ok"]
merged = dict(mans[0]); merged.pop("images")
merged.update(dict(created=max(m["created"] for m in mans), output_dir=str(out.resolve()), pattern="*.tif (in chunks)", n_images=len(rows),
                   n_errors=sum(r["status"] != "ok" for r in rows), n_warnings=sum("warning" in r for r in rows),
                   ridge_fraction_summary=dict(min=min(fr), median=float(np.median(fr)), max=max(fr)),
                   elapsed_seconds=round(sum(m["elapsed_seconds"] for m in mans), 1),
                   chunks=[dict(pattern=m["pattern"], n_images=m["n_images"], created=m["created"]) for m in mans], images=rows))
json.dump(merged, open(out / "preprocess_manifest.json", "w"), indent=1)
for c in A.chunks: shutil.rmtree(c)
print(f"merged {len(mans)} chunks -> {out}: {len(rows)} images, {merged['n_errors']} errors, overrides {merged['parameter_overrides']}")
