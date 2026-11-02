#!/usr/bin/env python3
"""
preprocess.py — enhanced binary ridge/valley images from raw fingerprint images.

Reconstruction (revision step 0.3) of the procedure that produced Data/FVC2000_binary/
in January 2026 (notebook Binarized_FVC2000.ipynb, commits 8cbd957 / da3b7b3 of the
TDA-Fingerprint-Analysis repository), reduced to its three operative lines:

    M   = cv2.imread(path, 0)                             # 8-bit grayscale
    out = fingerprint_enhancer.enhance_fingerprint(M)     # library defaults, resize=False
    cv2.imwrite(stem + "_bin.png", (out * 255).astype(np.uint8))

enhance_fingerprint is the Hong-Wan-Jain enhancement (Kovesi's implementation) as
packaged in Fingerprint-Enhancement-Python 0.0.14 (Utkarsh Deshmukh): block-wise ridge
segmentation (16 px blocks, std threshold 0.1), orientation field (sigmas 1/7/7), a single
mean ridge frequency per image (38 px blocks, wavelengths 5-15 px), oriented Gabor
filtering (kx = ky = 0.65, 3 deg steps) and the threshold response < -3 => ridge.
The output is boolean with True = ridge; it is written as 0/255 PNG, ridge = 255.

Deterministic: no random numbers anywhere; identical output across the library versions
0.0.11-0.0.14 and across numpy 1.23-2.5 / scipy 1.9-1.18 / OpenCV 4.6-5.0 (verified on
FVC2000 DB1, 11 Oct 2026).

Usage (from the project root):
    python Code/preprocess.py --input Data/FVC2000_raw     --output Data/FVC2000_binary_regen
    python Code/preprocess.py --input Data/FVC2002_DB3_raw --output Data/FVC2002_DB3_binary

Writes <output>/<stem>_bin.png for every image matching --pattern (default *.tif), plus
<output>/preprocess_manifest.json recording the environment (python, platform, library
versions), the exact parameters, and per image the input md5, output md5, ridge fraction
and any warning. Never overwrites an existing output directory unless --overwrite is given.
Requires: fingerprint-enhancer==0.0.14, numpy, scipy, opencv-python(-headless)
(see Code/requirements-preprocess.txt).
"""
import argparse
import hashlib
import json
import os
import platform
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from glob import glob

import numpy as np

LIBRARY_PARAMS = dict(          # the library's defaults, spelled out so the manifest is explicit
    resize=False, ridge_segment_blksze=16, ridge_segment_thresh=0.1, gradient_sigma=1,
    block_sigma=7, orient_smooth_sigma=7, ridge_freq_blksze=38, ridge_freq_windsze=5,
    min_wave_length=5, max_wave_length=15, relative_scale_factor_x=0.65,
    relative_scale_factor_y=0.65, angle_inc=3.0, ridge_filter_thresh=-3, invert_output=False)


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def process_one(args):
    """Enhance one image. Returns a manifest row; never raises (errors are recorded).
    args = (in_path, out_path) or (in_path, out_path, params); params defaults to LIBRARY_PARAMS.
    (Revision step 1.4: parameter overrides travel in the job tuple so that they reach the
    workers under every multiprocessing start method; the default call is unchanged.)"""
    in_path, out_path = args[0], args[1]
    params = args[2] if len(args) > 2 else LIBRARY_PARAMS
    import cv2                      # imported in the worker
    import fingerprint_enhancer
    row = {"input": os.path.basename(in_path), "output": os.path.basename(out_path),
           "input_md5": md5_file(in_path)}
    t0 = time.time()
    try:
        M = cv2.imread(in_path, 0)                                   # IMREAD_GRAYSCALE
        if M is None:
            raise ValueError("cv2.imread returned None (not an image?)")
        row["input_shape"] = list(M.shape)
        row["input_mean"] = float(M.mean())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out = fingerprint_enhancer.enhance_fingerprint(M, **params)
        out = np.asarray(out)
        if out.dtype != bool:
            raise ValueError(f"unexpected output dtype {out.dtype} (expected bool)")
        if out.shape != M.shape:
            raise ValueError(f"output shape {out.shape} != input shape {M.shape}")
        binary_img = (out * 255).astype(np.uint8)                    # exactly as the notebook did
        if not cv2.imwrite(out_path, binary_img):
            raise IOError("cv2.imwrite failed")
        frac = float(out.mean())
        # pixel_sha256 hashes the decoded pixel array (uint8, row-major); PNG file md5s can
        # differ between encoders while the pixels are identical, so both are recorded.
        row.update({"ridge_fraction": frac, "output_md5": md5_file(out_path),
                    "pixel_sha256": sha256_bytes(np.ascontiguousarray(binary_img).tobytes()),
                    "seconds": round(time.time() - t0, 3), "status": "ok"})
        if frac > 0.5:
            row["warning"] = "ridge fraction > 0.5: polarity suspect (the R loader will invert)"
        elif frac < 0.05:
            row["warning"] = "ridge fraction < 0.05: segmentation found almost no print"
    except Exception as e:                                           # noqa: BLE001
        row.update({"status": "error", "error": f"{type(e).__name__}: {e}",
                    "seconds": round(time.time() - t0, 3)})
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="directory of raw images")
    ap.add_argument("--output", required=True, help="directory for <stem>_bin.png outputs")
    ap.add_argument("--pattern", default="*.tif", help="glob for raw images (default *.tif)")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--overwrite", action="store_true", help="allow writing into an existing output directory")
    ap.add_argument("--param", action="append", default=[], metavar="KEY=VALUE",
                    help="override one enhancement parameter (repeatable; revision step 1.4 sensitivity Tier B). "
                         "Without --param the run is the canonical procedure with the library defaults.")
    a = ap.parse_args()

    params = dict(LIBRARY_PARAMS)
    for kv in a.param:
        k, _, v = kv.partition("=")
        if k not in params: sys.exit(f"--param: unknown parameter {k!r}; known: {sorted(params)}")
        d = params[k]
        params[k] = (v.lower() == "true") if isinstance(d, bool) else type(d)(v) if not isinstance(d, bool) else d
    overrides = {k: params[k] for k in params if params[k] != LIBRARY_PARAMS[k]}

    import cv2
    import scipy
    import fingerprint_enhancer
    try:
        from importlib.metadata import version as _v
        fe_version = _v("fingerprint-enhancer")
    except Exception:                                                # noqa: BLE001
        fe_version = "unknown"
    # Hash the enhancement source that is actually executed (a version number alone cannot
    # identify locally modified code) and this script itself.
    try:
        import fingerprint_enhancer.fingerprint_image_enhancer as _fie
        enhancer_source_file = _fie.__file__
        enhancer_source_sha256 = sha256_bytes(open(enhancer_source_file, "rb").read())
    except Exception as e:                                           # noqa: BLE001
        enhancer_source_file, enhancer_source_sha256 = None, f"unavailable: {e}"
    this_script_sha256 = sha256_bytes(open(os.path.abspath(__file__), "rb").read())

    files = sorted(f for f in glob(os.path.join(a.input, a.pattern)) if " 2." not in os.path.basename(f))
    if not files:
        sys.exit(f"no files matching {a.pattern} in {a.input}")
    if os.path.isdir(a.output) and os.listdir(a.output) and not a.overwrite:
        sys.exit(f"{a.output} exists and is not empty; use --overwrite to write into it")
    os.makedirs(a.output, exist_ok=True)

    jobs = [(f, os.path.join(a.output, os.path.splitext(os.path.basename(f))[0] + "_bin.png")) + ((params,) if overrides else ())
            for f in files]
    if overrides:
        print(f"preprocess.py: parameter overrides (sensitivity variant, NOT the canonical procedure): {overrides}")
    print(f"preprocess.py: {len(jobs)} images from {a.input} -> {a.output}  (fingerprint-enhancer {fe_version}, "
          f"numpy {np.__version__}, scipy {scipy.__version__}, opencv {cv2.__version__}, {a.workers} workers)")
    t0 = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for i, row in enumerate(ex.map(process_one, jobs, chunksize=4), 1):
            rows.append(row)
            if i % 100 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} done ({time.time() - t0:.0f} s)")

    n_err = sum(r["status"] != "ok" for r in rows)
    n_warn = sum("warning" in r for r in rows)
    fracs = [r["ridge_fraction"] for r in rows if r["status"] == "ok"]
    manifest = {
        "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "procedure": ("cv2.imread(path, 0) -> fingerprint_enhancer.enhance_fingerprint(M, defaults) -> (out*255).astype(uint8) -> PNG"
                      if not overrides else
                      "cv2.imread(path, 0) -> fingerprint_enhancer.enhance_fingerprint(M, parameters below; sensitivity variant) -> (out*255).astype(uint8) -> PNG"),
        "parameters": params,
        "parameter_overrides": overrides,
        "environment": {"python": sys.version.split()[0], "platform": platform.platform(),
                        "machine": platform.machine(), "fingerprint_enhancer": fe_version,
                        "fingerprint_enhancer_file": getattr(fingerprint_enhancer, "__file__", None),
                        "enhancer_source_file": enhancer_source_file,
                        "enhancer_source_sha256": enhancer_source_sha256,
                        "preprocess_py_sha256": this_script_sha256,
                        "numpy": np.__version__, "scipy": scipy.__version__, "opencv": cv2.__version__},
        "input_dir": os.path.abspath(a.input), "output_dir": os.path.abspath(a.output), "pattern": a.pattern,
        "n_images": len(rows), "n_errors": n_err, "n_warnings": n_warn,
        "ridge_fraction_summary": {"min": min(fracs), "median": float(np.median(fracs)), "max": max(fracs)} if fracs else None,
        "elapsed_seconds": round(time.time() - t0, 1),
        "images": rows,
    }
    with open(os.path.join(a.output, "preprocess_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    print(f"done: {len(rows) - n_err} ok, {n_err} errors, {n_warn} warnings; ridge fraction "
          f"{manifest['ridge_fraction_summary']}; manifest written to {a.output}/preprocess_manifest.json")
    if n_err:
        print("ERRORS:")
        for r in rows:
            if r["status"] != "ok":
                print("  ", r["input"], r["error"])
        sys.exit(1)


if __name__ == "__main__":
    main()
