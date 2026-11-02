#!/usr/bin/env python3
"""
Acquisition degradations of the sensitivity experiment (R4-10), applied to the RAW grayscale FVC2000 DB1
images BEFORE enhancement. The degraded raw images are then passed through
the canonical enhancer (Code/preprocess.py, library defaults), so the resulting binaries reflect how
the pipeline actually responds to acquisition degradation, cleanup included.

Input : Data/FVC2000_raw/<f>_<i>.tif        (300 x 300, 8-bit, dark ridges on a light background)
Output: Data/FVC2000_raw_degraded/<cond>/<f>_<i>.tif   (same names, 8-bit uncompressed TIFF)
        + degrade_manifest.json (family, strength, calibration record, per-image input/output md5,
          pixel sha256, seed, environment, script hash)

Families (strength s; I = raw image as float, bg = 99th-percentile intensity = background level):
  dry    light pressure / dry skin: local contrast loss.  I' = bg - (bg - I) * max(0, 1 - s * F(x)),
         F = smooth random field (Gaussian-filtered white noise, sigma 20 px, min-max scaled to [0, 1])
         -> ridge contrast attenuated patchily, fully where F(x) >= 1/s: for s <= 1 only the worst spot
         loses all contrast, for s in (1, 2] the fully faded patch grows (at s = 2, half of the field).
  wet    moisture / heavy pressure: ridges thicken and valleys fill.  Grayscale min-filter (3 x 3 cross;
         ridges are dark, so the minimum spreads ink) blended with the original: lambda = min(s, 1), then a
         second min-filter blend with lambda = s - 1 if s > 1; then Gaussian blur sigma = 0.5 * s.  s in [0, 2].
  noise  sensor noise: additive Gaussian noise with standard deviation s gray levels (0-255 scale).
  blur   defocus / motion smear: Gaussian blur with sigma = s pixels.
Each result is rounded and clipped to 0..255 (uint8). One RNG per (image, condition), seeded from
sha256(file name | condition), so results are reproducible and independent of the processing order.

Strengths are NOT free parameters: Code/sensitivity/calibrate_degradations.py chooses them by its rule
(the degraded image may differ from its clean version, in the monitored binary statistics, by no more than
two impressions of the same finger typically differ: the median within-finger difference; condition
<family>_typ) and writes Data/FVC2000_raw_degraded/calibration.json. This generator refuses to run without
that file and records its hash.

Usage (project root):
  python3 Code/sensitivity/degrade_raw.py --input Data/FVC2000_raw --out Data/FVC2000_raw_degraded \
      --calibration Data/FVC2000_raw_degraded/calibration.json [--conditions dry_typ,wet_typ,...]
Then, per condition:  python3 Code/preprocess.py --input Data/FVC2000_raw_degraded/<cond> --output Data/FVC2000_binary_degraded/<cond>
"""
import argparse, hashlib, json, os, platform, re, sys, time
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

FAMILIES = ("dry", "wet", "noise", "blur")
CROSS3 = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=bool)
FIELD_SIGMA = 20.0            # correlation scale of the dry-patch field (px)


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def smooth_field(shape, rng, sigma=FIELD_SIGMA):
    """Gaussian-filtered white noise, min-max scaled to [0, 1]."""
    f = ndi.gaussian_filter(rng.standard_normal(shape), sigma)
    return (f - f.min()) / max(f.max() - f.min(), 1e-12)


def to_uint8(x):
    return np.clip(np.rint(x), 0, 255).astype(np.uint8)


def degrade(img_u8, family, s, rng):
    """img_u8: raw grayscale uint8 (dark ridges, light background). Returns uint8 of the same shape."""
    I = img_u8.astype(np.float64)
    if family == "dry":
        if not 0.0 <= s <= 2.0: raise ValueError("dry strength must be in [0, 2]")
        bg = np.percentile(I, 99)
        F = smooth_field(I.shape, rng)
        out = bg - (bg - I) * np.maximum(0.0, 1.0 - s * F)
    elif family == "wet":
        if not 0.0 <= s <= 2.0: raise ValueError("wet strength must be in [0, 2]")
        lam1 = min(s, 1.0)
        out = (1 - lam1) * I + lam1 * ndi.minimum_filter(I, footprint=CROSS3, mode="nearest")
        if s > 1.0:
            lam2 = s - 1.0
            out = (1 - lam2) * out + lam2 * ndi.minimum_filter(out, footprint=CROSS3, mode="nearest")
        if s > 0: out = ndi.gaussian_filter(out, 0.5 * s)
    elif family == "noise":
        if s < 0: raise ValueError("noise strength must be >= 0")
        out = I + rng.normal(0.0, s, size=I.shape)
    elif family == "blur":
        if s < 0: raise ValueError("blur strength must be >= 0")
        out = ndi.gaussian_filter(I, s) if s > 0 else I
    else:
        raise ValueError(f"unknown family {family}")
    return to_uint8(out)


def seed_for(name, cond):
    return int.from_bytes(hashlib.sha256(f"{name}|{cond}".encode()).digest()[:8], "little")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--calibration", required=True, help="calibration.json written by calibrate_degradations.py")
    ap.add_argument("--conditions", default=None, help="comma-separated subset of the calibrated conditions (default: all)")
    ap.add_argument("--overwrite", action="store_true")
    A = ap.parse_args()
    cal = json.load(open(A.calibration)); cal_sha = hashlib.sha256(open(A.calibration, "rb").read()).hexdigest()
    chosen = cal["chosen"]                                  # {cond: {"family":..., "strength":..., "level":...}}
    conds = [c.strip() for c in A.conditions.split(",")] if A.conditions else list(chosen)
    for c in conds:
        if c not in chosen: sys.exit(f"condition {c} not in calibration file")
    files = sorted(p for p in Path(A.input).iterdir() if re.fullmatch(r"\d+_\d+\.tif", p.name))
    if len(files) != 800: sys.exit(f"expected 800 raw images in {A.input}, found {len(files)}")
    script_sha = hashlib.sha256(open(__file__, "rb").read()).hexdigest()
    t0 = time.time()
    for cond in conds:
        fam, s, level = chosen[cond]["family"], float(chosen[cond]["strength"]), chosen[cond]["level"]
        od = Path(A.out) / cond
        if od.exists() and any(od.iterdir()) and not A.overwrite: sys.exit(f"{od} is not empty; use --overwrite")
        od.mkdir(parents=True, exist_ok=True)
        rows = []
        for p in files:
            a = np.array(Image.open(p))
            if a.ndim != 2 or a.dtype != np.uint8 or a.shape != (300, 300): sys.exit(f"{p.name}: not a 300x300 8-bit grayscale image")
            seed = seed_for(p.name, cond)
            out = degrade(a, fam, s, np.random.default_rng(seed))
            Image.fromarray(out, mode="L").save(od / p.name, format="TIFF", compression=None)
            back = np.array(Image.open(od / p.name))
            if not np.array_equal(back, out): sys.exit(f"{p.name}: TIFF round trip changed pixels")
            rows.append(dict(file=p.name, input_md5=md5(p), output_md5=md5(od / p.name),
                             pixel_sha256=hashlib.sha256(np.ascontiguousarray(out).tobytes()).hexdigest(),
                             mean_abs_change=float(np.abs(out.astype(float) - a.astype(float)).mean()), seed=seed))
        manifest = dict(created=time.strftime("%Y-%m-%dT%H:%M:%S%z"), condition=cond, family=fam, strength=s, level=level,
                        calibration_file=os.path.abspath(A.calibration), calibration_sha256=cal_sha,
                        input_dir=str(Path(A.input).resolve()), output_dir=str(od.resolve()), n_images=len(rows),
                        definition=__doc__.split("Families")[1].split("Each result")[0].strip(),
                        seeding="numpy default_rng(int.from_bytes(sha256(file name + '|' + condition)[:8], 'little')) per image",
                        environment=dict(python=platform.python_version(), platform=platform.platform(), numpy=np.__version__,
                                         scipy=__import__("scipy").__version__, pillow=Image.__version__),
                        script_sha256=script_sha,
                        summary=dict(mean_abs_change_median=float(np.median([r["mean_abs_change"] for r in rows]))),
                        images=rows)
        json.dump(manifest, open(od / "degrade_manifest.json", "w"), indent=1)
        print(f"{cond:>10} ({fam}, s = {s:g}, {level}): {len(rows)} images; mean |delta gray| median {manifest['summary']['mean_abs_change_median']:.2f}  ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
