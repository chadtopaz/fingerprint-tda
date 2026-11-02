# Data

This directory is intentionally empty in the repository except for this file. The study uses two third-party benchmarks that we do not redistribute; obtain them, lay them out as below, and run the preprocessing script. Nothing under `Data/` other than this README is tracked (`.gitignore`); the released manifests under `Results/Data/` let you check your preprocessing against ours.

## Benchmarks

- **FVC2000 DB1** (set A: 100 fingers × 8 impressions = 800 images, optical sensor, 300 × 300 pixels, 500 dpi). Maio, D., Maltoni, D., Cappelli, R., Wayman, J. L., and Jain, A. K. (2002). FVC2000: Fingerprint verification competition. *IEEE TPAMI* 24(3), 402–412.
- **FVC2002 DB3** (set A: 100 × 8, capacitive sensor, 300 × 300 pixels). Maio, D., Maltoni, D., Cappelli, R., Wayman, J. L., and Jain, A. K. (2002). FVC2002: Second fingerprint verification competition. *Proc. ICPR*, 811–814.

The full 100-finger sets are distributed on the DVD accompanying the *Handbook of Fingerprint Recognition* (Maltoni, Maio, Jain, and Feng; Springer, 3rd ed., 2022); the competition sites (<http://bias.csr.unibo.it/fvc2000/>, <http://bias.csr.unibo.it/fvc2002/>) host the competition details and the small free subsets (set B, 10 fingers), which are *not* the sets used in the paper.

## Layout

```
Data/FVC2000_raw/<finger>_<impression>.tif       1_1.tif ... 100_8.tif   (800 raw 8-bit grayscale TIFFs)
Data/FVC2002_DB3_raw/<finger>_<impression>.tif   same naming
```

## Preprocessing (the binary maps the pipeline reads)

`Code/preprocess.py` applies the Hong–Wan–Jain enhancement as packaged in Fingerprint-Enhancement-Python 0.0.14 (library defaults) and writes 0/255 PNGs, ridge = 255, plus a manifest:

```bash
python3 -m venv .venv-preprocess && .venv-preprocess/bin/pip install -r Code/requirements-preprocess.txt
.venv-preprocess/bin/python Code/preprocess.py --input Data/FVC2000_raw     --output Data/FVC2000_binary_regen
.venv-preprocess/bin/python Code/preprocess.py --input Data/FVC2002_DB3_raw --output Data/FVC2002_DB3_binary
.venv-preprocess/bin/python Code/check_binary_set.py --binary Data/FVC2000_binary_regen --raw Data/FVC2000_raw --out Outputs_preprocess_check/db1
```

The procedure is deterministic and gave bit-identical output across the tested library and Python versions (see the requirement files). Each output directory's `preprocess_manifest.json` records the environment, the parameters and, per image, the input md5, the output md5 and the ridge fraction. Compare it with the released `Results/Data/FVC2000_binary_regen/preprocess_manifest.json` (and the DB3 one): equal output md5s for all 800 images mean your binary maps are the ones the paper used; the pipeline then reproduces the reported results exactly in the tested software configuration (R 4.5.3, the package versions recorded in each run manifest's session information, the seeds fixed in the code) and, across other R, BLAS or package versions, up to floating-point differences that the manifests let you measure.

The sensitivity experiment additionally needs the six enhancement variants (`preprocess.py --param KEY=VALUE`, overrides listed in `Code/sensitivity/score_perturbed.R`) under `Data/FVC2000_binary_variant/<condition>/`, and the four calibrated degradations of the raw DB1 images (`Code/sensitivity/calibrate_degradations.py`, `degrade_raw.py`, then `preprocess.py` on each `Data/FVC2000_raw_degraded/<condition>/` into `Data/FVC2000_binary_degraded/<condition>/`); the README's step 7 gives the order.

The loader reads every `*_bin.png` of a binary directory, infers the finger identity from the file-name prefix before the first underscore, and inverts an image whose white fraction exceeds one half (none of the canonical images does).
