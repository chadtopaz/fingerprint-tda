# Data

This directory is intentionally **empty** in the repository. The study uses a third-party benchmark that we do not redistribute; you must obtain it and prepare it locally.

## Dataset

**FVC2000 DB1**, full set ("DB1_A"): 100 fingers x 8 impressions = **800 images**, acquired at 500 dpi with a low-cost optical sensor (300 x 300 px). Cite as:

> Maio, D., Maltoni, D., Cappelli, R., Wayman, J. L., and Jain, A. K. (2002). FVC2000: Fingerprint verification competition. *IEEE Transactions on Pattern Analysis and Machine Intelligence*, 24(3), 402-412.

### Where to get it

- The **full 100-finger DB1** is distributed on the DVD accompanying the *Handbook of Fingerprint Recognition* (Maltoni, Maio, Jain, and Feng; Springer, 3rd ed., 2022).
- The FVC2000 competition site, <http://bias.csr.unibo.it/fvc2000/>, hosts competition details and a small free subset (**DB1_B**, 10 fingers x 8 = 80 images). DB1_B is *not* the 100-finger set used in the paper.

## Preparing the binary images the code reads

The pipeline does not read raw FVC2000 images directly; it reads **enhanced binary** ridge/valley segmentations.

1. Apply Hong-Wan-Jain enhancement, as implemented in the **Fingerprint-Enhancement-Python** library (Utkarsh Deshmukh, <https://github.com/Utkarsh-Deshmukh/Fingerprint-Enhancement-Python>), with default parameters, to each raw image.
2. Save the results as PNGs in this directory under `FVC2000_binary/`, named:

   ```
   Data/FVC2000_binary/<finger>_<impression>_bin.png
   ```

   for example `1_1_bin.png`, `1_2_bin.png`, ..., `100_8_bin.png` (800 files).

The loader (`read_binary_fingerprint`) reads every `*.png` here, infers the finger identity from the filename prefix before the first underscore, and treats ridge pixels as the minority class (inverting the labeling automatically if more than half the pixels are flagged as ridge).

Note: exact bit-for-bit reproduction of the binary images depends on the enhancement library version. Given a fixed set of binary images, the downstream topological pipeline is exactly reproducible (all seeds and iteration counts are fixed).
