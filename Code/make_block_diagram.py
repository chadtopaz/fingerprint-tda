# Block diagram of the end-to-end pipeline (Fig. 1, top panel). Run from any directory; writes fig_block_diagram.pdf
# (345 x 150 pt); copy to Pattern Recognition/.
# Audit corrective pass, step 4 (1 Nov 2026; memo item C3): the clamp-and-round step that precedes the persistence
# computation is shown as its own box; the geometry baselines branch off after downsampling (the unclamped DT L^2 and
# its aligned variant) and after clamping and rounding (the clamped DT L^2); the fusion box names the robust-z
# standardization that precedes the equal-weight combination.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams.update({"font.family": "serif", "font.size": 7.5, "pdf.fonttype": 42})
W, H = 345, 150
fig = plt.figure(figsize=(W/72, H/72)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
def box(x, y, w, h, text, fc="#f2f2f2", fs=7.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=3", fc=fc, ec="black", lw=0.7))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs, linespacing=1.1)
def line(pts, arrow=True, ls="-"):
    for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
        last = (x1, y1) == pts[-1]
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>" if (arrow and last) else "-", mutation_scale=6, lw=0.7, color="black", ls=ls, shrinkA=0, shrinkB=0))
# row 1 (y 106..146): six boxes
y1, h1 = 106, 40
R1 = [(2, 40, "raw\nimage", 7.5), (47, 52, "enhance,\nbinarize", 7.5), (104, 62, "ridge / valley\ndistance\ntransforms", 7.0),
      (171, 52, "downsample\n150 x 150", 7.2), (228, 40, "clamp\n[0, 20],\nround 0.1", 7.0), (273, 66, "$H_1$ persistence\nglobal + local\npatches", 7.0)]
for x, w, t, fs in R1:
    box(x, y1, w, h1, t, fs=fs)
for (xa, wa, _, _), (xb, _, _, _) in zip(R1[:-1], R1[1:]):
    line([(xa + wa, y1 + h1/2), (xb, y1 + h1/2)])
# branches in the gap (y 78..100): reference matcher from the raw image; geometry baselines from the downsampled
# fields (unclamped DT L^2, aligned) and from the clamped, rounded fields (clamped DT L^2)
ax.text(4, 88, "reference matcher\n(NBIS)", ha="left", va="center", fontsize=6.8, style="italic", linespacing=1.1)
line([(20, y1), (20, 97)], ls=(0, (2, 1.5)))
xg1, xg2 = 197, 248                      # centres of the downsample and clamp/round boxes
line([(xg1, y1), (xg1, 97)], ls=(0, (2, 1.5)))
line([(xg2, y1), (xg2, 97)], ls=(0, (2, 1.5)))
ax.text(xg1, 87, "geometry baselines\n(unclamped, aligned)", ha="center", va="center", fontsize=6.4, style="italic", linespacing=1.1)
ax.text(xg2 + 6, 87, "geometry\nbaseline\n(clamped)", ha="center", va="center", fontsize=6.4, style="italic", linespacing=1.05)
# row 2 (y 30..74), ends at x = 330 to leave the right channel
y2, h2 = 30, 44
box(14,  y2, 82, h2, "Betti curves,\npersistence\nimages", fs=7.2)
box(103, y2, 118, h2, "scorers: Betti $L^1$,\nPI $L^2$, TopoLR,\nLPOT", fs=7.2)
box(228, y2, 52, h2, "robust $z$,\nequal-weight\nfusion", fs=7.0)
box(287, y2, 43, h2, "score,\nthresholds", fs=7.0)
for a, b in [(96, 103), (221, 228), (280, 287)]:
    line([(a, y2 + h2/2), (b, y2 + h2/2)])
# connector: persistence box -> right channel -> bottom channel -> first row-2 box
line([(339, y1), (339, 16), (8, 16), (8, y2 + h2/2), (14, y2 + h2/2)])
fig.savefig("fig_block_diagram.pdf"); print("ok")
