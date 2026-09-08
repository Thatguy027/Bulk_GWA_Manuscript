#!/usr/bin/env python3
"""SID-2 ectodomain rendered by local electrostatic character, not confidence.

    cd <repo root> && python3 <this file>
        -> claude_science/plots/sid2_cartoon_charge.png   (whole ectodomain)
        -> claude_science/plots/sid2_zoom_charge.png      (the T96 pocket)

Writes only into claude_science/plots/; nothing elsewhere in the
repository is modified.

WHY CHARGE AND NOT pLDDT
McEwan, Weisman & Hunter 2012 (Mol Cell 47:746) show that SID-2 dsRNA uptake
needs an acidic extracellular environment and depends on three extracellular
histidines, which are protonated only at that pH -- and that a triple
histidine-to-arginine mutant, i.e. permanent rather than pH-dependent positive
charge, internalises MORE dsRNA than wild type. T96K adds a permanent positive
charge to the same lumenal domain. The panel therefore has to show charge.

WHAT IS COLOURED
Each residue is coloured by the NET CHARGE OF ITS 12 A NEIGHBOURHOOD, not by
its own identity -- a per-residue class colouring gives stripes, and the claim
is about a pocket, not a residue. Charges are Henderson-Hasselbalch side-chain
charges at pH 4.4, the gut-lumen pH, using the same pKa set as the manuscript's
electrostatics supplement (Asp 3.9, Glu 4.25, His 6.0, Lys 10.5, Arg 12.5,
Cys 8.3, Tyr 10.1). Diverging scale, white at neutral, saturating at +/-2 e.

Geometry, camera, ribbon construction, stick drawing and framing are all
imported from scripts/sid2_ribbon_render.py and scripts/sid2_zoom_render.py,
so these panels are the same figures in a different colour space.
"""

import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm, colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

sys.path.insert(0, "scripts")
import sid2_ribbon_render as rr
import sid2_zoom_render as zr

warnings.filterwarnings("ignore")

OUT = Path("plots/assets")
CHG = Path("supplemental_data/structure/sid2_local_charge.tsv")

FOCAL = 96
COL_T96 = "#F34C00"
COL_BASIC = "#0B4F9E"          # the two lysines that make the pocket
COL_HIS = "#1B7F79"            # the pH-titrating uptake-critical set
MEMBRANE = "#DCE3E8"
QLIM = 2.0                     # colour scale saturates here, in e

BASIC = {93: "K93", 132: "K132"}
HIS = {32: "H32", 168: "H168", 175: "H175"}

CMAP = matplotlib.colormaps["RdBu"]   # red negative, blue positive, white at 0
NORM = colors.Normalize(vmin=-QLIM, vmax=QLIM)


def qcol(q):
    return CMAP(NORM(q))[:3]


def draw_marker(ax, p, lab, col, focal=False, dz=3.0, side=None, xr=None):
    ax.scatter(*p, s=90 if focal else 55, color=col, edgecolors="white",
               linewidths=1.2 if focal else 0.9, depthshade=False, zorder=15)
    if side is None:
        return
    lx = xr[0] if side == "left" else xr[1]
    ax.plot([p[0], lx], [p[1], p[1]], [p[2], p[2] + dz], color=col,
            linewidth=0.7, alpha=0.85, solid_capstyle="round", zorder=15.5,
            path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])
    ax.text(lx, p[1], p[2] + dz, lab, color=col,
            fontsize=7.4 if focal else 6.8, fontweight="bold",
            ha="right" if side == "left" else "left", va="center", zorder=16,
            path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])


def overview(ids, CA, sse, q):
    quads, norms, ridx = zr.ribbon_quads(CA, sse)
    cols = [rr.shade(qcol(q[r]), n) for r, n in zip(ridx, norms)]

    fig = plt.figure(figsize=(3.9, 6.2), dpi=600)
    ax = fig.add_subplot(111, projection="3d")
    ax.add_collection3d(Poly3DCollection(quads, facecolors=cols,
                                         edgecolors="none", shade=False))

    x0, x1 = CA[:, 0].min() - 5, CA[:, 0].max() + 5
    ym = CA[:, 1].mean()
    ax.add_collection3d(Poly3DCollection(
        [[(x0, ym, -9.0), (x1, ym, -9.0), (x1, ym, 13.0), (x0, ym, 13.0)]],
        facecolors=MEMBRANE, edgecolors="none", alpha=0.85, zorder=0))

    ## At azimuth 108 the membrane-frame x axis maps onto screen-horizontal
    ## REVERSED, so side="left" (x minimum) puts a label on the screen right.
    ## The pocket trio -- T96, K93, K132 -- sit within ~7 A of one another, so
    ## they need a large vertical spread or their labels overlap; they take one
    ## side, the histidines the other.
    xr = (CA[:, 0].min() - 3.5, CA[:, 0].max() + 3.5)
    plan = [(FOCAL, "T96", COL_T96, "left", 9.0, True)]
    plan += [(p, l, COL_BASIC, "left", d, False)
             for (p, l), d in zip(sorted(BASIC.items()), (1.0, -7.0))]
    plan += [(p, l, COL_HIS, "right", d, False)
             for (p, l), d in zip(sorted(HIS.items()), (3.0, -1.5, -5.0))]
    for pos, lab, col, side, dz, foc in plan:
        w = np.where(ids == pos)[0]
        if not len(w):
            print(f"  WARNING residue {pos} absent, skipped")
            continue
        draw_marker(ax, CA[w[0]], lab, col, focal=foc, dz=dz, side=side, xr=xr)

    allp = np.vstack([np.array([p for quad in quads for p in quad]),
                      np.array([[x0, ym, -9.0], [x1, ym, 13.0]])])
    zr.frame_axes(ax, allp, pad=1.0, elev=zr.OVERVIEW_ELEV, azim=zr.OVERVIEW_AZIM)
    for txt, zc, va in (("lumen", CA[:, 2].max(), "bottom"),
                        ("membrane", 2.0, "center")):
        ax.text(CA[:, 0].mean(), ym, zc, txt, color="grey", fontsize=6.8,
                ha="center", va=va, zorder=17,
                path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])

    dest = OUT / "sid2_overview_charge.png"
    fig.subplots_adjust(0, 0, 1, 1)
    fig.savefig(dest, dpi=600, bbox_inches="tight", pad_inches=0.02,
                facecolor="white")
    plt.close(fig)
    rr.autocrop(dest)
    print(f"wrote {dest}")


def zoom(res, ids, CA, sse, q):
    p96 = CA[ids == FOCAL][0]
    keep = np.where(np.linalg.norm(CA - p96, axis=1) <= zr.ZOOM_RADIUS - 1)[0]
    quads, norms, ridx = zr.ribbon_quads(CA[keep], sse[keep], lw_scale=1.0)
    qk = q[keep]
    cols = [tuple(rr.shade(qcol(qk[r]), n)) + (0.72,) for r, n in zip(ridx, norms)]

    fig = plt.figure(figsize=(5.6, 5.0), dpi=600)
    ax = fig.add_subplot(111, projection="3d")
    ax.add_collection3d(Poly3DCollection(quads, facecolors=cols,
                                         edgecolors="none", shade=False))

    ## the pocket: T96 and the two lysines that make it basic
    print("  Ca-Ca and nearest heavy-atom distances from T96:")
    r96 = res[int(np.where(ids == FOCAL)[0][0])]
    for pos, lab in [(FOCAL, "T96")] + sorted(BASIC.items()):
        r = res[int(np.where(ids == pos)[0][0])]
        zr.sticks(r, ax)
        p = r["CA"].coord
        col = COL_T96 if pos == FOCAL else COL_BASIC
        ax.scatter(*p, s=30, color=col, edgecolors="white", linewidths=0.7,
                   depthshade=False, zorder=13)
        ## K93 and K132 are 6.6 and 6.8 A from T96 and close to each other, so
        ## their labels are pushed to opposite sides rather than both below
        off = {FOCAL: (0.0, 0.0, 3.4), 93: (-6.5, 0.0, -1.5),
               132: (6.5, 0.0, -3.0)}[pos]
        ax.text(p[0] + off[0], p[1] + off[1], p[2] + off[2], lab, color=col,
                fontsize=10, fontweight="bold",
                ha={FOCAL: "center", 93: "right", 132: "left"}[pos],
                va="bottom" if pos == FOCAL else "center", zorder=16,
                path_effects=[pe.withStroke(linewidth=2.6, foreground="white")])
        if pos == FOCAL:
            continue
        d_ca = float(np.linalg.norm(p - p96))
        d_min = min(float(np.linalg.norm(a.coord - b.coord))
                    for a in r96 for b in r
                    if a.element != "H" and b.element != "H")
        print(f"    {lab:5s} Ca-Ca {d_ca:5.1f} A   nearest atom {d_min:5.1f} A")
        ax.plot(*np.array([p96, p]).T, color=COL_BASIC, linewidth=1.0,
                linestyle=(0, (3, 2.5)), zorder=11,
                path_effects=[pe.withStroke(linewidth=2.4, foreground="white")])
        mid = p96 + (p - p96) * (0.42 if pos == 93 else 0.70)
        ax.text(mid[0], mid[1], mid[2], f"{d_ca:.1f} \u00c5", color=COL_BASIC,
                fontsize=8.6, ha="center", va="center", zorder=17,
                path_effects=[pe.withStroke(linewidth=2.8, foreground="white")])

    pts = np.array([p for quad in quads for p in quad])
    zr.frame_axes(ax, pts, pad=1.0, elev=zr.ZOOM_ELEV, azim=zr.ZOOM_AZIM)
    dest = OUT / "sid2_zoom_charge.png"
    fig.subplots_adjust(0, 0, 1, 1)
    fig.savefig(dest, dpi=600, bbox_inches="tight", pad_inches=0.02,
                facecolor="white")
    plt.close(fig)
    rr.autocrop(dest)
    print(f"wrote {dest}")


def main():
    res, ids, N, CA, C, O = zr.load("A")
    sse = rr.dssp_lite(N, CA, C, O)
    chg = pd.read_csv(CHG, sep="\t").set_index("resid")
    q = chg.loc[ids, "q_local_pH44"].to_numpy(float)
    assert np.isfinite(q).all(), "local charge missing for a modelled residue"
    print("local 12 A charge at pH 4.4: %.2f to %.2f e, T96 %+.2f"
          % (q.min(), q.max(), float(chg.loc[FOCAL, "q_local_pH44"])))
    overview(ids, CA, sse, q)
    zoom(res, ids, CA, sse, q)


if __name__ == "__main__":
    main()
