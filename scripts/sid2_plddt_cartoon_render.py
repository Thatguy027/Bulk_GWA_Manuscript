#!/usr/bin/env python3
"""Membrane-oriented SID-2 ectodomain cartoon, coloured by pLDDT.

    cd <repo root> && python3 <this file>

Writes only into claude_science/plots/; nothing elsewhere in the
repository is modified.

WHY THIS EXISTS
scripts/sid2_ribbon_render.py can already colour a cartoon by pLDDT, but it
orients the model by PCA, which is reproducible and anatomically meaningless.
scripts/sid2_zoom_render.py keeps the membrane frame -- z = 0 is the bilayer
centre, +z is the intestinal lumen -- but reads
supplemental_data/structure/sid2_membrane_oriented.pdb, whose B-factor column
holds assigned partial charges rather than pLDDT. So the released panel is a
correctly oriented cartoon that cannot be coloured by confidence.

This joins the two: geometry and camera from the zoom script, per-residue
pLDDT from supplemental_data/structure/sid2_per_residue.tsv keyed on residue
number, and the same Kabsch-Sander secondary structure, so strands still carry
their arrowheads. Every geometry helper is imported from the existing scripts
rather than reimplemented; only the colour source and the marker set differ.

MARKERS
T96 plus all four published uptake-critical residues. The released overview
marks T96 alone and the zoom draws three of the four; H175 at 37.8 A is the one
that never appears, and it is the reason the panel reads as stronger evidence
than the binomial supports.
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
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

sys.path.insert(0, "scripts")
import sid2_ribbon_render as rr
import sid2_zoom_render as zr

warnings.filterwarnings("ignore")

OUT = Path("plots/assets")
PERRES = Path("supplemental_data/structure/sid2_per_residue.tsv")

COL_T96 = "#F34C00"
COL_FUNC = "#16324A"
MEMBRANE = "#DCE3E8"
FOCAL = 96
ZOOM_RADIUS = 20.0

## Two distinct classes of published residue, which the released panel merges.
##
## The uptake-critical set from McEwan, Weisman & Hunter 2012 (Mol Cell 47:746,
## doi:10.1016/j.molcel.2012.07.014) is exactly three extracellular HISTIDINES:
## His32, His168, His175. They were targeted because imidazole is protonated
## only in the acidic conditions SID-2 requires; His->Ala and His->Glu each
## reduced dsRNA transport. That paper never mentions residue 34.
##
## D34 is a different kind of evidence -- the qt13 loss-of-function allele --
## so it is drawn in a separate colour rather than folded into the same set.
##
## The model's ectodomain contains exactly three histidines, at 32, 168 and
## 175, which independently confirms the residue numbering.
PUB = {32: "H32", 168: "H168", 175: "H175"}
QT13 = {34: "D34"}
COL_QT13 = "#7B5EA7"

## Labels are placed OUTSIDE the projected silhouette rather than at a small
## offset from the residue. At azimuth 108 the membrane-frame x axis maps
## almost exactly onto screen-horizontal, so pushing a label to just beyond
## min(x) or max(x) puts it clear of the ribbon whatever the residue's depth;
## a leader line carries the eye back. A hand-tuned offset per residue was the
## previous approach and put three labels on top of the beta-sandwich.
##
## Which side each residue takes, and the vertical nudge applied when two on
## the same side would otherwise sit within ~7 A of each other in z.
LABEL_SIDE = {96: "right", 34: "left", 32: "left", 168: "right", 175: "left"}
LABEL_DZ = {96: 3.0, 34: 1.5, 32: -3.0, 168: -1.0, 175: 0.0}

ELEV, AZIM = zr.OVERVIEW_ELEV, zr.OVERVIEW_AZIM


def main():
    res, ids, N, CA, C, O = zr.load("A")
    sse = rr.dssp_lite(N, CA, C, O)

    per = pd.read_csv(PERRES, sep="\t").set_index("resid")
    plddt = per.loc[ids, "plddt"].to_numpy(float)
    assert np.isfinite(plddt).all(), "pLDDT missing for a modelled residue"

    quads, norms, ridx = zr.ribbon_quads(CA, sse)
    cols = [rr.shade(rr.plddt_col(plddt[r]), n) for r, n in zip(ridx, norms)]

    fig = plt.figure(figsize=(3.9, 6.2), dpi=600)
    ax = fig.add_subplot(111, projection="3d")
    ax.add_collection3d(Poly3DCollection(quads, facecolors=cols,
                                         edgecolors="none", shade=False))

    ## the bilayer as one edge-on quad in the x-z plane, as in the released
    ## overview: two horizontal leaflets project as detached plates
    x0, x1 = CA[:, 0].min() - 5, CA[:, 0].max() + 5
    ym = CA[:, 1].mean()
    ## the slab is clipped at z = -9 rather than the full -13: the ectodomain
    ## starts at z = 11, so the lower leaflet is dead page in a panel that is
    ## already tall and narrow
    ax.add_collection3d(Poly3DCollection(
        [[(x0, ym, -9.0), (x1, ym, -9.0), (x1, ym, 13.0), (x0, ym, 13.0)]],
        facecolors=MEMBRANE, edgecolors="none", alpha=0.85, zorder=0))

    p96 = CA[ids == FOCAL][0]
    th = np.linspace(0, 2 * np.pi, 120)
    ax.plot(p96[0] + ZOOM_RADIUS * np.cos(th),
            p96[1] + ZOOM_RADIUS * np.sin(th) * 0.35,
            p96[2] + ZOOM_RADIUS * np.sin(th), color=COL_T96, linewidth=0.8,
            linestyle=(0, (3, 3)), alpha=0.75, zorder=14)

    ## markers: the focal residue, then the published set
    x_left, x_right = CA[:, 0].min() - 3.5, CA[:, 0].max() + 3.5
    marks = [(FOCAL, "T96")] + sorted(PUB.items()) + sorted(QT13.items())
    for pos, lab in marks:
        w = np.where(ids == pos)[0]
        if not len(w):
            print(f"  WARNING residue {pos} absent from the model, skipped")
            continue
        p = CA[w[0]]
        focal = pos == FOCAL
        col = COL_T96 if focal else (COL_QT13 if pos in QT13 else COL_FUNC)
        ax.scatter(*p, s=90 if focal else 55, color=col,
                   edgecolors="white", linewidths=1.2 if focal else 0.9,
                   depthshade=False, zorder=15)
        left = LABEL_SIDE[pos] == "left"
        lx = x_left if left else x_right
        lz = p[2] + LABEL_DZ[pos]
        ax.plot([p[0], lx], [p[1], ym], [p[2], lz], color=col,
                linewidth=0.7, solid_capstyle="round", zorder=15.5,
                alpha=0.85,
                path_effects=[pe.withStroke(linewidth=1.8,
                                            foreground="white")])
        ax.text(lx, ym, lz, lab, color=col,
                fontsize=7.4 if focal else 6.8, fontweight="bold",
                ha="right" if left else "left", va="center", zorder=16,
                path_effects=[pe.withStroke(linewidth=2.0,
                                            foreground="white")])

    allp = np.vstack([np.array([q for quad in quads for q in quad]),
                      np.array([[x0, ym, -9.0], [x1, ym, 13.0]])])
    zr.frame_axes(ax, allp, pad=1.0, elev=ELEV, azim=AZIM)

    for txt, zc, va in (("lumen", CA[:, 2].max(), "bottom"),
                        ("membrane", 2.0, "center")):
        ax.text(CA[:, 0].mean(), ym, zc, txt, color="grey", fontsize=6.8,
                ha="center", va=va, zorder=17,
                path_effects=[pe.withStroke(linewidth=2.0,
                                            foreground="white")])
    ## name the dashed ring, at its top-left where the page is empty, so it
    ## reads as the radius the companion histogram uses
    ax.text(p96[0] - ZOOM_RADIUS * 0.80, ym,
            p96[2] + ZOOM_RADIUS * 0.66, "20 \u00c5", color=COL_T96,
            fontsize=6.8, ha="right", va="center", zorder=17,
            path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / "sid2_cartoon_plddt.png"
    fig.subplots_adjust(0, 0, 1, 1)
    fig.savefig(dest, dpi=600, bbox_inches="tight", pad_inches=0.02,
                facecolor="white")
    plt.close(fig)
    rr.autocrop(dest)

    print(f"wrote {dest}")
    print("pLDDT bands (rr.PLDDT_BANDS):", rr.PLDDT_BANDS)
    print("T96 pLDDT %.1f | ectodomain n=%d | pLDDT>=70 %.3f"
          % (per.loc[FOCAL, "plddt"], len(ids), float((plddt >= 70).mean())))
    print("sse counts:", {s: int((sse == s).sum()) for s in set(sse)})


if __name__ == "__main__":
    main()
