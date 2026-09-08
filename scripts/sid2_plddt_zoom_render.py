#!/usr/bin/env python3
"""T96 environment zoom, with the scaffold coloured by pLDDT.

    cd <repo root> && python3 <this file>

Writes only into claude_science/plots/; nothing elsewhere in the
repository is modified.

WHAT CHANGES FROM scripts/sid2_zoom_render.py
Only the scaffold colour. The released zoom greys the ribbon deliberately --
it is scaffold, and grey lets the side chains carry the panel. That works, but
it also means the two halves of Figure 4C speak different languages: the
overview is about the model and the zoom is about the residues, and neither
tells the reader how much of either to trust.

Here the scaffold keeps the same role but carries pLDDT, at reduced opacity so
the sticks still read as the foreground. The reader can then see that T96's
immediate environment is confidently modelled (pLDDT 79 at T96) while the loop
that leaves the frame is not -- which is the honest version of the proximity
argument, and matches the overview panel beside it.

Camera, scaffold radius, stick geometry, label offsets and the Ca-Ca distance
annotations are all imported or copied unchanged from the released script, so
the panel is recognisably the same figure.

OPTIONAL
SHOW_BASIC draws K93 and K132, the two basic residues within 12 A of T96. They
are the reason the T96K charge argument is worth making -- T96's pocket is the
87th percentile of local net charge in an otherwise acidic ectodomain -- but
they are not part of the published uptake-critical set, so the default is off.
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

SHOW_BASIC = False
BASIC = {93: "K93", 132: "K132"}
## spread along the T96->target line, overriding zr.LABEL_FRAC
LABEL_FRAC = {34: 0.56, 32: 0.30, 168: 0.72}
COL_BASIC = "#2F6FB5"
COL_QT13 = "#7B5EA7"
SCAFFOLD_ALPHA = 0.60


def main():
    res, ids, N, CA, C, O = zr.load("A")
    sse = rr.dssp_lite(N, CA, C, O)
    per = pd.read_csv(PERRES, sep="\t").set_index("resid")
    plddt = per.loc[ids, "plddt"].to_numpy(float)

    p96 = CA[ids == zr.FOCAL][0]
    ## same scaffold extent as the released zoom
    keep = np.where(np.linalg.norm(CA - p96, axis=1) <= zr.ZOOM_RADIUS - 1)[0]
    quads, norms, ridx = zr.ribbon_quads(CA[keep], sse[keep], lw_scale=1.0)
    pl_keep = plddt[keep]

    ## Opacity carries confidence as well as hue: a very-low-pLDDT loop is
    ## drawn faintly, so the eye is not pulled to the part of the model that
    ## means least. A single flat alpha let the disordered loop leaving the
    ## frame dominate the panel purely because orange is a loud colour.
    def alpha_for(v):
        if v >= 90: return 0.78
        if v >= 70: return 0.66
        if v >= 50: return 0.42
        return 0.26

    cols = [tuple(rr.shade(rr.plddt_col(pl_keep[r]), n)) + (alpha_for(pl_keep[r]),)
            for r, n in zip(ridx, norms)]
    print("  scaffold: %d residues, pLDDT %.0f-%.0f"
          % (len(keep), pl_keep.min(), pl_keep.max()))

    fig = plt.figure(figsize=(5.6, 5.0), dpi=600)
    ax = fig.add_subplot(111, projection="3d")
    ## no alpha= kwarg: it would override the per-quad alpha above
    ax.add_collection3d(Poly3DCollection(quads, facecolors=cols,
                                         edgecolors="none", shade=False))

    ## zr.FUNC is {32: H32, 34: D34, 168: H168}. Only 32 and 168 are from the
    ## uptake-critical histidine set of McEwan et al. 2012; D34 is the qt13
    ## loss-of-function allele and takes the separate colour used in the
    ## cartoon and histogram panels. H175, the third histidine, is 37.8 A away
    ## and so falls outside this zoom's 19 A scaffold radius.
    show = {zr.FOCAL: ("T96", zr.COL_T96)}
    show.update({p: (l, COL_QT13 if p == 34 else zr.COL_FUNC)
                 for p, l in zr.FUNC.items()})
    if SHOW_BASIC:
        show.update({p: (l, COL_BASIC) for p, l in BASIC.items()})

    r96 = res[int(np.where(ids == zr.FOCAL)[0][0])]
    print("  Ca-Ca and nearest heavy-atom distances from T96:")
    for pos, (lab, col) in show.items():
        w = np.where(ids == pos)[0]
        if not len(w):
            print(f"    WARNING residue {pos} not in the model, skipped")
            continue
        r = res[int(w[0])]
        zr.sticks(r, ax)
        p = r["CA"].coord
        ax.scatter(*p, s=26, color=col, edgecolors="white", linewidths=0.7,
                   depthshade=False, zorder=13)
        off = zr.LABEL_OFF.get(pos, (0.0, 0.0, 2.6))
        ax.text(p[0] + off[0], p[1] + off[1], p[2] + off[2], lab, color=col,
                fontsize=10, fontweight="bold",
                ha=zr.LABEL_HA.get(pos, "center"),
                va=zr.LABEL_VA.get(pos, "bottom"), zorder=16,
                path_effects=[pe.withStroke(linewidth=2.6,
                                            foreground="white")])
        if pos == zr.FOCAL:
            continue
        d_ca = float(np.linalg.norm(p - p96))
        d_min = min(float(np.linalg.norm(a.coord - b.coord))
                    for a in r96 for b in r
                    if a.element != "H" and b.element != "H")
        print(f"    {lab:5s} Ca-Ca {d_ca:5.1f} A   nearest atom {d_min:5.1f} A")
        if pos in BASIC:
            continue          # no distance annotation for the optional pair
        ax.plot(*np.array([p96, p]).T, color=zr.COL_FUNC, linewidth=1.0,
                linestyle=(0, (3, 2.5)), zorder=11,
                path_effects=[pe.withStroke(linewidth=2.4,
                                            foreground="white")])
        ## the three dashed lines converge on T96, so the numbers are spread
        ## further along their own line than in the released panel, where two
        ## of them sat within a few points of each other
        f = LABEL_FRAC.get(pos, zr.LABEL_FRAC.get(pos, 0.5))
        mid = p96 + (p - p96) * f
        ax.text(mid[0], mid[1], mid[2], f"{d_ca:.1f} \u00c5",
                color=zr.COL_FUNC, fontsize=8.6, ha="center", va="center",
                zorder=17,
                path_effects=[pe.withStroke(linewidth=2.8,
                                            foreground="white")])

    pts = np.array([q for quad in quads for q in quad])
    zr.frame_axes(ax, pts, pad=1.0, elev=zr.ZOOM_ELEV, azim=zr.ZOOM_AZIM)

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / ("sid2_zoom_plddt_basic.png" if SHOW_BASIC
                  else "sid2_zoom_plddt.png")
    fig.subplots_adjust(0, 0, 1, 1)
    fig.savefig(dest, dpi=600, bbox_inches="tight", pad_inches=0.02,
                facecolor="white")
    plt.close(fig)
    rr.autocrop(dest)
    print(f"wrote {dest}")


if __name__ == "__main__":
    main()
