#!/usr/bin/env python3
"""Membrane-oriented overview and zoom of the SID-2 ectodomain, for Figure 4.

    python3 scripts/sid2_zoom_render.py

Writes to plots/assets/:
    sid2_overview_oriented.png   whole ectodomain, lumenal face up, with the
                                 zoom region ringed and the bilayer marked
    sid2_zoom_t96.png            zoom on T96 with side chains as sticks and
                                 dashed Ca-Ca distances to the published
                                 uptake-critical residues

WHY A DIFFERENT INPUT FILE
The AlphaFold output (data/structure/AF_sid2_wt_dimer/sid2wt.pdb) is in an
arbitrary frame, so scripts/sid2_ribbon_render.py orients it by PCA -- which
gives a reproducible view but not a meaningful one. This script instead reads
    data/structure_modeling/claude_docking/stage4_electrostatics/sid2_wt_elec.pdb
which is the same model already rotated onto the membrane normal by the
earlier stage3 work: z = 0 is the bilayer centre, +z is the intestinal lumen.
Keeping that frame and putting +z up the page means the picture is anatomically
right -- the lumenal face points into the gut, the way the protein sits.

That file covers residues 21-188, the UniProt ectodomain boundary, rather than
the 21-193 DeepTMHMM boundary used elsewhere. The difference is five residues
of linker and does not touch anything annotated here.

The B-factor column of that file holds assigned partial charges, NOT pLDDT, so
nothing here is coloured by confidence. Secondary structure is recomputed from
the backbone with the same Kabsch & Sander routine used by the ribbon script.

DISTANCES are Ca-Ca, matching the numbers quoted in the figure and the text.
Nearest heavy-atom distances are also printed, because they are what a reader
looking at sticks will try to estimate off the page, and they are shorter.
"""

import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBParser
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

sys.path.insert(0, "scripts")
import sid2_ribbon_render as rr          # geometry + Kabsch-Sander helpers

warnings.filterwarnings("ignore")

PDB = Path("data/structure_modeling/claude_docking/stage4_electrostatics/"
           "sid2_wt_elec.pdb")
OUT = Path("plots/assets")
OUT.mkdir(parents=True, exist_ok=True)

COL_T96  = "#F34C00"
COL_FUNC = "#16324A"
COL_STICK = {"C": "#4A5A66", "N": "#2F6FB5", "O": "#C0392B", "S": "#C9A227"}
MEMBRANE = "#DCE3E8"
ZOOM_RADIUS = 20.0          # A around T96, the radius used in the p-value
FOCAL = 96

## Camera angles, one place. Elevation is tilt out of the membrane plane and
## azimuth is rotation about the membrane normal, so the lumen stays up as long
## as the elevation is low. Pick replacements off the contact sheets built by
## scripts/sid2_orientation_sheet.py.
##
## OVERVIEW_AZIM is the previous view (-72) turned through 180 degrees, so the
## far face of the ectodomain now points at the reader.
OVERVIEW_ELEV, OVERVIEW_AZIM = 8, 108
## the zoom lifted and swung right from its previous (6, -104)
ZOOM_ELEV, ZOOM_AZIM = 18, -89
## label placement along the T96->target line, per residue: all three targets
## project below T96 in this view, so labels at a common midpoint stacked
FUNC = {32: "H32", 34: "D34", 168: "H168"}
LABEL_FRAC = {34: 0.40, 32: 0.58, 168: 0.76}
## per-residue label offsets in the membrane frame. A single +z offset for all
## four put H32's label on top of a dashed line at this camera angle.
LABEL_OFF = {96: (0.0, 0.0, 2.6), 34: (2.4, 0.0, 1.8),
             32: (-3.4, 0.0, 1.2), 168: (0.0, 0.0, -3.0)}
LABEL_HA = {96: "center", 34: "left", 32: "right", 168: "center"}
LABEL_VA = {96: "bottom", 34: "bottom", 32: "bottom", 168: "top"}


def load(cid="A"):
    st = PDBParser(QUIET=True).get_structure("sid2", PDB)
    ch = [c for c in list(st)[0] if c.id == cid][0]
    res = [r for r in ch if r.id[0] == " "
           and all(a in r for a in ("N", "CA", "C", "O"))]
    ids = np.array([r.id[1] for r in res])
    get = lambda nm: np.array([r[nm].coord for r in res], dtype=float)
    return res, ids, get("N"), get("CA"), get("C"), get("O")


def ribbon_quads(ca, sse, lw_scale=1.0):
    b = rr.frames(ca)
    w = rr.widths(sse) * lw_scale
    pts, side, wid, idx = rr.resample(ca, b, w)
    quads, norms, ridx = [], [], []
    for k in range(len(pts) - 1):
        p0, p1 = pts[k], pts[k + 1]
        s0, s1 = side[k] * wid[k], side[k + 1] * wid[k + 1]
        quads.append([p0 - s0, p0 + s0, p1 + s1, p1 - s1])
        n = np.cross(p1 - p0, side[k])
        nl = np.linalg.norm(n)
        norms.append(n / nl if nl > 1e-9 else np.array([0.0, 0.0, 1.0]))
        ridx.append(int(round(min(idx[k], len(sse) - 1))))
    return quads, norms, ridx


def sticks(res, ax, scale=1.0):
    """Side-chain bonds, inferred from interatomic distance.

    No connectivity table is needed: covalent bonds between heavy atoms are
    1.2-1.9 A and nothing else in a residue is that close.
    """
    atoms = [a for a in res if a.element != "H"]
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            p, q = atoms[i].coord, atoms[j].coord
            if 1.0 < np.linalg.norm(p - q) < 1.95:
                for a, xyz in ((atoms[i], p), (atoms[j], q)):
                    mid = (p + q) / 2
                    seg = np.array([xyz, mid])
                    ax.plot(seg[:, 0], seg[:, 1], seg[:, 2],
                            color=COL_STICK.get(a.element, "#4A5A66"),
                            linewidth=2.6 * scale, solid_capstyle="round",
                            zorder=12)


def frame_axes(ax, pts, pad=2.0, elev=8, azim=-72):
    lo, hi = pts.min(0) - pad, pts.max(0) + pad
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    ax.set_axis_off()
    ## keep the membrane normal (z) vertical on the page: a low elevation and
    ## no PCA rotation means +z reads as "up, into the lumen"
    ax.view_init(elev=elev, azim=azim)
    try:
        ax.set_box_aspect(tuple(hi - lo))
    except Exception:
        pass


def overview(res, ids, CA, sse, fname, elev=None, azim=None):
    quads, norms, ridx = ribbon_quads(CA, sse)
    cols = [rr.shade(rr.SSE_COL[sse[r]], n) for r, n in zip(ridx, norms)]

    fig = plt.figure(figsize=(3.6, 6.2), dpi=600)
    ax = fig.add_subplot(111, projection="3d")
    ax.add_collection3d(Poly3DCollection(quads, facecolors=cols,
                                         edgecolors="none", shade=False))

    ## the bilayer, from the stage3 orientation: z = 0 is its centre and the
    ## hydrophobic slab is about 32 A thick
    ## The bilayer as ONE quad in the x-z plane, not two horizontal leaflets.
    ## Viewed near edge-on it fills as a band, which is what a membrane
    ## cross-section should look like; two horizontal quads projected as
    ## detached plates with 26 A of empty page between them.
    x0, x1 = CA[:, 0].min() - 5, CA[:, 0].max() + 5
    ym = CA[:, 1].mean()
    ax.add_collection3d(Poly3DCollection(
        [[(x0, ym, -13.0), (x1, ym, -13.0), (x1, ym, 13.0), (x0, ym, 13.0)]],
        facecolors=MEMBRANE, edgecolors="none", alpha=0.85, zorder=0))

    p96 = CA[ids == FOCAL][0]
    ## ring the zoom region rather than boxing it: a wireframe box on an
    ## invisible-axes 3D plot reads as clutter
    th = np.linspace(0, 2 * np.pi, 120)
    for r in (ZOOM_RADIUS,):
        ax.plot(p96[0] + r * np.cos(th), p96[1] + r * np.sin(th) * 0.35,
                p96[2] + r * np.sin(th), color=COL_T96, linewidth=1.1,
                linestyle=(0, (4, 3)), zorder=14)
    ## the anchor of the whole panel: large enough to read after the rotation,
    ## with a white ring so it separates from the dark ribbon behind it
    ax.scatter(*p96, s=95, color=COL_T96, edgecolors="white", linewidths=1.4,
               depthshade=False, zorder=15)

    allp = np.vstack([np.array([q for quad in quads for q in quad]),
                      np.array([[x0, ym, -13.0], [x1, ym, 13.0]])])
    frame_axes(ax, allp, pad=1.0,
               elev=OVERVIEW_ELEV if elev is None else elev,
               azim=OVERVIEW_AZIM if azim is None else azim)
    for txt, zc, va in (("lumen", CA[:, 2].max(), "bottom"),
                        ("membrane", 0.0, "center")):
        ax.text(CA[:, 0].mean(), ym, zc, txt, color="grey", fontsize=8,
                ha="center", va=va, zorder=16,
                path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])

    fig.subplots_adjust(0, 0, 1, 1)
    fig.savefig(OUT / fname, dpi=600, bbox_inches="tight", pad_inches=0.02,
                facecolor="white")
    plt.close(fig)
    rr.autocrop(OUT / fname)
    print(f"  wrote {OUT / fname}")


def zoom(res, ids, CA, sse, fname, elev=None, azim=None):
    elev = ZOOM_ELEV if elev is None else elev
    azim = ZOOM_AZIM if azim is None else azim
    p96 = CA[ids == FOCAL][0]
    ## scaffold just wide enough to hold the farthest marked residue; the
    ## original +4 A pulled in the disordered linker, which dominated the frame
    near = np.linalg.norm(CA - p96, axis=1) <= ZOOM_RADIUS - 1
    keep = np.where(near)[0]

    quads, norms, ridx = ribbon_quads(CA[keep], sse[keep], lw_scale=1.0)
    cols = [rr.shade("#C7D2DA", n) for n in norms]     # faded, it is scaffold

    fig = plt.figure(figsize=(5.6, 5.0), dpi=600)
    ax = fig.add_subplot(111, projection="3d")
    ax.add_collection3d(Poly3DCollection(quads, facecolors=cols,
                                         edgecolors="none", shade=False,
                                         alpha=0.85))

    show = {FOCAL: ("T96", COL_T96)}
    show.update({p: (l, COL_FUNC) for p, l in FUNC.items()})
    print("  Ca-Ca and nearest heavy-atom distances from T96:")
    for pos, (lab, col) in show.items():
        r = res[int(np.where(ids == pos)[0][0])]
        sticks(r, ax)
        p = r["CA"].coord
        ax.scatter(*p, s=26, color=col, edgecolors="white", linewidths=0.7,
                   depthshade=False, zorder=13)
        off = LABEL_OFF.get(pos, (0.0, 0.0, 2.6))
        ax.text(p[0] + off[0], p[1] + off[1], p[2] + off[2], lab, color=col,
                fontsize=10, fontweight="bold",
                ha=LABEL_HA.get(pos, "center"), va=LABEL_VA.get(pos, "bottom"),
                zorder=16,
                path_effects=[pe.withStroke(linewidth=2.6, foreground="white")])
        if pos == FOCAL:
            continue
        r96 = res[int(np.where(ids == FOCAL)[0][0])]
        d_ca = float(np.linalg.norm(p - p96))
        d_min = min(float(np.linalg.norm(a.coord - b.coord))
                    for a in r96 for b in r
                    if a.element != "H" and b.element != "H")
        print(f"    {lab:5s} Ca-Ca {d_ca:5.1f} A   nearest atom {d_min:5.1f} A")
        ## dashed Ca-Ca line, labelled with the distance quoted in the text
        ax.plot(*np.array([p96, p]).T, color=COL_FUNC, linewidth=1.0,
                linestyle=(0, (3, 2.5)), zorder=11,
                path_effects=[pe.withStroke(linewidth=2.4,
                                            foreground="white")])
        f = LABEL_FRAC.get(pos, 0.5)
        mid = p96 + (p - p96) * f
        ax.text(mid[0], mid[1], mid[2], f"{d_ca:.1f} Å", color=COL_FUNC,
                fontsize=8.6, ha="center", va="center", zorder=17,
                path_effects=[pe.withStroke(linewidth=2.8,
                                            foreground="white")])

    pts = np.array([q for quad in quads for q in quad])
    frame_axes(ax, pts, pad=1.0, elev=elev, azim=azim)
    fig.subplots_adjust(0, 0, 1, 1)
    fig.savefig(OUT / fname, dpi=600, bbox_inches="tight", pad_inches=0.02,
                facecolor="white")
    plt.close(fig)
    rr.autocrop(OUT / fname)
    print(f"  wrote {OUT / fname}")


def main():
    res, ids, N, CA, C, O = load("A")
    sse = rr.dssp_lite(N, CA, C, O)
    print(f"chain A: {len(ids)} residues {ids.min()}-{ids.max()} | "
          f"z {CA[:, 2].min():.1f} to {CA[:, 2].max():.1f} A (+z = lumen)")
    print(f"  strand {np.mean(sse == 'E'):.0%} helix {np.mean(sse == 'H'):.0%}")
    print(f"overview, lumenal face up (elev {OVERVIEW_ELEV}, "
          f"azim {OVERVIEW_AZIM})")
    overview(res, ids, CA, sse, "sid2_overview_oriented.png")
    print(f"zoom on T96 (elev {ZOOM_ELEV}, azim {ZOOM_AZIM})")
    zoom(res, ids, CA, sse, "sid2_zoom_t96.png")
    ## the previous overview angle, kept so the 180-degree turn is checkable
    print("overview at the previous angle, for comparison")
    overview(res, ids, CA, sse, "sid2_overview_oriented_prev.png",
             elev=8, azim=-72)


if __name__ == "__main__":
    main()
