#!/usr/bin/env python3
"""Render the AlphaFold3 SID-2 model for Figure 4.

    python3 scripts/sid2_structure_render.py

Writes to plots/assets/:
    sid2_monomer_plddt.png   full chain A, C-alpha trace coloured by pLDDT
    sid2_ecd_plddt.png       extracellular domain only (residues 21-193)
    sid2_dimer_chains.png    both chains, for reference only -- see CAVEAT

Input: data/structure/AF_sid2_wt_dimer/sid2wt.pdb, the top-ranked of five
AlphaFold3 models. pLDDT is carried in the B-factor column.

WHY A C-ALPHA TRACE AND NOT A CARTOON
No PyMOL or ChimeraX on this machine, so secondary-structure cartoons are not
available. A pLDDT-coloured backbone trace is the honest alternative: it shows
fold, the position of the variant, and how much of the model can be trusted,
without implying a secondary-structure assignment we did not compute.

CAVEAT ON THE DIMER
The prediction is a dimer, but its interface should not be shown as evidence of
dimerisation: ipTM is 0.45-0.46 across all five models and pTM 0.46-0.48, well
below the ~0.8 usually required to believe an interface, with 43-46% of the
model called disordered. Per-domain pLDDT for chain A:

    signal peptide  1-20     56.4
    extracellular   21-193   72.3   <- the domain the variant sits in
    TM helix        194-211  46.0
    cytoplasmic     212-311  40.0

Only the extracellular domain reaches the "confident" band, which is why the
ECD-only render is the one intended for the figure. The variant neighbourhood
is locally confident: N94 70.5, C95 78.0, T96 79.3.
"""

import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBParser
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from PIL import Image
from scipy.interpolate import splev, splprep

warnings.filterwarnings("ignore")

PDB = Path("data/structure/AF_sid2_wt_dimer/sid2wt.pdb")
OUT = Path("plots/assets")
OUT.mkdir(parents=True, exist_ok=True)

# the AlphaFold pLDDT palette, as used by the EBI viewer
BANDS = [(90, 101, "#0053D6", "Very high (>90)"),
         (70, 90, "#65CBF3", "Confident (70-90)"),
         (50, 70, "#FFDB13", "Low (50-70)"),
         (0, 50, "#FF7D45", "Very low (<50)")]
COL_MARK = "#C4302B"
DOMAINS = {"signal": (1, 20), "extracellular": (21, 193),
           "TM helix": (194, 211), "cytoplasmic": (212, 311)}
# label offsets are per-site: the two are 2 residues apart and their
# labels collided when both were placed the same way
MARK = {94: ("N94", (0.0, 0.0, 4.5), "center", "bottom"),
        96: ("T96", (0.0, 0.0, -4.5), "center", "top")}


def plddt_colour(v):
    for lo, hi, col, _ in BANDS:
        if lo <= v < hi:
            return col
    return BANDS[0][2]


def load_chain(cid):
    st = PDBParser(QUIET=True).get_structure("sid2", PDB)
    ch = [c for c in list(st)[0] if c.id == cid][0]
    res = [r for r in ch if r.id[0] == " " and "CA" in r]
    ids = np.array([r.id[1] for r in res])
    xyz = np.array([r["CA"].coord for r in res], dtype=float)
    plddt = np.array([np.mean([a.get_bfactor() for a in r]) for r in res])
    return ids, xyz, plddt


def orient(xyz):
    """Rotate onto principal axes so the long axis is horizontal.

    A fixed view angle is meaningless for an arbitrary model orientation; this
    makes the framing reproducible instead of hand-tuned.
    """
    c = xyz - xyz.mean(0)
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    return c @ vt.T


def smooth(xyz, per_residue=6):
    """Spline through the C-alpha positions, for a trace that reads as a chain
    rather than as a polyline of straight hops."""
    if len(xyz) < 4:
        return xyz, np.arange(len(xyz), dtype=float)
    tck, u = splprep(xyz.T, s=0, k=3)
    un = np.linspace(0, 1, len(xyz) * per_residue)
    pts = np.array(splev(un, tck)).T
    # index back into residue space, so colours can follow the original pLDDT
    return pts, un * (len(xyz) - 1)


def autocrop(path, margin=8, bg=250):
    """Crop a saved render to its ink.

    mplot3d reserves margins inside the axes that no limit or box-aspect
    setting removes, and bbox_inches cannot recover them because the axes
    patch is invisible. Cropping the raster afterwards is the reliable fix.
    """
    im = Image.open(path).convert("RGB")
    a = np.asarray(im)
    ink = (a < bg).any(axis=2)
    if not ink.any():
        return
    rows, cols = np.where(ink.any(1))[0], np.where(ink.any(0))[0]
    box = (max(cols[0] - margin, 0), max(rows[0] - margin, 0),
           min(cols[-1] + margin + 1, a.shape[1]),
           min(rows[-1] + margin + 1, a.shape[0]))
    im.crop(box).save(path)
    print(f"    cropped to {box[2] - box[0]} x {box[3] - box[1]} px")


def render(ids, xyz, plddt, fname, colour_by="plddt", title=None,
           chain_break=None, elev=18, azim=-60, lw=3.2, figsize=(6.4, 5.4),
           legend=False):
    xyz = orient(xyz)
    fig = plt.figure(figsize=figsize, dpi=600)
    ax = fig.add_subplot(111, projection="3d")

    pieces = [(np.arange(len(ids)), None)] if chain_break is None else chain_break
    for idx, fixed in pieces:
        pts, ridx = smooth(xyz[idx])
        segs = np.stack([pts[:-1], pts[1:]], axis=1)
        if fixed is not None:
            cols = [fixed] * len(segs)
        else:
            src = plddt[idx]
            cols = [plddt_colour(src[int(round(min(i, len(src) - 1)))])
                    for i in ridx[:-1]]
        ax.add_collection3d(Line3DCollection(segs, colors=cols, linewidths=lw,
                                             capstyle="round"))

    # the variant neighbourhood
    for pos, (lab, off, ha, va) in MARK.items():
        w = np.where(ids == pos)[0]
        if not len(w):
            continue
        p = xyz[w[0]]
        ax.scatter(*p, s=110, color=COL_MARK, edgecolors="white",
                   linewidths=0.9, depthshade=False, zorder=10)
        ax.text(p[0] + off[0], p[1] + off[1], p[2] + off[2], lab,
                color=COL_MARK, fontsize=10, fontweight="bold",
                ha=ha, va=va, zorder=11)

    # Frame tightly. A cubic box around the largest dimension wastes most of
    # the canvas on an elongated molecule, and bbox_inches cannot crop it back
    # because the axes patch is invisible. Setting the limits per axis and the
    # box aspect to the same ratios keeps the scaling equal in all three
    # directions while filling the frame.
    pad = 2.0
    lo, hi = xyz.min(0) - pad, xyz.max(0) + pad
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    ax.set_zlim(lo[2], hi[2])
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)
    try:
        ax.set_box_aspect(tuple(hi - lo))
    except Exception:
        pass

    # The key is drawn in R instead, so the asset crops tightly to the
    # molecule and the legend typography matches the rest of the figure.
    if legend and colour_by == "plddt":
        handles = [Line2D([], [], color=c, lw=3, label=l)
                   for _, _, c, l in BANDS]
        handles.append(Line2D([], [], marker="o", color="none",
                              markerfacecolor=COL_MARK, markeredgecolor="white",
                              markersize=7, label="Variant site"))
        ax.legend(handles=handles, loc="upper left", frameon=False,
                  fontsize=7.5, handlelength=1.4, borderaxespad=0.0,
                  labelspacing=0.35, bbox_to_anchor=(0.0, 1.0))
    if title:
        ax.set_title(title, fontsize=10, color="#333333")

    fig.subplots_adjust(0, 0, 1, 1)
    fig.savefig(OUT / fname, dpi=600, bbox_inches="tight", pad_inches=0.02,
                transparent=False, facecolor="white")
    plt.close(fig)
    autocrop(OUT / fname)
    print(f"  wrote {OUT / fname}")


def main():
    ids, xyz, plddt = load_chain("A")
    print(f"chain A: {len(ids)} residues, pLDDT {plddt.min():.1f}-{plddt.max():.1f}")
    for name, (a, b) in DOMAINS.items():
        m = (ids >= a) & (ids <= b)
        print(f"  {name:15s} {a:3d}-{b:3d}  mean pLDDT {plddt[m].mean():5.1f}")

    print("\nfull monomer")
    render(ids, xyz, plddt, "sid2_monomer_plddt.png")

    print("extracellular domain only (21-193)")
    a, b = DOMAINS["extracellular"]
    m = (ids >= a) & (ids <= b)
    render(ids[m], xyz[m], plddt[m], "sid2_ecd_plddt.png",
           lw=3.6, figsize=(6.4, 5.2))

    print("dimer, reference only -- interface is NOT reliable (ipTM 0.46)")
    ia, xa, pa = load_chain("A")
    ib, xb, pb = load_chain("B")
    both = np.vstack([xa, xb])
    idx_a = np.arange(len(xa))
    idx_b = np.arange(len(xa), len(xa) + len(xb))
    render(np.concatenate([ia, ib]), both, np.concatenate([pa, pb]),
           "sid2_dimer_chains.png", colour_by="chain",
           chain_break=[(idx_a, "#2E4057"), (idx_b, "#9FB4C7")], lw=2.4)

    conf = sorted(Path("data/structure/AF_sid2_wt_dimer").glob(
        "*summary_confidences_*.json"))
    vals = [(json.load(open(f)).get("iptm"), json.load(open(f)).get("ptm"))
            for f in conf]
    print("\nmodel confidences (ipTM, pTM):", vals)


if __name__ == "__main__":
    main()
