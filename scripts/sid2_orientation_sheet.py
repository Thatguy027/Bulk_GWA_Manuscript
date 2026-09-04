#!/usr/bin/env python3
"""Contact sheets of candidate orientations for the Figure 4 structure panels.

    python3 scripts/sid2_orientation_sheet.py

Writes to plots/assets/:
    sid2_orient_overview_sheet.png   the ectodomain overview, 10 orientations
    sid2_orient_zoom_sheet.png       the T96 zoom, 12 orientations

Every tile is labelled with its elevation and azimuth. Pick one and pass it
back: scripts/sid2_zoom_render.py takes ELEV/AZIM for each panel at the top, so
choosing an orientation is a two-number edit and a re-render at full quality.

The membrane frame is preserved throughout -- z is the membrane normal and +z
is the lumen -- so elevation is tilt relative to the membrane plane and
azimuth is rotation about the membrane normal. Azimuth therefore changes which
face of the ectodomain points at the reader without ever tipping the lumen away
from up, which is the property that makes these comparable.

Tiles are rendered at low resolution on purpose; they are for choosing, not for
publication.
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
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

sys.path.insert(0, "scripts")
import sid2_ribbon_render as rr
import sid2_zoom_render as zm

warnings.filterwarnings("ignore")

OUT = Path("plots/assets")
COL_T96 = zm.COL_T96
COL_FUNC = zm.COL_FUNC
MEMBRANE = zm.MEMBRANE

## current settings, shown first for reference
OVERVIEW_NOW = (8, -72)
ZOOM_NOW = (6, -104)

## AZIMUTH CANNOT BRING T96 FORWARD, and it is worth knowing why before
## picking. T96's Ca is 3.0 A from the domain axis in the membrane plane -- the
## 5th percentile of the domain's in-plane radius -- and +15.0 A along the
## membrane normal. Rotating about the normal therefore moves it barely at
## all: 180 degrees of azimuth changes its depth toward the reader by about
## 1 A. Elevation is the control that works, taking it from +0.7 A at elev 8
## to +15.0 A at elev 88, because that is looking down onto the lumenal face.
## The cost is that "lumen up" becomes "lumen toward the reader".
OVERVIEW_GRID = ([(8, -72)]
                 + [(e, 108) for e in (8, 25, 40, 55, 70, 88)]
                 + [(e, 93) for e in (25, 45, 65)])
## "up and around to the right slightly" from the current zoom. Here the four
## residues are spread, so azimuth does separate them.
ZOOM_GRID = [(e, a) for e in (6, 18, 30, 45)
             for a in (-119, -104, -89, -74)]


def cam_dir(elev, azim):
    """Unit vector from the scene toward the matplotlib 3D camera."""
    e, a = np.radians(elev), np.radians(azim)
    return np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])


def depth_of(p, centre, elev, azim):
    """Signed depth of a point toward the viewer, in Angstroms.

    Needed because the residue markers are drawn with depthshade off and a
    forced zorder -- mplot3d will otherwise bury them under the ribbon -- so
    they LOOK frontmost in every orientation and the picture alone cannot tell
    you which side of the domain the residue is on. This is the honest cue.
    """
    return float(np.dot(p - centre, cam_dir(elev, azim)))


def draw_overview(ax, ids, CA, sse, elev, azim):
    quads, norms, ridx = zm.ribbon_quads(CA, sse)
    cols = [rr.shade(rr.SSE_COL[sse[r]], n) for r, n in zip(ridx, norms)]
    ax.add_collection3d(Poly3DCollection(quads, facecolors=cols,
                                         edgecolors="none", shade=False))
    x0, x1 = CA[:, 0].min() - 5, CA[:, 0].max() + 5
    ym = CA[:, 1].mean()
    ax.add_collection3d(Poly3DCollection(
        [[(x0, ym, -13.0), (x1, ym, -13.0), (x1, ym, 13.0), (x0, ym, 13.0)]],
        facecolors=MEMBRANE, edgecolors="none", alpha=0.85, zorder=0))
    p96 = CA[ids == zm.FOCAL][0]
    d = depth_of(p96, CA.mean(0), elev, azim)
    ## solid when T96 faces the reader, open when it is on the far side
    ax.scatter(*p96, s=40,
               color=COL_T96 if d > 0 else "white",
               edgecolors=COL_T96 if d <= 0 else "white",
               linewidths=1.1, depthshade=False, zorder=15)
    ax.text2D(0.5, -0.02, f"T96 {'front' if d > 0 else 'BACK'} ({d:+.0f} Å)",
              transform=ax.transAxes, ha="center", fontsize=7,
              color=COL_T96 if d > 0 else "#999999")
    allp = np.vstack([np.array([q for quad in quads for q in quad]),
                      np.array([[x0, ym, -13.0], [x1, ym, 13.0]])])
    zm.frame_axes(ax, allp, pad=1.0, elev=elev, azim=azim)


def draw_zoom(ax, res, ids, CA, sse, elev, azim):
    p96 = CA[ids == zm.FOCAL][0]
    keep = np.where(np.linalg.norm(CA - p96, axis=1) <= zm.ZOOM_RADIUS - 1)[0]
    quads, norms, _ = zm.ribbon_quads(CA[keep], sse[keep])
    cols = [rr.shade("#C7D2DA", n) for n in norms]
    ax.add_collection3d(Poly3DCollection(quads, facecolors=cols,
                                         edgecolors="none", shade=False,
                                         alpha=0.85))
    show = {zm.FOCAL: COL_T96}
    show.update({p: COL_FUNC for p in zm.FUNC})
    centre = CA[keep].mean(0)
    for pos, col in show.items():
        r = res[int(np.where(ids == pos)[0][0])]
        zm.sticks(r, ax, scale=0.7)
        p = r["CA"].coord
        d = depth_of(p, centre, elev, azim)
        ax.scatter(*p, s=24, color=col if d > 0 else "white",
                   edgecolors=col if d <= 0 else "white", linewidths=0.9,
                   depthshade=False, zorder=13)
        if pos != zm.FOCAL:
            ax.plot(*np.array([p96, p]).T, color=COL_FUNC, linewidth=0.8,
                    linestyle=(0, (3, 2.5)), zorder=11)
    d96 = depth_of(p96, centre, elev, azim)
    spread = max(depth_of(res[int(np.where(ids == p)[0][0])]["CA"].coord,
                          centre, elev, azim) for p in zm.FUNC)
    ax.text2D(0.5, -0.02,
              f"T96 {d96:+.0f} Å, nearest partner {spread:+.0f} Å",
              transform=ax.transAxes, ha="center", fontsize=7,
              color="#555555")
    pts = np.array([q for quad in quads for q in quad])
    zm.frame_axes(ax, pts, pad=1.0, elev=elev, azim=azim)


def sheet(grid, ncol, draw, fname, title, figsize):
    nrow = int(np.ceil(len(grid) / ncol))
    fig = plt.figure(figsize=figsize, dpi=170)
    for i, (elev, azim) in enumerate(grid):
        ax = fig.add_subplot(nrow, ncol, i + 1, projection="3d")
        draw(ax, elev, azim)
        tag = f"elev {elev}, azim {azim}"
        if (elev, azim) in (OVERVIEW_NOW, ZOOM_NOW):
            tag += "  (current)"
        ax.set_title(tag, fontsize=7.5, color="#333333", pad=1)
    fig.suptitle(title, fontsize=9.5, y=0.995, color="#111111")
    fig.subplots_adjust(0.01, 0.01, 0.99, 0.95, wspace=0.02, hspace=0.10)
    fig.savefig(OUT / fname, dpi=170, facecolor="white",
                bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"  wrote {OUT / fname}  ({len(grid)} orientations)")


def main():
    res, ids, N, CA, C, O = zm.load("A")
    sse = rr.dssp_lite(N, CA, C, O)
    print(f"chain A: {len(ids)} residues, membrane frame (+z = lumen)")

    print("overview sheet")
    sheet(OVERVIEW_GRID, 5,
          lambda ax, e, a: draw_overview(ax, ids, CA, sse, e, a),
          "sid2_orient_overview_sheet.png",
          "SID-2 ectodomain overview -- candidate orientations. T96 sits on "
          "the rotation axis, so azimuth barely moves it; elevation is what "
          "brings it forward, at the cost of tipping the lumen toward you.",
          figsize=(13.0, 6.2))

    print("zoom sheet")
    sheet(ZOOM_GRID, 4,
          lambda ax, e, a: draw_zoom(ax, res, ids, CA, sse, e, a),
          "sid2_orient_zoom_sheet.png",
          "T96 zoom -- candidate orientations "
          "(T96 orange, uptake-critical residues dark, dashed = Ca-Ca)",
          figsize=(12.0, 11.0))


if __name__ == "__main__":
    main()
