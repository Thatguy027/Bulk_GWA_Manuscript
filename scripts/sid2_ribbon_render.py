#!/usr/bin/env python3
"""Ribbon (cartoon) renders of the AlphaFold3 SID-2 model for Figure 4.

    python3 scripts/sid2_ribbon_render.py

Writes to plots/assets/:
    sid2_ecd_ribbon_plddt.png   extracellular domain, ribbon coloured by pLDDT
    sid2_ecd_ribbon_sse.png     the same, coloured by secondary structure
    sid2_monomer_ribbon_plddt.png   full chain A

Replaces the C-alpha trace in scripts/sid2_structure_render.py, which is kept
for the dimer reference image.

SECONDARY STRUCTURE
No DSSP binary and no biotite on this machine, so secondary structure is
assigned here from the backbone using the Kabsch & Sander (1983) electrostatic
hydrogen-bond criterion -- the same definition DSSP uses, not a geometric
approximation of it. The AlphaFold PDB carries N, CA, C and O for every
residue, which is all the criterion needs; the amide hydrogen is placed at
1 A from N along the C(i-1)->O(i-1) direction, as DSSP does.

    E = 332 * 0.42 * 0.20 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN)   kcal/mol
    hydrogen bond if E < -0.5

    alpha helix   two consecutive 4-turns
    strand        residues in a parallel or antiparallel bridge

Only H (helix) and E (strand) are separated from coil; 3-10 and pi helices are
folded into H. Isolated single-residue assignments are dropped, as DSSP also
does, because one residue cannot make a helix.

RIBBON GEOMETRY
A flat ribbon is swept along a spline through the C-alpha positions. The ribbon
plane at each residue is set by the local binormal, cross(CA(i+1)-CA(i),
CA(i)-CA(i-1)), sign-corrected against its predecessor so the band does not
flip; for a helix this is the vector that makes the familiar twisted band.
Width follows the assignment -- wide for helix and strand, thin for coil -- and
strands taper to an arrowhead over their last two residues. Faces are shaded by
the angle between the surface normal and a fixed light, which is what makes a
flat matplotlib surface read as a three-dimensional object.

CAVEAT ON THE MODEL: unchanged from the trace version -- ipTM 0.45-0.46 means
the dimer interface is not credible and is not drawn, and only the
extracellular domain (mean pLDDT 72.3) reaches the confident band.
"""

import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBParser
import matplotlib.patheffects as pe
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PIL import Image
from scipy.interpolate import splev, splprep

warnings.filterwarnings("ignore")

PDB = Path("supplemental_data/structure/sid2_alphafold_dimer.pdb")
OUT = Path("plots/assets")
OUT.mkdir(parents=True, exist_ok=True)

PLDDT_BANDS = [(90, 1e9, "#0053D6"), (70, 90, "#65CBF3"),
               (50, 70, "#FFDB13"), (-1e9, 50, "#FF7D45")]
SSE_COL = {"H": "#7E9BB5", "E": "#4E6E8A", "C": "#C4D0D8"}
## residue-class colours, for the charge-chemistry render. His is given its own
## colour because it is the class that titrates between neutral pH and the
## acidic intestinal lumen, which is the point of that figure.
CHARGE_COL = {"ARG": "#3B6FB6", "LYS": "#3B6FB6", "HIS": "#8E6BAF",
              "ASP": "#C0392B", "GLU": "#C0392B"}
CHARGE_NEUTRAL = "#D5DCE1"
COL_MARK = "#F34C00"   # the JU1793 strain colour: the model is the 96T allele
## white fill with a dark ring: a filled dark marker disappeared against the
## dark ribbon, and any single colour fails somewhere on a coloured surface
COL_FUNC = "#FFFFFF"
COL_FUNC_EDGE = "#16324A"
## default marks: the variant plus the neighbouring loop residue
MARK = {94: ("N94", (0, 0, 5.0), "center", "bottom", None),
        96: ("T96", (0, 0, -5.0), "center", "top", None)}
## for the main display: the variant against the published uptake-critical
## residues that share its face
## offsets are applied in the ORIENTED frame (x = long axis), so they are
## reproducible; the four labels overlapped when all were offset the same way
MARK_FUNC = {
    96:  ("T96",  (7.0, 0, -5.0), "left", "top", None),
    32:  ("H32",  (0.0, 0, 6.5), "center", "bottom", COL_FUNC),
    34:  ("D34",  (-7.0, 0, -6.0), "right", "top", COL_FUNC),
    168: ("H168", (-11.0, 0, 3.0), "right", "bottom", COL_FUNC)}
WIDTH = {"H": 1.90, "E": 1.60, "C": 0.32}
LIGHT = np.array([0.35, 0.35, 0.87])


# ---------------------------------------------------------------- structure
def backbone(cid="A"):
    st = PDBParser(QUIET=True).get_structure("sid2", PDB)
    ch = [c for c in list(st)[0] if c.id == cid][0]
    res = [r for r in ch if r.id[0] == " "
           and all(a in r for a in ("N", "CA", "C", "O"))]
    get = lambda nm: np.array([r[nm].coord for r in res], dtype=float)
    plddt = np.array([np.mean([a.get_bfactor() for a in r]) for r in res])
    ids = np.array([r.id[1] for r in res])
    return ids, get("N"), get("CA"), get("C"), get("O"), plddt


def dssp_lite(N, CA, C, O):
    """Kabsch & Sander hydrogen bonds -> H / E / C per residue."""
    n = len(CA)
    # amide H, 1 A from N along the C(i-1)->O(i-1) direction
    H = np.full_like(N, np.nan)
    d = C[:-1] - O[:-1]
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    H[1:] = N[1:] + d

    def dist(a, b):
        return np.linalg.norm(a - b, axis=-1)

    hb = np.zeros((n, n), bool)          # hb[donor, acceptor]
    for i in range(1, n):
        if np.isnan(H[i]).any():
            continue
        r_ON = dist(O, N[i]); r_CH = dist(C, H[i])
        r_OH = dist(O, H[i]); r_CN = dist(C, N[i])
        with np.errstate(divide="ignore", invalid="ignore"):
            E = 332.0 * 0.42 * 0.20 * (1 / r_ON + 1 / r_CH - 1 / r_OH - 1 / r_CN)
        E[np.abs(np.arange(n) - i) < 2] = 0.0   # no local self-bonds
        hb[i] = E < -0.5

    turn4 = np.zeros(n, bool)
    for i in range(n - 4):
        turn4[i] = hb[i + 4, i]

    sse = np.array(["C"] * n, dtype="<U1")
    for i in range(n - 5):
        if turn4[i] and turn4[i + 1]:
            sse[i + 1:i + 5] = "H"

    for i in range(1, n - 1):
        for j in range(i + 3, n - 1):
            anti = (hb[j, i] and hb[i, j]) or (hb[i + 1, j - 1] and hb[j + 1, i - 1])
            para = (hb[j, i - 1] and hb[i + 1, j]) or (hb[i, j - 1] and hb[j + 1, i])
            if anti or para:
                for k in (i, j):
                    if sse[k] == "C":
                        sse[k] = "E"

    # a single residue cannot be a helix or a strand
    out = sse.copy()
    k = 0
    while k < n:
        j = k
        while j + 1 < n and sse[j + 1] == sse[k]:
            j += 1
        if sse[k] in "HE" and (j - k + 1) < 2:
            out[k:j + 1] = "C"
        k = j + 1
    return out


# ------------------------------------------------------------------ geometry
def frames(ca):
    """Per-residue ribbon-plane vector, sign-corrected so the band cannot flip.

    Uses a CENTRAL difference for the tangent and a second difference for the
    curvature, both clamped at the ends. An earlier version took the forward
    difference, which is the zero vector at the last residue: the primary cross
    product and its fallback were then both zero and the function returned NaN
    there. That silently corrupted the final quad of every ribbon, and broke
    outright once axis limits were computed from the mesh.
    """
    n = len(ca)
    b = np.zeros_like(ca)
    for i in range(n):
        ip, im = min(i + 1, n - 1), max(i - 1, 0)
        t = ca[ip] - ca[im]
        if np.linalg.norm(t) < 1e-6:
            t = ca[min(i + 2, n - 1)] - ca[i]
        k = ca[ip] - 2 * ca[i] + ca[im]
        v = np.cross(t, k)
        if np.linalg.norm(v) < 1e-6:
            tn = t / max(np.linalg.norm(t), 1e-12)
            ref = np.array([0.0, 0.0, 1.0]) if abs(tn[2]) < 0.9 \
                else np.array([1.0, 0.0, 0.0])
            v = np.cross(t, ref)
        b[i] = v / np.linalg.norm(v)
    for i in range(1, n):
        if np.dot(b[i], b[i - 1]) < 0:
            b[i] = -b[i]
    assert np.isfinite(b).all(), "non-finite ribbon frame"
    return b


def widths(sse):
    """Half-width per residue, with strands tapering to an arrowhead."""
    w = np.array([WIDTH[s] for s in sse], dtype=float)
    n = len(sse)
    k = 0
    while k < n:
        j = k
        while j + 1 < n and sse[j + 1] == sse[k]:
            j += 1
        if sse[k] == "E" and (j - k + 1) >= 3:
            w[k:j + 1] = WIDTH["E"]
            w[j - 1] = WIDTH["E"] * 1.55      # arrow shoulder
            w[j] = WIDTH["E"] * 0.12          # tip
        k = j + 1
    return w


def resample(ca, b, w, per=10):
    tck, _ = splprep(ca.T, s=0, k=3)
    u = np.linspace(0, 1, (len(ca) - 1) * per + 1)
    pts = np.array(splev(u, tck)).T
    idx = u * (len(ca) - 1)
    lo = np.clip(np.floor(idx).astype(int), 0, len(ca) - 1)
    hi = np.clip(lo + 1, 0, len(ca) - 1)
    f = (idx - lo)[:, None]
    side = b[lo] * (1 - f) + b[hi] * f
    side /= np.linalg.norm(side, axis=1, keepdims=True)
    wid = w[lo] * (1 - f[:, 0]) + w[hi] * f[:, 0]
    return pts, side, wid, idx


def shade(rgb, normal):
    lam = LIGHT / np.linalg.norm(LIGHT)
    s = 0.55 + 0.45 * abs(float(np.dot(normal, lam)))
    return tuple(np.clip(np.array(matplotlib.colors.to_rgb(rgb)) * s, 0, 1))


def plddt_col(v):
    for lo, hi, c in PLDDT_BANDS:
        if lo <= v < hi:
            return c
    return PLDDT_BANDS[0][2]


def orient(xyz, ref):
    c = ref - ref.mean(0)
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    return (xyz - ref.mean(0)) @ vt.T


def autocrop(path, margin=8, bg=250):
    im = Image.open(path).convert("RGB")
    a = np.asarray(im)
    ink = (a < bg).any(axis=2)
    rows, cols = np.where(ink.any(1))[0], np.where(ink.any(0))[0]
    im.crop((max(cols[0] - margin, 0), max(rows[0] - margin, 0),
             min(cols[-1] + margin + 1, a.shape[1]),
             min(rows[-1] + margin + 1, a.shape[0]))).save(path)


# -------------------------------------------------------------------- render
def render(ids, ca, plddt, sse, fname, colour="plddt", elev=18, azim=-60,
           figsize=(7.2, 4.6), marks=None, resnames=None):
    b = frames(ca)
    w = widths(sse)
    pts, side, wid, idx = resample(ca, b, w)
    pts = orient(pts, ca)
    ca_o = orient(ca, ca)
    side = side @ np.linalg.svd(ca - ca.mean(0), full_matrices=False)[2].T

    quads, cols = [], []
    for k in range(len(pts) - 1):
        p0, p1 = pts[k], pts[k + 1]
        s0, s1 = side[k] * wid[k], side[k + 1] * wid[k + 1]
        quad = [p0 - s0, p0 + s0, p1 + s1, p1 - s1]
        tangent = p1 - p0
        nrm = np.cross(tangent, side[k])
        nl = np.linalg.norm(nrm)
        nrm = nrm / nl if nl > 1e-9 else np.array([0.0, 0.0, 1.0])
        r = int(round(min(idx[k], len(ids) - 1)))
        if colour == "plddt":
            base = plddt_col(plddt[r])
        elif colour == "charge":
            base = CHARGE_COL.get(resnames[r], CHARGE_NEUTRAL)
        else:
            base = SSE_COL[sse[r]]
        quads.append(quad)
        cols.append(shade(base, nrm))

    fig = plt.figure(figsize=figsize, dpi=600)
    ax = fig.add_subplot(111, projection="3d")
    ax.add_collection3d(Poly3DCollection(quads, facecolors=cols,
                                         edgecolors="none", shade=False))

    for pos, (lab, off, ha, va, col) in (marks or MARK).items():
        wq = np.where(ids == pos)[0]
        if not len(wq):
            continue
        p = ca_o[wq[0]]
        variant = col is None
        c = COL_MARK if variant else col
        ax.scatter(*p, s=115 if variant else 90, color=c,
                   edgecolors="white" if variant else COL_FUNC_EDGE,
                   linewidths=1.0 if variant else 1.4, depthshade=False,
                   zorder=10)
        ## leader line from the residue to its label. mplot3d depth-sorts the
        ## ribbon surface over the markers whatever the zorder, so an occluded
        ## marker plus an offset label is ambiguous without one.
        tip = (p[0] + off[0] * 0.82, p[1] + off[1] * 0.82,
               p[2] + off[2] * 0.82)
        ax.plot([p[0], tip[0]], [p[1], tip[1]], [p[2], tip[2]],
                color=COL_MARK if variant else COL_FUNC_EDGE,
                linewidth=0.9, solid_capstyle="round", zorder=10.5,
                path_effects=[pe.withStroke(linewidth=2.4,
                                            foreground="white")])
        ## a white stroke round the glyphs, so labels stay legible wherever
        ## they land on the ribbon
        ax.text(p[0] + off[0], p[1] + off[1], p[2] + off[2], lab,
                color=COL_MARK if variant else COL_FUNC_EDGE,
                fontsize=10.5 if variant else 9.5, fontweight="bold",
                ha=ha, va=va, zorder=11,
                path_effects=[pe.withStroke(linewidth=2.6,
                                            foreground="white")])

    allpts = np.vstack([pts, ca_o])
    pad = 2.5
    lo, hi = allpts.min(0) - pad, allpts.max(0) + pad
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)
    try:
        ax.set_box_aspect(tuple(hi - lo))
    except Exception:
        pass
    fig.subplots_adjust(0, 0, 1, 1)
    fig.savefig(OUT / fname, dpi=600, bbox_inches="tight", pad_inches=0.02,
                facecolor="white")
    plt.close(fig)
    autocrop(OUT / fname)
    print(f"  wrote {OUT / fname}")


def main():
    ids, N, CA, C, O, plddt = backbone("A")
    sse = dssp_lite(N, CA, C, O)
    frac = {s: float(np.mean(sse == s)) for s in "HEC"}
    print(f"chain A: {len(ids)} residues | helix {frac['H']:.0%} "
          f"strand {frac['E']:.0%} coil {frac['C']:.0%}")
    for a, b_, nm in [(1, 20, "signal"), (21, 193, "extracellular"),
                      (194, 211, "TM helix"), (212, 311, "cytoplasmic")]:
        m = (ids >= a) & (ids <= b_)
        print(f"  {nm:15s} {a:3d}-{b_:3d}  helix {np.mean(sse[m]=='H'):.0%}"
              f"  strand {np.mean(sse[m]=='E'):.0%}"
              f"  mean pLDDT {plddt[m].mean():.1f}")
    for p in (94, 96):
        i = int(np.where(ids == p)[0][0])
        print(f"  residue {p}: secondary structure {sse[i]}, pLDDT {plddt[i]:.1f}")

    # export the per-residue assignment so the R figure draws the same track
    topo = {}
    tl = Path("supplemental_data/structure/sid2_deeptmhmm_topology.3line")
    if tl.exists():
        code = {"S": "Signal peptide", "O": "Extracellular",
                "M": "TM helix", "I": "Cytoplasmic"}
        line = open(tl).read().split("\n")[2]
        topo = {i + 1: code[c] for i, c in enumerate(line)}
    tsv = Path("supplemental_data/structure/sid2_per_residue.tsv")
    with open(tsv, "w") as fh:
        fh.write("resid\tresname\tplddt\tsse\ttopology\tx\ty\tz\n")
        st = PDBParser(QUIET=True).get_structure("sid2", PDB)
        ch = [c for c in list(st)[0] if c.id == "A"][0]
        names = {r.id[1]: r.get_resname() for r in ch if r.id[0] == " "}
        for k, rid in enumerate(ids):
            x, y, zc = CA[k]
            fh.write(f"{rid}\t{names.get(rid,'')}\t{plddt[k]:.2f}\t{sse[k]}"
                     f"\t{topo.get(int(rid),'')}\t{x:.3f}\t{y:.3f}"
                     f"\t{zc:.3f}\n")
    print(f"\nwrote {tsv}")

    m = (ids >= 21) & (ids <= 193)
    print("extracellular domain, ribbon coloured by pLDDT")
    render(ids[m], CA[m], plddt[m], sse[m], "sid2_ecd_ribbon_plddt.png")
    print("extracellular domain, ribbon coloured by secondary structure")
    render(ids[m], CA[m], plddt[m], sse[m], "sid2_ecd_ribbon_sse.png",
           colour="sse")
    print("full monomer, ribbon coloured by pLDDT")
    render(ids, CA, plddt, sse, "sid2_monomer_ribbon_plddt.png",
           figsize=(7.6, 4.4))

    st = PDBParser(QUIET=True).get_structure("sid2", PDB)
    chA = [c for c in list(st)[0] if c.id == "A"][0]
    nm = {r.id[1]: r.get_resname() for r in chA if r.id[0] == " "}
    rn = np.array([nm[i] for i in ids])

    print("extracellular domain, variant against the uptake-critical residues")
    render(ids[m], CA[m], plddt[m], sse[m], "sid2_ecd_ribbon_func.png",
           colour="sse", marks=MARK_FUNC)
    print("extracellular domain, ribbon coloured by residue class")
    render(ids[m], CA[m], plddt[m], sse[m], "sid2_ecd_ribbon_charge.png",
           colour="charge", marks=MARK_FUNC, resnames=rn[m])


if __name__ == "__main__":
    main()
