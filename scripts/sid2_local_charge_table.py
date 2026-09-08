#!/usr/bin/env python3
"""Derive every table the rebuilt Figure 4C panels need.

    cd <repo root> && python3 claude_science/scripts/01_panelC_data.py

Reads only from supplemental_data/structure/ and writes only into
claude_science/data/. Nothing outside claude_science/ is modified.

OUTPUTS  (promoted from the claude_science session; writes into
         supplemental_data/structure/ so the figure scripts read the deposit,
         which is the convention every other figure follows)
  sid2_panelC_residues.tsv   per-residue: pLDDT, SSE, topology, coordinates,
                             Ca distance from T96, SASA, local charge
  sid2_panelC_published.tsv  the published residues with distances and class
  sid2_local_charge.tsv      per-residue net charge of the 12 A neighbourhood
                             at pH 4.4 (gut lumen) and pH 7.4

RESIDUE CLASSES
The uptake-critical set is three extracellular HISTIDINES -- His32, His168,
His175 -- from McEwan, Weisman & Hunter 2012, Mol Cell 47:746
(doi:10.1016/j.molcel.2012.07.014). That paper never mentions residue 34: D34
is the qt13 loss-of-function allele, a separate line of evidence, and is
labelled as such rather than folded into the uptake set. This model's
ectodomain contains exactly three histidines, at 32, 168 and 175, which
independently confirms the numbering.

CHARGE MODEL
Henderson-Hasselbalch side-chain charges, same pKa set as the manuscript's
electrostatics supplement: Asp 3.9, Glu 4.25, His 6.0, Lys 10.5, Arg 12.5,
Cys 8.3, Tyr 10.1. Formal, so it ignores pKa shifts from the local environment
and any glycan shielding; the supplement's Coulombic treatment is the more
careful one and this is used only to rank neighbourhoods against each other.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1 as three2one
from Bio.PDB.SASA import ShrakeRupley
from scipy.stats import binomtest

PDB = Path("supplemental_data/structure/sid2_membrane_oriented.pdb")
PERRES = Path("supplemental_data/structure/sid2_per_residue.tsv")
OUTDIR = Path("supplemental_data/structure")

FOCAL = 96
RADIUS = 12.0           # neighbourhood radius for the local-charge measure
NEAR_CUT = 20.0         # the radius drawn in the panels
CLASSES = {32: "uptake histidine", 168: "uptake histidine",
           175: "uptake histidine", 34: "qt13 allele"}

PKA = {"D": 3.9, "E": 4.25, "H": 6.0, "K": 10.5, "R": 12.5, "C": 8.3, "Y": 10.1}
NEG = set("DECY")


def aa1(resname):
    """Three-letter to one-letter. Biopython's table is keyed on UPPERCASE in
    current releases and Titlecase in older ones, so try both -- silently
    returning 'X' here would zero out every charge."""
    r = resname.strip()
    return three2one.get(r, three2one.get(r.capitalize(), "X"))


def side_chain_charge(aa, pH):
    if aa not in PKA:
        return 0.0
    if aa in NEG:
        return -1.0 / (1.0 + 10 ** (PKA[aa] - pH))
    return 1.0 / (1.0 + 10 ** (pH - PKA[aa]))


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    chain = PDBParser(QUIET=True).get_structure("sid2", str(PDB))[0]["A"]

    res = {r.id[1]: r for r in chain if r.id[0] == " " and "CA" in r}
    ids = np.array(sorted(res))
    one = {i: aa1(res[i].get_resname()) for i in ids}
    assert "X" not in one.values(), "unmapped residue name -- check three2one"
    ca = {i: res[i]["CA"].coord.astype(float) for i in ids}
    heavy = {i: np.array([a.coord for a in res[i] if a.element != "H"], float)
             for i in ids}

    per = pd.read_csv(PERRES, sep="\t").set_index("resid")

    sr = ShrakeRupley()
    sr.compute(chain, level="R")
    sasa = {i: float(res[i].sasa) for i in ids}

    xyz = np.array([ca[i] for i in ids])
    d_focal = np.linalg.norm(xyz - ca[FOCAL], axis=1)
    dmat = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)

    q44 = np.array([side_chain_charge(one[i], 4.4) for i in ids])
    q74 = np.array([side_chain_charge(one[i], 7.4) for i in ids])
    near = dmat <= RADIUS

    tab = pd.DataFrame(dict(
        resid=ids, aa=[one[i] for i in ids],
        plddt=per.loc[ids, "plddt"].to_numpy(float),
        sse=per.loc[ids, "sse"].to_numpy(),
        topology=per.loc[ids, "topology"].to_numpy(),
        x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2],
        d_ca_t96=d_focal, sasa=[sasa[i] for i in ids],
        n_near=near.sum(1), q_local_pH44=near @ q44, q_local_pH74=near @ q74))
    tab.to_csv(OUTDIR / "sid2_panelC_residues.tsv", sep="\t", index=False)
    tab[["resid", "aa", "q_local_pH44", "q_local_pH74", "n_near", "plddt",
         "x", "y", "z"]].to_csv(OUTDIR / "sid2_local_charge.tsv", sep="\t",
                                index=False)

    def nearest_heavy(i):
        return float(np.min(np.linalg.norm(
            heavy[i][:, None, :] - heavy[FOCAL][None, :, :], axis=2)))

    pub = pd.DataFrame([dict(
        resid=i, label=one[i] + str(i), cls=cls,
        d_ca=float(np.linalg.norm(ca[i] - ca[FOCAL])), d_heavy=nearest_heavy(i),
        x=float(ca[i][0]), y=float(ca[i][1]), z=float(ca[i][2]),
        plddt=float(per.loc[i, "plddt"]),
        q_local_pH44=float(tab.q_local_pH44[tab.resid == i].iloc[0]))
        for i, cls in CLASSES.items()]).sort_values("d_ca")
    pub.to_csv(OUTDIR / "sid2_panelC_published.tsv", sep="\t", index=False)

    ## the numbers quoted in the panel subtitles and in the memo
    ecd = tab[tab.resid != FOCAL]
    frac = float((ecd.d_ca_t96 <= NEAR_CUT).mean())
    his = pub[pub.cls == "uptake histidine"]
    k = int((his.d_ca <= NEAR_CUT).sum())
    q96 = float(tab.q_local_pH44[tab.resid == FOCAL].iloc[0])
    print("residues %d-%d (n=%d); histidines present: %s"
          % (ids.min(), ids.max(), len(ids),
             [int(i) for i in ids if one[i] == "H"]))
    print("median Ca distance from T96 %.1f A; %.1f%% within %.0f A"
          % (ecd.d_ca_t96.median(), 100 * frac, NEAR_CUT))
    print("uptake histidines within %.0f A: %d of %d, binomial p = %.2f"
          % (NEAR_CUT, k, len(his),
             binomtest(k, len(his), frac, alternative="greater").pvalue))
    print("ectodomain net charge: pH 4.4 %+.2f e, pH 7.4 %+.2f e"
          % (q44.sum(), q74.sum()))
    print("local charge pH 4.4: T96 %+.2f e (%.0fth pct) -> T96K %+.2f e (%.0fth pct)"
          % (q96, 100 * (tab.q_local_pH44 < q96).mean(),
             q96 + 1, 100 * (tab.q_local_pH44 < q96 + 1).mean()))
    print("T96 pLDDT %.1f, SASA %.1f A^2 (%.0fth pct)"
          % (per.loc[FOCAL, "plddt"], sasa[FOCAL],
             100 * np.mean([sasa[i] < sasa[FOCAL] for i in ids])))
    print("wrote 3 tables to", OUTDIR)


if __name__ == "__main__":
    main()
