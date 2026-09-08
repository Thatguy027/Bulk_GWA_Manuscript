#!/usr/bin/env python3
"""Is the SID-2 ectodomain a SID-1 beta-strand-rich domain? TM-align says no.

    cd <repo root> && python3 claude_science/scripts/04_fold_comparison.py

Needs network access on first run (fetches from files.rcsb.org into
claude_science/data/pdb_cache/) and the `tmtools` package:  pip install tmtools

Writes claude_science/data/sid2_vs_sid1_tmalign.csv; nothing elsewhere in the
repository is modified.

WHY
Wang, Cong, Qian, Yan & Gong 2024, Nucleic Acids Res 52:6718
(doi:10.1093/nar/gkae395) solved C. elegans SID-1 by cryo-EM at 2.21 A: PDB
8XBS apo and 8XC1 with dsRNA. Its extracellular region is two beta-strand-rich
domains, BRD1 (residues 18-178) and BRD2 (179-310), and the dsRNA binds at
their interface. SID-2's ectodomain is a similar size and similarly beta-rich,
which invites the guess that it is the same fold.

THE POINT OF THE CONTROLS
A TM-score against BRD1 alone is uninterpretable. Two references make it
readable: unrelated beta-sandwich domains of comparable size, which set the
background any compact beta domain scores; and BRD1 against BRD2, a genuine
structural relationship inside SID-1, which shows what a real match looks like
in this comparison set. The query is the well-modelled core of the SID-2
prediction (pLDDT >= 70), so the controls are scored against the same
imperfect query and the model's limitations do not bias the comparison.
"""

import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1 as three2one
from tmtools import tm_align

SID2_PDB = Path("supplemental_data/structure/sid2_membrane_oriented.pdb")
PERRES = Path("supplemental_data/structure/sid2_per_residue.tsv")
CACHE = Path("claude_science/data/pdb_cache")
OUT = Path("claude_science/data/sid2_vs_sid1_tmalign.csv")

PLDDT_MIN = 70.0
SID1 = "8XC1"                      # cSID1 + dsRNA; chain A is one subunit
SID1_DOMAINS = {"cSID1 BRD1 (18-178)": (18, 178),
                "cSID1 BRD2 (179-310)": (179, 310),
                "cSID1 TMD (311-776)": (311, 776)}
## unrelated beta-rich domains: fibronectin type III, an immunoglobulin
## domain, a galectin carbohydrate-recognition domain and a legume lectin
CONTROLS = {"fn3 1TEN": "1TEN", "Ig 12E8": "12E8",
            "galectin 2JJ6": "2JJ6", "lectin 1LOB": "1LOB"}


def aa1(resname):
    r = resname.strip()
    return three2one.get(r, three2one.get(r.capitalize(), "X"))


def ca_seq(chain, lo=None, hi=None):
    xyz, seq = [], []
    for r in chain:
        if r.id[0] != " " or "CA" not in r:
            continue
        if lo is not None and not (lo <= r.id[1] <= hi):
            continue
        xyz.append(r["CA"].coord)
        seq.append(aa1(r.get_resname()))
    return np.array(xyz, float), "".join(seq)


def fetch(pdb_id):
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / f"{pdb_id}.cif"
    if not dest.exists():
        urllib.request.urlretrieve(
            f"https://files.rcsb.org/download/{pdb_id}.cif", dest)
        print(f"  fetched {pdb_id}")
    return dest


def biggest_protein_chain(model):
    return max(model, key=lambda c: sum(
        1 for r in c if r.id[0] == " " and "CA" in r))


def score(Xq, Sq, Xt, St, label, role):
    r = tm_align(Xq, Xt, Sq, St)
    pairs = [(a, b) for a, b in zip(r.seqxA, r.seqyA) if a != "-" and b != "-"]
    return dict(target=label, n_target=len(St),
                tm_norm_sid2=round(r.tm_norm_chain1, 3),
                tm_norm_target=round(r.tm_norm_chain2, 3),
                rmsd=round(r.rmsd, 2), aligned=len(pairs),
                pct_id=round(100 * sum(a == b for a, b in pairs) / len(pairs), 1),
                role=role)


def main():
    sid2 = PDBParser(QUIET=True).get_structure("sid2", str(SID2_PDB))[0]["A"]
    per = pd.read_csv(PERRES, sep="\t").set_index("resid")
    keep = [r.id[1] for r in sid2
            if r.id[0] == " " and "CA" in r
            and per.loc[r.id[1], "plddt"] >= PLDDT_MIN]
    Xq = np.array([sid2[(" ", i, " ")]["CA"].coord for i in keep], float)
    Sq = "".join(aa1(sid2[(" ", i, " ")].get_resname()) for i in keep)
    print("query: SID-2 ectodomain core, %d residues at pLDDT >= %.0f"
          % (len(keep), PLDDT_MIN))

    mdl = MMCIFParser(QUIET=True).get_structure(SID1, str(fetch(SID1)))[0]
    rows = []
    for label, (lo, hi) in SID1_DOMAINS.items():
        Xt, St = ca_seq(mdl["A"], lo, hi)
        role = ("query vs SID-1 domain" if "BRD" in label
                else "negative control")
        rows.append(score(Xq, Sq, Xt, St, label, role))

    for label, pid in CONTROLS.items():
        m = MMCIFParser(QUIET=True).get_structure(pid, str(fetch(pid)))[0]
        ch = biggest_protein_chain(m)
        Xt, St = ca_seq(ch)
        rows.append(score(Xq, Sq, Xt, St, f"{label} chain {ch.id}",
                          "negative control"))

    ## the internal reference: the two genuine BRDs against each other
    Xb1, Sb1 = ca_seq(mdl["A"], 18, 178)
    Xb2, Sb2 = ca_seq(mdl["A"], 179, 310)
    rows.append(score(Xb1, Sb1, Xb2, Sb2,
                      "cSID1 BRD1 vs BRD2 (internal reference)",
                      "internal reference"))

    tab = pd.DataFrame(rows).sort_values("tm_norm_sid2", ascending=False)
    tab.to_csv(OUT, index=False)
    print(tab.to_string(index=False))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
