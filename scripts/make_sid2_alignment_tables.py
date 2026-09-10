#!/usr/bin/env python3
"""Stage the C. elegans / C. briggsae SID-2 alignment and its overlays.

  python3 scripts/make_sid2_alignment_tables.py
    -> supplemental_data/structure/sid2_species_alignment.tsv
       supplemental_data/structure/sid2_population_missense.tsv
       supplemental_data/structure/sid2_sequences.fa

WHY THIS EXISTS
C. briggsae is insensitive to environmental RNAi, and a C. elegans sid-2
transgene confers sensitivity on it (Winston et al. 2007), so the two species'
SID-2 proteins are the natural comparison for asking which residues matter.
Neither sequence is in this repository and the alignment is not either, so both
are staged here and the figure reads only these tables.

SOURCES, all outside the repository:
  UniProt G5EEV9  C. elegans SID-2, 311 aa
  UniProt A8XSB8  C. briggsae CBR-SID-2 (CBG18280), 314 aa
  the BCSQ-annotated CeNDR release (20231213) for population missense variants
  supplemental_data/structure/sid2_deeptmhmm_topology.3line for the topology
  supplemental_data/structure/sid2_variants_cendr.tsv for per-parent residues

ALIGNMENT. Global Needleman-Wunsch, BLOSUM62, gap open -11 extend -1
(Biopython PairwiseAligner defaults for protein). The two proteins are only
~47% identical, so the alignment carries real uncertainty in the low-complexity
stretches; the per-position table records the aligned Cb residue and a gap flag
so a reader can see where the comparison is weak rather than trusting it
uniformly.

ON THE "UPTAKE HISTIDINES". H32, H168 and H175 are the residues McEwan et al.
2012 mutated. That paper's triple His->Arg mutant internalised MORE dsRNA than
wild type, not less, so these positions are implicated in dsRNA handling
WITHOUT being required for it, and this file does not label them as critical.
"""
import subprocess, sys, os, datetime
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

OUT = "supplemental_data/structure"
ST  = OUT
VCF = os.environ.get("CENDR_BCSQ",
      "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213/bcsq.vcf.gz")
SID2_REGION = "III:13679445-13682450"
ACC = {"Ce_SID2_G5EEV9": "G5EEV9", "Cb_SID2_A8XSB8": "A8XSB8"}

def fetch(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    r = subprocess.run(["curl","-sS","--max-time","30",url],
                       capture_output=True, text=True)
    if r.returncode or not r.stdout.startswith(">"):
        sys.exit(f"could not fetch {acc}: {r.stderr.strip()[:200]}")
    return "".join(r.stdout.split("\n")[1:]).strip()

seqs = {name: fetch(acc) for name, acc in ACC.items()}
ce, cb = seqs["Ce_SID2_G5EEV9"], seqs["Cb_SID2_A8XSB8"]
assert len(ce) == 311 and len(cb) == 314, (len(ce), len(cb))

with open(f"{OUT}/sid2_sequences.fa","w") as fh:
    for name, s in seqs.items():
        fh.write(f">{name} {ACC[name]} len={len(s)} fetched={datetime.date.today()}\n")
        for i in range(0, len(s), 60): fh.write(s[i:i+60] + "\n")

al = Align.PairwiseAligner()
al.substitution_matrix = substitution_matrices.load("BLOSUM62")
al.open_gap_score, al.extend_gap_score = -11, -1
aln = al.align(ce, cb)[0]
A, B = aln[0], aln[1]

# topology, per C. elegans residue
topo_lines = open(f"{ST}/sid2_deeptmhmm_topology.3line").read().split("\n")
topo = topo_lines[2].strip()
assert len(topo) == len(ce), (len(topo), len(ce))
TOPO = {"S":"signal peptide","O":"extracellular","M":"transmembrane","I":"cytoplasmic"}

# per-parent residues from the curated table
import csv
par = {}
with open(f"{ST}/sid2_variants_cendr.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        par[int(row["residue"])] = row

rows = []
cepos = cbpos = 0
for x, y in zip(A, B):
    if x != "-": cepos += 1
    if y != "-": cbpos += 1
    if x == "-":  continue          # insertion in Cb, no Ce residue to key on
    p = par.get(cepos, {})
    rows.append(dict(
        ce_pos=cepos, ce_aa=x, cb_pos=(cbpos if y != "-" else ""), cb_aa=y,
        state=("gap" if y == "-" else "identical" if x == y else "different"),
        topology=TOPO.get(topo[cepos-1], topo[cepos-1]),
        n2_aa=p.get("n2_aa",""), xz1516_aa=p.get("xz1516_aa",""),
        ju1793_aa=p.get("ju1793_aa",""), ju2466_aa=p.get("ju2466_aa",""),
        xz1516_differs_from_n2=("TRUE" if p and p.get("n2_aa") != p.get("xz1516_aa") else ""),
    ))

with open(f"{OUT}/sid2_species_alignment.tsv","w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)

# population missense from the VCF
q = ["bcftools","query","-r",SID2_REGION,
     "-f","%POS\t%REF\t%ALT\t%AF\t%AC\t%AN\t%INFO/BCSQ\n", VCF]
out = subprocess.run(q, capture_output=True, text=True)
if out.returncode: sys.exit("bcftools failed: " + out.stderr[:300])
mis = []
for line in out.stdout.strip().split("\n"):
    if not line or "missense" not in line: continue
    pos, ref, alt, af, ac, an, bcsq = line.split("\t")
    entry = next((e for e in bcsq.split(",") if "missense" in e), None)
    if entry is None: continue
    f = entry.split("|")
    aa = f[5] if len(f) > 5 else ""
    import re
    m = re.match(r"(\d+)([A-Z])>\d+([A-Z])", aa)
    if not m: continue
    resid = int(m.group(1))
    mis.append(dict(pos=pos, ref=ref, alt=alt, resid=resid,
                    from_aa=m.group(2), to_aa=m.group(3), aa_change=aa,
                    af=f"{float(af):.5f}", ac=ac, an=an,
                    topology=TOPO.get(topo[resid-1], topo[resid-1]),
                    cb_aa=next((r["cb_aa"] for r in rows if r["ce_pos"]==resid), "")))
with open(f"{OUT}/sid2_population_missense.tsv","w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(mis[0]), delimiter="\t")
    w.writeheader(); w.writerows(sorted(mis, key=lambda r: r["resid"]))

ident = sum(1 for r in rows if r["state"]=="identical")
cov   = sum(1 for r in rows if r["state"]!="gap")
print(f"Ce {len(ce)} aa, Cb {len(cb)} aa; {ident}/{cov} aligned pairs identical "
      f"= {100*ident/cov:.1f}%")
print(f"population missense variants staged: {len(mis)}")
print("key positions (Ce -> Cb):")
for p, lab in [(32,"H32 uptake His"),(34,"D34 qt13"),(94,"N94 sequon"),
               (96,"T96 focal"),(168,"H168 uptake His"),(175,"H175 uptake His")]:
    r = rows[p-1]
    print(f"  {p:3d} {r['ce_aa']} -> {r['cb_aa']}  {r['state']:9s}  {lab}")
