#!/usr/bin/env python3
"""Stage the SID-2 ortholog survey across nematodes.

  python3 scripts/make_sid2_ortholog_tables.py
    -> supplemental_data/structure/sid2_ortholog_search.tsv
       supplemental_data/structure/sid2_ortholog_alignment.tsv
       supplemental_data/structure/sid2_ortholog_conservation.tsv
       supplemental_data/structure/sid2_ortholog_sequences.fa

WHY THIS EXISTS
sid-2 96K is a polymorphism inside C. elegans. Calling 96T the ancestral state
needs outgroups, and the existing two-species comparison
(make_sid2_alignment_tables.py, C. elegans against C. briggsae) is one
outgroup. This widens it to a ladder of 20 nematode proteomes and asks the
question the other way round: how far out does the comparison still work at
all?

THE SEARCH, WHICH THIS SCRIPT DOES NOT RE-RUN
Reciprocal-best-hit blastp of C. elegans SID-2 against twenty UniProt
reference proteomes, 488,718 proteins in total. Proteomes were downloaded whole
(~120 MB), one BLAST database built per species, and the search run with
NCBI blastp 2.x, `-evalue 10 -seg no -max_target_seqs 5`; a hit counts as an
ortholog at E < 1e-5. Forward best hits were confirmed by blasting the hit back
against the C. elegans proteome. Re-running that needs the proteomes, which are
far too large to deposit and are a public UniProt release rather than our data,
so the RESULT is pinned in SEARCH below and this script rebuilds everything
downstream of it from sequences fetched by accession. The proteome IDs are
recorded so the search can be repeated.

WHAT THE SEARCH FOUND, AND THE THREE PLACES THINGS BREAK DOWN
  1. The sequon breaks at the base of the genus. N94-x-[ST] is intact in all
     six Elegans-group orthologs; C. angaria, in the basal Angaria group, has
     Q-G-F at the same three columns and no sequon.
  2. The column stops being readable one step before that. C. japonica's
     ortholog is split across two proteome entries (A0A8R1DPT1 covers Ce
     167-311, A0A8R1IAX9 covers Ce 21-127), and inside the N-terminal fragment
     the alignment sits in a TTDT repeat and is offset by two to three
     residues, so the identity of the column aligned to Ce 96 is AMBIGUOUS
     there and is reported as such rather than as a substitution.
  3. Orthology itself stops being detectable immediately outside
     Caenorhabditis. No proteome outside the genus yields a hit at E < 1e-5 --
     not even Diploscapter pachys, the sister genus (best hit E = 4.1). Two
     independent resources agree: the UniRef50 cluster of G5EEV9 has exactly
     one member, and NCBI's ortholog set for sid-2 within Nematoda is empty.

SO THE ANCESTRAL-STATE CLAIM HAS A DEPTH, AND IT IS THE ELEGANS GROUP.
96T is ancestral at that depth and 96K is derived within C. elegans. It is NOT
"deeply conserved across nematodes", and this file is written so that the
figure cannot say otherwise.

WHY THE THREE-RESIDUE WINDOW IS THE STRONG FORM OF THE ARGUMENT
A single conserved column in a 42-48% identity alignment of a Thr-rich,
low-complexity region is weak evidence: Thr alone is 13.3% of the ectodomain,
so hitting one by chance is not unlikely. The window controls for that, because
all three columns come from the same alignment:

  N94  the constrained Asn of the sequon   6/6 conserved   92.5th percentile
  C95  the unconstrained X of the sequon   0/6 conserved   15.6th percentile
  T96  the constrained Ser/Thr             5/6 T, 6/6 ST   82.1st percentile

Against an ectodomain background of 2.23/6 (37.1%), with only 15.0% of
ectodomain positions conserved in all six. A bad alignment does not produce
100% / 0% / 83% across three adjacent columns.

CALIBRATION, INCLUDING THE PART THAT CUTS THE OTHER WAY. Three of the five
ectodomain sequons in C. elegans SID-2 (at 71, 81 and 94) are intact in all six
orthologs and two (at 100 and 148) are intact in one. So a fully conserved
sequon is common in this protein rather than remarkable, and the 94 sequon is
one of the conserved majority, not a standout.

NO FUNCTIONAL CLAIM IS MADE HERE. The glycosylation hypothesis for this sequon
was tested with N94A and failed -- see SUPP_FIG_XX_sid2_allele_swaps_full.R.
Conservation of the sequon is evidence about history, not about mechanism.
"""
import os, subprocess, sys, datetime, statistics

OUT = "supplemental_data/structure"
CE_ACC = "G5EEV9"
ECD = (21, 193)                 # DeepTMHMM extracellular span
E_ORTH = 1e-5                   # the significance an ortholog call required

# upid, species, group, depth rank, best-hit accession, %id, %query cov, E
SEARCH = [
 ("UP000001940", "Caenorhabditis elegans",        "Elegans group",                 0, "G5EEV9",     100.0, 100.0, 0.0),
 ("UP000008549", "Caenorhabditis briggsae",       "Elegans group",                 1, "A8XSB8",      41.9, 100.0, 4.43e-62),
 ("UP000230233", "Caenorhabditis nigoni",         "Elegans group",                 1, "A0A2G5UT38",  42.3, 100.0, 2.70e-61),
 ("UP000095282", "Caenorhabditis tropicalis",     "Elegans group",                 1, "A0A1I7TM24",  47.9,  95.0, 1.45e-87),
 ("UP000008281", "Caenorhabditis remanei",        "Elegans group",                 1, "E3LYG5",      43.1, 100.0, 8.06e-71),
 ("UP000216463", "Caenorhabditis latens",         "Elegans group",                 1, "A0A261BV85",  48.3, 100.0, 6.50e-81),
 ("UP000008068", "Caenorhabditis brenneri",       "Elegans group",                 1, "G0ML34",      42.8,  98.0, 2.22e-57),
 ("UP000005237", "Caenorhabditis japonica",       "Japonica group",                2, "A0A8R1DPT1",  64.3,  47.0, 2.55e-53),
 ("UP000494206", "Caenorhabditis bovis",          "basal Caenorhabditis",          3, "A0A8S1F7L6",  45.5,  14.0, 2.2),
 ("UP000835052", "Caenorhabditis auriculariae",   "Angaria group",                 3, "A0A8S1H7Y8",  32.6,  14.0, 0.15),
 ("UP001152747", "Caenorhabditis angaria",        "Angaria group",                 3, "A0A9P1N397",  28.6,  93.0, 7.14e-18),
 ("UP000218231", "Diploscapter pachys",           "Rhabditidae (sister genus)",    4, "A0A2A2JYI7",  24.7,  30.0, 4.1),
 ("UP001328107", "Pristionchus mayeri",           "Diplogastridae",                5, "A0AAN5CQ80",  30.0,  19.0, 0.83),
 ("UP000095283", "Heterorhabditis bacteriophora", "Heterorhabditidae",             5, "A0A1I7XUH2",  26.5,  30.0, 0.62),
 ("UP000025227", "Haemonchus contortus",          "Strongylida (clade V)",         6, "A0A7I5E9P0",  31.0,  41.0, 0.007),
 ("UP000492821", "Panagrellus redivivus",         "Panagrolaimomorpha (clade IV)", 7, "A0A7E4VDI2",  24.7,  42.0, 1.5),
 ("UP000035682", "Strongyloides ratti",           "Panagrolaimomorpha (clade IV)", 7, "A0A090L2T1",  31.1,  18.0, 1.9),
 ("UP000659654", "Bursaphelenchus xylophilus",    "Tylenchomorpha (clade IV)",     7, "A0A7I8XEM1",  26.4,  35.0, 2.1),
 ("UP000006672", "Brugia malayi",                 "Spiruromorpha (clade III)",     8, "A0A4E9FV15",  36.5,  16.0, 3.2),
 ("UP000054776", "Trichinella spiralis",          "Dorylaimia (clade I)",          9, "A0A0V1BU79",  25.6,  29.0, 0.33),
]

# The ortholog set the per-position tables are built from: one protein per
# species, covering Ce 96. C. japonica needs its N-terminal fragment, not its
# best hit, and its column is flagged ambiguous for the reason in the header.
ALIGN_SET = [
 ("Caenorhabditis briggsae",   "A8XSB8",     "Elegans group",  False),
 ("Caenorhabditis nigoni",     "A0A2G5UT38", "Elegans group",  False),
 ("Caenorhabditis tropicalis", "A0A1I7TM24", "Elegans group",  False),
 ("Caenorhabditis remanei",    "E3LYG5",     "Elegans group",  False),
 ("Caenorhabditis latens",     "A0A261BV85", "Elegans group",  False),
 ("Caenorhabditis brenneri",   "G0ML34",     "Elegans group",  False),
 ("Caenorhabditis japonica",   "A0A8R1IAX9", "Japonica group", True),
 ("Caenorhabditis angaria",    "A0A9P1N397", "Angaria group",  False),
]
# the group the conservation statistics are computed over: full-length
# orthologs with an unambiguous column, which is the Elegans group
CONS_GROUP = [sp for sp, _, grp, amb in ALIGN_SET if grp == "Elegans group" and not amb]


def fetch(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    r = subprocess.run(["curl", "-sS", "--max-time", "60", url],
                       capture_output=True, text=True)
    if r.returncode or not r.stdout.startswith(">"):
        sys.exit(f"could not fetch {acc}: {r.stderr.strip()[:200]}")
    return "".join(r.stdout.split("\n")[1:]).strip()


def ce_anchored(ce, other):
    """Global align and return one residue per C. elegans position.

    Needleman-Wunsch, BLOSUM62, gap open -11 extend -1: the same settings
    make_sid2_alignment_tables.py uses, so the two files agree on method.
    """
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score, al.mode = -11, -1, "global"
    a = al.align(ce, other)[0]
    s1, s2 = str(a[0]), str(a[1])
    return [c2 for c1, c2 in zip(s1, s2) if c1 != "-"]


def main():
    os.makedirs(OUT, exist_ok=True)
    ce = fetch(CE_ACC)
    assert len(ce) == 311, len(ce)
    assert ce[93] == "N" and ce[95] == "T", (ce[93], ce[95])

    seqs = {"Caenorhabditis elegans": (CE_ACC, ce)}
    for sp, acc, _grp, _amb in ALIGN_SET:
        seqs[sp] = (acc, fetch(acc))
        print(f"  fetched {acc:12s} {sp:30s} {len(seqs[sp][1]):4d} aa")

    with open(f"{OUT}/sid2_ortholog_sequences.fa", "w") as fh:
        for sp, (acc, s) in seqs.items():
            fh.write(f">{acc} {sp.replace(' ', '_')} len={len(s)} "
                     f"fetched={datetime.date.today()}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    # ---- the search table -------------------------------------------------
    with open(f"{OUT}/sid2_ortholog_search.tsv", "w") as fh:
        fh.write("proteome\tspecies\tgroup\tdepth_rank\taccession\tpercent_identity"
                 "\tquery_coverage\tevalue\tis_ortholog\n")
        for upid, sp, grp, rank, acc, pid, cov, ev in SEARCH:
            fh.write(f"{upid}\t{sp}\t{grp}\t{rank}\t{acc}\t{pid}\t{cov}\t{ev:.3g}"
                     f"\t{'TRUE' if ev < E_ORTH else 'FALSE'}\n")
    n_orth = sum(1 for *_, ev in [(r[-1],) for r in SEARCH] if ev < E_ORTH)
    n_orth = sum(1 for r in SEARCH if r[7] < E_ORTH)
    print(f"  search table: {len(SEARCH)} proteomes, {n_orth} with an ortholog "
          f"at E < {E_ORTH:g}")

    # ---- the per-position alignment --------------------------------------
    mapped = {sp: ce_anchored(ce, seqs[sp][1]) for sp, _, _, _ in ALIGN_SET}
    amb = {sp: a for sp, _, _, a in ALIGN_SET}
    grp = {sp: g for sp, _, g, _ in ALIGN_SET}
    with open(f"{OUT}/sid2_ortholog_alignment.tsv", "w") as fh:
        fh.write("ce_pos\tce_aa\tspecies\tgroup\taligned_aa\tidentical\tambiguous\n")
        for i in range(len(ce)):
            for sp in mapped:
                aa = mapped[sp][i]
                fh.write(f"{i+1}\t{ce[i]}\t{sp}\t{grp[sp]}\t{aa}\t"
                         f"{'TRUE' if aa == ce[i] else 'FALSE'}\t"
                         f"{'TRUE' if amb[sp] else 'FALSE'}\n")

    # ---- per-position conservation over the Elegans group ----------------
    n = len(CONS_GROUP)
    cons = [sum(1 for sp in CONS_GROUP if mapped[sp][i] == ce[i])
            for i in range(len(ce))]
    ecd_i = [i for i in range(len(ce)) if ECD[0] <= i + 1 <= ECD[1]]

    def pctile(i):
        below = sum(1 for j in ecd_i if cons[j] < cons[i])
        tied = sum(1 for j in ecd_i if cons[j] == cons[i])
        return 100.0 * (below + 0.5 * tied) / len(ecd_i)

    with open(f"{OUT}/sid2_ortholog_conservation.tsv", "w") as fh:
        fh.write("ce_pos\tce_aa\tn_conserved\tn_orthologs\tfraction"
                 "\tin_ectodomain\tecd_percentile\n")
        for i in range(len(ce)):
            ine = ECD[0] <= i + 1 <= ECD[1]
            fh.write(f"{i+1}\t{ce[i]}\t{cons[i]}\t{n}\t{cons[i]/n:.4f}\t"
                     f"{'TRUE' if ine else 'FALSE'}\t"
                     f"{pctile(i):.1f}\n" if ine else
                     f"{i+1}\t{ce[i]}\t{cons[i]}\t{n}\t{cons[i]/n:.4f}\tFALSE\t\n")

    # ---- the assertions that keep the figure honest ----------------------
    mean_cons = statistics.mean(cons[i] for i in ecd_i)
    all_cons = sum(1 for i in ecd_i if cons[i] == n)
    print(f"  conservation over {n} Elegans-group orthologs: "
          f"ectodomain mean {mean_cons:.2f}/{n} ({100*mean_cons/n:.1f}%), "
          f"{all_cons}/{len(ecd_i)} positions conserved in all "
          f"({100*all_cons/len(ecd_i):.1f}%)")
    for p in (94, 95, 96):
        i = p - 1
        print(f"  Ce {ce[i]}{p}: {cons[i]}/{n} conserved, "
              f"{pctile(i):.1f}th ectodomain percentile, "
              f"residues {''.join(mapped[sp][i] for sp in CONS_GROUP)}")

    assert cons[93] == n, "N94 is no longer conserved in every Elegans-group ortholog"
    assert cons[94] == 0, "C95 is no longer the freely varying position"
    assert cons[95] == n - 1, "T96 conservation has changed"
    assert all(mapped[sp][95] in "ST" for sp in CONS_GROUP), \
        "position 96 is no longer Ser or Thr in every Elegans-group ortholog"
    ang = mapped["Caenorhabditis angaria"]
    assert not (ang[93] == "N" and ang[95] in "ST"), \
        "C. angaria now carries the sequon, which is the claimed breakdown point"
    assert all(r[7] >= E_ORTH for r in SEARCH if r[3] >= 4), \
        "a proteome outside Caenorhabditis now yields a significant hit"
    seq_intact = sum(1 for sp in CONS_GROUP
                     if mapped[sp][93] == "N" and mapped[sp][95] in "ST"
                     and mapped[sp][94] != "P")
    assert seq_intact == n, "the N94 sequon is no longer intact in every ortholog"
    print(f"  N94-x-[ST] intact in {seq_intact}/{n}; C. angaria carries "
          f"{ang[93]}-{ang[94]}-{ang[95]} and no sequon")
    print("  all assertions passed")


if __name__ == "__main__":
    main()
