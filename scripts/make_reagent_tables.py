#!/usr/bin/env python3
"""Split the strain/oligo stock sheet into two clean reagent tables, and assign
the Kruglyak-series designations.

    python3 scripts/make_reagent_tables.py

    reads   supplemental_data/Strains_Oligos.csv      (the exported stock sheet)
    writes  supplemental_data/Table_S_strains.csv
            supplemental_data/Table_S_oligos.csv
            supplemental_data/reagent_table_issues.txt

WHY A SCRIPT. The source sheet is three tables stacked in one CSV with 26
trailing empty columns, and the designations have to be assigned in a fixed,
reproducible order or the mapping from wSZ to QX cannot be checked later. The
wSZ -> QX/allele map lives in DESIGNATIONS below and is the single source of
truth; rerunning reproduces the tables byte for byte.

THE SERIES, as supplied: QX2545 is the next free strain number, qq212 the next
free allele, qqIR59 the next free introgression. Assignment is by ascending wSZ
number, which is also roughly construction order, so the mapping stays legible.
  QX2545-QX2578   34 strains
  qqIR59-qqIR80   22 introgressions, one per NIL
  qq212-qq222     11 alleles over 12 edited strains -- wSZ203 and wSZ204 share
                  qq218, see ISSUES. qq223 is left free so that splitting them
                  later renumbers nothing.

WHAT AN ALLELE NAME MEANS HERE. An allele names a molecular lesion, not a
strain, so one name per independent repair event. Independent F1s from the same
injection are separate events and get separate names; clonal siblings share a
name. Four of the twelve edited strains carry a lesion that is NOT the intended
substitution -- two frameshifts, two deletions -- and are named as the indels
they are rather than as sid-2(T96K), because calling them T96K in a table would
be wrong.

INTERVALS. Six introgressions have WGS-confirmed breakpoints and are taken from
data/nil_ranges.bed. The rest are bounded by the genotyping markers recorded in
the sheet, and `breakpoint_source` says which is which. Do not quote a
marker-bounded interval as if it were sequenced.
"""

import csv
import os
import re
import sys

ISSUES = """Things the stock sheet leaves unresolved, found while building these tables.
Written by scripts/make_reagent_tables.py. Each one needs a person, not a script.

1. wSZ203 AND wSZ204 ARE RECORDED AS THE SAME ISOLATE.
   Both carry isolate "9_6" and wSZ204's note is "replicate from above ^^". If
   they are clonal siblings they share one repair event and therefore ONE allele,
   which is how the tables name them (both qq218). The manuscript says otherwise:
   METHODS.txt and the Results both call them "two independently edited lines",
   and the Figure 4B pooled 4.0% is built from them as two lines (4.4%, n=273 and
   3.7%, n=295). One of the two records is wrong.
     - if they ARE independent, wSZ204 becomes qq223, which is deliberately left
       free so that the split renumbers nothing else.
     - if they are NOT, the manuscript needs correcting, and wSZ205 (isolate 7_7,
       qq219) is the genuinely independent second N2 line already in hand.

2. THE N94A STRAIN NUMBERS DISAGREE WITH METHODS.
   METHODS.txt:36 lists "wSZ207  JU1793 edited to sid-2 N94A". There is no wSZ207
   in the stock sheet. The sheet's JU1793 N94A strain is wSZ209 (isolate 2_3),
   which METHODS does not mention. wSZ207 is absent from the sheet entirely.

3. METHODS MISDESCRIBES wSZ153.
   METHODS.txt:32 groups "wSZ192-wSZ195, wSZ153" as NILs whose "introgression
   breakpoints [are] not determined". wSZ153's breakpoint IS determined --
   12,166,600 by WGS, in the sheet and in data/nil_ranges.bed. wSZ192-wSZ195 are
   correctly described. METHODS also omits wSZ159 and wSZ167, which have WGS
   breakpoints and are in nil_ranges.bed.

4. FOUR "T96K" STRAINS DO NOT CARRY T96K.
   wSZ197 complex deletion, wSZ198 frameshift, wSZ199 frameshift, wSZ202
   deletion. They are named here as the lesions they are (qq212, qq213, qq214,
   qq217) rather than as sid-2(T96K). None is used in the manuscript. If they are
   not being kept, they need no allele names at all and the series can start at
   the first clean edit.

5. A GUIDE USED IN THE SHEET IS NOT IN THE gRNA TABLE.
   wSZ178 and wSZ179 have isolate IDs "oZ171_1_5" and "oZ171_8_6", whose prefix
   points at a guide 171, but the sheet's description column says gSZ169 for all
   four of wSZ176-wSZ179 and no gSZ171 exists in the gRNA block. The tables
   record gSZ169 as the sheet says. Confirm which guide made wSZ178 and wSZ179.

6. THREE PAIRING TYPOS IN THE SOURCE, CORRECTED HERE.
   oZ330 was listed as pairing with itself ("pair oZ330"), oZ331 likewise, and
   oZ371 was labelled "left" when it is the right-hand primer of the oZ370 pair.
   Corrected to oZ330<->oZ331 and oZ370(F)/oZ371(R).

7. THE oZ370/oZ371 COLUMNS WERE SWAPPED IN THE SOURCE.
   Every other primer has coordinates under "Prime site" and product sizes under
   "Expected side". oZ370 and oZ371 have them the other way round. Un-swapped
   here. The header itself reads "Expected side", which is presumably "size".

8. gSZ177 HAS A COORDINATE AT THE WRONG PRECISION.
   The sheet gives "chrIII:13.695" where the other guides carry base positions.
   Carried as 13,695,000 and flagged approximate. It matters because that cut
   defines the right edge of the 37 kb interval the paper's claim rests on.

9. SEQUENCE CASE IS INCONSISTENT AND MAY BE MEANINGFUL.
   oZ446 carries a lower-case "a" at the edited base and oZ456 an upper-case "C";
   gSZ169 is entirely lower case and gSZ177 begins with a capital G that is
   probably an appended transcription start rather than genomic sequence.
   Sequences are carried through verbatim. Confirm whether the leading G of
   gSZ177 is part of the target site before anyone orders from this table.

10. ONLY THE N94A TEMPLATES CARRY A PAM BLOCK.
    All four templates are 169 bp on one backbone, checked rather than assumed,
    and the script fails if any of this stops holding. In the +2 reading frame
    the backbone translates through the N94-C95-T96 sequon the paper names.
      oZ446/oZ456 differ from each other at exactly one base, 76, the
      III:13,680,248 site -- A in oZ446, C in oZ456, matching JU1793 C >
      JU2466 A. Neither carries any other change, so nothing stops Cas9
      re-cutting a correctly repaired allele.
      oZ490/oZ491 change bases 69-70 AAC>GCC, which is N94A, and base 83
      ACC>ACA, which is synonymous and is the PAM block. Each leaves base 76 at
      its own background's residue 96 -- C (96T) in oZ490, A (96K) in oZ491 --
      so they remove the glycosylation site without touching 96.
    That asymmetry is worth stating in METHODS.txt:42: the substitution-only
    templates had no block and four of twelve recovered lines came back as
    frameshifts or deletions; the N94A templates had one.

11. THREE STRAINS ARE LABELLED AS AN EDIT THEIR BACKGROUND ALREADY CARRIES.
    wSZ197, wSZ198 and wSZ199 are recorded as "JU2466 SID-2[T96K]", but JU2466
    already carries 96K, so T96K in that background is not a possible edit. The
    likely reading is that they are failed siblings of wSZ206, the JU2466 K96T
    line, which would make their repair template oZ456. That is a guess, so
    repair_template is left blank for those three rather than filled in.

12. SEVENTEEN OF THE 22 INTROGRESSIONS ARE MARKER-BOUNDED, NOT SEQUENCED.
    Only wSZ153, wSZ159, wSZ167, wSZ176, wSZ191 and wSZ196 have WGS breakpoints.
    The rest are bounded by the flanking genotyping markers, and interval_start /
    interval_end for those are marker positions, not breakpoints. The column
    breakpoint_source carries this. Three pairs are also indistinguishable at the
    markers used: wSZ158/wSZ160/wSZ161, and wSZ170/wSZ171.
"""

SRC   = "supplemental_data/Strains_Oligos.csv"
OUT_S = "supplemental_data/Table_S_strains.csv"
OUT_O = "supplemental_data/Table_S_oligos.csv"
OUT_I = "supplemental_data/reagent_table_issues.txt"

CHRIII_END = 13783801          # data/nil_ranges.bed
MARKERS = {                    # genotyping marker -> chrIII position
    "318-319": 11832608, "326-327": 12183696,
    "330-331": 12830894, "334-335": 13276309, "370-371": 13702324,
}
CUTS = {"gSZ169": 13657684, "gSZ177": 13695000, "gSZ179": 13736877,
        "gSZ182": 13680248}

# wSZ, QX, allele, class, background, donor, start, end, bp_source, derived_from,
# method, guide, lesion, verification
DESIGNATIONS = [
 ("wSZ153","QX2545","qqIR59","NIL","JU1793","JU2466",12166600,CHRIII_END,"WGS","JU1793 x JU2466","F2 recombinant, 6x backcross","", "",  "markers 326-327, 334-335 (JU2466)"),
 ("wSZ154","QX2546","qqIR60","NIL","JU1793","JU2466",11832608,CHRIII_END,"marker","JU1793 x JU2466","F2 recombinant, 6x backcross","", "", "markers 326-327, 334-335 (JU2466); not sequenced"),
 ("wSZ155","QX2547","qqIR61","NIL","JU1793","JU2466",12183696,CHRIII_END,"marker","wSZ154","recombinant arising during homozygosing","", "", "markers 326-327 (JU1793), 334-335 (JU2466); not sequenced"),
 ("wSZ158","QX2548","qqIR62","NIL","JU2466","JU1793",12183696,CHRIII_END,"marker","JU1793 x JU2466","F2 recombinant (B4 series)","", "", "markers 330-331, 334-335 (JU1793)"),
 ("wSZ159","QX2549","qqIR63","NIL","JU2466","JU1793",11945240,CHRIII_END,"WGS","JU1793 x JU2466","F2 recombinant (B4 series)","", "", "markers 326-327, 330-331, 334-335 (JU1793)"),
 ("wSZ160","QX2550","qqIR64","NIL","JU2466","JU1793",12183696,CHRIII_END,"marker","JU1793 x JU2466","F2 recombinant (B4 series)","", "", "markers 330-331, 334-335 (JU1793); genotype identical to wSZ158, wSZ161"),
 ("wSZ161","QX2551","qqIR65","NIL","JU2466","JU1793",12183696,CHRIII_END,"marker","JU1793 x JU2466","F2 recombinant (B4 series)","", "", "markers 330-331, 334-335 (JU1793); genotype identical to wSZ158, wSZ160"),
 ("wSZ167","QX2552","qqIR66","NIL","JU1793","JU2466",12560200,CHRIII_END,"WGS","wSZ153 x JU1793","F2 recombinant from 96-well screen","", "", "markers 330-331, 334-335 (JU2466)"),
 ("wSZ168","QX2553","qqIR67","NIL","JU1793","JU2466",12183696,CHRIII_END,"marker","wSZ153 x JU1793","F2 recombinant from 96-well screen","", "", "markers 330-331, 334-335 (JU2466); not sequenced"),
 ("wSZ169","QX2554","qqIR68","NIL","JU1793","JU2466",11832608,13276309,"marker","wSZ153 x JU1793","F2 recombinant from 96-well screen","", "", "internal segment: 326-327, 330-331 (JU2466), 334-335 (JU1793)"),
 ("wSZ170","QX2555","qqIR69","NIL","JU1793","JU2466",11832608,12830894,"marker","wSZ153 x JU1793","F2 recombinant from 96-well screen","", "", "internal segment: 326-327 (JU2466) only; genotype identical to wSZ171"),
 ("wSZ171","QX2556","qqIR70","NIL","JU1793","JU2466",11832608,12830894,"marker","wSZ153 x JU1793","F2 recombinant from 96-well screen","", "", "internal segment: 326-327 (JU2466) only; genotype identical to wSZ170"),
 ("wSZ176","QX2557","qqIR71","NIL","JU1793","JU2466",13657700,CHRIII_END,"WGS","wSZ167","Cas9-induced recombination","gSZ169","", "JU1793 left of cut (330-331), JU2466 right (370-371)"),
 ("wSZ177","QX2558","qqIR72","NIL","JU2466","JU1793",13657684,CHRIII_END,"marker","wSZ167","Cas9-induced recombination","gSZ169","", "reciprocal of wSZ176: JU2466 left of cut, JU1793 right; not sequenced"),
 ("wSZ178","QX2559","qqIR73","NIL","JU1793","JU2466",13657684,CHRIII_END,"marker","wSZ167","Cas9-induced recombination","gSZ169","", "as wSZ176; not sequenced"),
 ("wSZ179","QX2560","qqIR74","NIL","JU2466","JU1793",13657684,CHRIII_END,"marker","wSZ167","Cas9-induced recombination","gSZ169","", "as wSZ177; not sequenced"),
 ("wSZ191","QX2561","qqIR75","NIL","JU1793","JU2466",13657700,13695000,"WGS","wSZ176","Cas9-induced recombination","gSZ177","", "JU2466 left of cut (440-441/DraI), JU1793 right (370-371)"),
 ("wSZ192","QX2562","qqIR76","NIL","JU1793","JU2466",13657700,13695000,"marker","wSZ176","Cas9-induced recombination","gSZ177","", "as wSZ191; not sequenced"),
 ("wSZ193","QX2563","qqIR77","NIL","JU1793","JU2466",13657700,13736877,"marker","wSZ176","Cas9-induced recombination","gSZ179","", "JU2466 left of cut (370-371), JU1793 right (434-435/HpyAV); not sequenced"),
 ("wSZ194","QX2564","qqIR78","NIL","JU1793","JU2466",13657700,13736877,"marker","wSZ176","Cas9-induced recombination","gSZ179","", "as wSZ193; not sequenced"),
 ("wSZ195","QX2565","qqIR79","NIL","JU1793","JU2466",13736877,CHRIII_END,"marker","wSZ176","Cas9-induced recombination","gSZ179","", "JU1793 left of cut, JU2466 right; not sequenced"),
 ("wSZ196","QX2566","qqIR80","NIL","JU1793","JU2466",13695000,CHRIII_END,"WGS","wSZ176","Cas9-induced recombination","gSZ177","", "JU1793 left of cut (440-441/DraI), JU2466 right (370-371)"),
 ("wSZ197","QX2567","qq212","edit","JU2466","",0,0,"","JU2466","CRISPR-Cas9 HDR","gSZ182","complex deletion","Sanger: complex deletion, not the intended substitution"),
 ("wSZ198","QX2568","qq213","edit","JU2466","",0,0,"","JU2466","CRISPR-Cas9 HDR","gSZ182","frameshift","Sanger verified edit, but frameshifted"),
 ("wSZ199","QX2569","qq214","edit","JU2466","",0,0,"","JU2466","CRISPR-Cas9 HDR","gSZ182","frameshift","Sanger verified edit, but frameshifted"),
 ("wSZ200","QX2570","qq215","edit","JU1793","",0,0,"","JU1793","CRISPR-Cas9 HDR","gSZ182","T96K","Sanger verified"),
 ("wSZ201","QX2571","qq216","edit","JU1793","",0,0,"","JU1793","CRISPR-Cas9 HDR","gSZ182","T96K","Sanger verified"),
 ("wSZ202","QX2572","qq217","edit","JU1793","",0,0,"","JU1793","CRISPR-Cas9 HDR","gSZ182","deletion","Sanger: deletion present, not the intended substitution"),
 ("wSZ203","QX2573","qq218","edit","N2","",0,0,"","N2","CRISPR-Cas9 HDR","gSZ182","T96K","Sanger verified"),
 ("wSZ204","QX2574","qq218","edit","N2","",0,0,"","N2","CRISPR-Cas9 HDR","gSZ182","T96K","Sanger verified; sheet records the same isolate 9_6 as wSZ203 -- see ISSUES"),
 ("wSZ205","QX2575","qq219","edit","N2","",0,0,"","N2","CRISPR-Cas9 HDR","gSZ182","T96K","Sanger verified"),
 ("wSZ206","QX2576","qq220","edit","JU2466","",0,0,"","JU2466","CRISPR-Cas9 HDR","gSZ182","K96T","Sanger verified, no frameshift; non-clumping"),
 ("wSZ208","QX2577","qq221","edit","JU2466","",0,0,"","JU2466","CRISPR-Cas9 HDR","gSZ182","N94A","Sanger verified (PAM breaker)"),
 ("wSZ209","QX2578","qq222","edit","JU1793","",0,0,"","JU1793","CRISPR-Cas9 HDR","gSZ182","N94A","Sanger verified (PAM breaker)"),
]

# oligo, type, orientation, pair, target, product, assay, purpose
OLIGOS = [
 ("oZ318","genotyping primer","F","oZ319","III:11,832,608-11,832,695","JU1793 464 bp; JU2466 377 bp","length polymorphism","chrIII 11.83 Mb marker"),
 ("oZ319","genotyping primer","R","oZ318","III:11,832,608-11,832,695","JU1793 464 bp; JU2466 377 bp","length polymorphism","chrIII 11.83 Mb marker"),
 ("oZ326","genotyping primer","F","oZ327","III:12,183,696-12,183,821","JU1793 535 bp; JU2466 660 bp","length polymorphism","chrIII 12.18 Mb marker"),
 ("oZ327","genotyping primer","R","oZ326","III:12,183,696-12,183,821","JU1793 535 bp; JU2466 660 bp","length polymorphism","chrIII 12.18 Mb marker"),
 ("oZ330","genotyping primer","F","oZ331","III:12,830,894-12,831,147","JU1793 947 bp; JU2466 1177 bp","length polymorphism","chrIII 12.83 Mb marker"),
 ("oZ331","genotyping primer","R","oZ330","III:12,830,894-12,831,147","JU1793 947 bp; JU2466 1177 bp","length polymorphism","chrIII 12.83 Mb marker"),
 ("oZ334","genotyping primer","F","oZ335","III:13,276,309-13,276,443","JU1793 590 bp; JU2466 734 bp","length polymorphism","chrIII 13.28 Mb marker"),
 ("oZ335","genotyping primer","R","oZ334","III:13,276,309-13,276,443","JU1793 590 bp; JU2466 734 bp","length polymorphism","chrIII 13.28 Mb marker"),
 ("oZ370","genotyping primer","F","oZ371","III:13,702,324-13,702,495","JU1793 835 bp; JU2466 664 bp","length polymorphism","chrIII 13.70 Mb marker"),
 ("oZ371","genotyping primer","R","oZ370","III:13,702,324-13,702,495","JU1793 835 bp; JU2466 664 bp","length polymorphism","chrIII 13.70 Mb marker"),
 ("oZ434","genotyping primer","F","oZ435","III:13,773,217 G>A (JU1793)","","HpyAV RFLP","chrIII 13.77 Mb marker"),
 ("oZ435","genotyping primer","R","oZ434","III:13,773,217 G>A (JU1793)","","HpyAV RFLP","chrIII 13.77 Mb marker"),
 ("oZ440","genotyping primer","F","oZ441","III:13,686,395 A>G (JU1793)","","DraI RFLP","chrIII 13.69 Mb marker"),
 ("oZ441","genotyping primer","R","oZ440","III:13,686,395 A>G (JU1793)","","DraI RFLP","chrIII 13.69 Mb marker"),
 ("oZ461","genotyping primer","F","oZ468","III:13,680,248 (sid-2 edit site)","JU1793 440/140/78 bp; JU2466 440/218 bp","HpyCH4IV RFLP","screens the sid-2 96 edits"),
 ("oZ468","genotyping primer","R","oZ461","III:13,680,248 (sid-2 edit site)","JU1793 440/140/78 bp; JU2466 440/218 bp","HpyCH4IV RFLP","screens the sid-2 96 edits"),
 ("oZ446","repair template","","","III:13,680,248","","","installs the JU2466 allele (96K) on a 96T background; single change, no PAM block"),
 ("oZ456","repair template","","","III:13,680,248","","","installs the JU1793 allele (96T) on a 96K background; single change, no PAM block"),
 ("oZ490","repair template","","","III:13,680,231-13,680,248","","","N94A in JU1793, keeping 96T (A-x-T sequon); carries a synonymous PAM block"),
 ("oZ491","repair template","","","III:13,680,231-13,680,248","","","N94A in JU2466, keeping 96K (A-x-K sequon); carries a synonymous PAM block"),
 ("gSZ169","gRNA","","","III:13,657,684","","","Cas9-induced recombination, wSZ167 -> wSZ176 series"),
 ("gSZ177","gRNA","","","III:13,695,000 (approximate)","","","Cas9-induced recombination, wSZ176 -> wSZ191/wSZ192/wSZ196"),
 ("gSZ179","gRNA","","","III:13,736,877","","","Cas9-induced recombination, wSZ176 -> wSZ193/wSZ194/wSZ195"),
 ("gSZ182","gRNA","","","III:13,680,248","","","cuts at the sid-2 residue-96 site for all point edits"),
]


# Supplied after the stock sheet was exported, so they have no row in it. Both
# are 169 bp and share the oZ446/oZ456 backbone. Changes, in the +2 frame:
# bases 69-70 AAC>GCC is N94A, base 83 ACC>ACA is synonymous and is the PAM
# block, and base 76 is left at whichever residue-96 allele the background
# already carries -- 96T for JU1793, 96K for JU2466. That is what makes them a
# clean sequon test: they remove the glycosylation site without touching 96.
EXTRA_SEQ = {
 "oZ490": ("ACGGAACTGCCGCAATTTCGGACCTTAAAAATGTGACATTTATATTGGAGGTCACAACTGACA"
           "CTAAAgcCTGCACGTTTACaGCTAATTACACCGGATACTTCACTCCGGATCCCAAGAGCAAGC"
           "CATTTCAGTTAGGATTCGCAAGTGCCACGTTGAACCGAGATAT"),
 "oZ491": ("ACGGAACTGCCGCAATTTCGGACCTTAAAAATGTGACATTTATATTGGAGGTCACAACTGACA"
           "CTAAAgcCTGCAaGTTTACaGCTAATTACACCGGATACTTCACTCCGGATCCCAAGAGCAAGC"
           "CATTTCAGTTAGGATTCGCAAGTGCCACGTTGAACCGAGATAT"),
}

# which repair template built which edited strain. Blank where the sheet does
# not say and it cannot be inferred safely -- see ISSUES.
REPAIR = {
 "wSZ200": "oZ446", "wSZ201": "oZ446", "wSZ202": "oZ446",
 "wSZ203": "oZ446", "wSZ204": "oZ446", "wSZ205": "oZ446",
 "wSZ206": "oZ456",
 "wSZ208": "oZ491", "wSZ209": "oZ490",
}


def read_source(path):
    """Return {block: [row, ...]} with the trailing empty columns stripped."""
    blocks, cur = {}, None
    with open(path, newline="") as fh:
        for row in csv.reader(fh):
            while row and row[-1] == "":
                row.pop()
            if not row:
                continue
            if row[0] in ("STRAINS", "OLIGOS", "gRNAs"):
                cur = row[0]
                blocks[cur] = []
                continue
            if cur:
                blocks[cur].append(row)
    return blocks


def genotype(d):
    (wsz, qx, allele, cls, bg, donor, start, end, src,
     parent, method, guide, lesion, verif) = d
    if cls == "NIL":
        # no thousands separators: this is a CSV field and a comma here is a
        # trap for anything that splits on commas without a real parser
        return f"{allele}[{donor}>{bg} III:{start}-{end}]"
    # sid-2(qq215[T96K]) is the form for a designed substitution; an unintended
    # indel is just sid-2(qq212), because the bracket asserts what was built
    if lesion in ("T96K", "K96T", "N94A"):
        return f"sid-2({allele}[{lesion}]) {bg}"
    return f"sid-2({allele}) {bg}"


def main():
    if not os.path.exists(SRC):
        sys.exit(f"missing {SRC}")
    blocks = read_source(SRC)
    src_strains = {r[0] for r in blocks.get("STRAINS", [])}
    src_oligos = {r[0] for r in blocks.get("OLIGOS", [])} | \
                 {r[0] for r in blocks.get("gRNAs", [])}
    seq = {r[0]: r[1] for r in blocks.get("OLIGOS", []) + blocks.get("gRNAs", [])}
    seq.update(EXTRA_SEQ)
    isolate = {r[0]: r[1] for r in blocks.get("STRAINS", [])}

    # every source row must be carried over, and nothing invented
    named = {d[0] for d in DESIGNATIONS}
    assert named == src_strains, (
        f"strain mismatch: only in sheet {sorted(src_strains - named)}, "
        f"only in table {sorted(named - src_strains)}")
    named_o = {o[0] for o in OLIGOS}
    assert named_o == src_oligos | set(EXTRA_SEQ), (
        f"oligo mismatch: only in sheet {sorted(src_oligos - named_o)}, "
        f"only in table {sorted(named_o - src_oligos - set(EXTRA_SEQ))}")
    assert all(seq.get(o[0]) for o in OLIGOS), "an oligo has no sequence"

    with open(OUT_S, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["strain", "lab_id", "isolate", "class", "allele", "genotype",
                    "background", "donor", "chrom", "interval_start",
                    "interval_end", "interval_kb", "breakpoint_source",
                    "derived_from", "construction", "guide_rna",
                    "repair_template", "lesion", "verification"])
        for d in DESIGNATIONS:
            (wsz, qx, allele, cls, bg, donor, start, end, src_,
             parent, method, guide, lesion, verif) = d
            w.writerow([qx, wsz, isolate.get(wsz, ""), cls, allele, genotype(d),
                        bg, donor, "III" if cls == "NIL" else "",
                        start or "", end or "",
                        f"{(end - start) / 1000:.1f}" if cls == "NIL" else "",
                        src_, parent, method, guide, REPAIR.get(wsz, ""),
                        lesion, verif])

    with open(OUT_O, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["oligo", "type", "sequence", "orientation", "pairs_with",
                    "target", "expected_product", "assay", "purpose"])
        for o in OLIGOS:
            name, typ, orient, pair, target, product, assay, purpose = o
            w.writerow([name, typ, seq.get(name, ""), orient, pair, target,
                        product, assay, purpose])

    # The two repair templates are the whole sid-2 claim in reagent form, so
    # check them rather than trusting the labels: same length, differing at
    # exactly one base, in the direction the paper states for III:13,680,248
    # (JU1793 C > JU2466 A).
    a, b = seq["oZ446"], seq["oZ456"]
    assert len(a) == len(b) == 169, (len(a), len(b))
    d = [i for i, (x, y) in enumerate(zip(a, b)) if x.upper() != y.upper()]
    assert d == [75], f"repair templates differ at {d}, expected one site"
    assert a[75].upper() == "A" and b[75].upper() == "C", (a[75], b[75])
    # and the N94A pair: both knock out the sequon at bases 69-70 (AAC>GCC),
    # both carry the synonymous block at base 83, and each keeps its own
    # background's residue 96 -- C (96T) for JU1793, A (96K) for JU2466
    for name, r96 in (("oZ490", "C"), ("oZ491", "A")):
        t = seq[name]
        assert len(t) == 169, (name, len(t))
        assert t[68:70].upper() == "GC" and b[68:70].upper() == "AA", name
        assert t[82].upper() == "A" and b[82].upper() == "C", name
        assert t[75].upper() == r96, (name, t[75])

    with open(OUT_I, "w") as fh:
        fh.write(ISSUES)

    n_nil = sum(1 for d in DESIGNATIONS if d[3] == "NIL")
    n_ed = sum(1 for d in DESIGNATIONS if d[3] == "edit")
    n_al = len({d[2] for d in DESIGNATIONS if d[3] == "edit"})
    print(f"{OUT_S}: {len(DESIGNATIONS)} strains "
          f"({n_nil} NILs, {n_ed} edited), "
          f"{DESIGNATIONS[0][1]}-{DESIGNATIONS[-1][1]}")
    print(f"  introgressions: {n_nil} ({DESIGNATIONS[0][2]}-"
          f"{[d[2] for d in DESIGNATIONS if d[3]=='NIL'][-1]})")
    print(f"  alleles: {n_al} over {n_ed} edited strains")
    print(f"{OUT_O}: {len(OLIGOS)} oligos "
          f"({sum(1 for o in OLIGOS if o[1]=='genotyping primer')} primers, "
          f"{sum(1 for o in OLIGOS if o[1]=='repair template')} repair templates, "
          f"{sum(1 for o in OLIGOS if o[1]=='gRNA')} gRNAs)")
    print(f"{OUT_I}: written")


if __name__ == "__main__":
    main()
