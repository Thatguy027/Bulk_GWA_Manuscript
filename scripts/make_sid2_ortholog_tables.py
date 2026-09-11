#!/usr/bin/env python3
"""Stage the SID-2 ortholog survey across nematodes.

  python3 scripts/make_sid2_ortholog_tables.py
    -> supplemental_data/structure/sid2_ortholog_search.tsv
       supplemental_data/structure/sid2_ortholog_alignment.tsv
       supplemental_data/structure/sid2_ortholog_conservation.tsv
       supplemental_data/structure/sid2_ortholog_window_survey.tsv
       supplemental_data/structure/sid2_env_rnai_sensitivity.tsv
       supplemental_data/structure/sid2_species_name_map.tsv
       supplemental_data/structure/sid2_species_tree.nwk
       supplemental_data/structure/sid2_ortholog_sequences.fa

WHY THIS EXISTS
sid-2 96K is a polymorphism inside C. elegans. Calling 96T the ancestral state
needs outgroups, and the existing two-species comparison
(make_sid2_alignment_tables.py, C. elegans against C. briggsae) is one
outgroup. This widens it to a ladder of 52 nematode proteomes -- 20 UniProt
reference proteomes spanning the phylum plus 32 from the Caenorhabditis
Genomes Project -- and asks the question the other way round: how far out does
the comparison still work at all?

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

WHAT THE SEARCH FOUND
Orthology stops being detectable immediately outside Caenorhabditis. No
proteome outside the genus yields a hit at E < 1e-5 -- not even Diploscapter
pachys, the sister genus (best hit E = 4.1). Two independent resources agree:
the UniRef50 cluster of G5EEV9 has exactly one member, and NCBI's ortholog set
for sid-2 within Nematoda is empty. So the comparison is bounded by the genus.

A FIRST PASS OVERSTATED THE BREAKDOWN INSIDE THE GENUS, and this file records
the correction. With only the UniProt species it looked as though the sequon
was intact through the Elegans group and broken at C. angaria, the most distant
species with a detectable ortholog. Two things were wrong with that. The
C. angaria call came from a global alignment whose local block identity around
the window is 26%, far too low to read a column from; and the UniProt sample
contained no Elegans-group species that lacks the sequon, which the wider
sample shows do exist. Both are fixed below: the window is now read from the
BLAST HSP that spans it, and every call carries a confidence.

THE SURVEY WAS THEN WIDENED WITH THE CAENORHABDITIS GENOMES PROJECT
UniProt carries only a handful of Caenorhabditis proteomes. The Caenorhabditis
Genomes Project v2 release (caenorhabditis.org, Zenodo 10.5281/zenodo.12633738)
adds 32 proteomes, 750,673 proteins. Only one of the 32, C. auriculariae, is a
species the UniProt set already had, so the two sets together cover 51 distinct
species in 52 proteomes; 38 of the 51 comparators yield an ortholog and every
one of them is a Caenorhabditis. Those were searched the same way, and the
residues aligned to
C. elegans 94-96 were read from the BLAST HSP that spans them rather than from
a global alignment -- a global alignment of the full 311 aa query against a
partial ortholog misplaces the window, which is how a first pass produced
several calls this one contradicts.

EVERY WINDOW CALL CARRIES ITS OWN CONFIDENCE, and calls below the floor are
counted neither way. The measure is the local identity of the +/-10 residue
block around the window. The six UniProt Elegans-group orthologs span 37-47%,
so 37% is the floor: it is the level at which the original six were already
being trusted, not a threshold chosen to get an answer. Sequon-lacking species
sit in the same 37-47% band as sequon-carrying ones, so the losses cannot be
dismissed as bad alignment.

WHAT THE WIDER SAMPLE SAYS, over the 14 species that clear the floor
  residue 96 is Ser or Thr in 12 of 14 -- 11 Thr, 1 Ser -- and LYSINE IN NONE
  the sequon is intact in 11 of 14
  of the three losses, only C. doughertyi changes position 96 itself (to Ala);
  C. afra keeps Thr96 and loses Asn94, and C. sp54 changes both

So the constraint sits on residue 96 being small and hydroxylated rather than on
the sequon as a unit, and C. afra is a NATURAL AxT -- the same construct panel D
shows to be fully resistant. Two species below the floor (C. sp44, C. sp27)
appear to carry Lys96, but at 32% and 21% block identity those calls are not
trustworthy and are reported as unreliable rather than used or hidden.

ACKNOWLEDGEMENT REQUESTED BY THE RESOURCE: we thank members of the
Caenorhabditis Genomes Project for prepublication access to genome and
transcriptome data.

ENVIRONMENTAL RNAi COMPETENCE, LAYERED ON FROM THE LITERATURE
Nuez & Felix 2012 (PLoS ONE 7:e29811) scored Caenorhabditis species for the
response to ingested dsRNA. Their species are provisional numbers, and many
have since been named -- C. sp. 11 is C. tropicalis and C. sp. 10 is
C. doughertyi -- so every call here is mapped through Felix, Braendle & Cutter
2014 (PLoS ONE 9:e94723), whose type strains are the same isolates Nuez & Felix
tested. RNAI below records the provisional designation alongside the formal
name so the mapping can be checked rather than trusted.

ONLY CALLS THE PAPER STATES IN PROSE ARE USED. Its Table S1 is a workbook of
per-species Wilcoxon tests rather than a table of calls, and reading calls out
of it would mean deciding which p-value belongs to the strain and which to the
N2 control on each sheet. The prose is unambiguous, so that is the source, with
C. tropicalis added because Table S2 exists to rescue it and its own sheet
gives p = 0.93. Species the paper tested but does not call in prose are absent
from RNAI rather than guessed at: C. doughertyi and C. nigoni among them.

THE ISOLATE USUALLY DIFFERS FROM THE SEQUENCED ONE, which matters more here
than it normally would, because the same paper reports intraspecific variation
in C. elegans. Each row carries the tested strain so the mismatch is visible.

WHAT THE OVERLAY SHOWS: residue 96 does not predict the phenotype. Of the
species with an intact sequon and a published call, four are sensitive or
weakly so and four are insensitive. And C. afra, which has lost Asn94 and is a
natural AxT, IS sensitive -- the comparative mirror of the AxT construct being
fully resistant. That agrees with the paper's own conclusion that environmental
RNAi was gained or lost repeatedly in the genus, and with the N94A result: the
sequon is not the mechanism.

SO THE ANCESTRAL-STATE CLAIM HAS A DEPTH, AND IT IS THE GENUS.
Residue 96 is Ser or Thr across Caenorhabditis and Lys in none of the species
that can be read, so 96K is derived within C. elegans. It is NOT "deeply
conserved across nematodes": one genus out, the protein is not alignable at
all. This file is written so that the figures cannot say otherwise.

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

# The Caenorhabditis Genomes Project v2 window survey. Columns: species,
# strain, best-hit protein, %identity, %query coverage, E, the residues aligned
# to C. elegans 88-104 taken from the HSP that spans 94-96, and the local
# identity of the +/-10 block around the window. None means no significant hit,
# or none whose HSP spans the window.
WIN = list(range(88, 105))
CGP = [
 ("macrosperma", "JU2083", "CMACR.g10030.t1", 52.0, 47, 2.62e-42, None, None),
 ("sp54", "BRC20483", "CSP54.g24653.t2", 50.0, 100, 5.31e-80, "--TKNIHTNETANYTGF", 47.4),
 ("wallacei", "JU1898", "CWALL.g20938.t2", 47.8, 98, 4.67e-82, "ATSHVGNGTETNNYTGV", 47.4),
 ("kamaaina", "QG2077", "CKAMA.g23423.t2", 47.7, 94, 2.75e-73, "LSKDTNNQTVAKNYTGE", 52.6),
 ("sp33", "BRC20065", "CSP33.g16743.t2", 46.8, 95, 9.16e-75, "SKVNVLSDQKSYQFMGI", 15.8),
 ("oiwi", "ECA1100", "COIWI.g11269.t1", 46.1, 94, 5.37e-71, "LSDTTTNQTVAKNYTGD", 47.4),
 ("doughertyi", "JU1771", "CDOUG.g26961.t2", 46.0, 98, 9.78e-77, "KITDYSNSAATLNYTGV", 42.1),
 ("sp51", "QG2939", "CSP51.g5930.t1", 45.8, 94, 6.06e-71, "VFNQTTTISEVFNYTGV", 31.6),
 ("sp44", "BRC20300", "CSP44.g503.t1", 45.8, 98, 1.09e-79, "ASVSSANEKKTVIYTGN", 31.6),
 ("sp48", "BRC20454", "CSP48.g15529.t1", 44.5, 94, 1.78e-59, "VNQTLSNATETLTYNGT", 42.1),
 ("imperialis", "EG5942", "CIMPE.g15758.t1", 42.9, 99, 1.96e-59, "AEVKDLNGTVTSNFTGV", 36.8),
 ("sp49", "BRC20456", "CSP49.g20575.t1", 42.3, 100, 7.09e-68, "ANVKDFNGTVVANFTGV", 36.8),
 ("nouraguensis", "JU2079", "CNOUR.g1843.t1", 42.0, 94, 4.08e-57, "ITSGE-----ATNYTGI", 31.6),
 ("afra", "JU1286", "CAFRA.g13378.t1", 42.0, 94, 1.13e-60, "VASTTVVTTVVTNYTGV", 47.4),
 ("sp25", "ZF1457", "CSP25.g14051.t1", 41.1, 70, 6.27e-43, "ANKDYYNATAVVNFTGV", 36.8),
 ("sp8", "DF5173", "CSP08.g2846.t1", 40.4, 46, 2.8e-21, None, None),
 ("sp41", "BRC20276", "CSP41.g21785.t1", 40.1, 97, 1.67e-57, "TKVKILNVSIDEQFSGV", 15.8),
 ("yunquensis", "XZ1518", "CYUNQ.g1822.t1", 39.9, 95, 1.74e-52, "ITSGE-----ATNYTGV", 31.6),
 ("sp56", "JU2215", "CSP56.g4670.t1", 39.0, 46, 2.14e-20, None, None),
 ("sp24", "QG555", "CSP24.g7786.t1", 38.4, 46, 5.78e-20, None, None),
 ("drosophilae", "DF5077", "CDROS.g7470.t1", 38.0, 48, 2.57e-24, None, None),
 ("dolens", "NIC394", "CDOLE.g6252.t1", 37.7, 46, 3.97e-19, None, None),
 ("plicata", "SB355", "CPLIC.g4090.t1", 35.6, 47, 4.53e-12, None, None),
 ("sp27", "ECA211", "CSP27.g3101.t1", 33.2, 94, 4.54e-31, "VTGKGDNLKGVWEFNET", 21.1),
 ("virilis", "JU1968", "CVIRI.g4487.t1", 33.0, 96, 1.14e-29, "L--DPLNA--VTPFSGT", 21.1),
 ("portoensis", "EG5626", "CPORT.g619.t1", 32.6, 99, 1.14e-26, "ATNDSTNVTYSYSYTGV", 42.1),
 ("sp2", "DF5070", "CSP02.g1087.t1", 32.5, 93, 8.61e-31, "STE------FQTSFHGV", 21.1),
 ("sp30", "DF5174", "CSP30.g337.t1", 31.1, 99, 6.86e-34, "VTS------IETVLTGD", 26.3),
 ("castelli", "JU1956", "CCAST.g7759.t1", 30.4, 87, 7.7e-17, "WDNHYEQGVFVVNETA-", 26.3),
 ("astrocarya", "NIC1040", "CASTR.g11858.t1", 30.3, 99, 4.5e-24, "------NLNNADPFTGF", 21.1),
 ("auriculariae", "NKZ352", None, None, None, None, None, None),
 ("monodelphis", "JU1667", None, None, None, None, None, None),
]
CGP_FLOOR = 37.0   # block identity below which a window call is not counted
# CGP genome codes, needed for the two species with no hit at all (the rest are
# read off the accession prefix). C. auriculariae is the one species present in
# BOTH sets: UP000835052 and CAURI are independent proteomes of it, and neither
# yields an ortholog, so the 52 proteomes cover 51 distinct species.
CGP_CODE = {"auriculariae": "CAURI", "monodelphis": "CMONO"}
CGP_DUP = {"auriculariae"}          # also present as a UniProt proteome


def cgp_species(short):
    """CGP short name -> the species string used in the deposited tables."""
    return "Caenorhabditis " + short

# Response to ingested dsRNA, from Nuez & Felix 2012 (PLoS ONE 7:e29811).
# formal species, tested strain, provisional designation, call, evidence.
# "" for a strain the cited sentence does not name. Formal names follow Felix,
# Braendle & Cutter 2014 (PLoS ONE 9:e94723).
# Provisional-to-formal name mapping, so a reader can check the correspondence
# rather than trust it. Names from Felix, Braendle & Cutter 2014 (PLoS ONE
# 9:e94723); the provisional numbering is Kiontke et al. 2011 (BMC Evol Biol
# 11:339). "type strain" is that paper's type/reference strain, which for these
# species is the isolate Nuez & Felix 2012 tested.
NAME_MAP = [
 ("C. sp. 6",  "portoensis",     "EG4788", "Felix et al. 2014"),
 ("C. sp. 7",  "afra",           "JU1199", "Felix et al. 2014"),
 ("C. sp. 9",  "nigoni",         "JU1325", "Felix et al. 2014"),
 ("C. sp. 10", "doughertyi",     "JU1133", "Felix et al. 2014"),
 ("C. sp. 11", "tropicalis",     "JU1373", "Felix et al. 2014"),
 ("C. sp. 12", "castelli",       "JU1426", "Felix et al. 2014"),
 ("C. sp. 13", "virilis",        "JU1528", "Felix et al. 2014"),
 ("C. sp. 14", "imperialis",     "EG5716", "Felix et al. 2014"),
 ("C. sp. 15", "kamaaina",       "QG122",  "Felix et al. 2014"),
 ("C. sp. 16", "wallacei",       "JU1873", "Felix et al. 2014"),
 ("C. sp. 17", "nouraguensis",   "JU1825", "Felix et al. 2014"),
 ("C. sp. 18", "macrosperma",    "JU1857", "Felix et al. 2014"),
 ("C. sp. 19", "yunquensis",     "EG6142", "Felix et al. 2014"),
 ("C. sp. 20", "guadeloupensis", "NIC113", "Felix et al. 2014"),
 ("C. sp. 23", "latens",         "VX88",   "Felix et al. 2014"),
]

# Reference topology: the Open Tree of Life induced subtree for these species,
# synthetic tree opentree16.1 (2025-12-20), taxonomy 3.7draft3, retrieved from
# api.opentreeoflife.org/v3/tree_of_life/induced_subtree. Topology only, no
# branch lengths, so the figure draws it as a cladogram. Used ONLY to order and
# group rows; nothing here re-estimates a phylogeny, and a SID-2 gene tree would
# be the wrong thing to draw at 42-48% identity.
SPECIES_TREE = (
    "((((((((((((Caenorhabditis_brenneri_ott90647,(Caenorhabditis_dougherty"
    "i_ott624496,(Caenorhabditis_wallacei_ott624497,Caenorhabditis_tropical"
    "is_ott5701053)mrcaott624497ott5701053)mrcaott624496ott624497)mrcaott90"
    "647ott624496,(((Caenorhabditis_briggsae_ott395053,Caenorhabditis_nigon"
    "i_ott5701061)mrcaott395053ott5701061,(Caenorhabditis_sinica_ott571283)"
    "mrcaott571283ott7073354)mrcaott395053ott571283,(Caenorhabditis_remanei"
    "_ott396902,Caenorhabditis_latens_ott5490377)mrcaott396902ott5490377)mr"
    "caott395053ott396902)mrcaott90647ott395053,(Caenorhabditis_elegans_ott"
    "395048,Caenorhabditis_inopinata_ott7073336)mrcaott395048ott7073336)mrc"
    "aott90647ott395048,Caenorhabditis_kamaaina_ott624499)mrcaott90647ott62"
    "4499,(((Caenorhabditis_afra_ott102813)mrcaott102813ott7073353,((((Caen"
    "orhabditis_nouraguensis_ott454025)mrcaott454025ott7073333,Caenorhabdit"
    "is_yunquensis_ott624492)mrcaott454025ott624492,(Caenorhabditis_panamen"
    "sis_ott7073338)mrcaott7073338ott7073357)mrcaott454025ott7073338,Caenor"
    "habditis_macrosperma_ott624491)mrcaott454025ott624491)mrcaott102813ott"
    "454025,(Caenorhabditis_japonica_ott215930,Caenorhabditis_imperialis_ot"
    "t624498)mrcaott215930ott624498)mrcaott102813ott215930)mrcaott90647ott1"
    "02813)mrcaott90647ott7073355,Caenorhabditis_astrocarya_ott7073332)mrca"
    "ott90647ott7073332,((Caenorhabditis_angaria_ott94476,Caenorhabditis_ca"
    "stelli_ott5701050)mrcaott94476ott5701050,(Caenorhabditis_dolens_ott707"
    "3335)mrcaott7073335ott7073340)mrcaott94476ott7073335)mrcaott90647ott94"
    "476,Caenorhabditis_portoensis_ott102814,(Caenorhabditis_drosophilae_ot"
    "t362927,Caenorhabditis_virilis_ott624494)mrcaott215931ott362927)mrcaot"
    "t90647ott102814,Caenorhabditis_plicata_ott215924)mrcaott90647ott215924"
    ")mrcaott90647ott7073339,Caenorhabditis_guadeloupensis_ott102817,Caenor"
    "habditis_monodelphis_ott624495)Caenorhabditis_ott395055;"
)

RNAI = [
 ("elegans",     "N2",     "",          "sensitive",        "prose"),
 ("portoensis",  "EG4788", "C. sp. 6",  "sensitive",        "prose"),
 ("afra",        "JU1199", "C. sp. 7",  "sensitive",        "prose"),
 ("virilis",     "JU1528", "C. sp. 13", "sensitive",        "prose"),
 ("imperialis",  "EG5716", "C. sp. 14", "sensitive",        "prose"),
 ("kamaaina",    "QG122",  "C. sp. 15", "sensitive",        "prose"),
 ("wallacei",    "JU1873", "C. sp. 16", "weakly sensitive", "prose"),
 ("briggsae",    "AF16",   "",          "insensitive",      "prose"),
 ("remanei",     "",       "",          "insensitive",      "prose"),
 ("brenneri",    "",       "",          "insensitive",      "prose"),
 ("drosophilae", "DF5077", "",          "insensitive",      "prose"),
 ("tropicalis",  "JU1373", "C. sp. 11", "insensitive",      "Table S1 and S2"),
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
    # Both proteome sets, on one scale. The CGP species sit at depth_rank 3
    # with the other non-Elegans, non-Japonica Caenorhabditis: the CGP release
    # carries no group assignment we could cite, and the claim this table
    # supports is about the genus boundary, not about structure inside it.
    with open(f"{OUT}/sid2_ortholog_search.tsv", "w") as fh:
        fh.write("source\tproteome\tspecies\tstrain\tgroup\tdepth_rank"
                 "\taccession\tpercent_identity\tquery_coverage\tevalue"
                 "\tis_ortholog\n")
        for upid, sp, grp, rank, acc, pid, cov, ev in SEARCH:
            fh.write(f"UniProt\t{upid}\t{sp}\t\t{grp}\t{rank}\t{acc}\t"
                     f"{pid}\t{cov}\t{ev:.3g}"
                     f"\t{'TRUE' if ev < E_ORTH else 'FALSE'}\n")
        for (sp, strain, acc, pid, cov, ev, w, blk) in CGP:
            code = CGP_CODE.get(sp) or acc.split(".")[0]
            fh.write(f"CGP\t{code}\t{cgp_species(sp)}\t{strain}\t"
                     f"Caenorhabditis (CGP v2)\t3\t{acc or ''}\t"
                     f"{'' if pid is None else pid}\t"
                     f"{'' if cov is None else cov}\t"
                     f"{'' if ev is None else f'{ev:.3g}'}\t"
                     f"{'TRUE' if ev is not None and ev < E_ORTH else 'FALSE'}\n")
    n_up = sum(1 for r in SEARCH if r[7] < E_ORTH)
    n_cgp = sum(1 for r in CGP if r[5] is not None and r[5] < E_ORTH)
    n_prot = len(SEARCH) + len(CGP)
    print(f"  search table: {n_prot} proteomes "
          f"({len(SEARCH)} UniProt + {len(CGP)} CGP, "
          f"{n_prot - len(CGP_DUP)} distinct species), "
          f"{n_up + n_cgp} with an ortholog at E < {E_ORTH:g} "
          f"({n_up} UniProt incl. the query, {n_cgp} CGP)")
    # the claim the figure makes: every proteome that clears the threshold is a
    # Caenorhabditis, and no proteome outside the genus does
    out_genus = [r for r in SEARCH if not r[1].startswith("Caenorhabditis")]
    assert all(r[7] >= E_ORTH for r in out_genus), \
        "a proteome outside Caenorhabditis now clears the ortholog threshold"
    assert len(out_genus) == 9, "the outgroup ladder changed size"

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

    # ---- the CGP window survey ------------------------------------------
    with open(f"{OUT}/sid2_ortholog_window_survey.tsv", "w") as fh:
        fh.write("source\tspecies\tstrain\taccession\tpercent_identity"
                 "\tquery_coverage\tevalue\twindow_88_104\taa94\taa95\taa96"
                 "\tblock_identity\tconfident\tsequon\n")
        def row(src, sp, strain, acc, pid, cov, ev, w, blk, conf):
            a94 = a95 = a96 = ""
            if w:
                a94, a95, a96 = w[WIN.index(94)], w[WIN.index(95)], w[WIN.index(96)]
            seq = bool(w) and conf and a94 == "N" and a96 in "ST" and a95 != "P"
            fh.write(f"{src}\t{sp}\t{strain}\t{acc}\t{pid}\t{cov}\t{ev}\t"
                     f"{w or ''}\t{a94}\t{a95}\t{a96}\t{blk}\t"
                     f"{'TRUE' if conf else 'FALSE'}\t"
                     f"{'TRUE' if seq else 'FALSE'}\n")
            return seq, (a94, a95, a96) if w else None
        ## the UniProt orthologs, window taken from the same Ce-anchored map.
        ## Their identity and coverage come from SEARCH so that every row in
        ## this table sits on one scale and the figure can order by identity.
        by_sp = {r[1]: r for r in SEARCH}
        for sp, acc, grp, ambig in ALIGN_SET:
            short = sp.replace("Caenorhabditis ", "")
            w = "".join(mapped[sp][i - 1] for i in WIN)
            conf = (not ambig) and grp == "Elegans group"
            r = by_sp.get(sp)
            row("UniProt", short, "", acc,
                "" if r is None else r[5], "" if r is None else r[6],
                "" if r is None else f"{r[7]:.3g}", w, "", conf)
        for (sp, strain, acc, pid, cov, ev, w, blk) in CGP:
            conf = blk is not None and blk >= CGP_FLOOR
            row("CGP", sp, strain, acc or "", "" if pid is None else pid,
                "" if cov is None else cov, "" if ev is None else f"{ev:.3g}",
                w, "" if blk is None else blk, conf)

    # the combined confident set: the six full-length Elegans-group orthologs
    # plus every CGP species clearing the block-identity floor
    cgp_conf = [r for r in CGP if r[7] is not None and r[7] >= CGP_FLOOR]
    def trip(w):
        return w[WIN.index(94)], w[WIN.index(95)], w[WIN.index(96)]
    combined = [(mapped[sp][93], mapped[sp][94], mapped[sp][95])
                for sp in CONS_GROUP] + [trip(r[6]) for r in cgp_conf]
    n96 = sum(1 for a, b, c in combined if c in "ST")
    nK = sum(1 for a, b, c in combined if c == "K")
    nseq = sum(1 for a, b, c in combined if a == "N" and c in "ST" and b != "P")
    print(f"  CGP survey: {len(CGP)} species, {len(cgp_conf)} clear the "
          f"{CGP_FLOOR:.0f}% block-identity floor")
    print(f"  combined confident set ({len(combined)} species): residue 96 is "
          f"Ser/Thr in {n96}, Lys in {nK}; sequon intact in {nseq}")
    assert nK == 0, "a confidently aligned Caenorhabditis now carries Lys96"
    assert n96 >= len(combined) - 2, "residue 96 is no longer Ser/Thr in nearly all"
    assert nseq >= 11, "the sequon count has dropped"

    # ---- the environmental-RNAi layer -----------------------------------
    ## the sequenced isolate per species, so the mismatch can be shown
    ## the isolate each proteome was sequenced from. The UniProt records carry
    ## no strain field, so these were resolved through the assembly accession
    ## in each proteome record and the NCBI dataset report for it.
    seqd = {sp.replace("Caenorhabditis ", ""): "" for sp, *_ in ALIGN_SET}
    seqd.update({"elegans": "N2", "briggsae": "AF16", "remanei": "PB4641",
                 "brenneri": "PB2801", "nigoni": "JU1422",
                 "tropicalis": "JU1373"})
    for (sp, strain, acc, pid, cov, ev, w, blk) in CGP:
        seqd.setdefault(sp, strain)
    with open(f"{OUT}/sid2_env_rnai_sensitivity.tsv", "w") as fh:
        fh.write("species\tprovisional_name\ttested_strain\tsequenced_strain"
                 "\tsame_isolate\tresponse\tevidence\tsource\n")
        for sp, strain, prov, call, ev in RNAI:
            sq = seqd.get(sp, "")
            same = "TRUE" if (sq and strain and sq == strain) else \
                   ("FALSE" if (sq and strain) else "")
            fh.write(f"{sp}\t{prov}\t{strain}\t{sq}\t{same}\t{call}\t{ev}\t"
                     f"Nuez & Felix 2012 PLoS ONE 7:e29811\n")
    print(f"  environmental-RNAi layer: {len(RNAI)} species with a published "
          f"call")

    with open(f"{OUT}/sid2_species_name_map.tsv", "w") as fh:
        fh.write("provisional_name\tformal_name\ttype_strain\tsource\n")
        for prov, formal, strain, src in NAME_MAP:
            fh.write(f"{prov}\t{formal}\t{strain}\t{src}\n")
    print(f"  name map: {len(NAME_MAP)} provisional designations")

    with open(f"{OUT}/sid2_species_tree.nwk", "w") as fh:
        fh.write(SPECIES_TREE + "\n")
    n_tip = SPECIES_TREE.count("Caenorhabditis_")
    print(f"  reference topology: {n_tip} tips (Open Tree opentree16.1)")

    ## the overlay: species that are BOTH confidently aligned and called
    called = {sp: call for sp, _st, _p, call, _e in RNAI}
    conf_sp = [sp for sp in CONS_GROUP] + \
              [r[0] for r in CGP if r[7] is not None and r[7] >= CGP_FLOOR]
    conf_sp = [s2.replace("Caenorhabditis ", "") for s2 in conf_sp] + ["elegans"]
    overlay = [(sp, called[sp]) for sp in conf_sp if sp in called]
    ## "insensitive".endswith("sensitive") is True -- match the calls exactly
    sens = sum(1 for _s, c in overlay if c in ("sensitive", "weakly sensitive"))
    insens = sum(1 for _s, c in overlay if c == "insensitive")
    print(f"  overlay: {len(overlay)} confidently aligned species have a call; "
          f"{sens} sensitive or weakly so, {insens} insensitive")
    assert sens >= 3 and insens >= 3, \
        "the overlay no longer shows sensitive and insensitive species sharing " \
        "the same residue-96 state, which is the point it makes"
    assert called.get("afra") == "sensitive", \
        "C. afra is the natural AxT that is sensitive; the figure says so"

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
    ## C. angaria is NOT the claimed breakdown point any more: its window sits
    ## at 26% block identity, below the floor, so it is counted neither way.
    ## What must hold is that nothing outside the genus is alignable.
    assert all(r[7] >= E_ORTH for r in SEARCH if r[3] >= 4), \
        "a proteome outside Caenorhabditis now yields a significant hit"
    seq_intact = sum(1 for sp in CONS_GROUP
                     if mapped[sp][93] == "N" and mapped[sp][95] in "ST"
                     and mapped[sp][94] != "P")
    assert seq_intact == n, "the N94 sequon is no longer intact in every ortholog"
    print(f"  N94-x-[ST] intact in {seq_intact}/{n} full-length "
          f"Elegans-group orthologs")
    print("  all assertions passed")


if __name__ == "__main__":
    main()
