#!/usr/bin/env python3
"""Assemble manuscript_supplement/ -- the supplemental tables the manuscript
needs, numbered in order of first mention in the text.

    python3 scripts/make_manuscript_supplement.py

    writes  manuscript_supplement/Table_S*.{tsv,csv,...}
            manuscript_supplement/TABLE_INDEX.csv
            manuscript_supplement/README.md

WHY THE ORDER IS WHAT IT IS. supplemental_data/ is organised by analysis stage
(deconvolution, mapping, phenotypes, structure), which is how the pipeline
thinks and not how a reader moves through the paper. The numbering here is the
order of FIRST MENTION in the Results, so Table S1 is the first thing a reader
needs and Table Sn the last. Files that no manuscript figure or claim depends on
are deliberately left in supplemental_data/ and are listed at the end of
README.md as the deposit-only remainder.

EVERY ENTRY CARRIES AN ANCHOR: a short verbatim phrase from the manuscript,
checked against the extracted text at build time, marking where the citation
belongs. The point is that the anchor can be pasted into a find box.

THE ONE TABLE THE MANUSCRIPT ALREADY ASKS FOR is the cross QTL summary. The
Results say "Figure S7B-C, Supplemental table X" -- an unresolved placeholder --
and that is Table S16 here.
"""

import csv
import os
import re
import shutil
import sys

OUT = "manuscript_supplement"
TEXT = "/tmp/bulk_flat.txt"          # flattened Bulk Paper.pdf, see README
SD = "supplemental_data"

# (number, published stem, source path(s), description, manuscript anchor)
# Sources are relative to the repository root. Several entries carry more than
# one file because they are one table with more than one part.
TABLES = [
 (1, "pos1_pilot_strain_phenotypes", ["supplemental_data/phenotypes/pos1_2023_association_traits.csv.gz"],
  "Per-strain pos-1 RNAi response for the 231-isotype pilot pool. Columns: "
  "delta_ctrl_pos-1_T2 (RNAi minus control change in pool frequency), "
  "vst_ctrl_pos-1_T2 (variance-stabilised, the trait the association scan was "
  "run on), log2fc_ctrl_pos-1_T2 (the log2 ratio), and two negative-control "
  "columns, growth on HT115 against t0 and log10 control abundance.",
  "231 pooled wild isolates"),

 (2, "simulation_recovery_by_depth", ["supplemental_data/deconvolution/simulation_seeded_r2.tsv",
   "supplemental_data/deconvolution/simulation_fitness_traits.tsv"],
  "Deconvolution recovery in simulation. r-squared of estimated against known "
  "input frequency for each of seven fitness traits at eight sequencing depths, "
  "over ten seeded replicates, plus the seven trait vectors used as fitness and "
  "the strains carrying each.",
  "NNLS recovered strain frequencies"),

 (3, "known_mixture_design_and_recovery", ["supplemental_data/deconvolution/dilution_design.tsv",
   "supplemental_data/deconvolution/dilution_strain_sets.tsv",
   "supplemental_data/deconvolution/dilution_predictions_poolref.tsv.gz"],
  "The designed DNA mixture, in three parts: (a) the titration design, giving "
  "each sample's set B and set C volumes, water, total DNA mass and nominal set "
  "fractions; (b) which of the four sets each strain belongs to, with its "
  "isotype; (c) the recovered per-strain frequency in every library.",
  "root mean squared error"),

 (4, "nnls_vs_mipseq_frequencies", ["supplemental_data/deconvolution/baugh_nnls_dep103_with_mipseq.tsv.gz",
   "supplemental_data/deconvolution/mipseq_frequencies.txt.gz"],
  "Per-strain, per-sample frequencies from this platform (NNLS on whole-genome "
  "sequence) beside the published MIP-seq frequencies for the same 23 samples of "
  "the L1 starvation time course. The basis for the platform comparison.",
  "highly correlated with those derived from MIP-seq"),

 (5, "downsampling_recovery", ["supplemental_data/deconvolution/baugh_downsample_trait_recovery.tsv"],
  "Recovery of the two published traits at each subsampled depth (0.25, 0.5, 1, "
  "3, 5 and 10x): Spearman correlation against the published PC1 and Slope, and "
  "against Slope computed on the difference-based definition used here.",
  "as we downsampled reads"),

 (6, "pos1_pilot_replicate_frequencies", ["supplemental_data/phenotypes/pos1_2023_sample_frequencies.csv.gz"],
  "Sample-level strain frequencies for the pilot pool: four pos-1 replicate pools "
  "and two control pools at three read-depth cutoffs. Supports the "
  "replicate-agreement comparison and the per-isolate responsiveness counts.",
  "across all six pairwise replicate comparisons"),

 (7, "pos1_pilot_gwas", ["supplemental_data/mapping/pos1_2023_gemma_loco.csv.gz"],
  "Genome-wide association results for the pilot pos-1 response, 464,045 markers "
  "over 231 strains, GEMMA with leave-one-chromosome-out kinship.",
  "genome-wide association scan on the pos-1 responses"),

 (8, "gwas_significance_thresholds", ["supplemental_data/mapping/eigen_independent_tests.tsv"],
  "Effective number of independent tests by the Li and Ji (2005) eigenvalue "
  "method, per panel and per chromosome, with the marker counts, the imputed "
  "fraction, the trace check, and the resulting Bonferroni and eigenvalue "
  "thresholds drawn on every association panel.",
  "passed the Bonferroni significance threshold"),

 (9, "manual_plate_scores", ["supplemental_data/phenotypes/plate_scores_pos1.tsv"],
  "Manual pos-1 RNAi score per wild isolate, on a six-level ordinal scale where "
  "0 is a complete RNAi response and 5 is none. Two columns, strain and score; "
  "which strains were carried into the 93-strain pool is not recorded here.",
  "191 wild strains"),

 (10, "paaby2015_comparison", ["supplemental_data/phenotypes/paaby2015_embryonic_lethality.txt.gz"],
  "Well-level embryo and larva counts from Paaby et al. 2015, restricted to the "
  "pos-1 clone, as used for the independent comparison against the manual plate "
  "scores.",
  "previously published evaluation of wild isolate RNAi responses"),

 (11, "pooled_ten_target_phenotypes", ["supplemental_data/phenotypes/pooled_vst_traits.csv.gz"],
  "Per-strain responses for the 93-strain pool across all ten RNAi targets "
  "against the HT115 empty-vector control, one column per target in each of the "
  "three parameterisations (delta_ctrl, log2fc, vst). The vst columns are what "
  "the association scans were run on.",
  "each of ten target genes"),

 (12, "pooled_gwas_and_cross_scans", ["supplemental_data/mapping/pooled_cross_bundle_thinned.rds"],
  "The thinned mapping bundle behind Figure 2: pooled association scans for all "
  "ten targets and the cross contrast scans, one marker per bin retained on "
  "maximum LOD so every peak survives. Peak positions and heights are identical "
  "to the full data; interval widths may be 0-9 kb narrower.",
  "exceeded the Bonferroni threshold for pos-1 and mig-6"),

 (13, "cross_allele_frequencies", ["supplemental_data/mapping/cross_af_JU1793xJU2466.tsv.gz",
   "supplemental_data/mapping/cross_af_N2xXZ1516.tsv.gz",
   "supplemental_data/mapping/cross_af_samples.tsv"],
  "Genome-wide parental allele frequencies for both advanced-intercross crosses, "
  "per sample, with the sample sheet naming the condition and timepoint of each.",
  "allele frequency deviations between conditions"),

 (14, "cross_qtl_summary", ["plots/TABLE_cross_qtl_full.tsv"],
  "THE TABLE THE TEXT ALREADY CITES AS \"Supplemental table X\". Every cross QTL "
  "peak above the genome-wide threshold in every contrast (495 rows): position, "
  "peak LOD, support interval, parental frequencies either side, rank on its "
  "chromosome, whether the trough separates it from a taller peak, and whether "
  "the other cross calls it by interval overlap or by position.",
  "Supplemental table X"),

 (15, "cross_locus_classification", ["plots/TABLE_cross_qtl_locus_classification.tsv",
   "plots/TABLE_cross_qtl_condition_dfreq.tsv"],
  "General against target-specific classification of each cross locus: how many "
  "RNAi targets respond at it and in which direction, with the per-condition "
  "frequency differences the call is made on, swept over four response cutoffs.",
  "QTL that are shared between conditions"),

 (16, "nil_introgression_ranges", ["supplemental_data/hatching_assays/nil_introgression_ranges.bed"],
  "Introgression boundaries for the NIL series as BED: chromosome, start, end, "
  "strain and donor parent, with the two parents carried as whole-chromosome "
  "rows. These are the sequence-confirmed lines, and the difference between the "
  "wSZ191 and wSZ196 rows is the 37 kb interval. Which boundaries are "
  "sequence-confirmed as against marker-bounded is in Table S25, not here.",
  "37 kb interval"),

 (17, "nil_hatching_counts", ["supplemental_data/hatching_assays/nil_series_hatching.tsv"],
  "Embryo hatching for every NIL and both parents on pos-1 RNAi and on control "
  "food. Columns: strain, condition, embryos plated, number unhatched, and "
  "hatched fraction. One plate per strain per condition, so the binomial "
  "intervals quoted in the text are computed from these counts rather than "
  "stored here.",
  "distinct pos-1 RNAi phenotypes"),

 (18, "nil_interval_content", ["supplemental_data/mapping/nil_interval_genes.tsv",
   "supplemental_data/mapping/nil_interval_parent_variants.tsv",
   "supplemental_data/mapping/nil_interval_exons.tsv"],
  "What the 37 kb interval contains: every gene, every difference between the "
  "cross parents with its annotated consequence and impact class, and the exon "
  "models. The basis for the two protein-altering differences in sid-2.",
  "27 parental differences across twelve genes"),

 (19, "sid2_allele_swap_hatching_ju", ["supplemental_data/hatching_assays/ju_allele_swaps_hatching.csv"],
  "Embryo hatching for the JU1793 and JU2466 sid-2 allele swaps and their "
  "parents on pos-1 and control food. Columns: experiment, strain, genotype, "
  "glycosylation motif, condition, embryos plated, number unhatched, hatched "
  "fraction. The motif column is what distinguishes the residue-94 and "
  "residue-96 states (NxT, NxK, AxT, AxK) and records that JU2466 appears as two "
  "isolates, A and B, with the 96T edit made in A. Note that this file calls the "
  "JU1793 N94A strain wSZ207 where the stock sheet calls it wSZ209.",
  "reciprocal allele-swap strains"),

 (20, "sid2_allele_swap_hatching_n2", ["supplemental_data/hatching_assays/n2_allele_swaps_hatching.tsv"],
  "Embryo hatching for N2 and the two N2 sid-2 96K lines (wSZ203, wSZ204) across "
  "the 0, 25, 50, 75 and 100% pos-1 food dose series. One row per strain per "
  "dose; the pooled 96K figure quoted in the text is computed from the two lines' "
  "counts rather than stored as a row.",
  "lowered the N2 hatching rate"),

 (21, "sid2_ortholog_conservation", ["supplemental_data/structure/sid2_ortholog_conservation.tsv",
   "supplemental_data/structure/sid2_ortholog_search.tsv",
   "supplemental_data/structure/sid2_ortholog_window_survey.tsv",
   "supplemental_data/structure/sid2_species_name_map.tsv"],
  "Conservation of the N94-C95-T96 sequon across the Caenorhabditis species "
  "surveyed, with the ortholog search results, the alignable-window survey "
  "bounding how far SID-2 can be compared, and the species name map.",
  "Elegans supergroup"),

 (22, "sid2_environmental_competence", ["supplemental_data/structure/sid2_env_rnai_sensitivity.tsv"],
  "Published environmental-RNAi competence per Caenorhabditis species, curated "
  "from the literature: the species, the strain tested for RNAi response, the "
  "strain whose genome was sequenced, whether those are the same isolate, the "
  "reported response, the kind of evidence, and the source. Layered onto the "
  "conservation figure.",
  "sensitive to ingested double stranded RNA"),

 (23, "sid2_local_charge", ["supplemental_data/structure/sid2_local_charge.tsv",
   "supplemental_data/structure/sid2_per_residue.tsv"],
  "Local net charge in 12 A windows across the SID-2 ectodomain, per residue, at "
  "two pH values (4.4 and 7.4), with the number of residues inside each window, "
  "the model pLDDT and the residue coordinates. The second file carries the "
  "per-residue model quantities behind it: pLDDT, assigned secondary structure, "
  "topology region and coordinates.",
  "local net charge in 12"),

 (24, "sid2_population_variants", ["supplemental_data/structure/sid2_variants_cendr.tsv",
   "supplemental_data/structure/sid2_population_missense.tsv",
   "supplemental_data/structure/sid2_parental_variants.tsv"],
  "sid-2 variation across the wild population from CaeNDR, the missense variants "
  "among it, and the differences between the parents of each cross, including the "
  "T96K site and its allele assignment.",
  "high-frequency variant in the dsRNA transporter SID-2"),

 (25, "strains", ["supplemental_data/Table_S_strains.csv"],
  "Every strain constructed for this work: QX designation, lab identifier, "
  "introgression or allele designation, genotype, background, construction route, "
  "guide and repair template used, and verification status.",
  "Methods: NIL Construction"),

 (26, "oligonucleotides", ["supplemental_data/Table_S_oligos.csv"],
  "Every oligonucleotide: genotyping primers with their pairings and expected "
  "products, restriction assays, the four sid-2 repair templates and the four "
  "guide RNAs.",
  "CRISPR Design to edit sid-2"),
]


def main():
    if not os.path.exists(TEXT):
        sys.exit(f"missing {TEXT}: re-extract Bulk Paper.pdf first (see README.md)")
    text = open(TEXT).read()

    # anchors must actually occur in the manuscript, or they are useless
    # An anchor has to be present AND unique, or it cannot be used to find the
    # sentence. Both are build-time errors rather than something to eyeball.
    missing = [(n, a) for n, _, _, _, a in TABLES
               if not a.startswith("Methods:") and a not in text]
    if missing:
        sys.exit("anchors not found in the manuscript text: " + repr(missing))
    ambiguous = [(n, a, text.count(a)) for n, _, _, _, a in TABLES
                 if not a.startswith("Methods:") and text.count(a) > 1]
    if ambiguous:
        sys.exit("anchors that match more than once: " + repr(ambiguous))

    nums = [t[0] for t in TABLES]
    assert nums == list(range(1, len(TABLES) + 1)), "table numbers must be 1..n"

    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)

    rows = []
    for num, stem, srcs, desc, anchor in TABLES:
        published = []
        for src in srcs:
            if not os.path.exists(src):
                sys.exit(f"Table S{num}: missing source {src}")
            ext = "".join(
                s for s in [os.path.splitext(src)[1]] if s)
            if src.endswith(".tsv.gz"): ext = ".tsv.gz"
            elif src.endswith(".csv.gz"): ext = ".csv.gz"
            elif src.endswith(".txt.gz"): ext = ".txt.gz"
            part = "" if len(srcs) == 1 else f"_{chr(ord('a') + srcs.index(src))}"
            name = f"Table_S{num:02d}_{stem}{part}{ext}"
            shutil.copy2(src, os.path.join(OUT, name))
            published.append(name)
        nbytes = sum(os.path.getsize(os.path.join(OUT, f)) for f in published)
        # Journals cap supplemental files and expect something a reader can
        # open. Anything large, or in an R binary format, is a deposit item and
        # belongs in the data repository rather than the article supplement.
        binary = any(f.endswith((".rds", ".rda", ".RData")) for f in published)
        kind = "dataset" if (nbytes > 5e6 or binary) else "table"
        rows.append({
            "table": f"Table S{num}",
            "kind": kind,
            "file": "; ".join(published),
            "size_MB": f"{nbytes / 1e6:.2f}",
            "description": desc,
            "manuscript_anchor": anchor,
            "source_in_repo": "; ".join(srcs),
        })

    with open(os.path.join(OUT, "TABLE_INDEX.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["table", "kind", "file", "size_MB",
                                           "description", "manuscript_anchor",
                                           "source_in_repo"])
        w.writeheader()
        w.writerows(rows)

    # the remainder: present in supplemental_data but not promoted
    promoted = {s for _, _, srcs, _, _ in TABLES for s in srcs}
    remainder = []
    for root, _, fs in os.walk(SD):
        for f in fs:
            if f.startswith("."): continue
            p = os.path.join(root, f)
            if p not in promoted:
                remainder.append(p)

    with open(os.path.join(OUT, "README.md"), "w") as fh:
        fh.write(README_HEAD)
        for r in rows:
            tag = "" if r["kind"] == "table" else "  *(dataset -- deposit, not an article supplement)*"
            fh.write(f"### {r['table']} — `{r['file']}`{tag}\n\n")
            fh.write(f"{r['description']}\n\n")
            fh.write(f"**Cite near:** “{r['manuscript_anchor']}”\n\n")
            fh.write(f"*Source in repository:* `{r['source_in_repo']}`  "
                     f"({r['size_MB']} MB)\n\n")
        fh.write("\n## Not promoted\n\n")
        fh.write(f"{len(remainder)} files in `supplemental_data/` support no "
                 "manuscript figure or claim and stay where they are -- "
                 "alternative deconvolution references, the mig-6 locus census "
                 "that no cited figure draws, intermediate caches, and the "
                 "structure files behind figures the manuscript does not cite.\n\n")
        for p in sorted(remainder):
            fh.write(f"- `{p}`\n")

    n_ds = sum(1 for r in rows if r["kind"] == "dataset")
    print(f"{OUT}/: {len(TABLES)} entries over "
          f"{sum(len(t[2]) for t in TABLES)} files -- "
          f"{len(rows) - n_ds} reader-facing tables, {n_ds} datasets")
    print(f"  index: {OUT}/TABLE_INDEX.csv")
    print(f"  not promoted: {len(remainder)} files left in {SD}/")


README_HEAD = """# Manuscript supplement

Supplemental tables for *Flexible pooled phenotyping enables population-scale
mapping of RNAi-sensitivity modifiers*, numbered in order of first mention in
the Results.

Built by `scripts/make_manuscript_supplement.py` from `supplemental_data/` and
`plots/`. Rerunning it rebuilds this directory from scratch, so edit the script
rather than the files here. The script reads a flattened text extraction of
`Bulk Paper.pdf` at `/tmp/bulk_flat.txt` in order to check every anchor against
the manuscript; regenerate it with:

```sh
python3 -c "
import re, pypdf
r = pypdf.PdfReader('Bulk Paper.pdf')
t = '\\n'.join((p.extract_text() or '') for p in r.pages)
open('/tmp/bulk_flat.txt','w').write(re.sub(r'\\s+',' ',t))"
```

**Cite near** gives a short verbatim phrase from the manuscript, checked against
the text at build time. Paste it into a find box to land on the sentence where
the table reference belongs.

`TABLE_INDEX.csv` is the same list as a spreadsheet.

`kind` separates the two things in here. A **table** is small and opens in a
spreadsheet, and is what a journal expects as a supplemental table. A
**dataset** is either large or in an R binary format; those belong in the data
repository the paper deposits to, and should be cited under data availability
rather than attached to the article. The datasets are most of the bulk here.

"""


if __name__ == "__main__":
    main()
