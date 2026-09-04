# Bulk GWA manuscript — RNAi sensitivity in wild *Caenorhabditis elegans*

Analysis code and figures for a manuscript on natural variation in RNAi
sensitivity: pooled competition assays across wild isolates, genome-wide
association mapping, F2 bulk-segregant crosses, near-isogenic line
fine-mapping, and characterisation of a *sid-2* coding variant.

Everything here runs from the repository root with no arguments.

## What is in this repository

| path | contents |
|---|---|
| `scripts/` | 16 figure scripts, one per figure and named for it, plus 12 shared, data-prep and rendering scripts |
| `scripts/legacy/` | 30 superseded and orphaned scripts, kept so earlier figures can be reproduced; they write to `plots/legacy/` |
| `plots/` | every figure as PDF and PNG, plus the structure renders in `plots/assets/` |
| `data/` | the **small** input tables the figure scripts read — trait tables, GWAS result exports, plate phenotyping, NIL breakpoints, structure models, *sid-2* variant tables |
| `METHODS.txt` | draft methods text, with every number traceable to the script that produces it |
| `FIGURE_CAPTIONS.txt` | a caption for every figure and supplement, with the numbers and the caveats |
| `DATA_AVAILABILITY.md` | what is archived externally, and which figures need it |
| `main_displays.pdf` | the assembled main display |

Tracked size is about 72 MB.

## What is excluded

The working tree is roughly **12.8 GB**; about 12.7 GB of it is excluded by
`.gitignore` and archived separately. Excluded, in descending size:

- `data/genotypes/` — CeNDR PLINK panel and processed genotype matrix (7.1 GB)
- `data/structure_modeling/` — exploratory AlphaFold, HADDOCK and
  electrostatics work (2.0 GB); one PDB and the text reports are kept
- `data/cross_experiments/*/` — xQTL exports: allele counts, AFD tables and all
  contrast tables (1.4 GB); one contrast table is kept
- `data/pos1_original/gemma_output/`, `data/pos1_rnai_sensitive/gemma/` —
  per-chromosome GEMMA output from earlier runs (0.7 GB)
- `data/pooled_RNAi_expt/` — pooled competition raw data (0.7 GB)
- `data/pooled_cross_intersection/` — regenerable cache (281 MB)
- `data/cross_experiments/Nov2024_JU_cross_pos_mig/` — superseded analysis (245 MB)
- `data/LD/` — LD matrices (56 MB)
- `plots/pooled_cross_intersection/`, `plots/legacy/` — exploratory and
  superseded figure sets (96 MB); the variant tables in the first are kept
- large R objects (`*.Rdata`, `*.rda`, `*.rds`), alignment and variant files
  (`*.bam`, `*.vcf.gz`), PLINK binaries, GATK count tables and MSAs, wherever
  they appear under `data/` — with named exceptions for the few small ones the
  figures read

`.gitignore` excludes by rule, not by listing files, so unpacking the archive
over the repository root restores the tree and git keeps ignoring it.

`*.key` is excluded: a Keynote bundle is rewritten in full on every save and
would grow history quickly. `main_displays.pdf` is tracked instead.

## The curated figure set

`plots/` holds exactly the figures intended for the manuscript. Superseded
variants live in `plots/legacy/`, which is not tracked; the scripts that build
them now write there too, so re-running anything will not put them back in
`plots/`.

| figure | script |
|---|---|
| `Figure1_pos1` | `Figure1_pos1.R` |
| `Figure2` | `Figure2.R` |
| `Figure3_quad` | `Figure3_quad.R` |
| `Figure4_sid2` | `Figure4_sid2.R` |
| `SUPP_FIG_plate_vs_paaby_vs_pos1original` | `SUPP_FIG_plate_vs_paaby_vs_pos1original.R` |
| `SUPP_FIG_XX_cross_contrast_panels` | `SUPP_FIG_XX_cross_contrast_panels.R` |
| `SUPP_FIG_XX_pooled_phenotype_ranks` | `SUPP_FIG_XX_pooled_phenotype_ranks.R` |
| `SUPP_FIG_XX_original_pos1_dfreq_rep_correlation` | `SUPP_FIG_XX_original_pos1_dfreq_rep_correlation.R` |
| `SUPP_FIG_XX_baugh_per_sample_frequencies` | `SUPP_FIG_XX_baugh_per_sample_frequencies.R` |
| `SUPP_FIG_XX_bootstrap_propagation_checks` | `SUPP_FIG_XX_bootstrap_propagation_checks.R` |
| `SUPP_FIG_XX_downsample_per_sample` | `SUPP_FIG_XX_downsample_per_sample.R` |
| `SUPP_FIG_XX_sid2_electrostatics` | `SUPP_FIG_XX_sid2_electrostatics.R` |
| `SUPP_FIG_XX_sid2_allele_in_panel` | `SUPP_FIG_XX_sid2_allele_in_panel.R` |
| `SUPP_FIG_XX_n2_swap_dose` | `SUPP_FIG_XX_n2_swap_dose.R` |
| `SUPP_FIG_XX_nil_hatching_full` | `SUPP_FIG_XX_nil_hatching_full.R` |
| `SUPP_FIG_XX_sid2_allele_swaps_full` | `SUPP_FIG_XX_sid2_allele_swaps_full.R` |

Each is present as both `.pdf` and `.png`, and each script is named for the
figure it produces. Verified by deleting every file in `plots/` and rebuilding:
all sixteen regenerate, with nothing extra.

`FIGURE_CAPTIONS.txt` still carries captions for the superseded variants,
marked where they are; they describe the same measurements in a different
arrangement and are useful when choosing between layouts.

## Regenerating the figures

Every figure is `Rscript scripts/<figure name>.R`. Two figures need their
assets built first:

```sh
python3 scripts/sid2_ribbon_render.py    # once, Figure4_sid2 structure panel
python3 scripts/sid2_zoom_render.py      # once, Figure4_sid2 structure panel
Rscript scripts/Figure4_sid2.R
```

The other twelve scripts in `scripts/` are shared panel builders
(`Figure1_common.R`, `Figure3_common.R`, `n2_swap_panels.R`), the
multiple-testing machinery (`gwas_thresholds.R`, `eigen_independent_tests.R`),
data preparation (`pooled_cross_intersection_prep.R`,
`pooled_cross_candidate_variation.R`, `baugh_L1_DownSample_Counts.R`,
`sid2_variant_table.R`) and an orientation-picking tool
(`sid2_orientation_sheet.py`).

**Three things need the external archive**, because they read the genotype
panel or the intersection cache:

- Figure 2 and the cross-QTL supplements
- `scripts/eigen_independent_tests.R`
- `scripts/SUPP_FIG_XX_sid2_allele_in_panel.R`

Figure 1 rebuilds from a committed NNLS result; only `FIG1_REFRESH=1`, which
redoes the deconvolution from the 31 MB genotype matrix, needs the archive.
Everything else builds from a clone. See `DATA_AVAILABILITY.md` for the table.

## Software

R 4.5.2 with tidyverse, data.table, patchwork, ggtext, RcppML, png.
Python 3.13 with NumPy, SciPy, Biopython, Matplotlib, Pillow.
PLINK 2.00a3, bcftools 1.11. GEMMA was run outside this repository; only its
exported association tables are here.

Two scripts, `scripts/plate_pheno_comparisons.R` and `scripts/pos1_script1.R`,
still resolve their working directory through the RStudio API and run only
inside RStudio. Every other script runs under `Rscript`.
