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
| `supplemental_data/` | **every input the figure scripts read**, 37.0 MB. Self-contained: all eighteen figures rebuild from this directory alone. Documented file by file, column by column, in `supplemental_data/SUPPLEMENTAL_DATA_OVERVIEW.md` |
| `METHODS.txt` | draft methods text, with every number traceable to the script that produces it |
| `FIGURE_CAPTIONS.txt` | a caption for every figure and supplement, with the numbers and the caveats |
| `DATA_AVAILABILITY.md` | what is archived externally |
| `SUPPLEMENTAL_DATA_SURVEY.md` | how the deposit set was chosen, and which files have no recipe here |
| `main_displays.pdf` | the assembled main display |

Tracked size is about 58 MB.

## What is excluded

`data/` is **not** in git. It is the provenance tree — 12.8 GB of raw sequence,
genotype panels, alignment intermediates, exports and caches — and it is the
*input* to the deposit rather than the deposit itself. Everything the figures
actually read has been staged into `supplemental_data/`, which is tracked.
`DATA_AVAILABILITY.md` lists the archive contents and
`SUPPLEMENTAL_DATA_SURVEY.md` explains how the 12.8 GB was reduced to 37.0 MB.

Two reductions are worth knowing about when using the deposit:

- **The pooled/cross bundle is thinned** to one marker per 5 kb (cross scans)
  and 10 kb (pooled GWAS), keeping the maximum-LOD marker so every peak
  survives exactly: 9.6 MB instead of 258.7 MB. Peak positions and peak LODs
  are identical to the full data; **QTL support interval widths are 0–9 kb
  narrower**. `scripts/compare_full_vs_thinned_bundle.R` rebuilds Figure 2 from
  both bundles and reports the difference.
- **The genotypes are a region subset**, chromosome III 13,679,000–13,682,000,
  because the two scripts that read genotypes need the *sid-2* region.

Also excluded: `plots/legacy/` and `plots/pooled_cross_intersection/`
(superseded and exploratory figures, regenerable from `scripts/legacy/`), and
`*.key` — a Keynote bundle is rewritten whole on every save and would grow
history quickly, so `main_displays.pdf` is tracked instead.

## The figure report

[`FIGURE_REPORT.md`](FIGURE_REPORT.md) presents all eighteen figures with their
captions, ordered by the argument the manuscript makes rather than by build
order: the assay, the map, the interval, the residue. It renders inline on
GitHub — click the link above.

`FIGURE_REPORT.html` is the same report as a single self-contained file, with
the figures embedded and click-to-zoom, for reading offline or on another
machine. Download it and open it in any browser; there is nothing to install.

Both come from `FIGURE_REPORT.Rmd`, which reads the figures from `plots/` and
recomputes its data tables from `supplemental_data/` on every knit, so the
report cannot drift from the deposit:

```sh
Rscript -e 'rmarkdown::render("FIGURE_REPORT.Rmd", output_format = "all")'
```

## The curated figure set

`plots/` holds exactly the figures intended for the manuscript. Superseded
variants live in `plots/legacy/`, which is not tracked; the scripts that build
them now write there too, so re-running anything will not put them back in
`plots/`.

| figure | script |
|---|---|
| `SUPP_FIG_XX_simulation_depth` | `SUPP_FIG_XX_simulation_depth.R` |
| `SUPP_FIG_XX_dilution_validation` | `SUPP_FIG_XX_dilution_validation.R` |
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
all eighteen regenerate, with nothing extra.

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
