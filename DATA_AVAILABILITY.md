# Data availability

This repository holds the analysis code, the figures, and `supplemental_data/`
— every input the figure scripts read, staged and compressed, 38.2 MB. All
eighteen figures rebuild from it alone, verified by deleting `data/` and running
them. See `supplemental_data/SUPPLEMENTAL_DATA_OVERVIEW.md` for a file-by-file,
column-by-column description, and `SUPPLEMENTAL_DATA_SURVEY.md` for how the
deposit set was chosen.

`data/` is the provenance tree and is archived here rather than tracked. The bulk data are archived
separately because they total roughly **12.8 GB across ~9,700 files**, which is
far beyond what a git repository should carry.

`.gitignore` excludes them by rule rather than by listing every file, so
dropping the archive back into place restores the tree exactly and nothing has
to be re-pointed.

## What is in the archive

| size | files | path | contents |
|---|---|---|---|
| 7.08 GB | 78 | `data/genotypes/` | CeNDR 20210121 PLINK panel; processed genotype matrix |
| 1.99 GB | 9200 | `data/structure_modeling/` | exploratory AlphaFold runs, HADDOCK docking, electrostatics |
| 1.27 GB | 112 | `data/cross_experiments/N2-XZ_export/` | N2 × XZ1516 xQTL: counts, AFD, all contrast tables |
| 751 MB | 35 | `data/pos1_original/` | 2023 *pos-1* experiment (small trait tables stay in git) |
| 735 MB | 63 | `data/pooled_RNAi_expt/` | pooled competition experiment (small trait tables stay in git) |
| 281 MB | 3 | `data/pooled_cross_intersection/` | cache built by `pooled_cross_intersection_prep.R` |
| 245 MB | 22 | `data/cross_experiments/Nov2024_JU_cross_pos_mig/` | superseded JU cross analysis |
| 134 MB | 36 | `data/cross_experiments/JU1793-JU2466_export/` | JU cross xQTL (one contrast table stays in git) |
| 117 MB | 26 | `data/pos1_rnai_sensitive/` | earlier *pos-1* GEMMA run |
| 82 MB | 70 | `plots/pooled_cross_intersection/` | exploratory figure set |
| 56 MB | 4 | `data/LD/` | LD matrices |
| 31 MB | 5 | `data/baugh/` | Baugh L1 validation inputs (small caches stay in git) |
| 23 MB | 24 | `data/structure/` | AlphaFold bundle (the PDB and confidence JSONs stay in git) |
| 14 MB | 11 | `plots/legacy/` | superseded figures |

Two of these are **regenerable rather than raw**, and it may be better to
document the command than to archive them:

- `data/pooled_cross_intersection/` (281 MB) — rebuild with
  `Rscript scripts/pooled_cross_intersection_prep.R`, given the cross exports.
- `plots/pooled_cross_intersection/` (82 MB) and `plots/legacy/` (14 MB) —
  rebuild by running the scripts in `scripts/` and `scripts/legacy/`.

## What still works from a clone alone

Since `supplemental_data/` is tracked, **all eighteen figures now build from a
clone**, with no archive needed. The table below is about rebuilding the
*inputs* rather than the figures.

| figure | runs from a clone? |
|---|---|
| Figure 1 and its supplements | yes |
| Figure 1 with `FIG1_REFRESH=1` | **no** — needs `data/baugh/2024bootstrapINPUT.Rdata` (31 MB) |
| Figure 2 and cross-QTL supplements | **no** — needs `data/pooled_cross_intersection/bundle.rds`, or the cross exports plus the prep script |
| Figure 3 and all four variants | yes |
| Figure 4 and its supplements | yes |
| 2023 *pos-1* supplements | yes |
| `SUPP_FIG_XX_sid2_allele_in_panel` | **no** — needs `data/genotypes/CeNDR20210121_Plink/` |
| `SUPP_FIG_XX_sid2_electrostatics` | yes |
| `eigen_independent_tests.R` | **no** — needs `data/genotypes/` |
| `baugh_strain_similarity.R` | **no** — needs `data/baugh/2024bootstrapINPUT.Rdata` (31 MB) |
| `baugh_leakage_vs_similarity.R` | yes, on a borrowed predictor — see below |
| plate-phenotype comparisons | yes |

So Figure 2, the eigen-threshold computation, the allele-frequency supplement
and the Baugh similarity table are the four things that need the archive.
Everything else builds from the repository as it stands.

`scripts/baugh_leakage_vs_similarity.R` is the one script that runs from a
clone but runs *better* with the archive. It tests whether the strains the
deconvolution resolves badly are the genetically similar ones, and needs a
relatedness measure to do it. From a clone the only one available is
`dilution_strain_similarity.tsv`, whose distances are to the closest of the 170
dilution strains rather than of the 98-strain Baugh panel — a borrowed
predictor that covers 66 of the 98 and is well posed for 24. Running
`scripts/baugh_strain_similarity.R` once against the archive writes a 102-row
table computed inside the Baugh reference itself; it is small enough to commit,
after which the analysis runs from a clone on the right predictor. The script
prints how much the substitution mattered.

## Kept in git on purpose

These are inputs a figure script reads, small enough to carry, and awkward to
regenerate:

- `data/cross_experiments/JU1793-JU2466_export/plot_data/*HT115g-POS1g*_plot_DF.tsv.gz`
  (15.5 MB) — panel A of every Figure 3 variant
- `data/pooled_RNAi_expt/reanalysis/mapping/*_loco_results.csv.gz` (8 MB) — GWAS panels
- `data/pos1_original/updated_analysis/vst_ctrl_pos-1_T2_loco_results.csv.gz` (7.5 MB)
- `data/baugh/2024_processedBOOTs_with_MIP.RData` (68 KB) — the committed NNLS
  result, which is why Figure 1 builds without the 31 MB input
- `data/baugh/2024baugh_bootstrap_prediction.rda` (3.2 MB) — the bootstrap array
- `data/structure/AF_sid2_wt_dimer/sid2wt.pdb` and the `summary_confidences`
  JSONs; the five CIF models, MSAs and templates are archived
- `data/structure_modeling/.../stage4_electrostatics/sid2_wt_elec.pdb` — the
  membrane-oriented model used by the Figure 4 structure panel
- `plots/assets/` (4.4 MB) — the structure renders, so Figure 4 builds without
  running the Python renderers first

## Restoring the archive

Unpack it over the repository root; the paths in the archive match the paths in
`.gitignore`, so files land where the scripts expect them and git continues to
ignore them.

## To deposit

[TO FILL: Dryad DOI and the date of deposit.]

Suggested exclusions from the deposit, since they are neither raw data nor
results: `data/structure_modeling/` is 1.99 GB of exploratory modelling of
which only one PDB is used and the docking is explicitly unusable (see
`METHODS.txt`), and `data/cross_experiments/Nov2024_JU_cross_pos_mig/` is a
superseded analysis. Dropping both would take the deposit from 12.8 GB to
about 10.6 GB. Depositing the raw counts under
`data/cross_experiments/*/counts/` matters more than either.
