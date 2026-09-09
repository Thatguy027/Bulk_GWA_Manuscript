# Files to keep available across computers

`DATA_AVAILABILITY.md` inventories the whole 12.8 GB archive and is the
reference for a *deposit*. This file is the working document for a *second
machine*: the subset that has to be reachable to keep analysing, ordered so the
small, high-value pieces move first. Sizes for `data/` come from the archive
inventory; the two marked (est.) are computed from marker and sample counts
because the files are not on this machine.

## Tier 1 — needed to define QTL intervals (~120 MB)

Nothing in the repository can produce a linkage-disequilibrium interval, which
is the interval type an association scan in a wild-isolate panel calls for.
These four items are what unblock it, and they are small enough to sync
routinely rather than archive.

| item | size | why |
|---|---|---|
| `data/LD/` | 56 MB, 4 files | already-computed r² of every chrIII marker against the *sid-2* focal variant III:13,680,248, from `scripts/legacy/LD_calc.sh`. Answers directly whether the right-arm association cluster is in LD with *sid-2* |
| `data/genotypes/CeNDR20210121_Plink/III.{bed,bim,fam}` | ~56 MB (est.) | 341,971 chrIII variants × 540 samples. Needed to run `plink --r2 --ld-snp III:12718465` and put an r² ≥ 0.5 interval around the GWA peak itself, not just around *sid-2* |
| `supplemental_data/mapping/pos1_2023_gemma_loco.csv.gz` | in git | already available everywhere |
| `supplemental_data/mapping/ju_cross_ht115_vs_pos1_scan.tsv.gz` | in git | already available everywhere |
| `supplemental_data/mapping/cross_af_*.tsv.gz` | 11.2 MB, in git | genome-wide parental allele counts for every pool of both crosses, so a windowed parental frequency at any position is available from a clone. See below |

The rest of `data/genotypes/CeNDR20210121_Plink/` (all six chromosomes,
~400 MB est.) is only needed if intervals are wanted for QTL off chromosome III.
The full `data/genotypes/` is 7.08 GB; **do not sync it wholesale for this.**

## Tier 2 — needed to re-run analyses that do not build from a clone

`DATA_AVAILABILITY.md` establishes that four things need the archive. These are
their inputs.

| item | size | unblocks |
|---|---|---|
| `data/pooled_cross_intersection/bundle.rds` | part of 281 MB | `compare_full_vs_thinned_bundle.R`, `make_thinned_bundle.R`, `pooled_cross_candidate_variation.R`. **Not** Figure 2 or `SUPP_FIG_XX_cross_contrast_panels.R`: both read the 9.6 MB thinned bundle in git, so neither is blocked on a second machine |
| `data/baugh/2024bootstrapINPUT.Rdata` | 31 MB | Figure 1 with `FIG1_REFRESH=1`; `baugh_strain_similarity.R` |
| `data/genotypes/processed_genotype_matrix.Rda` | part of 7.08 GB | `baugh_leakage_vs_similarity.R` on the right predictor; `make_experiments_deposit.R` |
| `data/genotypes/CeNDR20210121_Plink/` | ~400 MB (est.) | `eigen_independent_tests.R` (the Bonferroni and eigen thresholds); `SUPP_FIG_XX_sid2_allele_in_panel.R`; `simulation_deconvolution.R` (all six chromosomes, via `CENDR_PLINK`) |

Running `scripts/baugh_strain_similarity.R` once against Tier 2 writes a
102-row table small enough to commit, after which that analysis runs from a
clone permanently. Worth doing on whichever machine gets the archive first.

## Tier 3 — needed only to rebuild the deposit

`make_supplemental_data.R` and `make_experiments_deposit.R` regenerate
`supplemental_data/`. They read, beyond Tiers 1–2:

- `data/baugh/` — `2024_processedBOOTs_with_MIP.RData`,
  `2024baugh_bootstrap_prediction.rda`, `2026downsampled_dataset.rda`,
  `cache_boot_freq.rds`, `cache_boot_slopes.rds`, `MIPseq_frequencies.txt`
- `data/cross_experiments/JU1793-JU2466_export/` (134 MB) — the cross exports
- `data/pos1_original/updated_analysis/` — `association_traits.csv`,
  `final_dataset.csv`, `vst_ctrl_pos-1_T2_loco_results.csv.gz`
- `data/pooled_RNAi_expt/` — `reanalysis/vst_association_traits.csv`,
  `reanalysis/mapping/vst_ctrl_mig-6_T2_loco_results.csv.gz`,
  `meta_files/strain_isotype_lookup.tsv`, `paaby2015/emb_leth_data.txt`
- `data/plate_rnai_phenotyping/` — three small phenotype tables
- `data/pos1_plate_phenotyping/pos1_phenotypes_first2rounds.tsv`
- `data/nil_ranges.bed`, `data/eigen_independent_tests.tsv`
- `data/structure/` and one PDB under `data/structure_modeling/`
- `data/experiments/20211111_BulkGWA/`, `data/experiments/initial_sims/`
- `data/experiments/baugh_wgs/` — read as `../data/...` by
  `baugh_L1_DownSample_Counts.R`, i.e. from a sibling checkout, not this tree

Deposit rebuilds are rare. Archive these rather than sync them.

## The cross allele frequencies no longer need the exports

`scripts/make_cross_af_tables.R` reduces both cross exports to two committed
tables, `supplemental_data/mapping/cross_af_{N2xXZ1516,JU1793xJU2466}.tsv.gz`,
11.2 MB together at full marker resolution — 522,357 markers × 10 samples and
153,963 × 4. Everything a frequency question needs is in them:

| what the exports hold | size | what is kept |
|---|---|---|
| `plot_data/*_plot_DF.tsv.gz` | 780 + 45 MB | nothing: each sample's 11 columns are repeated once per contrast it appears in, five times over for the N2 × XZ1516 pools |
| `afd/*.afd.tsv.gz` | 158 + 18 MB | `chrom`, `physical.position`, and each sample's `p1`/`p2` — the counts every other column is derived from |

So the exports are needed only to *rebuild* those tables, which is a one-time
job on the machine that has the archive. `plots/TABLE_cross_qtl_full.tsv` and
`SUPP_FIG_XX_cross_contrast_panels.R` both read the committed tables, not the
exports, and `Figure3_common.R` reads the staged JU scan as before.

## Tier 4 — do not sync

Large, and nothing in the current figure set reads them:
`data/structure_modeling/` (1.99 GB, exploratory docking that `METHODS.txt`
calls unusable), `data/cross_experiments/N2-XZ_export/` (1.27 GB),
`data/cross_experiments/Nov2024_JU_cross_pos_mig/` (245 MB, superseded),
`plots/pooled_cross_intersection/` and `plots/legacy/` (96 MB, regenerable).

## The R environment is not the same on both machines

`METHODS.txt` records the environment the deposited figures were built in. A
second machine will not match it by default, and one of the differences changes
how figures render rather than only which numbers come out:

| package | METHODS.txt (build machine) | this MacBook, 2026-09-08 |
|---|---|---|
| R | 4.5.2 | 4.5.1 |
| ggplot2 | 4.0.3 | **3.5.2** |
| data.table | 1.18.2.1 | 1.17.6 |
| patchwork | 1.3.2 | 1.3.1 |
| ggtext | 0.1.2 | **0.2.0** |
| ggrepel | 0.9.7 | 0.9.6 |
| png | 0.1.8 | 0.1.9 |

The ggplot2 gap is a major version. `scripts/check_repo_invariants.sh` passes
here, but it checks pinned NUMBERS, which are mostly renderer-independent — it
does not compare the rendered figure against the committed PNG. Rebuilding a
figure on this machine and committing the result would therefore quietly mix
two renderers into `plots/`. Rebuild figures on the machine whose versions
match `METHODS.txt`, or update that section deliberately once both agree.

`ggtext` is a hard dependency of `SUPP_FIG_XX_dilution_validation.R`, which
fails to load without it; `chrIII_association_support.R` degrades gracefully via
`requireNamespace`. Install it before running the invariant checks on a new
machine.

## Outside the archive entirely — two hard-coded absolute paths

These point at a location `DATA_AVAILABILITY.md` does not cover, because it is
not under `data/`. They will break on any machine that is not the iMac:

- `scripts/pooled_cross_candidate_variation.R:32`
  `VCF <- "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213/bcsq.vcf.gz"`
- `scripts/pooled_cross_intersection_prep.R:31`
  `CENDR <- "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213"`

The CeNDR 20231213 release is public and re-downloadable, so it does not need
to be synced — but the paths should read an environment variable with a
documented default before either script is run elsewhere. `bcftools` is also
required on PATH by the first of the two.
