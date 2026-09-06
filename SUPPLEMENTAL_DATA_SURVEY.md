# Supplemental data survey

What it would cost to stage a self-contained supplemental data folder that
regenerates all sixteen figures, and which of those files have no recipe in
this repository.

Survey only — nothing has been copied or moved.

> **Historical.** This is the survey as written *before* the deposit was
> staged, and its numbers describe that moment: sixteen figures, 33 files,
> 376 MB raw. The deposit was subsequently built along the lines recommended
> at the end, and the two deconvolution-validation supplements (Figures S1 and
> S2) were added afterwards, bringing their own inputs with them. For the
> deposit as it now stands — **eighteen figures, 46 files, 38.2 MB** — see
> `supplemental_data/SUPPLEMENTAL_DATA_OVERVIEW.md`. The numbers below are
> left as they were, because rewriting them piecemeal would make the
> cost argument they support incoherent.

## Headline

**33 files, 376 MB raw, 324 MB gzipped.** Since these files already exist under
`data/`, staging a compressed copy costs about **324 MB of extra disk** on top
of the current tree.

Compression buys little (14%) because most of the set is *already* compressed:
`.rds`, `.rda`, `.RData`, `.csv.gz` and `.tsv.gz` account for 316 of the
376 MB. The plain-text files compress 46–90%, but they are small.

One file dominates everything:

| file | raw | gzipped | share |
|---|---|---|---|
| `data/pooled_cross_intersection/bundle.rds` | 258.7 MB | 258.7 MB | **80%** |
| `data/baugh/2024bootstrapINPUT.Rdata` | 31.0 MB | 31.0 MB | 10% |
| JU cross contrast scan (`…HT115g-POS1g…tsv.gz`) | 15.5 MB | 15.5 MB | 5% |
| `data/genotypes/CeNDR20210121_Plink/III.bed` + `.bim` | 56.5 MB | 6.8 MB | 2% |
| 2023 *pos-1* GEMMA scan (`…loco_results.csv.gz`) | 7.5 MB | 7.5 MB | 2% |
| Baugh bootstrap array (`…bootstrap_prediction.rda`) | 3.2 MB | 3.2 MB | 1% |
| everything else (27 files) | 3.8 MB | 1.2 MB | <1% |

## Four ways to make it much smaller

**1. Thin the bundle — the big one.** `bundle.rds` is 259 MB of which the
`scans` element is essentially all of it: 18 contrast tables, 8,297,244 rows of
per-marker LOD. Every figure that reads it thins to 5 kb bins for display, and
at 5 kb only **342,768 rows survive — 4.1% of the original**. A pre-thinned
bundle would be roughly **11 MB instead of 259**, taking the whole deposit to
about **45 MB**.

  The caveat, and it is a real one: `Figure2.R` recomputes QTL support
  intervals from the full-resolution scans (5% LOD drop). Thinning first would
  shift interval bounds slightly. The fix is to ship the intervals that the
  bundle *already contains* (`ivs`, 108 × 13) and have `Figure2.R` read those
  instead of recomputing — a small change, and arguably better practice for a
  deposit, since the published intervals then cannot drift from the figure.

**2. Drop the `FIG1_REFRESH` input.** `2024bootstrapINPUT.Rdata` (31 MB) is
needed only to redo the NNLS deconvolution from scratch. Every figure builds
from the 68 KB committed result. Saves 31 MB and costs only the ability to
re-run the slow path.

**3. Subset the genotypes.** `III.bed` + `.bim` are 56.5 MB raw because they
carry all 341,971 chromosome III markers, and the two scripts that read them
want **one variant**, `III:13680248`. A region subset would be a few kB. Note
PLINK needs these uncompressed, so the 6.8 MB gzipped figure is misleading —
staged uncompressed they cost 56.5 MB.

**4. Recompress the `.rds`/`.rda` with xz.** R's default is gzip; `xz` typically
gains a further 15–30% on this kind of numeric data. Worth doing only if the
bundle is kept at full resolution.

Taking 1–3 together: **about 45 MB**, from 324 MB.

## Files with no recipe in this repository

You asked which files were pre-processed elsewhere. Of the 33, **11 are
produced by scripts here** and **22 are not**.

### Produced by scripts in this repository

| file | script |
|---|---|
| `data/pooled_cross_intersection/bundle.rds` | `pooled_cross_intersection_prep.R` |
| `data/baugh/2026downsampled_dataset.rda` | `baugh_L1_DownSample_Counts.R` |
| `data/baugh/2024_processedBOOTs_with_MIP.RData` | `Figure1_common.R` (`baugh_frequencies(refresh = TRUE)`) |
| `data/structure/sid2_per_residue.tsv` | `sid2_ribbon_render.py` |
| `data/structure/sid2_variants_cendr.tsv` | `sid2_variant_table.R` |
| `data/structure/sid2_variant_ld.tsv` | `sid2_variant_table.R` |
| `data/eigen_independent_tests.tsv` | `eigen_independent_tests.R` |
| `plots/pooled_cross_intersection/TABLE_sid2_parental_variants.tsv` | `pooled_cross_candidate_variation.R` |
| `plots/assets/*.png` (3 renders) | `sid2_ribbon_render.py`, `sid2_zoom_render.py` |

### No recipe here — pre-processed elsewhere or externally sourced

**Analysis done in another folder.** The upstream pipeline is not in this
repository, so these arrive as finished products:

- `data/pooled_RNAi_expt/reanalysis/vst_association_traits.csv` — the pooled
  phenotypes. The variance-stabilising transform is *not* reproducible from
  anything here; this is already flagged in `METHODS.txt`.
- `data/pos1_original/updated_analysis/final_dataset.csv`
- `data/pos1_original/updated_analysis/association_traits.csv`
- `data/pos1_original/updated_analysis/vst_ctrl_pos-1_T2_loco_results.csv.gz`
  — all three from `/UCLA/Projects/bulkGWAS/lipid_RNAi/2023_pos1`
- `data/baugh/2024bootstrapINPUT.Rdata` — genotype and allele-count matrices
- `data/baugh/2024baugh_bootstrap_prediction.rda` — the bootstrap array. Only
  its *outputs* were archived, never the resampling code, which is why
  `METHODS.txt` cannot say what was resampled.

**Reproducible, but by a pipeline shipped alongside the data rather than in
`scripts/`:**

- `data/cross_experiments/JU1793-JU2466_export/plot_data/…_plot_DF.tsv.gz` —
  the export carries its own `scripts/run_all.sh`, `xqtl_pipeline.R` and
  `config.yml`, and its README states that regenerating from
  `counts/aser_combined/` reproduces every peak and interval exactly, with LOD
  agreeing to 1e-14. Those inputs are 1.4 GB and are not in the deposit set.

**Externally sourced — no recipe possible:**

- `data/pooled_RNAi_expt/paaby2015/emb_leth_data.txt` — published data
- `data/baugh/MIPseq_frequencies.txt` — published data
- `data/genotypes/CeNDR20210121_Plink/III.{bed,bim,fam}` — CeNDR 20210121 release
- `data/structure/AF_sid2_wt_dimer/sid2wt.pdb` and the five
  `summary_confidences_*.json` — AlphaFold3 output
- `data/structure/SID2deeptmhmm/predicted_topologies.3line` — DeepTMHMM output
- `data/structure_modeling/claude_docking/stage4_electrostatics/sid2_wt_elec.pdb`
  — membrane-oriented model. A recipe *does* exist, in
  `data/structure_modeling/claude_docking/stage3_membrane/*.py`, but it is not
  in `scripts/`.

**Measured by hand — the primary data:**

- `data/nil_ranges.bed` — NIL introgression breakpoints
- `data/plate_rnai_phenotyping/RNAi Sensitivity - JU1793 - November2025 gSZ177 - gSZ179 NILs.tsv`
- `data/plate_rnai_phenotyping/20260409_ju2466swap_plus_N2A_swap.csv`
- `data/plate_rnai_phenotyping/n2_swap.tsv`
- `data/pos1_plate_phenotyping/pos1_phenotypes_first2rounds.tsv`

These five are irreplaceable — they are the hatching counts and plate scores
themselves, and nothing in the repository or any archive can regenerate them.
They total under 20 kB.

## What I would do

Stage the deposit as `supplemental_data/` with the thinned bundle, no
`FIG1_REFRESH` input and a region-subset genotype file: about **45 MB**, and
every figure still regenerates. Keep the full-resolution bundle and the cross
exports in the Dryad archive described in `DATA_AVAILABILITY.md`, so the
deposit is the reproduction set and Dryad is the provenance set.

Awaiting a decision on the bundle before restructuring the scripts.
