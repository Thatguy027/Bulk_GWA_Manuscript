# Per-paragraph review of the manuscript draft

One file per drafted paragraph. Each records what the repository's own record says about the
numbers, claims and figure citations in that paragraph, the recomputed values, and a rewritten
version of the paragraph. Written 2026-09-09.

These are **working notes, not a deliverable**. They exist here rather than in a chat log so
they travel between machines with the repository.

## The paragraphs

| file | paragraph | figure it leans on |
|---|---|---|
| `para1_simulation.md` | NNLS recovers simulated strain frequencies as a function of depth | `SUPP_FIG_XX_simulation_depth` |
| `para2_dilution.md` | the designed DNA mixture, sets B and C | `SUPP_FIG_XX_dilution_validation` |
| `para3_mipseq.md` | comparison against published MIP-seq (Webster et al. 2022) | `Figure1_pos1` panel A, `SUPP_FIG_XX_downsample_per_sample` |
| `para4_pos1_pilot.md` | the pooled *pos-1* pilot, plate re-evaluation, and the GWA | `Figure1_pos1` panels B and C, `SUPP_FIG_XX_original_pos1_dfreq_rep_correlation`, `SUPP_FIG_plate_vs_paaby_vs_pos1original` |

## Recomputed tables

Everything under `data/` was computed from `supplemental_data/` on this repository, not
transcribed from a caption. Where a figure script pins its plotted values as literals, the
recomputation reproduces those literals exactly and the file says so.

| file | what it holds | computed from |
|---|---|---|
| `simulation_trait_strain_counts.csv` | per-trait non-NA strain counts, both deposit files side by side | `deconvolution/simulation_gwas_traits.tsv.gz`, `simulation_nnls_frequencies.tsv.gz`, `simulation_reported_r2.tsv` |
| `panelB_r2_padded_vs_scanned.csv` | pooled r² against the 500x estimate, full grid vs non-NA pairs only | `deconvolution/simulation_nnls_frequencies.tsv.gz` |
| `dilution_design_vs_recovered.csv` | designed vs recovered set-B fraction per titration step | `deconvolution/dilution_predictions_bcref.tsv.gz`, `dilution_design.tsv`, `dilution_strain_sets.tsv` |
| `mipseq_depth_series.csv` | slope-level and per-sample agreement with MIP-seq by depth | `FIGURE_CAPTIONS.txt` / `METHODS.txt` (the downsampled slopes themselves are in `deconvolution/baugh_downsampled_slopes.rda`) |
| `pos1_replicate_correlations.csv` | all six pairwise Spearman rho between the four *pos-1* replicates | `phenotypes/pos1_2023_sample_frequencies.csv.gz` |
| `pos1_response_counts.csv` | responsive-strain counts under five different definitions | `phenotypes/pos1_2023_sample_frequencies.csv.gz`, `pos1_2023_association_traits.csv.gz` |
| `pos1_gwa_bonferroni_markers.csv` | the ten markers clearing Bonferroni, with effect sizes | `mapping/pos1_2023_gemma_loco.csv.gz` |
| `pos1_gwa_eigen_summary.csv` | per-chromosome eigen-threshold counts and spans | `mapping/pos1_2023_gemma_loco.csv.gz` |

## Relationship to the number checker

`scripts/check_manuscript_numbers.py` attests draft numbers against `FIGURE_REPORT.md` and
nothing else. Its default globs are `manuscript/*.txt`, `manuscript/*.md` and `MANUSCRIPT.md`,
which do **not** reach this subdirectory — so these files are not scanned, and they should not
be. They deliberately quote values that are absent from the report (that is often the finding).

If the pre-push hook is ever widened to `manuscript/**`, exclude `manuscript/review/` rather
than adding these numbers to `manuscript_number_exceptions.txt`.

Three findings in these files are about the checker itself and are worth keeping in view: it
passed a paragraph claiming an "average discrepancy of 2.85%" because 2.85 appears in the report
as an unrelated per-mille leakage value; it passed "~48 strains" because 48 appears as a
chromosome IV strain count; and it passed a broad-sense heritability of 0.32 because 0.32 rounds
onto the leakage Spearman rho of 0.326. It checks drift, not meaning.

## Open items the repository cannot settle

- **The *pos-1* pool size.** The draft says 224 strains. That is not reproducible from the
  deposit, and `pos1_2023_association_traits.csv.gz` carries 231 strains with a response value —
  seven more than the stated pool. Needs the pool composition list from the lab record.
- **Broad-sense heritability.** No estimate exists anywhere in this repository.
- **QTL intervals.** No linkage-disequilibrium interval can be computed here; see
  `SYNC_MANIFEST.md`, Tier 1.
- **The simulation's provenance.** ~~Per `METHODS.txt`, the simulation script, drawn fitness
  values and expected input frequencies were never archived~~ — **partly settled.** The script
  was recovered and is archived as `scripts/simulation_deconvolution.R`, with the original
  working file at `scripts/legacy/haploReg_original.R`; it pins the inverse chi-squared
  parameters, the fitness-to-frequency mapping and the bootstrap unit. Two things remain open
  and are marked `[TO FILL]` in `METHODS.txt`: the fitness draw was never seeded, so the exact
  simulated populations cannot be reproduced, and the seven-trait simulation used the trait
  values themselves as fitness rather than a draw, which is not what that paragraph currently
  says.
