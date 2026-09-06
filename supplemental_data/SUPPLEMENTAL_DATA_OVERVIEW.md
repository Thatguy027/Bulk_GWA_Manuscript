# Supplemental data — overview

Every file needed to regenerate all eighteen manuscript figures, and nothing
else. **38.2 MB in 46 files.** Verified by deleting `data/` entirely and
rebuilding: all eighteen figures and the three asset builders run from this
directory alone.

Reproduce with the scripts in `scripts/`, from the repository root:

```sh
Rscript scripts/<figure name>.R          # e.g. scripts/Figure2.R
python3 scripts/sid2_ribbon_render.py    # structure assets, run once
python3 scripts/sid2_zoom_render.py
```

## Two things to know before using these files

**The pooled/cross bundle is thinned.** `mapping/pooled_cross_bundle_thinned.rds`
keeps one marker per 5 kb bin (cross scans) and per 10 kb bin (pooled GWAS) —
the maximum-LOD marker, so every peak survives exactly. It is 9.6 MB instead of
258.7 MB. Peak positions and peak LODs are identical to the full data, but
**QTL support interval widths are 0–9 kb narrower**, because an interval edge is
the outermost marker within a 5% LOD drop of the peak and thinning can discard
it. Widths quoted from this deposit are therefore not interchangeable to the
kilobase with widths from the full bundle. The full bundle is in the Dryad
archive, and `scripts/compare_full_vs_thinned_bundle.R` rebuilds Figure 2 from
both and reports the difference.

**The genotypes are a region subset.** `genotypes/sid2_region.*` covers
chromosome III 13,679,000–13,682,000 only — 81 variants of the 341,971 on that
chromosome — because the two scripts that read genotypes need the *sid-2*
region. The full CeNDR 20210121 panel is a public release, not our data.

## File index

### `phenotypes/`

| file | rows | what it is |
|---|---|---|
| `pooled_vst_traits.csv.gz` | 93 | per-strain pooled RNAi response for ten targets, three parameterisations |
| `pos1_2023_association_traits.csv.gz` | 366 | per-strain 2023 *pos-1* response |
| `pos1_2023_sample_frequencies.csv.gz` | 6,606 | sample-level strain frequencies behind those traits |
| `plate_scores_pos1.tsv` | 191 | manual plate scores on *pos-1* RNAi |
| `paaby2015_embryonic_lethality.txt.gz` | 27,360 | Paaby et al. 2015 well-level counts (published data) |

**`pooled_vst_traits.csv.gz`** — one row per strain, one column set per RNAi
target (`fog-2`, `mig-6`, `pos-1`, `rde-3`, `ric-3`, `rpn-12`, `spe-19`,
`spe-43`, `unc-39`, `vha-5`), all at timepoint 2.

| column | description |
|---|---|
| `strain` | CeNDR isotype name |
| `delta_ctrl_<gene>_T2` | change in pool frequency under RNAi minus control. Positive = gained frequency = resistant |
| `vst_ctrl_<gene>_T2` | the same contrast, variance-stabilised. **This is the trait the association mapping used and the one every figure plots** |
| `log2fc_ctrl_<gene>_T2` | log2 ratio of the same contrast |
| `PC1_delta_ctrl`, `PC1_vst_ctrl` | first principal component across targets, on each scale |
| `negctrl_growth_HT115_delta_t0` | growth on empty-vector control relative to t0; a per-strain quality covariate |
| `negctrl_abundance_log10` | log10 starting abundance in the pool |

Nine of the 93 strains have no `vst` value and are dropped by the figures,
which is why the pooled panel is n = 84.

**`pos1_2023_association_traits.csv.gz`** — as above, for the 2023 *pos-1*
experiment only: `strain`, `delta_ctrl_pos-1_T2`, `vst_ctrl_pos-1_T2`,
`log2fc_ctrl_pos-1_T2`, `negctrl_growth_HT115_delta_t0`,
`negctrl_abundance_log10`. 231 of the 366 rows carry a `vst` value.

**`pos1_2023_sample_frequencies.csv.gz`** — the sample-level data those traits
were reduced from.

| column | description |
|---|---|
| `sample` | sequencing sample identifier |
| `sample_info` | free-text sample label as recorded |
| `time` | timepoint (`t0`, `t2`) |
| `rnai` | food condition: `ctrl` (HT115) or `pos-1` |
| `replicate` | replicate within condition; control A–B, *pos-1* 1–4 |
| `depth_cutoff` | minimum read depth applied when calling frequencies (3, 5 or 10). **Use 5**: it is the cutoff the shipped traits were built from, reproducing `delta_ctrl_pos-1_T2` to 2.5e-16 |
| `strain` | isotype name |
| `frq` | deconvolved frequency of that strain in that sample |
| `t0_frq`, `ctrl_frq` | that strain's frequency at t0, and mean across control replicates |
| `delta_t0`, `delta_ctrl` | `frq` minus each baseline |
| `log2fc_t0`, `log2fc_ctrl` | log2 ratio against each baseline |

CAVEAT: JU1793, and only JU1793, appears **twice per sample** at every depth
cutoff — one row carries the frequency, the other ~1e-21. This is non-negative
least squares splitting one strain across two identical entries in the genotype
reference, not two measurements. Sum rows sharing a `(strain, sample)` key
before taking a delta; the shipped traits average them instead, which halves
JU1793. This matters because JU1793 is a cross parent. See
`scripts/SUPP_FIG_XX_original_pos1_dfreq_rep_correlation.R`.

**`plate_scores_pos1.tsv`** — `strain`, and `trait`: an ordinal plate score,
0 = complete RNAi response (sensitive) to 5 = no response (resistant). Six
levels, so use rank correlation with it.

**`paaby2015_embryonic_lethality.txt.gz`** — published well-level data. The
columns used here are `strain`, `vector` (RNAi clone; filter to `pos-1`),
`eggs` and `larvae` (automated counts). Lethality per well is
`eggs / (eggs + larvae)`, averaged within a strain. `adults_man`, `eggs_man`,
`larvae_man` are the manual counts; `*_area`, `date`, `well_row`, `well_col`,
`image`, `fol_*` are acquisition metadata.

### `hatching_assays/`

The primary hand-scored data. Left as plain text so they can be opened without
tooling. **One plate per strain per condition — no strain is replicated**, so
intervals computed from them describe counting uncertainty, not between-plate
variability.

**`nil_series_hatching.tsv`** (20 rows) — the NIL fine-mapping assay at
**50% *pos-1* RNAi**.

| column | description |
|---|---|
| `Strain` | strain name; parents JU1793 and JU2466, NILs `wSZ*` |
| `condition` | `ht115` (empty vector) or `pos1` |
| `plated embryo` | embryos plated |
| `unhatched` | embryos that failed to hatch |
| `fraction` | fraction unhatched as recorded |

**`ju_allele_swaps_hatching.csv`** (14 rows) — reciprocal *sid-2* residue-96
edits in JU1793 and JU2466, at **50% *pos-1* RNAi**.

| column | description |
|---|---|
| `experiment` | experiment identifier |
| `strain` | strain assayed; `wSZ200`, `wSZ206`, `wSZ207`, `wSZ208` are edits |
| `genotype` | background and residue-96 allele, e.g. `JU1793[96T]` |
| `glycosylation motif` | the N94-x-96 sequon state, e.g. `JU1793[NxT]` |
| `condition` | `pos` or `ht115` |
| `n_plated`, `n_unhatched`, `fraction_hatched` | counts and the hatched fraction |

**`n2_allele_swaps_hatching.tsv`** (15 rows) — the same substitution in N2
across a dose series. `condition` here is the **percentage of *pos-1* RNAi
bacteria** in the lawn (0, 25, 50, 75, 100), not a category. Only the 25% dose
has dynamic range: at 0% every genotype hatches, from 50% up every genotype is
dead. `wSZ203` and `wSZ204` are independent lines of the same edit — the
closest thing to a biological replicate in these assays.

**`nil_introgression_ranges.bed`** — headerless BED: chromosome, start, end,
strain, genotype (which parent the interval came from). The interval the series
resolves is chromosome III 13,657,700–13,695,000.

### `mapping/`

| file | rows | what it is |
|---|---|---|
| `pooled_cross_bundle_thinned.rds` | — | thinned pooled GWAS + 18 cross contrast scans, plus intervals and metadata |
| `ju_cross_ht115_vs_pos1_scan.tsv.gz` | 153,963 | full-resolution JU1793 × JU2466 HT115-vs-*pos-1* scan |
| `pos1_2023_gemma_loco.csv.gz` | 464,045 | GEMMA LOCO association for `vst_ctrl_pos-1_T2` |
| `eigen_independent_tests.tsv` | 14 | effective test counts and significance thresholds |

**`pooled_cross_bundle_thinned.rds`** — an R list. Read with `readRDS()`.

| element | description |
|---|---|
| `gwas` | pooled GWAS, thinned to one marker per 10 kb: `gene`, `chr`, `rs`, `ps`, `af`, `beta`, `se`, `p_wald`, `neglog10p` |
| `gwas_peaks`, `gwas_sig` | per-trait peak markers, and those clearing significance |
| `GWAS_BF`, `N_MARKER` | Bonferroni threshold and marker count for the pooled panel |
| `scans` | named list of 18 cross contrast scans, thinned to one marker per 5 kb. Names are `<cross> \| <condA> vs <condB>`. Columns: `chrom`, `physical.position`, `contrast.beta`, `z`, `LOD`, `neglog10p` |
| `scan_meta` | which cross and conditions each scan name refers to |
| `ivs` | QTL support intervals per chromosome per scan, **at full resolution** — use these rather than recomputing from the thinned scans if you need published interval bounds |
| `CROSS_THR` | genome-wide LOD threshold for the cross scans (3.57) |
| `CROSSES` | the two crosses and their parents |
| `genes`, `rnai_genes` | WBcel235 gene coordinates, and the RNAi targets among them |
| `panel_gt`, `parent_gt_qtl`, `inform` | panel and parental genotypes at QTL, and cross informativeness |
| `isect_A`, `isect_B`, `isect_C` | GWAS-to-cross-interval intersections at three stringencies |
| `thinned` | the thinning parameters and before/after row counts |
| `built` | when the bundle was made |

**`ju_cross_ht115_vs_pos1_scan.tsv.gz`** — the scan behind panel B of
Figure 3, at full marker resolution. `chrom`, `physical.position`, `ID`, then
per-sample blocks suffixed `_S32_..._HT115g` and `_S34_..._POS1g`:

| column stem | description |
|---|---|
| `ref_`, `alt_` | raw reference and alternate allele counts |
| `p1_`, `p2_` | counts assigned to the JU1793 and JU2466 haplotypes after phasing |
| `nindv_`, `ndepth_`, `n_` | contributing individuals, read depth, total counts |
| `afd_`, `afd.se_` | allele-frequency deviation and its standard error |

then the contrast: `contrast.beta` (effect, positive = the phased parent's
allele rose on the left side), `contrast.beta.se`, `z`, `p`, `LOD`, `log10p`,
`neglog10p`.

**Use `LOD` or `neglog10p`, never `p`.** `p` is formed on the linear scale and
underflows to exactly 0 above |z| ≈ 38.5, which happens often here. `LOD` is
computed in log space and cannot saturate.

**`pos1_2023_gemma_loco.csv.gz`** — standard GEMMA output: `chr`, `rs`, `ps`
(position), `n_miss`, `allele1` (alternate), `allele0` (reference), `af`
(alternate frequency), `beta`, `se`, `logl_H1`, `l_remle`, `p_wald`, `trait`.

**`eigen_independent_tests.tsv`** — one row per panel per scope.

| column | description |
|---|---|
| `panel` | `pooled_RNAi_expt` or `pos1_2023` |
| `chrom`, `scope` | chromosome, and whether the row is genome-wide or per chromosome |
| `n_strain`, `n_tested`, `n_marker` | strains, markers tested, markers in the panel |
| `frac_calls_imputed` | fraction of genotype calls mean-imputed before the eigendecomposition |
| `rank_nonzero` | non-zero eigenvalues of the marker correlation matrix |
| `M_eff_liji` | effective independent tests, Li & Ji (2005) |
| `M_eff_var995` | alternative: components explaining 99.5% of variance |
| `trace_err` | relative error on `trace(A) = M`, a correctness check |
| `thr_bonferroni` | −log10 p threshold, Bonferroni over every marker |
| `thr_eigen_liji`, `thr_eigen_var995` | −log10 p thresholds from each effective count |

### `deconvolution/`

Validation of the NNLS strain-frequency deconvolution. Three experiments, in
the order the manuscript reports them — a simulation with a known input and
synthetic counts, a designed DNA mixture with a known input and real counts,
and the Baugh L1 starvation time course compared against published MIP-seq.

#### The simulation (`simulation_*`)

Seven published *C. elegans* traits with validated QTL each seeded a simulated
pooled population of 327 strains; counts were simulated by binomial sampling at
eight depths and deconvolved back by NNLS.

| file | what it is |
|---|---|
| `simulation_nnls_frequencies.tsv.gz` | the archived NNLS output, 18,312 rows = 7 traits × 8 depths × 327 strains |
| `simulation_reported_r2.tsv` | 56 rows: the r² of estimate against **known input**, per trait and depth |
| `simulation_gwas_traits.tsv.gz` | 327 strains × 56 trait-depth columns — the simulated frequencies as a GWAS trait file |

`simulation_nnls_frequencies.tsv.gz` columns: `trait` (which trait seeded the
population), `depth` (simulated coverage, 1–500), `strain` (isotype),
`coefficient` (the archived NNLS coefficient — **not** a frequency: each
trait-and-depth's coefficients sum to the depth itself, so they are on an
expected-count scale), `frequency` (`coefficient` divided by that sum, which is
the quantity to use).

`simulation_reported_r2.tsv` columns: `trait`, `depth`, `r2`. **Read the
provenance before quoting these.** The expected input frequencies and the
simulation script were never archived, so this r² cannot be recomputed; the
values are recovered from text embedded in the original per-trait PDFs by
`scripts/extract_sim_reported_r2.py`, which is the only surviving record of the
comparison against the known input. Everything else in the deposit is
recomputable; this one file is not.

Two caveats live with these data. **139 of the 18,312 coefficients (0.8%) are
negative**, at depths 5 through 100 — a strict non-negative solver cannot return
those, so the archive carries noise of order 1e-2 on the count scale. And at 1×
the reported r² spans **0.52 to 0.91** across the seven traits (median 0.79), so
"accurate at 1×" is trait-dependent; 30× is the lowest depth at which all seven
reach 0.95.

#### The designed DNA mixture (`dilution_*`)

174 wild isolates in four sets (A–D); each set pooled and sequenced pure in
triplicate, plus a seven-step titration of set B against set C (BC1–BC7).

| file | what it is |
|---|---|
| `dilution_strain_sets.tsv` | 174 rows: which set each strain went into |
| `dilution_predictions_poolref.tsv.gz` | 170 strains × 19 samples, reference restricted to the pooled strains |
| `dilution_predictions_bcref.tsv.gz` | 84 strains × 7 samples, reference restricted to sets B+C — the original analysis |
| `dilution_predictions_fullref.tsv.gz` | 540 strains × 19 samples, full CeNDR reference |
| `dilution_predictions_regenotype.tsv.gz` | as above, from the regenotyped VCF |
| `dilution_strain_similarity.tsv` | 170 rows: how genetically close each pooled strain is to its nearest neighbour |

`dilution_strain_sets.tsv` columns: `strain` (name as pooled), `isotype` (CeNDR
isotype; `NA` for ECA252 and LSJ1, which have no isotype in the lookup and are
dropped by the analysis), `set` (A, B, C or D). Kept at **strain** level on
purpose — see the ambiguity note below.

Prediction files share the columns `strain` (isotype), `sample`, `frequency`,
and the three wide-reference files also carry `set`. Frequencies sum to 1 within
a sample.

**One isotype is assigned to two sets, and it matters.** Strains JU1580 and
JU1793 are the same isotype (JU1793) and the experiment put them in sets B and
D. The deconvolution returns one frequency for the pair, so
`dilution_predictions_poolref.tsv.gz` carries that isotype **twice**, once
labelled `B` and once labelled `D` — 171 rows per sample for 170 distinct
strains. Summing by `set` without resolving this adds the single estimate to
both, where it is 5–19% of set D's apparent frequency in BC1–BC7 and *rises
across the titration*, making the untitrated set D look like it drifts.
`scripts/SUPP_FIG_XX_dilution_validation.R` resolves it to set B — what the
original analysis did, and what the estimate's own behaviour indicates — and
asserts that no other isotype is ambiguous. This is the same strain whose
duplicated genotype-reference entry is documented under `phenotypes/`, and it is
a cross parent in Figures 2 and 3.

`dilution_strain_similarity.tsv` columns: `strain` (isotype), `nn_ibs`
(identity-by-state — the fraction of markers carrying the same allele — to the
closest *other* strain in the reference), `nn_partner` (which strain that is),
`mean_ibs` (mean identity to all other strains). Computed over 231,615 variable
markers sampled under a fixed seed from the 12 GB CeNDR genotype matrix, which
is Dryad-hosted; this 170-row reduction is what the figure reads, so the
deposit-only rule holds. Used by panel D of the dilution figure, where leakage
into pools a strain is absent from rises with `nn_ibs` (Spearman ρ = 0.326,
p = 1.4e-05) while own-pool recovery does not (ρ = 0.046, p = 0.55).

**The nominal mixing ratios are not recorded.** Nothing in the archived
experiment states the intended B:C proportions for BC1–BC7, so these data show
monotonic and complementary recovery, not quantitative accuracy.

Both groups are rebuilt by `scripts/make_experiments_deposit.R` from the raw
folders under `data/experiments/` (Dryad-hosted), except
`simulation_reported_r2.tsv`, which comes from the Python extractor above.

#### MIP-seq comparison (`baugh_*`)

Validation against published MIP-seq, using the Baugh L1 starvation time
course.

| file | what it is |
|---|---|
| `baugh_nnls_with_mipseq.RData` | the deconvolution result joined to MIP-seq — `wgs_mip_results`, 2,346 rows |
| `baugh_bootstrap_array.rda` | 100 bootstrap replicates of the deconvolution |
| `baugh_downsampled_slopes.rda` | the same deconvolution at six sequencing depths |
| `baugh_strain_order.txt` | 102 strain names, in the bootstrap array's row order |
| `mipseq_frequencies.txt.gz` | published MIP-seq frequencies |
| `cache_boot_slopes.rds`, `cache_boot_freq.rds` | precomputed bootstrap reductions |

**`baugh_nnls_with_mipseq.RData`** — `load()` gives `wgs_mip_results`:
`sample`, `replicate`, `day`, `baseline` (is this the day-1 baseline sample),
`strain`, `bootmean` and `bootse` (mean and SD across bootstrap replicates),
`frq` (the point-estimate NNLS frequency) and `published_frq` (MIP-seq).

**`baugh_bootstrap_array.rda`** — `ab`, a 102 × 23 × 100 array (strain ×
sample × bootstrap replicate) of deconvolved frequencies, normalised to sum to
1 within each sample; `bootse`, the 102 × 23 standard errors; `blist`, the same
content reshaped. **The first dimension is unnamed** — its order is
`baugh_strain_order.txt`. Only the outputs were archived, never the resampling
code, so what was resampled (reads, markers or strains) is not recoverable.

**`baugh_downsampled_slopes.rda`** — `ds_predictions_df`: `strain`, `sample`,
`ds_frq` (frequency at the subsampled depth), `ds_n` (depth: 0.25, 0.5, 1, 3, 5
or 10×), `frq` (full-depth frequency).

**`baugh_strain_order.txt`** — one column, `strain`, with a header line.

**`mipseq_frequencies.txt.gz`** — strains in the first column (unnamed), then
one column per sample: `rep<N>_d<day>` for the time course and `rep<N>_BL` for
the baselines, which the scripts rename to `rep<N>_d1_baseline`.

**`cache_boot_slopes.rds`, `cache_boot_freq.rds`** — precomputed reductions of
the bootstrap array: a 102 × 100 matrix of per-strain slopes, and per-sample
frequency percentile intervals. Shipped so the figures build quickly; delete
them and they rebuild from the array.

### `genotypes/`

`sid2_region.{bed,bim,fam}` — PLINK binary set, chromosome III
13,679,000–13,682,000, 81 variants across 540 CeNDR isotypes. Standard PLINK
format: `.bed` is binary genotypes, `.bim` is chromosome / variant ID /
centimorgan / position / allele 1 / allele 2 (headerless), `.fam` is family ID /
individual ID / paternal / maternal / sex / phenotype (headerless). `.log` is
the PLINK log from the subsetting.

The variant of interest is **`III:13680248`**, C>A: REF C encodes 96T and ALT A
encodes 96K. Orientation was verified against four strains whose allele is
known independently — JU1793 and N2 are 96T, JU2466 and XZ1516 are 96K.

### `structure/`

| file | what it is |
|---|---|
| `sid2_alphafold_dimer.pdb` | AlphaFold3 SID-2 dimer, top-ranked model, arbitrary frame |
| `sid2_membrane_oriented.pdb` | the same model rotated onto the membrane normal: z = 0 is the bilayer centre, +z the intestinal lumen. Residues 21–188, both chains |
| `sid2_alphafold_confidences/*.json` | ipTM, pTM, ranking score and disorder fraction for all five models |
| `sid2_deeptmhmm_topology.3line` | DeepTMHMM output: header, sequence, per-residue topology code |
| `sid2_per_residue.tsv` | per-residue annotation derived from the model |
| `sid2_variants_cendr.tsv` | *sid-2* protein-altering variants with population frequencies |
| `sid2_variant_ld.tsv` | pairwise LD among those variants |
| `sid2_parental_variants.tsv` | annotated variants distinguishing the cross parents |

**Model confidence limits use.** ipTM is 0.45–0.46 and pTM 0.46–0.48 across all
five models, with 43–46% called disordered, so the **dimer interface is not
credible** and is not shown in any figure. Mean pLDDT by region for chain A:
signal peptide 56.4, extracellular 72.3, transmembrane helix 46.0, cytoplasmic
40.0. Only the extracellular domain reaches the confident band.

**`sid2_deeptmhmm_topology.3line`** — three lines: `>SID2 | SP+TM`, the 311-aa
sequence, then one topology character per residue: `S` signal peptide,
`O` extracellular, `M` transmembrane, `I` cytoplasmic. Gives signal 1–20,
extracellular 21–193, TM helix 194–211, cytoplasmic 212–311.

**`sid2_per_residue.tsv`** (311 rows)

| column | description |
|---|---|
| `resid` | residue number, 1–311 |
| `resname` | three-letter amino acid |
| `plddt` | AlphaFold per-residue confidence, mean over the residue's atoms |
| `sse` | secondary structure: `H` helix, `E` strand, `C` coil. Assigned from the backbone with the Kabsch & Sander (1983) hydrogen-bond criterion, not DSSP |
| `topology` | DeepTMHMM region name |
| `x`, `y`, `z` | Cα coordinates, in the model's own frame |

**`sid2_variants_cendr.tsv`** (8 rows) — provenance and the coverage caveat are
in `sid2_variants_cendr_README.txt` beside it.

| column | description |
|---|---|
| `chrom`, `pos` | chromosome and position |
| `ref`, `alt` | reference and alternate allele |
| `consequence` | BCSQ consequence; all eight rows are `missense` |
| `aa_change` | BCSQ notation, e.g. `96T>96K`. Multiple transcripts are `;`-separated |
| `residue` | protein residue number |
| `label` | display form, e.g. `T96K`; alternates joined with `/` |
| `af` | alternate allele frequency across the CeNDR isotype set |
| `n_alleles`, `n_isotypes` | alleles and isotypes observed at that site |
| `ju1793_aa`, `ju2466_aa`, `n2_aa`, `xz1516_aa` | the amino acid each strain carries |
| `parents_differ` | do JU1793 and JU2466 differ here? TRUE for exactly two of the eight, `V5L` and `T96K` |

**NOT an exhaustive catalogue.** 81 variants of any kind segregate in the 3 kb
*sid-2* span in CeNDR; 11 are annotated here and 8 are protein-altering. The
annotation covers variants segregating among the four cross parents, so a
variant private to another isolate is absent.

**`sid2_variant_ld.tsv`** (28 rows) — `variant_a`, `variant_b` (labels as
above), `r` (Pearson correlation of alternate-allele dosage across CeNDR
isotypes) and `r2`. Note **T96K and P153T are r² = 0.935**: 96K never occurs
without 153T, though 153T occurs without 96K in nine isotypes — JU1793 among
them, which is why the JU1793 × JU2466 cross can isolate residue 96.

**`sid2_parental_variants.tsv`** (30 rows) — `cross`, `pos`, `chrom`, `ref`,
`alt`, `gt_parent1` / `gt_parent2` (VCF genotypes), `parent1` / `parent2`
(names), `alt_carrier`, `gene`, `consequence`, `aa_change`,
`protein_altering`, `in_sid2`, `focal` (is this the T96K variant).

## Provenance

Files here that **no script in the repository creates** — pre-processed
elsewhere or externally sourced — are listed in `../SUPPLEMENTAL_DATA_SURVEY.md`
with the reason for each. In short: the pooled and 2023 *pos-1* trait tables
come from an analysis folder outside this repository and the
variance-stabilising transform is not reproducible from anything here; the
Baugh inputs, the CeNDR panel, the AlphaFold and DeepTMHMM outputs and the
Paaby data are external; and the hatching assays and NIL breakpoints are
primary measurements that nothing can regenerate.
