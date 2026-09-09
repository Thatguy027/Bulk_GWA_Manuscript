# Review of manuscript paragraph 1 — NNLS simulation

Source of truth: `supplemental_data/deconvolution/simulation_reported_r2.tsv`,
`FIGURE_CAPTIONS.txt` (SUPP_FIG_XX_simulation_depth section), `METHODS.txt` (Simulation).
All values below recomputed from the TSV, not transcribed.

## 1. The 1x claim is not supported

At 1x, **no trait reaches r2 >= 0.95 (0 of 7)**. Spread: 0.52 (mtDNA_ratio) to
0.91 (PC1), median 0.79. `FIGURE_CAPTIONS.txt` flags the sentence
"NNLS can accurately infer strain frequencies with as little as 1X sequencing depth"
verbatim as CAVEAT 2.

Traits with r2 >= 0.95 by depth: 1x: 0/7, 3x: 2/7, 5x: 3/7, 10x: 5/7, 30x: 7/7, 50x: 7/7, 100x: 7/7, 500x: 7/7.
30x is the lowest depth at which all seven sit at or above 0.95.
All seven reach r2 = 1.00 at 500x.

### Attested r2, estimated vs known input (panel A)

| trait | 1x | 3x | 5x | 10x | 30x | 50x | 100x | 500x |
|---|---|---|---|---|---|---|---|---|
| Albendazole_q75.TOF | 0.56 | 0.77 | 0.84 | 0.90 | 0.95 | 0.97 | 0.97 | 1.00 |
| PC1 | 0.91 | 0.98 | 0.98 | 0.99 | 1.00 | 1.00 | 1.00 | 1.00 |
| amsacrine_f.L1 | 0.79 | 0.93 | 0.95 | 0.97 | 0.99 | 1.00 | 1.00 | 1.00 |
| assay_norm | 0.81 | 0.89 | 0.93 | 0.96 | 0.99 | 1.00 | 1.00 | 1.00 |
| etoposide_median.TOF | 0.72 | 0.89 | 0.93 | 0.97 | 0.99 | 1.00 | 1.00 | 1.00 |
| mtDNA_ratio | 0.52 | 0.69 | 0.78 | 0.89 | 0.95 | 0.96 | 0.99 | 1.00 |
| value | 0.86 | 0.96 | 0.98 | 0.99 | 0.99 | 1.00 | 1.00 | 1.00 |

Lowest depth from which each trait stays >= 0.95: PC1 3x, value 3x, amsacrine_f.L1 5x,
assay_norm 10x, etoposide_median.TOF 10x, Albendazole_q75.TOF 30x, mtDNA_ratio 30x.

## 2. Depths are eight discrete values, not a range

Binomial sampling was done at 1, 3, 5, 10, 30, 50, 100 and 500x. "A range of sequencing
depths (1-500x)" reads as a continuous sweep. List the eight values, in Methods at minimum.

## 3. The paragraph is missing its denominator

Seven independent simulated populations of 327 strains; fitness structure supplied by seven
published C. elegans traits with validated QTL (Albendazole_q75.TOF, PC1, assay_norm,
mtDNA_ratio, value, amsacrine_f.L1, etoposide_median.TOF). Without this, "accurately" has no
denominator and the r2 spread across traits is invisible to the reader.

## 4. Figure citation: cite the panel, not the figure

Target is `SUPP_FIG_XX_simulation_depth`, the earliest supplement (S1 if the curated set is
numbered in `FIGURE_CAPTIONS.txt` order; the `XX` placeholders mean no supplement is numbered
yet). Panel A is the comparison against known input. Panel B plots estimates against the 500x
estimate -- a convergence panel with 500x standing in for truth. An accuracy claim must cite
S1A specifically.

## 5. Reproducibility gap behind "(Methods)"

`METHODS.txt` carries a `[TO FILL]`: the inverse-chi-squared parameters (df, scale) and the
fitness -> expected-frequency mapping were never recorded, and neither the simulation script nor
the drawn fitness values were archived. Only the NNLS output survives; panel A's r2 values are
read back out of text embedded in the original per-trait PDFs by
`scripts/extract_sim_reported_r2.py`. "Compared the estimates to the known input values"
describes the original analysis correctly but cannot currently be re-run. Either recover the
script plus the 327 x 7 expected input frequencies, or state in Methods that these r2 are
reported from the original analysis.

Related, if the archive is described as clean NNLS output: 139 of 18,312 archived coefficients
(0.8%) are negative, which a strict non-negative solver cannot produce.

## 6. What the number checker does not catch

`python3 scripts/check_manuscript_numbers.py` passes this paragraph (exit 0, every number
attested): bare `1` and small integers are ignored and `500` is attested. The guard catches
drift in quoted values, not overstated claims.

## Suggested revision

To determine whether non-negative least squares (NNLS) regression can infer strain frequencies
from pooled populations, we established a simulation framework (Methods). Briefly, each wild
isolate was assigned a fitness value drawn from an inverse-chi-squared distribution, with the
fitness structure taken from seven published C. elegans traits with validated QTL, giving seven
independent simulated populations of 327 strains. For each, we computed the expected pooled
allele frequencies such a population would produce and simulated observed alt-allele counts by
binomial sampling at 1, 3, 5, 10, 30, 50, 100 and 500x. We then deconvolved these counts back to
per-strain frequencies by NNLS and compared the estimates with the known input as a function of
depth. Recovery was near-perfect at high depth -- every trait reached r2 = 1.00 by 500x -- and
degraded gracefully as depth fell: 10x was sufficient for five of the seven traits (r2 >= 0.95),
30x was the lowest depth at which all seven met that bar, and even 1x recovered most of the
signal for most traits (median r2 0.79, range 0.52-0.91) (Figure S1A).

## 7. Per-trait strain counts: the GWAS trait file is not 327 for six of seven traits

`simulation_gwas_traits.tsv.gz` (327 rows x 56 trait-depth columns) carries NAs. The non-NA
strain count per trait, constant across all eight depths within a trait:

| trait | strains | NA | % NA | r2 @1x | r2 @500x |
|---|---|---|---|---|---|
| mtDNA_ratio | 327 | 0 | 0.0 | 0.52 | 1.00 |
| Albendazole_q75.TOF | 199 | 128 | 39.1 | 0.56 | 1.00 |
| assay_norm | 152 | 175 | 53.5 | 0.81 | 1.00 |
| etoposide_median.TOF | 137 | 190 | 58.1 | 0.72 | 1.00 |
| amsacrine_f.L1 | 135 | 192 | 58.7 | 0.79 | 1.00 |
| value | 131 | 196 | 59.9 | 0.86 | 1.00 |
| PC1 | 84 | 243 | 74.3 | 0.91 | 1.00 |

Total non-NA strain-trait pairs: 1165 of the 2,289 (7 x 327) full grid.
Only 31 strains carry a value for all seven traits. mtDNA_ratio is the sole complete trait.

WHAT THE TWO FILES DISAGREE ABOUT. `simulation_nnls_frequencies.tsv.gz` is a complete
7 x 8 x 327 grid (18,312 rows, zero NAs) -- the NNLS itself ran on all 327 strains, which is why
`SUPPLEMENTAL_DATA_OVERVIEW.md` describes it that way. The NAs appear only in the GWAS trait
file built downstream from those estimates (`processed_simFreq_traits.tsv`, read by
`make_experiments_deposit.R:121`). So "seven independent simulated populations of 327 strains"
is accurate for the deconvolution and wrong for anything describing the simulated GWAS.

The NA set is NOT a magnitude filter: within a trait it is identical at every depth, and absent
strains routinely carry larger coefficients than present ones (e.g. assay_norm at 500x, max
absent coefficient 0.231 vs min present 0.000). It behaves like a fixed per-trait strain set --
most likely the isotype set of the original published trait -- but the upstream
`processed_simFreq_traits.tsv` is in the unarchived simulation directory, so the filter cannot
be confirmed from this clone.

CONSEQUENCE FOR PANEL B. See section 8 -- the pad does inflate the pooled r2, but not by the
zero-against-zero mechanism first suspected here.

## 8. Do the zeros correspond to the NA strains?

Not symmetrically. Tested at the strain-trait pair level (2,289 pairs) and the row level
(18,312 rows) of `simulation_nnls_frequencies.tsv.gz`.

ONE DIRECTION HOLDS EXACTLY. Every pair that is zero at all eight depths is an NA strain:
136 of 136. No strain the GWAS trait file keeps is zero at all eight depths (0 of 1,165).
So "zero everywhere" implies "excluded".

THE CONVERSE FAILS. Only 136 of the 1,124 NA pairs (12.1%) are zero at every depth; the other
988 carry a non-zero NNLS estimate at some depth, up to 0.2308 on the coefficient scale. Per
trait, NA-but-non-zero-somewhere: Albendazole_q75.TOF 101, PC1 206, assay_norm 162, value 170,
amsacrine_f.L1 175, etoposide_median.TOF 174.

ROW LEVEL. 95.2% of zero rows are NA strains (6,356 of 6,677); 70.7% of NA rows are zero
(6,356 of 8,992). The 321 zero rows that fall on kept strains concentrate at low depth --
128 at 1x, 68 at 3x, falling to 4 at 500x -- which is what binomial sampling zeroing out a
genuinely present strain looks like.

READING. The zeros are not the mechanism behind the NA set. NNLS assigned most of the excluded
strains non-zero frequency, i.e. the deconvolution ran against the full 327-strain reference
whatever the trait, and the GWAS trait file subset afterwards.

### Panel B, recomputed

Pooled r2 of the frequency estimate against the 500x estimate, all seven traits pooled:

| depth | all pairs (n=2,289) | non-NA pairs only (n=1,165) | delta |
|---|---|---|---|
| 1x | 0.787 | 0.699 | +0.088 |
| 3x | 0.913 | 0.865 | +0.048 |
| 5x | 0.941 | 0.903 | +0.038 |
| 10x | 0.967 | 0.944 | +0.023 |
| 30x | 0.987 | 0.978 | +0.009 |
| 50x | 0.993 | 0.987 | +0.006 |
| 100x | 0.995 | 0.991 | +0.004 |

The all-pairs column reproduces the caption (0.79 at 1x, 0.91 at 3x, 0.94 at 5x, 0.97 at 10x,
0.99 at 30x and 50x, 1.00 at 100x), so the panel B pipeline is confirmed. Restricting to the
1,165 scanned pairs costs 0.088 at 1x and 0.023 at 10x, converging above 30x. The inflation is
real but modest, and it comes from the excluded strains sitting at low frequency and widening
the range rather than from zero-against-zero pairs.
