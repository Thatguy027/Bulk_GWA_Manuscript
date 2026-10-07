# 2023 *pos-1*: deconvolution against the pool, not the whole panel

Kept separate from `supplemental_data/` on purpose. Nothing here replaces the
deposited analysis; it exists so the two can be compared before anything does.

## The problem

The deposited analysis deconvolved every sample against a **367-column**
genotype reference covering the whole strain panel. The 2023 *pos-1* pool was
not the whole panel. It is sets **B, C, E, F and G** of
`meta_files/BulkCe_strainsets.tsv` — 226 strains, **224** of them with a CeNDR
isotype.

Offering NNLS 142 strains that were never in the tube lets it place mass on
them, and that mass comes off the strains that were. Measured against the
deposited traits:

| | |
|---|---|
| strains in the deposited output | 366 |
| …actually in the pool | 224 |
| …**never** in the pool | **142** |
| strains given a phenotype | 231 |
| …in the pool | 166 |
| …**not** in the pool, so leakage | **65** |
| pool strains left **without** a phenotype | **58** |

A second, smaller defect falls out of the same place: the reference carries
`JU1793_JU1793` **twice**. Two identical columns, so NNLS splits that strain's
mass between them arbitrarily — which is the duplicated-row behaviour
`METHODS.txt` documents downstream. One column per isotype removes it at source.

## What was done

Same markers, same counts, same solver. The marker filter is not re-derived: it
is inherited by reusing the `gt` matrix the original run stored for each depth
cutoff, so the only thing that changes is which columns the design matrix has.

Before trusting the corrected fit, `01` reproduces the original 367-column fit
from those same inputs. It agrees at *r* = 0.997 but not exactly, and that is a
finding rather than a defect: with 367 near-collinear columns the NNLS solution
is **not uniquely determined**. Refitting lands on a different near-optimal
member of the same set, with a slightly *lower* residual. The 367-column fit put
mass on 305 of 367 strains; the pool fit puts it on 211 of 224.

`02` builds the traits, using the upstream pipeline's own transforms
(`scripts/07_make_phenotypes.R` and `04_control_corrections.R` in
`/Users/Stefan/UCLA/Projects/bulkGWAS/lipid_RNAi`, documented in
`TRAIT_SPEC.md`):

```
delta_ctrl  = f - p
vst_ctrl    = asin(sqrt(f)) - asin(sqrt(p))
log2fc_ctrl = log2((f + c) / (p + c)),  c = half the smallest non-zero frequency
```

with `f` the mean frequency over the four *pos-1* replicates, `p` the control
frequency, and all three NA where `p == 0` — the pipeline's own `usable` rule.
Fed the *deposited* frequencies, this reproduces the deposited
`delta_ctrl_pos-1_T2` to **2.5e-16** and `vst_ctrl_pos-1_T2` to **7.8e-16**, so
the formulae are the originals and any difference in the new traits comes from
the deconvolution. JU1793 is the single exception, for the reason above.

## What came out

| | deposited | pool reference |
|---|---|---|
| strains with a phenotype | 231 | **184** of 224 |
| of those, not in the pool | 65 | **0** |
| declining under *pos-1* | 183 (79%) | 140 (**76%**) |
| Spearman, old vs new, shared strains | — | **0.968** (n = 164) |

The phenotype **values** are stable for strains that were genuinely in the pool
(Pearson 0.997). What was wrong was *which* strains had one.

Replicate reproducibility, on the corrected frequencies and the 184 strains with
a phenotype:

```
rep1 vs rep2  0.901     rep1 vs rep3  0.875
rep3 vs rep4  0.887     rep1 vs rep4  0.837
rep2 vs rep3  0.884     rep2 vs rep4  0.829
```

## The association scan

Run on `vst_ctrl_pos1_T2`, GEMMA 0.98.5 LOCO, kinship per chromosome. GEMMA's
own logs record the panel: **224 individuals supplied, 184 analysed** — it drops
the 40 with no phenotype, so the scan is on the strains that have one.

Thresholds are recomputed for this panel rather than carried over, because the
effective number of independent tests falls with the number of strains:
Bonferroni **6.96** over 457,571 markers, and the Li & Ji eigenvalue threshold
**4.50** from M_eff **1590** (against 4.60 from M_eff 1972 on the deposited
231-strain panel). Method as in `scripts/eigen_independent_tests.R`.

| chr | deposited peak | corrected peak | position |
|---|---|---|---|
| I | 3.59 | 4.00 | — |
| II | 3.99 | 3.20 | — |
| **III** | **8.68** | **5.90** | III:5,965,738 (same marker) |
| **IV** | **8.84** | **8.22** | IV:15,323,414 (same marker) |
| V | 5.40 | 5.78 | — |
| **X** | **7.83** | **8.05** | X:4,875,969 (same marker) |

**Chromosome IV and chromosome X survive.** Both peak on the same marker as
before, both clear Bonferroni, and the X locus is *stronger* on the corrected
phenotype than on the contaminated one.

**Chromosome III does not.** The same marker falls from 8.68 to 5.90 — below
Bonferroni, though still above the eigenvalue line. It was the second strongest
association in the deposited scan.

Chromosome V now clears the eigenvalue threshold (5.78) where it did not before,
so it is the locus to look at next if the eigenvalue line is the one being used.

## Not done

- **Figure 1C** is not redrawn in the manuscript's own style; the comparison
  figure here is `plots/mapping_old_vs_pool_reference.png`.

`METHODS.txt` in the manuscript repository carries a `[TO FILL]` saying the vst
definition is not reproducible from the deposit. It is reproducible — just not
from the deposit. That note can be closed with the formula above.

## Files

```
scripts/01_deconvolve_pool_reference.R   the 224-column deconvolution
scripts/02_build_traits.R                traits, with the formula check
scripts/03_figures.R                     old-vs-new comparison figures
data/pool_reference_frequencies.csv.gz   224 strains x 7 samples x 3 cutoffs
data/pool_reference_traits_dp5.csv       all trait columns
data/mapping_traits_dp5.csv              strain + delta, vst and log2fc
plots/phenotype_old_vs_pool_reference.*  distribution old vs new, and a scatter
plots/replicate_reproducibility_*.*      all six replicate pairs
scripts/04_mapping_figure.R              thresholds and the two scans
mapping/vst_...loco_results.csv.gz       GEMMA output, 457,571 markers
mapping/gemma_LOCOmapping.*.log.txt      one GEMMA log per chromosome
mapping/plink/traits.{log,sample}        how the genotypes were prepared
data/eigen_independent_tests_pool.tsv    M_eff per chromosome, this panel
plots/mapping_old_vs_pool_reference.png  both scans, with both thresholds
```

`mapping/plink/traits.gen` (627 MB) is not tracked; `traits.log` records the
exact plink invocation that rebuilds it from the 20210121 VCF.

Run in order from the repository root. `01` reads the upstream analysis
directory `/Users/Stefan/UCLA/Projects/bulkGWAS/lipid_RNAi/2023_pos1` and the
metadata at `.../2023_original_pos1/meta_files`, neither of which is in the
repository, so this does not run from a clone.

## If this is adopted

It invalidates more than the figures here: Figure 1B's distribution, Figure S4's
reproducibility, and the Results numbers built on the 231-strain table — 79%
declining, 81 absent from every *pos-1* pool, 141 of 231 responsive, and the
binomial *p*.

For Figure 1C and the text around it, the claim that changes is the number of
QTL. Two of the three survive on the same markers; the chromosome III locus does
not clear Bonferroni on the corrected phenotype.
