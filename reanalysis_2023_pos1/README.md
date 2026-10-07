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

`02` builds the traits. Fed the *deposited* frequencies, the same trait code
reproduces the deposited `delta_ctrl_pos-1_T2` to **2.5e-16** — so any
difference in the new traits comes from the deconvolution, not the formula.
JU1793 is the single exception, for the reason above.

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

## Not done

- **`vst` is not rebuilt.** The deposited mapping trait was
  `vst_ctrl_pos-1_T2`. It is not a function of the delta alone — against the
  deposited delta it is Spearman 0.870, not 1 — so it cannot be recovered from
  what is here, and its definition is in neither this repository (`METHODS.txt`
  carries a `[TO FILL]` saying exactly that) nor the upstream directory, which
  holds no code. The traits written here are the two that are reproducible.
- **The association scan.** Run it on `data/mapping_traits_dp5.csv`.
- **Figure 1C** cannot be redone until that scan exists.

## Files

```
scripts/01_deconvolve_pool_reference.R   the 224-column deconvolution
scripts/02_build_traits.R                traits, with the formula check
scripts/03_figures.R                     old-vs-new comparison figures
data/pool_reference_frequencies.csv.gz   224 strains x 7 samples x 3 cutoffs
data/pool_reference_traits_dp5.csv       all trait columns
data/mapping_traits_dp5.csv              strain + the two mapping traits
plots/phenotype_old_vs_pool_reference.*  distribution old vs new, and a scatter
plots/replicate_reproducibility_*.*      all six replicate pairs
```

Run in order from the repository root. `01` reads the upstream analysis
directory `/Users/Stefan/UCLA/Projects/bulkGWAS/lipid_RNAi/2023_pos1` and the
metadata at `.../2023_original_pos1/meta_files`, neither of which is in the
repository, so this does not run from a clone.

## If this is adopted

It invalidates more than the figures here: Figure 1B's distribution, Figure 1C's
scan and its three QTL, Figure S4's reproducibility, and the Results numbers
built on the 231-strain table — 79% declining, 81 absent from every *pos-1*
pool, 141 of 231 responsive, and the binomial *p*.
