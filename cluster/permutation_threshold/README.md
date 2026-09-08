# Permutation threshold for the pooled GWA scan

Genome-wide significance measured in this panel rather than approximated. Drops
in beside `gemma_nf` on Hoffman2 and follows its conventions: SGE executor,
`-V -l h_data=...,h_rt=...`, the same plink/gemma binaries, and the same
R-environment fix-up.

## Transfer and run

```sh
scp -r permutation_threshold thatguy0@hoffman2.idre.ucla.edu:/u/project/kruglyak/thatguy0/nextflow_pipes/
ssh thatguy0@hoffman2.idre.ucla.edu
cd /u/project/kruglyak/thatguy0/nextflow_pipes/permutation_threshold

# smoke test first -- 20 permutations, ~15 minutes, proves the plumbing
nextflow run main.nf -profile hoffman2 \
  --pheno traits/2023_pos1_association_traits.csv \
  --traits vst_ctrl_pos-1_T2 \
  --n_perm 20 --perm_batch 10 --name smoke

# then the real thing
nextflow run main.nf -profile hoffman2 \
  --pheno traits/2023_pos1_association_traits.csv \
  --traits vst_ctrl_pos-1_T2 \
  --n_perm 1000 --perm_batch 25 --name pos1_perm1000
```

`-resume` works: the conversions and the six kinship matrices are cached, so a
re-run after a queue eviction only redoes the mapping jobs that failed.

## What comes out

| file | contents |
|---|---|
| `permutation_thresholds.tsv` | the threshold per trait per alpha, the observed maximum, and its empirical p |
| `permutation_maxima.tsv` | genome-wide maximum for every permutation — the null distribution itself |
| `permutation_threshold.{pdf,png}` | that null, with the threshold, the observed scan, and the eigen/Bonferroni lines for comparison |
| `kinship/` | the six LOCO kinship matrices, reusable by other runs |

Read the threshold off `permutation_thresholds.tsv` at `alpha == 0.05`. For the
manuscript it replaces the eigen value of 4.60 in
`scripts/gwas_qtl_intervals.R` — set `THRESHOLD` there and nothing else changes.

## How it works, and the one thing that makes it affordable

Genotypes are fixed and the phenotype is shuffled across strains. Each shuffle
is mapped genome-wide with the same leave-one-chromosome-out kinship matrices as
the real scan, and only the genome-wide **maximum** −log10 p is kept. The
(1 − alpha) quantile of those maxima is the threshold.

Shuffling the phenotype does not change the genotypes, so **the kinship matrices
are computed once and reused by every permutation**. Recomputing them per
permutation would multiply the cost by `n_perm` for an identical answer.

The observed phenotype rides along as permutation 0, through the same code path,
so the observed maximum and the permutations cannot diverge for reasons
unrelated to the null.

## Cost

Jobs are `(n_perm / perm_batch) × 6` mapping jobs, plus 6 kinship matrices and
the conversions. At `n_perm 1000, perm_batch 25` that is **240 mapping jobs**,
each running 25 GEMMA calls in series. Batching exists because 6,000
single-call jobs would spend most of the wall clock in the SGE queue. With
`queueSize 50` and roughly a minute per GEMMA call on this panel, expect
**2–3 hours**.

Raise `--perm_batch` to cut scheduler pressure, but raise `h_rt` for
`GEMMA_PERM` in `nextflow.config` with it — that job's wall clock scales with
the batch.

## What permuting a label does and does not measure

**Read this before quoting the number.** Shuffling phenotype labels destroys the
relatedness structure the mixed model is fitted to. The null it samples is
therefore *no association and no population structure*, while the model assumes
structure is present and corrects for it. This is the standard permutation
threshold in the *C. elegans* GWAS literature and is what cegwas/NemaScan
report, but it is not exact — it tends to be slightly **anti-conservative**
where structure inflates the real scan, because the permuted scans have no
structure to inflate.

The defensible claim is "a threshold calibrated to this panel's linkage
disequilibrium and marker density", not "an exact family-wise error rate". The
structure-preserving alternative is to permute in the space rotated by the
eigenvectors of the kinship matrix. That is not implemented here, and it is the
thing to do if a reviewer presses on it.

## Ordering, which is the failure mode to watch

BIMBAM files carry no sample IDs — dosage columns are positional, in `.fam`
order. A phenotype vector in a different order silently maps the wrong values to
the wrong strains and yields a plausible, meaningless scan. `bin/make_permutations.R`
joins by strain ID against the `.fam` and refuses to continue if nothing matches
or if fewer than 20 individuals carry a value. If the trait file ever stops using
the same strain names as the VCF, that assertion is what will catch it.

Strains with a missing phenotype stay missing and are never shuffled into
phenotyped positions — otherwise the sample size would change between
permutations and the maxima would not be comparable.

## The smoke run found a real problem — read this before the long run

The smoke run completed but returned **`observed_max = 8.6894`** where the
shipped scan's genome-wide maximum is **8.84**. That mismatch is the whole point
of carrying the observed phenotype as permutation 0, and it meant the pipeline
was not mapping the same data as the scan it was meant to threshold.

The plink log said why: **519,341 markers passed filters** against the scan's
**464,045**. Two causes, both now fixed.

**MAF was computed on all 540 isotypes, not the 231 phenotyped ones.** A marker
at 5% across the whole collection can sit below 5% among the strains that
actually carry a *pos-1* value. `PREP_PANEL` now writes the phenotyped panel and
`--keep` restricts the conversion to it. On chromosome III alone that changes
the count from 80,639 to 78,507.

**The scan's filter chain is not recoverable, so it is no longer guessed.**
Instead the permutation scan tests **exactly the markers the scan tested**,
supplied as an id list (`markers/pos1_2023_scan_markers.txt.gz`, 464,045 ids,
all of which resolve in the CeNDR plink set) and applied with `--extract`. No
`--maf` or `--geno`: the list defines the set, and any further filter would
silently shrink it. `PLINK_CONVERT` now **fails the run** if the retained count
is not `expect_markers` (464045), rather than discovering the problem after a
thousand permutations.

**One trait per panel, enforced.** Traits in this file do not share a panel —
the three *pos-1* traits have 231 strains each, `negctrl_growth_HT115_delta_t0`
has all 366. Requesting traits with different panels together would compute MAF
on a superset for at least one of them, reintroducing the first bug by the back
door, so `prep_panel.R` refuses and tells you to split the run.

After re-running, **`observed_max` must be 8.84**. If it is not, stop and work
out why before trusting the threshold.

## Fixed after the first cluster attempt

Three things, recorded because two of them would have produced a *wrong number*
rather than an error.

**Params must be declared in `nextflow.config`, not `main.nf`.** The config is
parsed before the script, so a `${params.x}` interpolated inside a `process` or
`profiles` block can only resolve against params defined in the config file.
Declaring them in the script gave

```
Unknown config attribute `process.withName:MAKE_PERMS|COLLECT_THRESHOLD.params.r_env_bin`
```

They now live in a `params { }` block at the top of the config, above everything
that interpolates them — config is evaluated top-down, and a reference above its
definition silently yields `[:]/...` instead of failing.

**GEMMA now runs `-lmm 1`, not `-lmm 4`, and `p_wald` is found by name.**
`-lmm 4` emits `p_wald`, `p_lrt` *and* `p_score`, so reading the last column
positionally picked up **p_score** — a different test statistic from the one the
shipped scan reports, which would have thresholded the wrong thing. The column
is now located from the header, so a GEMMA version that reorders its output
cannot silently change the answer.

**BIMBAM allele order.** `.traw` dosages count the `COUNTED` allele (column 5),
and BIMBAM's dosages count the allele listed *first*. The first version emitted
column 6 then 5, flipping the coding. That changes the sign of `beta` and leaves
`p_wald` alone, so it would not have broken this threshold — but it would have
quietly corrupted any effect size read from those files.

## Notes on this copy

`traits/2023_pos1_association_traits.csv` is the same file the shipped scan was
run on: 366 strains, 231 with a `vst_ctrl_pos-1_T2` value, the trait behind the
464,045-marker scan in `supplemental_data/mapping/pos1_2023_gemma_loco.csv.gz`.
Replace it or point `--pheno` elsewhere for another experiment.

Paths to `plink`, `gemma`, the conda bin and the R env are the defaults from
`gemma_nf`'s config. Override on the command line if any of them moves.

## Two assertions, at both ends

`PLINK_CONVERT` fails unless it retains exactly `expect_markers` (464,045)
markers, and `COLLECT_THRESHOLD` fails unless the observed genome-wide maximum
equals `expect_observed_max` (8.8361) within `observed_tol`.

Both are needed, and the history says why. The marker assertion was added after
a first run computed MAF on all 540 strains instead of the 231 phenotyped ones
and tested 519,341 markers, giving 8.6894. With the marker set pinned, the next
run matched 464,045 markers and 231 strains exactly -- and still returned 8.5700,
because the kinship was built with `-gk 1` (centered) where the scan used
`-gk 2` (standardized). A correct panel and a correct marker set are not
sufficient; the model has to match too, and only the observed maximum tests that.

To threshold a trait with no shipped scan to compare against, pass
`--expect_observed_max 0`.

## The observed maximum: what has been ruled out

Target 8.8361, the maximum of `supplemental_data/mapping/pos1_2023_gemma_loco.csv.gz`.

| Run | observed_max | cause |
|---|---|---|
| smoke  | 8.6894 | MAF computed on all 540 strains, 519,341 markers tested |
| smoke2 | 8.5700 | kinship built with `-gk 1` (centered) where the scan used `-gk 2` |
| smoke3 | 8.8759 | open -- 0.4%, cause not yet identified |

Eliminated as causes of the remaining 0.4%, each by measurement rather than
by argument:

* **Marker set.** `--extract` of the scan's own 464,045 IDs; PLINK_CONVERT
  asserts the count and reported exactly 464045.
* **Panel.** 231 strains, and `panel.txt` is the same set as the trait's
  non-missing strains and as trait 1's, so the GRM's individuals match however
  GEMMA selects them.
* **Phenotype.** `traits/2023_pos1_association_traits.csv` is byte-identical
  (md5 b35a18aeeb0acfbc741c33c4159d12ae) to the `association_traits.csv` that
  produced the scan, at full precision; the 231 mapped values differ by 0.
* **GEMMA version.** 0.98.5 on the cluster, 0.98.5 in the scan's archived log.
* **Missing genotypes.** `n_miss` is 0 for all 464,045 markers in the scan, so
  the oxford (`dosage 0`) versus `.traw` (`NA`, mean-imputed) difference in
  missing handling has nothing to act on.
* **p-value column.** `-lmm 1`, and `p_wald` located by header name.

What remains, in order of suspicion:

1. **The kinship's marker set.** gemma_nf hands GEMMA the whole genotype file
   with `-loco ${chrom}`; this pipeline pre-splits with plink `--not-chr` and
   passes no `-loco`. GEMMA_GRM now publishes its log, so its analysed count can
   be compared against the expected kinship size for each chromosome:
   I 412322, II 388154, III 399622, IV 387572, V 341788, X 390767.
2. **Individuals in the genotype file.** gemma_nf passes 366 columns and lets
   GEMMA select the 231 by phenotype missingness; this passes 231 columns.
   Equivalent unless GEMMA's internal MAF filter uses the file rather than the
   analysed subset -- the scan retains all 464,045 markers, so it dropped none,
   and a difference here would show as a marker-count difference.
3. **`-loco` at the mapping step**, which gemma_nf passes and this does not,
   having already split the genotypes by chromosome.

`scripts/compare_observed_scan.R` in the main repo distinguishes 1 from 3: a
near-constant ratio across the whole range of the statistic is the signature of
a variance-component difference, disagreement confined to particular markers is
not.
