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

## Notes on this copy

`traits/2023_pos1_association_traits.csv` is the same file the shipped scan was
run on: 366 strains, 231 with a `vst_ctrl_pos-1_T2` value, the trait behind the
464,045-marker scan in `supplemental_data/mapping/pos1_2023_gemma_loco.csv.gz`.
Replace it or point `--pheno` elsewhere for another experiment.

Paths to `plink`, `gemma`, the conda bin and the R env are the defaults from
`gemma_nf`'s config. Override on the command line if any of them moves.
