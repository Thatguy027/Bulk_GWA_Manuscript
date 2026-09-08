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

## Reproducing the scan, and why every flag matters

A threshold is only meaningful for the scan it is computed against, so this
pipeline reproduces the 2023 pos-1 scan rather than doing the same thing in
spirit. Four attempts were needed. Each returned a different observed
genome-wide maximum against the scan's **8.8361**, and each difference had a
single cause:

| Run | observed_max | cause |
|---|---|---|
| smoke  | 8.6894 | MAF computed on the 231 phenotyped strains, so 519,341 markers were tested |
| smoke2 | 8.5700 | kinship built with `-gk 1` (centered); the scan used `-gk 2` (standardized) |
| smoke3 | 8.8759 | markers, panel and kinship type all correct — and still 34,316 markers short |
| current | — | the missing-genotype encoding below |

### The panel is 366 strains, not 231

The scan's plink step kept every strain in the phenotype file, so `--maf 0.05`
and `--geno` were computed across all 366. A marker carried by two of the 231
phenotyped strains can therefore be in the scan. GEMMA then drops the strains
with no value for the trait and analyses 231. Filtering on the 231 instead is
the first bug in the table.

### The kinship is `-gk 2`

`-gk 1` is the centered relatedness matrix, `-gk 2` the standardized one, where
each marker is divided by its own standard deviation before the cross-product.
Different matrices, different p-values. The scan's archived
`gemmeGRM.*.sXX.txt` names which.

### Missing genotypes become dosage 0, deliberately

This is the subtle one, and it is why the pipeline goes through plink's oxford
format rather than `.traw`.

The scan's plink log removed **339,834** variants at >10% missingness and kept
464,209 that each carry up to 10%. GEMMA then dropped **none** of them — but its
`-miss` default is 0.05, which those markers plainly exceed. A marker above 5%
missing exceeds that threshold under any individual set, so GEMMA cannot have
seen a missing value at all.

The reason is the encoding. `--recode oxford` writes a missing call as `0 0 0`,
and the dosage expression `2*P(AA) + P(AB)` turns that into **0** — a homozygous
call for the second allele, indistinguishable from a real one. A `.traw` route
passes `NA` through instead, GEMMA counts it missing, and 34,316 markers fall
out. That was the whole of the smoke3 discrepancy.

**This is not the better choice on its own terms.** It substitutes a fabricated
genotype for an absent one, and does so most often on the chromosome arms where
the hyper-divergent regions are and where calls fail: 41.8% of the markers in
the V:17 Mb window, 24% at V:18, 20% on the II left arm. It is here because the
published scan was built that way. `scripts/compare_observed_scan.R` in the main
repo quantifies what it costs: 94.75% of shared markers differ, the largest by
3.58 in -log10 p, and the biggest disagreements are a chrX 10.2-10.9 Mb block
that moves in both directions. The top of the scan is stable — the same peak
marker, `IV:15323414` — but a scan with missingness handled honestly is a
different scan and needs its own threshold. Emitting `NA` from that one awk
expression in `BUILD_BIMBAM` is the change.

### The marker list is shipped, not derived

`snps/pos1_2023_gemma_snps.tsv.gz` is the 464,208 ids GEMMA actually received,
from the scan's own `traits_gemmaSnps.tsv`. It is one marker short of plink's
464,209 because the scan built its annotation with `awk 'NR!=1'` on a `.gen`
file, which has no header, silently dropping the first variant. Shipping the
list reproduces the scan exactly without reproducing the bug in code:

```
464,209   plink retained (snps/pos1_2023_plink_traits.log)
     -1   the `awk 'NR!=1'` on a headerless .gen
464,208   the -snps list GEMMA received
   -163   MtDNA, which -loco never maps
464,045   shipped scan rows  = supplemental_data/mapping/pos1_2023_gemma_loco.csv.gz
```

`snps/pos1_2023_traits.sample` is the scan's oxford `.sample`. BIMBAM has no
sample ids — dosage columns are positional — so that file *is* the definition of
which column is which strain, and it is checked against plink's own output
rather than trusted.

## Three assertions, at three points

* `PLINK_CONVERT` fails unless it retains exactly `expect_variants` (464,209)
  and `expect_individuals` (366).
* `PREP_ORDER` fails if a strain in the scan's order is absent from the
  phenotype file, or if several requested traits have different phenotyped-strain
  sets — GEMMA's `-gk` takes its individuals from phenotype column 1, so one
  kinship cannot serve two different sets.
* `COLLECT_THRESHOLD` fails unless the observed genome-wide maximum equals
  `expect_observed_max` (8.8361) within `observed_tol`.

All three are needed, and the table above is why. The marker count passed while
the kinship was wrong; the marker count and the kinship together passed while
the encoding was wrong. Only the observed maximum tests the whole chain.

To threshold a trait with no shipped scan to compare against, pass
`--expect_observed_max 0`.

## Ruled out as causes of the smoke3 discrepancy

Each by measurement, not argument, and recorded so they are not re-litigated:

* **GEMMA version** — 0.98.5 on the cluster, 0.98.5 in the scan's archived log.
* **Phenotype file** — byte-identical to the `association_traits.csv` that
  produced the scan (md5 `b35a18aeeb0acfbc741c33c4159d12ae`), full precision,
  and the 231 mapped values differ by 0.
* **The strain set** — `panel.txt` was the same 231 as the trait's non-missing
  set and as trait 1's, so the kinship's individuals matched however GEMMA
  selects them.
* **MAF** — of the 34,316 markers GEMMA dropped, *none* has af below 0.01 or
  above 0.99; their frequencies are unremarkable (median 0.139). Not the MAF
  filter, which left `-miss` as the only remaining default that discards a
  marker.
* **The LOCO split** — GEMMA's own logs reported `total SNPs` matching the
  expected per-chromosome kinship sizes exactly (I 412322, II 388154,
  III 399622, IV 387572, V 341788, X 390767).
