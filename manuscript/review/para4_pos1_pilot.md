# Review of manuscript paragraph 4 -- the pos-1 pilot and the plate re-evaluation

The paragraph is written against an EARLIER version of this experiment. Three of its numbers
come from a two-replicate dataset and a figure that no longer exists in the repository; the
current export has four pos-1 replicates and the figures were rebuilt on a different scale.
Everything below is recomputed from `supplemental_data/phenotypes/`.

## THE HEADLINE DISCREPANCY: the dataset was replaced under the text

`scripts/SUPP_FIG_XX_original_pos1_dfreq_rep_correlation.R`, in its own words:

  "The updated dataset has FOUR pos-1 replicates where the earlier version had two, so the
  replicate figure now shows all six pairwise comparisons rather than a single scatter."

`pos1_2023_sample_frequencies.csv.gz` confirms it: 366 strains, T2 only, **control replicates
A and B, pos-1 replicates 1-4**, each at read-depth cutoffs 3, 5 and 10.

So "we grew these populations for two generations in these conditions across two replicates"
and a single "Spearman's rho = 0.72" describe the superseded two-replicate scatter. The current
figure has six panels, each with its own rho and n printed on it.

### 1. rho = 0.72 matches no pair in the current data

`scripts/SUPP_FIG_XX_original_pos1_dfreq_rep_correlation.R` applies no strain filter: it
collapses duplicate (strain, sample) rows by summing, takes the control baseline as the mean
over both control replicates, and drops only NA pairs. Reproducing that recipe exactly gives
what the six facets print:

| pair | Spearman rho | as printed | n |
|---|---|---|---|
| rep1 vs rep2 | 0.867 | 0.87 | 366 |
| rep1 vs rep3 | 0.771 | 0.77 | 366 |
| rep1 vs rep4 | 0.781 | 0.78 | 366 |
| rep2 vs rep3 | 0.804 | 0.80 | 366 |
| rep2 vs rep4 | 0.781 | 0.78 | 366 |
| rep3 vs rep4 | 0.813 | 0.81 | 366 |

**The range to quote is 0.77 to 0.87, n = 366 on every facet.** Nothing is 0.72, and there is
no longer a single number to quote. (Excluding the 89 strains whose delta is exactly zero in
all four replicates raises the band to 0.81-0.90, but the figure does not exclude them, so
quoting that band would not match the panel.)

### 2. "146 of 224 (65%)" is one replicate of a figure that is not in the repository

The source is the caption for `SUPP_FIG_XX_original_pos1_strain_freq_change`
(`scripts/pos1_script1.R`), which reads: "224 strains: in replicate A, 146 decrease and 78
increase; in replicate B, 144 decrease and 80 increase."

Two problems. First, 146 is **replicate A alone** -- replicate B gives 144 -- so the sentence
presents one replicate's count as the experiment's result. Second, that figure is **not in
`plots/`, not in `plots/legacy/`, and not in the curated set of twenty**; it is a caption for
output of a legacy script on the two-replicate dataset. There is nothing to cite.

What the current data support, by definition of "responsive":

| definition | responsive | total | % |
|---|---|---|---|
| replicate A (superseded dataset, per caption) | 146 | 224 | 65.2 |
| replicate B (superseded dataset, per caption) | 144 | 224 | 64.3 |
| shipped `delta_ctrl_pos-1_T2` < 0 | 183 | 231 | 79.2 |
| negative mean delta across the four pos-1 replicates | 183 | 277 | 66.1 |
| negative in ALL four pos-1 replicates | 151 | 277 | 54.5 |

65% matches none of them on the current data. The strictest and most defensible statement is
the last row: 151 of 277 strains (55%) decrease in frequency under pos-1 RNAi in every one of
the four replicates. If the denominator should be the pool as constructed rather than the
reference, note that the trait file carries 366 strains of which 231 have a value -- 224 is not
reproducible from the deposit and needs to come from the lab record.

### 3. rho = 0.42, n = 106, p = 6e-6 is the superseded panel A

`SUPP_FIG_plate_vs_paaby_vs_pos1original` was rebuilt. From its caption: "The earlier version
of panel A used a bootstrap frequency change from a different export (rho = 0.42, n = 106),
which was on no scale used elsewhere in the manuscript; the VST version is on the same scale as
every other phenotype panel."

Recomputed from the deposit, panel A is now plate score against `vst_ctrl_pos-1_T2`:
**rho = 0.410, n = 111, p = 7.79e-06** -- reproducing the caption's 0.41 / 111 / 7.8e-06
exactly. Against the shipped `delta_ctrl_pos-1_T2` instead it is rho = 0.340, n = 111,
p = 2.6e-04, which is neither the old value nor the new one.

## 4. The Paaby comparison is correct, but it is panel B

rho = -0.55, n = 19, p = 0.014 reproduces exactly (recomputed as unhatched eggs over eggs plus
larvae per well from the machine-scored counts, averaged within strain, pos-1 clone only; the
manual count columns give -0.592, p = 0.008, so use the machine columns to match the figure).

Two corrections to the sentence: the text cites both comparisons to `SUPP FIGURE XX_A`, but the
Paaby comparison is **panel B**. And the caption's own instruction: with only 19 shared strains
it "should be described as consistent rather than as independent confirmation." The published
measurement is embryonic lethality for the pos-1 clone, which is worth naming rather than
calling it an "evaluation of wild isolate RNAi responses".

## 5. Broad-sense heritability 0.32 is not computed anywhere in this repository

No heritability estimate appears in `METHODS.txt`, `FIGURE_CAPTIONS.txt`, `FIGURE_REPORT.md` or
any script. The number checker passes `0.32` only because it rounds onto the leakage Spearman
rho (0.323/0.326) from an unrelated analysis in the report -- the same bag-of-numbers accident
as the earlier `2.85`. Either compute it and deposit it (the four pos-1 and two control
replicates support a repeatability/variance-component estimate) or cite its source.

## 6. Numbers that are right but unattested

`191` is exact: `plate_scores_pos1.tsv` holds 191 strains. `93` matches `METHODS.txt:64` ("Nine
of the 93 strains in the pool have no vst value and were dropped"). The checker fails both 191
and -0.55 anyway, because it attests only against `FIGURE_REPORT.md` and neither number is in
it. Add them to the report or to `manuscript_number_exceptions.txt`.

For reference, the checker's verdict on this paragraph: FAIL on 146, 224, 191 and -0.55, and it
PASSES 0.72, 0.32, 0.42 and 6e-6 -- i.e. it flags two correct numbers and clears four wrong
ones.

## The paragraph, reworded

Numbers here are figure-exact: the replicate range is what the six facets print, 0.41/111/7.8e-06
reproduces panel A, -0.55/19/0.014 reproduces panel B. Two slots still need you (marked).

> C. elegans RNAi screens are routinely performed in the N2 strain. However, extending these
> screens to multiple wild isolates is challenging because each isolate needs to get assayed
> independently. We reasoned that the pooled phenotyping approach would facilitate fast
> evaluation of wild isolate RNAi responses. N2 strains exposed to pos-1 RNAi exhibit an
> embryonic lethal phenotype, making this treatment amenable to the pooled phenotyping because
> RNAi-responsive strains will fail to contribute progeny to subsequent generations. To explore
> if we could quickly evaluate wild isolate RNAi responses, we performed a pilot experiment
> where we exposed 224 pooled wild isolates to pos-1 RNAi-expressing HT115 bacteria or control
> HT115 bacteria (Methods). We grew these populations for two generations in these conditions,
> sequenced the resulting F3 L1 populations, and inferred the individual strain frequencies in
> each population, across four pos-1 replicates and two control replicates. We calculated
> individual strain frequency differences between the control and pos-1 RNAi conditions and
> found good agreement between replicates (Spearman's rho = 0.77-0.87 across all six pairwise
> comparisons of the pos-1 replicates, n = 366; SUPP FIG XX). As expected from previous reports
> that have identified substantial RNAi-response variation across wild C. elegans isolates, we
> observed that 183 of the 231 strains with a pos-1 response value (79%) were responsive to
> RNAi, as indicated by these strains having a lower frequency in pos-1 RNAi conditions than
> they did in the control condition (Figure 1B). We observed substantial variation within
> RNAi-responsive and -insensitive strain groups and estimated the broad-sense heritability to
> be [HERITABILITY], suggesting that genetic factors influence strain responses to pos-1 RNAi.
> These results motivated us to construct a pooled population of RNAi-responsive strains that we
> could expose to different RNAi conditions. To construct this population, we manually
> re-evaluated pos-1 RNAi responses of 191 wild strains on agar plates and identified 93 strains
> with robust RNAi responses that we pooled to evaluate additional RNAi responses (Methods). We
> compared the manual pos-1 RNAi phenotypes we collected to the pooled pos-1 response on the
> same VST scale used for association mapping and found good agreement between the methods
> (Spearman's rho = 0.41, p = 7.8x10-6; n = 111, SUPP FIGURE XX_A). Our manual phenotypes are
> also consistent with a previously published measurement of pos-1 embryonic lethality in wild
> isolates (Spearman's rho = -0.55, p = 0.014; n = 19; SUPP FIGURE XX_B) (Paaby et al. 2015).

### What changed and why

| was | now | reason |
|---|---|---|
| "across two replicates" | four pos-1 and two control replicates | the export has replicates 1-4 and A/B |
| rho = 0.72, p < 0.001 | rho = 0.77-0.87, n = 366, six comparisons | 0.72 was the old single scatter |
| 146 of 224 (65%), SUPP FIG | 183 of 231 (79%), Figure 1B | the 224/146 figure is not in the repository; Figure 1B is the curated panel that shows this distribution |
| rho = 0.42, p = 6e-6, n = 106 | rho = 0.41, p = 7.8e-06, n = 111 | panel A was rebuilt on the VST scale |
| Paaby cited to panel A | panel B | it is panel B of that supplement |
| "a previously published evaluation of wild isolate RNAi responses" | "a previously published measurement of pos-1 embryonic lethality" | that is what the Paaby data are |
| "good agreement" for Paaby | "consistent with" | the caption asks for this given 19 shared strains |

### The two slots

1. **[HERITABILITY]** -- 0.32 is not computed anywhere in the repository. I can estimate
   repeatability across the four pos-1 replicates if you want a number derived from these data,
   but it would be a new estimate, not a recovered one.
2. **224** -- kept as you wrote it, but it is not reproducible from the deposit, and it sits
   awkwardly with the sentence that follows: the trait file gives 231 strains with a pos-1
   response value, seven MORE than the stated pool size. If 224 is the pool as constructed, then
   231 includes strains the deconvolution assigned frequency to without their being pooled --
   exactly the leakage the dilution experiment quantifies -- and the responsive fraction should
   be computed on pool members only. Send me the pool composition list and I will recompute it
   restricted to those strains.

## 7. Replacing the heritability sentence with the QTL

Recomputed from `supplemental_data/mapping/pos1_2023_gemma_loco.csv.gz` (464,045 markers,
n = 231), thresholds as given in `FIGURE_CAPTIONS.txt`: Bonferroni -log10 p = 6.97, eigen 4.60
over 1,972 effective independent tests. 10 markers clear Bonferroni, 465 clear eigen.

### The ten Bonferroni markers

`support` is the number of markers in the same 100 kb window that also clear the eigen
threshold, from `plots/diagnostics/TABLE_gwas_local_support.tsv`; `window` is how many markers
that region contains at all.

| marker | Mb | -log10 p | allele freq | beta | support | window |
|---|---|---|---|---|---|---|
| IV:15323414 | 15.323 | 8.84 | 0.208 | +0.0258 | 95 | 1094 |
| **III:5965738** | **5.966** | **8.68** | **0.091** | **+0.0608** | **0** | **628** |
| IV:15323794 | 15.324 | 8.52 | 0.208 | +0.0252 | 95 | 1093 |
| X:4875969 | 4.876 | 7.83 | 0.087 | +0.0363 | 2 | 926 |
| IV:13408895 | 13.409 | 7.49 | 0.100 | +0.0314 | 20 | 703 |
| IV:13413079 | 13.413 | 7.49 | 0.100 | +0.0314 | 20 | 697 |
| IV:13404290 | 13.404 | 7.42 | 0.100 | +0.0311 | 20 | 701 |
| IV:15986656 | 15.987 | 7.29 | 0.450 | +0.0184 | 17 | 1093 |
| IV:13405229 | 13.405 | 7.11 | 0.095 | +0.0311 | 20 | 694 |
| IV:13405230 | 13.405 | 7.11 | 0.095 | +0.0311 | 20 | 694 |

Every beta is positive: the minor allele raises the VST response, i.e. toward resistance. All
but IV:15986656 (af 0.450) are rare in the panel, af 0.087-0.208.

**III:5965738 is the only one of the ten with no local support at all**, in a window holding
628 markers. Genome-wide only 14 of the 465 eigen-significant markers have support 0, and this
is the only Bonferroni marker among them. See the correction below.

### CORRECTION to the framing, and then a correction to the correction

The first pass through this file said: "Chromosome III is not eigen-only -- III:5965738 clears
Bonferroni at 8.68, the second strongest marker genome-wide." That is arithmetically true and
substantively wrong, and it should not go into the draft in that form.

**III:5,965,738 does not behave like an association.** Two independent lines in this
repository agree, having been written without reference to each other:

- `scripts/chrIII_association_support.R` (commit 97057ad) tests every Bonferroni marker for
  local support, on the reasoning that a causal variant is detected through the markers in
  linkage disequilibrium with it, so a real association is a CLUSTER. Of the 628 markers within
  100 kb of this one, **zero** clear the eigen threshold and their median -log10 p is 0.11.
- `TABLE_gwas_local_support.tsv`, from the interval pipeline, records the same marker at
  `support = 0, n_window = 628` -- the identical numbers, and the only one of the ten
  Bonferroni markers with no support. Every other has 2 to 95 supporting markers.

The recombination domain is what makes this decisive rather than suggestive. 5.97 Mb sits in
the chromosome III low-recombination CENTRE, where LD blocks are long and support is easiest to
come by; a real association there should drag up markers across hundreds of kilobases. It drags
up none. A lone significant marker in a well-covered region is the signature of something that
tracks no haplotype -- genotyping error, an alignment artefact, a mismapped duplication.

The interval pipeline reaches the same verdict structurally: **III:5,965,738 appears in neither
`TABLE_gwas_qtl_intervals.tsv` nor its eigen variant.** The three chromosome III loci that do
are 3.818, 4.407 and 12.718 Mb. It was not dropped as a bad marker upstream -- it is absent from
`TABLE_scan_markers_dropped.txt` -- it was excluded for having no supporting markers to build
an interval from.

**What is left on chromosome III is the right-arm cluster**, 15 markers over 12.702-12.800 Mb
peaking at III:12,718,465 (6.31), which clears the eigen threshold but not Bonferroni (6.97).
That cluster is the tightest-localising locus in the whole scan: an r2 >= 0.8 interval of
**14 kb**, the smallest of the eleven, and one of the six flagged `localises = TRUE`. The next
tightest is IV:13.409 Mb at 36 kb. The two other chromosome III loci sprawl over 3,765 kb each
and are flagged FALSE.

Distances to the *sid-2* focal variant III:13,680,248: the right-arm cluster peak is 0.962 Mb
proximal, the discarded singleton 7.715 Mb. The Figure 1C caption quotes the 7.715 Mb figure
when it says no Bonferroni marker is near the NIL-resolved interval -- which remains true, and
is now true for a better reason than distance.

None of this puts the association on *sid-2*. The 37 kb NIL interval holds 88 markers whose
best -log10 p is 0.98, and the T96K variant itself is 0.62, ranking 18,662 of 64,423 on the
chromosome. The claim the scan supports is the right ARM of chromosome III, at the resolution
231 phenotyped strains buy.

### Eigen-level support by chromosome

Counts below are markers at or above the eigen threshold, Bonferroni markers included.

| chr | Bonferroni | eigen total | max -log10 p | eigen span (Mb) |
|---|---|---|---|---|
| III | 1 | 39 | 8.68 † | 3.086-12.800 |
| IV | 8 | 413 | 8.84 | 0.773-17.098 (concentrated in two blocks) |
| V | 0 | 5 | 5.40 | 0.605-17.857 |
| X | 1 | 8 | 7.83 | 4.876-14.753 |

† chromosome III's maximum is the unsupported singleton at 5.97 Mb. Excluding it, the
chromosome III maximum is 6.31 at the right-arm cluster, which is below Bonferroni. Read the
"1" in the Bonferroni column with that in mind: on the local-support criterion chromosome III
carries no Bonferroni-supported association.

Chromosome IV carries two distinct blocks: 13.379-13.415 Mb (21 eigen markers, max 7.49) and
15.103-16.120 Mb (250 eigen markers, max 8.84). The caption's "two clusters near 13.41 and
15.32 Mb" glosses over IV:15986656, which is a third position 0.66 Mb from the 15.32 Mb cluster.

### The replacement sentences

> To ask whether this variation has a genetic basis, we performed a genome-wide association
> scan on the variance-stabilized pos-1 response (464,045 markers, n = 231 strains; Figure 1C).
> Ten markers exceeded a Bonferroni threshold of -log10 p = 6.97 and 465 exceeded an
> eigenvalue-based threshold of 4.60. Eight of the ten Bonferroni markers fall on chromosome
> IV, in a cluster at 13.40-13.41 Mb and a stronger cluster at 15.32 Mb that carries the
> genome-wide maximum (-log10 p = 8.84), together with a single marker at 15.99 Mb; a ninth
> lies on chromosome X at 4.88 Mb (7.83). The tenth, on chromosome III at 5.97 Mb, has no
> supporting marker within 100 kb of it although 628 are genotyped there, and it falls in the
> low-recombination centre of the chromosome, where an association should be flanked by the
> markers in linkage disequilibrium with it; we therefore do not carry it forward. The
> chromosome III signal we do carry forward is a cluster of 15 markers spanning 12.70-12.80 Mb
> on the right arm, which clears the eigenvalue-based threshold but not Bonferroni (peak
> -log10 p = 6.31) and is the most tightly localised locus in the scan. Each of these
> associations acts in the same direction, with the minor allele conferring resistance, and all
> but the chromosome IV marker at 15.99 Mb are rare in the panel (allele frequency 0.09 to
> 0.21).

If you want it shorter, the first two sentences plus "Eight of the ten Bonferroni markers fall
on chromosome IV and a ninth on chromosome X at 4.88 Mb; the tenth, on chromosome III at
5.97 Mb, has no supporting marker within 100 kb and is not treated as an association. A cluster
on the right arm of chromosome III (12.70-12.80 Mb) clears the eigenvalue-based threshold
(Figure 1C)" carries the result.

Whether to connect the right-arm cluster to *sid-2* is still unsettled and is deliberately not
written into the draft above. It is 0.962 Mb from the focal variant, the NIL interval itself is
flat (88 markers, best -log10 p 0.98; T96K 0.62, rank 18,662 of 64,423), and the Figure 1C
caption stresses that no Bonferroni marker is near the NIL-resolved interval. The defensible
claim is the right ARM, which is also what the JU1793 x JU2466 cross supports independently --
though the cross peaks at 13.78 Mb, roughly 1 Mb distal to the GWA cluster, so the two associate
the same arm without overlapping. Present them as two independent associations to one arm, not
as an intersection.

Intervals ARE now available for the eigen loci -- `TABLE_gwas_qtl_intervals_eigen.tsv` carries
r2 intervals at four cutoffs -- but not for *sid-2* itself, which needs Tier 1 of
`SYNC_MANIFEST.md`.
