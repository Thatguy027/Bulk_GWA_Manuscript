# The chromosome V centre in N2 x XZ1516: two loci, not one, and not five

Analysis date 2026-09-11. Data: `data/cross_experiments/N2-XZ_export`
(10 pools, 6 conditions, 249,106 usable chrV markers per pool) and
`data/cross_experiments/JU1793-JU2466_export` for the cross comparison.
Nothing in this report modifies the pipeline or the shipped tables; it
reinterprets them.

---

## 1. Verdict

The chromosome V centre carries **two separable response regions**, and the
thing that separates them is **pos-1**:

| | peak | who responds | who does not |
|---|---|---|---|
| **Left locus** | **V:10.6-10.9 Mb** | mig-6, pos-1, par-1, rpn-12, vha-5 | — |
| **Right locus** | **V:13.6-14.6 Mb** | mig-6, par-1, rpn-12, vha-5 | **pos-1** |

Peak separation is **5.4 cM** (V:10.6 Mb = 24.5 cM, V:14.6 Mb = 29.9 cM on this
cross's map). The bootstrap 95% CIs for the two argmaxes do not overlap
(pos-1 10.0-11.1 Mb; par-1 14.4-15.8 Mb).

**The right locus is not mig-6-specific.** It looks that way because the
repo's "specific" column is defined as the *pos-1 vs mig-6 contrast*
(`scripts/SUPP_FIG_XX_cross_contrast_panels.R:53`), which is large at 13.8 Mb
(LOD 730) precisely because pos-1 does **not** respond there. par-1, rpn-12 and
vha-5 all do, with dfreq 0.10-0.20 and LOD 39-91. mig-6 is the *strongest*
responder at every point on chrV, not the only one.

**Direction is uniform: the N2 allele is enriched and the XZ1516 allele
depleted under every RNAi condition, at both loci.** This is the opposite sign
from the chrI:1.87 Mb mig-6 locus, where XZ1516 is the enriched parent.

**What that sign means (section 6.0): the enriched allele is the RNAi-RESISTANT
one.** So the N2 chrV-centre allele confers resistance and the XZ1516 allele is
the more RNAi-competent one, being purged by the selection.

---

## 2. The baseline has to be subtracted before anything else

The HT115 control is nowhere near 0.5. N2 allele frequency in S29 (HT115, t2),
1 Mb bins:

```
chrV   0.76 0.73 0.71 0.73 0.69 0.65 0.59 0.62 0.54 0.49 0.42 0.36 0.35 0.35 0.28 0.32 0.31 0.33 0.18 0.18 0.01
chrI   0.64 0.69 0.80 0.89 0.90 0.92 0.93 0.93 0.93 0.93 0.91 0.89 0.86 0.74 0.68 0.59
chrIV  0.48 0.43 0.36 0.28 0.20 0.15 0.14 0.12 0.12 0.11 0.12 0.13 0.16 0.23 0.26 0.32 0.41 0.40
```

Every chromosome carries a centre-wide sweep in the control, in different
directions (N2 near-fixed in the I, III and X centres; XZ1516 near-fixed in the
IV centre; a full 0.76 -> 0.01 cline on V). A "deviation from 0.5" test is
therefore unavailable in this cross, and the only usable statistic is
condition-minus-control at the same timepoint. All numbers below are that.

---

## 3. Why two components rather than one

### 3.1 The profile ratios are not constant

If all conditions responded to one shared locus at different strengths, the
ratio of any two control-corrected profiles would be flat across the region. It
is not. Ratio to mig-6 (t2), control-corrected, 200 kb bins:

| Mb | mig-6 | vha-5 | rpn-12 | par-1 | pos-1 | mig-6 t1 |
|---|---|---|---|---|---|---|
| 6.0 | 0.126 | 0.94 | 0.68 | 0.00 | -0.27 | 0.87 |
| 8.4 | 0.203 | 0.81 | 0.57 | 0.29 | 0.21 | 0.94 |
| 10.8 | 0.394 | 0.57 | 0.39 | 0.25 | 0.28 | 0.79 |
| 12.4 | 0.405 | 0.55 | 0.32 | 0.29 | 0.15 | 0.67 |
| 13.8 | 0.365 | 0.57 | 0.35 | 0.33 | -0.07 | 0.88 |
| 15.4 | 0.314 | 0.49 | 0.40 | 0.40 | -0.05 | 0.97 |
| 16.8 | 0.212 | 0.38 | 0.55 | 0.46 | -0.11 | 1.02 |

vha-5 declines monotonically left to right, par-1 rises monotonically, pos-1 is
sharply peaked at 10.6 and negative beyond 13 Mb. **mig-6 t1 / mig-6 t2 is flat
at 0.87 +/- 0.08 across the whole region** — the replicate behaves as one
rescaled profile, the other conditions do not. Monotone trends cannot be
produced by a symmetric saturation effect, so the left-right asymmetry between
vha-5 and par-1 is real.

### 3.2 NMF, k = 2

Non-negative decomposition of the six control-corrected profiles (200 kb bins,
V:5-18 Mb) puts the components at **10.6 Mb** and **13.6 Mb**, residual 0.10 of
total norm:

```
                mig6.2  vha5.2  rpn12.2  par1.2  pos1.2  mig6.1
comp1 (10.6 Mb)  0.170   0.149    0.070   0.006   0.104   0.122
comp2 (13.6 Mb)  0.339   0.173    0.130   0.132   0.000   0.314
```

Peak positions are stable at 10.2-10.6 / 13.6-14.6 Mb when the low-marker bin
at 7.8 Mb is dropped, when bins with < 350 markers are dropped, when mig-6 t1
is dropped, and across five random seeds. The *magnitudes* for
mig-6/vha-5/rpn-12 on comp1 are seed-unstable (identical residual, loadings
ranging 0.19-0.38) because the components overlap — read the loadings as an
ordering, not as values. What is stable in every run is the qualitative
anchoring: **pos-1 loads 0.000 on comp2 and par-1 loads lowest of all
conditions on comp1.**

k = 3 lowers the residual to 0.065 with components at 10.6 / 11.8 / 14.4 Mb
(the 11.8 component being vha-5 + rpn-12). The first SVD PC already explains
97.4% of the matrix, so the profiles are highly collinear and k = 3 is not
defensible from these data alone.

### 3.3 Block bootstrap on the argmax

500 resamples over 200 kb blocks, 1 Mb rolling mean, V:8-17 Mb:

| condition | peak | bootstrap median | 95% CI |
|---|---|---|---|
| pos-1 | 10.8 Mb | 10.6 | **10.0 - 11.1** |
| mig-6 | 11.0 Mb | 11.4 | 10.8 - 12.0 |
| rpn-12 | 11.2 Mb | 11.2 | 10.6 - 11.8 |
| vha-5 | 11.6 Mb | 11.4 | 10.6 - 12.0 |
| par-1 | 14.6 Mb | 14.8 | **14.4 - 15.8** |

pos-1 and par-1 do not overlap. The 20-LOD-drop intervals agree: pos-1
10.22-11.26 Mb, rpn-12 10.34-11.20, mig-6 10.59-10.92, vha-5 10.18-12.43,
**par-1 14.26-16.14** — par-1's interval is disjoint from all the others.

### 3.4 The repo's own tables already say this

`plots/TABLE_cross_qtl_locus_classification.tsv`:

| locus | class | responds |
|---|---|---|
| N2xXZ1516 V:9.28 | general:3 | mig6, rpn12, vha5 |
| **N2xXZ1516 V:10.76** | **general:5** | **mig6, par1, pos1, rpn12, vha5** |
| N2xXZ1516 V:12.21 | general:4 | mig6, par1, rpn12, vha5 |
| N2xXZ1516 V:13.81 | general:4 | mig6, par1, rpn12, vha5 |
| N2xXZ1516 V:16.26 | general:4 | mig6, par1, rpn12, vha5 |

V:10.76 is the only chrV locus where pos-1 responds, and everything from
12.2 Mb rightward is the same pos-1-absent set. The five rows are two
components: the pipeline's peak-splitting rule cuts one broad plateau into
12.21 / 13.81 / 16.26, and V:9.28 is the left tail of the 10.76 locus before
par-1 and pos-1 have risen.

pos-1 dfreq along the chromosome makes the split explicit:
**+0.083 (9.28) -> +0.129 (10.76) -> +0.048 (12.21) -> -0.021 (13.81) -> -0.033 (16.26).**

---

## 4. Replication: only mig-6 replicates

`dfreq.rep1` (timepoint 1) from the repo's own condition table:

| locus | mig-6 t2 / t1 | par-1 t2 / t1 | rpn-12 t2 / t1 |
|---|---|---|---|
| V:10.76 | 0.384 / **0.325** | 0.113 / **0.005** | 0.142 / **0.035** |
| V:12.21 | 0.387 / **0.332** | 0.118 / 0.043 | 0.154 / 0.059 |
| V:13.81 | 0.374 / **0.295** | 0.105 / **-0.012** | 0.139 / 0.022 |
| V:16.26 | 0.323 / **0.252** | 0.124 / 0.014 | 0.112 / 0.022 |

mig-6 reproduces at ~85% of its t2 effect at every locus. par-1 and rpn-12 are
at or below zero at t1. There is no pos-1 or vha-5 pool at t1.

**This is the main weakness of the two-locus conclusion.** The split rests on
pos-1 (present only at t2) being absent where par-1 (present only convincingly
at t2) is strongest. Both halves of the discriminating comparison are
single-pool observations. The *left* locus is replicated (mig-6 t1 dfreq +0.325
there); the claim that a *second, distinct* locus exists at 13.6-14.6 Mb is
not.

---

## 5. Cross comparison

**Left locus: flatly absent in JU1793 x JU2466.** Max LOD per Mb in
HT115g-MIG6g on chrV: 30.7 at 0 Mb, then **0.1-1.5 across 9-13 Mb**, against
806 in N2 x XZ1516. `TABLE_general_vs_specific.tsv` records the same thing at
the locus level.

**Right locus: not cleanly cross-private.** HT115g-POS1g in the JU cross peaks
on chrV at **V:14,534,981, LOD 5.5** — within 10 kb of the N2 x XZ1516
HT115-par1 peak at V:14,545,025. LOD 5.5 over an 884 kb interval is weak and
could be coincidence, but the right locus cannot be treated as absent in the
second cross the way the left one can.

So the cross filter is decisive for the left locus only: the causal variant
there must differ between **N2 and XZ1516** while **JU1793 and JU2466 share**
an allele.

---

## 6. Candidate variants after the cross filter

Filter: `N2 != XZ1516` AND `JU1793 == JU2466` AND all four called AND position
outside a called XZ1516 divergent region
(`/Users/Stefan/UCLA/Genomics_Data/div_regions/FINAL_RECALLED_DIVERGENT_fixMissingHaps.tsv`).
Source VCF `/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213/bcsq.vcf.gz`.

| window | records | N2 != XZ | ...and JU1793 == JU2466 | HIGH-impact, non-divergent |
|---|---|---|---|---|
| V:9.9-11.2 Mb | 23,097 | 3,183 | 2,882 | **10** |
| V:14.4-15.8 Mb | 92,600 | 22,026 | 20,243 | **21** |

### 6.0 The sign of the shift, calibrated

The selection regime itself is **not recorded** -- `METHODS.txt:430` leaves
"selection or competition regime" as `[TO FILL]`. The direction is therefore
established from two anchors that agree:

1. **Allele swap.** On 25% pos-1 RNAi food, unedited N2 (sid-2 96T) hatches
   **32.3%** of plated embryos; two independent N2 lines edited to 96K hatch
   **4.0%** pooled (`METHODS.txt:958-962`). Hatching on pos-1 RNAi is escape
   from the knockdown, so **96K is the RNAi-sensitive allele and 96T the
   resistant one**.
2. **The crosses.** N2 and JU1793 carry 96T; XZ1516 and JU2466 carry 96K
   (`METHODS.txt:545`). At the chromosome III right-arm locus the 96T parent
   rises to near fixation under RNAi in both crosses -- N2 0.677 -> 0.996
   (mig-6), JU1793 0.349 -> 0.964. `METHODS.txt:559` states the same reading
   directly: "the XZ1516 haplotype is the one removed wherever RNAi selects".

**Therefore: allele enriched in an RNAi pool = allele conferring resistance /
a weaker RNAi response. The pools are survivors; responders are purged.**

Applied to chrV: the N2 allele rises 0.39 -> 0.77 under mig-6, so **the N2
chromosome V centre allele is the resistant one and the XZ1516 allele is the
more RNAi-competent one being removed.**

### 6.0.1 This inverts the candidate logic, and rules out the obvious hits

A loss-of-function allele in RNAi machinery carried by **XZ1516** predicts the
XZ1516 haplotype is **protected and enriched**. It is depleted. So
XZ1516-private damaging variants in *required* RNAi genes are the **wrong
sign** and cannot explain either locus. Two candidates are retracted on this
basis and kept below only as documented exclusions:

- **rde-1 G75E** (V:9,991,336, XZ1516, AF 0.030) -- wrong sign. rde-1 is
  required for exogenous RNAi; an XZ1516 hypomorph would confer resistance.
- **set-5 Y221\*** (V:14,701,477, XZ1516-private nonsense) -- wrong sign for
  the same reason.

The right-sign candidate classes are:

- XZ1516 **loss of a negative regulator** of exogenous RNAi -- ERI pathway,
  `rrf-3`, `lin-15b`, `adr-1`/`adr-2`, `eri-6`/`eri-7`, `ergo-1`. Losing a
  repressor raises silencing, which raises sensitivity, which gets the haplotype
  purged.
- an XZ1516 **hypermorph**, or *cis*-regulatory variation raising expression of
  silencing machinery.
- an **N2 hypomorph**. This class is **structurally invisible** to the census
  below: N2 is the reference, so N2 is 0/0 at essentially every site in a
  reference-based VCF. This is a limit on the whole approach, not only on chrV.

Scanning both intervals for right-sign genes: `ergo-1` (V:1.01 Mb) and `tofu-2`
(V:7.05 Mb) are far outside. The only one inside is **`tofu-1`,
V:9,981,147-9,982,766** -- 6.8 kb from rde-1, piRNA biogenesis, and loss of the
piRNA/26G branch can enhance exogenous RNAi. Right sign, weaker mechanistic
link, and at the same boundary position ~780 kb outside the support region.

**Net: after the sign correction neither chrV interval contains a
well-directed RNAi-machinery candidate.** That removes the interval-content
support for an RNAi-pathway explanation and favours the alternatives in
section 7 -- a general fitness locus visible only under strong selection, or
*cis* regulatory variation. The mapping result is untouched: the two-locus
split, the 5.4 cM separation, the disjoint bootstrap CIs and the pos-1 absence
do not depend on sign interpretation.

### Left locus

None of the 10 HIGH-impact survivors is in an RNA-silencing gene
(`C50H2.7`, `fipr-13`, `C08B6.3`, `F15H10.9`, `F15H10.10`, `Y32F6A.4`,
`Y32F6A.5`, `ugt-35`, `C45B11.2`).

Of the canonical RNAi-pathway genes on chrV (`sid-1` 5.12, `mut-14` 5.66,
`sago-1` 6.30, `rde-1` 9.99, `rde-12` 13.67, `mut-15` 15.00, `hrde-2` 15.01,
`prde-1` 15.54, `ergo-1` 1.01 Mb), the nearest to the comp1 peak is **rde-1**:

**rde-1, V:9,987,916-9,991,642, 0% divergent in XZ1516.** XZ1516 carries three
missense changes with N2, JU1793 and JU2466 all reference:

| site | change | CeNDR AF | carriers |
|---|---|---|---|
| V:9,991,336 | **G75E** | **0.0295** (AC 36/1222) | 17 ECA isolates + XZ1516 |
| V:9,991,040 | I160V | 0.2253 | common |
| V:9,990,221 | D369G | 0.2253 | common |

plus a 12 bp XZ1516-private 3'UTR insertion at V:9,987,947.

**RETRACTED as a candidate -- wrong sign (section 6.0.1).** rde-1 is the
Argonaute absolutely required for exogenous RNAi, so an XZ1516 hypomorph
predicts XZ1516 is *enriched*; it is depleted. Two further problems stand
independently of sign: rde-1 sits 0.77 Mb (about 0.6 cM) left of the peak
marker, outside the left bound of the pos-1 bootstrap CI and 230 kb outside the
20-LOD bound; and there is no null allele in it -- the only protein changes are
missense and two of the three are at 22% species frequency. Recorded here
because the variants are real and the exclusion is on direction, not data
quality.

### Right locus

21 HIGH-impact survivors. The only one in a chromatin- or silencing-related
gene is:

**set-5 Y221\*, V:14,701,477, T>G, stop_gained. RETRACTED as a candidate --
wrong sign (section 6.0.1).** N2 = JU1793 = JU2466
reference, XZ1516 homozygous alternate; **CeNDR AF 0.0016 (AC 2/1222), XZ1516 the only carrier in the
species panel**; gene 0% divergent; 156 kb from the HT115-par1
peak; truncates a >1200-residue protein at residue 221. A second XZ1516
frameshift in the same gene (V:14,703,806) is flagged `*frameshift` with an
upstream-dependency pointer.

Other protein changes in silencing genes in the right window, all outside
divergent regions:

| gene | position | XZ1516 change | CeNDR AF |
|---|---|---|---|
| rde-12 | V:13,665,389-13,669,292 | 7 missense (Y48D, D63E, H217Q, K261N, P313S, K779M, ...) + splice_region | — |
| mut-15 | V:15,002,834 / V:15,003,678 | N25D / F178V | 0.218 / 0.041 |
| hrde-2 | V:15,005,999 | R245H | 0.217 |

`mut-15` N25D and `hrde-2` R245H sit at the same ~22% frequency as the common
rde-1 pair, i.e. on one widespread haplotype rather than an XZ1516-private
lesion.

**Excluded as uncallable:** `prde-1` (V:15,538,244-15,547,086) and `hda-1`
(V:14,525,389-14,527,411) are **100% inside called XZ1516 divergent regions**.
prde-1's apparent frameshift and stop_gained at V:15,546,362/15,546,398 are
`./.` in XZ1516 — no call, not a variant. Do not use them.

---

## 7. Limitations

1. **LOD magnitudes are not interpretable in absolute terms.** `scripts/config.yml`
   sets `sample.size: 10000`, `sel.strength: 0.95`, so every statistic assumes
   9,500 independent individuals. The 2-LOD intervals that follow (80-250 kb on
   a chromosome-centre bulk-segregant peak) are not credible; use the peak plus
   a bootstrap CI.
2. **The chrV centre is ~12 cM across 10 Mb** (V:6 Mb = 20.0 cM, V:16 Mb =
   32.4 cM), about 7x lower recombination than either arm. Any selected site in
   the centre drags the whole centre with it, which is why one locus can produce
   a 6 Mb plateau.
3. **Marker density varies 10-fold between 200 kb bins** (346 to 6,958), and
   density is correlated with divergence — hyper-divergent haplotypes carry more
   callable variants, not fewer.
4. **V:15.2-17.6 Mb is 75-100% hyper-divergent in XZ1516**, with 15,386 of
   92,600 records uncalled in the right window. The right locus's CI extends
   into that zone.
5. **Only mig-6 has a timepoint-1 replicate** (see section 4). The two-locus
   split is a t2-only result.
6. **Neither locus is resolved to a gene, and the variant census is
   sign-blind in one direction.** 575 genes lie in the left window and 601 in
   the right. Both of the function-based nominations (rde-1, set-5) are
   retracted on direction (section 6.0.1). Because N2 is the reference, an
   N2-hypomorph explanation cannot be seen in this VCF at all; testing it needs
   parental expression data or an XZ1516-anchored assembly, not a
   reference-based variant call.
7. **Response versus fitness is not separated.** mig-6 is the strongest-selecting
   condition in this cross at every locus; a fitness locus visible only under
   strong selection would look like a graded RNAi-response locus.

---

## 8. Reproducing

Scratch scripts (not tracked):
`/private/tmp/claude-501/-Users-Stefan-UCLA-Projects-RNAi-Manuscript-Bulk-GWA-Manuscript/914ad822-1756-45db-a32c-68329e5f52a4/scratchpad/`
— `chrV.R`, `allchrom.R`, `raw3.R` (200 kb profiles), `ratio.R`, `nmf.R`,
`robust.R`, `fine.R` (50 kb profiles), `boot.R`.

Repo tables used unchanged: `plots/TABLE_cross_qtl_locus_classification.tsv`,
`plots/TABLE_cross_qtl_condition_dfreq.tsv`, `plots/TABLE_general_vs_specific.tsv`.
