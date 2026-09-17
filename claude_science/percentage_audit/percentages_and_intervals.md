# Every percentage and confidence interval in the draft, recomputed

Audited 2026-09-17 against repository HEAD. Sources: `supplemental_data/hatching_assays/`
(`ju_allele_swaps_hatching.csv`, `n2_allele_swaps_hatching.tsv`, `nil_series_hatching.tsv`).
Recomputation in `plots/diagnostics/CACHE_hatching_ci_audit.tsv`.

## The headline: the draft mixes two interval conventions

`METHODS.txt:1352` declares **"Binomial proportions are reported with Wilson score 95%
intervals"**, and seven figure scripts implement Wilson (`Figure4_sid2.R`,
`n2_swap_panels.R`, `Figure3_common.R`, `SUPP_FIG_XX_nil_hatching_full.R`,
`SUPP_FIG_XX_sid2_allele_swaps_full.R`, `SUPP_FIG_XX_sid2_ortholog_conservation.R`,
`DIAG_cross_lod_conventions.R`).

The draft text does not follow it consistently:

- **allele-swap intervals are Clopper-Pearson** (every one matches the exact column to
  the decimal, and none matches Wilson)
- **NIL / plate intervals are Wilson** (these match)

METHODS itself is clean — every value in `EMBRYO HATCHING ASSAYS` reproduces exactly on
Wilson. Only the draft prose diverges, so the fix is to restate the allele-swap intervals
as Wilson rather than to change the convention.

## Full check

| # | draft text | draft value | Wilson | Clopper-Pearson | verdict |
|---|---|---|---|---|---|
| 1 | JU1793[96T] | 94.8% (90.9-97.4) | **91.0-97.1** | 90.9-97.4 | value right; restate CI as Wilson |
| 2 | JU1793[96K] (wSZ200) | 53.1% (46.1-60.1) | **46.3-59.8** | 46.1-60.1 | value right; restate CI as Wilson |
| 3 | JU2466[96K] | 4.5% (2.8-7.0) | 2.9-7.0 (pooled) | 2.8-7.0 (pooled) | **WRONG DENOMINATOR — see below** |
| 4 | JU2466[96T] (wSZ206) | 18.4% (13.5-24.2) | **13.8-24.1** | 13.5-24.2 | value right; restate CI as Wilson |
| 5 | N2[96T], 25% dose | 32.3% (2.61-38.9) | **26.4-38.7** | 26.1-38.9 | **lower bound is 26.1 with the decimal shifted** |
| 6 | N2[96K] two lines pooled | 4.0% (2.6-6.0) | **2.7-6.0** | 2.6-6.0 | value right; restate CI as Wilson |
| 7 | N2 at 50% dose | "6% unhatched (3.2-10.0)" | **3.5-10.0** | 3.2-10.0 | **6.0% is the HATCHED fraction, not unhatched** |
| 8 | NIL JU2466 | 35.6% (29.7-44.0) | 35.7%, **29.7-42.0** | 29.5-42.2 | **upper bound wrong: 44.0 should be 42.0**; value 35.7% |
| 9 | NIL JU1793 | 99.4% (96.9-99.9) | 96.9-99.9 | 96.9-100.0 | correct, Wilson |
| 10 | NIL wSZ153 | 39% (32.8-45.6) | 32.8-45.6 | 32.5-45.8 | correct, Wilson |
| 11 | wSZ208 N94A in JU2466 | 42.5% | 42.51% (88/207) | — | correct |
| 12 | N2 96K lines separately | 4.4% and 3.7% | 4.40, 3.73 | — | correct |
| 13 | wSZ206 vs N94A contrast | 42.5% vs 18.4%, OR 3.3 | — | — | percentages correct |

## Item 3 in detail — JU2466 is 5.4%, not 4.5%

Not a digit transposition. `results_review` recorded it as one (`N4`) and that has
propagated through several sessions; it is wrong.

| | n | hatched | rate | Wilson | exact |
|---|---|---|---|---|---|
| JU2466_A[96K] | 204 | 11 | **5.4%** | 3.0-9.4 | 2.7-9.4 |
| JU2466_B[96K] | 215 | 8 | 3.7% | 1.9-7.2 | 1.6-7.2 |
| A + B pooled | 419 | 19 | 4.5% | 2.9-7.0 | **2.8-7.0** |

4.5% is the pooled rate and its Clopper-Pearson interval is 2.8-7.0, matching the draft
exactly. It is arithmetically correct and is the **wrong comparator**:

- `Figure4_sid2.R:233` filters `strain %in% c("JU1793","JU2466_A","wSZ200","wSZ206")`
- the Fisher test backing the sentence (`:258`) is `JU2466_A` vs `wSZ206`
- `MANUSCRIPT_CAPTIONS.txt` reports 5.4% (n = 204)
- wSZ206's `glycosylation motif` field reads **`JU2466_A[NxT]`** — the edit was made in A

So the text compares an A-derived edit against an A+B average while the figure and p value
compare A against A, and the two isolates genuinely differ. **Use 5.4% (11/204,
Wilson 3.0-9.4%).** Recorded in `METHODS.txt` under EMBRYO HATCHING ASSAYS and in the
Figure 4 caption.

## Replacements for the draft

| find | replace |
|---|---|
| `from 94.8% (95% CI 90.9-97.4) to 53.1% (95% CI 46.1-60.1)` | `from 94.8% (95% CI 91.0-97.1) to 53.1% (95% CI 46.3-59.8)` |
| `raised its hatching rate from 4.5% (95% CI 2.8-7.0) to 18.4% (95% CI 13.5-24.2)` | `raised its hatching rate from 5.4% (95% CI 3.0-9.4) to 18.4% (95% CI 13.8-24.1)` |
| `from 32.3% (95% CI 2.61-38.9) to 4.0% (95% CI 2.6-6.0)` | `from 32.3% (95% CI 26.4-38.7) to 4.0% (95% CI 2.7-6.0)` |
| `(6% unhatched, 95% CI 3.2-10.0)` | `(6.0% hatched, 95% CI 3.5-10.0)` |
| `JU2466 (35.6% hatched; 95% CI 29.7-44.0%)` | `JU2466 (35.7% hatched; 95% CI 29.7-42.0%)` |

Items 9-13 need no change.

## The odds ratios, recomputed

The draft states the 96T allele "raises hatching odds 16.1-fold in JU1793
(95% CI 8.2-34.8) but only 4.8-fold in JU2466 (95% CI 2.6-9.1)".

**Method identified:** Fisher conditional MLE with the exact conditional interval
(R `fisher.test`). JU1793 reproduces to the decimal — conditional OR 16.08 -> 16.1,
CI 8.16-34.75 -> 8.2-34.8. Neither a sample odds ratio nor a Woolf interval matches.

| comparison (96T vs 96K, pos-1 food) | OR | exact 95% CI | Fisher p |
|---|---|---|---|
| JU1793 (202/11 vs 110/97) | 16.08 | 8.16-34.75 | 1.7e-24 |
| **JU2466_A alone (40/177 vs 11/193)** | **3.95** | **1.92-8.82** | **3.9e-05** |
| JU2466 A+B pooled (40/177 vs 19/400) | 4.74 | 2.60-8.93 | 3.3e-08 |
| JU2466_B alone (40/177 vs 8/207) | 5.83 | 2.60-14.79 | 9.1e-07 |

**The caption already settles which is right.** Figure 4A reports "p = 1.7e-24 and
p = 3.9e-05". The second is the **JU2466_A** test; the pooled comparison gives 3.3e-08.
So the figure uses A and the text's odds ratio uses A+B — the two currently describe
different comparisons.

The draft's stated JU2466 interval (2.6-9.1) does not reproduce exactly under any method
tried: the pooled conditional interval is 2.60-8.93 (lower bound matches, upper does not)
and its point estimate 4.74 rounds to 4.7, not 4.8, though the *sample* OR pooled is 4.76.
Unlike the JU1793 figure, which reproduces perfectly, this one appears to mix estimators
or predates a data revision. Either way it is superseded.

**Replacement:** `only 4.8-fold in JU2466 (95% CI 2.6-9.1)` ->
`only 4.0-fold in JU2466 (95% CI 1.9-8.8)`.

Note the contrast gets *stronger*, not weaker: 16.1 against 4.0 rather than 16.1 against 4.8.

## The background-dependence claim now has a test

The draft infers background dependence by comparing two separately estimated ratios. Their
intervals overlap (8.16-34.75 against 1.92-8.82), so that comparison does not establish it.
A logistic model of hatched/unhatched on background, allele and their interaction does:

- ratio of odds ratios **4.08 (95% CI 1.56-10.71)**
- Wald p = 0.0042, likelihood-ratio p = 0.0049

Recorded in `METHODS.txt` under EMBRYO HATCHING ASSAYS. Worth stating in the Results in
place of, or alongside, the two-interval comparison.

---

# The allele-swap paragraph, rewritten

Every interval below is Wilson, per `METHODS.txt:1352`. Odds ratios are the Fisher
conditional MLE with its exact interval, per the convention identified above.

## Tracked changes

We generated reciprocal allele-swap strains to test whether this variant was contributing
to *pos-1* RNAi responses. Introducing the 96K allele into the JU1793 strain lowered the
~~strains~~ **[strain's]** hatching rate from 94.8% (95% CI ~~90.9-97.4~~ **[91.0-97.1]**)
to 53.1% (95% CI ~~46.1-60.1~~ **[46.3-59.8]**), while introducing the 96T allele into the
JU2466 strain raised its hatching rate from ~~4.5%~~ **[5.4%]** (95% CI ~~2.8-7.0~~
**[3.0-9.4]**) to 18.4% (95% CI ~~13.5-24.2~~ **[13.8-24.1]**) (Figure ~~4a~~ **[4A]**).
Therefore, the 96T reference allele produces a weaker RNAi response than the alternate 96K
allele. These results suggest that this allele interacts with other variants in the parental
strain's background to produce background-dependent effects: the 96T allele raises hatching
odds 16.1-fold in JU1793 (95% CI 8.2-34.8~~)~~ **[; p = 1.7e-24]**) but only ~~4.8~~
**[4.0]**-fold in JU2466 (95% CI ~~2.6-9.1~~ **[1.9-8.8; p = 3.9e-05]**). **[Fitting
hatching to background, allele and their interaction confirms that the two backgrounds
differ in how much residue 96 matters: the ratio of the two odds ratios is 4.08
(95% CI 1.56-10.71), likelihood-ratio p = 0.0049.]**

## Clean text to paste

> We generated reciprocal allele-swap strains to test whether this variant was contributing
> to *pos-1* RNAi responses. Introducing the 96K allele into the JU1793 strain lowered the
> strain's hatching rate from 94.8% (95% CI 91.0-97.1) to 53.1% (95% CI 46.3-59.8), while
> introducing the 96T allele into the JU2466 strain raised its hatching rate from 5.4%
> (95% CI 3.0-9.4) to 18.4% (95% CI 13.8-24.1) (Figure 4A). Therefore, the 96T reference
> allele produces a weaker RNAi response than the alternate 96K allele. These results
> suggest that this allele interacts with other variants in the parental strain's background
> to produce background-dependent effects: the 96T allele raises hatching odds 16.1-fold in
> JU1793 (95% CI 8.2-34.8; p = 1.7e-24) but only 4.0-fold in JU2466 (95% CI 1.9-8.8;
> p = 3.9e-05). Fitting hatching to background, allele and their interaction confirms that
> the two backgrounds differ in how much residue 96 matters: the ratio of the two odds
> ratios is 4.08 (95% CI 1.56-10.71), likelihood-ratio p = 0.0049.

## Why each change

| change | reason |
|---|---|
| four intervals restated | the draft used Clopper-Pearson; `METHODS.txt:1352` and all seven figure scripts use Wilson |
| 4.5% -> 5.4% | wrong denominator: 4.5% pools JU2466_A and JU2466_B, but the 96T edit was made in A and the caption's p value is the A test |
| 4.8 -> 4.0-fold | follows from the same denominator change |
| p values added | they are already in the Figure 4A caption; giving them in the text lets a reader tie the ratio to the test |
| interaction test added | the two odds-ratio intervals **overlap** (8.16-34.75 against 1.92-8.82), so comparing them does not establish background dependence; the interaction term does |
| "strains" -> "strain's" | possessive, logged as `T3` in the earlier results review |

## If you would rather not add a sentence

The interaction result can ride inside the existing one:

> ... the 96T allele raises hatching odds 16.1-fold in JU1793 (95% CI 8.2-34.8) but only
> 4.0-fold in JU2466 (95% CI 1.9-8.8), a 4.08-fold difference between backgrounds
> (95% CI 1.56-10.71, likelihood-ratio p = 0.0049).

This is the minimum that makes the claim testable rather than inferred from two intervals.
