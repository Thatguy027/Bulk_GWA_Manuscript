# Review of manuscript paragraph 3 -- the MIP-seq comparison

Sources: `METHODS.txt` (accuracy against MIP-seq; sequencing-depth requirements),
`FIGURE_CAPTIONS.txt` (FIGURE 1, FIGURE 1 VARIANT Figure1_pos1,
SUPP_FIG_XX_downsample_per_sample), `plots/`, and PubMed for the citation.

## 1. The citation is right; the expansion of MIP-seq is not

Verified: Webster AK, Chitrakar R, Powell M, Chen J, Fisher K, Tanny RE, Stevens L, Evans K,
Wei A, Antoshechkin I, Andersen EC, Baugh LR (2022) "Using population selection and sequencing
to characterize natural variation of starvation resistance in Caenorhabditis elegans." eLife
11:e80204, doi:10.7554/eLife.80204. So "Webster et al. 2022" is correct.

But MIP-seq is **molecular inversion probe** sequencing, which is what the paper's abstract
says. `METHODS.txt:192` expands it as "multiplexed inverse PCR sequencing", which is wrong and
should be fixed there before it propagates into the draft.

Note also that the paper applied MIP-seq to 100 wild strains, while the repo's dataset is
described as 102 strains and 23 samples and the comparison itself runs on n = 98. Worth
reconciling if the text quotes a strain count.

## 2. The slope is measured DURING starvation, not on recovery afterwards

`METHODS.txt` and the Figure 1A caption both define it as the per-strain rate of change in pool
frequency during L1 starvation: the slope of frequency change regressed on day, with day 1 as
the baseline and day 17 excluded, averaged over five replicate arms, computed in closed form.
"Recovery slopes after L1 arrest" describes a different measurement.

N2 is excluded throughout, which is why n = 98 rather than 102: as the reference strain its
genotype shares sites with every other strain in the pool, so its NNLS frequency is not
identified on the same footing. That exclusion is a deliberate choice and belongs in the
sentence.

## 3. "The same analysis": what the code actually does

**VERIFIED (2026-09-08), and METHODS.txt now states it.** Re-read line by line at author
request. `platform_slopes()` builds one frame, applies `filter(!baseline, day != 17)` and the
`day == 1` baseline join once to both columns, forms `delta_d1_wgs` and `delta_d1_mip` in a
single `mutate()`, and passes each through the same `fit()` -> same `ols_slope(day, .)` under
the same `group_by(replicate, strain)`. The only asymmetry between the two paths is which
column is named. `ols_slope` centres x before the closed-form fit, so the intercept is free.
Confirmed: the two experiments are analysed identically, and `METHODS.txt` now says so
explicitly rather than leaving it to be inferred from the code.

CORRECTED after reading `scripts/Figure1_common.R:213`. `platform_slopes()` puts both platforms
through one identical path -- `delta_d1_wgs = frq - base_frq` and
`delta_d1_mip = published_frq - base_pubfrq`, then the same `ols_slope(day, ...)` per replicate
arm. One transformation applied to both datasets, which is what makes the comparison fair. The
earlier claim in this file that parity "is not established" was wrong.

The only live issue is what "the same" modifies. Webster et al.'s published Slope trait is
defined differently from this one in two respects (their Methods): days 1, 9 and 13 recovery
samples normalized by day-1 frequencies and LOG2 TRANSFORMED, and a line fit with INTERCEPT
FORCED TO 0. This repository uses a raw difference from day 1 and a free-intercept OLS slope
(`ols_slope` centres x, so the intercept is not constrained). Day 1 as baseline, day 13 as the
final point and day 17 excluded all match -- they dropped day 17 because mortality by then left
recovery cultures with few progeny, and their PC1 metric is the one that retains it.

Consequence: none for rho = 0.97, which is a rank correlation between two identically treated
columns. But the slope values are NOT on the scale of the published Slope trait, so they should
never be compared numerically against it. Phrase the sentence as "the same slope-based
transformation to both datasets" so no reader checks it against the log2 definition and finds a
mismatch that is not there.

## 4. Figure 1 no longer contains the panels this paragraph needs

`plots/` holds `Figure1_pos1` only -- there is no `plots/Figure1.pdf` -- so the manuscript's
Figure 1 is the pos-1 variant, whose panels are:

| panel | content |
|---|---|
| A | MIP-seq against NNLS slope, n = 98, rho = 0.97, p < 1e-4 |
| B | distribution of the 2023 pooled pos-1 vst phenotype, n = 231 |
| C | GEMMA LOCO scan for that phenotype, 464,045 markers, n = 231 |

The per-sample agreement panel (23 samples, median 0.84, minimum 0.517) and the
sequencing-depth panel were panels B and C of the OLDER Figure 1 and are not in the curated
figure. So:

- rho = 0.97 -> **Figure 1A**, correct as cited.
- the downsampling claim -> **SUPP_FIG_XX_downsample_per_sample**, panel B (the slope-level
  curve with 0.5x restored). It cannot be cited to Figure 1.

## 5. The correlation does degrade at 1x

Slope-level Spearman rho against the MIP-seq slopes:

| depth | rho |
|---|---|
| 0.25x | 0.833 |
| 0.5x | 0.870 |
| 1x | 0.851 |
| 3x | 0.903 |
| 5x | 0.906 |
| 10x | 0.904 |
| full | 0.97 |

1x gives 0.851 against 0.97 at full depth -- a drop of about 0.12 -- so "did not degrade
substantially" overstates it. What the series supports: agreement is already 0.83-0.87 at a
quarter to half an x, reaches 0.90 by 3x and saturates there, and never recovers the full-depth
0.97 within the range tested.

Two things the caption flags that the text should respect:

- The 0.5x value (0.870) sits ABOVE the 1x value (0.851). That inversion is why 0.5x is omitted
  from the Figure 1C variant, and the caption reads it as noise in one downsampling draw.
  Quoting 1x alone as the floor is therefore slightly arbitrary.
- "Quoting the depth requirement from B alone overstates how well any single sample is
  measured." At the per-sample level the median rises 0.724 at 0.25x to 0.827 at 5x, against a
  full-depth median of 0.835, and samples span roughly 0.43 to 0.88 at any fixed depth.

## 6. The closing claim needs a boundary

"Accurately infer pooled strain frequencies" is carried by three results that each qualify it:
median per-sample agreement 0.84 with a minimum of 0.517 (and rep5_d13 stuck at 0.428-0.478 at
every depth, so its disagreement is not a depth problem); the simulation's 1x spread of
0.52-0.91; and pure-pool recovery of 0.775-0.853 in the dilution experiment. The defensible
version scopes the claim to what the figures measure -- per-strain slopes averaged over
replicate arms -- rather than to per-sample frequencies.

## 7. The number checker passes this paragraph too

Third in a row. `0.97` is attested, `2022` is read as a year and ignored, `1X` is a bare 1 and
ignored. Nothing in the paragraph's numbers is wrong, but the guard would not have caught the
figure mis-citation or the depth overstatement either.

## Suggested revision

As a final test of the NNLS approach, we compared it with a previously published
pooled-phenotyping experiment that used molecular inversion probe sequencing (MIP-seq) to infer
strain frequencies (Webster et al. 2022). To make the comparison, we performed whole-genome
sequencing on the same samples. For each strain, we took the rate of change in pool frequency
across days of L1 starvation as the slope of frequency change regressed on day, using day 1 as
the baseline and excluding day 17, and averaged those slopes over the five replicate arms; N2
was excluded throughout, because as the reference strain its genotype shares sites with every
other strain in the pool. We found the NNLS-derived slopes to be highly correlated with those
derived from MIP-seq (Spearman's rho = 0.97, n = 98 strains, p < 1e-4) (Figure 1A). Agreement
was largely retained as we downsampled reads, reaching rho = 0.85 at 1x and saturating at 0.90
from 3x upward (SUPP FIGURE). Taken together, these three validation experiments gave us
confidence that we could infer per-strain frequency dynamics accurately enough to treat them as
a quantitative phenotype.

If you would rather keep the original closing sentence, the honest hedge is to add the
per-sample number: agreement is lower for individual samples than for the averaged slopes
(median rho = 0.84 across 23 samples), which the supplement shows.
