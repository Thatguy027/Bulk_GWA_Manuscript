# Review of manuscript paragraph 2 -- the designed DNA mixture

Sources: `METHODS.txt` (Designed DNA mixture), `FIGURE_CAPTIONS.txt`
(SUPP_FIG_XX_dilution_validation), `scripts/SUPP_FIG_XX_dilution_validation.R`,
`supplemental_data/deconvolution/{dilution_design,dilution_strain_sets,dilution_predictions_bcref,dilution_predictions_poolref}`.
Every statistic below was recomputed from the deposit; the recomputation reproduces the
script's pinned panel-C fractions to 0.0000.

## 1. The sizes are wrong; "two pools" is right

CORRECTED after author input. The first pass through this file said the experiment "was four
pools, not two". That conflates two things. The DILUTION experiment is two sets titrated into
each other, exactly as the paragraph says; sets A and D were sequenced in the same run but are
not part of it. Nothing needs changing on that point.

What does need changing is the sizes. 174 wild isolates were divided into four sets and genomic
DNA pooled within each; the titration used sets B and C.

| set | strains | isotypes | in the dilution experiment? |
|---|---|---|---|
| A | 46 | 46 | no -- sequenced pure only |
| B | 46 | 46 | **yes** |
| C | 40 | 38 | **yes** |
| D | 42 | 42 | no -- sequenced pure only |

So "two pooled populations composed of ~48 strains each" should be **sets B and C at 46 and 40
strains**, 84 isotypes combined, which is the reference panel C uses. The "~48 each" is wrong
for both and hides that the two pools are of different sizes.

## 2. DNA was pooled, not populations

`METHODS.txt` and the script both say genomic DNA from each set was pooled. "We extracted DNA
from two pooled populations" reads as pooled worms -- one culture per set, grown together --
which would introduce competition and differential growth that this experiment deliberately
avoids. If the strains were in fact grown separately and their DNA combined, say so: it is the
reason this validation isolates inference error from biology.

## 3. "Summed ... to establish an input pooled population ratio" inverts the logic

Summing per-strain frequencies within a set gives the INFERRED share of that set. The INPUT
ratio comes from the design: a two-fold doubling series of set B against a fixed 1 uL of set C,
made up to 10 uL with water, giving designed B fractions of 0.091, 0.167, 0.286, 0.444, 0.616,
0.762 and 0.865. Stocks measured 100 ng/uL (B1) and 99.9 ng/uL (C1), so mass fraction equals
volume fraction to within 2.5e-4 -- the equal-concentration assumption is verified, not
assumed. The sentence should say that because the per-strain input frequencies within each pool
were unknown, recovery was assessed at the level of the pools.

## 4. The 2.85% is not one of the attested statistics

Deviation of recovered from designed set-B fraction, panel C reference (recomputed):

| sample | designed B | recovered B | deviation |
|---|---|---|---|
| BC1 | 0.091 | 0.1713 | +0.0803 |
| BC2 | 0.167 | 0.2095 | +0.0427 |
| BC3 | 0.286 | 0.2847 | -0.0012 |
| BC4 | 0.445 | 0.4775 | +0.0328 |
| BC5 | 0.616 | 0.6200 | +0.0044 |
| BC6 | 0.762 | 0.7706 | +0.0085 |
| BC7 | 0.865 | 0.8433 | -0.0217 |

- mean absolute deviation **2.74%**
- RMSE **0.0376** in fraction units (the figure and Methods quote 0.038)
- largest deviation **+0.080** at BC1, whose 0.1 uL of set B is the smallest volume in the series
- dropping BC1: RMSE 0.024
- mean signed bias **+0.0208**, toward set B
- Pearson r **0.99716**, Spearman +1

No computation on the deposit returns 2.85%. The closest quantity to "average discrepancy" is
the mean absolute deviation, 2.74%. If the intended statistic is the one the figure and
Methods already quote, it is RMSE 0.038 (3.8%), not 2.85%.

NOTE that `check_manuscript_numbers.py` PASSES the paragraph as written. `2.85` is attested in
`FIGURE_REPORT.md` -- as a per-mille leakage value in an unrelated analysis -- and `48` is
attested as a strain count at the chromosome IV peak. The checker matches numbers, not their
meaning, so a wrong value that happens to appear elsewhere in the report sails through. Second
paragraph in a row where the guard passes and the claim is still wrong.

## 5. The assumption-free statistic the script prefers

The script names metric 1 "the metric to quote": only B and C are titrated against each other,
so the pair's combined share of the pool must be constant across BC1-BC7 whatever the intended
ratios were. Recomputed on the pool reference: B+C mean 0.8226, sd 0.0066 (0.80% of the mean),
maximum deviation 0.0106 (1.29%). This bounds inference error without invoking the design at
all, which suits a paragraph whose premise is that the strains were not pooled at known
frequencies. Worth a second sentence.

## 6. Figure citation

Second supplement, `SUPP_FIG_XX_dilution_validation` -- S2 if numbered in `FIGURE_CAPTIONS.txt`
order. The accuracy comparison is **panel E** (recovered against designed B fraction). Panel C
is the same series renormalised to the B+C reference; panel A is the pure pools; panel B is the
pool-wide reference. Cite S2E for this claim, not the whole figure.

## 7. Two caveats the paragraph should not carry but Methods must

- There are no replicate dilutions, so pipetting error is unreplicated and enters the
  comparison in full.
- Total input DNA was not constant: 11 ng at BC1 to 74 ng at BC7, because only the B volume was
  varied. Input mass and B fraction are perfectly confounded across the series.

## Suggested revision

Next, we wanted to establish an effective strain pooling and frequency inference strategy for
genetically diverse C. elegans isolates. To this end, we performed a small-scale experiment in
which we pooled genomic DNA from two sets of wild isolates, sets B (46 strains) and C (40
strains), and combined them across a seven-step titration at known mass ratios, sequenced the
resulting libraries, and inferred individual strain frequencies. Because the strains within
each pool were not combined at known frequencies, we assessed recovery at the level of the
pools, summing the inferred strain frequencies within each set to obtain its share of the
mixture. We found that the inferred shares tracked the designed ratios closely (Pearson
r = 0.997, root mean squared error 0.038 in fraction units), with the largest deviation at the
smallest-volume step of the series (+0.080) (SUPP FIGURE, panel E).

If you want the design-free bound as well, append: Consistent with this, the combined share of
the two titrated pools stayed constant across the series (mean 0.823, s.d. 0.007, maximum
deviation 1.3%), which bounds the inference error without reference to the intended mixing
ratios.
