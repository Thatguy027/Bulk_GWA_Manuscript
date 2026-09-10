# Results section review

Every checkable claim in the drafted Results section, recomputed from
`supplemental_data/`. Written 2026-09-10 against repository HEAD `8c140d8`.

**42 claims checked: 34 verified, 5 wrong, 3 needing an author decision.**
Nothing at the repository root was modified.

## Files

| file | what it is |
|---|---|
| `results_draft_v1.txt` | the draft as supplied, unmodified |
| `verify_claims.R` | recomputes all 42 claims from the deposit |
| `claims_verified.tsv` | one row per claim: drafted value, computed value, verdict, source file |
| `make_tracked_changes.py` | applies the edit list; asserts every target string matches exactly once |
| `results_tracked_changes.md` | full text, `~~struck~~` -> **[replacement]**, with the reason for each |
| `results_corrected.txt` | clean text with all edits applied |
| `edits_table.tsv` | find/replace pairs for editing the Google Doc |

Reproduce:

```sh
Rscript claude_science/results_review/verify_claims.R
python3 claude_science/results_review/make_tracked_changes.py
```

## The five errors

**1. 1X sequencing depth is not "accurate" (`S3`).** The draft says NNLS "can
accurately infer strain frequencies with as little as 1X sequencing depth". At
1x the mean r² per trait runs 0.51-0.91 with a median of 0.76; the mean first
reaches 0.95 at **10x**. `SUPP_FIG_XX_simulation_depth.R` prints these under
the heading "accuracy at 1x, the depth the text claims", so the script was
already flagging this sentence. Note the simulation was reseeded and replicated
in commit `651a0e1`, after the draft was written — this number moved.

**2. The dilution pools are 46 and 40 strains, not ~48 each (`D1`).**
`dilution_strain_sets.tsv` puts 174 isolates in four sets: A=46, B=46, C=40,
D=42. The titration is set B against set C.

**3. Four pos-1 replicate pools, not two (`P2a`).**
`pos1_2023_sample_frequencies.csv.gz` has pos-1 replicates 1, 2, 3, 4 and
control pools A, B. This one is self-evident from the draft itself: two
replicates give one pairwise comparison, and the next sentence reports six.

**4. Ten RNAi targets, not nine (`C1`).** The pooled panel was phenotyped on
fog-2, mig-6, pos-1, rde-3, ric-3, rpn-12, spe-19, spe-43, unc-39 and vha-5,
93 strains each.

**5. JU2466 hatches at 5.4%, not 4.5% (`N4`).** `ju_allele_swaps_hatching.csv`
on the pos-1 condition: JU2466_A[96K] is 11 hatched of 204 = 5.4%, and
`Figure4_sid2.R` prints 0.054. The 18.4% for the 96T swap is correct. The two
digits appear to be transposed.

## The three author decisions

**Chromosome X is the left arm, not the centre (`P4b`).** The single chromosome
X marker clearing Bonferroni is at 4,875,969 — 28% along a 17.72 Mb chromosome.

**The chromosome III Bonferroni marker goes unmentioned (`P4d`).** The draft
reports chromosome III only as an eigen-threshold locus on the right arm, which
follows this repository's admission rule: `FIGURE_REPORT.md` ("GWAS interval
admission") shows the 5.966 Mb marker clears Bonferroni at 8.68 with **zero**
eigen-passing neighbours within 100 kb, while the 12.70-12.80 Mb right-arm
cluster peaks below Bonferroni at 6.31 with **14**. The framing is defensible,
but a reader comparing the text against Figure 1C sees a red chromosome III
marker the text never accounts for. `manuscript/review/para4_pos1_pilot.md`
already drafts a sentence that disposes of it.

**The N2 rescue is two independent lines (`N6`).** 32.3% -> 4.0% is right; at
the 25% dose the 4.0% pools wSZ203 (4.4%, n=273) and wSZ204 (3.7%, n=295),
which sum to the pooled n=568. Naming both is the
stronger claim and matches the caption. Separately, the draft's "no difference"
at the 50:50 dose is supported: 6.0% vs 3.4%, Fisher p = 0.108.

One further wording note (`M1`): the draft quotes the MIP-seq agreement as
"ρ = 0.97 … p < 1e-4". Recomputed ρ is 0.974 over 98 strains, but Figure 1A
deliberately reports a bootstrap interval instead of a p value —
`Figure1_common.R` notes that a p value against ρ = 0 is not the question being
asked. Quoting the p value reintroduces what the figure dropped.

## Supplement numbering does not match the repository

The draft's supplement numbers are its own sequence. Apparent mapping:

| draft | repository | content |
|---|---|---|
| S1 | S1 | simulation depth |
| S2 | S2 | designed DNA mixture |
| S3 | S5 | sequencing-depth requirement (downsampling) |
| S4 | S6 | replicate reproducibility |
| S5A/B | S7 | plate assay vs pooled response and vs Paaby 2015 |
| S6 | S8 | why these cross parents |
| S7B-C | S9 | expanded view behind Figure 2's tracks |
| S8 | S11 | full N2 dose series |
| S9 | S15 | where T96's pocket sits in the charge distribution |

The draft also cites "Supplemental table X" for the cross QTL. `TABLE_cross_qtl_intervals`
and `scripts/cross_qtl_table.R` exist for this; the script is staged but not
yet committed in the working tree.

## Four numbers the report cannot attest

`scripts/check_manuscript_numbers.py` on the corrected text still flags four
values: `X:4,875,969`, `III:12,353,680`, `V:14,647,434` and `-log10p 7.71`.
All four are correct — this review recomputed each from
`pos1_2023_gemma_loco.csv.gz` and the pooled bundle — but `FIGURE_REPORT.md`
never prints them, so the checker has nothing to match against. The fix is to
surface the pooled peak coordinates in `FIGURE_REPORT.Rmd` rather than to add
them to `manuscript_number_exceptions.txt`, which is for numbers the code does
*not* produce.

## What could not be checked here

Husbandry and protocol claims with no data in the deposit: two generations of
growth, F3 L1 sequencing, ten rounds of intercrossing, "we performed
whole-genome sequencing on the same samples", the Webster et al. 2022 and
Zdraljevic et al. 2025 attributions, and the McEwan et al. 2012 arginine
argument. `METHODS.txt` carries `[TO FILL]` markers for several of these.

## Codebase state

All 23 R figure scripts and all three Python asset builders were run from a
scratch copy of the deposit. **21 of 23 R scripts and 3 of 3 Python builders
completed clean.** The two failures are the documented archive-dependent ones:
`gwas_qtl_intervals.R` (Figure S18) exits with "need
data/genotypes/CeNDR20210121_Plink for LD", and
`SUPP_FIG_XX_sid2_allele_in_panel.R` fails in the `system2` call that shells out
to the same panel. Both are declared in `DATA_AVAILABILITY.md` and the first is
the declared `ARCHIVE_EXEMPT` entry in `scripts/check_repo_invariants.sh`.

`bash scripts/check_repo_invariants.sh` reports **0 failed checks**.
