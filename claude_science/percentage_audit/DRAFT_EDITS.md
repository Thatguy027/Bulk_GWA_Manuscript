# Draft edits — outstanding percentage and interval corrections

Three sentences. Everything else from the earlier audit is already applied.

Every replacement number was recomputed from
`supplemental_data/hatching_assays/` on 2026-09-21 by
`scripts/DIAG_hatching_interval_audit.R`. CURRENT text is quoted verbatim from
`Bulk Paper.pdf` as exported 21 Sep, so the spacing oddities (`( Figure 4a )`)
are the PDF extractor's, not necessarily your document's.

Why the intervals move at all: `METHODS.txt:1352` declares Wilson score
intervals and all seven figure scripts compute Wilson, but these four
allele-swap intervals were Clopper-Pearson. The NIL intervals elsewhere in the
draft were already Wilson and are correct.

---

## Sentence 1 — allele swaps (Results, SID-2 section)

Five changes: four intervals restated as Wilson, and **4.5% becomes 5.4%**.

**CURRENT**

> Introducing the 96K allele into the JU1793 strain lowered the strain’s hatching rate from 94.8% (95% CI 90.9-97.4) to 53.1% (95% CI 46.1-60.1), while introducing the 96T allele into the JU2466 strain raised its hatching rate from 4.5% (95% CI 2.8-7.0) to 18.4% (95% CI 13.5-24.2) ( Figure 4a ).

**CORRECTED**

> Introducing the 96K allele into the JU1793 strain lowered the strain’s hatching rate from 94.8% (95% CI 91.0-97.1) to 53.1% (95% CI 46.3-59.8), while introducing the 96T allele into the JU2466 strain raised its hatching rate from 5.4% (95% CI 3.0-9.4) to 18.4% (95% CI 13.8-24.1) (Figure 4A).

| what | from | to | why |
|---|---|---|---|
| JU1793[96T] CI | 90.9-97.4 | **91.0-97.1** | Wilson, 202/213 |
| JU1793[96K] CI | 46.1-60.1 | **46.3-59.8** | Wilson, 110/207 |
| JU2466 value and CI | 4.5% (2.8-7.0) | **5.4% (3.0-9.4)** | wrong denominator — see below; Wilson, 11/204 |
| JU2466[96T] CI | 13.5-24.2 | **13.8-24.1** | Wilson, 40/217 |
| figure callout | Figure 4a | **Figure 4A** | casing |

**The 4.5% to 5.4% change is the only one that alters a claim.** 4.5% is the
pooled rate of JU2466_A and JU2466_B (19/419). The 96T edit was made in
JU2466_A — wSZ206's motif field reads `JU2466_A[NxT]` — and the p value the
Figure 4A caption reports for this comparison, 3.9e-05, is the JU2466_A test;
pooled would give 3.3e-08. So the sentence currently compares an A-derived edit
against an A+B average while the figure compares A against A. JU2466_A alone is
5.4% (11 of 204). The two isolates genuinely differ: B is 3.7% (8 of 215).

The odds ratios later in the same paragraph already use JU2466_A (4.0-fold,
95% CI 1.9-8.8), so this change makes the percentage agree with them.

---

## Sentence 2 — N2 at the standard 50% dose (Results, before Figure S9)

**CURRENT**

> This is likely because at this dose N2 fails to hatch 94% of its embryos (6% hatched, 95% CI 3.2-10.0) and there is little room for further sensitization ( Figure S9 ).

**CORRECTED**

> This is likely because at this dose N2 fails to hatch 94% of its embryos (6.0% hatched, 95% CI 3.5-10.0) and there is little room for further sensitization (Figure S9).

| what | from | to | why |
|---|---|---|---|
| CI | 3.2-10.0 | **3.5-10.0** | Wilson, 13/217 |
| value | 6% | **6.0%** | one decimal, to match every other rate in the paragraph |

---

## Sentence 3 — N2 at the 25% dose (Results, Figure 4B)

**CURRENT**

> we found that the 96K allele lowered the N2 hatching rate from 32.3% (95% CI 26.4-38.7) to 4.0% (95% CI 2.6-6.0) (two independently edited lines, 4.4% and 3.7%; Figure 4B ),

**CORRECTED**

> we found that the 96K allele lowered the N2 hatching rate from 32.3% (95% CI 26.4-38.7) to 4.0% (95% CI 2.7-6.0) (two independently edited lines, 4.4% and 3.7%; Figure 4B ),

| what | from | to | why |
|---|---|---|---|
| pooled 96K CI | 2.6-6.0 | **2.7-6.0** | Wilson, 23/568 |

The other numbers in this sentence are already correct: 32.3% (95% CI
26.4-38.7) is Wilson on 71/220, and 4.4% and 3.7% are the two lines separately.

---

## If you would rather use find-and-replace

Seven fragments, each unique in the draft:

```
90.9-97.4        ->  91.0-97.1
46.1-60.1        ->  46.3-59.8
4.5% (95% CI 2.8-7.0)  ->  5.4% (95% CI 3.0-9.4)
13.5-24.2        ->  13.8-24.1
3.2-10.0         ->  3.5-10.0
2.6-6.0          ->  2.7-6.0
Figure 4a        ->  Figure 4A
```

Search `2.6-6.0` rather than `4.0% (95% CI 2.6-6.0)`; the longer string may not
match if your document spaces the parenthesis differently from the PDF export.
