# Draft edits — outstanding percentage and interval corrections

Four sentences. Everything else from the earlier audit is already applied.

Sentences 1-3 are percentage and interval corrections from the 2026-09-21
hatching audit. Sentence 4 is separate: a clause lost its verb somewhere in
editing, and checking the numbers turned up a p value that was never the
computed one.

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

## Sentence 4 — the panel-level decline (Results, Figure 1B)

A broken clause, a wrong p value, and a term that does not match the figure.

**CURRENT**

> We found that 79% of wild isolates declined in frequency in pos-1 RNAi relative to the control condition, which is significantly if frequency changes reflected noise alone (binomial test, p < 2x10-16) (Figure 1B).

**CORRECTED**

> We found that 79% of wild isotypes (183 of 231) declined in frequency under pos-1 RNAi relative to the control condition, more than would be expected if frequency changes reflected noise alone (binomial test, p = 9.3 × 10⁻²⁰) (Figure 1B).

| what | from | to | why |
|---|---|---|---|
| the clause | "which is significantly if" | **"more than would be expected if"** | the comparative went missing; as written the sentence has no predicate |
| p value | p < 2x10-16 | **p = 9.3 × 10⁻²⁰** | 2.2e-16 is not the result, it is R's printing floor — see below |
| denominator | — | **(183 of 231)** | 231 is the subset carrying a variance-stabilised value, not the 366 isolates in the pool |
| term | wild isolates | **wild isotypes** | matches the Figure 1B caption and METHODS |

**Where `p < 2x10-16` came from.** `scripts/DIAG_pos1_responsiveness_test.R`
printed the test through `format.pval()`, whose default `eps` is
`.Machine$double.eps`; anything smaller prints as `< 2.22e-16` no matter what it
is. The test itself is an exact binomial tail sum and evaluates far below that:
**9.32e-20** for 183 of 231 against 0.5. The script now prints the computed
value, and METHODS.txt quotes it.

**One thing to decide before this goes in.** The sign test is against 0.5, and
the sentence reads that null as "noise alone". Pool frequencies are
compositional and these behave that way — the 231 retained isolates' changes sum
to -0.062, which is zero to within rounding. The two sides are not symmetric:

| | isolates | total mass | mean per isolate |
|---|---|---|---|
| declined | 183 | 0.664 | 0.0036 |
| gained | 48 | 0.603 | 0.0126, largest +0.128 |

A minority gaining that much forces the majority to lose whether or not those
isolates responded at all, so a majority declining is not on its own evidence
against noise. METHODS.txt already says this —

> Pool frequencies are compositional, so this proportion measures decline
> relative to the pool rather than an absolute per-strain response

— and the per-isolate tests in the same section are the claim that does carry
individual responsiveness: **141 of 231, 61%**, by the floor construction.

So the corrected sentence above is safe if the claim is "declines outnumber
gains". If the intended claim is that most isolates responded, use the
per-isolate number instead, something like:

> Under pos-1 RNAi, 183 of 231 wild isotypes declined in pool frequency while 48 gained, with the gains concentrated in a few isolates (largest +0.13); 141 of the 231 (61%) can be shown to have responded individually (Figure 1B).

---

## If you would rather use find-and-replace

Eight fragments, each unique in the draft:

```
90.9-97.4        ->  91.0-97.1
46.1-60.1        ->  46.3-59.8
4.5% (95% CI 2.8-7.0)  ->  5.4% (95% CI 3.0-9.4)
13.5-24.2        ->  13.8-24.1
3.2-10.0         ->  3.5-10.0
2.6-6.0          ->  2.7-6.0
Figure 4a        ->  Figure 4A
which is significantly if frequency changes reflected noise alone (binomial test, p < 2x10-16)
                 ->  more than would be expected if frequency changes reflected noise alone (binomial test, p = 9.3 × 10⁻²⁰)
```

Search `2.6-6.0` rather than `4.0% (95% CI 2.6-6.0)`; the longer string may not
match if your document spaces the parenthesis differently from the PDF export.
The sentence 4 fragment leaves `79% of wild isolates` alone — change that to
`79% of wild isotypes (183 of 231)` separately, since the surrounding wording
may differ from the PDF export.
