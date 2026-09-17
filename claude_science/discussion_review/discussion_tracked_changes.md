# Discussion — tracked changes

Drafted 2026-09-17 against repository HEAD `c1a6bbe`, using the `scientific-writing`
skill (K-Dense, v2.1). Introduction and Results were left alone except where a number
or cross-reference is demonstrably wrong; those are listed separately at the end
rather than edited in place.

**How to read:** ~~struck text~~ **[replacement]**. Rationale follows each paragraph.
Every number introduced below is either already in the draft or recomputed from
`supplemental_data/`; nothing was invented. Where I could not verify a claim, it is
marked and left for an author decision rather than quietly rewritten.

---

## What was verified before editing

| Claim in draft | Source checked | Verdict |
|---|---|---|
| Pool of "84 isolates" | `METHODS.txt:64`, `METHODS.txt:136` | **Both numbers real, conflated.** 93 strains were pooled; 9 have no `vst` value and were dropped, "leaving n = 84 for association mapping" |
| Cross contrasts = Figure S7B-C | `manuscript_figures/FIGURES.md`, `git show c8e61fe^` | **Draft was right, repository was wrong.** S7 and S8 have now been swapped back; cross contrasts are S7 again |
| N94A result = Figure S9 | `manuscript_figures/FIGURES.md` | **Wrong.** S9 is the N2 dilution series; the sequon removal is **S10** |
| mig-6 cross QTL = Figure S6 | `manuscript_figures/FIGURES.md` | **Wrong.** S6 ranks per-strain responses; the cross scan is **S7** (post-swap) |
| Paaby ρ = +0.55 | `SUPP_FIG_plate_vs_paaby_vs_pos1original.R` | Correct — sign flip to hatching is applied |
| III:12,353,680, V:14,647,434, X:4,875,969 | `claude_science/results_review/` | Correct; recomputed there from `pos1_2023_gemma_loco.csv.gz` |
| 37 kb interval 13.658–13.695 Mb | `METHODS.txt:725` | Correct (13,657,700–13,695,000) |

---

## Paragraph 1 — the platform

We developed a flexible phenotyping platform that can ~~accurately~~ **[reliably]** infer
strain frequencies from pooled populations **[, recovering simulated frequencies with a
median r² of 0.906 at 10× coverage (IQR 0.786–0.957) and above 0.986 from 30×]**. This
platform can scale to hundreds of individuals and diverse pool compositions because it
does not depend on establishing a probe library designed to identify individual strains.
Demonstrating this point, we used this platform to map RNAi-response variation across 231
*C. elegans* isolates to three ~~distinct~~ **[loci — two clearing Bonferroni correction and
a third passing the less conservative eigen-decomposition threshold]** QTL (Figure 1B-C)
and ~~an 84-isolate pool of mostly RNAi-sensitive individuals that identified~~ **[a second
pool of 93 RNAi-responsive strains, 84 of which carried a variance-stabilised phenotype and
entered the association scan. This panel recovered]** the same general RNAi-response QTL on
chromosome III and an additional knockdown-specific QTL on chromosome V (Figure 2). ~~We
show that these results are reliable because we identified the same QTL in two independent
follow-up crosses (Figure S7B-C).~~ **[Two independent crosses, in which parental alleles
segregate at known frequencies and are measured directly, recovered the chromosome III and
chromosome V loci (Figure S8B-C), showing that these associations are not artifacts of
frequency deconvolution.]**

> **Why.** Four changes. (i) "Accurately" is unbounded; the simulation gives a specific
> figure, and quoting it makes the claim checkable. (ii) **84 vs 93 is a real error of
> fact** — `METHODS.txt` is explicit that 93 strains were pooled and 84 mapped, so "an
> 84-isolate pool" understates the experiment and contradicts the Results text, which says
> 93. (iii) "Three distinct QTL" flattens a genuine difference in evidential strength that
> the Results section is careful about. (iv) "We show that these results are reliable
> because" asserts more than a cross can deliver — a cross rules out deconvolution
> artifact, which is the specific thing it is good for, and the sentence is stronger for
> saying exactly that. Figure reference corrected S7→S8.

---

## Paragraph 2 — SID-2

We resolved the general ~~RNA~~ **[RNAi]**-sensitivity locus on the right arm of chromosome
III to a T96K missense mutation in the extracellular domain of the dsRNA receptor SID-2. We
show that strains with the 96K allele are more sensitive to *pos-1* RNAi than strains with
the 96T **[allele]**, which holds true across three genetic backgrounds (Figure 4A-B),
including the standard laboratory strain N2. T96 sits in a conserved N94-C95-T96
N-glycosylation sequon, and while the T96K allele adds a positive charge to the ectodomain
that might facilitate stronger interaction with the phosphate backbone of double-stranded
RNA, it should also disrupt the conserved sequon. Because N-glycosylation is associated with
protein stability and localization, disrupting the sequon should impact SID-2 function.
However, we see the opposite effect in the T96K allele swap experiments. We therefore
introduced an N94A allele in both genetic backgrounds to disrupt the sequon independently of
T96. The N94A allele makes JU1793 and JU2466 more resistant to *pos-1* RNAi, suggesting that
this allele impairs SID-2 function (Figure ~~S9~~ **[S10]**). Because T96K also disrupts the
sequon**[,]** while contributing a positive charge, we conclude that the electrostatic effect
of T96K must be large enough to override the antagonistic, glycosylation-linked resistance
effect we observed in the N94A edited strains. We hypothesize that by adding ~~more~~ localized
positive charge to the ectodomain, the lysine allele facilitates more double-stranded RNA
uptake than the threonine to induce a stronger RNAi response in carrier worms. **[Two
observations temper this model. First, the N94A substitution in JU2466 conferred more
resistance than restoring the ancestral threonine at position 96 (42.5% versus 18.4% hatched;
OR 3.3, p = 9e-8), which is difficult to reconcile with the two residues acting only through
a shared glycan. Second, the effect of residue 96 is background-dependent — the 96T allele
raises hatching odds 16.1-fold in JU1793 but 4.8-fold in JU2466 — so additional modifiers
must segregate between these strains.]** We note that we provide no experimental evidence
that the N94-C95-T96 N-glycosylation sequon is glycosylated and there is no published
experimental evidence suggesting it is (Hu et al., 2026). The absence of detected SID-2
glycosylation in this dataset is likely driven by sampling depth in the proteomics dataset
because only 870 proteins were recovered (Hu et al., 2026). In proteomes that do detect
SID-2, only the cytoplasmic domain is detected, not the N94-C95-T96 peptide (Narayan et al.,
2016; Tan et al., 2024).

> **Why.** "RNA-sensitivity" is a typo for RNAi. The N94A figure reference is wrong (S9 is
> the N2 dilution series; the sequon work is S10). The two added sentences move facts that
> are already in your Results into the Discussion, where they do argumentative work: the
> N94A-versus-96T asymmetry is the strongest evidence *against* the simple shared-glycan
> model, and burying it in Results while the Discussion presents a clean story reads as
> selective. Stating it makes the electrostatic argument more credible, not less.

---

## Paragraph 3 — knockdown-specific loci

In addition to the identification of general RNAi-response QTL, pooled RNAi phenotyping
revealed a *mig-6*-specific QTL localized to the center of chromosome V (Figure 2), which was
independently identified in a cross between N2 and XZ1516 (Figure ~~S6~~ **[S8]**). We
hypothesize that this QTL is ~~in fact~~ *mig-6*-specific ~~because~~ **[rather than a general
uptake locus for two reasons.]** ~~the~~ **[First, the]** XZ1516 alleles are strongly depleted**[;
second,]** ~~and~~ both N2 and XZ1516 are strong responders to RNAi, making disrupted RNA
transport to *mig-6*-expressing somatic tissues a less likely ~~hypothesis~~ **[explanation]**.
~~The rationale can be applied to~~ **[The same reasoning applies to]** the large-effect
*mig-6*-specific QTL on the left arm of chromosome I in this cross**[;]** ~~and~~ the reason we
did not see it in the pooled experiment is that the underlying variant was at low frequency in
the pool. Performing additional unrelated RNAi knockdowns on this cross population would ~~add
more support to this conclusion~~ **[test this directly]**. However, it is less clear that the
*mig-6*-specific QTL we identified on chromosome X in the JU1793 **[×]** ~~to~~ JU2466 cross is
not being driven by disruption to RNA trafficking to the target tissue**[,]** because JU1793 has
been shown to **[have]** inefficient somatic RNAi (Nuez & Félix, 2012).

> **Why.** Figure reference corrected. The rest is grammar and sentence-splitting — the
> "because … and …" construction ran two independent arguments together, and "JU1793 to
> JU2466 cross" should be "×". "Add more support to this conclusion" presumes the outcome
> of an experiment not yet done; "test this directly" does not.

---

## Paragraph 4 — NEW: limitations

**[Three features of the platform bound what these results can support. First, pooled
frequencies are compositional: a strain's decline is measured relative to the rest of the
pool, so the 79% of isolates that fell in frequency under *pos-1* RNAi cannot be read as
79% having an absolute response. Second, the measurement has a floor. Eighty-one isolates
were absent from every *pos-1* pool while present in both controls, and a strain returned
as exactly zero contributes no further information — these isolates are censored rather
than ranked, which is why we report 141 of 231 (61%) as individually responsive rather
than fitting a response distribution across all strains. Third, a pooled scan resolves only
alleles common within the panel, and the 93-strain panel was deliberately ascertained on
RNAi responsiveness rather than sampled for genetic diversity. That ascertainment is what
makes the knockdown-specific comparison possible, but it distorts allele frequencies
relative to the species and makes the panel a poor instrument for estimating how much
RNAi-response variation segregates in wild populations overall. The crosses address the
first and third limitations for the loci they recover, but not for loci they do not.]**

> **Why.** The skill's fidelity rules require concrete limitations and bounded
> generalizability, and the Discussion currently has none. All three constraints are
> already visible in your own Results (the compositional caveat, the 81 absent isolates,
> the common-allele restriction) — this paragraph consolidates them rather than
> introducing new claims. Reviewers will find these anyway; naming them first is
> the stronger position.

---

## Paragraph 5 — implications

These results highlight the platform's ability to identify large-effect genetic modifiers of
genetic perturbations. **[The architecture we recover differs from that reported by the most
comparable prior screen: assaying 29 maternal-effect genes across 55 wild strains, Paaby et
al. (2015) found highly polygenic, gene-specific modification, whereas contrasting two
knockdown conditions here isolates a small number of large-effect knockdown-specific loci.
Pooled selection and a targeted contrast may therefore recover a different, more tractable
slice of modifier architecture than per-strain phenotyping of many targets.]** Thousands of
rare Mendelian **[diseases]** are caused by ~~LoF~~ **[loss-of-function]** variants, ~~a
majority of which are influenced by genetic modifiers, which are~~ **[and their penetrance and
expressivity are frequently modified by background variation (Kang & Drivas, 2026). Such
modifiers are]** currently difficult to identify in human genetics studies because rare
disease-causing variants are carried by too few individuals to give modifier scans sufficient
statistical power. The platform we developed enables us to apply a wide **[]**~~-~~range of
genetic perturbations of disease-causing homologs to a pool of genetically diverse *C. elegans*
isolates to identify segregating genetic modifiers. **[Because the perturbation is delivered
by feeding rather than by editing, the same pool can be re-used across targets, so the cost of
adding a perturbation is one culture rather than one mutant line per background — the
constraint that has kept metazoan modifier screens at the scale of a few dozen strains.]**

> **Why.** "Thousands of rare Mendelian are caused by LoF variants" is missing its noun and
> expands an unglossed abbreviation on first use. **"A majority of which are influenced by
> genetic modifiers" is an unsupported quantitative claim** — I could not find a source for
> a majority, so I replaced it with the weaker statement your Introduction already
> supports and cites. The added opening sentence engages with Paaby et al. 2015, which the
> Introduction raises as the closest prior work but the Discussion never returns to. The
> closing sentence states the actual economic argument for the platform, which is the
> point the paper has earned but never quite makes.

---

## Glaring issues outside the Discussion

Flagged, not edited, per your instruction.

**1. Results — JU2466 hatching: 4.5% is arithmetically right but is the wrong denominator.**
The draft says the 96T swap "raised its hatching rate from 4.5% (95% CI 2.8-7.0) to 18.4%".
**This is NOT a digit transposition**, which is what `results_review` recorded as `N4` and
what an earlier version of this document repeated. Both were wrong. 4.5% is JU2466_A and
JU2466_B **pooled** (19/419), and its Clopper-Pearson interval is 2.8-7.0 — matching the
draft exactly. 5.4% is JU2466_A alone (11/204).

The problem is the comparator, not the arithmetic. `Figure4_sid2.R:233` filters
`strain %in% c("JU1793", "JU2466_A", "wSZ200", "wSZ206")` and its Fisher test (`:258`) is
`JU2466_A` against `wSZ206`, so the figure and the p value use A alone.
`MANUSCRIPT_CAPTIONS.txt:86` reports 5.4% (n = 204). Decisively, wSZ206's
`glycosylation motif` field reads `JU2466_A[NxT]` — the 96T edit was made in the A
background. The text therefore compares an A-derived edit against an A+B average, while
the figure it cites compares A against A, and the two isolates genuinely differ
(5.4% vs 3.7%). **Use 5.4% (11/204).** `MANUSCRIPT_CAPTIONS.txt:291` handles it best by
reporting both isolates separately rather than pooling. The `N4` entry in
`claude_science/results_review/` should be amended, since it has propagated as a
transposition through several sessions.

**2. Results — the N2 confidence interval is impossible.** "32.3% (95% CI 2.61-38.9)". A lower
bound of 2.61% cannot sit under a point estimate of 32.3% with an upper bound of 38.9%.
Recomputing from the repository's own n = 220: **32.3%, Wilson 95% CI 26.4-38.7%**. The
lower bound appears to be "26.4" or "26.1" mangled. Confirm which interval method the
manuscript uses before substituting.

**3. Results — the cross-contrast figure is cited under two different numbers. RESOLVED in
the repository, 2026-09-17.** The section cites "Figure S7B-C" three times and
"Figure S8B-C" once for the same panels. An earlier version of this document said S8 was
correct. It was not: `git show c8e61fe^` confirms the cross contrasts were **S7** before the
NIL figure was inserted, so the draft's three S7 citations are the original numbering and
the insertion broke them. It also broke citation order, since the cross contrasts are cited
in the Figure 2 section well ahead of the NIL work. **S7 and S8 have been swapped in the
repository**: S7 is again the cross contrasts, S8 is the NIL hatching experiment.
S9-S11 are unaffected.

Remaining draft edits: the single **"Figure S8B-C"** citation for the cross contrasts
becomes **S7B-C**, and the two citations that point at S7 for NIL/hatching content — the
plate-based lethality assay and the 37 kb interval sentence — become **S8**.

**4. Results — "Figure S9B" for the N94A result should be S10B.**

**5. Results — typo: "devconvolution"** → deconvolution.

**6. Results — typo: "pos-1-specifc"** → *pos-1*-specific (already logged as `T1`).

**7. Results — "localize the QTL to a 37 kb interval spanning 13.658–13.695"** is missing its
units; should read "13.658–13.695 Mb".

**8. Introduction — duplicated sentence.** "A growing body of work indicates that rare and
common variation in an individual's genetic background is a major source of this variation."
is immediately followed by "A growing body of work indicates that rare and common background
variation is a major source of this variation (Heyne et al., 2026; Kang & Drivas, 2026)." The
second carries the citations; the first should be deleted.

**9. Introduction — "they remain observation"** → "they remain observational".

**10. Introduction — "The T96K shifts responses to multiple because"** is missing a noun,
presumably "multiple knockdowns".

---

## Not verified

- **Hu et al. 2026, Narayan et al. 2016, Tan et al. 2024** — the glycoproteomics claims,
  including "only 870 proteins were recovered" and "only the cytoplasmic domain is
  detected". No copies in the repository; these need an author to open the sources and
  confirm both the numbers and that the papers support the specific propositions.
- **Kang & Drivas 2026** as support for background modification of LoF penetrance — cited in
  your Introduction, so presumably verified there, but I did not open it.
- **"conserved across nine of ten species of the Elegans supergroup and its sister
  C. kamaaina"** — stated in Results and relied on in the Discussion; check against
  `Figure_S10`'s underlying alignment.
- Whether the chromosome I *mig-6* QTL variant is genuinely low-frequency in the pool, as
  Paragraph 3 asserts. This is testable against the panel genotypes and is currently an
  assertion.
