# Review notes for the pooled-pool / cross paragraphs

Sources: `FIGURE_CAPTIONS.txt` (FIGURE 2, SUPP_FIG_XX_pooled_phenotype_ranks,
SUPP_FIG_XX_pooled_phenotype_heatmap, SUPP_FIG_XX_cross_contrast_panels,
SUPP_FIG_XX_cross_qtl_all), `METHODS.txt` (pooled RNAi competition; F2 bulk-segregant crosses),
and `supplemental_data/mapping/pooled_cross_bundle_thinned.rds`. Values recomputed from the
bundle unless noted.

## What the deposit supports

**The pooled experiment.** The pool was grown against ten targets (fog-2, mig-6, pos-1, rde-3,
ric-3, rpn-12, spe-19, spe-43, unc-39, vha-5) plus the HT115 empty-vector control, sequenced at
two timepoints. Nine of the 93 strains have no vst value, so n = 84 for mapping. Association
used the vst parameterisation, 322,010 markers, Bonferroni -log10 p = 6.81, eigen 4.17 over 732
effective tests.

**Exactly two markers clear Bonferroni across both scans, one per treatment:**

| treatment | peak marker | Mb | -log10 p | allele freq | beta |
|---|---|---|---|---|---|
| mig-6 | V:14647434 | 14.647 | 7.71 | 0.131 | +0.084 |
| pos-1 | III:12353680 | 12.354 | 7.04 | 0.107 | +0.046 |

Both betas positive: the minor allele is the resistant one. The next-best peaks are well below
threshold (mig-6 II:8.16 Mb at 5.54; pos-1 IV:16.26 Mb at 5.51), so each treatment gives one
locus and they are on different chromosomes.

**The chromosome III peak is the pilot's locus.** III:12,353,680 is 0.365 Mb from the pilot
scan's chromosome III eigen-level peak at III:12,718,465. Locus-level, the same signal.

**Of the four cross parents, only JU1793 carries the minor allele** -- at BOTH peaks (1/1 at
V:14647434 and III:12353680; N2, XZ1516 and JU2466 are all 0/0).

**The crosses.** N2 x XZ1516 and JU1793 x JU2466, parents taken from opposite extremes of the
pooled assay. 18 contrast scans, 108 genome-wide significant intervals at LOD 3.57, of which 46
have peak LOD > 100; Figure 2 draws the nine with peak LOD > 100 in the two contrasts it uses.
The relevant ones:

| contrast | cross | chrom | peak Mb | LOD | interval kb |
|---|---|---|---|---|---|
| mig-6 vs pos-1 | N2xXZ1516 | V | 13.81 | 730 | 278 |
| HT115 vs pos-1 | N2xXZ1516 | III | 13.31 | 710 | 394 |
| HT115 vs pos-1 | JU1793xJU2466 | III | 13.78 | 140 | 51 |

## Three things to get right in the text

### 0. LD makes the chromosome III correspondence defensible, not just close

`plots/diagnostics/TABLE_gwas_qtl_intervals_eigen.tsv` (built since these notes were started,
and covering the 2023 pilot scan) gives the pilot's chromosome III locus at III:12,718,465 an
r2 >= 0.5 interval of **11.436-13.784 Mb** (2,347 kb, 178 markers), and the interval is
unchanged at r2 >= 0.6 and r2 >= 0.7 (27 and 13 markers). That interval contains the pooled
pos-1 peak at III:12.354 Mb AND both cross chromosome III peaks, 13.31 Mb (N2 x XZ1516) and
13.78 Mb (JU1793 x JU2466). In other words the "they are megabases apart" objection is answered
by linkage disequilibrium on chromosome III: markers out to 13.78 Mb are in LD with the
association peak at r2 >= 0.7, so all three signals sit inside one LD block.

Two limits on that argument. At r2 >= 0.8 the interval collapses to 12.718-12.732 Mb (14 kb),
so the claim is cutoff-dependent and should name its cutoff. And the table covers the PILOT
scan only -- its single chromosome V locus is V:605248, not the pooled mig-6 peak at
V:14.647 Mb -- so the chromosome V correspondence (0.84 Mb from the cross peak) still rests on
proximity alone.

### 1. Co-localisation is locus-level, and the caption asks you to say so

`FIGURE_CAPTIONS.txt` marks this "CAVEAT to state in the text rather than leave a reader to
notice": no cross interval contains its matching pooled GWAS peak marker. Nearest
correspondences: mig-6 V:14.65 Mb against the N2 x XZ1516 V interval at 13.81 Mb (0.84 Mb), and
pos-1 III:12.35 Mb against the N2 x XZ1516 III interval at 13.31 Mb (0.96 Mb) and the
JU1793 x JU2466 interval at 13.78 Mb (1.43 Mb). The caption's own resolution is the sentence to
borrow: "A GWAS on 84 strains is not expected to localise to a cross interval; the claim is
concordance of locus, not of marker." One clause covers it and forecloses the obvious referee
question.

### 2. Chromosome III is the GENERAL locus, not the pos-1-specific one -- and the contrast design is why

CONFIRMED BY THE AUTHOR. The chromosome III locus is not *pos-1*-specific: it shows up under
multiple RNAi conditions, and the contrast design is what makes that legible.

**Read the contrast design before reading the QTL.** Each condition is contrasted against the
frequency changes under *pos-1*, so *pos-1* is the common reference arm rather than one
treatment among equals. A knockdown-SPECIFIC QTL is therefore one that appears in an
X-versus-*pos-1* contrast: the two arms differ at that locus because the locus acts on X and
not on *pos-1*. A locus affecting the RNAi machinery generally does the opposite -- it moves
the same way in both arms, cancels in the contrast between them, and shows up instead against
the HT115 empty-vector control. Chromosome III behaves the second way: it rises in
HT115-vs-*pos-1* and is flat between *mig-6* and *pos-1* (cross_LOD 0.0 at the GWAS peak
marker). Under the contrast design that is the signature of a general locus, not a *pos-1* one,
and the flatness is the evidence rather than an absence of it.

The repository already says so in two places. `SUPP_FIG_XX_cross_contrast_panels` exists to
"justify calling the chromosome III locus general and the chromosome V, I and X loci specific",
and Figure 3 fine-maps "the chromosome III RNAi-response QTL" rather than a *pos-1* QTL.

Calling it *pos-1*-specific would collide with both, and would also read the contrast
backwards -- treating a locus's absence from the *mig-6*-vs-*pos-1* contrast as evidence that
it is *pos-1*-specific, when under this design that absence is exactly what a general locus
produces. It also throws away the payoff set up in the previous paragraph: from the pilot alone
it was unclear whether these were *pos-1*-specific or general RNAi-response QTL, and the
contrast design is precisely what answers that. Suggested split -- the ASSOCIATION detected one
locus per treatment; the CROSSES then show chromosome III is general and chromosome V is
*mig-6*-specific.

### 3. The crossing scheme is a [TO FILL]

`METHODS.txt`, F2 BULK-SEGREGANT (xQTL) CROSSES: "[TO FILL: crossing scheme, F2 population
size, selection or competition regime, generations, DNA extraction, library preparation and
sequencing platform.]" Ten rounds of intercrossing followed by two rounds of selection on pos-1
and on mig-6 RNAi is not recorded anywhere in the repository, so it comes from you and needs to
land in Methods before the results text leans on it.

## Do not quote cross interval WIDTHS from the thinned bundle

`pooled_cross_bundle_thinned.rds` reproduces every peak position (to the kb) and every peak LOD
in the caption's table, but its LOD-drop interval bounds are inflated by thinning: chromosome
III HT115-vs-pos-1 comes out 547 kb against the caption's 394 kb, and the JU1793 x JU2466 III
interval 104 kb against 51 kb. Peaks and LODs from the bundle, widths from the caption table
(which was built from the full bundle in `data/`).

## Draft

> Having assembled a pool of RNAi-responsive strains, we next asked whether we could identify
> QTL that shape the response to individual RNAi treatments. We grew this pool on bacteria
> expressing double-stranded RNA against each of ten target genes and on the HT115 empty-vector
> control, sequenced the resulting populations, and inferred strain frequencies as above
> (Methods). We then performed genome-wide association scans on the variance-stabilized
> responses to two of these treatments, pos-1 and mig-6 RNAi (322,010 markers, n = 84 strains;
> Figure 2). Two markers exceeded the Bonferroni threshold, and each was significant for only
> one treatment: the mig-6 response mapped to V:14,647,434 (-log10 p = 7.71, allele frequency
> 0.13), a locus we did not detect in the pilot experiment, and the pos-1 response mapped to
> III:12,353,680 (-log10 p = 7.04, allele frequency 0.11), 0.4 Mb from the chromosome III
> association we identified in the pilot. At both loci the minor allele conferred resistance,
> and of the four wild strains we went on to use as cross parents, only JU1793 carried it.
>
> To test these loci independently and to ask whether they act generally or on specific
> knockdowns, we constructed two F2 bulk-segregant populations from strains at opposite extremes
> of the pooled assay, N2 x XZ1516 and JU1793 x JU2466. We intercrossed each population for ten
> rounds to break up parental haplotypes and then selected for two rounds on pos-1 and on mig-6
> RNAi (Methods), sequenced the selected populations, and contrasted allele frequencies between
> conditions. The strongest QTL fell on the same chromosomes as the pooled associations: the
> mig-6-versus-pos-1 contrast peaked at V:13.81 Mb (LOD 730) in the N2 x XZ1516 cross, and the
> HT115-versus-pos-1 contrast peaked at III:13.31 Mb (LOD 710) in the same cross and at
> III:13.78 Mb (LOD 140) in the JU1793 x JU2466 cross (Figure 2). No cross interval contains its
> matching association marker -- the nearest correspondences are 0.84 Mb apart on chromosome V
> and 0.96 Mb on chromosome III -- so the concordance is one of locus rather than of marker,
> which is what an association scan on 84 strains can support. On chromosome III that
> correspondence is supported by linkage disequilibrium: markers spanning 11.44 to 13.78 Mb
> remain in LD with the association peak at r2 >= 0.7, placing both cross peaks and the
> association peak within a single LD block. Comparing the two contrasts
> within each cross then separated the loci by mechanism. Because each condition is contrasted
> against the frequency changes under pos-1, a locus acting on one knockdown alone appears in
> that contrast, whereas a locus acting on the RNAi machinery generally moves both arms together
> and cancels. The chromosome III locus rose in the HT115-versus-pos-1 contrast but was flat
> between mig-6 and pos-1, marking it as a general RNAi-response locus, whereas the chromosome V
> locus appeared only in the mig-6 comparisons and is therefore specific to that knockdown
> (SUPP FIGURE XX).

If you want the specific/general verdict to arrive later instead -- with the NILs in Figure 3 --
cut the last sentence and keep the cross paragraph purely about concordance.
