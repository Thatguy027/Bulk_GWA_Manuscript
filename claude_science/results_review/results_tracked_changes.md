# Results section, tracked changes

Struck text is the current draft; **[bracketed]** text is the proposed
replacement. Every change is keyed to a claim id in `claims_verified.tsv`,
which carries the recomputed value and the file it came from.

5 corrections of fact, 4 wording calls for the author, 8 typos.

---

Non-negative least squares regression enables accurate strain frequency inference
To determine if non-negative least squares (NNLS) regression can accurately infer strain frequencies from pooled populations, we established a simulation framework (Methods). Briefly, each wild isolate was assigned a fitness value drawn from an inverse-χ² distribution, with the fitness structure taken from seven published C. elegans traits with validated QTL. For each trait, we computed the expected pooled allele frequencies such a population would produce and simulated observed alt-allele counts via binomial sampling across a range of sequencing depths (1–500×). We then deconvolved these simulated counts back to per-strain frequencies using NNLS regression and compared the estimates with the known input as a function of depth. ~~This simulation framework revealed that NNLS regression can accurately infer strain frequencies with as little as 1X sequencing depth (Figure S1).~~ **[This simulation framework revealed that NNLS regression recovers strain frequencies with a mean r² of 0.76 at 1X sequencing depth (range 0.51-0.91 across the seven traits), rising above 0.95 from 10X (Figure S1).]**  

Next, we wanted to establish an effective strain pooling and frequency inference strategy for genetically diverse C. elegans isolates. To this end, we performed a small-scale experiment where we extracted DNA from ~~two pooled populations composed of ~48 strains each~~ **[two pooled populations composed of 46 and 40 strains]**. We combined the DNA extracted from these populations at known ratios, sequenced the resulting libraries, and inferred individual strain frequencies. Because we did not pool the strains at known frequencies, we summed the strain frequencies derived from each pool to establish an input pooled population ratio. We found that the inferred shares tracked the designed ratios closely (Pearson r = 0.997, root mean squared error 0.038 in fraction units), with the largest deviations occurring at the smallest-volume step of the dilution (Figure S2).

As a final test for the NNLS approach, we compared it with a previously published pooled-phenotyping experiment that used molecular inversion probe sequencing to infer strain frequencies (Webster et al. 2022). To make the comparison, we performed whole-genome sequencing on the same samples. For each strain, we took the rate of change in pool frequency across days of L1 starvation as the slope of frequency change regressed on day, using day 1 as the baseline, and averaged those slopes over the five replicate arms (Methods). We found the NNLS-derived slopes to be highly correlated with those derived from MIP-seq ~~(Spearman’s ρ = 0.97, n = 98 strains, p < 1e-4)~~ **[(Spearman’s ρ = 0.974, n = 98 strains)]** (Figure 1A). Agreement was largely retained as we downsampled reads, reaching ρ = 0.85 at 1x (Figure S3). Taken together, these three validation experiments gave us confidence that we could infer per-strain frequency dynamics accurately enough to treat them as a quantitative phenotype.
Pooled phenotyping enables evaluation of wild isolate RNAi responses
C. elegans RNAi screens are routinely performed in the N2 strain. However, extending these screens to multiple wild isolates is challenging because each isolate needs to get assayed independently. We reasoned that the pooled phenotyping approach would facilitate fast evaluation of wild isolate RNAi responses. N2 strains exposed to pos-1 RNAi exhibit an embryonic lethal phenotype, making this treatment amenable to the pooled phenotyping because RNAi-responsive strains will fail to contribute progeny to subsequent generations. To explore if we could quickly evaluate wild isolate RNAi responses, we performed a pilot experiment where we exposed 231 pooled wild isolates to pos-1 RNAi-expressing HT115 bacteria or control HT115 bacteria (Methods). ~~We grew these populations for two generations in these conditions across two replicates,~~ **[We grew these populations for two generations in these conditions across four pos-1 replicate pools and two control pools,]** sequenced the resulting F3 L1 populations, and inferred the individual strain frequencies in each population. We calculated individual strain frequency differences between the control and pos-1 RNAi conditions and found good agreement between the replicate conditions (Figure S4; Spearman’s ρ = 0.77-0.87 across all six pairwise replicate comparisons). As expected from previous reports that have identified substantial RNAi-response variation across wild C. elegans isolates, we observed that 183 of 231 (79%) of the wild isolates were responsive to RNAi, as indicated by these strains having a lower frequency in pos-1 RNAi conditions than they did in the control condition (Figure 1B). We performed a genome-wide association scan on the variance-stabilized pos-1 responses, which identified QTL on the right arm of chromosome IV and ~~the center of chromosome X that passed the bonferroni significance threshold~~ **[the left arm of chromosome X (X:4,875,969) that passed the Bonferroni significance threshold]** ~~and a third QTL on the right arm of chromosome III above the eigen-decomposition threshold (Figure 1C).~~ **[and a third QTL on the right arm of chromosome III above the eigen-decomposition threshold (Figure 1C). A tenth Bonferroni-passing marker, on chromosome III at 5.97 Mb, has no supporting marker within 100 kb of it and we do not carry it forward.]** From this analysis alone, ~~it is unclear that these QTL are pos-1-specifc QTL~~ **[it is unclear whether these QTL are pos-1-specific QTL]** or general RNAi-response QTL. We reasoned that if they are general RNAi-response QTL, we would want to minimize their contribution to future experiments. Therefore, we decided to construct a pooled population of RNAi-responsive strains. To construct this population, we manually re-evaluated pos-1 RNAi responses of 191 wild strains on agar plates and identified 93 strains with robust RNAi responses that we pooled to evaluate additional RNAi responses (Methods). We did not include N2 in this panel because it is the reference strain and by definition does not have alternate genotype calls that can be used for strain inference. We compared the manual pos-1 RNAi phenotypes we collected to the pooled RNAi responses based on strain frequency changes above and found good agreement between the methods (Spearman’s rho = 0.41, p = 7.8×10⁻⁶; n =111, Figure S5A). We also found good agreement between our manual phenotypes and a previously published evaluation of wild isolate RNAi responses (Spearman’s rho = −0.55, p = 0.014; n=19;  Figure S5B) (Paaby et al. 2015).
A large-effect RNAi sensitivity QTL maps to the right arm of chromosome III
Having assembled a pool of RNAi-responsive strains, we next asked whether we could identify QTL that modify responses to individual RNAi treatments. We grew this pool on bacteria expressing double-stranded RNA ~~against each of nine target genes~~ **[against each of ten target genes]** and on the HT115 empty-vector control, sequenced the resulting populations, and inferred strain frequencies (Methods). We performed genome-wide association scans on the variance-stabilized responses to these knockdowns and identified QTL that exceeded the Bonferroni threshold for pos-1 and mig-6 RNAi (Figure 2). The pos-1 response mapped to the right arm of chromosome III (III:12,353,680, -log10p = 7.04, allele frequency 0.11), just 0.4 Mb from the association we identified in our pilot experiment (Figure 1C). The mig-6 response mapped to the center of chromosome V (V:14,647,434, -log10p = 7.71, allele frequency 0.13), a locus we did not identify in either pos-1 experiment, suggesting this locus is mig-6-specific. 

To attempt to reproduce these QTL we exposed two crosses to mig-6 and pos-1 RNAi after intercrossing for ten generations (Methods). One cross was constructed between JU1793, which was resistant to both mig-6 and pos-1 RNAi and JU2466, which was sensitive to pos-1 RNAi and had an average mig-6 response (Figure S6). We made use of a previously constructed cross between N2 and XZ1516 because XZ1516 had among the strongest responses to both RNAi conditions (Figure S6) (Zdraljevic et al. 2025). We sequenced the cross populations after two generations on the RNAi conditions and used allele frequency deviations between conditions to perform QTL mapping. By comparing parental frequencies after RNAi exposure to the control HT115 condition, we identified multiple QTL shared between the crosses and cross-dependent QTL (Figure S7B-C, Supplemental table X). ~~Additionally, we identified  QTL that are shared between conditions~~ **[Additionally, we identified QTL that are shared between conditions]**, which we refer to as general RNAi-response QTL. By contrasting allele frequencies after pos-1 and mig-6 RNAi exposure, ~~we identified RNA-specific QTL~~ **[we identified RNAi-specific QTL]**. One major general RNAi-response QTL we identified in both crosses after mig-6 and pos-1 exposure localized to the right arm of chromosome III (Figure S7B-C), in close proximity to the QTL we observed from the pooled association mapping (Figure 2). We also identified a large mig-6-specific QTL in the XZ1516 cross that co-localized with the mig-6-specific QTL we identified in the pooled experiment. In addition to these re-discovered QTL, we identified two novel large-effect mig-6-specific QTL on chromosome I in the XZ1516 cross (LOD = 751.1, β = 0.391) and on chromosome X in the JU1793 cross (LOD = 361.4, β = 0.713). Taken together, these results highlight the power of the pooled-phenotyping approach to accurately evaluate RNAi responses across the C. elegans population, point to substantial RNAi-response variation segregating in the population, and highlight several loci that modify gene-specific knockdown responses. 
A high-frequency variant in the dsRNA transporter SID-2 confers increased sensitivity to ingested RNAi

We next sought to uncover the variant underlying the general RNAi-response QTL on the right arm of chromosome III. We used JU1793 and JU2466 as parents to construct near-isogenic lines. The pooled assay showed that JU1793 is resistant to pos-1 RNAi (Figure 3A), which was confirmed by allele-frequency distortions in the cross experiment described above (Figure 3B). To track the NIL phenotypes for resolving this locus, we turned to plate-based embryonic lethality assays (Methods). Through iterative NIL construction we were able to ~~localize the QTL to a 37 kb interval spanning 13.658 - 13.695.~~ **[localize the QTL to a 37 kb interval spanning 13.658-13.695 Mb.]** There are 27 variants that distinguish the two parental strains within this region, including three synonymous, 22 intergenic, intronic, or in UTRs, and two missense variants located in the sid-2 gene - V5L and T96K. sid-2 encodes ~~a single-pass transmembrane protein localized apical membrane of intestinal cells~~ **[a single-pass transmembrane protein localized to the apical membrane of intestinal cells]** that facilitates uptake of dsRNA from the gut lumen, making it an excellent candidate to harbor variation that contributes to RNAi sensitivity across wild isolates. We focused on the T96K variant because it is predicted to reside in the extracellular domain of SID-2 and therefore might alter interactions with double-stranded RNA. We generated reciprocal allele-swap strains to test whether this variant was contributing to pos-1 RNAi responses. Introducing the 96K allele into the JU1793 strain ~~lowered the strains hatching rate~~ **[lowered the strain’s hatching rate]** from 94.8% to 53.1%, while introducing the 96T allele 
into the JU2466 strain ~~raised its hatching rate from 4.5% to 18.4%~~ **[raised its hatching rate from 5.4% to 18.4%]** (Figure 4a). Therefore, the 96T reference allele produces a weaker RNAi response than the alternate 96K allele.  These results suggest that this allele interacts with other variants in the parental strain’s background to produce background-dependent effects. 

Given that N2 is the standard laboratory reference strain that has been used in countless RNAi studies to great success, we were surprised to find that the 96T reference allele caused a weaker RNAi response in the JU2466 background. We therefore asked if we could increase the RNAi sensitivity in N2 by introducing the 96K allele. Introducing this allele into the N2 strain produced no difference in the embryo hatching rate in the same experimental condition we used for all plate-based embryo hatching experiments (a 50:50 mixture of pos-1 and HT115 bacteria) (Figure S8). However, when we lowered the pos-1 bacteria concentration on the assay plates, ~~we found that the 96K allele lowered the N2 hatching rate from 32.3% to 4% (Figure 4B)~~ **[we found that the 96K allele lowered the N2 hatching rate from 32.3% to 4.0% (two independently edited lines, 4.4% and 3.7%; Figure 4B)]**, therefore increasing RNAi sensitivity in this background as well. 

We hypothesized that by introducing a positively charged residue the T96K allele might facilitate a stronger interaction with the phosphate backbone of double-stranded RNA and therefore better uptake of the double-stranded RNA. To add support to this hypothesis ~~we calculated the per-reside local net charge~~ **[we calculated the per-residue local net charge]** in 12 Å windows across the SID-2 ectodomain and found that the T96K allele sits in a positively charged region of the ectodomain (Figure 4C). The local 12 Å region of the 96T is in the 82nd percentile of all residues with a net charge of +1.24 e, largely driven by K93 and K132 sitting at 6.6 and 6.8 Å away, respectively (Figure S9). Introducing the 96K variant to this local region shifts the net charge to +2.24 e, which corresponds to the 98th percentile of the ectodomain's local net charge distribution (Figure S9). While these observations don’t prove that the T96K allele is causing a stronger interaction with double-stranded RNA by ~~increasing the net charge of it’s local environment~~ **[increasing the net charge of its local environment]**, it does agree with previous reports that modifying critical histidines in the ectodomain to positive arginines is less disruptive to uptake of double-stranded RNA than modifying the histidines to alanines (McEwan et al. 2012).


---

## Why each change

**S3 - ERROR**

- current: This simulation framework revealed that NNLS regression can accurately infer strain frequencies with as little as 1X sequencing depth (Figure S1).
- proposed: This simulation framework revealed that NNLS regression recovers strain frequencies with a mean r² of 0.76 at 1X sequencing depth (range 0.51-0.91 across the seven traits), rising above 0.95 from 10X (Figure S1).
- reason: simulation_seeded_r2.tsv: mean r2 per trait at 1x is 0.51-0.91, median 0.76; the mean first reaches 0.95 at 10x. SUPP_FIG_XX_simulation_depth.R prints this under the heading 'accuracy at 1x, the depth the text claims'.

**D1 - ERROR**

- current: two pooled populations composed of ~48 strains each
- proposed: two pooled populations composed of 46 and 40 strains
- reason: dilution_strain_sets.tsv: 174 isolates in four sets, A=46, B=46, C=40, D=42. The titration is set B against set C, so the two populations are 46 and 40, not ~48 each.

**P2a - ERROR**

- current: We grew these populations for two generations in these conditions across two replicates,
- proposed: We grew these populations for two generations in these conditions across four pos-1 replicate pools and two control pools,
- reason: pos1_2023_sample_frequencies.csv.gz: pos-1 T2 has replicates 1,2,3,4 and ctrl T2 has A,B. Four pos-1 pools is also what makes the six pairwise comparisons cited in the next sentence possible - two replicates give only one pair.

**C1 - ERROR**

- current: against each of nine target genes
- proposed: against each of ten target genes
- reason: pooled_vst_traits.csv.gz and bundle $pheno both carry ten RNAi targets at 93 strains each: fog-2, mig-6, pos-1, rde-3, ric-3, rpn-12, spe-19, spe-43, unc-39, vha-5.

**N4 - ERROR**

- current: raised its hatching rate from 4.5% to 18.4%
- proposed: raised its hatching rate from 5.4% to 18.4%
- reason: ju_allele_swaps_hatching.csv, pos-1 condition: JU2466_A[96K] is 11 hatched of 204 = 5.4%. Figure4_sid2.R prints 0.054. The 18.4% for wSZ206 is correct. Looks like a digit transposition.

**P4b - JUDGE**

- current: the center of chromosome X that passed the bonferroni significance threshold
- proposed: the left arm of chromosome X (X:4,875,969) that passed the Bonferroni significance threshold
- reason: The single chromosome X marker clearing Bonferroni is at 4,875,969, which is 28% along a 17.72 Mb chromosome. Also fixes the lower-case 'bonferroni'.

**P4d - JUDGE**

- current: and a third QTL on the right arm of chromosome III above the eigen-decomposition threshold (Figure 1C).
- proposed: and a third QTL on the right arm of chromosome III above the eigen-decomposition threshold (Figure 1C). A tenth Bonferroni-passing marker, on chromosome III at 5.97 Mb, has no supporting marker within 100 kb of it and we do not carry it forward.
- reason: FIGURE_REPORT.md 'GWAS interval admission' settles this: the 5.966 Mb marker clears Bonferroni at 8.68 with zero eigen-passing neighbours in 100 kb, while the 12.70-12.80 Mb cluster peaks below Bonferroni at 6.31 with 14. The draft's framing follows the repository's admission rule, but a reader comparing against Figure 1C will see a red marker on chromosome III that the text does not explain.

**N6 - JUDGE**

- current: we found that the 96K allele lowered the N2 hatching rate from 32.3% to 4% (Figure 4B)
- proposed: we found that the 96K allele lowered the N2 hatching rate from 32.3% to 4.0% (two independently edited lines, 4.4% and 3.7%; Figure 4B)
- reason: n2_allele_swaps_hatching.tsv at the 25% dose: N2[96T] 32.3% (n=220), and the pooled 96K value 4.0% is two independent lines wSZ203 4.4% (n=273) and wSZ204 3.7% (n=295); 273+295=568, the pooled n. The pooled figure is right; naming the two lines matches the caption and is the stronger claim.

**M1 - JUDGE**

- current: (Spearman’s ρ = 0.97, n = 98 strains, p < 1e-4)
- proposed: (Spearman’s ρ = 0.974, n = 98 strains)
- reason: Recomputed rho is 0.974 over 98 strains. Figure 1A deliberately reports a bootstrap interval rather than a p value - Figure1_common.R notes that a p value against rho = 0 is not the question. Quoting p < 1e-4 in the text reintroduces what the figure dropped.

**T1 - STYLE**

- current: it is unclear that these QTL are pos-1-specifc QTL
- proposed: it is unclear whether these QTL are pos-1-specific QTL
- reason: typo: specifc -> specific; 'unclear that' -> 'unclear whether'.

**T2 - STYLE**

- current: we identified RNA-specific QTL
- proposed: we identified RNAi-specific QTL
- reason: RNA-specific should be RNAi-specific; the contrast is between RNAi treatments.

**T3 - STYLE**

- current: lowered the strains hatching rate
- proposed: lowered the strain’s hatching rate
- reason: possessive apostrophe.

**T4 - STYLE**

- current: increasing the net charge of it’s local environment
- proposed: increasing the net charge of its local environment
- reason: it’s -> its.

**T5 - STYLE**

- current: a single-pass transmembrane protein localized apical membrane of intestinal cells
- proposed: a single-pass transmembrane protein localized to the apical membrane of intestinal cells
- reason: missing preposition.

**T6 - STYLE**

- current: localize the QTL to a 37 kb interval spanning 13.658 - 13.695.
- proposed: localize the QTL to a 37 kb interval spanning 13.658-13.695 Mb.
- reason: missing units. The wSZ191 introgression in nil_introgression_ranges.bed is 13,657,700-13,695,000 = 37.3 kb, so both the width and the bounds are right.

**T7 - STYLE**

- current: Additionally, we identified  QTL that are shared between conditions
- proposed: Additionally, we identified QTL that are shared between conditions
- reason: double space.

**T8 - STYLE**

- current: we calculated the per-reside local net charge
- proposed: we calculated the per-residue local net charge
- reason: typo: reside -> residue.
