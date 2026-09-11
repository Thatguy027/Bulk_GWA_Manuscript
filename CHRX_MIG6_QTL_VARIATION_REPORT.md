# The chromosome X *mig-6* QTL in JU1793 × JU2466: variation analysis

Analysis run 2026-09-11. **Standalone reference note, not a draft export.** Keep
it at the repository root: `check_manuscript_numbers.py` globs `manuscript/*.md`
and `manuscript/*.txt`, and these numbers are finer-grained than the ones
`FIGURE_REPORT.md` attests.

## The question

The chrX QTL in the *mig-6*-vs-*pos-1* contrast has an enormous effect and
enriches JU1793 alleles. JU1793 and JU1580 are the same isotype, and JU1580 has
been described as competent for germline RNAi but not for somatic RNAi. If that
is right, this QTL might be a **somatic dsRNA-trafficking locus** rather than a
***mig-6*-specific** one. This note asks whether the sequence in the region
supports that.

**Answer in one line:** it does not. No RNA-silencing gene in the enriched region
carries a coding difference between the parents, while a gene in the *mig-6*
glycosylation pathway sits 6.2 kb from the peak with a large parental expression
difference. The two hypotheses remain formally confounded by the experimental
design, and the discriminating experiment is a second somatic RNAi target.

---

## 1. The QTL

Peak **X:5,955,109**, from
`data/cross_experiments/JU1793-JU2466_export/plot_data/JU1793_JU2466_F2-2_contrast_MIG6g-POS1g_10000_plot_DF.tsv.gz`.

| contrast | LOD at X:5,955,109 | chrX peak LOD | chrX peak position |
| --- | --- | --- | --- |
| HT115g vs MIG6g | **687.4** | 798.0 | X:7,151,301 |
| HT115g vs POS1g | **5.8** | 9.4 | X:6,904,461 |
| MIG6g vs POS1g | **361.4** | 361.4 | X:5,955,109 |

JU1793 allele frequency at the peak marker:

| condition | JU1793 frequency |
| --- | --- |
| HT115 control | **0.091** |
| *mig-6* RNAi | **0.894** |
| *pos-1* RNAi | **0.181** |

The locus has a very large effect on the *mig-6* response and essentially none
on the *pos-1* response. **This rules out a general dsRNA-uptake locus,** which
would affect both conditions — the contrast with the chrIII/*sid-2* locus, whose
LOD collapses 94% in the RNAi-vs-RNAi contrast, is the relevant comparison.

**Two things it does not rule out,** and cannot: a *mig-6*-pathway modifier and a
somatic-RNAi-competence locus predict identical data here, because *mig-6* is the
only somatic target in the design and *pos-1* the only germline one.

**A confound worth stating.** The JU1793 allele sits at 0.091 in the untreated
control, so the haplotype is already at roughly 1:10 before any RNAi. Part of
this locus may be viability or fertility distortion rather than RNAi response,
and the "enrichment" is measured against that skew.

## 2. Interval definition

The shipped support interval uses a LOD drop of 10% of the chromosome peak
(36.14 LOD units). Because the effect is so large, that interval is narrow
relative to where signal actually extends, so the search was widened and
ultimately run across the whole chromosome.

| basis | interval | width |
| --- | --- | --- |
| 1.5-LOD drop | X:5,904,993–6,004,204 | 0.10 Mb |
| **shipped support interval (10% drop)** | **X:5,554,409–7,987,415** | **2.43 Mb** |
| 25% drop | X:5,277,492–9,409,907 | 4.13 Mb |
| 50% drop | X:4,930,672–13,036,529 | 8.11 Mb |
| JU1793 frequency ≥ 0.80 under *mig-6* | X:5,239,020–10,537,682 | 5.30 Mb |
| JU1793 frequency ≥ 0.70 under *mig-6* | X:4,916,851–16,472,067 | 11.56 Mb |

## 3. Data sources

| source | what it gave |
| --- | --- |
| `/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213/bcsq.vcf.gz` | BCSQ-annotated CeNDR isotype VCF, 611 isotypes, JU1793 and JU2466 genotypes and CeNDR allele frequencies |
| `Annotations/c_elegans.PRJNA13758.WS283.csq.gff3.gz` | chrX gene coordinates and locus names (6,050 genes; 1,613 with a public name) |
| `CeNDR/expression/Ce207expression.csv` | the 207-isolate expression matrix, 25,849 transcripts |
| `CeNDR/expression/ce207_qtl.tsv` | mapped eQTL with classification, interval and variance explained |
| `div_regions/FINAL_RECALLED_DIVERGENT_fixMissingHaps.tsv` | hyper-divergent region calls per isotype |
| GO Consortium `wb.gaf.gz` + `go-basic.obo` | unbiased gene-set selection by GO term name |

## 4. Parental variation on chromosome X

496,453 chrX sites, all PASS. **18,628 sites are homozygous-discordant between
JU1793 and JU2466** (11,915 with the alternate allele in JU1793, 6,713 in
JU2466), an average of **1.05 discordant sites per kb**. 1,159 distinct
(site, gene) pairs carry a protein-altering consequence, across **712 genes**.

**Positional overlap therefore carries no weight here.** The 2.43 Mb support
interval alone contains 113 genes with protein-altering parental differences;
X:4.9–13.1 Mb contains 318.

## 5. RNA-silencing genes: a clean negative

Genes were selected without a hand-curated list: GO terms were matched by name
against the ontology, then genes carrying those terms were taken from the
WormBase GAF. **Ten chrX genes carry an RNA-silencing GO term.**

| gene | position | LOD | JU1793 freq, control → *mig-6* | coding difference between parents |
| --- | --- | --- | --- | --- |
| **sid-5** | X:6.90 Mb | **767** | 0.09 → 0.91 | **none** |
| **chup-1** | X:8.79 Mb | 531 | 0.10 → 0.86 | **none** |
| **pgp-4** | X:11.36 Mb | 288 | 0.12 → 0.77 | **none** |
| alg-1 | X:13.95 Mb | 220 | 0.16 → 0.76 | none |
| sid-3 | X:17.19 Mb | 70 | 0.23 → 0.62 | none |
| stau-1 | X:1.03 Mb | 42 | 0.13 → 0.41 | none |
| nrde-3 | X:0.37 Mb | 31 | 0.13 → **0.37** | 4 missense, all JU1793 |
| mrp-1 | X:0.58 Mb | 33 | 0.13 → 0.38 | 1 missense, JU1793 |
| C30E1.4 | X:16.92 Mb | 76 | 0.22 → 0.62 | 1 missense, JU2466 |
| C08A9.3 | X:17.09 Mb | 71 | 0.22 → 0.62 | 1 missense, JU2466 |

**Every gene positioned where the JU1793 allele is enriched has no coding
difference. Every gene with a coding difference is positioned where it is not.**
*nrde-3*, the somatic nuclear Argonaute, carries the most attractive variant set
(H128Q, R467Q, E757A, E856D, all JU1793, two at CeNDR frequency 0.072) but sits
at X:374 kb where the JU1793 allele frequency is **0.372** — depleted, not
enriched. It fails the direction test.

**A structural point that constrains the whole hypothesis:** none of the
canonical somatic-RNAi genes is on chromosome X. *rde-1* (V), *rde-4* (III),
*rde-10* (I), *rde-11* (IV), *rde-12* (V), *ppw-1* (I), *rrf-1* (I), *mut-16*
(I), *sago-1* (V), *sago-2* (I), *sid-1* (V), *sid-2* (III), *drh-1/2* (IV). If
this is a somatic-RNAi-competence locus, it is not one of the known machinery
genes.

### *sid-5*

*sid-5* (X:6,904,721–6,905,849) is the best-positioned silencing gene on the
chromosome, at LOD 767 with a 0.09 → 0.91 allele-frequency swing, and it is a
genuine systemic-RNAi gene (Hinas et al. 2012). It is **not** a candidate here:

- 19 parent-discordant variants within ±5 kb, **all non-coding**. This
  reproduces the conclusion already recorded in
  `scripts/pooled_cross_candidate_variation.R`, which found zero protein-altering
  differences within ±2 kb and judged a chrX peak landing 260 bp away to be
  coincidence in a ~900 kb interval.
- **Expression:** JU1793 3.204 (36th percentile), JU2466 3.445 (71st),
  gap −0.64 SD — ordinary, at the 44.6th percentile of random strain pairs.
  h² = 0.000, H² = 0.255. **No mapped eQTL.**

The direction is the one the hypothesis wants (JU1793 lower), but the magnitude
is unremarkable and there is no additive heritability to support a segregating
*cis*-regulatory variant.

## 6. Deleterious JU1793 variants in the support interval

**29 JU1793-carried protein-altering sites across 26 genes. None is truncating** —
no stop-gained, frameshift, splice-acceptor/donor, start-lost or stop-lost.
28 missense plus one inframe deletion.

| gene | variant | CeNDR AF | isotypes | position |
| --- | --- | --- | --- | --- |
| F14B8.5 | Q24E | **0.0033** | 2/611 | X:6,928,580 |
| T22B7.7 | inframe del 17MS>17I | 0.0098 | 6/611 | X:5,653,693 |
| C38C5.1 | S28L | 0.0115 | 7/611 | X:5,558,140 |
| set-19 | Y861N | 0.0115 | 7/611 | X:5,676,663 |
| F49E10.4 | Y399F | 0.0115 | 7/611 | X:5,878,867 |
| W05H9.4 | T706A | 0.0115 | 7/611 | X:6,266,406 |
| lev-9 | T208I | 0.0115 | 7/611 | X:6,290,170 |
| fah-1 | R98G | 0.0115 | 7/611 | X:6,448,565 |
| C15B12.3 | E51G | 0.0115 | 7/611 | X:6,467,411 |
| C15B12.4 | G205R | 0.0115 | 7/611 | X:6,472,425 |
| gnrr-4 | R399Q | 0.0164 | 10/611 | X:5,740,067 |
| F43E12.1 | E80D, H81Y | 0.0168 | 10/611 | X:7,177,927–8 |
| C39D10.7 | E408K | 0.0230 | 14/611 | X:7,911,300 |
| hog-1 | I15V | 0.0443 | 27/611 | X:5,853,327 |

### The shared rare haplotype

**13 of the 29 sites sit at exactly CeNDR frequency 0.0115 — 7 of 611 isotypes —
spanning X:5,558,140–6,938,107 (1.38 Mb).** That is one rare haplotype block
carried by seven isotypes including JU1793. These variants are in near-complete
linkage and **cannot be separated genetically in this cross**; choosing among
them by sequence is not possible and would need recombinants or transgenics.

Two are worth naming. **F14B8.5 Q24E** is the rarest variant in the interval
(2 of 611) and lies 23 kb from *sid-5* (F14B8.2) in the same cosmid — though
"adjacent to a candidate gene" is precisely the positional coincidence this
repository's own analysis warns against. **set-19 Y861N** is a SET-domain
methyltransferase, the family that contains the nuclear-RNAi effectors *set-25*
and *set-32*, with a parental expression gap of −0.90 SD.

## 7. *mig-23*: the strongest positive, and it points the other way

***mig-23* (R07E4.4) is at X:5,945,813–5,948,869 — 6.2 kb from the peak.**

| strain | ce207 value | z | percentile |
| --- | --- | --- | --- |
| **JU1793** | 2.712 | −0.71 | **19.8th** |
| **JU2466** | 3.070 | +0.93 | **90.3rd** |
| N2 | 2.810 | −0.26 | 40.6th |
| XZ1516 | 2.920 | +0.25 | 68.1st |

Gap **−1.64 SD**, at the **85th percentile** of random strain pairs.
h² = 0.045, H² = 0.446, **no mapped eQTL**. Only 5 parent-discordant variants
within ±2 kb, all intronic or intergenic — no coding, promoter or UTR candidate.
Local divergence is 0.82/kb against a chromosome average of 1.05/kb, so this is
not a divergence desert either.

**The functional link is to *mig-6*, not to RNAi.** MIG-23 is a Golgi-resident
NDPase that links ADAM protease glycosylation to organ morphogenesis
(Nishiwaki et al. 2004); MIG-17 secretion and localisation to the distal tip
cell surface requires N-glycosylation and the basement-membrane protein
MIG-6/papilin, and *mig-6* interacts genetically with *mig-17* and collagen IV
(Kawano et al. 2009). *mig-23* and *mig-6* are therefore in the same distal-tip-
cell migration pathway.

**Caveats.** *mig-23* is not the largest gap in the neighbourhood: *klp-13*
(29.9 kb away) is −2.04 SD, and *R07E4.1* (4.0 kb away, 2 coding variants)
appears to be −2.71 SD but JU1793 sits exactly at the matrix detection floor
(−1.000 is the global minimum, 1.87% of all values; 7 of 207 strains are floored
for that gene), so its z is not trustworthy. And the direction is not the naive
prediction: JU1793 has *less* MIG-23 yet resists *mig-6* RNAi, which requires a
suppressor rather than a synthetic relationship. Parallel routes exist that could
buffer such a loss (Rahman et al. 2025).

## 8. *trans*-eQTL mapping into the region

3,360 distant eQTL are mapped in ce207, 459 of them to chromosome X. **77 peak
within the support interval, on 75 target genes** — about 1.3× the 64 expected
for an interval of that width, so a mild enrichment rather than a hotspot.

Two targets are RNAi-related, and **only these two eQTL intervals cover the QTL
peak**:

| target | eQTL peak | eQTL interval | logP | var_exp | covers X:5,955,109 | JU1793 − JU2466 |
| --- | --- | --- | --- | --- | --- | --- |
| **haf-6** (Y48G8AL.11a.1) | X:6,769,644 | X:303,285–7,975,338 | 11.6 | 0.21 | **yes** | **+1.34 SD** (73rd pct) |
| **haf-6** (Y48G8AL.11a.2) | X:1,612,592 | X:305,933–7,975,338 | 15.6 | 0.28 | **yes** | −0.06 SD (12th pct) |
| **mut-15** (T01C3.8a.1) | X:6,723,875 | X:5,590,940–7,798,884 | 6.5 | 0.14 | **yes** | +0.14 SD (9th pct) |

**HAF-6** is a *bona fide* RNAi gene — an ABC transporter required for efficient
RNAi that interacts genetically with *rde-2* and *mut-7* (Sundaram et al. 2006,
2008). It is the only result in this analysis that connects the locus to RNAi
machinery. Three things weaken it:

1. The eQTL interval is **7.7 Mb wide** (X:0.30–7.98 Mb), so "covers the peak"
   is close to uninformative on a 17.7 Mb chromosome.
2. **The direction is backwards.** HAF-6 is *required* for efficient RNAi;
   JU1793 has *more* of it, yet JU1793 is the resistant allele.
3. The three *haf-6* isoforms disagree (+1.34, −0.06, +0.88 SD).

*mut-15*'s eQTL interval is much tighter and essentially coincides with the
support interval, but the two parents barely differ at the target (+0.14 SD,
9th percentile of random pairs), so it carries no information about this cross.

No *mig-6*-pathway gene has a *trans*-eQTL covering the peak. *pat-3* has the
largest parental expression gap found anywhere in this analysis
(ZK1058.2.2, −2.61 SD, 95th percentile) but its eQTL maps to X:10.4–12.8 Mb, a
different locus.

## 9. What was ruled out

- **Hyper-divergence is not hiding the answer.** JU1793's divergent regions on
  chrX total roughly 45 kb (X:1.76, 5.05, 12.61, 14.22, 14.37 Mb) and none
  overlaps the peak or *sid-5*, so the absence of coding variants is real rather
  than an ascertainment artifact.
- **The known JU1580 RNAi lesion is on chromosome IV, not X — and JU1793
  carries it.** An earlier draft of this note said the lesion was "elsewhere"
  and gave IV:6,607,376–6,613,353. Those are the *drh-1* **gene** coordinates
  from the WS283 annotation, not the deletion's, and the statement was written
  from background knowledge with only the citations verified. The correction:
  the deletion is *niDf250*, 159 bp at **IV:6,607,635–6,607,793**, removing most
  of exon 19 and part of exon 20 and truncating the RIG-I C-terminal domain
  (Ashe et al. 2013); it is at 22/97 (23%) in wild isolates. In
  `VCFs/WI.MANTAsv.soft-filter.vcf.gz` it is called at IV:6,607,644–6,607,803
  (PASS) with **JU1793 1/1 and JU2466 not carrying it** (60 of 328 samples are
  1/1, 1 is 0/1, 267 are `./.`; 18.6% carries, matching the published 23%).
  **So *drh-1* segregates in this cross.** It falls inside the HT115g-POS1g
  chrIV support interval (IV:4,691,159–7,592,922, 189 kb from that peak, LOD
  56.6 against a chrIV maximum of 56.7) and 63 kb outside the HT115g-MIG6g
  interval, where the local LOD is still 150.8. The JU1793 allele frequency at
  *drh-1* runs 0.25 (control) → 0.78 (*mig-6*) → 0.61 (*pos-1*), i.e. enriched
  under **both** targets — the opposite of the chrX pattern and the signature of
  a general RNAi-response locus. Of the 7 homozygous JU1793 deletions among 16
  PASS structural variants differing between the parents in that interval, only
  two overlap coding sequence: *drh-1* and *srx-50*.

  **A methodological warning this episode earned.** An intermediate draft also
  claimed, from the SNV VCF, that JU1793 does *not* carry the deletion, on the
  basis that JU1793 is reference at all 215 variant sites in the gene with zero
  missing calls. That inference is invalid: **there are zero SNV records inside
  the 159 bp deletion window**, so call patterns outside it carry no information
  about the deletion. Deletion status must be read from the structural-variant
  call set, never inferred from SNV genotypes — which is the same limitation
  section 10 states and which the SNV-based reasoning ignored.

  **Still open:** DRH-1 is characterised as antiviral, dicing viral RNA. Whether
  *niDf250* affects the response to *exogenous* dsRNA is a separate question and
  is not assumed here.
- **JU1580 and JU1793 are the same isotype**, which this repository already
  documents (`METHODS.txt`, `FIGURE_CAPTIONS.txt`) and which is consistent with
  JU1580 being absent from the 611-isotype CeNDR set while JU1793 is present.

## 10. Limitations

1. **Structural variants are not ascertained in the chrX analysis.** Sections
   4-8 rest on a short-read SNV and indel call set. The *drh-1* result in
   section 9 shows what that misses: a 159 bp deletion with a published
   loss-of-function phenotype is entirely absent from the SNV VCF and is visible
   only in `VCFs/WI.MANTAsv.soft-filter.vcf.gz`. **The chrX interval has not
   been searched that way**, and it should be before the negative in section 5
   is treated as final.
2. **ce207 is whole-animal, one value per strain, unchallenged.** A difference
   that is intestine-specific, or that only appears on exposure to dsRNA, is
   invisible. There is no within-strain replication, so no p-value on any
   parental difference; gaps are reported as z-scores against the 207-strain
   distribution.
3. **The design confounds the two hypotheses.** One somatic and one germline
   target cannot separate "*mig-6*-specific" from "somatic-RNAi-specific".
4. **The control is skewed** (JU1793 at 0.091 before treatment).

## 11. What would settle it

- **A second somatic RNAi target.** If the locus is somatic-RNAi trafficking it
  appears for any somatic target; if it is *mig-6* pathway it does not. No amount
  of further sequence mining separates these.
- **Structural-variant calling** across X:4.9–13.0 Mb in JU1793.
- **An LD check on the 7/611 haplotype block**, to confirm formally that the
  13 variants are inseparable.
- **Direct measurement of *mig-23* and *sid-5* transcript levels** in JU1793 and
  JU2466, in the intestine and under dsRNA exposure, which is the condition ce207
  does not cover.

## 12. Citations

- **Ashe A, Bélicard T, Le Pen J, Sarkies P, Frézal L, Lehrbach NJ, Félix MA, Miska EA** (2013) A deletion polymorphism in the *Caenorhabditis elegans* RIG-I homolog disables viral RNA dicing and antiviral immunity. *eLife* 2:e00994. doi:10.7554/eLife.00994. PMID 24137537.
- **Félix MA, Ashe A, Piffaretti J, Wu G, Nuez I, Bélicard T, et al.** (2011) Natural and experimental infection of *Caenorhabditis* nematodes by novel viruses related to nodaviruses. *PLoS Biol* 9(1):e1000586. doi:10.1371/journal.pbio.1000586. PMID 21283608.
- **Hinas A, Wright AJ, Hunter CP** (2012) SID-5 is an endosome-associated protein required for efficient systemic RNAi in *C. elegans*. *Curr Biol* 22(20):1938–1943. doi:10.1016/j.cub.2012.08.020. PMID 22981770.
- **Kawano T, Zheng H, Merz DC, Kohara Y, Tamai KK, Nishiwaki K, Culotti JG** (2009) *C. elegans mig-6* encodes papilin isoforms that affect distinct aspects of DTC migration, and interacts genetically with *mig-17* and collagen IV. *Development* 136(9):1433–1442. doi:10.1242/dev.028472. PMID 19297413.
- **McEwan DL, Weisman AS, Hunter CP** (2012) Uptake of extracellular double-stranded RNA by SID-2. *Mol Cell* 47(5):746–754. doi:10.1016/j.molcel.2012.07.014. PMID 22902558.
- **Nishiwaki K, Kubota Y, Chigira Y, Roy SK, Suzuki M, Schvarzstein M, Jigami Y, Hisamoto N, Matsumoto K** (2004) An NDPase links ADAM protease glycosylation with organ morphogenesis in *C. elegans*. *Nat Cell Biol* 6(1):31–37. doi:10.1038/ncb1079. PMID 14688791.
- **Rahman MM, et al.** (2025) Distal tip cell migration mutants of *Caenorhabditis elegans* are rescued by bioequivalent outputs from chondroitin and N-glycosylation pathways. PMID 41197717. *(Author list and journal not verified beyond PubMed metadata — confirm before citing.)*
- **Sundaram P, Echalier B, Han W, Hull D, Timmons L** (2006) ATP-binding cassette transporters are required for efficient RNA interference in *Caenorhabditis elegans*. *Mol Biol Cell* 17(8):3678–3688. doi:10.1091/mbc.e06-03-0192. PMID 16723499.
- **Sundaram P, Han W, Cohen N, Echalier B, Albin J, Timmons L** (2008) *Caenorhabditis elegans* ABCRNAi transporters interact genetically with *rde-2* and *mut-7*. *Genetics* 178(2):801–814. doi:10.1534/genetics.107.081588. PMID 18245353.
- **Winston WM, Sutherlin M, Wright AJ, Feinberg EH, Hunter CP** (2007) *Caenorhabditis elegans* SID-2 is required for environmental RNA interference. *PNAS* 104(25):10565–10570. doi:10.1073/pnas.0611282104. PMID 17563372.
- **Wong MC, Schwarzbauer JE** (2012) Gonad morphogenesis and distal tip cell migration in the *C. elegans* hermaphrodite. *WIREs Dev Biol* 1(4):519–531. PMC3614366.

Every citation above was verified against PubMed except where noted.

## 13. Reproducing this

The analysis ran outside the figure pipeline and depends on files that are not
in the repository (the CeNDR BCSQ VCF, the WS283 annotation, the ce207
expression matrix and eQTL table, the divergent-region calls). Intermediate
tables were written to the session scratchpad:

| file | contents |
| --- | --- |
| `chrX_per_gene.tsv` | 712 chrX genes with protein-altering parental differences, with interval tier and carrier counts |
| `core_ju1793_variants.tsv` | the 29 JU1793-carried protein-altering sites in the support interval |
| `core_ju1793_genes.tsv` | the 26 genes those sit in, with ce207 expression gap and eQTL status |

If any of this reaches the manuscript, it should be rebuilt by a committed
script with pinned inputs, in the style of
`scripts/pooled_cross_candidate_variation.R`, rather than cited from this note.
