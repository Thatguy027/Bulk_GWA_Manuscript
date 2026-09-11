# Is the N94-C95-T96 sequon glycosylated? — source material for the Discussion

Drafted 2026-09-11. **This is not a draft export.** It is reference material to
be pulled into the Google Doc, in the same category as `METHODS.txt`. Keep it
at the repository root: `check_manuscript_numbers.py` globs `manuscript/*.md`
and `manuscript/*.txt`, so a file placed there is checked as though it were an
export of the draft, and these numbers are deliberately finer-grained than the
ones `FIGURE_REPORT.md` attests.

The corresponding METHODS passage is the block beginning `IS THE N94 SEQUON
ACTUALLY GLYCOSYLATED?` in `METHODS.txt`.

---

## 1. Citation list

### Datasets queried directly

**Hu X, Xiang W, Qin S, Zhang X, Tian Z (2026)** Integrative quantitative
proteomics and N-glycoproteomics reveal age-processing of N-glycosylation in
*Caenorhabditis elegans*. *Proteomics* 26(9):39–48.
doi:10.1002/pmic.70133. PMID 42050943; PMC13519370.
*Queried:* Supplementary Tables 1–7 (`pmic70133-sup-0002-tables.rar`), 22
sheets. 369 glycoforms at 196 N-glycosites on 161 glycoproteins; N2 under FUdR
at days 1, 5 and 10 of adulthood; ZIC-HILIC enrichment with TMT labelling.
*Raw data:* ProteomeXchange **PXD067150** — reserved but **not yet released**
as of 11 Sep 2026 (absent from ProteomeXchange, PRIDE, iProX and MassIVE).

**Tan CH, Wang TY, Park H, Lomenick B, Chou TF, Sternberg PW (2024)**
Single-tissue proteomics in *Caenorhabditis elegans* reveals proteins resident
in intestinal lysosome-related organelles. *PNAS* 121(25):e2322588121.
doi:10.1073/pnas.2322588121. PMID 38861598; PMC11194598.
*Queried:* Datasets S1–S7 (`pnas.2322588121.sd01`–`sd07.xlsx`); and PRIDE
**PXD047792**, search-result file `James_20220518-G-I_chimerys20230616.msf`
(9.53 GB), read at PSM level.
*Search settings:* oxidation (M) dynamic, carbamidomethylation (C) fixed — **no
glycan modifications**.

**Narayan V, Ly T, Pourkarimi E, Murillo AB, Gartner A, Lamond AI, Kenyon C
(2016)** Deep proteome analysis identifies age-related processes in
*C. elegans*. *Cell Systems* 3(2):144–159. doi:10.1016/j.cels.2016.06.011.
PMID 27453442; PMC5003814.
*Queried:* Supplementary Tables 1–4 (`mmc2`–`mmc5.xlsx`). 9,398 proteins
identified, 7,380 quantified; N2 under FUdR at days 1, 5 and 10.
*Raw data:* ProteomeXchange **PXD004584** / MassIVE **MSV000079263** — the
ProteomeXchange record is annotated "No PTMs are included in the dataset".
*Search settings:* deamidation (N/Q), oxidation (M), pyro-Glu and protein
N-terminal acetylation variable; N-ethylmaleimide (C) fixed — **no glycan
modifications**.

### Annotation databases queried

- UniProt Knowledgebase, reviewed entry **G5EEV9** (`SID2_CAEEL`) — no
  `CARBOHYD` feature.
- **GlyGen** — no record for G5EEV9.
- **GlycoProtDB** / IGOT *C. elegans* N-glycoproteome (`jcggdb.jp`) — resource
  withdrawn; all known URLs return 404.
- **WormBase** — gene **ZK520.2** (*sid-2*).

### Primary literature full-text searched (not datasets)

- **Winston WM, Sutherlin M, Wright AJ, Feinberg EH, Hunter CP (2007)**
  *Caenorhabditis elegans* SID-2 is required for environmental RNA
  interference. *PNAS* 104(25):10565–10570. doi:10.1073/pnas.0611282104.
  PMID 17563372.
- **McEwan DL, Weisman AS, Hunter CP (2012)** Uptake of extracellular
  double-stranded RNA by SID-2. *Mol Cell* 47(5):746–754.
  doi:10.1016/j.molcel.2012.07.014. PMID 22902558.
- **Braukmann F, Jordan D, Jenkins B, Koulman A, Miska EA (2021)** SID-2
  negatively regulates development likely independent of nutritional dsRNA
  uptake. *RNA Biol* 18(6):888–899. doi:10.1080/15476286.2020.1827619.
  PMID 33044912.

None of the three contains glycosylation terminology across eleven search
terms; the single "glycan" hit in Braukmann et al. is a title in its reference
list.

### Cited via Tan et al., not queried directly

- **Han S, Schroeder EA, Silva-García CG, Hebestreit K, Mair WB, Brunet A
  (2017)** Mono-unsaturated fatty acids link H3K4me3 modifiers to *C. elegans*
  lifespan. *Nature* 544(7649):185–190. doi:10.1038/nature21686.
  PMID 28379943. — the dissected intestine-against-gonad transcriptome Tan
  et al. use as their intestine-enriched comparison set.

---

## 2. Statement for the Discussion — full version

> We note that we provide no experimental evidence that the N94-C95-T96
> N-glycosylation sequon is glycosylated, and we are aware of no published
> evidence that it is; the sequon is a sequence prediction. SID-2 carries no
> carbohydrate annotation in its reviewed UniProt entry (G5EEV9), has no record
> in GlyGen, and glycosylation is not mentioned in any of the primary SID-2
> papers (Winston et al. 2007; McEwan et al. 2012; Braukmann et al. 2021). The
> only *C. elegans* N-glycoproteome we were able to obtain and search — 369
> glycoforms at 196 N-glycosites on 161 glycoproteins from N2 adults at days 1,
> 5 and 10 (Hu et al. 2026) — does not contain SID-2, but that absence is
> uninformative rather than negative: the accompanying unenriched proteome
> comprises only 3,966 peptides from 870 protein accessions, and just 53 of
> those 161 glycoproteins appear in it, so SID-2 was never a realistic
> candidate for glycopeptide enrichment. Deeper proteomes of the same strain at
> the same ages do detect SID-2 (Narayan et al. 2016: 6, 5 and 5 peptides
> across three replicates, 20.3% coverage, median abundance among 9,398
> proteins), and single-tissue proteomics assigns it specifically to the
> intestine, detected in all six dissected intestines and absent from all six
> gonads (Tan et al. 2024). Neither of those studies searched for glycan
> modifications, however, so neither can address occupancy. In the deposited
> spectra of the intestinal dataset (PXD047792), SID-2's confident
> identification rests on two cytoplasmic-tail peptides covering 7.1% of the
> protein, and the tryptic peptide spanning N94-C95-T96 is not identified in
> any form. The ectodomain is therefore effectively unsampled by existing data,
> and settling occupancy would require targeted glycopeptide-enriched mass
> spectrometry or PNGase F sensitivity of a tagged protein rather than
> reanalysis.
>
> Two points bound the interpretation. First, cysteine at the central position
> of a sequon is unusual but not prohibitive: 7 of the 182 N-glycosites in Hu
> et al. for which the sequon can be read within the identified peptide are
> N-C-[ST], five of them N-C-T, so this sequon context is demonstrably occupied
> *in vivo* in *C. elegans*. Second, and more importantly, glycan occupancy is
> not the mechanism at issue here: removing the acceptor asparagine (N94A)
> leaves animals fully resistant to ingested *pos-1* dsRNA (99.3% hatching
> against 94.8% for wild type, *p* = 0.004) rather than phenocopying the
> sensitivity conferred by 96K (53.1%), and *C. afra*, which naturally carries
> the same AxT configuration, responds to ingested dsRNA.

## 3. Statement for the Discussion — compact version

> We note that we provide no experimental evidence that the N94-C95-T96 sequon
> is glycosylated, and we are aware of no published evidence that it is. SID-2
> has no carbohydrate annotation in UniProt (G5EEV9) or GlyGen, and
> glycosylation is absent from the primary SID-2 literature (Winston et al.
> 2007; McEwan et al. 2012; Braukmann et al. 2021). The one searchable
> *C. elegans* N-glycoproteome (Hu et al. 2026) lacks SID-2, but its unenriched
> proteome covers only 870 protein accessions, so the absence reflects sampling
> depth: proteomes that do detect SID-2 (Narayan et al. 2016; Tan et al. 2024)
> searched no glycan modifications, and in the latter's deposited spectra the
> peptide spanning N94-C95-T96 is never identified. Cysteine at the central
> position is not prohibitive — 7 of 182 assignable *C. elegans* glycosites are
> N-C-[ST] — but occupancy is not the mechanism at issue, since N94A leaves
> animals fully resistant rather than phenocopying 96K, and *C. afra* is a
> naturally occurring AxT that responds to ingested dsRNA.

---

## 4. Provenance of every number above

Each was measured from the source named, not taken from the papers' prose.

| Number | Where it comes from |
| --- | --- |
| 369 glycoforms, 196 glycosites, 161 glycoproteins | Hu et al. Table S1, `Content` sheet |
| 3,966 peptides from 870 accessions | Hu et al. Table S3, `proteome_matrix`, counted |
| 53 of 161 glycoproteins in their own unenriched proteome | intersection of Hu et al. S1 `161proteins` with S3 `proteome_matrix` accessions |
| 7 of 182 sites are N-C-[ST], 5 of them N-C-T | Hu et al. S1 `196glycosites`; the modified Asn's in-peptide index read from `ID_Comp`, not from the first sequon match. The 14 unassignable sites place the modified Asn at the penultimate residue, so the +2 position falls outside the peptide |
| named N-C-T sites: *mig-6* (O76840) 445 and 638, Y45F10C.4 (O45944) 62, O01454 84, *irg-7* (A0A131MBU3) 279; N-C-S: Q22720 360, Q09967 250 | same sheet; gene names via Hu et al. Table S5 `ProteinID`→`Gene` |
| 9,398 identified / 7,380 quantified | Narayan et al. `mmc3` and `mmc4` row counts |
| SID-2: 6, 5, 5 peptides; 20.3% coverage; rank 4,387 of 9,398 (mean intensity 2.09e9 against a median of 1.74e9) | Narayan et al. `mmc3` and `mmc2` |
| 50.7% membrane / 40.0% cytosol / 9.0% nucleus; >75% needed to assign a compartment | Narayan et al. `mmc2`, `Localization Profiles` sheet and the `Membrane` sheet's own description |
| Narayan and Tan search settings; "No PTMs are included in the dataset" | Narayan Supplemental Experimental Procedures (`mmc1.pdf`) and the ProteomeXchange record for PXD004584; Tan et al. Materials and Methods |
| SID-2 in all 6 single intestines, 0 of 6 gonads; log2 intestine/gonad 6.64, adjusted *p* = 2.5e-17; unchanged in *glo-1(lf)* (1.15x, *p* = 0.74) | Tan et al. `sd02` `Fig. 3B` column A; `sd07` `Batch#1`–`Batch#3-3`; the PD export in PXD047792 |
| 5 confident PSMs, 2 peptides, 22 residues, 7.1% coverage, all from intestine runs, both peptides cytoplasmic; `NCTFTANYTGYFTPDPK` absent in any form | `scripts/query_sid2_psms_pride.py` — reads the PSM table out of the 9.53 GB `.msf` in PXD047792 over HTTP range requests (538 requests, 2.20 MB, 0.023% of the file) and asserts each of these |
| N94A 99.3% hatching against 94.8% wild type, *p* = 0.004; 96K 53.1% | the editing series in `scripts/SUPP_FIG_XX_sid2_ortholog_conservation.R` |
| *C. afra* is a natural AxT and responds to ingested dsRNA | `supplemental_data/structure/sid2_ortholog_window_survey.tsv` and `sid2_env_rnai_sensitivity.tsv` |

---

## 5. Check before submission

1. **Hu et al. pagination.** PubMed gives `26(9):39–48`; Wiley's e-locator is
   `e70133`. Confirm which the target journal wants. The DOI is unambiguous
   either way.
2. **PXD067150 may be released.** The "not yet accessible" check is dated
   11 Sep 2026. If those spectra open, a targeted extraction of the four SID-2
   peptides spanning N94-C95-T96 against their reported glycan compositions
   becomes possible, and the sentence about needing a targeted experiment
   rather than reanalysis would need softening.
3. **SIDT1/SIDT2 are not SID-2 orthologs.** They are relatives of SID-1. Their
   N-glycosylation is published but is not evidence about SID-2, and should not
   be cited as though it were.
