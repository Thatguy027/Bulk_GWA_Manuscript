# Does the SID-2 ectodomain resemble the SID-1 dsRNA-binding domains?

Reference: Wang R, Cong Y, Qian D, Yan C, Gong D. *Structural basis for
double-stranded RNA recognition by SID1.* Nucleic Acids Res 2024;52(11):6718.
doi:10.1093/nar/gkae395. PMID 38742627. Coordinates: PDB **8XBS** (apo cSID1,
cryo-EM 2.21 Å) and **8XC1** (cSID1 + dsRNA, 2.21 Å).

Note: that paper does not mention SID-2 anywhere. Any comparison here is ours.

## What the SID-1 structure establishes

cSID1 is a dimer; each subunit has an extracellular region built from two
β-strand-rich domains, **BRD1 (residues 18–178)** and **BRD2 (179–310)**,
followed by a transmembrane domain. BRD1 is 2 short α-helices plus 11
β-strands; BRD2 is 2 short α-helices plus 8 β-strands. The dsRNA sits at the
BRD1/BRD2 interface, nearly parallel to the membrane plane, and each dimer
engages two dsRNA duplexes. Recognition is **electrostatic and
sequence-independent**: ionic contacts between basic side chains and the
phosphate backbone, hydrogen bonds to 2′-OH groups, and three basic surface
regions seating into three successive major grooves.

The size coincidence that prompted the comparison is real: SID-2's ectodomain
is residues 21–188 (168 residues, 72 assigned to β-strand and none to helix by
the Kabsch–Sander call used in Figure 4), against BRD1's 18–178.

## The fold comparison does not support homology

TM-align, query = the well-modelled core of the SID-2 ectodomain model
(pLDDT ≥ 70, 117 of 168 residues), TM-score normalised by that query.

| Target | TM-score | RMSD (Å) | % identity |
|---|---|---|---|
| cSID1 BRD1 vs BRD2 — *internal reference* | **0.525** | 2.75 | 10.4 |
| Immunoglobulin domain, 12E8 chain H — *control* | 0.427 | 4.03 | 5.0 |
| cSID1 BRD1 (18–178) | 0.404 | 4.23 | 5.2 |
| cSID1 BRD2 (179–310) | 0.386 | 4.46 | 5.0 |
| Fibronectin type III, 1TEN — *control* | 0.385 | 3.72 | 5.7 |
| Galectin CRD, 2JJ6 — *control* | 0.337 | 4.61 | 2.9 |
| cSID1 TMD (311–776) — *control* | 0.333 | 4.99 | 6.7 |
| Legume lectin, 1LOB — *control* | 0.317 | 4.35 | 1.6 |

Two comparisons decide it. An **unrelated immunoglobulin domain scores higher
against SID-2 (0.427) than BRD1 does (0.404)**, and a fibronectin type III
domain scores the same as BRD2. Meanwhile the two genuine BRDs of cSID1 — a
real structural relationship — score 0.525 against each other at 2.75 Å. So
SID-2's resemblance to the SID-1 BRDs is the background level for any compact
β-sandwich of this size, not a specific relationship. Sequence identity across
every alignment is 1.6–10.4%.

**Conclusion:** the visual similarity is a fold-class resemblance — both are
compact β-rich extracellular domains — and it does **not** license transferring
SID-1's dsRNA-binding geometry onto SID-2 by homology. Do not superpose the
8XC1 dsRNA onto the SID-2 model.

## The inference that does survive

SID-1's mechanism is electrostatic, not shape-specific: basic patches against
the phosphate backbone, sequence-independent. That is transferable as
*physics* without any fold relationship, and it makes the charge state of
SID-2's lumenal surface the relevant question.

Measured on the membrane-oriented ectodomain (chain A, 21–188):

- The ectodomain is **net acidic: 10 basic vs 18 acidic residues, net −8**;
  the lumenal cap (z > 60 Å) is net −4. SID-2's ectodomain is therefore *not*
  a SID-1-like polybasic dsRNA-binding surface overall.
- But **T96's own neighbourhood is one of the few basic pockets in it**: within
  12 Å there are 2 basic residues (K93, K132) and 1 acidic, net +1, which is
  the **87th percentile** of local net charge across the ectodomain (domain
  median −1).
- T96 is solvent-exposed at the median level for the domain (SASA 62.4 Å²,
  48th percentile), so the side chain is available to solvent.
- **T96K takes that pocket from +1 to +2**, and the phenotype runs the
  consistent direction: 96K increases pos-1 RNAi sensitivity in all three
  backgrounds tested (JU1793 0.95 → 0.53 hatched, JU2466 0.19 → 0.05,
  N2 0.32 → 0.04).

So the defensible statement is a charge hypothesis, not a homology one: *in a
pathway where dsRNA recognition is known to be electrostatic and
sequence-independent, the T96K substitution adds a positive charge at the most
basic solvent-exposed pocket of an otherwise acidic lumenal domain, and it
increases dsRNA sensitivity.* This is consistent with — and independent of —
the Coulombic surface analysis already in
`data/structure_modeling/claude_docking/stage4_electrostatics/`.

## Addendum: the primary SID-2 paper resolves the residue set

McEwan DL, Weisman AS, Hunter CP. *Uptake of extracellular double-stranded RNA
by SID-2.* Mol Cell 2012;47(5):746–754. doi:10.1016/j.molcel.2012.07.014.
PMID 22902558.

**The uptake-critical set is three histidines, not four residues.** The paper
states that SID-2 contains three extracellular histidines — **His32, His168 and
His175** — and tests exactly those. It never mentions residue 34. Our model's
ectodomain (21–188) contains **exactly three histidine residues, at 32, 168 and
175**, which independently confirms the numbering and closes the
"verify residue positions" item on the pre-submission list. UniProt is not the
right source for these; this paper is.

Consequences for Figure 4C and METHODS:

- D34 is the **qt13 loss-of-function allele**, a separate line of evidence from
  the histidine mutagenesis. It should be cited separately and drawn as a
  separate class, not folded into "residues with a published effect on dsRNA
  uptake".
- With the corrected set the proximity statistic gets **weaker**: 2 of 3
  histidines within 20 Å, binomial **p = 0.38** (was 3 of 4, p = 0.19–0.20).
  The proximity argument should be dropped rather than defended.
- `scripts/sid2_zoom_render.py` line 76 has `FUNC = {32, 34, 168}` — it omits
  H175 and includes D34, so the released zoom shows neither the full histidine
  set nor a clean one.

**Why histidines, and why this matters for T96K.** The paper targeted
histidines because the imidazole group is protonated only in the acidic
conditions SID-2 requires; SID-2-dependent transport needs an acidic
extracellular environment and is selective for dsRNA of at least 50 bp. His→Ala
and His→Glu each reduced transport. That is the same pH story already built in
the supplementary electrostatics panel (net −8.4 e at pH 7.4, −0.2 e at pH 4.4,
pI 4.38).

**The precedent that makes the T96K charge argument strong.** Individual His→Arg
substitutions retained more function than alanine but did not reach wild-type
levels — and the **triple His→Arg mutant internalised more dsRNA than wild type
(p < 0.05)**. Replacing pH-dependent positive charge with permanent positive
charge *increased* uptake. T96K adds a permanent positive charge to the same
lumenal domain and increases RNAi sensitivity in all three backgrounds. The
published experiment and ours point the same way, which is a far better footing
for Figure 4C than proximity. The paper also notes these histidines are not
conserved in Cbr-SID-2, and cites the TLR3 ectodomain, where multiple
histidines contact the RNA phosphate backbone.

## Numbers behind the charge version of Figure 4C

Local net charge = Henderson–Hasselbalch side-chain charge summed over all
ectodomain residues with a Cα within 12 Å, same pKa set as the electrostatics
supplement (Asp 3.9, Glu 4.25, His 6.0, Lys 10.5, Arg 12.5, Cys 8.3, Tyr 10.1).

At **pH 4.4**, the gut-lumen pH at which SID-2 functions:

- Whole ectodomain (21–188, the modelled span): **+0.47 e**, against −8.34 e at
  pH 7.4. The supplement's −0.2 e is for 21–193, so the two differ by the five
  unmodelled residues; quote the supplement's figure in the text.
- **T96's pocket: +1.24 e — the 82nd percentile** of the domain (median 0.00).
- **T96K takes it to +2.24 e — the 98th percentile.**
- The pocket is made by **K93 at 6.6 Å and K132 at 6.8 Å** (4.5 and 4.4 Å
  nearest heavy atom), the two nearest basic residues.
- The three uptake-critical histidines sit at −0.15 (H32), +0.85 (H168) and
  +1.46 (H175) — i.e. H175, the one 37.8 Å from T96, is in the most positive
  environment of the three.
- 81 of 168 ectodomain residues have positive local charge at this pH, so a
  positive pocket per se is unremarkable; being at the 98th percentile is the
  claim.

## Caveats

- The SID-2 side is a **prediction**, not a structure. A low-confidence model
  depresses TM-score, so the negative fold result is partly confounded by model
  quality; the Ig and fn3 controls are what make it interpretable, since they
  are scored against the same imperfect query.
- Charge counting is formal (Asp/Glu −1, Lys/Arg +1) and ignores pKa shifts,
  local dielectric and glycan shielding. Three of the nine N-glycosylation
  sequons lie in the ectodomain.
- Whether SID-2 contacts dsRNA at all is not established by anything here.
  The charge hypothesis is a reason to test, not a result.
- Directionality is a correlation across three backgrounds, not a
  demonstration that the charge is what causes it. The obvious discriminating
  experiment is a charge-matched control: T96R should behave like T96K if
  charge is the mechanism, and T96Q (isosteric, neutral) should behave like the
  wild type.
