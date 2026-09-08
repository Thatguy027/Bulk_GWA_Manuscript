# claude_science — Figure 4C rebuild and the SID-1 comparison

Everything in this directory was produced in a Claude Science session and is
self-contained: it reads from `supplemental_data/structure/` and writes only
inside `claude_science/`. Nothing under `scripts/`, `plots/`, `data/` or
`supplemental_data/` at the repository root has been modified.

Run every script **from the repository root**.

## Two versions of Figure 4C

The panel exists in two variants that make different claims. Pick one; they
are not compatible in a single panel.

**`plots/Figure4C_charge.pdf` — the charge argument (recommended).**
Ectodomain coloured by the net charge of each residue's 12 Å neighbourhood at
pH 4.4, the T96 pocket with the two lysines that make it, and where T96K moves
that pocket in the domain-wide distribution. Rests on McEwan et al. 2012
showing that permanent positive charge at this surface *increases* uptake, and
does not depend on the model's global accuracy — only on which residues
neighbour T96 and on T96 being solvent-exposed, both of which the model gets
right.

**`plots/Figure4C_rebuild.pdf` — the confidence/proximity argument.**
Cartoon coloured by pLDDT, the T96 zoom in the same colouring, and the
proximity null as a histogram. Makes the weaker, defensive claim. Useful as a
supplementary panel showing the model is sound where it matters, since the
charge version carries no confidence encoding.

## Pipeline

    python3 claude_science/scripts/01_panelC_data.py          # tables
    python3 claude_science/scripts/02_render_charge.py         # charge cartoons
    Rscript  claude_science/scripts/03_figure4C_charge.R       # -> Figure4C_charge
    python3 claude_science/scripts/04_fold_comparison.py       # TM-align table
    Rscript  claude_science/scripts/05_figure_fold_comparison.R

    # the pLDDT variant
    python3 claude_science/scripts/alt_render_plddt_cartoon.py
    python3 claude_science/scripts/alt_render_plddt_zoom.py
    Rscript  claude_science/scripts/alt_figure4C_plddt.R

Step 01 must run before 02/03 and before the `alt_*` scripts. Step 04 needs
network access on first run and caches its downloads in `data/pdb_cache/`
(git-ignored, re-downloadable).

Every render script imports its geometry, camera, ribbon construction, stick
drawing and framing from `scripts/sid2_ribbon_render.py` and
`scripts/sid2_zoom_render.py` via `sys.path.insert(0, "scripts")`. Only the
colour source and the marker set differ from the released panel, so these are
the same figures in a different colour space. If those two modules change,
these panels follow.

Requirements: `numpy`, `pandas`, `scipy`, `matplotlib`, `biopython`, and
`tmtools` (step 04 only; `pip install tmtools`). R side: `tidyverse`,
`patchwork`, `ggtext`, `png`, `ggrepel`, `RColorBrewer`. Running the R scripts
through `Rscript` in a non-UTF-8 locale prints harmless "unable to translate"
warnings for Å and β; the glyphs still render. `LC_ALL=en_US.UTF-8 Rscript …`
silences them.

## The two corrections this work produced

**The uptake-critical residue set is three histidines.** McEwan, Weisman &
Hunter 2012 (*Mol Cell* 47:746, doi:10.1016/j.molcel.2012.07.014) identify
**His32, His168 and His175** and never mention residue 34. This model's
ectodomain contains exactly three histidines, at those three positions, which
independently confirms the numbering — so the "verify residue positions
against UniProt" item on the pre-submission list is closed, and UniProt was
never the right source. Residue 199 from that UniProt-derived list is in the
transmembrane helix and does not belong.

Consequences still to apply at the repository root:

- `FIGURE_CAPTIONS.txt` (~line 327) describes H32, D34 and H168 together as
  "residues with a published effect on dsRNA uptake". D34 is the qt13
  loss-of-function allele — a separate line of evidence needing its own
  citation, not part of the histidine set.
- `scripts/sid2_zoom_render.py` line 76 has `FUNC = {32, 34, 168}`: it omits
  H175 and includes D34.
- `METHODS.txt` (~line 453) reports three of four annotated residues within
  20 Å, binomial p = 0.19, permutation p = 0.30. With the corrected
  three-histidine set it is **2 of 3 within 20 Å, binomial p = 0.38**. The
  proximity claim is better dropped than defended.

**SID-2's ectodomain is not a SID-1 β-strand-rich domain.** See
`plots/sid2_vs_sid1_fold.pdf` and `SID2_vs_SID1_memo.md`. An unrelated
immunoglobulin domain scores higher against SID-2 (TM 0.427) than SID-1's BRD1
does (0.404), while the two genuine BRDs score 0.525 against each other. Do
not superpose the 8XC1 dsRNA onto the SID-2 model.

## Contents

    SID2_vs_SID1_memo.md    the full write-up: methods, numbers, caveats,
                            and the charge-matched T96R / T96Q test proposed
    scripts/                see the pipeline above
    data/                   derived tables (regenerate with 01 and 04)
    plots/                  figures
    articles/               empty; a place for the two source PDFs if wanted

`data/sid2_panelC_residues.tsv` is the per-residue table behind every panel:
pLDDT, secondary structure, topology, coordinates, Cα distance from T96, SASA
and local charge at both pH values.
