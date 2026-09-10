# Reproduction environment for Bulk_GWA_Manuscript

Two conda environments, pinned to the versions `METHODS.txt` records for the
build machine. Verified against the committed figures rather than assumed:
the R stack reproduces a figure PNG byte-for-byte, and the Python stack
reproduces 12 of the 17 structure assets byte-for-byte.

Built 2026-09-09 against repository HEAD `8c140d8`.

## `bulkgwa-r` — figures, deposit builders, invariant checks

R 4.5.2, matching `METHODS.txt` exactly.

| package | installed | METHODS.txt | |
|---|---|---|---|
| R | 4.5.2 | 4.5.2 | match |
| ggplot2 | 4.0.3 | 4.0.3 | match |
| ggtext | 0.1.2 | 0.1.2 | match |
| patchwork | 1.3.2 | 1.3.2 | match |
| scales | 1.4.0 | 1.4.0 | match |
| png | 0.1.8 | 0.1.8 | match |
| jsonlite | 2.0.0 | 2.0.0 | match |
| tidyverse | 2.0.0 | 2.0.0 | match |
| BEDMatrix | 2.0.4 | 2.0.4 | match |
| abind | 1.4.8 | 1.4.8 | match |
| data.table | 1.18.6.1 | 1.18.2.1 | patch differs |
| ggrepel | 0.9.8 | 0.9.7 | patch differs |
| RcppML | 0.3.7.1 | 0.3.7 | patch differs |
| extraDistr | 1.10.0.5 | 1.10.0 | patch differs |

The four deltas are the versions conda-forge carries; the exact patch releases
`METHODS.txt` names are not packaged there. Every renderer-relevant package —
ggplot2, ggtext, patchwork, scales, png — is exact, which is what the
byte-identity result below turns on.

Also present, unpinned because `METHODS.txt` does not pin them: GGally 2.4.0,
broom 1.0.13, lme4 2.0.6, rmarkdown 2.32 (for knitting `FIGURE_REPORT.Rmd`).

## `bulkgwa-py` — structure analysis and rendering

Every package matches `METHODS.txt` exactly.

| package | installed | METHODS.txt | |
|---|---|---|---|
| Python | 3.13.2 | 3.13.2 | match |
| NumPy | 2.4.3 | 2.4.3 | match |
| SciPy | 1.17.1 | 1.17.1 | match |
| Biopython | 1.86 | 1.86 | match |
| Matplotlib | 3.10.8 | 3.10.8 | match |
| Pillow | 12.1.1 | 12.1.1 | match |

Plus `tmtools` (pip) for the TM-align step in `claude_science/scripts/04_fold_comparison.py`.

Pinning mattered here. An unpinned build (NumPy 2.5.3, SciPy 1.18.0,
Matplotlib 3.11.1, Pillow 12.3.0) failed to reproduce three assets that the
pinned build reproduces exactly.

## Not installed

PLINK 2.00a3 / 1.90b6.21 and bcftools 1.11. Nothing in the curated figure set
needs them: the two scripts that call `bcftools` resolve hard-coded absolute
paths to a CeNDR release outside the repository, and the PLINK panel lives in
the Dryad archive rather than a clone.

## Verification performed

1. **`SUPP_FIG_XX_dilution_validation.R`** — the repository's pinned-literal
   canary — runs clean, all asserted numbers agree, and the PNG it writes is
   **byte-identical** to the committed one.
2. **`bash scripts/check_repo_invariants.sh`** — full fast suite, **0 checks
   failed**. Both interpreters must be on PATH; see below.
3. **`sid2_ribbon_render.py` and `sid2_zoom_render.py`** — 12 of 17 assets
   byte-identical. The 5 that differ are content differences, not environment
   artifacts; see the note below.

Consequence: this machine can now rebuild figures that match `plots/`, which is
the hazard `SYNC_MANIFEST.md` raises under "The R environment is not the same on
both machines". That section is now out of date for this machine.

## Running the invariant suite

The suite shells out to both `Rscript` and `python3`, so it needs both
environments visible at once:

```sh
ENVS=~/.claude-science/conda/envs
PATH="$ENVS/bulkgwa-r/bin:$PATH:$ENVS/bulkgwa-py/bin" \
  LC_ALL=en_US.UTF-8 bash scripts/check_repo_invariants.sh
```

`LC_ALL=en_US.UTF-8` silences the harmless "unable to translate" warnings for
the Å and β glyphs.

## The five assets that do not reproduce

`sid2_ecd_ribbon_charge.png`, `sid2_ecd_ribbon_func.png`,
`sid2_overview_oriented.png`, `sid2_overview_oriented_prev.png`,
`sid2_zoom_t96.png`.

These are not a version problem — 12 sibling assets from the same two scripts
in the same environment reproduce byte-for-byte. Two contributing causes:

- **The committed rasters predate their script.**
  `sid2_ecd_ribbon_charge.png` and `sid2_ecd_ribbon_func.png` were last
  committed in `45e9c39` (09-03); their producer `sid2_ribbon_render.py` was
  last changed in `fc922f9` (09-04).
- **The residue-set correction was applied to one renderer, not both.**
  Commit `29c7075`, "Correct the Figure 4C residue set: three histidines, not
  four residues", rewrote `sid2_zoom_render.py` to split `HIS = {32, 168}` from
  `ALLELE = {34}` and documented why. `sid2_ribbon_render.py` still carries the
  uncorrected grouping at lines 87-91:

  ```python
  MARK_FUNC = {
      96:  ("T96",  ...),
      32:  ("H32",  ..., COL_FUNC),
      34:  ("D34",  ..., COL_FUNC),
      168: ("H168", ..., COL_FUNC)}
  ```

  D34 is drawn in the same class as the histidines and H175 is absent — the
  exact error `29c7075` fixed elsewhere.

### Why this reaches a figure

`scripts/SUPP_FIG_XX_sid2_electrostatics.R:58` embeds
`plots/assets/sid2_ecd_ribbon_charge.png` as an `annotation_raster`. That raster
is rendered with `marks=MARK_FUNC`, so the supplement's structure panel labels
**H32, D34 and H168, and omits H175**.

The same script's header (lines 78-89) states the corrected position — three
extracellular histidines, D34 a separate line of evidence, "the statistic in
panel C is computed over the HISTIDINES only". So the script's prose and the
image it embeds disagree, which is invariant 4 in
`.claude/agents/repo-reviewer.md`.

The separate charge renders from `sid2_charge_render.py`
(`sid2_overview_charge.png`, `sid2_zoom_charge.png`) are correct — they use
`HIS = {32, 168, 175}` and `BASIC = {93, 132}`, they match their caption, and
both reproduce byte-identically. Only the ribbon renderer is affected.

Not fixed here, only reported.
