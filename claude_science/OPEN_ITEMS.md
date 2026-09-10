# Open items found by Claude Science, 2026-09-09 and 2026-09-10

Findings from two sessions that were reported in conversation but had no home
in the repository. Nothing here has been fixed. Each item names the file and
line so it can be checked before it is acted on.

Related: `ENVIRONMENT.md` (how to rebuild figures faithfully on this machine),
`results_review/` (the Results-section claim audit).

## 1. The sid-2 residue-set correction was applied to one renderer, not both

Commit `29c7075` ("Correct the Figure 4C residue set: three histidines, not
four residues") rewrote `scripts/sid2_zoom_render.py` to separate `HIS = {32,
168, 175}` from the *qt13* allele at D34. The same correction was never applied
to `scripts/sid2_ribbon_render.py`, whose `MARK_FUNC` (lines 87-91) still draws
D34 in the same class as the histidines and omits H175 entirely.

`scripts/SUPP_FIG_XX_sid2_electrostatics.R:58` embeds that renderer's
`plots/assets/sid2_ecd_ribbon_charge.png` as an `annotation_raster`, so the
supplement's structure panel labels **H32, D34 and H168 and shows no H175** —
while the same script's own header (lines 78-89) states the corrected
three-histidine position and says the statistic is computed over the histidines
only. Script prose contradicts the image it embeds, which is invariant 4 in
`.claude/agents/repo-reviewer.md`.

The separate `scripts/sid2_charge_render.py` outputs
(`sid2_overview_charge.png`, `sid2_zoom_charge.png`) are correct, match their
caption, and rebuild byte-identically. Only the ribbon renderer is affected.
`sid2_ecd_ribbon_func.png` comes from the same script and is consumed by
nothing.

Fix is a one-line change to `MARK_FUNC` plus a rebuild of that script's assets,
which changes tracked files under `plots/assets/`.

## 2. Five committed structure assets do not reproduce

`sid2_ecd_ribbon_charge.png`, `sid2_ecd_ribbon_func.png`,
`sid2_overview_oriented.png`, `sid2_overview_oriented_prev.png`,
`sid2_zoom_t96.png`. Twelve sibling assets from the same two scripts in the
same pinned environment reproduce byte-for-byte, so this is content drift, not
a renderer-version artifact. Two contributing causes: item 1 above, and the
committed rasters predating their producer — the two ribbon PNGs were last
committed in `45e9c39` (09-03) while `sid2_ribbon_render.py` was last changed
in `fc922f9` (09-04).

## 3. The "verify against UniProt" item is settled but still posted as open

`METHODS.txt:513-520` now records that McEwan et al. 2012 test exactly H32,
H168 and H175; that the modelled ectodomain contains exactly three histidines
at those positions, which confirms the numbering independently of any sequence
database; and that residue 199 came from the UniProt-derived list and sits in
the transmembrane helix rather than the ectodomain.

The instruction to check UniProt nevertheless survives in three places
(line numbers at `e6b1bbd`; they move as the report is re-knit, so grep for
"against UniProt" rather than trusting them):

- `FIGURE_REPORT.Rmd:2011-2013` — Figure 4 caveat 8
- `FIGURE_REPORT.Rmd:3481-3482` — the pre-submission list
- `FIGURE_CAPTIONS.txt:504-506` — the same caveat 8 in the caption file

The report's copy also says those residues "carry the proximity argument in
Figure 4C", which is no longer true: Figure 4C is the charge panel, and
`METHODS.txt` recommends dropping the proximity claim rather than defending it
at p = 0.38.

## 4. The curated figure count disagrees three ways

At `e6b1bbd`, `plots/` holds **23** curated figures. But `README.md` says
"twenty-one" in three places (lines 17, 53, 100) and its curated table lists
18 rows; `DATA_AVAILABILITY.md` says "twenty-one" twice (lines 5, 51);
`supplemental_data/SUPPLEMENTAL_DATA_OVERVIEW.md` says "eighteen" twice
(lines 3, 5). It was 22 on 09-09, so the count is also still moving.

This matters more than a stale count
usually would, because the deposit's self-containment claim is stated as a
*verified* count — so that sentence currently attests something narrower than
what is actually in the tree.

## 5. FIGURE_REPORT.md cannot attest the pooled peak coordinates

`scripts/check_manuscript_numbers.py` flags `X:4,875,969`, `III:12,353,680`,
`V:14,647,434` and `-log10p 7.71` in the Results draft. All four are correct —
`results_review/verify_claims.R` recomputes each from
`pos1_2023_gemma_loco.csv.gz` and `pooled_cross_bundle_thinned.rds` — but
`FIGURE_REPORT.md` never prints them, so the checker has nothing to match
against.

The fix is to surface the pooled peak coordinates in `FIGURE_REPORT.Rmd`, not
to add them to `manuscript_number_exceptions.txt`, which is reserved for
numbers the code does not produce.

## 6. Supplement numbering in the draft does not match the repository

The draft's supplement numbers are their own sequence. Apparent mapping, which
needs the author's confirmation: draft S1 = repo S1, S2 = S2, S3 = S5, S4 = S6,
S5 = S7, S6 = S8, S7 = S9, S8 = S11, S9 = S15.

## Concurrency note

Claude Code commits to this repository while Claude Science sessions are
running — HEAD moved from `86410ad` to `8c140d8` mid-session on 09-09, and
again to `e6b1bbd` overnight. Do not snapshot and restore files inside the
working tree to test a rebuild; copy the deposit to a scratch directory and run
there. Re-check `git log` before assuming the tree is where it was left.
