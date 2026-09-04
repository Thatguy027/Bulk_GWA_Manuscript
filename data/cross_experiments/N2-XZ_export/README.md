# N2 x XZ1516 — counts, plotting data and figures

Exported from `NJX_rnai` on 2026-09-02. Everything here concerns the N2 x XZ1516 cross only.

- **10 samples** (S11, S17, S20, S23, S29, S32, S35, S38, S41, S44), conditions: HT115, mig6, par1, pos1, rpn12, vha5
- **15 contrasts**, all pairwise among the timepoint-2 samples
- allele-frequency deviations are phased to **N2**, so a positive `contrast.beta` means the N2 allele is at higher frequency on the left side of the contrast

## Layout

```
sample_sheet.tsv              the samples and which contrasts were run
counts/
  aser_combined/<s>.table.gz  GATK ASEReadCounter counts, sequencing runs summed
  phased/<label>.counts.tsv.gz  the same counts phased onto p1/p2 at genetic-map markers
afd/<label>.afd.tsv.gz        per-sample allele frequency deviation and its standard error
plot_data/<contrast>_plot_DF.tsv.gz   full per-marker statistics behind every figure
intervals/<contrast>_intervals.tsv    support intervals per chromosome
plots/
  per_sample/<label>_cuts.png         allele frequency across the genome, one per sample
  contrast/<contrast>_LOD.{pdf,png}   clean LOD trace with the support interval shaded
scripts/                      everything needed to regenerate the above
MANIFEST.tsv                  every file with its size
```

## Columns

`counts/phased/*.counts.tsv.gz` — `id`, `chr`, `pos`, `ref`, `alt`, `p1`, `p2`.
`p1` and `p2` are the counts assigned to the N2 and XZ1516 haplotypes; `ref`/`alt` are the raw reference and alternate counts before phasing.

`afd/*.afd.tsv.gz` — `chrom`, `physical.position`, `ID`, then per-sample
`ref_`, `alt_`, `p1_`, `p2_`, `nindv_`, `ndepth_`, `n_`, `afd_`, `afd.se_`.

`plot_data/*_plot_DF.tsv.gz` — the two samples' AFD columns side by side, then
`contrast.beta`, `contrast.beta.se`, `z`, `p`, `log10p`, `neglog10p`, `LOD`.

Use `neglog10p` rather than `p` for anything quantitative: `p` is computed on the
linear scale and underflows to exactly 0 above |z| ~ 38.5, which happens often here.
`LOD` is computed in log space and does not saturate.

## Support intervals

`lod.drop.used` in each interval table is the LOD drop that defined it. The drop is
10% of each chromosome's own peak, calibrated so that a known causal gene falls inside
the interval rather than just outside it. It is proportional rather than a fixed
number of LOD units because an absolute drop large enough for a strong peak spans
an entire chromosome on a weak one.

## Regenerating

```sh
bash scripts/run_all.sh
```

Outputs go to `regenerated/`, alongside a `cache/` of intermediate per-sample
results and an `af_counts.rda`, so the shipped copies are left alone. Expect
roughly 3 GB and half an hour.

This was checked: regenerating from `counts/aser_combined/` reproduces every peak
position and interval bound **exactly**, with LOD values agreeing to 1e-14 relative
(floating-point accumulation order, not a real difference).

Requires R with
`xQTLStats`, `tidyverse`, `vcfR`, `AlphaSimR`, `Rfast`, `ggpubr`, `data.table`,
`yaml`, `qs2`, and the reference objects named in `scripts/config.yml`
(`/Users/Stefan/github_repos/xQTLSims/`).

`scripts/NJX_combine_aser.R` produced `counts/aser_combined/` by summing the two
sequencing runs per library. The split per-run tables are not bundled (7.5 GB in the
source experiment); the combined tables here are its output.

No report-generation code is included.
