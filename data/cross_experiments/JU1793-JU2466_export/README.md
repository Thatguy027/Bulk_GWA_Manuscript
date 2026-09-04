# JU1793 x JU2466 — counts, plotting data and figures

Exported from `Nov2024_JU_cross_pos_mig` on 2026-09-02. Everything here concerns the JU1793 x JU2466 cross only.

- **4 samples** (S31, S32, S33, S34), conditions: HT115g, MIG6g, OP50, POS1g
- **3 contrasts** between the samples listed in `contrast_to`; S31, S34 carry no contrast and are included for their allele frequencies only
- allele-frequency deviations are phased to **JU1793**, so a positive `contrast.beta` means the JU1793 allele is at higher frequency on the left side of the contrast

## Layout

```
sample_sheet.tsv              the samples and which contrasts were run
counts/
  aser/<s>.table.gz           GATK ASEReadCounter counts
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
`p1` and `p2` are the counts assigned to the JU1793 and JU2466 haplotypes; `ref`/`alt` are the raw reference and alternate counts before phasing.

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
results and an `af_counts.rda`, so the shipped copies are left alone.
Budget roughly 2 GB and ten minutes.

This was checked: regenerating from `counts/aser/` reproduces every peak position and
interval bound in `intervals/` exactly, with LOD values matching to floating-point
accumulation order.

Requires R with
`xQTLStats`, `tidyverse`, `vcfR`, `AlphaSimR`, `Rfast`, `ggpubr`, `data.table`,
`yaml`, `qs2`, and the reference objects named in `scripts/config.yml`
(`/Users/Stefan/github_repos/xQTLSims/`).

No report-generation code is included.
