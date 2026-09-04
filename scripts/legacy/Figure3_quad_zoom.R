## Figure 3, version 5 -- four panels, chromosome III cut at 8 Mb ------------
##
##   Rscript scripts/Figure3_quad_zoom.R  ->  plots/Figure3_quad_zoom.{pdf,png}
##
##   A  the pooled pos-1 phenotype distribution, with both cross parents marked
##   B  the cross scan over chromosome III from 8 Mb to the right telomere
##   C  embryo hatching under pos-1 RNAi for the NIL series
##   D  the NIL introgressions on the right arm of chromosome III
##
## Everything shared with the other variants is in Figure3_common.R.
##
## The cut
## -------
## Figure3_quad.R draws all 13.8 Mb of chromosome III, of which the first 10
## are flat: the whole left arm buys empty space. Starting the axis at 8 Mb
## spends that space on the rise instead, and each panel then covers roughly an
## order of magnitude less sequence than the one before it -- 5.8 Mb in B,
## 184 kb in D, 37 kb resolved -- so the figure reads as a zoom series.
##
## The cut is honest for this QTL because the trace is at LOD < 2 across the
## whole excluded interval, but it is a cut: a reader cannot confirm from panel
## B alone that nothing else on chromosome III is significant. Figure3_quad.R
## keeps the full chromosome and Figure3.R keeps the whole genome, for that.
## The line is clipped rather than filtered, so it enters the panel at its true
## height instead of appearing to start from zero at 8 Mb.
##
## Layout and the panel C / D row are as in Figure3_quad.R.
##
## CAVEAT: one plate per strain per condition -- see Figure3_common.R.
## ---------------------------------------------------------------------------

source("scripts/Figure3_common.R")

## SUPERSEDED. The curated figure set is listed in README.md; this variant was
## kept for comparison but is no longer a manuscript figure, so it writes to
## plots/legacy/ rather than cluttering plots/ on every run.
OUT <- "plots/legacy"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)


CUT <- 8   # Mb

msg("panel A: pooled pos-1 phenotype")
pA <- pheno_inset(type = "hist", base_size = 11.5, letter = "A")

msg("panel B: cross scan, chromosome III from ", CUT, " Mb")
scan_t <- thin_scan(load_scan(), bin = 1e3)
pB <- panel_A_chr3(scan_t, letter = "B", from_mb = CUT)

## what the cut hides, for the record and for the caption
excl <- scan_t %>% filter(chrom == "III", pos.mb < CUT)
msg("  excluded 0-", CUT, " Mb: max LOD ", round(max(excl$LOD), 2),
    " over ", nrow(excl), " bins")

msg("panel C: hatching")
pC <- panel_C(letter = "C")

msg("panel D: NIL introgressions")
pD <- panel_B(letter = "D")

fig <- four_panel(pA, pB, pC, pD)

ggsave(file.path(OUT, "Figure3_quad_zoom.pdf"), fig, width = 9.6, height = 6.0,
       device = cairo_pdf)
ggsave(file.path(OUT, "Figure3_quad_zoom.png"), fig, width = 9.6, height = 6.0,
       dpi = 300, bg = "white")
msg("wrote Figure3_quad_zoom.{pdf,png}")
