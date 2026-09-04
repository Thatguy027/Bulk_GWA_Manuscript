## Figure 3, version 2 -- chromosome III only in panel A ---------------------
##
##   Rscript scripts/Figure3_chrIII.R  ->  plots/Figure3_chrIII.{pdf,png}
##
##   A  the JU1793 x JU2466 HT115-vs-pos-1 cross scan on chromosome III, with
##      the pooled pos-1 phenotype distribution inset over the left arm, where
##      the scan is flat
##   B  the NIL introgressions on the right arm of chromosome III
##   C  embryo hatching under pos-1 RNAi for the same strains
##
## Everything shared with the other two variants is in Figure3_common.R.
##
## Why this version exists
## -----------------------
## All three panels are then about one chromosome, so the figure reads as a
## single continuous zoom: chromosome III, then its right arm, then the
## phenotype. The genome-wide context is lost, which is the trade -- a reader
## cannot see from this panel that chromosome III is the only strong signal.
## That is what version 1 is for.
##
## Because the axis is 13.8 Mb rather than 98 Mb, the 37 kb mapped region is
## 0.3% of the panel width and can be drawn as a band with dotted edges, which
## is not possible at genome scale.
##
## The QTL is at the right telomere and the trace is flat across 12 of 13.8 Mb,
## so the left two thirds of the panel is empty. The pooled phenotype inset
## goes there: both cross parents were phenotyped in the pooled panel and sit
## at opposite ends of it, which is why they were chosen, and putting that in
## the same frame connects the pooled experiment, the cross and the NILs. It is
## drawn as a histogram here rather than a kernel density: 84 strains are few
## enough to show as countable bins, and the density version invents small
## bumps in both tails that no strain supports. See Figure3_common.R for the
## trait and the equal-height lollipop reasoning.
## ---------------------------------------------------------------------------

source("scripts/Figure3_common.R")

msg("panel A: cross scan, chromosome III")
## finer bins than the genome-wide version: one chromosome across the same
## page width is ~7x the resolution, so 1 kb bins are still ~1 point per pixel
scan_t <- thin_scan(load_scan(), bin = 1e3)
pA  <- panel_A_chr3(scan_t)
pIn <- pheno_inset(type = "hist")

## left/right in panel-relative units, so 0.035-0.530 is 0.5-7.3 Mb -- entirely
## within the flat left arm, clear of the rise that starts near 10 Mb
pA_in <- pA + inset_element(pIn, left = 0.035, bottom = 0.13,
                            right = 0.530, top = 0.98, align_to = "panel")

msg("panel B: NIL introgressions")
pB <- panel_B()

msg("panel C: hatching")
pC <- panel_C()

fig <- three_panel(pA_in, pB, pC, heights = c(0.8, 1))

ggsave(file.path(OUT, "Figure3_chrIII.pdf"), fig, width = 9.6, height = 6.4,
       device = cairo_pdf)
ggsave(file.path(OUT, "Figure3_chrIII.png"), fig, width = 9.6, height = 6.4,
       dpi = 300, bg = "white")
msg("wrote Figure3_chrIII.{pdf,png}")
