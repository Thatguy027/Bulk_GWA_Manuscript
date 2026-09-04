## Figure 3, version 1 -- genome-wide cross scan in panel A ------------------
##
##   Rscript scripts/Figure3.R      ->  plots/Figure3.{pdf,png}
##
##   A  the JU1793 x JU2466 HT115-vs-pos-1 cross scan, all six chromosomes
##   B  the NIL introgressions on the right arm of chromosome III
##   C  embryo hatching under pos-1 RNAi for the same strains
##
## Everything shared with the other two variants is in Figure3_common.R; see
## that file's header for the data sources, the LOD-vs--log10p reasoning and
## the panel C caveat.
##
## Two things specific to this version
## -----------------------------------
## 1  No marker on the mapped region in panel A. The region is 37 kb, which is
##    a fraction of a pixel at genome scale; drawn at that size it reads as
##    noise or as part of the chrIII peak rather than as an annotation. Panel B
##    shades it, and scripts/Figure3_chrIII.R marks it on an axis where 37 kb
##    is large enough to place honestly.
##
## 2  The whole figure is deliberately compact. The chrIII peak is LOD 140 and
##    the next highest chromosome reaches 57, so a panel A tall enough to be
##    square is mostly empty air; the height ratio below and the small upper y
##    expansion pull that back without shortening the peak. The canvas is also
##    smaller than a full page, which shrinks the panel B and C bars -- they
##    were larger than anything read off them requires -- and, because the
##    type size is fixed, tightens the empty regions relative to the ink.
## ---------------------------------------------------------------------------

source("scripts/Figure3_common.R")

msg("panel A: cross scan, genome-wide")
scan_t <- thin_scan(load_scan())
pA <- panel_A_genome(scan_t)

msg("panel B: NIL introgressions")
pB <- panel_B()

msg("panel C: hatching")
pC <- panel_C()

fig <- three_panel(pA, pB, pC, heights = c(0.72, 1))

ggsave(file.path(OUT, "Figure3.pdf"), fig, width = 9.6, height = 5.9,
       device = cairo_pdf)
ggsave(file.path(OUT, "Figure3.png"), fig, width = 9.6, height = 5.9,
       dpi = 300, bg = "white")
msg("wrote Figure3.{pdf,png}")
