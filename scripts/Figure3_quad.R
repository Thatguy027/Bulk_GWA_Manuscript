## Figure 3, version 4 -- four panels, phenotype first -----------------------
##
##   Rscript scripts/Figure3_quad.R  ->  plots/Figure3_quad.{pdf,png}
##
##   A  the pooled pos-1 phenotype distribution, with both cross parents marked
##   B  the JU1793 x JU2466 HT115-vs-pos-1 cross scan on chromosome III
##   C  embryo hatching under pos-1 RNAi for the NIL series
##   D  the NIL introgressions on the right arm of chromosome III
##
## Everything shared with the other three variants is in Figure3_common.R.
##
## How this differs from the other three
## -------------------------------------
## The pooled phenotype is a panel in its own right rather than an inset, so
## the figure reads as the argument runs: these two strains sit at opposite
## ends of the pooled panel (A), crossing them maps a QTL (B), NILs carrying
## pieces of that interval give an allelic series of hatching (C), and the
## pieces they carry are these (D).
##
## The cost is that panel A takes real space instead of occupying dead space,
## so the chromosome III scan is narrower here than in Figure3_chrIII.R. The
## benefit is that nothing is hidden inside another panel, which matters if
## the figure is ever printed at column width.
##
## Layout: column 1 holds the two narrow panels (A, C), column 2 the two wide
## ones (B, D). The scan and the introgressions both want horizontal room; a
## distribution of 84 values and five percentage bars do not.
##
## The bottom row is C then D, phenotype then genotype, as specified. Rows in
## C and D are the same five strains in the same order, JU1793 at the bottom to
## JU2466 at the top, and share y limits so the rows line up across the pair.
##
## CAVEAT: one plate per strain per condition -- see Figure3_common.R.
## ---------------------------------------------------------------------------

source("scripts/Figure3_common.R")

msg("panel A: pooled pos-1 phenotype")
## a histogram, not a kernel density: 84 strains are few enough to show as
## countable bins, and the density version invents small bumps in both tails
## that one strain each supports
pA <- pheno_inset(type = "hist", base_size = 11.5, letter = "A")

msg("panel B: cross scan, chromosome III")
scan_t <- thin_scan(load_scan(), bin = 1e3)
pB <- panel_A_chr3(scan_t, letter = "B")

msg("panel C: hatching")
pC <- panel_C(letter = "C")

msg("panel D: NIL introgressions")
pD <- panel_B(letter = "D")

fig <- four_panel(pA, pB, pC, pD)

ggsave(file.path(OUT, "Figure3_quad.pdf"), fig, width = 9.6, height = 6.2,
       device = cairo_pdf)
ggsave(file.path(OUT, "Figure3_quad.png"), fig, width = 9.6, height = 6.2,
       dpi = 300, bg = "white")
msg("wrote Figure3_quad.{pdf,png}")
