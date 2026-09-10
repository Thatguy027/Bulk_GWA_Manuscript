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
## Panel C is the same five strains in one row each, JU1793 at the bottom to
## JU2466 at the top: genotype on the left, hatching on the right. The file
## keeps the name Figure3_quad because that stem is registered in
## FIGURE_REPORT.Rmd, FIGURE_CAPTIONS.txt and the figure-list check; it is
## three panels now, not four.
##
## CAVEAT: one plate per strain per condition -- see Figure3_common.R.
## ---------------------------------------------------------------------------

source("scripts/Figure3_common.R")

msg("panel A: pooled pos-1 phenotype")
## a histogram, not a kernel density: 84 strains are few enough to show as
## countable bins, and the density version invents small bumps in both tails
## that one strain each supports
pA <- pheno_inset(type = "hist", base_size = 11.5, letter = "A")

msg("panel B: parental allele frequency, chromosome III right arm")
freq <- load_parent_freq(chrom_keep = "III", from_mb = 8)
pB <- panel_parent_freq_chr3(freq, letter = "B")

msg("panel C: NIL genotypes and hatching, merged")
## Was two panels, C (hatching) and D (genotypes), drawing the same five
## strains as rows twice with two x axes. panel_nil_geno_hatch() puts each
## strain's genotype immediately left of its own hatching bar, so the pairing
## no longer has to be carried across a panel boundary. panel_C() and panel_B()
## are left in Figure3_common.R: Figure3_chrIII.R and the supplements still use
## them, and they are what the two-panel version was.
## labels = FALSE: no strain names on the rows. The row order is fixed and
## stated in the caption -- JU1793 at the bottom, then wSZ196, wSZ191,
## wSZ176, JU2466 at the top -- and Figure S10 shows the same series with
## the names on, so nothing is only knowable from this panel.
pC <- panel_nil_geno_hatch(letter = "C", labels = FALSE)

## C is now as wide as A and B together, so it takes the whole bottom row
fig <- (pA + pB) / pC + plot_layout(heights = c(1, 1.08))

ggsave(file.path(OUT, "Figure3_quad.pdf"), fig, width = 9.6, height = 6.2,
       device = cairo_pdf)
ggsave(file.path(OUT, "Figure3_quad.png"), fig, width = 9.6, height = 6.2,
       dpi = 300, bg = "white")
msg("wrote Figure3_quad.{pdf,png}")
