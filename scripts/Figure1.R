## Figure 1 -- the pooled assay measures what it claims to -------------------
##
##   Rscript scripts/Figure1.R      ->  plots/Figure1.{pdf,png}
##
##   A  MIP-seq slope against NNLS slope, one point per wild isolate
##   B  the distribution of per-sample agreement
##   C  how much sequencing depth the deconvolution needs
##
## Three panels in one row.
##
## Panels and the shared theme are in Figure1_common.R; see that file's header
## for what changed from scripts/legacy/Figure1_legacy_*.R and why.
##
## This is the arrangement the legacy script called figure1_wcor -- the most
## complete of the three it produced. The other two were subsets: the committed
## figure1.pdf was A over C, and the uncommitted edit had reduced it to A alone
## with C moved to a supplement. Nothing here is new data; the panels are the
## same three plots with real letters, a shared theme and a cached input.
##
## Set FIG1_REFRESH=1 in the environment to redo the NNLS from the genotype
## matrix instead of reading the committed result.
## ---------------------------------------------------------------------------

source("scripts/Figure1_common.R")

REFRESH <- nzchar(Sys.getenv("FIG1_REFRESH"))

msg("frequencies")
freq   <- baugh_frequencies(refresh = REFRESH)
slopes <- platform_slopes(freq)

msg("panels")
pA <- panel_slope(slopes,      letter = "A")
pB <- panel_sample_rho(freq,   letter = "B")
pC <- panel_downsample(slopes, letter = "C")

## One row of three. Panel A uses coord_equal -- the two axes are the same
## quantity in the same units, so an unequal aspect would misrepresent the
## scatter about y = x -- which means a double-width panel A would sit as a
## square in a wide box with dead space on both sides. Three equal columns give
## the square panel a square box.
fig <- pA + pB + pC + plot_layout(nrow = 1)

write_fig(fig, "Figure1", width = 11, height = 3.9)
