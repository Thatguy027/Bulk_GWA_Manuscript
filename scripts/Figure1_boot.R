## Figure 1, bootstrap variant -----------------------------------------------
## The platform-agreement panel with uncertainty on the NNLS slopes.
##
##   Rscript scripts/Figure1_boot.R  ->  plots/Figure1_boot.{pdf,png}
##
##   A  MIP-seq slope against NNLS slope, with bootstrap intervals on the NNLS
##      slope and a bootstrap interval on the correlation itself
##   B  the distribution of per-sample agreement
##   C  how much sequencing depth the deconvolution needs
##
## Identical to Figure1.R apart from panel A. See Figure1_common.R for how the
## intervals are computed -- briefly, the whole slope calculation is redone
## inside each of the 100 bootstrap replicates of the deconvolution, so the
## interval is on the slope rather than on a frequency, and the reported rho
## carries the percentile interval of rho across those replicates instead of a
## p value against rho = 0.
##
## First run opens the 30 MB genotype matrix once, to recover the strain order
## of the bootstrap array, and caches the result to
## data/baugh/cache_boot_slopes.rds. Set FIG1_REFRESH=1 to rebuild both caches.
## ---------------------------------------------------------------------------

source("scripts/Figure1_common.R")

REFRESH <- nzchar(Sys.getenv("FIG1_REFRESH"))

msg("frequencies")
freq   <- baugh_frequencies(refresh = REFRESH)
slopes <- platform_slopes(freq)

msg("bootstrap slopes")
S <- boot_slope_matrix(refresh = REFRESH)

msg("panels")
pA <- panel_slope(slopes,      letter = "A", boot_mat = S)
pB <- panel_sample_rho(freq,   letter = "B")
pC <- panel_downsample(slopes, letter = "C")

## one row of three; see Figure1.R for why panel A is not double width
fig <- pA + pB + pC + plot_layout(nrow = 1)

write_fig(fig, "Figure1_boot", width = 11, height = 3.9)
