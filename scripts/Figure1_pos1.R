## Figure 1, pos-1 variant ---------------------------------------------------
## Validation, then the phenotype the validated assay produces, then the map.
##
##   Rscript scripts/Figure1_pos1.R  ->  plots/Figure1_pos1.{pdf,png}
##
##   A  MIP-seq slope against NNLS slope -- the assay agrees with the published
##      platform on the same samples
##   B  the 2023 pos-1 phenotype distribution across wild isotypes
##   C  the association scan for that phenotype
##
## Panels and the shared theme are in Figure1_common.R.
##
## B and C are panels A and B of SUPP_FIG_XX_2023_pos1_vst_and_mapping, moved
## into the main figure. The argument becomes one continuous claim: the
## deconvolution is accurate (A), so the phenotype it yields is a phenotype (B),
## and it maps (C). The supplement then only has to carry the replicate
## reproducibility and the depth work.
##
## The trait and the scan are used AS SHIPPED. scripts/2023_pos1_analysis.R
## documents the read-depth cutoff, the JU1793 duplicate-entry correction, and
## the caveat that the shipped vst still carries the averaged JU1793 value --
## which matters here because JU1793 is a cross parent in Figures 2 and 3.
##
## MARK names strains in panel B; NULL is the plain distribution, which is what
## is used. It was briefly set to the two cross parents, but only JU1793 draws
## -- JU2466 is not in this 2023 panel -- and a single marked strain reads as a
## claim about that strain rather than as context for Figures 2 and 3.
## ---------------------------------------------------------------------------

source("scripts/Figure1_common.R")

REFRESH <- nzchar(Sys.getenv("FIG1_REFRESH"))
MARK    <- NULL

msg("frequencies")
freq   <- baugh_frequencies(refresh = REFRESH)
slopes <- platform_slopes(freq)

msg("panels")
pA <- panel_slope(slopes, letter = "A")
pB <- panel_pos1_dist(letter = "B", mark = MARK)
pC <- panel_pos1_manhattan(letter = "C")

## A and B side by side on top, the genome-wide scan across the bottom, which
## is the only panel that needs the full page width.
##
## NESTED, deliberately. A flat plot_layout(design = ...) with C spanning both
## columns makes C wider than A and B together, by exactly A's y-axis gutter
## (103 px): the spanning panel keeps its own gutter and starts inside the
## column-1 gutter track. Nesting puts both gutters in the same outer track and
## the panel edges coincide -- verified at x = 117 and x = 1415 for both rows.
## See the note in Figure1_common.R; padding the gutter does not fix it.
fig <- (pA | pB) / pC + plot_layout(heights = c(1, 1))

write_fig(fig, "Figure1_pos1", width = 9.6, height = 7.0)
