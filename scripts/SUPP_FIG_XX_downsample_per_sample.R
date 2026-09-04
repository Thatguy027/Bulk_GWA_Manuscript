## Supplement -- the depth requirement, per sample ---------------------------
##
##   Rscript scripts/SUPP_FIG_XX_downsample_per_sample.R
##     -> plots/SUPP_FIG_XX_downsample_per_sample.{pdf,png}
##
##   A  per-sample agreement against sequencing depth: one line per sample
##   B  the aggregate, slope-level curve for comparison -- panel C of
##      Figure1_boot.R, with the 0.5x depth restored
##
## Panels and the shared theme are in Figure1_common.R.
##
## Why the two levels differ
## -------------------------
## Panel B correlates SLOPES, fitted across three timepoints and averaged over
## five replicate arms, so it answers "how deep must we sequence before the
## per-strain rate of change is recovered". Panel A correlates individual
## sample frequencies, so it answers "how deep before any one measurement is
## right". A slope aggregates fifteen measurements, so B sits higher than A at
## every depth and saturates sooner -- averaging removes noise the individual
## samples still carry. Both are true; they are answers to different questions,
## and B alone would overstate how well a single sample is measured.
##
## Two things panel A shows that panel B cannot
## --------------------------------------------
## 1  The spread between samples is larger than the effect of depth. The
##    median moves about 0.10 across a fortyfold range of depth, while at any
##    fixed depth the samples span roughly 0.43 to 0.88.
## 2  The worst sample is worst at every depth, by a wide margin and without
##    improving. That is not a depth problem, and no amount of sequencing
##    would have fixed it -- see also
##    scripts/SUPP_FIG_XX_baugh_per_sample_frequencies.R, where the same
##    sample is the one visibly off the diagonal.
##
## The 0.5x depth is restored in both panels here. Figure1_boot.R drops it,
## inherited from the legacy script, and the aggregate is why: at the slope
## level 0.5x scores above the 1x that follows it, so the curve is not
## monotone. No such inversion appears per sample, so the wobble looks like
## noise in one downsampling draw rather than a property of that depth.
##
## No error bars: the bootstrap array holds replicates of the full-depth
## deconvolution only, so there is nothing to resample at a reduced depth. The
## dashed line in A is the full-depth median, which is the ceiling.
## ---------------------------------------------------------------------------

source("scripts/Figure1_common.R")

REFRESH <- nzchar(Sys.getenv("FIG1_REFRESH"))

msg("frequencies")
freq   <- baugh_frequencies(refresh = REFRESH)
slopes <- platform_slopes(freq)

msg("panels")

## A SHARED y AXIS. The panels exist to be compared, and on independent scales
## they cannot be: the aggregate spans 0.833-0.906 and would fill its panel
## with a curve that looks as dramatic as A's while covering a seventh of the
## range. On one scale the actual finding is visible -- B sits above A and is
## nearly flat, because averaging fifteen samples removes most of what depth
## costs an individual one.
YLIM <- c(0.40, 0.95)

## the panel subtitles are wrapped for a HALF-width panel. wrap_md()'s default
## is set for a full-page panel and the two subtitles overlapped each other.
pA <- panel_downsample_samples(freq, letter = "A", bare = FALSE) +
  labs(subtitle = wrap_md(paste0(
    "One grey line per sample, 23 samples; dark line is the median across ",
    "samples. Dashed line is the median full-depth agreement."),
    width = 62)) +
  coord_cartesian(ylim = YLIM)

## drop_depth = NULL restores the 0.5x that Figure1_boot.R panel C omits
pB <- panel_downsample(slopes, letter = "B", bare = FALSE,
                       drop_depth = NULL) +
  labs(subtitle = wrap_md(paste0(
    "Slopes recomputed at each depth, then correlated with the MIP-seq ",
    "slopes. Same y axis as A."), width = 62)) +
  coord_cartesian(ylim = YLIM)

fig <- (pA | pB) +
  plot_annotation(
    title = "**Sequencing depth: per sample and in aggregate**",
    subtitle = wrap_md(paste0(
      "Reads subsampled to each depth and the deconvolution rerun, then ",
      "compared with the published MIP-seq frequencies. A compares individual ",
      "sample frequencies; B compares the per-strain slopes that Figure 1 ",
      "plots, which average fifteen samples and so sit higher and saturate ",
      "sooner. Spread between samples exceeds the effect of depth, and the ",
      "worst sample does not improve with depth.")),
    theme = theme(plot.title = element_markdown(size = 13),
                  plot.subtitle = element_markdown(size = 8.4,
                                                   colour = "grey30"),
                  plot.title.position = "plot"))

write_fig(fig, "SUPP_FIG_XX_downsample_per_sample", width = 10.4, height = 5.2)
