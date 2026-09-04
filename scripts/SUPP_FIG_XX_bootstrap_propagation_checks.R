## Supplement -- does the bootstrap propagate correctly? ---------------------
##
##   Rscript scripts/SUPP_FIG_XX_bootstrap_propagation_checks.R
##     -> plots/SUPP_FIG_XX_bootstrap_propagation_checks.{pdf,png}
##
## The intervals in Figure1_boot.R are obtained by redoing the whole slope
## calculation inside each of the 100 bootstrap replicates of the
## deconvolution -- delta against the day-1 sample of the same replicate arm,
## regression on day, average over arms -- rather than by pushing the stored
## bootse through that chain. This figure is the audit of that claim. Every
## panel is a check that could fail, drawn so it can be seen to have passed.
##
## Panels and the shared theme are in Figure1_common.R.
##
## Only panel F uses coord_equal, because only there is the comparison against
## y = x the actual claim. In A and B the correlation is 0.999984 and 0.999994,
## so nothing is being judged by eye off the diagonal and an equal aspect only
## left dead height in the row.
##
##   A  frequency level. Does `ab` resample the same estimator that Figure 1
##      plots? If the bootstrap were centred on some other quantity, or
##      normalised differently, this would not lie on y = x.
##   B  slope level. Same question after the transformation. This is the panel
##      that would break if the day-1 baseline were taken from the point
##      estimate instead of from within each replicate, or if the arms were
##      averaged in the wrong order.
##   C  standardised bias, (point estimate - bootstrap mean) / bootstrap SD.
##      A transformation applied inconsistently shows up here as a shift away
##      from zero measured in units the intervals themselves set.
##   D  coverage. Each strain's interval with its point estimate marked. An
##      interval that does not contain its own point estimate is misplaced.
##   E  the correlation's own bootstrap distribution, which is what
##      Figure1_boot.R reports in place of a p value against rho = 0.
##   F  scale. Bootstrap interval width against distance from y = x, per
##      strain. This is the panel that explains why the error bars in
##      Figure1_boot.R are invisible: the uncertainty is genuinely smaller
##      than the platform disagreement, so most points fall below the line.
##
## CAVEAT, and it is not resolvable from this repository: only the bootstrap
## OUTPUTS were saved. data/baugh/2024baugh_bootstrap_prediction.rda holds `ab`
## (102 strains x 23 samples x 100 replicates) and `blist`, which is the same
## content reshaped -- not the resampling scheme. So these panels establish
## that the array is resampling the plotted estimator and that the
## transformation is carried through it consistently. They cannot establish
## WHAT was resampled -- reads, markers or strains -- because the generating
## code is not here. Any statement stronger than "100 bootstrap replicates of
## the deconvolution" is not supported by what is on disk.
## ---------------------------------------------------------------------------

source("scripts/Figure1_common.R")

REFRESH <- nzchar(Sys.getenv("FIG1_REFRESH"))

msg("inputs")
freq   <- baugh_frequencies(refresh = REFRESH)
slopes <- platform_slopes(freq)
S      <- boot_slope_matrix(refresh = REFRESH)
bs     <- boot_slope_summary(S)

## --- A: frequency level ----------------------------------------------------
fq <- freq %>% filter(is.finite(bootmean), is.finite(frq))
r_fq <- cor(fq$bootmean, fq$frq)
msg("  A frequency level: n = ", nrow(fq), " r = ", signif(r_fq, 6),
    " median diff = ", signif(median(fq$bootmean - fq$frq), 3))

pA <- ggplot(fq, aes(frq, bootmean)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.35, colour = COL_FIT) +
  geom_point(size = 0.7, alpha = 0.35, colour = COL_PT) +
  annotate("richtext", x = -Inf, y = Inf,
           label = sprintf("r = %.6f<br>n = %s", r_fq,
                           format(nrow(fq), big.mark = ",")),
           hjust = -0.08, vjust = 1.25, size = 2.7, colour = COL_FIT,
           fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt")) +
  labs(x = "Point-estimate frequency", y = "Bootstrap mean frequency",
       title = titled("A", "**Same estimator, frequency level**")) +
  theme_pub(10.5)

## --- B: slope level --------------------------------------------------------
pt <- slopes %>% filter(strain != "N2") %>% group_by(strain) %>%
  summarise(wgs = mean(wgs_slope, na.rm = TRUE), .groups = "drop")
j <- inner_join(pt, bs, by = "strain") %>% filter(is.finite(wgs))
r_sl <- cor(j$wgs, j$boot_mean)
cover <- sum(j$wgs >= j$lo & j$wgs <= j$hi)
msg("  B slope level: n = ", nrow(j), " r = ", signif(r_sl, 6),
    " | coverage ", cover, "/", nrow(j))

pB <- ggplot(j, aes(wgs, boot_mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.35, colour = COL_FIT) +
  geom_point(size = 1.4, alpha = 0.6, colour = COL_PT) +
  annotate("richtext", x = -Inf, y = Inf,
           label = sprintf("r = %.6f<br>n = %d", r_sl, nrow(j)),
           hjust = -0.08, vjust = 1.25, size = 2.7, colour = COL_FIT,
           fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt")) +
  labs(x = "Point-estimate slope", y = "Bootstrap mean slope",
       title = titled("B", "**Same estimator, after the transformation**")) +
  theme_pub(10.5)

## --- C: standardised bias --------------------------------------------------
j <- j %>% mutate(z = (wgs - boot_mean) / boot_sd)
msg("  C standardised bias: median ", signif(median(j$z), 3),
    " | mean ", signif(mean(j$z), 3), " | max |z| ", signif(max(abs(j$z)), 3))

pC <- ggplot(j, aes(z)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = COL_FIT) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dotted",
             linewidth = 0.35, colour = "grey55") +
  geom_histogram(bins = 26, fill = COL_HIST, colour = NA) +
  annotate("richtext", x = Inf, y = Inf,
           label = sprintf("median = %.3f<br>largest = %.2f SD",
                           median(j$z), max(abs(j$z))),
           hjust = 1.05, vjust = 1.25, size = 2.7, colour = COL_FIT,
           fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt")) +
  labs(x = "(point &minus; bootstrap mean) / bootstrap SD", y = "Strains",
       title = titled("C", "**Bias in units of the interval**")) +
  theme_pub(10.5) +
  theme(axis.title.x = element_markdown())

## --- D: coverage caterpillar ----------------------------------------------
cat_d <- j %>% arrange(wgs) %>% mutate(rank = row_number(),
                                       inside = wgs >= lo & wgs <= hi)
pD <- ggplot(cat_d, aes(rank, wgs)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.4,
                 colour = COL_PT, alpha = 0.5) +
  geom_point(aes(colour = inside), size = 0.7) +
  scale_colour_manual(values = c(`TRUE` = COL_PT, `FALSE` = COL_FIT),
                      guide = "none") +
  labs(x = "Strain, ordered by slope", y = "NNLS slope",
       title = titled("D", sprintf("**Coverage: %d of %d intervals contain their point estimate**",
                                   cover, nrow(j)))) +
  theme_pub(10.5)

## --- E: the correlation's bootstrap distribution --------------------------
d98 <- slopes %>% filter(strain != "N2") %>% group_by(strain) %>%
  summarise(mip = mean(mip_slope, na.rm = TRUE),
            wgs = mean(wgs_slope, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(mip), is.finite(wgs)) %>%
  left_join(bs, by = "strain")
br <- boot_rho(S, setNames(d98$mip, d98$strain))
r_full <- cor(d98$mip, d98$wgs, method = "spearman")
msg("  E bootstrap rho: median ", signif(br$med, 4),
    " [", signif(br$lo, 4), ", ", signif(br$hi, 4),
    "] | full-data rho ", signif(r_full, 4))

pE <- ggplot(tibble(rho = br$rho), aes(rho)) +
  geom_histogram(bins = 24, fill = COL_HIST, colour = NA) +
  geom_vline(xintercept = c(br$lo, br$hi), linetype = "dotted",
             linewidth = 0.4, colour = "grey45") +
  geom_vline(xintercept = r_full, linewidth = 0.5, colour = COL_FIT) +
  annotate("richtext", x = -Inf, y = Inf,
           label = sprintf("median %.3f<br>95%% [%.3f, %.3f]<br>full data %.3f",
                           br$med, br$lo, br$hi, r_full),
           hjust = -0.08, vjust = 1.25, size = 2.7, colour = COL_FIT,
           fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt")) +
  labs(x = "Spearman's &rho; vs MIP-seq, per bootstrap replicate",
       y = "Replicates",
       title = titled("E", "**The correlation is stable under resampling**")) +
  theme_pub(10.5) +
  theme(axis.title.x = element_markdown())

## --- F: uncertainty against platform disagreement -------------------------
d98 <- d98 %>% mutate(width = hi - lo, resid = abs(wgs - mip))
below <- sum(d98$width < d98$resid, na.rm = TRUE)
msg("  F scale: median width ", signif(median(d98$width), 3),
    " | RMS deviation from y=x ", signif(sqrt(mean((d98$wgs - d98$mip)^2)), 3),
    " | width < |residual| for ", below, "/", nrow(d98))

pF <- ggplot(d98, aes(resid, width)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.35, colour = COL_FIT) +
  geom_point(size = 1.4, alpha = 0.6, colour = COL_PT) +
  annotate("richtext", x = Inf, y = -Inf,
           label = sprintf(paste0("below the line: %d of %d<br>",
                                  "uncertainty smaller than disagreement"),
                           below, nrow(d98)),
           hjust = 1.05, vjust = -0.45, size = 2.7, colour = COL_FIT,
           fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt")) +
  coord_equal() +
  labs(x = "Distance from *y* = *x* (absolute slope difference)",
       ## short: a longer rotated title reaches into the title row and
       ## overlaps the panel letter, which sits at the plot's left edge
       y = "95% interval width",
       title = titled("F", "**Why the error bars are invisible**")) +
  theme_pub(10.5) +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown())

## ===========================================================================
fig <- (pA + pB) / (pC + pD) / (pE + pF) +
  plot_annotation(
    title = "**Bootstrap propagation checks for Figure 1**",
    subtitle = wrap_md(paste0(
      "The slope intervals in Figure1_boot.R are obtained by redoing the ",
      "whole slope calculation inside each of the 100 bootstrap replicates of ",
      "the deconvolution, not by propagating a stored standard error through ",
      "it. A-D check that the array resamples the estimator the figure plots ",
      "and that the transformation is carried through it consistently; E gives ",
      "the correlation its own interval; F shows that the resulting ",
      "uncertainty is smaller than the disagreement between platforms, which ",
      "is why the bars cannot be seen at full scale. Only the bootstrap ",
      "outputs were saved, so what was resampled is not recoverable here.")),
    theme = theme(plot.title = element_markdown(size = 13),
                  plot.subtitle = element_markdown(size = 8.2,
                                                   colour = "grey30"),
                  plot.title.position = "plot"))

write_fig(fig, "SUPP_FIG_XX_bootstrap_propagation_checks",
          width = 9.4, height = 10.0)
