## Supplement -- every sample that goes into the Figure 1 slopes -------------
##
##   Rscript scripts/SUPP_FIG_XX_baugh_per_sample_frequencies.R
##     -> plots/SUPP_FIG_XX_baugh_per_sample_frequencies.{pdf,png}
##
## Figure 1A plots one point per strain: a slope, fitted across three
## timepoints and averaged over five replicate arms. This unpacks it. Each
## facet is one sample, each point one strain, x the published MIP-seq
## frequency and y the NNLS frequency, with a bootstrap percentile interval on
## the NNLS estimate. These 15 samples are exactly the ones the slope fits use
## -- non-baseline, day 17 excluded -- so nothing here is aggregated.
##
## Panels and the shared theme are in Figure1_common.R.
##
## Error bars
## ----------
## Vertical only, and they are percentile intervals over the 100 bootstrap
## replicates of the deconvolution, taken at the sample level rather than
## propagated from a summary. There is no interval on the x axis: the published
## MIP-seq frequencies are point values with no uncertainty distributed
## alongside them, so drawing one would mean inventing it.
##
## Two things the aggregate figure cannot show
## -------------------------------------------
## 1  NNLS assigns an EXACT zero to a substantial share of strain-sample
##    combinations -- the non-negativity constraint clamping a strain out of a
##    pool -- while the MIP-seq table never reports zero. Those points sit on
##    y = 0 and are counted in the console output. They depress the per-sample
##    rank correlation in Figure 1B relative to the slope correlation in 1A,
##    because a run of ties at zero cannot be ranked.
## 2  Which samples carry the disagreement. Figure 1B reduces each of these
##    facets to a single rho and shows the distribution; here the sample with
##    the low rho is identifiable.
##
## The facet order is replicate then day, so a row is one arm's time course.
## ---------------------------------------------------------------------------

source("scripts/Figure1_common.R")

REFRESH <- nzchar(Sys.getenv("FIG1_REFRESH"))
LEVEL   <- 0.95

msg("frequencies")
freq <- baugh_frequencies(refresh = REFRESH)
bfi  <- boot_freq_intervals(refresh = REFRESH, level = LEVEL)

## the samples the slope fits actually use
d <- freq %>%
  filter(!baseline, day != 17, strain != "N2") %>%
  left_join(bfi, by = c("strain", "sample")) %>%
  filter(is.finite(frq), is.finite(published_frq))

excl <- freq %>% filter(baseline | day == 17) %>% distinct(sample) %>% pull(sample)
msg("  ", n_distinct(d$sample), " samples feeding the slopes | ",
    length(excl), " excluded: ", paste(excl, collapse = ", "))

n_zero <- sum(d$frq == 0)
msg("  ", nrow(d), " strain-sample points | NNLS exact zeros: ", n_zero,
    " (", round(100 * n_zero / nrow(d), 1), "%) | MIP zeros: ",
    sum(d$published_frq == 0))
msg("  median bootstrap interval width: ",
    signif(median(d$hi - d$lo), 3), " | median frequency: ",
    signif(median(d$frq), 3))

## per-sample Spearman, the quantity Figure 1B turns into a histogram
st <- d %>% group_by(sample) %>%
  summarise(rho = cor(frq, published_frq, method = "spearman"),
            n = n(), .groups = "drop") %>%
  mutate(lab = sprintf("rho = %.2f", rho))

cat("\n== per-sample agreement, slope-feeding samples only ==\n")
print(as.data.frame(st %>% transmute(sample, n, rho = round(rho, 3)) %>%
                      arrange(rho)), row.names = FALSE)

## facet order: replicate, then day, so each row is one arm's time course
ord <- d %>% distinct(sample, replicate, day) %>% arrange(replicate, day)
lv  <- ord$sample
lb  <- setNames(sprintf("%s · day %d", ord$replicate, ord$day), lv)
d  <- d  %>% mutate(sample = factor(sample, levels = lv))
st <- st %>% mutate(sample = factor(sample, levels = lv))

lim <- c(0, max(c(d$published_frq, d$hi), na.rm = TRUE))

fig <- ggplot(d, aes(published_frq, frq)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.35, colour = COL_FIT) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.3,
                colour = COL_PT, alpha = 0.5) +
  geom_point(size = 0.9, alpha = 0.6, colour = COL_PT) +
  geom_richtext(data = st, aes(x = -Inf, y = Inf, label = lab),
                inherit.aes = FALSE, hjust = -0.08, vjust = 1.4, size = 2.7,
                colour = COL_FIT, fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  facet_wrap(~sample, ncol = 3,
             labeller = labeller(sample = lb)) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "MIP-seq frequency", y = "NNLS frequency",
       title = titled(NULL, "**Sample-level frequencies behind the Figure 1 slopes**"),
       subtitle = wrap_md(sprintf(paste0(
         "One facet per sample, one point per wild isolate, for the %d samples ",
         "the slope fits use (non-baseline, day 17 excluded). Vertical bars are ",
         "%g%% percentile intervals over %d bootstrap replicates of the ",
         "deconvolution, taken at the sample level; the MIP-seq frequencies are ",
         "published point values and carry no interval. Dashed line is y = x, ",
         "not a fit. N2 excluded. %s of %s points (%s%%) are exact NNLS zeros, ",
         "where the non-negativity constraint drops a strain from the pool; the ",
         "MIP-seq table reports no zeros. Rows are one replicate arm's time ",
         "course."),
         n_distinct(d$sample), LEVEL * 100, 100,
         format(n_zero, big.mark = ","), format(nrow(d), big.mark = ","),
         round(100 * n_zero / nrow(d), 1)), width = 132)) +
  theme_pub(10.5) +
  theme(strip.text = element_text(face = "bold", size = 9.5),
        panel.spacing = grid::unit(6, "pt"),
        axis.text = element_text(size = 7.5),
        plot.title = element_markdown(size = 12),
        plot.subtitle = element_markdown(size = 8, colour = "grey30"))

write_fig(fig, "SUPP_FIG_XX_baugh_per_sample_frequencies",
          width = 8.6, height = 11)
