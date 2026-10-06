## Supplement -- the Baugh L1 comparisons, and why they give different numbers --
##
##   Rscript scripts/SUPP_FIG_XX_baugh_analyses.R
##     -> plots/SUPP_FIG_XX_baugh_analyses.{pdf,png}
##
## Figure N1 of SUPPLEMENTARY_NOTE_BAUGH.md, numbered within the note and
## outside the supplementary series. The Baugh data set is scored
## several ways in this repository, and two of those numbers -- 0.974 and 0.89
## -- reach the manuscript text side by side. They answer different questions:
##
##   A  PLATFORM agreement. Pooled WGS frequencies against MIP-seq frequencies
##      for the same samples, both turned into the SAME trait (the difference
##      slope). Only the frequency column differs. This is Figure 1A, rebuilt
##      by the same function.
##   B  TRAIT agreement. Our pooled WGS trait against the PUBLISHED eLife Slope,
##      our trait built on the published recipe (log2 f / baseline, floored at
##      1/(4n), day 17 excluded).
##   C  As B, with our trait built as the difference slope of panel A.
##   D  Every comparison on one axis, so the gap between A and B can be
##      attributed: MIP-seq's own difference slope scores no better against the
##      published Slope than ours does, so the drop from A to B is the change of
##      trait definition, not the deconvolution.
##
## All on one frequency set: the cached dep102 deconvolution Figure 1 reads,
## N2 excluded, n = 98. The dep103 mapping table (baugh_mapping_traits.csv)
## scores 0.884 at n = 99 on the same recipe; that is a different reference
## panel, not a different method.
##
## Every value the note quotes is printed below and pinned. Depth is not
## repeated here; that is Figure S6 (SUPP_FIG_XX_baugh_downsample_traits.R).
##
## Panels and the shared theme are in Figure1_common.R.
## ---------------------------------------------------------------------------

source("scripts/Figure1_common.R")

REFRESH <- nzchar(Sys.getenv("FIG1_REFRESH"))
PUB     <- file.path(POS1, "baugh_published_traits.txt")

msg("frequencies")
freq   <- baugh_frequencies(refresh = REFRESH)
slopes <- platform_slopes(freq)
pub    <- read_tsv(PUB, show_col_types = FALSE)

n_str  <- n_distinct(freq$strain)
FLOOR  <- 1 / (4 * n_str)
## the conventional alternative: half the smallest positive value
HALFMIN <- min(freq$frq[freq$frq > 0]) / 2

wgs   <- published_recipe_traits(freq, "frq", floor = FLOOR)
wgs_h <- published_recipe_traits(freq, "frq", floor = HALFMIN)
mip   <- published_recipe_traits(freq, "published_frq")

score <- function(t) {
  j <- pub %>% inner_join(t, by = "strain") %>% filter(strain != "N2")
  sp <- function(a, b) abs(cor(a, b, method = "spearman", use = "complete.obs"))
  list(j = j, n = nrow(j), cols = t$n_cols[1],
       slope = sp(j$Slope, j$slope), pc1 = sp(j$PC1, j$pc1),
       delta = sp(j$Slope, j$delta_slope))
}
sw <- score(wgs); sh <- score(wgs_h); sm <- score(mip)

## platform agreement, exactly as panel_slope() computes it
plat <- slopes %>% filter(strain != "N2") %>% group_by(strain) %>%
  summarise(mip = mean(mip_slope, na.rm = TRUE),
            wgs = mean(wgs_slope, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(mip), is.finite(wgs))
rho_plat <- cor(plat$mip, plat$wgs, method = "spearman")
ps <- full_depth_per_sample(freq)

n_zero <- sum(freq$frq == 0)

V <- tibble::tribble(
  ~key,              ~value,
  "platform_slope",  rho_plat,
  "per_sample_med",  median(ps$rho),
  "wgs_recipe_slope", sw$slope,
  "wgs_recipe_pc1",  sw$pc1,
  "wgs_delta_slope", sw$delta,
  "mip_delta_slope", sm$delta,
  "mip_recipe_slope", sm$slope,
  "mip_recipe_pc1",  sm$pc1,
  "halfmin_slope",   sh$slope,
  "halfmin_pc1",     sh$pc1)

cat("\n== Baugh L1: every comparison, one frequency set ==\n")
cat(sprintf("  strains %d (N2 excluded) | samples %d | recipe columns %d\n",
            sw$n, nrow(ps), sw$cols))
cat(sprintf("  floor 1/(4n) = %.5f | half-min floor = %.2e | NNLS exact zeros %d of %d\n",
            FLOOR, HALFMIN, n_zero, nrow(freq)))
print(as.data.frame(V %>% mutate(value = round(value, 3))), row.names = FALSE)

## ---- pins: every number SUPPLEMENTARY_NOTE_BAUGH.md quotes ----------------
PINNED <- c(platform_slope = 0.974, per_sample_med = 0.835,
            wgs_recipe_slope = 0.890, wgs_recipe_pc1 = 0.822,
            wgs_delta_slope = 0.888, mip_delta_slope = 0.890,
            mip_recipe_slope = 0.984, mip_recipe_pc1 = 0.961,
            halfmin_slope = 0.663, halfmin_pc1 = 0.638)
got <- setNames(V$value, V$key)[names(PINNED)]
off <- abs(got - PINNED) > 0.0005
if (any(off))
  stop("pinned values moved: ",
       paste(sprintf("%s %.4f (pinned %.3f)", names(PINNED)[off], got[off],
                     PINNED[off]), collapse = "; "))
stopifnot(sw$n == 98, sm$n == 98, sw$cols == 15, nrow(ps) == 23,
          n_zero == 331)
msg("  pins agree")

## ---- panels ----------------------------------------------------------------
pA <- panel_slope(slopes, letter = "A", bare = FALSE) +
  labs(title = titled("A", "**Platform agreement**"),
       subtitle = wrap_md(paste0(
         "Same trait, two platforms: the difference slope (f &minus; f<sub>day 1</sub> ",
         "on day) from pooled WGS against MIP-seq. As Figure 1A."), 62))

trait_panel <- function(y, ylab, rho, letter, ttl, sub) {
  lab <- sprintf("rho = %.3f, n = %d", rho, sw$n)
  ggplot(sw$j, aes(Slope, .data[[y]])) +
    geom_point(size = 1.9, alpha = 0.65, colour = COL_PT) +
    geom_richtext(data = tibble(lab = lab), aes(x = -Inf, y = Inf, label = lab),
                  inherit.aes = FALSE, hjust = -0.06, vjust = 1.5, size = 3.1,
                  colour = COL_FIT, fill = NA, label.color = NA,
                  label.padding = grid::unit(rep(0, 4), "pt")) +
    labs(x = "Published Slope (Webster et al.)", y = ylab,
         title = titled(letter, ttl), subtitle = wrap_md(sub, 62)) +
    theme_pub(11.5)
}
pB <- trait_panel("slope", "Pooled WGS Slope, published recipe", sw$slope, "B",
  "**Trait agreement, published recipe**",
  paste0("log<sub>2</sub>(f / f<sub>baseline</sub>) on day, frequencies floored ",
         "at 1/(4n). Axes are on the same scale."))
pC <- trait_panel("delta_slope", "Pooled WGS difference slope", sw$delta, "C",
  "**Trait agreement, difference slope**",
  paste0("The panel A trait against the published Slope. Ranks only: the axes ",
         "are on different scales."))

## D: every comparison on one axis
lad <- tibble::tribble(
  ~group,                                ~what,                                   ~key,
  "Same trait,\ntwo platforms",           "Difference slope (panel A)",            "platform_slope",
  "Same trait,\ntwo platforms",           "Per-sample frequency, median",          "per_sample_med",
  "Pooled WGS vs\npublished trait",       "Slope, published recipe (panel B)",     "wgs_recipe_slope",
  "Pooled WGS vs\npublished trait",       "Difference slope (panel C)",            "wgs_delta_slope",
  "Pooled WGS vs\npublished trait",       "PC1, published recipe",                 "wgs_recipe_pc1",
  "MIP-seq vs\npublished trait",          "Slope, published recipe (ceiling)",     "mip_recipe_slope",
  "MIP-seq vs\npublished trait",          "PC1, published recipe (ceiling)",       "mip_recipe_pc1",
  "MIP-seq vs\npublished trait",          "Difference slope",                      "mip_delta_slope",
  "Pooled WGS,\nhalf-min floor",          "Slope, published recipe",               "halfmin_slope",
  "Pooled WGS,\nhalf-min floor",          "PC1, published recipe",                 "halfmin_pc1") %>%
  left_join(V, by = "key") %>%
  mutate(group = factor(group, levels = unique(group)),
         what = factor(what, levels = rev(unique(what))))

pD <- ggplot(lad, aes(value, what)) +
  geom_segment(aes(x = 0.6, xend = value, yend = what), colour = "grey80",
               linewidth = 0.5) +
  geom_point(size = 2.6, colour = COL_PT) +
  geom_text(aes(label = sprintf("%.3f", value)), hjust = -0.35, size = 3.1,
            colour = COL_FIT) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(limits = c(0.6, 1.04), breaks = seq(0.6, 1, 0.1),
                     expand = expansion(0)) +
  labs(x = "Spearman's &rho;", y = NULL,
       title = titled("D", "**Every comparison on one axis**"),
       subtitle = wrap_md(paste0(
         "MIP-seq's own difference slope agrees with the published Slope no ",
         "better than ours, so the drop from A to B is the change of trait. ",
         "Ceiling: the published recipe on MIP-seq, restricted to the 15 ",
         "columns the pooled set covers."), 62)) +
  theme_pub(11.5) +
  theme(axis.title.x = element_markdown(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8.5),
        panel.spacing.y = grid::unit(4, "pt"))

fig <- (pA | pB) / (pC | pD)
write_fig(fig, "SUPP_FIG_XX_baugh_analyses", width = 12, height = 11)
