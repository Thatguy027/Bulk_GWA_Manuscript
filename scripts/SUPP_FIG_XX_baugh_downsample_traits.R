## Sequencing depth against the traits that actually get mapped ---------------
##
## The existing depth supplement asks how well downsampled DECONVOLUTION
## FREQUENCIES track full-depth ones, and how well the difference-based slope
## tracks MIP-seq. This asks the question the mapping actually depends on: at
## each depth, how well do the traits recover the PUBLISHED eLife phenotypes?
##
## Traits are rebuilt from scratch at every depth on the published recipe --
## log2(f_day / f_baseline), prcomp(scale, center) for PC1 and the regression of
## that ratio on day excluding day 17 for Slope -- with the low-frequency floor
## at 1/(4n), a quarter of an equal share. The floor matters more at low depth,
## because subsampling produces more exact zeros: 1,897 of 14,076 cells are zero
## across the six depths, against 339 at full depth.
##
## Three lines are drawn, and all three are pooled-WGS traits scored against the
## PUBLISHED eLife phenotypes -- they differ in how the trait is built, not in
## what they are compared to:
##
##   PC1           published recipe, scored against the published PC1
##   Slope         published recipe -- the slope of log2(f/baseline) on day --
##                 scored against the published Slope
##   Slope, delta  the difference-based slope, f minus its day-1 value regressed
##                 on day, on RAW frequencies with no floor, scored against the
##                 same published Slope
##
## The last exists because it takes no logarithm and so cannot be broken by a
## zero, which is what the depth axis stresses. It is robust where the log-ratio
## traits are not.
##
## Uses the deposited downsampling output, supplemental_data/deconvolution/
## baugh_downsampled_slopes.rda (102 strains x 23 samples x 6 depths), produced
## by scripts/baugh_L1_DownSample_Counts.R.
##
## Writes plots/SUPP_FIG_XX_baugh_downsample_traits.{pdf,png} and
## supplemental_data/deconvolution/baugh_downsample_trait_recovery.tsv
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(ggrepel)
})

DEC <- "supplemental_data/deconvolution"
PH  <- "supplemental_data/phenotypes"
OUT <- "plots"

e <- new.env(); load(file.path(DEC, "baugh_downsampled_slopes.rda"), e)
ds <- as_tibble(e$ds_predictions_df)
pub <- read_tsv(file.path(PH, "baugh_published_traits.txt"), show_col_types = FALSE)

ols <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); if (sum(ok) < 2) return(NA_real_)
  x <- x[ok]; y <- y[ok]; xc <- x - mean(x); s <- sum(xc^2)
  if (s == 0) NA_real_ else sum(xc * y) / s
}

## traits from one frequency column, published recipe with the 1/(4n) floor
traits <- function(d, col) {
  n_str <- n_distinct(d$strain)
  fl <- 1 / (4 * n_str)
  l <- d %>% rename(f_raw = all_of(col)) %>%
    mutate(f = pmax(f_raw, fl)) %>%
    separate(sample, into = c("rep", "day"), sep = "_", remove = FALSE, extra = "merge") %>%
    mutate(is_bl = grepl("baseline", sample),
           dnum = as.numeric(gsub("d", "", sub("_baseline", "", day))),
           day2 = ifelse(is_bl, "BL", paste0("d", dnum)))
  bl <- l %>% filter(day2 == "BL") %>% select(strain, rep, base = f)
  w  <- l %>% filter(day2 != "BL") %>% left_join(bl, by = c("strain", "rep")) %>%
    mutate(l2 = log2(f / base), colk = paste0(rep, "_", day2)) %>%
    filter(is.finite(l2))
  m <- w %>% select(strain, colk, l2) %>%
    pivot_wider(names_from = colk, values_from = l2) %>%
    column_to_rownames("strain") %>% as.matrix()
  m <- m[, colSums(is.na(m)) < 0.1 * nrow(m), drop = FALSE]
  m <- m[complete.cases(m), , drop = FALSE]
  p <- prcomp(m, scale. = TRUE, center = TRUE)
  ## The difference-based slope, for contrast. Computed on the RAW frequencies:
  ## it takes no logarithm, so the floor is neither needed nor appropriate, and
  ## applying it would clip real low-frequency variation. (Flooring it anyway
  ## moves the result by 0.001 to 0.017, so this is a correctness point rather
  ## than a material one.) Fitted pooled across replicate arms, which for this
  ## balanced design with a common set of days is identical to fitting each arm
  ## and averaging, as platform_slopes() does.
  d1 <- l %>% filter(!is_bl, dnum == 1) %>% select(strain, rep, first = f_raw)
  dl <- l %>% filter(!is_bl, dnum != 17) %>% left_join(d1, by = c("strain", "rep")) %>%
    mutate(dv = f_raw - first) %>% group_by(strain) %>%
    summarise(delta_slope = ols(dnum, dv), .groups = "drop")
  tibble(strain = rownames(m), pc1 = p$x[, 1]) %>%
    left_join(w %>% filter(day2 != "d17") %>% group_by(strain) %>%
                summarise(slope = ols(dnum, l2), .groups = "drop"), by = "strain") %>%
    left_join(dl, by = "strain")
}

sgn <- function(a, b) if (isTRUE(cor(a, b, method = "spearman") < 0)) -a else a

rows <- map_dfr(sort(unique(ds$ds_n)), function(dep) {
  t <- traits(ds %>% filter(ds_n == dep), "ds_frq")
  j <- pub %>% inner_join(t, by = "strain") %>% filter(strain != "N2")
  tibble(depth = dep, n = nrow(j),
         PC1         = abs(cor(j$PC1,   j$pc1,         method = "spearman")),
         Slope       = abs(cor(j$Slope, j$slope,       method = "spearman")),
         `Slope, delta` = abs(cor(j$Slope, j$delta_slope, method = "spearman")))
})
tf <- traits(ds %>% filter(ds_n == max(ds_n)) %>% distinct(strain, sample, .keep_all = TRUE) %>%
               mutate(ds_frq = frq), "ds_frq")
jf <- pub %>% inner_join(tf, by = "strain") %>% filter(strain != "N2")
full <- tibble(depth = Inf, n = nrow(jf),
  PC1 = abs(cor(jf$PC1, jf$pc1, method = "spearman")),
  Slope = abs(cor(jf$Slope, jf$slope, method = "spearman")),
  `Slope, delta` = abs(cor(jf$Slope, jf$delta_slope, method = "spearman")))

res <- bind_rows(rows, full)
write_tsv(res, file.path(DEC, "baugh_downsample_trait_recovery.tsv"))
cat("recovery of the published traits by sequencing depth\n")
print(as.data.frame(res %>% mutate(across(where(is.numeric), ~round(.x, 3)))), row.names = FALSE)

long <- res %>% filter(is.finite(depth)) %>%
  pivot_longer(c(PC1, Slope, `Slope, delta`), names_to = "trait", values_to = "rho")
fl <- res %>% filter(!is.finite(depth)) %>%
  pivot_longer(c(PC1, Slope, `Slope, delta`), names_to = "trait", values_to = "rho")

COL <- c(PC1 = "#C4302B", Slope = "#2E4057", `Slope, delta` = "#1A7F5A")
theme_set(theme_bw(9) + theme(
  plot.title = element_text(face = "bold", size = 9.5),
  plot.subtitle = element_text(size = 7.6, colour = "grey30"),
  panel.grid.minor = element_blank(), legend.position = "top",
  legend.title = element_blank()))

pA <- ggplot(long, aes(depth, rho, colour = trait)) +
  geom_hline(data = fl, aes(yintercept = rho, colour = trait),
             linetype = 2, linewidth = 0.4, show.legend = FALSE) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  scale_colour_manual(values = COL) +
  scale_x_log10(breaks = sort(unique(long$depth)),
                labels = function(x) paste0(x, "x")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "A  Recovery of the published eLife traits against depth",
       subtitle = "Spearman against the published Slope and PC1. Dashed lines are full depth.",
       x = "Sequencing depth", y = "Spearman rho vs published trait")

zer <- ds %>% group_by(ds_n) %>%
  summarise(pz = mean(ds_frq == 0), .groups = "drop")
pB <- ggplot(zer, aes(ds_n, pz)) +
  geom_line(colour = "grey40", linewidth = 0.7) +
  geom_point(size = 2, colour = "grey25") +
  scale_x_log10(breaks = zer$ds_n, labels = function(x) paste0(x, "x")) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "B  Why the log-ratio traits suffer at low depth",
       subtitle = "Share of strain x sample cells the deconvolution sets to exactly zero",
       x = "Sequencing depth", y = "cells at exactly zero")

fig <- pA / pB + plot_layout(heights = c(1.5, 1))
ggsave(file.path(OUT, "SUPP_FIG_XX_baugh_downsample_traits.pdf"), fig, width = 7.2, height = 6.6)
ggsave(file.path(OUT, "SUPP_FIG_XX_baugh_downsample_traits.png"), fig, width = 7.2, height = 6.6, dpi = 200)
cat(sprintf("\nwrote %s/SUPP_FIG_XX_baugh_downsample_traits.{pdf,png}\n", OUT))
