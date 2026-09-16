## Figure S1, with the strains that were never in the pool dropped ------------
##
##   Rscript scripts/DIAG_simulation_depth_present_only.R
##     -> plots/diagnostics/DIAG_simulation_depth_present_only.{pdf,png}
##
## A COMPANION to SUPP_FIG_XX_simulation_depth.R, not a replacement. That figure
## is computed over all 327 strains, which is how the simulation was always
## reported and what the 2021 record can be asserted against. This one asks the
## narrower question: among the strains that were actually IN a given trait's
## pool, how well are their frequencies recovered?
##
## WHY IT IS A DIFFERENT QUESTION. The original transform sends a strain with no
## published value for a trait to fitness exactly 0, so it is absent from that
## pool, and NNLS returns exactly 0 for an absent strain because zero is a hard
## boundary of the feasible set. Those points sit on the origin. They are not
## wrong and they are not padding -- correctly identifying an absent strain is a
## real thing the method does -- but they answer "did it notice the strain was
## missing", and the question behind a pooled experiment is "how precisely did
## it measure the strains that were there".
##
## The gap is smaller than the number of zeros suggests, and it does NOT scale
## with how many there are: see DIAG_pool_simulation_optimism.R, where the loss
## on dropping them is uncorrelated with the fraction absent (rho = -0.04,
## p = 0.96) and strongly correlated with how dispersed the present strains are
## (rho = -0.86, p = 0.02). A trait whose present strains sit in a tight cluster
## leans on the distant zeros for its dynamic range; a well-spread one does not.
##
## Exploratory, on the pool_optimization branch. Nothing in the manuscript reads
## this, and Figure S1 is unchanged.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(ggtext)
})

DEC  <- "supplemental_data/deconvolution"
DIAG <- "plots/diagnostics"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)
SFQ  <- file.path(DEC, "simulation_seeded_frequencies.tsv.gz")
stopifnot(file.exists(SFQ))

TRAIT_COL <- c(
  `PC1`                  = "#1B6C7A",
  `value`                = "#3E8E5A",
  `amsacrine_f.L1`       = "#7A9A3B",
  `assay_norm`           = "#C08A2E",
  `etoposide_median.TOF` = "#B5623C",
  `Albendazole_q75.TOF`  = "#9E4257",
  `mtDNA_ratio`          = "#6A5A8C")
panel_title <- function(l)
  paste0("<span style='font-size:13pt;color:#111111'>**", l, "**</span>")
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = base_size - 1),
          plot.title = element_markdown(size = base_size + 0.5),
          plot.subtitle = element_markdown(size = base_size - 2.5, colour = "grey30"),
          plot.title.position = "plot",
          legend.key.size = grid::unit(9, "pt"))
}
DEPTHS <- c(1, 3, 5, 10, 30, 50, 100, 500)

sfq <- read_tsv(SFQ, show_col_types = FALSE)
NREP <- n_distinct(sfq$replicate)

## r-squared both ways, so the companion carries its own comparison
r2 <- sfq %>% group_by(trait, depth, replicate) %>%
  summarise(all     = cor(frequency, input)^2,
            present = cor(frequency[input > 0], input[input > 0])^2,
            n_present = sum(input > 0), .groups = "drop") %>%
  group_by(trait, depth) %>%
  summarise(across(c(all, present), mean), n_present = first(n_present), .groups = "drop")

np <- r2 %>% distinct(trait, n_present) %>% arrange(desc(n_present))
cat("strains present per trait, of 327:\n")
print(as.data.frame(np), row.names = FALSE)
cat(sprintf("\nmean r-squared over the seven traits, %d replicates:\n", NREP))
print(as.data.frame(r2 %>% group_by(depth) %>%
  summarise(all = round(mean(all), 3), present = round(mean(present), 3),
            drop = round(mean(all - present), 3))), row.names = FALSE)

L <- r2 %>% pivot_longer(c(all, present), names_to = "set", values_to = "r2") %>%
  mutate(set = factor(set, c("all", "present"),
                      c("all 327 strains (Figure S1)", "pool strains only (this figure)")))

pA <- ggplot(L, aes(depth, r2, colour = trait)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.1) +
  scale_x_log10(breaks = DEPTHS, labels = paste0(DEPTHS, "x")) +
  scale_y_continuous(limits = c(NA, 1)) +
  scale_colour_manual(values = TRAIT_COL, name = NULL) +
  facet_wrap(~ set) +
  labs(x = "simulated sequencing depth", y = "r-squared against known input",
       title = paste0(panel_title("A"),
         "  Recovery against depth, with and without the absent strains"),
       subtitle = paste0("Mean over ", NREP, " replicates. Dropping the strains that were never in the pool costs ",
                         "0.19 of r-squared at 1x and 0.07 at 10x.")) +
  theme_pub() + theme(legend.position = "bottom")

## FALSE NEGATIVES. A strain that WAS in the pool but is returned as exactly
## zero. On log axes these cannot be drawn, so they are counted and annotated
## rather than silently dropped -- they are the failure mode that matters for a
## pooled experiment, since a zeroed strain contributes no phenotype at all.
fn <- sfq %>% filter(input > 0) %>% group_by(depth) %>%
  summarise(pct_zero = 100 * mean(frequency == 0), .groups = "drop")
cat("\npresent strains returned as exactly zero (false negatives):\n")
print(as.data.frame(fn %>% mutate(pct_zero = round(pct_zero, 2))), row.names = FALSE)

B <- sfq %>% filter(replicate == 1, input > 0) %>%
  mutate(depth_lab = factor(paste0(depth, "x"), paste0(DEPTHS, "x")))
## -Inf is NaN once the axis is log, so the label gets a real coordinate
FN <- fn %>% mutate(depth_lab = factor(paste0(depth, "x"), paste0(DEPTHS, "x")),
                    lab = sprintf("%.1f%% zeroed", pct_zero),
                    x = min(B$input), y = max(B$frequency))
pB <- ggplot(B, aes(input, frequency, colour = trait)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_point(size = 0.5, alpha = 0.45) +
  scale_x_log10() + scale_y_log10() +
  geom_text(data = FN, aes(x = x, y = y, label = lab), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 2.6, colour = "grey25") +
  scale_colour_manual(values = TRAIT_COL, guide = "none") +
  facet_wrap(~ depth_lab, nrow = 2) +
  labs(x = "known input frequency", y = "NNLS estimate",
       title = paste0(panel_title("B"), "  Estimate against truth, pool strains only"),
       subtitle = paste0("Replicate 1. Both axes log, which Figure S1 cannot do while the absent strains are in: ",
                         "their true value is exactly 0.<br>The annotation counts strains that WERE in the pool ",
                         "and were returned as exactly zero. They cannot be drawn on a log axis, and they are ",
                         "the<br>failure mode that matters -- a zeroed strain contributes no phenotype. At 10x, ",
                         "87 of 327 strains are zeroed in at least one trait; it tracks both<br>abundance ",
                         "(rho = -0.20) and privateness (rho = -0.25).")) +
  theme_pub()

fig <- pA / pB + plot_layout(heights = c(1, 1.25))
ggsave(file.path(DIAG, "DIAG_simulation_depth_present_only.pdf"), fig,
       width = 9, height = 9.5, device = cairo_pdf)
ggsave(file.path(DIAG, "DIAG_simulation_depth_present_only.png"), fig,
       width = 9, height = 9.5, dpi = 200, bg = "white")
cat("\nwrote DIAG_simulation_depth_present_only.{pdf,png}\n")
