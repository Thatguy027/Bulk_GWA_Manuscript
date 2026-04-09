library(tidyverse)

# Read data
df <- read_csv("../data/plate_rnai_phenotyping/20260409_ju2466swap_plus_N2A_swap.csv")

# Compute 95% Wilson score confidence intervals for binomial proportions
# n_hatched = n_plated - n_unhatched, n = n_plated
z <- qnorm(0.975)

df <- df %>%
  mutate(
    n_hatched = n_plated - n_unhatched,
    p = fraction_hatched,
    # Wilson score interval
    ci_lo = (p + z^2 / (2 * n_plated) - z * sqrt(p * (1 - p) / n_plated + z^2 / (4 * n_plated^2))) /
      (1 + z^2 / n_plated),
    ci_hi = (p + z^2 / (2 * n_plated) + z * sqrt(p * (1 - p) / n_plated + z^2 / (4 * n_plated^2))) /
      (1 + z^2 / n_plated),
    ci_lo = pmax(ci_lo, 0),
    ci_hi = pmin(ci_hi, 1)
  )

# Set factor levels to preserve row order for x-axis
motif_order <- unique(df$`glycosylation motif`)
df <- df %>%
  mutate(
    `glycosylation motif` = factor(`glycosylation motif`, levels = motif_order),
    condition = factor(condition, levels = c("pos", "ht115"))
  )

# Plot
p <- ggplot(df, aes(x = `glycosylation motif`, y = fraction_hatched,
                    fill = condition, group = condition)) +
  geom_col(
    aes(color = condition),
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = ci_lo, ymax = ci_hi),
    position = position_dodge(width = 0.8),
    width = 0.25,
    color = "black",
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = c("pos" = "firebrick3", "ht115" = "black"),
    labels = c("pos" = expression(italic("pos-1")), "ht115" = "HT115"),
    name = "Condition"
  ) +
  scale_y_continuous(
    limits = c(0, 1.05),
    expand = c(0, 0),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Glycosylation motif",
    y = "Fraction hatched"
  ) +
  theme_bw(18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

print(p)

ggsave("../plots/Figure4.pdf", p, width = 10, height = 6)
ggsave("../plots/Figure4.png", p, width = 10, height = 6, dpi = 300)

# --- Figure 4 clean: pos-1 RNAi only, pairwise significance tests ---

df_pos <- df %>%
  filter(condition == "pos", strain != "JU2466_B") %>%
  mutate(`glycosylation motif` = str_replace(`glycosylation motif`, "_A", ""))

# Group JU1793 backgrounds together, then JU2466 backgrounds
motif_levels_clean <- c("JU1793[NxT]", "JU1793[NxK]", "JU1793[AxT]",
                        "JU2466[NxK]", "JU2466[NxT]", "JU2466[AxK]")

df_pos <- df_pos %>%
  mutate(`glycosylation motif` = factor(`glycosylation motif`, levels = motif_levels_clean))

# Fisher's exact test for two strains (count data: hatched vs unhatched)
run_fisher <- function(s1, s2) {
  d1 <- df_pos %>% filter(strain == s1)
  d2 <- df_pos %>% filter(strain == s2)
  mat <- matrix(c(d1$n_hatched, d1$n_unhatched,
                  d2$n_hatched, d2$n_unhatched),
                nrow = 2, byrow = TRUE)
  fisher.test(mat)$p.value
}

fmt_p <- function(p) {
  if (p < 0.001) "p < 0.001"
  else if (p < 0.01) sprintf("p = %.3f", p)
  else sprintf("p = %.3f", p)
}

# Comparisons: strain names map to x-axis positions in motif_levels_clean
comparisons <- tibble(
  s1      = c("JU1793",   "JU1793",   "JU2466_A", "JU2466_A"),
  s2      = c("wSZ200",   "wSZ207",   "wSZ206",   "wSZ208"),
  x1      = c(1,          1,          4,          4),
  x2      = c(2,          3,          5,          6),
  y_level = c(1,          2,          1,          2)
) %>%
  rowwise() %>%
  mutate(pval  = run_fisher(s1, s2),
         label = fmt_p(pval)) %>%
  ungroup()

# Bracket geometry parameters
y_base  <- 1.05
y_step  <- 0.09
tip_len <- 0.025

bracket_segs <- comparisons %>%
  mutate(yb = y_base + (y_level - 1) * y_step,
         yt = yb + tip_len) %>%
  rowwise() %>%
  reframe(
    x    = c(x1, x1, x2),
    xend = c(x2, x1, x2),
    y    = c(yb, yb, yb),
    yend = c(yb, yt, yt)
  )

bracket_text <- comparisons %>%
  mutate(
    x     = (x1 + x2) / 2,
    y     = y_base + (y_level - 1) * y_step + tip_len + 0.018
  )

p_clean <- ggplot(df_pos, aes(x = `glycosylation motif`, y = fraction_hatched)) +
  geom_col(fill = "firebrick3", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.25, color = "black", linewidth = 0.5) +
  geom_segment(data = bracket_segs,
               aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = bracket_text,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 5) +
  scale_y_continuous(
    expand = c(0, 0),
    labels = scales::percent_format(accuracy = 1)
  ) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  labs(x = "Strain genotype", y = "Fraction hatched") +
  theme_bw(18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 90, r = 10, b = 10, l = 10, unit = "pt")
  )

print(p_clean)

ggsave("../plots/Figure4_clean.pdf", p_clean, width = 10, height = 7)
ggsave("../plots/Figure4_clean.png", p_clean, width = 10, height = 7, dpi = 300)


