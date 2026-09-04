library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggtext)
library(patchwork)
library(glue)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))

dir.create("../plots/plate_pheno_comparisons", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Compare manual plate phenotyping of wild isolates on pos-1 RNAi
# (data/pos1_plate_phenotyping/pos1_phenotypes_first2rounds.tsv; 0 = complete
# RNAi response/sensitive, 5 = no response/resistant) against:
#   1) Paaby et al. 2015 embryonic lethality for the pos-1 RNAi clone
#      (data/pooled_RNAi_expt/paaby2015/emb_leth_data.txt)
#   2) the pos1_original t2 (RNAi - ctrl) bootstrap frequency change
#      (data/pos1_original/pos1_strain_freq_change_mean.tsv, from
#      scripts/pos1_strain_freq_change_export.R)
#
# Plate score is a 6-level ordinal scale, not a continuous measurement, so
# Spearman rank correlation is used for both comparisons (same rationale as
# scripts/paaby_comparison.R) -- it doesn't require the two variables to
# share units or a linear relationship, only a monotonic one.
#
# Expected direction of each relationship (for interpreting the plots):
#   - plate score vs Paaby lethality: NEGATIVE. Low score = sensitive =
#     strong RNAi response; high lethality = sensitive. So low score should
#     pair with high lethality.
#   - plate score vs pos1_original delta freq: POSITIVE. Sensitive strains
#     lose pool frequency under pos-1 RNAi relative to control (delta_freq
#     = pos1_frq - ctrl_frq, more negative = more sensitive), so low plate
#     score should pair with more negative delta_freq.
# ---------------------------------------------------------------------------

plate <- readr::read_tsv("../data/pos1_plate_phenotyping/pos1_phenotypes_first2rounds.tsv", show_col_types = FALSE) %>%
  dplyr::rename(plate_score = trait)

# ---------------------------------------------------------------------------
# Paaby et al. 2015: pos-1 embryonic lethality per strain (mean across wells)
# ---------------------------------------------------------------------------

paaby_raw <- data.table::fread("../data/pooled_RNAi_expt/paaby2015/emb_leth_data.txt")

paaby_pos1 <- paaby_raw %>%
  dplyr::filter(vector == "pos-1") %>%
  dplyr::filter(!is.na(eggs), !is.na(larvae), (eggs + larvae) > 0) %>%
  dplyr::mutate(leth = eggs / (eggs + larvae)) %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(mean_leth = mean(leth), n_wells = dplyr::n(), .groups = "drop")

# ---------------------------------------------------------------------------
# pos1_original: t2 (pos-1 RNAi - ctrl) bootstrap frequency change per strain
# ---------------------------------------------------------------------------

pos1orig <- readr::read_tsv("../data/pos1_original/pos1_strain_freq_change_mean.tsv", show_col_types = FALSE)

# ---------------------------------------------------------------------------
# comparison 1: plate score vs Paaby lethality
# ---------------------------------------------------------------------------

cmp_paaby <- dplyr::inner_join(plate, paaby_pos1, by = "strain")

ct_paaby <- suppressWarnings(cor.test(cmp_paaby$plate_score, cmp_paaby$mean_leth, method = "spearman"))

p_paaby <- ggplot(cmp_paaby, aes(x = factor(plate_score), y = mean_leth)) +
  geom_boxplot(outlier.shape = NA, fill = "grey95") +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.7) +
  ggrepel::geom_text_repel(aes(label = strain), size = 2.8, min.segment.length = 0, seed = 42) +
  labs(
    title = glue::glue("Plate pos-1 RNAi score vs Paaby et al. 2015 embryonic lethality (n = {nrow(cmp_paaby)})"),
    subtitle = glue::glue("Spearman rho = {round(unname(ct_paaby$estimate), 3)}, p = {signif(ct_paaby$p.value, 3)}"),
    x = "Plate score (0 = complete response, 5 = no response)",
    y = "Paaby et al. 2015 mean embryonic lethality [eggs / (eggs+larvae)]"
  ) +
  theme_bw(13) +
  theme(plot.title = element_text(size = 12))

ggsave(p_paaby, filename = "../plots/plate_pheno_comparisons/plate_vs_paaby.pdf", height = 6, width = 8)
ggsave(p_paaby, filename = "../plots/plate_pheno_comparisons/plate_vs_paaby.png", height = 6, width = 8, dpi = 300)

# ---------------------------------------------------------------------------
# comparison 2: plate score vs pos1_original t2 (pos-1 - ctrl) delta freq
# ---------------------------------------------------------------------------

cmp_pos1orig <- dplyr::inner_join(plate, pos1orig, by = "strain")

ct_pos1orig <- suppressWarnings(cor.test(cmp_pos1orig$plate_score, cmp_pos1orig$mean_delta_freq, method = "spearman"))

p_pos1orig <- ggplot(cmp_pos1orig, aes(x = factor(plate_score), y = mean_delta_freq)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_boxplot(outlier.shape = NA, fill = "grey95") +
  geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.5) +
  labs(
    title = glue::glue("Plate pos-1 RNAi score vs pos1_original t2 (pos-1 - ctrl) frequency change (n = {nrow(cmp_pos1orig)})"),
    subtitle = glue::glue("Spearman rho = {round(unname(ct_pos1orig$estimate), 3)}, p = {signif(ct_pos1orig$p.value, 3)}"),
    x = "Plate score (0 = complete response, 5 = no response)",
    y = "pos1_original: mean (pos-1 RNAi - ctrl) t2 frequency"
  ) +
  theme_bw(13) +
  theme(plot.title = element_text(size = 12))

ggsave(p_pos1orig, filename = "../plots/plate_pheno_comparisons/plate_vs_pos1original.pdf", height = 6, width = 8)
ggsave(p_pos1orig, filename = "../plots/plate_pheno_comparisons/plate_vs_pos1original.png", height = 6, width = 8, dpi = 300)

# ---------------------------------------------------------------------------
# summary table + merged data (for the visual report)
# ---------------------------------------------------------------------------

comparison_summary <- tibble::tibble(
  comparison = c("plate_score_vs_paaby_lethality", "plate_score_vs_pos1original_delta_freq"),
  n_strains = c(nrow(cmp_paaby), nrow(cmp_pos1orig)),
  spearman_rho = c(unname(ct_paaby$estimate), unname(ct_pos1orig$estimate)),
  p_value = c(ct_paaby$p.value, ct_pos1orig$p.value)
)

print(comparison_summary)

dir.create("../data/pos1_plate_phenotyping/comparisons", recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(comparison_summary, "../data/pos1_plate_phenotyping/comparisons/comparison_summary.tsv")
readr::write_tsv(cmp_paaby, "../data/pos1_plate_phenotyping/comparisons/plate_vs_paaby_merged.tsv")
readr::write_tsv(cmp_pos1orig, "../data/pos1_plate_phenotyping/comparisons/plate_vs_pos1original_merged.tsv")

# ---------------------------------------------------------------------------
# 2-panel supplemental figure: panel A = pooled bulk-competition phenotyping
# (pos1_original) on top, panel B = Paaby et al. 2015 embryonic lethality on
# bottom, sharing an x-axis (plate score) via patchwork stacking. rho and n
# are embedded as an in-plot annotation (corner text) rather than a
# title/subtitle, per the figure's intended use as a standalone supplemental
# panel.
# ---------------------------------------------------------------------------

rho_label <- function(rho, n) glue::glue("Spearman's ρ = {sprintf('%.2f', rho)}\nn = {n}")

panel_A <- ggplot(cmp_pos1orig, aes(x = factor(plate_score), y = mean_delta_freq)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_boxplot(outlier.shape = NA, fill = "grey95", width = 0.6) +
  geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.5) +
  annotate("text", x = -Inf, y = Inf, label = rho_label(unname(ct_pos1orig$estimate), nrow(cmp_pos1orig)),
           hjust = -0.08, vjust = 1.3, size = 3.6, lineheight = 0.95) +
  scale_x_discrete(limits = as.character(0:5)) +
  labs(x = NULL, y = "Pooled *pos-1* phenotypes<br>(*pos-1* RNAi − HT115 control)") +
  theme_bw(13) +
  theme(
    axis.title.y = ggtext::element_markdown(lineheight = 1.1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 10, r = 12, b = 4, l = 10)
  )

panel_B <- ggplot(cmp_paaby, aes(x = factor(plate_score), y = mean_leth)) +
  geom_boxplot(outlier.shape = NA, fill = "grey95", width = 0.6) +
  geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.5) +
  annotate("text", x = Inf, y = Inf, label = rho_label(unname(ct_paaby$estimate), nrow(cmp_paaby)),
           hjust = 1.08, vjust = 1.3, size = 3.6, lineheight = 0.95) +
  scale_x_discrete(limits = as.character(0:5)) +
  labs(
    x = "Plate score (0 = complete response, 5 = no response)",
    y = "Paaby et al. 2015 mean\nembryonic lethality"
  ) +
  theme_bw(13) +
  theme(plot.margin = margin(t = 4, r = 12, b = 10, l = 10))

supp_figure <- panel_A / panel_B +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(supp_figure, filename = "../plots/SUPP_FIG_plate_vs_paaby_vs_pos1original.pdf", height = 9, width = 6.5, device = cairo_pdf)
ggsave(supp_figure, filename = "../plots/SUPP_FIG_plate_vs_paaby_vs_pos1original.png", height = 9, width = 6.5, dpi = 300, type = "cairo")
