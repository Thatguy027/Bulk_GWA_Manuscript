## Supplement -- the plate assay against two independent measurements --------
##
##   Rscript scripts/SUPP_FIG_plate_vs_paaby_vs_pos1original.R
##     -> plots/SUPP_FIG_plate_vs_paaby_vs_pos1original.{pdf,png}
##
##   A  manual plate score against the pooled 2023 pos-1 phenotype, VST scale
##   B  manual plate score against Paaby et al. 2015 embryonic lethality
##
## Replaces the panel A that scripts/plate_pheno_comparisons.R produced. Two
## changes:
##
## 1  PANEL A IS NOW THE VST TRAIT, vst_ctrl_pos-1_T2 from
##    supplemental_data/phenotypes/pos1_2023_association_traits.csv.gz -- the same
##    trait the association mapping was run on, so this figure is on the same
##    scale as every other phenotype panel in the manuscript. The old version
##    used a bootstrap frequency change from a different export, which was on
##    no scale used elsewhere.
##
## 2  IT RUNS UNDER Rscript. The old script resolved its working directory
##    through the RStudio API and ran only inside RStudio.
##
## DIRECTIONS, so the signs can be checked rather than trusted
## The plate score is ordinal, 0 = complete RNAi response (sensitive) to
## 5 = no response (resistant).
##   - VST is positive for strains that GAINED pool frequency under pos-1 RNAi,
##     i.e. resistant, so the correlation with plate score should be POSITIVE.
##   - Paaby lethality is high for sensitive strains, so the correlation with
##     plate score should be NEGATIVE.
## Both hold; the console prints the values.
##
## Spearman throughout: the plate score is a six-level ordinal scale, so a rank
## correlation is the only defensible choice -- it needs a monotonic
## relationship, not shared units or linearity.
##
## CAVEAT: the Paaby comparison rests on the strains shared with that study and
## is far smaller than panel A. It should be described as consistent, not as
## independent confirmation.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(ggtext)
})

OUT   <- "plots"
PLATE <- "supplemental_data/phenotypes/plate_scores_pos1.tsv"
TRAIT <- "supplemental_data/phenotypes/pos1_2023_association_traits.csv.gz"
PAABY <- "supplemental_data/phenotypes/paaby2015_embryonic_lethality.txt.gz"

COL_PT <- "#2E4057"

panel_title <- function(letter) {
  paste0("<span style='font-size:13pt;color:#111111'>**", letter, "**</span>")
}
theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_markdown(size = base_size + 0.5),
          plot.title.position = "plot")
}
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

plate <- read_tsv(PLATE, show_col_types = FALSE) %>% rename(plate_score = trait)
msg("plate scores: ", nrow(plate), " strains")

## --- A: the VST pooled phenotype -----------------------------------------
vst <- read_csv(TRAIT, show_col_types = FALSE) %>%
  transmute(strain, vst = `vst_ctrl_pos-1_T2`) %>% filter(is.finite(vst))
a <- inner_join(plate, vst, by = "strain")
ct_a <- suppressWarnings(cor.test(a$plate_score, a$vst, method = "spearman"))
msg("panel A: ", nrow(a), " strains | rho = ", sprintf("%+.3f", ct_a$estimate),
    " p = ", signif(ct_a$p.value, 3))

## --- B: Paaby et al. 2015 embryonic lethality ----------------------------
## lethality per well as unhatched eggs over eggs plus larvae, averaged over
## wells within a strain, for the pos-1 clone only
paaby <- fread(PAABY) %>% as_tibble() %>%
  filter(vector == "pos-1", !is.na(eggs), !is.na(larvae), (eggs + larvae) > 0) %>%
  mutate(leth = eggs / (eggs + larvae)) %>%
  group_by(strain) %>%
  summarise(mean_leth = mean(leth), n_wells = n(), .groups = "drop")
b <- inner_join(plate, paaby, by = "strain")
ct_b <- suppressWarnings(cor.test(b$plate_score, b$mean_leth, method = "spearman"))
msg("panel B: ", nrow(b), " strains | rho = ", sprintf("%+.3f", ct_b$estimate),
    " p = ", signif(ct_b$p.value, 3))

cat("\n== plate score vs each measurement (Spearman) ==\n")
print(data.frame(
  comparison = c("pooled pos-1 response (VST)", "Paaby 2015 embryonic lethality"),
  n_strains = c(nrow(a), nrow(b)),
  rho = round(c(ct_a$estimate, ct_b$estimate), 3),
  p_value = signif(c(ct_a$p.value, ct_b$p.value), 3),
  expected_sign = c("positive", "negative")), row.names = FALSE)
cat("\n== strains per plate score, panel A ==\n")
print(as.data.frame(a %>% count(plate_score)), row.names = FALSE)

## the comparison sign lives in the string, so "p < 1e-4" does not come out
## as "p = < 1e-4"
lab <- function(ct, n) sprintf("Spearman's &rho; = %.2f<br>%s<br>n = %d",
                               ct$estimate,
                               if (ct$p.value < 1e-4) "*p* &lt; 1e-4"
                               else sprintf("*p* = %s", signif(ct$p.value, 2)),
                               n)

pA <- ggplot(a, aes(factor(plate_score), vst)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4,
             colour = "grey60") +
  geom_boxplot(outlier.shape = NA, fill = "grey95", width = 0.6,
               linewidth = 0.35) +
  ## seeded, so the figure is reproducible from the deposit
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1),
             size = 1.5, alpha = 0.5, colour = COL_PT) +
  geom_richtext(data = tibble(x = -Inf, y = Inf, l = lab(ct_a, nrow(a))),
                aes(x, y, label = l), inherit.aes = FALSE, hjust = -0.08,
                vjust = 1.15, size = 3.2, fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  scale_x_discrete(limits = as.character(0:5)) +
  labs(x = NULL, y = "Pooled *pos-1* response (VST)",
       title = panel_title("A")) +
  theme_pub() +
  theme(axis.title.y = element_markdown(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(t = 8, r = 12, b = 4, l = 10))

pB <- ggplot(b, aes(factor(plate_score), mean_leth)) +
  geom_boxplot(outlier.shape = NA, fill = "grey95", width = 0.6,
               linewidth = 0.35) +
  ## seeded, so the figure is reproducible from the deposit
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1),
             size = 1.5, alpha = 0.5, colour = COL_PT) +
  geom_richtext(data = tibble(x = Inf, y = Inf, l = lab(ct_b, nrow(b))),
                aes(x, y, label = l), inherit.aes = FALSE, hjust = 1.08,
                vjust = 1.15, size = 3.2, fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  scale_x_discrete(limits = as.character(0:5)) +
  labs(x = "Plate score (0 = complete response, 5 = no response)",
       y = "Paaby et al. 2015 mean<br>embryonic lethality",
       title = panel_title("B")) +
  theme_pub() +
  theme(axis.title.y = element_markdown(lineheight = 1.1),
        plot.margin = margin(t = 4, r = 12, b = 8, l = 10))

fig <- pA / pB
ggsave(file.path(OUT, "SUPP_FIG_plate_vs_paaby_vs_pos1original.pdf"), fig,
       width = 6.5, height = 8.4, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_plate_vs_paaby_vs_pos1original.png"), fig,
       width = 6.5, height = 8.4, dpi = 300, bg = "white")
msg("wrote SUPP_FIG_plate_vs_paaby_vs_pos1original.{pdf,png}")
