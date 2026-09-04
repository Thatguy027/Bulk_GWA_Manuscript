## Supplement -- pooled phenotype across all targets, vst scale --------------
## Every strain against every RNAi target on the vst scale, so a strain that is
## an outlier for one target can be told from one that is an outlier generally.
## XZ1516 sits at the sensitive extreme for nearly every target, which is what
## makes it a general rather than target-specific outlier.
##
##   Rscript scripts/SUPP_FIG_XX_pooled_phenotype_heatmap.R
##
## Phenotype
## ---------
## Uses the vst_ctrl_<gene>_T2 columns of
## data/pooled_RNAi_expt/reanalysis/vst_association_traits.csv -- the exact
## traits the GEMMA association was run on, so this figure is on the same scale
## as the Manhattan in Figure 2.
##
## 9 of the 93 panel strains have no vst value and are dropped, so n = 84.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggtext)
})

TRAITS <- "data/pooled_RNAi_expt/reanalysis/vst_association_traits.csv"
OUT    <- "plots/pooled_cross_intersection"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PARENTS <- c("XZ1516", "JU1793", "JU2466")
COL_HI  <- "#C4302B"

v <- read_csv(TRAITS, show_col_types = FALSE)

hm <- v %>%
  select(strain, starts_with("vst_ctrl_")) %>%
  pivot_longer(-strain, names_to = "gene", values_to = "vst") %>%
  mutate(gene = str_replace(gene, "^vst_ctrl_", "") %>% str_replace("_T2$", "")) %>%
  group_by(strain) %>% filter(any(!is.na(vst))) %>% ungroup()

n_tot  <- n_distinct(v$strain)
n_used <- n_distinct(hm$strain)
n_gene <- n_distinct(hm$gene)
cat("panel strains: ", n_tot, " | with a vst value: ", n_used,
    " | targets: ", n_gene, "\n\n", sep = "")

## strains most sensitive at the left, targets by how much they vary
ord  <- hm %>% group_by(strain) %>% summarise(m = mean(vst, na.rm = TRUE)) %>%
  arrange(m) %>% pull(strain)
gord <- hm %>% group_by(gene) %>% summarise(s = sd(vst, na.rm = TRUE)) %>%
  arrange(desc(s)) %>% pull(gene)

lim <- max(abs(hm$vst), na.rm = TRUE)

fig <- hm %>%
  mutate(strain = factor(strain, levels = ord),
         gene   = factor(gene, levels = gord)) %>%
  ggplot(aes(strain, gene, fill = vst)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#8C2D04", mid = "grey96", high = "#08519C",
                       midpoint = 0, limits = c(-lim, lim),
                       na.value = "grey88",
                       name = "pooled\nresponse\n(vst)") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6,
                                   face = ifelse(ord %in% PARENTS, "bold", "plain"),
                                   colour = ifelse(ord %in% PARENTS, COL_HI, "grey35")),
        axis.text.y = element_text(size = 9, face = "italic"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5),
        legend.key.height = grid::unit(24, "pt"),
        plot.margin = margin(6, 8, 4, 6))

ggsave(file.path(OUT, "SUPP_FIG_XX_pooled_phenotype_heatmap.pdf"), fig,
       width = 11, height = 3.6, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_pooled_phenotype_heatmap.png"), fig,
       width = 11, height = 3.6, dpi = 300, bg = "white")
cat("wrote SUPP_FIG_XX_pooled_phenotype_heatmap.{pdf,png}\n")

## ---------------------------------------------------------------------------
cat("\n== most sensitive and most resistant strains, by mean vst ==\n")
mm <- hm %>% group_by(strain) %>% summarise(mean_vst = mean(vst, na.rm = TRUE)) %>%
  arrange(mean_vst)
print(as.data.frame(bind_rows(
  mm %>% slice_head(n = 5) %>% mutate(end = "most sensitive"),
  mm %>% slice_tail(n = 5) %>% mutate(end = "most resistant")) %>%
  transmute(end, strain, mean_vst = round(mean_vst, 3))), row.names = FALSE)

cat("\n== cross parents ==\n")
print(as.data.frame(mm %>% filter(strain %in% PARENTS) %>%
  mutate(rank = rank(mean_vst)[match(strain, strain)]) %>%
  transmute(strain, mean_vst = round(mean_vst, 3),
            rank_of = paste0(match(strain, mm$strain), " of ", nrow(mm)))),
  row.names = FALSE)
