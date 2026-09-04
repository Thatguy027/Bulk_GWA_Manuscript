## Supplement -- pooled phenotype, ranked, on the vst scale -------------------
## Per-strain pooled RNAi response for mig-6 and pos-1, ranked, with the cross
## parents highlighted. This is the panel that shows the crosses were built from
## the extremes of the pooled assay.
##
##   Rscript scripts/SUPP_FIG_XX_pooled_phenotype_ranks.R
##
## Phenotype
## ---------
## Uses the vst_ctrl_<gene>_T2 columns of
## supplemental_data/phenotypes/pooled_vst_traits.csv.gz -- the exact
## traits the GEMMA association was run on, so this figure and the Manhattan in
## Figure 2 are on the same scale. The file also carries delta_ctrl_ (the raw
## RNAi-minus-control frequency change) and log2fc_ columns; the vst column is
## the one used here.
##
## 9 of the 93 panel strains have no vst value, so n = 84. Ranks below are out
## of those 84, not 93.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggtext)
  library(ggrepel)
  library(patchwork)
})

TRAITS <- "supplemental_data/phenotypes/pooled_vst_traits.csv.gz"
## a kept supplement (see README.md), so it writes to plots/ -- it used to
## write into plots/pooled_cross_intersection/, which is not tracked, so
## re-running it left the curated copy stale
OUT    <- "plots"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PARENTS <- c("XZ1516", "JU1793", "JU2466")
COL_HI  <- "#C4302B"
GENES   <- c("mig-6", "pos-1")

v <- read_csv(TRAITS, show_col_types = FALSE)

ph <- v %>%
  select(strain, all_of(paste0("vst_ctrl_", GENES, "_T2"))) %>%
  pivot_longer(-strain, names_to = "gene", values_to = "vst") %>%
  mutate(gene = str_replace(gene, "^vst_ctrl_", "") %>% str_replace("_T2$", "")) %>%
  filter(!is.na(vst))

n_tot  <- n_distinct(v$strain)
n_used <- n_distinct(ph$strain)
cat("panel strains: ", n_tot, " | with a vst value: ", n_used, "\n\n", sep = "")

rank_panel <- function(g) {
  d <- ph %>% filter(gene == g) %>%
    arrange(vst) %>%
    mutate(i = row_number(), is_par = strain %in% PARENTS)
  ggplot(d, aes(i, vst)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_segment(aes(xend = i, yend = 0), linewidth = 0.28, colour = "grey72") +
    geom_point(aes(colour = is_par, size = is_par)) +
    geom_text_repel(data = d %>% filter(is_par), aes(label = strain),
                    size = 3.2, fontface = "bold", seed = 42,
                    min.segment.length = 0, box.padding = 0.5,
                    segment.size = 0.3, colour = "grey15") +
    scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = COL_HI), guide = "none") +
    scale_size_manual(values = c(`FALSE` = 1.2, `TRUE` = 2.6), guide = "none") +
    labs(x = sprintf("wild isotype, ranked (n = %d)", nrow(d)),
         y = "pooled response (vst)",
         title = paste0("*", g, "* RNAi")) +
    theme_classic(base_size = 11.5) +
    theme(plot.title = element_markdown(face = "bold", size = 11.5),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          plot.title.position = "plot")
}

fig <- (rank_panel("mig-6") | rank_panel("pos-1")) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 14)))

ggsave(file.path(OUT, "SUPP_FIG_XX_pooled_phenotype_ranks.pdf"), fig,
       width = 10, height = 4.4, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_pooled_phenotype_ranks.png"), fig,
       width = 10, height = 4.4, dpi = 300, bg = "white")
cat("wrote SUPP_FIG_XX_pooled_phenotype_ranks.{pdf,png}\n")

## ---------------------------------------------------------------------------
cat("\n== cross parents on the vst scale (rank 1 = most resistant) ==\n")
print(as.data.frame(ph %>% group_by(gene) %>%
  mutate(rank = rank(-vst), n = n()) %>% ungroup() %>%
  filter(strain %in% PARENTS) %>%
  transmute(gene, strain, vst = round(vst, 3), rank, n) %>%
  arrange(gene, rank)), row.names = FALSE)
