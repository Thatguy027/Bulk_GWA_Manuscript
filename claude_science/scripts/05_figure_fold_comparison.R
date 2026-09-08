## Is the SID-2 ectodomain a SID-1 beta-strand-rich domain? -----------------
##
##   Rscript claude_science/scripts/05_figure_fold_comparison.R -> claude_science/plots/sid2_vs_sid1_fold.{pdf,png}
##
## TM-align of the SID-2 ectodomain model (well-modelled core, pLDDT >= 70,
## 117 residues) against the two extracellular beta-strand-rich domains of
## cSID1 (PDB 8XC1, Wang et al. 2024 NAR, doi:10.1093/nar/gkae395), with
## controls. TM-score is normalised by the SID-2 query throughout, so the bars
## are comparable; 0.5 is the conventional same-fold threshold.
##
## The comparison that matters is not SID-2 vs BRD1 in isolation but SID-2 vs
## BRD1 against (a) what an unrelated beta-sandwich scores and (b) what the two
## genuine BRDs score against each other. Both are on the plot.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggtext)
})

TAB <- "claude_science/data/sid2_vs_sid1_tmalign.csv"
OUT <- "claude_science/plots"

ROLE_COL <- c(`internal reference` = "#37474F",
              `query vs SID-1 domain` = "#F34C00",
              `negative control` = "grey72")

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_markdown(size = base_size + 0.5),
          plot.subtitle = element_markdown(size = base_size - 2.5,
                                           colour = "grey30"),
          plot.title.position = "plot",
          legend.key.size = grid::unit(8, "pt"))
}

d <- read_csv(TAB, show_col_types = FALSE) |>
  mutate(role = factor(role, levels = names(ROLE_COL)),
         target = fct_reorder(target, tm_norm_sid2))

p <- ggplot(d, aes(tm_norm_sid2, target, fill = role)) +
  annotate("rect", xmin = 0.5, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#37474F", alpha = 0.05) +
  geom_col(width = 0.66) +
  geom_vline(xintercept = 0.5, linewidth = 0.4, colour = "#37474F",
             linetype = "22") +
  geom_text(aes(label = sprintf("%.3f", tm_norm_sid2)), hjust = -0.18,
            size = 2.9, colour = "grey25") +
  annotate("text", x = 0.505, y = 0.62, hjust = 0, size = 2.7,
           colour = "#37474F", label = "same-fold threshold") +
  scale_fill_manual(values = ROLE_COL, name = NULL) +
  scale_x_continuous("TM-score, normalised by the SID-2 query",
                     limits = c(0, 0.68), breaks = seq(0, 0.6, 0.1),
                     expand = expansion(c(0, 0))) +
  labs(y = NULL,
       title = "The SID-2 ectodomain does not match a SID-1 &beta;-strand-rich domain",
       subtitle = paste0("An unrelated immunoglobulin domain scores higher against SID-2 than BRD1 does, ",
                         "while the two<br>genuine BRDs of cSID1 score 0.525 against each other. ",
                         "Sequence identity across every<br>alignment is 1.6-10.4%, i.e. background.")) +
  theme_pub() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7.6),
        legend.margin = margin(t = -4),
        axis.title.x = element_text(size = 9),
        axis.text.y = element_text(size = 8.4))

ggsave(file.path(OUT, "sid2_vs_sid1_fold.pdf"), p, width = 8.2, height = 4.2,
       device = cairo_pdf)
ggsave(file.path(OUT, "sid2_vs_sid1_fold.png"), p, width = 8.2, height = 4.2,
       dpi = 300, bg = "white")
cat("rows:", nrow(d), "| max TM:", max(d$tm_norm_sid2), "\n")
