## Figure 4C, charge version ---------------------------------------------------
##
##   Rscript claude_science/scripts/03_figure4C_charge.R -> claude_science/plots/Figure4C_charge.{pdf,png}
##
## Requires claude_science/scripts/02_render_charge.py to have been run first (it writes the
## two cartoon rasters and reads claude_science/data/sid2_local_charge.tsv).
##
## THE CLAIM THIS PANEL MAKES
## Not "T96 is near the published residues" -- it is not, and with the corrected
## three-histidine set of McEwan et al. 2012 the proximity statistic is
## p = 0.38. Instead:
##
##   SID-2 uptake is electrostatic and works at gut-lumen pH. T96 sits in the
##   most positive solvent-exposed pocket of a domain that is close to neutral
##   overall at that pH, and T96K adds a further permanent positive charge
##   there. McEwan et al. showed that swapping the pH-dependent histidines for
##   permanent arginines INCREASED dsRNA uptake; T96K increases RNAi
##   sensitivity in all three backgrounds. Same direction, same physics.
##
##   i    ectodomain cartoon coloured by local (12 A) net charge at pH 4.4
##   ii   the T96 pocket: T96 with K93 and K132, the two lysines that make it
##   iii  where that pocket sits in the domain-wide distribution, and where
##        T96K moves it
##
## Charges are Henderson-Hasselbalch side-chain charges using the same pKa set
## as the manuscript's electrostatics supplement.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggtext)
})

CHG     <- "supplemental_data/structure/sid2_local_charge.tsv"
OUT     <- "plots"

COL_T96   <- "#F34C00"
COL_BASIC <- "#0B4F9E"
COL_HIS   <- "#1B7F79"
QLIM      <- 2                     # matches the renderer's colour saturation
FOCAL     <- 96
HIS_IDS   <- c(32, 168, 175)

panel_title <- function(letter)
  paste0("<span style='font-size:13pt;color:#111111'>**", letter, "**</span>")

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

chg <- read_tsv(CHG, show_col_types = FALSE)
q96 <- chg$q_local_pH44[chg$resid == FOCAL]
q96k <- q96 + 1
pct <- function(v) 100 * mean(chg$q_local_pH44 < v)
## "%.0fth" gives "82th". Ordinal suffixes, with the 11/12/13 exceptions.
ord <- function(n) {
  n <- round(n)
  suff <- if (n %% 100 %in% 11:13) "th"
          else switch(as.character(n %% 10), "1" = "st", "2" = "nd",
                      "3" = "rd", "th")
  paste0(n, suff)
}

## ---- iii. where the pocket sits, and where T96K moves it ------------------
his <- chg |> filter(resid %in% HIS_IDS) |>
  mutate(label = paste0("H", resid))

p_iii <- ggplot(chg, aes(q_local_pH44)) +
  geom_histogram(binwidth = 0.25, boundary = 0, fill = "grey78",
                 colour = "white", linewidth = 0.25) +
  geom_point(data = his, aes(q_local_pH44, 21.5), inherit.aes = FALSE,
             size = 2.2, shape = 21, fill = "white", colour = COL_HIS,
             stroke = 0.7) +
  ggrepel::geom_text_repel(data = his, aes(q_local_pH44, 21.5, label = label),
                           inherit.aes = FALSE, size = 2.6, colour = COL_HIS,
                           fontface = "bold", direction = "x", nudge_y = 2.6,
                           seed = 3, segment.size = 0.25,
                           min.segment.length = 0) +
  geom_vline(xintercept = q96, colour = COL_T96, linewidth = 0.45) +
  geom_vline(xintercept = q96k, colour = COL_BASIC, linewidth = 0.45,
             linetype = "22") +
  annotate("segment", x = q96, xend = q96k, y = 15.5, yend = 15.5,
           colour = COL_BASIC, linewidth = 0.4,
           arrow = arrow(length = grid::unit(5, "pt"), type = "closed")) +
  annotate("text", x = (q96 + q96k) / 2, y = 17.2, size = 2.7,
           colour = COL_BASIC, fontface = "bold", label = "T96K") +
  annotate("text", x = q96 - 0.12, y = 11.5, hjust = 1, size = 2.7,
           colour = COL_T96,
           label = sprintf("T96  %+.2f e\n%s pct", q96, ord(pct(q96)))) +
  annotate("text", x = q96k + 0.12, y = 11.5, hjust = 0, size = 2.7,
           colour = COL_BASIC,
           label = sprintf("T96K  %+.2f e\n%s pct", q96k, ord(pct(q96k)))) +
  ## The extra room on the right is for the T96K annotation, and it is made
  ## with coord_cartesian rather than scale limits ON PURPOSE. Setting
  ## limits = c(NA, 4.1) on the scale extends the BINNING range too, adding six
  ## empty bins past the data and emitting "Removed 1 row containing missing
  ## values (geom_bar)". No residue was ever lost -- all 168 are counted either
  ## way -- but a warning that reads like data loss in a figure script is worth
  ## not shipping. coord_cartesian zooms the view and leaves the bins alone.
  scale_x_continuous("Net charge within 12 \u00c5 at pH 4.4 (e)",
                     breaks = seq(-3, 3, 1),
                     expand = expansion(c(0.02, 0))) +
  coord_cartesian(xlim = c(NA, 4.1)) +
  scale_y_continuous("Ectodomain residues", expand = expansion(c(0, 0.12))) +
  labs(subtitle = paste0(
    "Permanent positive charge at this surface increases uptake: the<br>",
    "triple His&rarr;Arg mutant of McEwan *et al.* internalised more dsRNA<br>",
    "than wild type, and T96K increases RNAi sensitivity in all three backgrounds")) +
  theme_pub() +
  theme(axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        plot.subtitle = element_markdown(size = 8.3, colour = "grey30"))


ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_local_charge.pdf"), p_iii,
       width = 7.4, height = 4.6, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_local_charge.png"), p_iii,
       width = 7.4, height = 4.6, dpi = 300, bg = "white")
cat("wrote SUPP_FIG_XX_sid2_local_charge.{pdf,png}\n")
