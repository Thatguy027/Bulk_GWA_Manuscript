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

CHG     <- "claude_science/data/sid2_local_charge.tsv"
CARTOON <- "claude_science/plots/sid2_cartoon_charge.png"
ZOOM    <- "claude_science/plots/sid2_zoom_charge.png"
OUT     <- "claude_science/plots"

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

## ---- i. the charge map ----------------------------------------------------
im <- png::readPNG(CARTOON)
H_IM <- 4.05
W_IM <- H_IM * (dim(im)[2] / dim(im)[1])

## diverging key drawn as a gradient strip, same RdBu ramp as the renderer
## NOT reversed: matplotlib's RdBu runs red at the low end to blue at the high
## end, and the renderer maps charge through Normalize(-QLIM, +QLIM), so red is
## negative and blue is positive. Reversing here silently inverted the key.
ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(64)
key <- tibble(i = seq_along(ramp),
              xmin = (i - 1) / length(ramp) * W_IM * 0.72,
              xmax = i / length(ramp) * W_IM * 0.72,
              col = ramp)
KEY_Y <- -0.34; KEY_H <- 0.15

p_i <- ggplot() +
  annotation_raster(im, xmin = 0, xmax = W_IM, ymin = 0, ymax = H_IM,
                    interpolate = TRUE) +
  geom_rect(data = key, aes(xmin = xmin, xmax = xmax, ymin = KEY_Y,
                            ymax = KEY_Y + KEY_H), fill = key$col) +
  annotate("text", x = 0, y = KEY_Y - 0.10, hjust = 0, size = 2.3,
           colour = "grey30", label = paste0("\u2212", QLIM)) +
  annotate("text", x = W_IM * 0.72, y = KEY_Y - 0.10, hjust = 1, size = 2.3,
           colour = "grey30", label = paste0("+", QLIM)) +
  annotate("text", x = 0, y = KEY_Y + KEY_H + 0.16, hjust = 0, size = 2.5,
           colour = "grey30",
           label = "net charge within 12 \u00c5 (e), pH 4.4") +
  coord_fixed(ratio = 1, xlim = c(-0.05, W_IM + 0.05),
              ylim = c(KEY_Y - 0.24, H_IM), expand = FALSE, clip = "off") +
  labs(title = panel_title("C"),
       subtitle = paste0("Local charge at gut-lumen pH; T96<br>",
                         "sits in the most positive lumenal pocket")) +
  theme_void(base_size = 11) +
  theme(plot.title = element_markdown(size = 11.5),
        plot.subtitle = element_markdown(size = 8.3, colour = "grey30"),
        plot.title.position = "plot",
        plot.margin = margin(2, 6, 2, 6))

## ---- ii. the pocket -------------------------------------------------------
im_z <- png::readPNG(ZOOM)
W_Z <- H_IM * (dim(im_z)[2] / dim(im_z)[1])

p_zoom <- ggplot() +
  annotation_raster(im_z, xmin = 0, xmax = W_Z, ymin = 0, ymax = H_IM,
                    interpolate = TRUE) +
  coord_fixed(ratio = 1, xlim = c(-0.05, W_Z + 0.05), ylim = c(0, H_IM),
              expand = FALSE, clip = "off") +
  labs(subtitle = paste0("The pocket: K93 and K132 lie 6.6 and 6.8 &#197; from<br>",
                         "T96 (4.5 and 4.4 &#197; heavy-atom)")) +
  theme_void(base_size = 11) +
  theme(plot.subtitle = element_markdown(size = 8.3, colour = "grey30"),
        plot.title.position = "plot",
        plot.margin = margin(2, 6, 2, 6))

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
           label = sprintf("T96  %+.2f e\n%.0fth pct", q96, pct(q96))) +
  annotate("text", x = q96k + 0.12, y = 11.5, hjust = 0, size = 2.7,
           colour = COL_BASIC,
           label = sprintf("T96K  %+.2f e\n%.0fth pct", q96k, pct(q96k))) +
  scale_x_continuous("Net charge within 12 \u00c5 at pH 4.4 (e)",
                     breaks = seq(-3, 3, 1), limits = c(NA, 4.1),
                     expand = expansion(c(0.02, 0))) +
  scale_y_continuous("Ectodomain residues", expand = expansion(c(0, 0.12))) +
  labs(subtitle = paste0(
    "Permanent positive charge at this surface increases uptake: the<br>",
    "triple His&rarr;Arg mutant of McEwan *et al.* internalised more dsRNA<br>",
    "than wild type, and T96K increases RNAi sensitivity in all three backgrounds")) +
  theme_pub() +
  theme(axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        plot.subtitle = element_markdown(size = 8.3, colour = "grey30"))

fig <- (p_i | p_zoom | p_iii) + plot_layout(widths = c(0.60, 0.78, 1.30))

dir.create(OUT, showWarnings = FALSE)
ggsave(file.path(OUT, "Figure4C_charge.pdf"), fig, width = 12.4, height = 5.0,
       device = cairo_pdf)
ggsave(file.path(OUT, "Figure4C_charge.png"), fig, width = 12.4, height = 5.0,
       dpi = 300, bg = "white")

cat(sprintf("T96 %+.2f e (%.0fth pct) -> T96K %+.2f e (%.0fth pct); domain median %+.2f\n",
            q96, pct(q96), q96k, pct(q96k), median(chg$q_local_pH44)))
cat(sprintf("ectodomain residues with positive local charge: %d of %d\n",
            sum(chg$q_local_pH44 > 0), nrow(chg)))
