## The sid-2 variant panel -- shared by Figure 4D and its full-catalogue supplement
##
## Sourced, not run. Figure4_sid2.R draws the sites that differ among the three
## parents it names; SUPP_FIG_XX_sid2_variants_all.R draws every protein-altering
## site in the latest CeNDR release. The drawing is one function so the two cannot
## drift apart in geometry, colour or the charge strip.
##
## The protein runs vertically with residue 1 at the top, the topology as a
## filled bar, and each variant labelled to the right with a bar for its
## population frequency (or count). Callout rows are evenly spaced down the
## panel and joined to the residue by a leader, rather than sitting at the
## residue's own height: variants cluster between residues 141 and 153 and their
## labels overlapped completely at true scale.
##
## `vr` needs: residue, label, focal, ju1793_aa, ju2466_aa, xz1516_aa,
## parents_differ, and the column named by `bar`. Callers must have sourced
## figure_palette.R and figure_theme.R and defined TOPO_COL2 and panel_title().
## ---------------------------------------------------------------------------

SID2_PERRES <- "supplemental_data/structure/sid2_per_residue.tsv"
SID2_LOCALQ <- "supplemental_data/structure/sid2_local_charge.tsv"

sid2_variant_panel <- function(vr, letter = "D", base_size, ramp, qlim,
                               bar = "af", bar_max = 1,
                               bar_text = function(x) sprintf("%.0f%%", 100 * x),
                               bar_header = "CeNDR frequency",
                               band_half = 11) {
  pr2 <- read_tsv(SID2_PERRES, show_col_types = FALSE) %>%
    mutate(topology = factor(topology, levels = names(TOPO_COL2)))
  LEN <- max(pr2$resid)

  dom2 <- (function(v) { r <- rle(as.character(v))
    tibble(value = factor(r$values, levels = names(TOPO_COL2)),
           end = cumsum(r$lengths),
           start = cumsum(r$lengths) - r$lengths + 1) })(pr2$topology)

  ## x geometry, in arbitrary units: the topology bar, then the labels, then the
  ## frequency bars. Residue is on y, reversed so residue 1 is at the top.
  X_BAR <- c(0, 0.5); X_LAB <- 0.78; X_FRQ <- c(1.36, 2.55)
  ## A local-net-charge strip immediately left of the topology bar, on the SAME
  ## ramp and the SAME limits as Figure 4C, so the two can be read against each
  ## other. It needs no key of its own; panel C's key serves both.
  X_CHG <- c(-0.86, -0.16)
  ## Colours are precomputed to hex rather than mapped through a second fill
  ## scale: the panel already uses fill for the topology, and ggnewscale is not a
  ## dependency of this repository.
  chg <- read_tsv(SID2_LOCALQ, show_col_types = FALSE) %>%
    transmute(resid, q = q_local_pH44,
              col = ramp[pmax(1, pmin(length(ramp),
                       round((pmin(pmax(q, -qlim), qlim) + qlim) /
                             (2 * qlim) * (length(ramp) - 1)) + 1))])
  ## The charge is defined only where there is a model to measure it in: the
  ## AlphaFold ectodomain, residues 21-188. The rest is left blank rather than
  ## filled with a zero that would read as "neutral here".
  CHG_RANGE <- range(chg$resid)
  stopifnot(nrow(chg) == 168, CHG_RANGE[1] == 21, CHG_RANGE[2] == 188)
  ## Three allele columns. X_P1/X_P2 are the parents of the JU1793 x JU2466
  ## cross; X_P3 is XZ1516, the parent of the OTHER mapping cross, carried for
  ## reference. It sits outside the highlight band on purpose -- the band marks
  ## sites that segregate in JU1793 x JU2466, and XZ1516 is not in that cross.
  X_P1 <- 3.16; X_P2 <- 3.62; X_P3 <- 4.08
  vr <- vr %>% arrange(residue) %>%
    mutate(row = seq(14, LEN - 14, length.out = n()),
           value = .data[[bar]],
           xend = X_FRQ[1] + value / bar_max * diff(X_FRQ))
  ## the band may not be taller than the row pitch, or adjacent bands merge
  if (nrow(vr) > 1) band_half <- min(band_half, 0.45 * diff(vr$row[1:2]))

  ggplot() +
    ## the charge strip, and an outline showing how far the model reaches.
    ## height slightly over 1 so adjacent residues abut with no hairline seam:
    ## a white seam would read as neutral charge or as missing data
    geom_tile(data = chg, aes(x = mean(X_CHG), y = resid),
              fill = chg$col, width = diff(X_CHG), height = 1.02) +
    annotate("rect", xmin = X_CHG[1], xmax = X_CHG[2],
             ymin = CHG_RANGE[1] - 0.5, ymax = CHG_RANGE[2] + 0.5,
             fill = NA, colour = "grey45", linewidth = 0.25) +

    geom_rect(data = dom2,
              aes(xmin = X_BAR[1], xmax = X_BAR[2],
                  ymin = start - 0.5, ymax = end + 0.5, fill = value),
              colour = "grey30", linewidth = 0.25) +
    ## tick across the topology bar at the true residue, then a leader out to
    ## the evenly spaced callout row
    geom_segment(data = vr,
                 aes(x = X_BAR[1], xend = X_BAR[2], y = residue, yend = residue),
                 linewidth = 0.5,
                 colour = ifelse(vr$focal, COL_JU1793, "grey35")) +
    geom_segment(data = vr,
                 aes(x = X_BAR[2], xend = X_LAB - 0.04, y = residue, yend = row),
                 linewidth = 0.3,
                 colour = ifelse(vr$focal, COL_JU1793, "grey55")) +
    geom_richtext(data = vr, aes(X_LAB, row, label = label),
                  colour = ifelse(vr$focal, COL_JU1793, "grey15"),
                  fontface = ifelse(vr$focal, "bold", "plain"),
                  size = 2.9, hjust = 0, vjust = 0.5, fill = NA,
                  label.color = NA,
                  label.padding = grid::unit(rep(0, 4), "pt")) +
    ## frequency bars
    geom_segment(data = vr, aes(x = X_FRQ[1], xend = X_FRQ[2], y = row,
                                yend = row),
                 linewidth = 3.0, colour = "grey92", lineend = "butt") +
    geom_segment(data = vr, aes(x = X_FRQ[1], xend = xend, y = row,
                                yend = row),
                 linewidth = 3.0, lineend = "butt",
                 colour = ifelse(vr$focal, COL_JU1793, "grey55")) +
    geom_richtext(data = vr, aes(X_FRQ[2] + 0.06, row,
                                 label = bar_text(value)),
                  colour = ifelse(vr$focal, COL_JU1793, "grey25"),
                  size = 2.7, hjust = 0, vjust = 0.5, fill = NA,
                  label.color = NA,
                  label.padding = grid::unit(rep(0, 4), "pt")) +
    ## the allele each strain carries. Rows where the two CROSS PARENTS differ
    ## are banded; XZ1516 is drawn plain, outside the band, because it is
    ## reference rather than a term in that comparison.
    geom_rect(data = vr %>% filter(parents_differ),
              aes(xmin = X_P1 - 0.22, xmax = X_P2 + 0.22,
                  ymin = row - band_half, ymax = row + band_half),
              fill = COL_JU1793, alpha = 0.10) +
    geom_richtext(data = vr, aes(X_P1, row, label = ju1793_aa),
                  colour = ifelse(vr$parents_differ, COL_JU1793, "grey45"),
                  fontface = ifelse(vr$parents_differ, "bold", "plain"),
                  size = 2.9, hjust = 0.5, vjust = 0.5, fill = NA,
                  label.color = NA,
                  label.padding = grid::unit(rep(0, 4), "pt")) +
    geom_richtext(data = vr, aes(X_P2, row, label = ju2466_aa),
                  colour = ifelse(vr$parents_differ, COL_JU2466, "grey45"),
                  fontface = ifelse(vr$parents_differ, "bold", "plain"),
                  size = 2.9, hjust = 0.5, vjust = 0.5, fill = NA,
                  label.color = NA,
                  label.padding = grid::unit(rep(0, 4), "pt")) +
    geom_richtext(data = vr, aes(X_P3, row, label = xz1516_aa),
                  colour = "grey45", size = 2.9, hjust = 0.5, vjust = 0.5,
                  fill = NA, label.color = NA,
                  label.padding = grid::unit(rep(0, 4), "pt")) +
    ## column headers, in the margin reserved above residue 1
    geom_richtext(data = tibble(
        x = c(mean(X_CHG), mean(X_FRQ), X_P1, X_P2, X_P3), y = -7,
        lab = c("Net charge", bar_header, "JU1793", "JU2466", "XZ1516")),
        aes(x, y, label = lab),
        colour = c("grey30", "grey30", COL_JU1793, COL_JU2466, COL_XZ),
        size = 2.6, hjust = 0.5, vjust = 0.5, fill = NA, label.color = NA,
        label.padding = grid::unit(rep(0, 4), "pt")) +
    ## The ends of the charge column, so it reads as a scaled quantity. geom_text,
    ## NOT geom_richtext: gridtext parses a bare "+" as a markdown list item.
    geom_text(data = tibble(x = X_CHG, y = c(-2.5, -2.5),
                            lab = c("−", "+")),
              aes(x, y, label = lab), size = 2.6, colour = "grey45",
              hjust = 0.5, vjust = 0.5) +
    scale_fill_manual(values = TOPO_COL2, name = NULL) +
    scale_x_continuous(limits = c(X_CHG[1] - 0.04, X_P3 + 0.30),
                       expand = expansion(0)) +
    ## 311 dropped from the breaks: it collided with the 300 tick
    scale_y_reverse(breaks = c(1, seq(50, 300, 50)),
                    limits = c(LEN + 6, -16), expand = expansion(0)) +
    labs(x = NULL, y = "SID-2 residue", title = panel_title(letter)) +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
    theme_pub(base_size) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          legend.position = "bottom", legend.margin = margin(t = -4),
          legend.text = element_text(size = 7.6),
          plot.margin = margin(4, 4, 4, 4))
}
