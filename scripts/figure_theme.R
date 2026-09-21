## One theme and one type scale for the manuscript figures -------------------
##
##   source("scripts/figure_theme.R")
##
## Every main display figure takes its theme from here, so text is the same
## size in Figure 1 as in Figure 4. Before this there were twenty-two separate
## copies of theme_pub() across scripts/, at base sizes 10, 11, 11.5 and 12,
## and the main figures alone split 11.5 (Figures 1-3) against 11 (Figure 4).
##
## THE UNIT TRAP THIS EXISTS TO CLOSE. ggplot sizes theme text in POINTS and
## geom_text/annotate text in MILLIMETRES, and the two are silently four-tenths
## apart: size = 2.2 in an annotate() is 6.3 pt, not 6.3 mm and not 2.2 pt.
## That is how Figure 4C's colourbar labels ended up at 6.3 pt beside axis text
## at 9.2 pt in the same figure. Annotation sizes are therefore written as
## pt_mm(9.2) or one of the named constants below, never as a bare number.
##
## A panel that draws its own axis -- Figure 3C puts two scales side by side in
## one coordinate system, so it has to -- uses TXT_AXIS for the tick labels and
## TXT_TITLE for the axis title. That is what makes a hand-drawn axis the same
## size as a real one in the panel beside it.

BASE_SIZE <- 11.5

## millimetres for a given point size -- the conversion geom_text wants
pt_mm <- function(pt) pt / (72.27 / 25.4)

TXT_TITLE <- pt_mm(BASE_SIZE)         # matches axis.title
TXT_AXIS  <- pt_mm(BASE_SIZE * 0.8)   # matches axis.text, theme_classic's rel(0.8)
TXT_NOTE  <- pt_mm(BASE_SIZE * 0.7)   # in-panel annotations that sit below the axis
TXT_SMALL <- pt_mm(BASE_SIZE * 0.6)   # crowded labels: residue names, tick rows

theme_pub <- function(base_size = BASE_SIZE) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold", size = base_size),
      panel.spacing.x  = grid::unit(8, "pt"),
      axis.line        = ggplot2::element_line(linewidth = 0.3),
      axis.ticks       = ggplot2::element_line(linewidth = 0.3),
      plot.title       = ggtext::element_markdown(size = base_size),
      plot.subtitle    = ggtext::element_markdown(size = base_size - 3,
                                                  colour = "grey30"),
      plot.title.position = "plot",
      legend.key.size  = grid::unit(9, "pt"))
}
