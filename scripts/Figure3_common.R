## Figure 3 -- shared vocabulary and panel builders --------------------------
## Sourced by the three Figure 3 variants. Not run on its own.
##
##   scripts/Figure3.R          genome-wide cross scan in panel A
##   scripts/Figure3_chrIII.R   chromosome III only in panel A
##   scripts/Figure3_pheno.R    genome-wide, with the pooled pos-1 phenotype
##                              distribution inset over the empty left half
##
## The three differ only in panel A, so panels B and C, the colours, the row
## geometry and the theme live here once. Editing a colour or a row order here
## changes all three, which is the point: they are three renderings of one
## figure, not three figures.
##
## Fine-mapping the chromosome III QTL with near-isogenic lines.
##   A  the JU1793 x JU2466 HT115-vs-pos-1 cross scan that defines the QTL
##   B  the NIL introgressions on the right arm of chromosome III
##   C  embryo hatching under pos-1 RNAi for the same strains
##
## Styled to match Figure 2: theme_classic, chromosome facets pinned to the
## WBcel235 lengths, LOD rather than -log10 p for a cross scan, cairo output,
## gene names in italics, panel letters only. Panel A carries no threshold
## line, no LOD value and no peak marker.
##
## Notes of substance
## ------------------
## 1  Panel A is read from data/cross_experiments/JU1793-JU2466_export, not the
##    legacy Nov2024 plots directory, and plots LOD. The legacy file's `p`
##    column is formed on the linear scale and underflows to exactly 0 above
##    |z| ~ 38.5; it happens not to bite for this contrast (peak LOD 140, |z|
##    ~ 25) but -log10(p) is not safe to use from that file in general. LOD is
##    computed in log space and cannot saturate. The threshold is unchanged in
##    meaning: alpha = 0.05 over 2000 effective tests is -log10 p = 4.60, which
##    is LOD 3.57 in the package convention the exports use. It is reported to
##    the console and deliberately not drawn.
##
## 2  Strain order and the JU1793 / JU2466 colours are exactly as previously
##    defined for this figure, and neither B nor C carries strain labels, per
##    the existing convention. Rows run JU1793 at the bottom to JU2466 at the
##    top in both panels.
##
## 3  The mapped region, 13.658-13.695 Mb, is shaded in panel B. It is 37 kb,
##    which is sub-pixel at genome scale, so it is NOT marked in the
##    genome-wide panel A -- a mark that small reads as noise or, worse, as
##    part of the chrIII peak. It is marked in the chromosome III variant,
##    where 37 kb is 0.3% of the axis and can be placed honestly.
##
## Panel C shows pos-1 only. The HT115 control arm, and the four further NILs
## measured in the same experiment, are in SUPP_FIG_XX_nil_hatching_full.R.
##
## RNAI DOSE: 50% pos-1 RNAi bacteria, diluted with HT115. The dose is NOT
## recorded in the source data file -- the condition column carries only
## "pos"/"ht115" -- so it comes from the lab record. Figure 4B uses 25% and its
## hatching percentages are therefore not comparable with these; see
## FIGURE_CAPTIONS.txt.
##
## CAVEAT: one plate per strain per condition. The intervals in panel C are
## Wilson binomial intervals on that single plate's embryo count, so they
## describe counting uncertainty, not between-plate variability.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(ggtext)
})

OUT  <- "plots"
NILB <- "supplemental_data/hatching_assays/nil_introgression_ranges.bed"
PHEN <- "supplemental_data/hatching_assays/nil_series_hatching.tsv"
SCAN <- paste0("supplemental_data/mapping/",
               "ju_cross_ht115_vs_pos1_scan.tsv.gz")
POOL <- "supplemental_data/phenotypes/pooled_vst_traits.csv.gz"

## the strain colours as defined for this figure; do not re-map these
COL_JU1793 <- "#F34C00"
COL_JU2466 <- "#40B4AB"
COL_THR    <- "grey50"
COL_REGION <- "#4D4D4D"   # the interval the NIL series resolves

CHROMS  <- c("I", "II", "III", "IV", "V", "X")
ALL_LEN <- c(I = 15072434, II = 15279421, III = 13783801,
             IV = 17493829, V = 20924180, X = 17718942)
fct_chr <- function(x) factor(as.character(x), levels = CHROMS)

CROSS_THR <- stats::qchisq(2 * 0.05 / 2000, df = 1, lower.tail = FALSE) / (2 * log(10))

WIN      <- c(13.60e6, ALL_LEN[["III"]])          # panel B window
RESOLVED <- tibble::tibble(xmin = 13657700, xmax = 13695000)

## strain order as originally defined
LEVELS <- c("JU1793", "wSZ196", "wSZ191", "wSZ176", "JU2466")

## B and C sit side by side and must share row positions exactly, so both use a
## numeric y with the same limits. Levels map straight to y, so JU1793 is row 1
## at the bottom and JU2466 row 5 at the top.
ROW   <- setNames(seq_along(LEVELS), LEVELS)
Y_LIM <- c(0.45, length(LEVELS) + 1.05)
y_axis <- function(labels = TRUE) list(
  scale_y_continuous(limits = Y_LIM, breaks = seq_along(LEVELS),
                     labels = if (labels) LEVELS else NULL,
                     expand = expansion(mult = 0)))

## panel letters only; the panels are described in the caption, not on the plot
panel_title <- function(letter) {
  paste0("<span style='font-size:14pt;color:#111111'>**", letter, "**</span>")
}

theme_pub <- function(base_size = 11.5) {
  theme_classic(base_size = base_size) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = base_size),
          ## 4pt was enough on an 11in canvas; at 9.6in the last tick label of
          ## one chromosome collides with the first of the next
          panel.spacing.x = grid::unit(8, "pt"),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_markdown(size = base_size),
          plot.title.position = "plot",
          legend.key.size = grid::unit(9, "pt"))
}

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## ===========================================================================
## the cross scan, shared by all three panel A variants
## ===========================================================================
load_scan <- function() {
  s <- fread(SCAN, select = c("chrom", "physical.position", "LOD")) %>%
    as_tibble() %>% filter(chrom %in% CHROMS, is.finite(LOD)) %>%
    transmute(chrom = fct_chr(chrom), pos = physical.position,
              pos.mb = physical.position / 1e6, LOD)
  pk <- s %>% filter(chrom == "III") %>% slice_max(LOD, n = 1)
  msg("  chrIII peak: ", round(pk$pos.mb, 3), " Mb, LOD ", round(pk$LOD, 1),
      "  (genome-wide threshold LOD ", round(CROSS_THR, 2), ", not drawn)")
  s
}

## peak-preserving thinning: bin, keep the maximum. Preserves peak heights and
## positions, which slice-by-row or every-nth thinning does not.
thin_scan <- function(s, bin = 5e3) {
  s %>% mutate(bin = floor(pos / bin)) %>%
    group_by(chrom, bin) %>% slice_max(LOD, n = 1, with_ties = FALSE) %>%
    ungroup()
}

## ---------------------------------------------------------------------------
## A, genome-wide. facet_grid(space = "free_x") would otherwise size each panel
## by its own data range, so every panel is pinned to the WBcel235 length.
## ---------------------------------------------------------------------------
panel_A_genome <- function(scan_t, letter = "A") {
  span <- bind_rows(
    tibble(chrom = fct_chr(names(ALL_LEN)), pos.mb = unname(ALL_LEN) / 1e6),
    tibble(chrom = fct_chr(names(ALL_LEN)), pos.mb = 0))

  ggplot(scan_t, aes(pos.mb, LOD)) +
    geom_blank(data = span, aes(x = pos.mb, y = 0), inherit.aes = FALSE) +
    geom_line(linewidth = 0.4, colour = "grey15") +
    facet_grid(. ~ chrom, scales = "free_x", space = "free_x") +
    scale_x_continuous(breaks = seq(0, 25, 5), expand = expansion(mult = 0.02),
                       guide = guide_axis(check.overlap = TRUE)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.07))) +
    labs(x = "Genomic Position (Mb)", y = "LOD", title = panel_title(letter)) +
    theme_pub() +
    theme(axis.text.x = element_text(size = 8.5))
}

## ---------------------------------------------------------------------------
## A, chromosome III only. Here 37 kb is 0.3% of the axis, so the mapped region
## can be drawn as a band with dotted edges rather than left unmarked.
##
## The panel is deliberately squat. This QTL sits at the right telomere, so the
## trace is flat across 12 of 13.8 Mb; a panel tall enough to be square is
## almost entirely empty. A short panel reads as a baseline with a rise at the
## end, which is what the data are. An earlier version bracketed the panel B
## window under the axis, but the window is the last 1.3% of the chromosome and
## the bracket collided with the axis text -- panel B's own axis labels make
## the zoom obvious without it.
## ---------------------------------------------------------------------------
panel_A_chr3 <- function(scan_t, letter = "A", from_mb = 0, strip = FALSE) {
  s3  <- scan_t %>% filter(chrom == "III")
  lim <- c(from_mb, ALL_LEN[["III"]] / 1e6)

  ggplot(s3, aes(pos.mb, LOD)) +
    geom_rect(data = RESOLVED, inherit.aes = FALSE,
              aes(xmin = xmin / 1e6, xmax = xmax / 1e6, ymin = -Inf, ymax = Inf),
              fill = COL_REGION, alpha = 0.18) +
    geom_vline(xintercept = c(RESOLVED$xmin, RESOLVED$xmax) / 1e6,
               linetype = "dotted", linewidth = 0.4, colour = COL_REGION) +
    geom_line(linewidth = 0.45, colour = "grey15") +
    facet_grid(. ~ chrom) +
    scale_x_continuous(breaks = scales::pretty_breaks(6),
                       expand = expansion(mult = 0.01)) +
    scale_y_continuous(breaks = scales::pretty_breaks(4),
                       expand = expansion(mult = c(0.02, 0.08))) +
    ## clip rather than filter, so the line is drawn from the full scan and
    ## does not appear to start from zero at the left edge of a cut window
    coord_cartesian(xlim = lim) +
    labs(x = "Chromosome III (Mb)", y = "LOD", title = panel_title(letter)) +
    theme_pub() +
    ## the strip would read "III" beside an axis already titled
    ## "Chromosome III (Mb)"; element_blank collapses the strip row entirely
    theme(strip.text = if (strip) NULL else element_blank())
}

## ---------------------------------------------------------------------------
## the pooled pos-1 phenotype, for the inset in the third variant
##
## This is the VST-transformed pos-1 response from the pooled panel, the same
## trait that was mapped by GWAS in Figure 2. Positive = the strain gained
## frequency under pos-1 RNAi = resistant. Both cross parents were in the
## pooled panel, at opposite ends of it, which is why they were chosen as
## parents; the inset makes that selection visible rather than asserted.
## ---------------------------------------------------------------------------
pooled_pos1 <- function() {
  read_csv(POOL, show_col_types = FALSE) %>%
    transmute(strain, v = `vst_ctrl_pos-1_T2`) %>%
    filter(is.finite(v))
}

pheno_inset <- function(type = c("density", "hist"), base_size = 8.5,
                        bins = 22, letter = NULL) {
  type <- match.arg(type)
  ph <- pooled_pos1()

  ## Both lollipops are drawn to one common height rather than to the height of
  ## the distribution at each strain's value. Height at the curve encodes how
  ## COMMON a phenotype is, which is not the point being made, and it made
  ## JU1793 -- the more extreme of the two parents -- the shorter marker,
  ## reading backwards. Equal heights leave position on the x axis as the only
  ## signal, which is the one that matters.
  if (type == "density") {
    d  <- stats::density(ph$v, adjust = 0.9)
    dn <- tibble(x = d$x, y = d$y)
    H  <- max(dn$y) * 0.75
    body <- list(
      geom_area(data = dn, aes(x, y), fill = "grey88", colour = NA),
      geom_line(data = dn, aes(x, y), linewidth = 0.3, colour = "grey40"),
      scale_y_continuous(expand = expansion(mult = c(0, 0.16))),
      labs(y = NULL))
    y_theme <- theme(axis.line.y = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
  } else {
    ## a histogram shows the sample directly -- 84 strains in countable bins --
    ## where a kernel density is a smooth of it, and invents the small bumps in
    ## the tails that the density version shows near -0.25 and +0.23
    bw <- diff(range(ph$v)) / bins
    H  <- max(graphics::hist(ph$v, breaks = seq(min(ph$v),
                                                max(ph$v) + bw, by = bw),
                             plot = FALSE)$counts) * 0.75
    body <- list(
      geom_histogram(data = ph, aes(x = v), binwidth = bw,
                     boundary = 0, fill = "grey85", colour = "grey40",
                     linewidth = 0.25),
      scale_y_continuous(breaks = scales::pretty_breaks(4),
                         expand = expansion(mult = c(0, 0.16))),
      labs(y = "Strains"))
    y_theme <- theme(axis.line.y = element_line(linewidth = 0.3),
                     axis.ticks.y = element_line(linewidth = 0.3))
  }

  par2 <- ph %>% filter(strain %in% c("JU1793", "JU2466")) %>%
    mutate(y = H,
           ## label away from the centre of the distribution on each side
           hj = ifelse(v > stats::median(ph$v), 0, 1),
           nx = ifelse(hj == 0, 1, -1) * diff(range(ph$v)) * 0.035)

  msg("  inset (", type, "): pooled pos-1 VST, n = ", nrow(ph),
      " | ", paste(sprintf("%s %.3f (rank %d/%d)", par2$strain, par2$v,
                           rank(ph$v)[match(par2$strain, ph$strain)], nrow(ph)),
                   collapse = " | "))

  ggplot() +
    body +
    ## every strain as a tick, so the reader sees the sample, not just its shape
    geom_rug(data = ph, aes(x = v), sides = "b", length = grid::unit(3.5, "pt"),
             linewidth = 0.25, colour = "grey45") +
    geom_segment(data = par2, aes(x = v, xend = v, y = 0, yend = y,
                                  colour = strain), linewidth = 0.7) +
    geom_point(data = par2, aes(v, y, colour = strain), size = 1.9) +
    geom_richtext(data = par2,
                  aes(x = v + nx, y = y, label = strain, colour = strain,
                      hjust = hj),
                  size = 2.7, vjust = 0.35, fill = NA, label.color = NA,
                  label.padding = grid::unit(rep(0, 4), "pt")) +
    scale_colour_manual(values = c(JU1793 = COL_JU1793, JU2466 = COL_JU2466),
                        guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0.10)) +
    labs(x = "Pooled *pos-1* response (VST)",
         ## a letter only when this stands as a panel rather than an inset
         title = if (is.null(letter)) NULL else panel_title(letter)) +
    theme_classic(base_size = base_size) +
    y_theme +
    theme(axis.line.x = element_line(linewidth = 0.3),
          axis.ticks.x = element_line(linewidth = 0.3),
          ## must name axis.title.x, not the parent axis.title: theme_classic
          ## sets axis.title.x explicitly and a parent element_markdown loses
          ## to it, which silently renders the asterisks as literal text
          axis.title.x = element_markdown(size = base_size - 0.5),
          axis.title.y = element_text(size = base_size - 0.5),
          plot.title = element_markdown(size = base_size),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(2, 4, 1, 2))
}

## ---------------------------------------------------------------------------
## layout
##
## Both layouts are FLAT -- one patchwork with an explicit design, never
## pA / (pB | pC). Nesting puts the bottom row on its own alignment level, so
## patchwork aligns the row as a block and panel A ends up wider than B and C
## together, by the difference between A's y-axis furniture and B's (which has
## none). A single design region per panel puts every panel in one grid, and
## column 1's axis gutter is then the widest of the panels in it, so the edges
## line up.
##
## C is given less width than B: the hatching bars only need to be long enough
## to read a percentage off, whereas B carries 180 kb of breakpoint detail.
## ---------------------------------------------------------------------------
three_panel <- function(pA, pB, pC, heights = c(0.72, 1), widths = c(1.45, 1)) {
  pA + pB + pC +
    plot_layout(design = c(area(1, 1, 1, 2), area(2, 1, 2, 1),
                           area(2, 2, 2, 2)),
                heights = heights, widths = widths)
}

## the four-panel arrangement: phenotype and scan on top, phenotype-by-genotype
## below. Column 1 holds the two narrow panels, column 2 the two wide ones.
four_panel <- function(pA, pB, pC, pD, heights = c(1, 1), widths = c(1, 1.45)) {
  pA + pB + pC + pD +
    plot_layout(ncol = 2, heights = heights, widths = widths)
}

## ===========================================================================
## B -- the NIL introgressions
## ===========================================================================
panel_B <- function(letter = "B") {
  bed <- fread(NILB, header = FALSE,
               col.names = c("chr", "start", "end", "strain", "geno")) %>%
    as_tibble() %>% filter(strain %in% LEVELS) %>%
    mutate(y = ROW[strain],
           other = ifelse(geno == "JU1793", "JU2466", "JU1793"),
           ## the parents are uniform across the window; the NILs carry a
           ## JU2466 introgression on a JU1793 background
           seg_start = pmax(start, WIN[1]),
           is_parent = strain %in% c("JU1793", "JU2466"))

  ## background bar spanning the whole window, then the introgression on top
  bg  <- bed %>% transmute(y, xmin = WIN[1] / 1e6, xmax = WIN[2] / 1e6,
                           geno = ifelse(is_parent, geno, other))
  seg <- bed %>% filter(!is_parent) %>%
    transmute(y, xmin = seg_start / 1e6, xmax = end / 1e6, geno)

  ## The shaded region is drawn AFTER the genotype bars, so it tints them
  ## rather than sitting behind them, and it stops at the top and bottom bar
  ## edges rather than running to the panel limits -- it used to extend up
  ## behind the coordinate label.
  BAR_H <- 0.32
  band <- RESOLVED %>%
    mutate(ymin = min(ROW) - BAR_H, ymax = max(ROW) + BAR_H)

  ggplot() +
    geom_rect(data = bg, aes(xmin = xmin, xmax = xmax, fill = geno,
                             ymin = y - BAR_H, ymax = y + BAR_H),
              colour = "grey25", linewidth = 0.3) +
    geom_rect(data = seg, aes(xmin = xmin, xmax = xmax, fill = geno,
                              ymin = y - BAR_H, ymax = y + BAR_H),
              colour = "grey25", linewidth = 0.3) +
    geom_rect(data = band, inherit.aes = FALSE,
              aes(xmin = xmin / 1e6, xmax = xmax / 1e6,
                  ymin = ymin, ymax = ymax),
              fill = COL_REGION, alpha = 0.22) +
    geom_segment(data = band, inherit.aes = FALSE,
                 aes(x = xmin / 1e6, xend = xmin / 1e6,
                     y = ymin, yend = ymax),
                 linetype = "dotted", linewidth = 0.4, colour = COL_REGION) +
    geom_segment(data = band, inherit.aes = FALSE,
                 aes(x = xmax / 1e6, xend = xmax / 1e6,
                     y = ymin, yend = ymax),
                 linetype = "dotted", linewidth = 0.4, colour = COL_REGION) +
    geom_richtext(data = RESOLVED %>%
                    mutate(mid = (xmin + xmax) / 2e6,
                           lab = sprintf("%.3f&ndash;%.3f Mb",
                                         xmin / 1e6, xmax / 1e6)),
                  inherit.aes = FALSE,
                  aes(x = mid, y = length(LEVELS) + 0.62, label = lab),
                  colour = COL_REGION, size = 3, hjust = 0.5, vjust = 0.5,
                  ## nothing is drawn behind the label any more, so it needs
                  ## no white plate
                  fill = NA, label.color = NA,
                  label.padding = grid::unit(rep(0, 4), "pt")) +
    scale_fill_manual(values = c(JU1793 = COL_JU1793, JU2466 = COL_JU2466),
                      name = NULL) +
    y_axis(labels = FALSE) +
    coord_cartesian(xlim = WIN / 1e6) +
    labs(x = "Chromosome III (Mb)", y = NULL, title = panel_title(letter)) +
    theme_pub() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.position = "top",
          legend.margin = margin(b = -6))
}

## ===========================================================================
## C -- hatching under pos-1 RNAi
## ===========================================================================
panel_C <- function(verbose = TRUE, letter = "C") {
  z <- qnorm(0.975)
  phen <- fread(PHEN) %>% as_tibble() %>%
    filter(Strain %in% LEVELS) %>%
    transmute(strain = factor(Strain, levels = rev(LEVELS)), row = ROW[Strain],
              cond = ifelse(condition == "ht115", "HT115", "*pos-1*"),
              n = `plated embryo`, hatched = `plated embryo` - unhatched) %>%
    mutate(p = hatched / n,
           ## Wilson interval, which behaves near 0 and 1 where Wald does not
           lo = pmax(0, (p + z^2 / (2 * n) -
                           z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
                       (1 + z^2 / n)),
           hi = pmin(1, (p + z^2 / (2 * n) +
                           z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
                       (1 + z^2 / n)),
           cond = factor(cond, levels = c("HT115", "*pos-1*"))) %>%
    ## main display shows pos-1 only; the HT115 control arm is
    ## SUPP_FIG_XX_nil_hatching_full.R, which shows it is ~97-100% throughout
    filter(cond == "*pos-1*") %>%
    mutate(yc = row)

  if (verbose) {
    cat("\n== hatching, fraction (Wilson 95% CI) ==\n")
    print(as.data.frame(phen %>%
      transmute(strain, condition = gsub("\\*", "", cond), embryos = n,
                hatched = round(p, 3),
                CI = sprintf("%.3f-%.3f", lo, hi)) %>%
      arrange(condition, desc(strain))), row.names = FALSE)
  }

  ggplot(phen) +
    geom_rect(aes(xmin = 0, xmax = p, ymin = yc - 0.3, ymax = yc + 0.3),
              fill = "grey70", colour = "grey25", linewidth = 0.3) +
    geom_segment(aes(x = lo, xend = hi, y = yc, yend = yc),
                 linewidth = 0.4, colour = "grey15") +
    geom_segment(aes(x = lo, xend = lo, y = yc - 0.11, yend = yc + 0.11),
                 linewidth = 0.4, colour = "grey15") +
    geom_segment(aes(x = hi, xend = hi, y = yc - 0.11, yend = yc + 0.11),
                 linewidth = 0.4, colour = "grey15") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1.04), expand = expansion(mult = c(0, 0))) +
    y_axis(labels = FALSE) +
    labs(x = "Embryos hatched", y = NULL, title = panel_title(letter)) +
    theme_pub() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.position = "none")
}
