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

## --- chromosome III parental allele frequency -----------------------------
##
## Replaces the LOD trace that used to be Figure 3B. A LOD trace answers "is
## there a QTL here", which Figure 2 already answers; it says nothing about
## WHICH parent's allele the pos-1 selection favoured, which is the point being
## made -- that the resistant parent's haplotype sweeps toward the right end of
## chromosome III, where sid-2 sits at 13.68 Mb.
##
## PARENT ASSIGNMENT IS FROM THE EXPORT, NOT INFERRED. The p1/p2 columns are,
## per data/cross_experiments/JU1793-JU2466_export/README.md, "the counts
## assigned to the JU1793 and JU2466 haplotypes". So JU1793 frequency is
## p1/(p1+p2). Getting this backwards would invert the whole panel and still
## look plausible, so the direction is asserted at the bottom of this function
## against the known result: the pos-1 pool must end up MORE JU1793 than the
## HT115 control at the right end.
##
## FREQUENCIES ARE COUNT-WEIGHTED, not averages of per-marker frequencies.
## Counts are summed within each physical bin and the frequency taken from the
## sums, so a marker with 200 reads counts for more than one with 4. Averaging
## per-marker frequencies would let the shallowest markers pull the trace
## around. The rolling mean is then applied to the binned frequencies for
## display only.
FREQ_BIN_KB  <- 50    # physical bin for the count-weighted frequency
FREQ_ROLL_N  <- 5     # centred rolling mean, in bins: 5 x 50 kb = 250 kb

load_parent_freq <- function(chrom_keep = "III", from_mb = 8) {
  nm <- names(fread(SCAN, nrows = 0))
  p1_ht <- grep("^p1_.*HT115g$", nm, value = TRUE)
  p2_ht <- grep("^p2_.*HT115g$", nm, value = TRUE)
  p1_ps <- grep("^p1_.*POS1g$",  nm, value = TRUE)
  p2_ps <- grep("^p2_.*POS1g$",  nm, value = TRUE)
  if (length(c(p1_ht, p2_ht, p1_ps, p2_ps)) != 4)
    stop("expected one p1/p2 column per sample in ", SCAN, ", found: ",
         paste(c(p1_ht, p2_ht, p1_ps, p2_ps), collapse = ", "), call. = FALSE)

  d <- fread(SCAN, select = c("chrom", "physical.position",
                              p1_ht, p2_ht, p1_ps, p2_ps)) %>% as_tibble()
  names(d) <- c("chrom", "pos", "p1_ht", "p2_ht", "p1_ps", "p2_ps")
  d <- d %>% filter(chrom == chrom_keep, pos >= from_mb * 1e6)
  msg("    ", nrow(d), " markers on ", chrom_keep, " from ", from_mb, " Mb")

  roll <- function(x, n) {
    if (n <= 1 || length(x) < n) return(x)
    stats::filter(x, rep(1 / n, n), sides = 2) %>% as.numeric() %>%
      ## the filter leaves NA at both ends; hold the nearest fitted value so
      ## the fill reaches the panel edges instead of stopping short
      (function(v) { k <- which(!is.na(v))
                     v[seq_len(min(k) - 1)] <- v[min(k)]
                     v[seq(max(k) + 1, length(v), length.out =
                             max(0, length(v) - max(k)))] <- v[max(k)]
                     v })()
  }

  b <- d %>%
    mutate(bin = floor(pos / (FREQ_BIN_KB * 1e3))) %>%
    group_by(bin) %>%
    summarise(pos.mb = mean(pos) / 1e6,
              f_ps = sum(p1_ps) / sum(p1_ps + p2_ps),
              f_ht = sum(p1_ht) / sum(p1_ht + p2_ht),
              n = n(), .groups = "drop") %>%
    arrange(pos.mb) %>%
    mutate(f_ps_s = roll(f_ps, FREQ_ROLL_N), f_ht_s = roll(f_ht, FREQ_ROLL_N))

  ## the assertion that catches a swapped parent assignment
  tail_n <- max(3, round(nrow(b) * 0.1))
  end_ps <- mean(tail(b$f_ps_s, tail_n)); end_ht <- mean(tail(b$f_ht_s, tail_n))
  msg("    JU1793 frequency over the last ", tail_n, " bins: pos-1 ",
      sprintf("%.3f", end_ps), " vs HT115 ", sprintf("%.3f", end_ht))
  if (!(end_ps > end_ht))
    stop("the pos-1 pool is not enriched for JU1793 at the right end of ",
         chrom_keep, " (", sprintf("%.3f vs %.3f", end_ps, end_ht), ").\n",
         "  JU1793 is the RESISTANT parent, so pos-1 selection must raise its\n",
         "  frequency relative to the HT115 control. If this fails, p1/p2 have\n",
         "  most likely been swapped -- see the README quoted above.",
         call. = FALSE)
  b
}

panel_parent_freq_chr3 <- function(b, letter = "B") {
  xr <- range(b$pos.mb)
  ggplot(b, aes(pos.mb)) +
    ## JU2466 above the trace, JU1793 below it: the two fills sum to 1, so the
    ## panel reads as "which parent occupies this interval of the pool"
    geom_ribbon(aes(ymin = f_ps_s, ymax = 1), fill = COL_JU2466, alpha = 0.85) +
    geom_ribbon(aes(ymin = 0, ymax = f_ps_s), fill = COL_JU1793, alpha = 0.85) +
    geom_hline(yintercept = 0.5, linewidth = 0.3, linetype = "dashed",
               colour = "grey30") +
    ## the HT115 control, as a line only. It is the reason this panel shows
    ## selection rather than a segregation artefact: the control runs the other
    ## way over the same interval.
    geom_line(aes(y = f_ht_s), linewidth = 0.5, colour = "grey20") +
    annotate("richtext", x = xr[1], y = 0.5, hjust = -0.03, vjust = -0.45,
             label = "HT115 control", size = 2.6, colour = "grey20",
             fill = NA, label.color = NA,
             label.padding = grid::unit(rep(0, 4), "pt")) +
    annotate("richtext", x = xr[2], y = 0.02, hjust = 1.02, vjust = 0,
             label = "JU1793", size = 2.9, colour = "white",
             fill = NA, label.color = NA, fontface = "bold",
             label.padding = grid::unit(rep(0, 4), "pt")) +
    annotate("richtext", x = xr[2], y = 0.98, hjust = 1.02, vjust = 1,
             label = "JU2466", size = 2.9, colour = "white",
             fill = NA, label.color = NA, fontface = "bold",
             label.padding = grid::unit(rep(0, 4), "pt")) +
    scale_x_continuous("Chromosome III (Mb)", expand = expansion(0)) +
    scale_y_continuous("Parental allele frequency",
                       labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1), expand = expansion(0)) +
    labs(title = panel_title(letter)) +
    theme_pub() +
    theme(panel.grid = element_blank())
}

## --- C: NIL genotype and hatching, in ONE panel ---------------------------
##
## Replaces the old panels C and D, which drew the same five strains as rows
## twice, in two panels with two x axes, so a reader had to carry a row
## position across a panel boundary to pair a genotype with its phenotype.
## Here each row is one strain and the genotype sits immediately left of its
## own hatching bar.
##
## TWO X SCALES IN ONE PANEL, so both are drawn on an abstract coordinate: the
## genotype occupies [0, GEN_W], the hatching [GEN_W + GAP, GEN_W + GAP +
## HAT_W], and each gets its own tick labels annotated beneath the rows. A
## single scale_x_continuous cannot carry two units, and patchwork would put
## them back into two panels, which is what this change is undoing.
##
## THE INTROGRESSIONS RUN TO THE END OF THE CHROMOSOME, and the panel has to
## say so, because that is what distinguishes the strains. Five of the six NILs
## carry a JU2466 segment that starts at a breakpoint and continues to the
## chromosome III terminus at 13,783,801; wSZ191 is the exception, an internal
## 13.658-13.695 Mb segment that stops short. Compressing the axis to the
## informative span alone would crop exactly that contrast, so the window keeps
## its right edge at the terminus and marks it, and the compression is in the
## left flank and the row height instead.
## The genotype track is deliberately the SMALLER half. It carries two
## breakpoints; the hatching bars carry a five-level ordinal series with
## intervals, so they earn the width. An equal split gave the genotypes half the
## panel to show 37 kb of real information.
GEN_W <- 0.58      # genotype track width, abstract x units
GAP   <- 0.13
HAT_W <- 0.92      # hatching axis width
## Row half-height. The bars are THICK vertically and short horizontally: the
## horizontal extent is set by GEN_W/HAT_W above, which is where "smaller"
## belonged, and thinning them vertically as well just made them hard to read.
BAR_H2 <- 0.30

## The window trims dead flank on the LEFT only. Its right edge must stay at
## the chromosome terminus: wSZ191 stopping short of it, where the others run
## to it, is the fine-mapping contrast, and cropping there would remove the
## point of the panel. 13.635 Mb leaves ~23 kb of proximal flank, enough to
## show the NILs are JU1793 proximal to their breakpoints.
GWIN <- c(13.635e6, 13783801)

panel_nil_geno_hatch <- function(verbose = TRUE, letter = "C",
                                 labels = TRUE) {
  gx <- function(pos) (pos - GWIN[1]) / diff(GWIN) * GEN_W
  hx <- function(p)   GEN_W + GAP + p * HAT_W

  bed <- fread(NILB, header = FALSE,
               col.names = c("chr", "start", "end", "strain", "geno")) %>%
    as_tibble() %>% filter(strain %in% LEVELS) %>%
    mutate(y = ROW[strain],
           other = ifelse(geno == "JU1793", "JU2466", "JU1793"),
           seg_start = pmax(start, GWIN[1]),
           is_parent = strain %in% c("JU1793", "JU2466"))
  ## every NIL segment in the file must actually reach the terminus, except
  ## wSZ191; if that changes, the sentence in the caption stops being true
  reach <- bed %>% filter(!is_parent) %>%
    transmute(strain, to_end = end >= ALL_LEN[["III"]])
  if (verbose)
    msg("    NIL segments reaching the chrIII terminus: ",
        paste(reach$strain[reach$to_end], collapse = ", "),
        " | stopping short: ",
        paste(reach$strain[!reach$to_end], collapse = ", "))

  bg  <- bed %>% transmute(y, xmin = gx(GWIN[1]), xmax = gx(GWIN[2]),
                           geno = ifelse(is_parent, geno, other))
  seg <- bed %>% filter(!is_parent) %>%
    transmute(y, xmin = gx(seg_start), xmax = gx(end), geno)

  z <- qnorm(0.975)
  phen <- fread(PHEN) %>% as_tibble() %>%
    filter(Strain %in% LEVELS, condition != "ht115") %>%
    transmute(strain = Strain, y = ROW[Strain],
              n = `plated embryo`, hatched = `plated embryo` - unhatched) %>%
    mutate(p = hatched / n,
           lo = pmax(0, (p + z^2 / (2 * n) -
                           z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
                       (1 + z^2 / n)),
           hi = pmin(1, (p + z^2 / (2 * n) +
                           z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
                       (1 + z^2 / n)))
  if (verbose) {
    cat("\n== hatching under pos-1 RNAi (Wilson 95% CI) ==\n")
    print(as.data.frame(phen %>%
      transmute(strain, embryos = n, hatched = round(p, 3),
                CI = sprintf("%.3f-%.3f", lo, hi)) %>%
      arrange(desc(hatched))), row.names = FALSE)
  }

  band <- RESOLVED %>% mutate(ymin = min(ROW) - BAR_H2,
                              ymax = max(ROW) + BAR_H2)
  y0   <- min(ROW) - 0.72          # where the two tick rows sit
  gen_ticks <- tibble(pos = c(13.65, 13.70, 13.75) * 1e6) %>%
    mutate(x = gx(pos), lab = sprintf("%.2f", pos / 1e6))
  hat_ticks <- tibble(p = c(0, 0.5, 1)) %>%
    mutate(x = hx(p), lab = scales::percent(p, accuracy = 1))

  ggplot() +
    ## genotype: background haplotype, then the introgression over it
    geom_rect(data = bg, aes(xmin = xmin, xmax = xmax, ymin = y - BAR_H2,
                             ymax = y + BAR_H2, fill = geno),
              colour = "grey25", linewidth = 0.3) +
    geom_rect(data = seg, aes(xmin = xmin, xmax = xmax, ymin = y - BAR_H2,
                              ymax = y + BAR_H2, fill = geno),
              colour = "grey25", linewidth = 0.3) +
    ## the resolved interval, tinting the bars rather than sitting behind them
    geom_rect(data = band, inherit.aes = FALSE,
              aes(xmin = gx(xmin), xmax = gx(xmax), ymin = ymin, ymax = ymax),
              fill = COL_REGION, alpha = 0.22) +
    geom_segment(data = band, inherit.aes = FALSE,
                 aes(x = gx(xmin), xend = gx(xmin), y = ymin, yend = ymax),
                 linetype = "dotted", linewidth = 0.4, colour = COL_REGION) +
    geom_segment(data = band, inherit.aes = FALSE,
                 aes(x = gx(xmax), xend = gx(xmax), y = ymin, yend = ymax),
                 linetype = "dotted", linewidth = 0.4, colour = COL_REGION) +
    geom_richtext(data = band %>% mutate(mid = gx((xmin + xmax) / 2)),
                  inherit.aes = FALSE,
                  aes(x = mid, y = max(ROW) + 0.62,
                      label = sprintf("%.3f&ndash;%.3f Mb",
                                      RESOLVED$xmin / 1e6, RESOLVED$xmax / 1e6)),
                  colour = COL_REGION, size = 2.9, hjust = 0.5,
                  fill = NA, label.color = NA,
                  label.padding = grid::unit(rep(0, 4), "pt")) +
    ## the chromosome end, named rather than implied
    annotate("segment", x = gx(GWIN[2]), xend = gx(GWIN[2]),
             y = min(ROW) - BAR_H2, yend = max(ROW) + BAR_H2,
             linewidth = 0.5, colour = "grey20") +
    annotate("richtext", x = gx(GWIN[2]), y = min(ROW) - 0.42,
             label = "end of III", size = 2.5, colour = "grey30",
             hjust = 1, vjust = 1, fill = NA, label.color = NA,
             label.padding = grid::unit(rep(0, 4), "pt")) +
    ## hatching bars, on their own stretch of the same abstract axis
    geom_rect(data = phen, aes(xmin = hx(0), xmax = hx(p),
                               ymin = y - BAR_H2, ymax = y + BAR_H2),
              fill = "grey70", colour = "grey25", linewidth = 0.3) +
    geom_segment(data = phen, aes(x = hx(lo), xend = hx(hi), y = y, yend = y),
                 linewidth = 0.4, colour = "grey15") +
    geom_segment(data = phen, aes(x = hx(lo), xend = hx(lo),
                                  y = y - 0.11, yend = y + 0.11),
                 linewidth = 0.4, colour = "grey15") +
    geom_segment(data = phen, aes(x = hx(hi), xend = hx(hi),
                                  y = y - 0.11, yend = y + 0.11),
                 linewidth = 0.4, colour = "grey15") +
    ## two tick rows and two sub-axis titles, annotated because the panel
    ## carries two units on one coordinate
    annotate("segment", x = gx(GWIN[1]), xend = gx(GWIN[2]), y = y0, yend = y0,
             linewidth = 0.3, colour = "grey30") +
    annotate("segment", x = hx(0), xend = hx(1), y = y0, yend = y0,
             linewidth = 0.3, colour = "grey30") +
    geom_text(data = gen_ticks, aes(x = x, y = y0, label = lab),
              vjust = 1.6, size = 2.6, colour = "grey25") +
    geom_text(data = hat_ticks, aes(x = x, y = y0, label = lab),
              vjust = 1.6, size = 2.6, colour = "grey25") +
    annotate("richtext", x = gx(mean(GWIN)), y = y0 - 0.42,
             label = "Chromosome III (Mb)", size = 3.1, colour = "grey15",
             vjust = 1, fill = NA, label.color = NA,
             label.padding = grid::unit(rep(0, 4), "pt")) +
    annotate("richtext", x = hx(0.5), y = y0 - 0.42,
             label = "Embryos hatched", size = 3.1, colour = "grey15",
             vjust = 1, fill = NA, label.color = NA,
             label.padding = grid::unit(rep(0, 4), "pt")) +
    scale_fill_manual(values = c(JU1793 = COL_JU1793, JU2466 = COL_JU2466),
                      name = NULL) +
    scale_y_continuous(breaks = seq_along(LEVELS),
                       labels = if (labels) LEVELS else NULL,
                       limits = c(y0 - 1.05, max(ROW) + 1.05),
                       expand = expansion(mult = 0)) +
    coord_cartesian(xlim = c(0, GEN_W + GAP + HAT_W), clip = "off") +
    labs(x = NULL, y = NULL, title = panel_title(letter)) +
    theme_pub() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          axis.text.y = element_text(size = 8.6),
          panel.grid = element_blank(),
          legend.position = "top", legend.margin = margin(b = -6))
}
