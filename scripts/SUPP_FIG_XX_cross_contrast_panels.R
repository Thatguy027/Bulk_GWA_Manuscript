## Supplement -- the two cross contrasts, panel per cross --------------------
## Pooled GWAS (mirrored) above, then one panel per cross with the two contrasts
## overlaid: HT115 vs pos-1 in purple, mig-6 vs pos-1 in blue. Purple high and
## blue flat marks a general RNAi-response locus, blue high a knockdown-specific
## one. This is the expanded view behind Figure 2.
##
##   Rscript scripts/pooled_cross_intersection_prep.R   # once, builds the cache
##   Rscript scripts/SUPP_FIG_XX_cross_contrast_panels.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggtext)
})

b   <- readRDS("supplemental_data/mapping/pooled_cross_bundle_thinned.rds")
## a kept supplement (see README.md), so it writes to plots/ -- it used to
## write into plots/pooled_cross_intersection/, which is not tracked, so
## re-running it left the curated copy stale
OUT <- "plots"

## Bonferroni assumes independent markers, which this panel is not. The eigen
## threshold uses the effective number of independent tests from the eigenvalues
## of the marker correlation matrix (Li & Ji 2005); see
## scripts/eigen_independent_tests.R.
source("scripts/gwas_thresholds.R")
TH <- gwas_thresholds("pooled_RNAi_expt")
GWAS_BF    <- TH$bonferroni
GWAS_EIGEN <- TH$eigen
COL_THR    <- "grey50"
COL_EIG    <- "#1A7F5A"

## -- colours -----------------------------------------------------------------
COL_MIG  <- "#2166AC"   # mig-6 trait, and the mig-6 vs pos-1 contrast
COL_POS  <- "#7C6A9C"   # pos-1 trait, and the HT115 vs pos-1 contrast
COL_PEAK <- "#C4302B"   # the significant pooled GWAS peaks and their window

CROSS_LAB <- c(N2xXZ1516 = "N2 × XZ1516", JU1793xJU2466 = "JU1793 × JU2466")
CROSS_COL <- c(N2xXZ1516 = "#0B7A75", JU1793xJU2466 = "#D57A00")

## the two contrasts, given one meaning each regardless of the order the two
## exports happen to write them in
ROLE_OF  <- c(`ht115 vs pos1` = "response",
              `pos1 vs mig6`  = "specific",
              `mig6 vs pos1`  = "specific")
ROLE_LAB <- c(response = "HT115 ; *pos-1*",
              specific = "*mig-6* ; *pos-1*")
ROLE_COL <- c(response = COL_POS, specific = COL_MIG)

KEYS <- c("N2xXZ1516 | ht115 vs pos1", "N2xXZ1516 | pos1 vs mig6",
          "JU1793xJU2466 | ht115 vs pos1", "JU1793xJU2466 | mig6 vs pos1")

CHROMS  <- c("I", "II", "III", "IV", "V", "X")
ALL_LEN <- c(I = 15072434, II = 15279421, III = 13783801,
             IV = 17493829, V = 20924180, X = 17718942)

fct_chr <- function(x) factor(as.character(x), levels = CHROMS)

## "A  Pooled GWAS ..." -- the letter is set in the title rather than as a
## patchwork tag so it is guaranteed flush top-left of its own panel and cannot
## overlap the title text
panel_title <- function(letter, txt, colour) {
  paste0("<span style='font-size:14pt;color:#111111'>**", letter, "**</span>",
         "<span style='color:", colour, "'>  ", txt, "</span>")
}
wrap_md <- function(x, n = 128) paste(strwrap(x, width = n), collapse = "<br>")

## facet_grid(space = "free_x") sizes each panel by its own data range, so every
## panel is pinned to the full chromosome to keep the rows aligned
chrom_span <- bind_rows(
  tibble(chrom = fct_chr(names(ALL_LEN)), pos.mb = unname(ALL_LEN) / 1e6),
  tibble(chrom = fct_chr(names(ALL_LEN)), pos.mb = 0))

## ===========================================================================
## the only annotation: 1 Mb either side of each significant pooled GWAS peak
## ===========================================================================
FLANK <- 1e6
peaks <- b$gwas_sig %>%
  transmute(gene, chrom = fct_chr(chr), ps, pos.mb = ps / 1e6, neglog10p,
            y = ifelse(gene == "mig-6", neglog10p, -neglog10p),
            xmin = pmax(ps - FLANK, 0) / 1e6,
            xmax = pmin(ps + FLANK, ALL_LEN[as.character(chr)]) / 1e6)

peak_band <- function() list(
  geom_rect(data = peaks, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = COL_PEAK, alpha = 0.10),
  geom_vline(data = peaks, inherit.aes = FALSE, aes(xintercept = xmin),
             linetype = "dotted", linewidth = 0.25, colour = COL_PEAK),
  geom_vline(data = peaks, inherit.aes = FALSE, aes(xintercept = xmax),
             linetype = "dotted", linewidth = 0.25, colour = COL_PEAK))

## ===========================================================================
## theme / facets
## ===========================================================================
theme_pub <- function(base_size = 11.5) {
  theme_classic(base_size = base_size) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = base_size),
          panel.spacing.x = grid::unit(5, "pt"),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          plot.title.position = "plot",
          legend.key.size = grid::unit(9, "pt"))
}

## no `switch`, so the chromosome strip sits at the TOP of whichever panel shows it
genome_facets <- function() list(
  geom_blank(data = chrom_span, aes(x = pos.mb, y = 0), inherit.aes = FALSE),
  facet_grid(. ~ chrom, scales = "free_x", space = "free_x"),
  scale_x_continuous(breaks = seq(0, 25, 5), expand = expansion(mult = 0.02),
                     guide = guide_axis(check.overlap = TRUE)))

## ticks stay on every panel so the reader can register position in A and B;
## only the bottom panel carries the tick labels
axis_ctl <- function(strip = FALSE, x = FALSE) theme(
  strip.text   = if (strip) NULL else element_blank(),
  axis.ticks.x = element_line(linewidth = 0.3),
  axis.text.x  = if (x) element_text(size = 9) else element_blank())

## thin, so a 500k-marker trace is not written at 150 vertices per pixel
thin_extreme <- function(d, y, bp = 5e3, by = NULL) {
  g <- c(by, "chrom", ".bin")
  d %>% mutate(.bin = floor(pos / bp)) %>%
    group_by(across(all_of(g))) %>%
    slice_max(abs(.data[[y]]), n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(-.bin) %>%
    arrange(across(all_of(c(by, "chrom", "pos"))))
}

## ===========================================================================
## data
## ===========================================================================
scans <- imap_dfr(b$scans[KEYS], function(s, k) {
  m <- b$scan_meta %>% filter(key == k)
  s %>% filter(is.finite(LOD), chrom %in% CHROMS) %>%
    transmute(key = k, cross = m$cross, role = ROLE_OF[[m$label]],
              chrom = fct_chr(chrom), pos = physical.position,
              pos.mb = physical.position / 1e6, LOD)
}) %>% thin_extreme("LOD", by = "key") %>%
  mutate(role = factor(role, levels = names(ROLE_LAB)))

gw <- b$gwas %>% filter(chr %in% CHROMS) %>%
  transmute(gene, chrom = fct_chr(chr), pos = ps, pos.mb = ps / 1e6, neglog10p) %>%
  thin_extreme("neglog10p", bp = 1e4, by = "gene") %>%
  mutate(y = ifelse(gene == "mig-6", neglog10p, -neglog10p))

## which trait is which, stated against the threshold line on chromosome I
## rather than in the axis title
trait_lab <- tibble(
  chrom  = fct_chr(c("I", "I")),
  pos.mb = 0.3,
  y      = c(GWAS_BF + 1.1, -(GWAS_BF + 1.1)),
  label  = c("*mig-6* RNAi", "*pos-1* RNAi"),
  col    = c(COL_MIG, COL_POS))

## ===========================================================================
## row 1: pooled GWAS, mirrored, chromosome names on top
## ===========================================================================
p_gwas <- ggplot(gw, aes(pos.mb, y)) +
  peak_band() +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey30") +
  geom_hline(yintercept = c(-1, 1) * GWAS_BF, linetype = "dashed",
             linewidth = 0.3, colour = COL_THR) +
  geom_hline(yintercept = c(-1, 1) * GWAS_EIGEN, linetype = "dotted",
             linewidth = 0.42, colour = COL_EIG) +
  geom_point(data = gw %>% filter(gene == "mig-6"), size = 0.4, alpha = 0.5,
             colour = COL_MIG) +
  geom_point(data = gw %>% filter(gene == "pos-1"), size = 0.4, alpha = 0.5,
             colour = COL_POS) +
  geom_point(data = peaks, aes(y = y), colour = COL_PEAK, size = 2) +
  geom_text(data = peaks, aes(y = y, label = sprintf("%.2f", neglog10p)),
            vjust = ifelse(peaks$gene == "mig-6", -0.9, 1.8), size = 3.1,
            colour = COL_PEAK, fontface = "bold") +
  geom_richtext(data = trait_lab, aes(y = y, label = label),
                colour = trait_lab$col, size = 3.3, hjust = 0, vjust = 0.5,
                fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  genome_facets() +
  scale_y_continuous(labels = abs, expand = expansion(mult = 0.2)) +
  labs(x = NULL, y = "−log<sub>10</sub>*p*",
       title = panel_title("A", "Pooled GWAS, 93 wild isotypes", "grey20")) +
  theme_pub(11.5) +
  theme(plot.title = element_markdown(size = 11.5),
        axis.title.y = element_markdown(size = 10),
        legend.position = "none") +
  axis_ctl(strip = TRUE, x = FALSE)

## ===========================================================================
## rows 2-3: one panel per cross, the two contrasts overlaid
## ===========================================================================
cross_panel <- function(cr, letter, xax = FALSE, legend = FALSE) {
  ggplot(scans %>% filter(cross == cr), aes(pos.mb, LOD, colour = role)) +
    peak_band() +
    geom_hline(yintercept = b$CROSS_THR, linetype = "dashed",
               linewidth = 0.3, colour = "grey50") +
    geom_line(aes(group = key), linewidth = 0.55) +
    scale_colour_manual(values = ROLE_COL, breaks = names(ROLE_LAB),
                        labels = ROLE_LAB, name = NULL) +
    genome_facets() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = if (xax) "Genomic Position (Mb)" else NULL, y = "LOD",
         title = panel_title(letter, CROSS_LAB[[cr]], CROSS_COL[[cr]])) +
    theme_pub(11.5) +
    theme(plot.title = element_markdown(size = 11.5),
          axis.title.y = element_text(size = 10),
          legend.position = if (legend) "top" else "none",
          legend.justification = "left",
          legend.direction = "horizontal",
          legend.text = element_markdown(size = 9.5),
          legend.key.width = grid::unit(17, "pt"),
          legend.margin = margin(t = -2, b = -6)) +
    guides(colour = guide_legend(nrow = 1)) +
    axis_ctl(strip = FALSE, x = xax)
}

## ===========================================================================
## assemble
## ===========================================================================
fig <- p_gwas / cross_panel("N2xXZ1516", "B", legend = TRUE) /
  cross_panel("JU1793xJU2466", "C", xax = TRUE) +
  plot_layout(heights = c(1, 1.3, 1.3)) +
  plot_annotation(theme = theme(plot.margin = margin(6, 8, 4, 8)))

ggsave(file.path(OUT, "SUPP_FIG_XX_cross_contrast_panels.pdf"), fig,
       width = 13, height = 8.2, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_cross_contrast_panels.png"), fig,
       width = 13, height = 8.2, dpi = 300, bg = "white")
cat("wrote SUPP_FIG_XX_cross_contrast_panels.{pdf,png}\n")

## ---------------------------------------------------------------------------
## the same argument as numbers
## ---------------------------------------------------------------------------
LOCI <- tribble(
  ~locus,     ~chrom, ~pos,
  "I:1.9",    "I",     1876106,
  "III:12.4", "III",  12353680,
  "III:13.7", "III",  13679445,
  "IV:9.0",   "IV",    8975323,
  "V:13.8",   "V",    13805216,
  "V:14.6",   "V",    14647434,
  "X:6.9",    "X",     6904721,
  "X:9.7",    "X",     9705058)

tab <- pmap_dfr(LOCI, function(locus, chrom, pos) {
  scans %>% filter(as.character(chrom) == !!chrom) %>%
    group_by(cross, role) %>%
    slice_min(abs(pos - !!pos), n = 1, with_ties = FALSE) %>%
    ungroup() %>% transmute(locus, cross, role, LOD = round(LOD, 1))
}) %>% pivot_wider(names_from = c(cross, role), values_from = LOD)

write_tsv(tab, file.path(OUT, "TABLE_general_vs_specific.tsv"))
cat("\n== LOD at each locus, by cross and contrast ==\n")
print(as.data.frame(tab), row.names = FALSE)
