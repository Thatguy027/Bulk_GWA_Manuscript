## Figure 2 -------------------------------------------------------------------
## Pooled GWAS Manhattan with the cross QTL compressed onto interval tracks:
## mig-6-specific above the axis, the pos-1 response below, coloured by cross,
## opacity scaled to peak LOD. Only QTL with peak LOD > 100 are drawn; the
## complete set is SUPP_FIG_XX_cross_qtl_all.R.
##
##   Rscript scripts/pooled_cross_intersection_prep.R   # once, builds the cache
##   Rscript scripts/Figure2.R
##
## The GWAS traits are the vst-transformed pooled phenotypes (vst_ctrl_<gene>_T2),
## so this figure is already on that scale; nothing here re-derives a phenotype.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggtext)
})

b   <- readRDS("supplemental_data/mapping/pooled_cross_bundle_thinned.rds")
## Figure 2 is a manuscript figure, so it writes to plots/. It used to write
## into plots/pooled_cross_intersection/, which is untracked, so the copy in
## plots/ was a stale hand-placed duplicate that re-running never refreshed.
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

DROP_FRAC <- 0.05

COL_MIG  <- "#2166AC"   # mig-6 trait
COL_POS  <- "#7C6A9C"   # pos-1 trait
COL_PEAK <- "#C4302B"   # significant pooled GWAS peaks
CROSS_LAB <- c(N2xXZ1516 = "N2 × XZ1516", JU1793xJU2466 = "JU1793 × JU2466")
CROSS_COL <- c(N2xXZ1516 = "#0B7A75", JU1793xJU2466 = "#D57A00")

CHROMS  <- c("I", "II", "III", "IV", "V", "X")
ALL_LEN <- c(I = 15072434, II = 15279421, III = 13783801,
             IV = 17493829, V = 20924180, X = 17718942)
fct_chr <- function(x) factor(as.character(x), levels = CHROMS)

## the four scans the classification needs, per cross
SCAN <- tribble(
  ~cross,          ~cond,   ~key,
  "N2xXZ1516",     "mig-6", "N2xXZ1516 | ht115 vs mig6",
  "N2xXZ1516",     "pos-1", "N2xXZ1516 | ht115 vs pos1",
  "N2xXZ1516",     "contr", "N2xXZ1516 | pos1 vs mig6",
  "JU1793xJU2466", "mig-6", "JU1793xJU2466 | ht115 vs mig6",
  "JU1793xJU2466", "pos-1", "JU1793xJU2466 | ht115 vs pos1",
  "JU1793xJU2466", "contr", "JU1793xJU2466 | mig6 vs pos1")

## ===========================================================================
## intervals at a 5% LOD drop
## ===========================================================================
## Faithful to peak_intervals() in the export pipeline: take the chromosome
## peak, walk outward, and keep the contiguous run of markers staying within
## `frac` of that peak.
MIN_GAP       <- 1e6    # buffer masked either side of a called interval
PEAK_MIN_FRAC <- 0.30   # a second peak must reach this share of the chromosome max
MAX_PEAKS     <- 1      # see note above: one interval per chromosome per scan
MAX_WIDTH     <- Inf

peak_intervals_frac <- function(d, frac = DROP_FRAC, threshold = 0) {
  d <- d[is.finite(d$LOD), ]
  if (!nrow(d)) return(NULL)
  bind_rows(lapply(split(d, as.character(d$chrom)), function(x) {
    x <- x[order(x$pos), ]
    chrom_max <- max(x$LOD)
    avail <- rep(TRUE, nrow(x))
    out <- list()
    for (k in seq_len(MAX_PEAKS)) {
      if (!any(avail)) break
      i <- which(avail)[which.max(x$LOD[avail])]
      if (x$LOD[i] <= threshold || x$LOD[i] < PEAK_MIN_FRAC * chrom_max) break
      ## The run must stay inside the markers still available. Measuring it on
      ## the full trace let a secondary peak's interval bleed through a taller
      ## neighbour: the N2 x XZ1516 pos-1 chromosome III peak at 12.11 Mb (LOD
      ## 245) ran to 13.78 Mb because the 709-LOD primary peak at 13.31 sits
      ## above its 5% floor. Masking bounds it correctly.
      above <- avail & (x$LOD > (x$LOD[i] - x$LOD[i] * frac))
      run <- cumsum(c(TRUE, diff(above) != 0))
      sel <- which(run == run[i])
      lo <- min(sel); hi <- max(sel)
      out[[k]] <- tibble(chrom = as.character(x$chrom[1]),
                         peak.position = x$pos[i], peak.LOD = x$LOD[i],
                         lcon = x$pos[lo], rcon = x$pos[hi], peak.rank = k)
      avail[x$pos >= x$pos[lo] - MIN_GAP & x$pos <= x$pos[hi] + MIN_GAP] <- FALSE
    }
    bind_rows(out)
  })) %>%
    mutate(width.kb = (rcon - lcon) / 1e3,
           significant = peak.LOD > threshold & (rcon - lcon) <= MAX_WIDTH)
}

## full-resolution scans for the interval maths (thinning happens only to plot)
scan_of <- function(key) {
  b$scans[[key]] %>% filter(is.finite(LOD), chrom %in% CHROMS) %>%
    transmute(chrom = as.character(chrom), pos = physical.position, LOD)
}
scans <- set_names(map(SCAN$key, scan_of), paste(SCAN$cross, SCAN$cond))

lod_at <- function(nm, ch, ps) {
  d <- scans[[nm]]; d <- d[d$chrom == ch, ]
  if (!nrow(d)) return(NA_real_)
  d$LOD[which.min(abs(d$pos - ps))]
}

## ===========================================================================
## intervals: every significant QTL from each HT115 contrast
## ===========================================================================
## upper track from the contrast, lower from HT115 vs pos-1
TRACK_SRC <- tribble(
  ~track,   ~cond,
  "mig-6",  "contr",
  "pos-1",  "pos-1")

ivs <- pmap_dfr(SCAN %>% inner_join(TRACK_SRC, by = "cond") %>%
                  select(cross, cond, key, track),
                function(cross, cond, key, track) {
  peak_intervals_frac(scan_of(key), threshold = b$CROSS_THR) %>%
    mutate(cross = cross, cond = cond, track = track, .before = 1)
}) %>% filter(significant) %>%
  rowwise() %>%
  mutate(lod_mig6 = lod_at(paste(cross, "mig-6"), chrom, peak.position),
         lod_pos1 = lod_at(paste(cross, "pos-1"), chrom, peak.position)) %>%
  ungroup()

cat("cross genome-wide LOD threshold: ", round(b$CROSS_THR, 2), "\n",
    "intervals: ", 100 * DROP_FRAC,
    "% LOD drop from the peak marker of each HT115 contrast\n\n", sep = "")
cat("== every interval on the tracks ==\n")
print(as.data.frame(ivs %>%
  transmute(track = ifelse(track == "mig-6", "UP  mig-6 vs pos-1",
                           "DN  HT115 vs pos-1"),
            cross, chrom, peak_Mb = round(peak.position / 1e6, 2),
            interval = sprintf("%.2f-%.2f", lcon / 1e6, rcon / 1e6),
            width_kb = round(width.kb), peak_LOD = round(peak.LOD, 1),
            HT115_mig6 = round(lod_mig6, 1), HT115_pos1 = round(lod_pos1, 1),
            mig6_dominant = lod_mig6 > lod_pos1) %>%
  arrange(track, chrom, desc(peak_LOD))), row.names = FALSE)

spec <- ivs

## how close does each pooled GWAS peak come to a track interval of its own trait?
gp <- b$gwas_sig %>% transmute(gene, chrom = chr, ps)
cat("\n== each pooled GWAS peak vs the matching-trait cross intervals ==\n")
print(as.data.frame(pmap_dfr(gp, function(gene, chrom, ps) {
  s <- ivs %>% filter(chrom == !!chrom, track == gene)
  if (!nrow(s)) return(tibble(gwas = gene, at = paste0(chrom, ":", round(ps / 1e6, 2)),
                              cross = "-", cross_peak_Mb = NA_real_,
                              dist_Mb = NA_real_, contains_peak = NA))
  s %>% transmute(gwas = gene, at = paste0(chrom, ":", round(ps / 1e6, 2)), cross,
                  cross_peak_Mb = round(peak.position / 1e6, 2),
                  dist_Mb = round(abs(peak.position - ps) / 1e6, 2),
                  contains_peak = ps >= lcon & ps <= rcon)
})), row.names = FALSE)

write_tsv(ivs, file.path(OUT, "TABLE_cross_qtl_tracks_5pct.tsv"))

## ===========================================================================
## the Manhattan, with the interval tracks over and under it
## ===========================================================================
thin_extreme <- function(d, y, bp = 1e4, by = NULL) {
  g <- c(by, "chrom", ".bin")
  d %>% mutate(.bin = floor(pos / bp)) %>%
    group_by(across(all_of(g))) %>%
    slice_max(abs(.data[[y]]), n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(-.bin) %>% arrange(across(all_of(c(by, "chrom", "pos"))))
}

gw <- b$gwas %>% filter(chr %in% CHROMS) %>%
  transmute(gene, chrom = fct_chr(chr), pos = ps, pos.mb = ps / 1e6, neglog10p) %>%
  thin_extreme("neglog10p", by = "gene") %>%
  mutate(y = ifelse(gene == "mig-6", neglog10p, -neglog10p))

peaks <- b$gwas_sig %>%
  transmute(gene, chrom = fct_chr(chr), pos.mb = ps / 1e6, neglog10p,
            y = ifelse(gene == "mig-6", neglog10p, -neglog10p))

## 1 Mb either side of each bulk GWAS peak. Each window is drawn only on its own
## side of the axis -- the mig-6 peak upward, the pos-1 peak downward -- so a
## window is never read against the trait it does not belong to.
FLANK <- 1e6
bulk_band <- b$gwas_sig %>%
  transmute(gene, chrom = fct_chr(chr),
            xmin = pmax(ps - FLANK, 0) / 1e6,
            xmax = pmin(ps + FLANK, ALL_LEN[as.character(chr)]) / 1e6)

peak_band <- function(top) {
  d <- bulk_band %>%
    mutate(ylo = ifelse(gene == "mig-6", 0, -top),
           yhi = ifelse(gene == "mig-6", top, 0))
  list(
    geom_rect(data = d, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = ylo, ymax = yhi),
              fill = COL_PEAK, alpha = 0.09),
    geom_segment(data = d, inherit.aes = FALSE,
                 aes(x = xmin, xend = xmin, y = ylo, yend = yhi),
                 linetype = "dotted", linewidth = 0.25, colour = COL_PEAK),
    geom_segment(data = d, inherit.aes = FALSE,
                 aes(x = xmax, xend = xmax, y = ylo, yend = yhi),
                 linetype = "dotted", linewidth = 0.25, colour = COL_PEAK))
}

## One lane per cross in each track, so two crosses overlapping the same region
## stay legible. Lanes sit outside the data range, which the y expansion opens up.
LANE <- tibble(cross = names(CROSS_LAB), lane = c(0, 1))
BASE <- 8.7        # first lane, just clear of the thresholds
STEP <- 1.15       # lane spacing
BAR  <- 0.34       # bar half-height

build <- function(lod_min, name, alpha_breaks) {
  alpha_lim <- range(c(alpha_breaks, spec$peak.LOD[spec$peak.LOD > lod_min]))

bars <- spec %>%
  filter(peak.LOD > lod_min) %>%
  left_join(LANE, by = "cross") %>%
  mutate(chrom = fct_chr(chrom),
         sign = ifelse(track == "mig-6", 1, -1),
         yc = sign * (BASE + lane * STEP),
         ymin = yc - BAR, ymax = yc + BAR,
         xmin = lcon / 1e6, xmax = rcon / 1e6,
         cross = factor(cross, levels = names(CROSS_LAB)))

## a bar can be narrower than a pixel, so give every one a visible minimum
MINW <- 0.09
bars <- bars %>% mutate(
  xmid = (xmin + xmax) / 2,
  xmin = pmin(xmin, xmid - MINW / 2),
  xmax = pmax(xmax, xmid + MINW / 2))

ylim_hi <- BASE + STEP + BAR + 0.55
track_lab <- tibble(
  chrom = fct_chr(c("I", "I")),
  pos.mb = 0.25,
  y = c(ylim_hi - 0.25, -(ylim_hi - 0.25)),
  label = c("*mig-6*-specific", "*pos-1* response"))

## sit on the significance line itself, each label just outside it
## name the two thresholds at the far right, clear of the trait labels
## chromosome I, left-aligned. The trait labels already sit just above the
## Bonferroni line there, so this one goes just below it and the eigen label
## just above its own line.
thr_lab <- tibble(
  chrom  = fct_chr(rep("I", 2)),
  pos.mb = 0.3,
  vj     = c(1.35, -0.35),
  y      = c(GWAS_BF, GWAS_EIGEN),
  label  = c(sprintf("Bonferroni %.2f", GWAS_BF),
             sprintf("eigen %.2f (M<sub>eff</sub> = %d)", GWAS_EIGEN, round(TH$m_eff))),
  col    = c(COL_THR, COL_EIG))

trait_lab <- tibble(
  chrom  = fct_chr(c("I", "I")),
  pos.mb = 0.25,
  y      = c(GWAS_BF, -GWAS_BF),
  vj     = c(-0.35, 1.35),
  label  = c("*mig-6* RNAi", "*pos-1* RNAi"),
  col    = c(COL_MIG, COL_POS))

fig <- ggplot(gw, aes(pos.mb, y)) +
  ## chromosome-length blanks keep the panels aligned and full width
  geom_blank(data = bind_rows(
    tibble(chrom = fct_chr(names(ALL_LEN)), pos.mb = unname(ALL_LEN) / 1e6, y = 0),
    tibble(chrom = fct_chr(names(ALL_LEN)), pos.mb = 0, y = 0))) +
  ## the bulk GWAS QTL windows, behind everything
  peak_band(ylim_hi) +
  ## separators between the Manhattan and each track
  geom_hline(yintercept = c(-1, 1) * (BASE - BAR - 0.35),
             linewidth = 0.25, colour = "grey85") +
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
            vjust = ifelse(peaks$gene == "mig-6", 1.8, -0.9), hjust = 1.15,
            size = 3.1, colour = COL_PEAK, fontface = "bold") +
  ## the interval tracks
  geom_rect(data = bars, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = cross, alpha = peak.LOD),
            colour = NA) +
  geom_richtext(data = trait_lab, aes(y = y, label = label, vjust = vj),
                colour = trait_lab$col, size = 3.3, hjust = 0,
                fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  geom_richtext(data = thr_lab, aes(y = y, label = label, vjust = vj),
                colour = thr_lab$col, size = 2.6, hjust = 0,
                fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  geom_richtext(data = track_lab, aes(y = y, label = label),
                colour = "grey35", size = 3.1, hjust = 0, vjust = 0.5,
                fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  scale_fill_manual(values = CROSS_COL, labels = CROSS_LAB, name = NULL) +
  ## LOD spans a wide range, so opacity is on a square-root scale
  scale_alpha_continuous(range = c(0.30, 1), trans = "sqrt",
                         breaks = alpha_breaks, limits = alpha_lim,
                         name = "QTL peak LOD  ") +
  facet_grid(. ~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, 25, 5), expand = expansion(mult = 0.02),
                     guide = guide_axis(check.overlap = TRUE)) +
  scale_y_continuous(breaks = seq(-6, 6, 3), labels = abs,
                     limits = c(-ylim_hi, ylim_hi), expand = expansion(mult = 0)) +
  labs(x = "Genomic Position (Mb)", y = "−log<sub>10</sub>*p*") +
  theme_classic(base_size = 11.5) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),
        panel.spacing.x = grid::unit(5, "pt"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.text.x = element_text(size = 9),
        axis.title.y = element_markdown(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 9.5),
        legend.key.height = grid::unit(8, "pt"),
        legend.title = element_text(size = 9.5),
        legend.box = "horizontal",
        legend.margin = margin(t = -4),
        plot.margin = margin(8, 10, 4, 8)) +
  guides(fill  = guide_legend(order = 1, override.aes = list(alpha = 1)),
         alpha = guide_legend(order = 2, override.aes = list(fill = "grey25")))

ggsave(file.path(OUT, paste0(name, ".pdf")), fig,
       width = 13, height = 5.4, device = cairo_pdf)
ggsave(file.path(OUT, paste0(name, ".png")), fig,
       width = 13, height = 5.4, dpi = 300, bg = "white")
cat("wrote ", name, ".{pdf,png}  (", nrow(bars), " intervals, LOD > ", lod_min,
    ")\n", sep = "")
invisible(bars)
}


## ---------------------------------------------------------------------------
LOD_MIN <- 100
kept <- build(LOD_MIN, "Figure2", c(100, 300, 600, 900))

cat("\n== intervals on Figure 2 (peak LOD > ", LOD_MIN, ") ==\n", sep = "")
print(as.data.frame(kept %>%
  transmute(track = ifelse(track == "mig-6", "UP  mig-6 vs pos-1",
                           "DN  HT115 vs pos-1"), cross, chrom,
            peak_Mb = round(peak.position / 1e6, 2),
            interval = sprintf("%.2f-%.2f", lcon / 1e6, rcon / 1e6),
            width_kb = round(width.kb), peak_LOD = round(peak.LOD, 1)) %>%
  arrange(track, chrom, desc(peak_LOD))), row.names = FALSE)
