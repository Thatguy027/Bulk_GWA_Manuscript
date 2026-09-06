## Which chromosome III association is real? ---------------------------------
##
##   Rscript scripts/chrIII_association_support.R
##     -> plots/diagnostics/chrIII_association_support.{pdf,png}
##     -> plots/diagnostics/TABLE_association_clusters.tsv
##
##   A  the chromosome III scan, with the two candidate signals marked
##   B  III:5,965,738 in close-up: 628 markers around it, all null
##   C  III:12.70-12.80 Mb in close-up: 15 markers above the eigen threshold
##   D  every eigen-significant cluster genome-wide, by width and peak
##
## THE QUESTION. The 2023 pos-1 scan's strongest chromosome III marker is a
## single site at 5,965,738 (-log10 p = 8.68), the only chromosome III marker to
## clear Bonferroni. The report has been treating that as "the chromosome III
## association peak" and noting that it lies 7.7 Mb from sid-2. This script asks
## whether that marker is an association at all, and what else on the
## chromosome would qualify if it is not.
##
## THE TEST IS LOCAL SUPPORT, AND IT NEEDS NO EXTRA DATA. A causal variant is
## detected through the markers in linkage disequilibrium with it, so a real
## association is a CLUSTER: the peak marker sits in a run of neighbours that
## are also elevated. A lone significant marker in a well-covered region is the
## signature of something that is not inherited with its neighbourhood --
## genotyping error, alignment artefact, a mismapped duplication -- because
## whatever it is tracks no haplotype.
##
## THE RECOMBINATION DOMAIN IS WHAT MAKES THIS DECISIVE. C. elegans chromosomes
## have low-recombination centres and high-recombination arms, so LD blocks are
## long in the centre and short on the arms. III:5.97 Mb is in the CENTRE, where
## support is easiest to come by and a real association should drag up markers
## across hundreds of kilobases. It drags up none: of the 628 markers within
## 100 kb, zero exceed the eigen threshold and their median -log10 p is 0.11.
## The right-arm signal at 12.70-12.80 Mb is on the ARM, where LD decays fast,
## and there a real association should look like a tight cluster -- which is
## exactly what 15 markers over 97 kb is.
##
## Domain boundaries are the Rockman & Kruglyak (2009) centre definitions,
## hard-coded below. The conclusion does not depend on them being exact:
## 5.97 Mb sits well inside the chromosome III centre on any version of these
## boundaries and 12.72 Mb well outside it.
##
## WHAT THIS DOES NOT SHOW. It does not put the association on sid-2. The 37 kb
## NIL interval contains 88 markers whose best -log10 p is 0.98, and the T96K
## variant itself is 0.62, ranking 18,662 of 64,423 on the chromosome. The
## right-arm cluster is 0.86-0.96 Mb proximal to the interval. The claim this
## supports is that the scan found the RIGHT ARM of chromosome III, at the
## resolution 231 phenotyped strains buy -- not that it found the gene.
##
## THE CLUSTER IS BELOW BONFERRONI. At -log10 p = 6.31 against a Bonferroni
## threshold of 6.97, the right-arm cluster clears only the eigen threshold.
## That is the threshold gwas_thresholds.R already calls "appropriate to the LD
## in the panel", and the argument here is explicitly a linkage-disequilibrium
## argument, so it is the consistent one to use -- but it must be stated, not
## slipped past.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
})

source("scripts/gwas_thresholds.R")

SCAN <- "supplemental_data/mapping/pos1_2023_gemma_loco.csv.gz"
OUT  <- "plots/diagnostics"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(SCAN))

HAS_GGTEXT <- requireNamespace("ggtext", quietly = TRUE)
panel_title <- function(letter) {
  if (HAS_GGTEXT)
    paste0("<span style='font-size:13pt;color:#111111'>**", letter, "**</span>")
  else letter
}
theme_pub <- function(base_size = 11) {
  title_el <- if (HAS_GGTEXT) ggtext::element_markdown(size = base_size + 0.5)
              else element_text(size = base_size + 2, face = "bold")
  theme_classic(base_size = base_size) +
    theme(axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          plot.title = title_el, plot.title.position = "plot",
          plot.subtitle = element_text(size = 8, colour = "grey35"),
          legend.key.size = grid::unit(9, "pt"))
}

## the numbers this script asserts, pinned as everywhere else in scripts/
PINNED <- list(singleton_ps = 5965738L, singleton_lp = 8.68,
               n_within_100kb = 628L, n_sig_within_100kb = 0L,
               median_lp_neighbour = 0.11,
               arm_n = 15L, arm_lp = 6.31,
               nil_n = 88L, nil_max_lp = 0.98,
               sid2_lp = 0.62, sid2_rank = 18662L)

th <- gwas_thresholds("pos1_2023")
d  <- fread(cmd = paste("gzcat", shQuote(SCAN)))
d[, lp := -log10(p_wald)]

## Rockman & Kruglyak (2009) low-recombination centre boundaries, Mb
CEN <- list(I = c(3.86, 11.04), II = c(4.88, 12.02), III = c(3.72, 10.34),
            IV = c(3.90, 12.79), V = c(5.20, 16.47), X = c(6.34, 12.72))
d[, domain := ifelse(ps / 1e6 >= sapply(chr, function(c) CEN[[c]][1]) &
                     ps / 1e6 <= sapply(chr, function(c) CEN[[c]][2]),
                     "centre", "arm")]

cat("== thresholds ==\n")
cat(sprintf("  Bonferroni %.2f | eigen (Li & Ji) %.2f | %d markers, M_eff %.0f\n\n",
            th$bonferroni, th$eigen, th$n_marker, th$m_eff))

## ---------------------------------------------------------------------------
## 1. local support for every marker that clears Bonferroni
## ---------------------------------------------------------------------------
support <- function(chrom, pos, win = 1e5) {
  w <- d[chr == chrom & ps >= pos - win & ps <= pos + win]
  o <- w[ps != pos]
  data.table(chr = chrom, Mb = round(pos / 1e6, 3),
             domain = d[chr == chrom & ps == pos, domain][1],
             lp = round(d[chr == chrom & ps == pos, max(lp)], 2),
             n_win = nrow(w), n_sig = sum(o$lp > th$eigen),
             max_nb = round(max(o$lp), 2), med_nb = round(median(o$lp), 2))
}
tops <- d[lp > th$bonferroni][order(-lp)]
sup  <- rbindlist(lapply(seq_len(nrow(tops)), function(i)
          support(tops$chr[i], tops$ps[i])))

cat("== every Bonferroni-significant marker, with its 100 kb neighbourhood ==\n")
cat("   n_win = markers within 100 kb | n_sig = of those, above the eigen line\n")
cat("   med_nb = median -log10 p of the neighbours\n\n")
print(as.data.frame(sup), row.names = FALSE)

lone <- sup[n_sig == 0]
cat("\n  ", nrow(lone), " of ", nrow(sup),
    " has no supporting marker at all: III:", format(PINNED$singleton_ps,
    big.mark = ","), "\n", sep = "")
stopifnot(nrow(lone) == 1, lone$Mb == 5.966)

## ---------------------------------------------------------------------------
## 2. is that unusual? every eigen-significant cluster genome-wide
## ---------------------------------------------------------------------------
sig <- d[lp > th$eigen][order(chr, ps)]
sig[, cl := cumsum(c(1, (diff(ps) > 1e5) |
                        (head(as.character(chr), -1) != tail(as.character(chr), -1))))]
cls <- sig[, .(chr = chr[1], domain = domain[1], n = .N,
               start_Mb = min(ps) / 1e6, end_Mb = max(ps) / 1e6,
               span_kb = (max(ps) - min(ps)) / 1e3, max_lp = max(lp)), by = cl]

cat("\n== the strongest signals genome-wide, ranked ==\n")
print(as.data.frame(cls[order(-max_lp)][1:10, .(chr, domain, n_markers = n,
      start_Mb = round(start_Mb, 3), end_Mb = round(end_Mb, 3),
      span_kb = round(span_kb, 1), max_lp = round(max_lp, 2))]), row.names = FALSE)
cat("\n  Of the seven signals above -log10 p = 6, six are clusters of 3 to 260\n",
    "  markers. One is a single marker: III:5.97 Mb. It is also the only one\n",
    "  of the seven in a chromosome CENTRE, where LD is long and support is\n",
    "  easiest to obtain.\n", sep = "")

cat(sprintf("\n  singletons among all %d eigen-significant markers: %d (%.1f%%)\n",
            nrow(sig), sum(cls$n == 1), 100 * sum(cls$n == 1) / nrow(sig)))

## ---------------------------------------------------------------------------
## 3. the right-arm cluster, and what it is not
## ---------------------------------------------------------------------------
## three small clusters sit between 12 and 13 Mb; the one meant here is the
## strongest, at 12.70-12.80 Mb -- select it by peak rather than by window
arm <- cls[chr == "III" & start_Mb > 12 & start_Mb < 13][which.max(max_lp)]
arm_peak_ps <- d[chr == "III" & ps >= arm$start_Mb * 1e6 &
                 ps <= arm$end_Mb * 1e6][which.max(lp), ps]
cat(sprintf("\n== the chromosome III right-arm cluster ==\n"))
cat(sprintf("  %d markers, %.3f-%.3f Mb (%.1f kb), peak -log10 p %.2f\n",
            arm$n, arm$start_Mb, arm$end_Mb, arm$span_kb, arm$max_lp))
cat(sprintf("  clears the eigen threshold (%.2f) but not Bonferroni (%.2f)\n",
            th$eigen, th$bonferroni))

NIL <- c(13.658e6, 13.695e6); SID2 <- 13680248
nil <- d[chr == "III" & ps >= NIL[1] & ps <= NIL[2]]
v   <- d[chr == "III" & ps == SID2]
c3  <- d[chr == "III"]
sid2_rank <- c3[, rank(-lp, ties.method = "min")][c3$ps == SID2]

cat(sprintf("\n  cluster peak marker: III:%s\n", format(arm_peak_ps, big.mark = ",")))
cat(sprintf("  distance from that peak to sid-2          : %.2f Mb\n",
            (SID2 - arm_peak_ps) / 1e6))
cat(sprintf("  distance from the cluster edge to the NIL interval: %.2f Mb\n",
            (NIL[1] - arm$end_Mb * 1e6) / 1e6))
cat(sprintf("  the 37 kb NIL interval holds %d markers, best -log10 p %.2f\n",
            nrow(nil), max(nil$lp)))
cat(sprintf("  sid-2 T96K (III:%s): -log10 p %.2f, rank %s of %s on chrIII\n",
            format(SID2, big.mark = ","), v$lp,
            format(sid2_rank, big.mark = ","), format(nrow(c3), big.mark = ",")))
cat("\n  So the scan reaches the right ARM, not the gene. The NIL interval\n",
    "  itself is flat.\n", sep = "")

stopifnot(abs(lone$lp - PINNED$singleton_lp) < 0.01,
          lone$n_win == PINNED$n_within_100kb,
          lone$n_sig == PINNED$n_sig_within_100kb,
          abs(lone$med_nb - PINNED$median_lp_neighbour) < 0.01,
          arm$n == PINNED$arm_n, abs(arm$max_lp - PINNED$arm_lp) < 0.01,
          nrow(nil) == PINNED$nil_n,
          abs(max(nil$lp) - PINNED$nil_max_lp) < 0.01,
          abs(v$lp - PINNED$sid2_lp) < 0.01,
          sid2_rank == PINNED$sid2_rank)

write_tsv(cls[order(-max_lp)], file.path(OUT, "TABLE_association_clusters.tsv"))

## ===========================================================================
## panels
## ===========================================================================
c3d <- d[chr == "III"]
mark <- tibble(Mb = c(5.965738, 12.718465),
               lab = c("III:5.97 Mb\n1 marker, no support",
                       "III:12.70–12.80 Mb\n15 markers"),
               lp = c(8.68, 6.31))

pA <- ggplot(c3d, aes(ps / 1e6, lp)) +
  annotate("rect", xmin = CEN$III[1], xmax = CEN$III[2], ymin = -Inf, ymax = Inf,
           fill = "grey92") +
  annotate("text", x = mean(CEN$III), y = Inf, label = "low-recombination centre",
           vjust = 1.5, size = 2.5, colour = "grey45") +
  geom_point(colour = "grey55", size = 0.5, alpha = 0.55) +
  geom_hline(yintercept = th$bonferroni, linetype = "dashed",
             linewidth = 0.35, colour = "grey25") +
  geom_hline(yintercept = th$eigen, linetype = "dotted",
             linewidth = 0.35, colour = "#B5446E") +
  geom_point(data = c3d[ps == 5965738], colour = "#B5446E", size = 2.2) +
  geom_point(data = c3d[ps >= 12.702e6 & ps <= 12.800e6 & lp > th$eigen],
             colour = "#3F7CAC", size = 1.6) +
  annotate("segment", x = 13.676, xend = 13.676, y = 3.3, yend = 1.3,
           linewidth = 0.3, colour = "grey30",
           arrow = arrow(length = unit(4, "pt"), type = "closed")) +
  annotate("text", x = 13.676, y = 3.7, label = "sid-2", size = 2.6,
           fontface = "italic", colour = "grey20") +
  labs(x = "Chromosome III position (Mb)", y = expression(-log[10]~italic(p)),
       title = panel_title("A"),
       subtitle = "dashed = Bonferroni, dotted = eigen; pink = the singleton, blue = the arm cluster") +
  theme_pub()

zoom <- function(lo, hi, letter, sub, col) {
  w <- c3d[ps >= lo & ps <= hi]
  ggplot(w, aes(ps / 1e6, lp)) +
    geom_hline(yintercept = th$eigen, linetype = "dotted",
               linewidth = 0.35, colour = "grey40") +
    geom_point(aes(colour = lp > th$eigen), size = 1.3, alpha = 0.85) +
    scale_colour_manual(values = c(`TRUE` = col, `FALSE` = "grey65"),
                        guide = "none") +
    labs(x = "Position (Mb)", y = expression(-log[10]~italic(p)),
         title = panel_title(letter), subtitle = sub) +
    theme_pub()
}
pB <- zoom(5965738 - 1e5, 5965738 + 1e5, "B",
           sprintf("%d markers within 100 kb, %d significant, median %.2f",
                   lone$n_win, lone$n_sig, lone$med_nb), "#B5446E")
pC <- zoom(12.65e6, 12.85e6, "C",
           sprintf("%d markers above the eigen line across %.0f kb",
                   arm$n, arm$span_kb), "#3F7CAC")

## width 0 would vanish on a log axis, so singletons are drawn at 0.3 kb and
## the axis is opened past it -- otherwise the very point the panel is about
## sits on the panel edge and is clipped
cls[, hi := fifelse(chr == "III" & start_Mb == 5.965738, "singleton",
             fifelse(chr == "III" & start_Mb == arm$start_Mb, "arm cluster",
                     domain))]
pD <- ggplot(cls, aes(pmax(span_kb, 0.3), max_lp)) +
  geom_hline(yintercept = th$bonferroni, linetype = "dashed",
             linewidth = 0.3, colour = "grey50") +
  geom_point(aes(size = n, fill = hi), shape = 21, stroke = 0.3,
             colour = "grey25", alpha = 0.9) +
  scale_x_log10(limits = c(0.2, 2000)) +
  scale_fill_manual(values = c(singleton = "#B5446E", `arm cluster` = "#3F7CAC",
                               centre = "#E8A33D", arm = "grey75"), name = NULL) +
  scale_size_continuous(range = c(1.4, 6), name = "markers") +
  annotate("text", x = 0.3, y = 8.68, label = "III:5.97 Mb  ", hjust = 0,
           vjust = -1.1, size = 2.6, colour = "#B5446E") +
  annotate("text", x = 97, y = 6.31, label = "III:12.7 Mb", hjust = 1.15,
           size = 2.6, colour = "#3F7CAC") +
  labs(x = "Cluster width (kb, log scale)", y = expression(peak~-log[10]~italic(p)),
       title = panel_title("D"),
       subtitle = "every eigen-significant cluster: the strong ones are wide, except one") +
  theme_pub() + theme(legend.position = "right")

fig <- (pA | pD) / (pB | pC)
pdf_path <- file.path(OUT, "chrIII_association_support.pdf")
dev <- if (isTRUE(capabilities("cairo"))) cairo_pdf else pdf
ggsave(pdf_path, fig, width = 9.4, height = 7.2, device = dev)
if (!file.exists(pdf_path)) ggsave(pdf_path, fig, width = 9.4, height = 7.2, device = pdf)
ggsave(file.path(OUT, "chrIII_association_support.png"), fig,
       width = 9.4, height = 7.2, dpi = 300, bg = "white")
stopifnot(file.exists(pdf_path))
cat("\nwrote chrIII_association_support.{pdf,png} and TABLE_association_clusters.tsv to ",
    OUT, "\n", sep = "")
