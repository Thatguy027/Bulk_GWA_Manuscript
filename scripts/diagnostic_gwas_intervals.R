## Diagnostic -- what a GWAS QTL interval should and should not capture -------
##
##   Rscript scripts/diagnostic_gwas_intervals.R
##     -> plots/diagnostics/gwas_interval_diagnostic.{pdf,png}
##     -> plots/diagnostics/TABLE_gwas_local_support.tsv
##
## NOT A MANUSCRIPT FIGURE. This exists to settle one decision: how to define a
## QTL interval from the pooled pos-1 association scan so that an unsupported
## single marker is excluded while a supported but sub-Bonferroni region is
## kept.
##
##   A  chromosome III, every marker, with both thresholds
##   B  the isolated marker at 5.966 Mb, +/- 150 kb -- nothing around it
##   C  the supported cluster at 12.70-12.80 Mb, +/- 150 kb
##   D  local support against significance for every marker above the eigen
##      line, genome-wide: the separating statistic
##
## THE PROBLEM. On chromosome III one marker clears Bonferroni (5.966 Mb,
## -log10 p 8.68, allele frequency 0.091) and has ZERO other eigen-passing
## markers within 100 kb, out of 628 markers present in that window. A cluster
## at 12.70-12.80 Mb peaks at 6.31 -- short of Bonferroni, above eigen -- and
## carries 14 eigen-passing markers within 100 kb. The eye reads the second as
## a QTL and the first as noise, and a threshold alone cannot tell them apart.
##
## THE SEPARATING STATISTIC is local support: how many OTHER markers within a
## window also clear the threshold. In a panel with linkage disequilibrium a
## true association is tagged by several correlated markers; a lone spike with
## flat neighbours is genotyping error, an unshared rare haplotype, or chance.
## This is the clumping logic used in human GWAS, applied here to a panel whose
## Bonferroni threshold is known to be over-conservative because the markers are
## not independent.
##
## WHAT THIS SCRIPT DOES NOT DO. It does not draw intervals. Interval ENDS need
## linkage disequilibrium to the peak marker, which needs the genotype matrix;
## the recommendation is to take the span of markers with r-squared >= 0.5 to the
## peak, which is what makes a GWAS interval commensurable with a linkage
## interval. This script establishes only that local support is the right
## admission criterion, and what cutoff separates the two cases cleanly.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(patchwork); library(ggtext)
})

OUT   <- "plots/diagnostics"
SCAN  <- "supplemental_data/mapping/pos1_2023_gemma_loco.csv.gz"
EIG   <- "supplemental_data/mapping/eigen_independent_tests.tsv"
stopifnot(file.exists(SCAN))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

BONF  <- 6.97      # Bonferroni over all 464,045 markers
EIGEN <- 4.60      # alpha / 1,972 effective independent tests (Li & Ji)
WIN   <- 1e5       # the support window, +/- 100 kb

COL_SPUR <- "#B23A48"   # the isolated marker
COL_REAL <- "#1B6C7A"   # the supported cluster
COL_BG   <- "grey72"

panel_title <- function(l)
  paste0("<span style='font-size:13pt;color:#111111'>**", l, "**</span>")
theme_pub <- function(base = 11) {
  theme_classic(base_size = base) +
    theme(axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_markdown(size = base + 0.5),
          plot.title.position = "plot",
          legend.key.size = grid::unit(9, "pt"))
}

d <- fread(SCAN)[, .(chr, ps, af, beta, p = p_wald)]
d[, lp := -log10(p)]
setorder(d, chr, ps)
cat("markers:", nrow(d), "| above eigen:", sum(d$lp > EIGEN),
    "| above Bonferroni:", sum(d$lp > BONF), "\n")

## ---- local support, genome-wide -------------------------------------------
## for each marker above the eigen line, how many OTHER markers within WIN of
## it also clear that line. Done per chromosome so windows never span
## chromosomes.
hi <- d[lp > EIGEN]
hi[, support := {
  s <- integer(.N)
  for (i in seq_len(.N)) s[i] <- sum(abs(ps - ps[i]) <= WIN) - 1L
  s
}, by = chr]
## and how many markers of any kind sit in that window, so a low support count
## cannot be blamed on thin coverage
hi[, n_window := {
  dd <- d[chr == .BY$chr]
  sapply(ps, function(p) sum(abs(dd$ps - p) <= WIN))
}, by = chr]

fwrite(hi[order(-lp), .(chr, ps, Mb = round(ps / 1e6, 4), lp = round(lp, 3),
                        af = round(af, 4), beta = signif(beta, 4),
                        support, n_window)],
       file.path(OUT, "TABLE_gwas_local_support.tsv"), sep = "\t")

cat("\n== every marker above the eigen line, by local support ==\n")
print(as.data.frame(hi[, .(markers = .N, median_lp = round(median(lp), 2),
                           median_af = round(median(af), 3)),
                       by = .(isolated = support == 0)][order(isolated)]),
      row.names = FALSE)
cat("\nisolated markers (support == 0), all chromosomes:\n")
print(as.data.frame(hi[support == 0][order(-lp),
      .(chr, Mb = round(ps / 1e6, 3), lp = round(lp, 2), af = round(af, 3),
        n_window)]), row.names = FALSE)

cat("\n== the two chromosome III cases side by side ==\n")
cases <- rbindlist(list(
  hi[chr == "III"][which.max(lp), .(case = "isolated 5.966 Mb", Mb = ps / 1e6,
                                    lp, af, support, n_window)],
  hi[chr == "III" & ps / 1e6 > 12.6 & ps / 1e6 < 12.9][which.max(lp),
     .(case = "cluster 12.7 Mb", Mb = ps / 1e6, lp, af, support, n_window)]))
print(as.data.frame(cases[, .(case, Mb = round(Mb, 3), lp = round(lp, 2),
                              af = round(af, 3), support, n_window)]),
      row.names = FALSE)

## does a support cutoff separate them, and is it sharp?
cat("\n== admission at support >= k, chromosome III ==\n")
for (k in c(1, 2, 3, 5)) {
  keep <- hi[chr == "III" & support >= k]
  cat(sprintf("  k = %d: %2d markers admitted | 5.966 in? %-5s | 12.7 cluster in? %s\n",
      k, nrow(keep), any(abs(keep$ps / 1e6 - 5.966) < 0.001),
      any(keep$ps / 1e6 > 12.6 & keep$ps / 1e6 < 12.9)))
}

## ---- A: chromosome III --------------------------------------------------
c3 <- d[chr == "III"][, Mb := ps / 1e6]
mark <- data.frame(Mb = c(5.966, 12.718), lp = c(8.68, 6.31),
                   lab = c("isolated", "supported cluster"),
                   col = c(COL_SPUR, COL_REAL))
pA <- ggplot(c3, aes(Mb, lp)) +
  geom_point(shape = 16, size = 0.35, colour = COL_BG, alpha = 0.55) +
  geom_hline(yintercept = BONF, linetype = "dashed", linewidth = 0.35,
             colour = "grey35") +
  geom_hline(yintercept = EIGEN, linetype = "dotted", linewidth = 0.4,
             colour = "#2E7D32") +
  geom_point(data = c3[lp > EIGEN], shape = 16, size = 0.9, colour = "grey25") +
  geom_point(data = mark, aes(Mb, lp), colour = mark$col, size = 2.2) +
  annotate("text", x = 0.15, y = BONF + 0.28, label = "Bonferroni 6.97",
           hjust = 0, size = 2.5, colour = "grey35") +
  annotate("text", x = 0.15, y = EIGEN + 0.28, label = "eigen 4.60",
           hjust = 0, size = 2.5, colour = "#2E7D32") +
  annotate("text", x = 5.966, y = 8.68, label = "isolated", vjust = -0.9,
           size = 2.6, colour = COL_SPUR, fontface = "bold") +
  annotate("text", x = 12.718, y = 6.31, label = "supported", vjust = -0.9,
           hjust = 0.9, size = 2.6, colour = COL_REAL, fontface = "bold") +
  scale_y_continuous(limits = c(0, 9.8), breaks = 0:9) +
  labs(x = "Chromosome III (Mb)", y = "−log10 p", title = panel_title("A")) +
  theme_pub()

## ---- B and C: the two windows ------------------------------------------
zoom <- function(centre, letter, col, lab) {
  w <- c3[Mb >= centre - 0.15 & Mb <= centre + 0.15]
  ggplot(w, aes(Mb, lp)) +
    geom_hline(yintercept = BONF, linetype = "dashed", linewidth = 0.3,
               colour = "grey35") +
    geom_hline(yintercept = EIGEN, linetype = "dotted", linewidth = 0.4,
               colour = "#2E7D32") +
    geom_point(shape = 16, size = 0.7, colour = COL_BG) +
    geom_point(data = w[lp > EIGEN], shape = 16, size = 1.5, colour = col) +
    annotate("text", x = centre - 0.148, y = 9.2,
             label = sprintf("%s\n%d markers in window, %d above eigen",
                             lab, nrow(w), sum(w$lp > EIGEN)),
             hjust = 0, vjust = 1, size = 2.5, colour = "grey20", lineheight = 1.1) +
    scale_y_continuous(limits = c(0, 9.8), breaks = 0:9) +
    labs(x = "Chromosome III (Mb)", y = "−log10 p",
         title = panel_title(letter)) +
    theme_pub()
}
pB <- zoom(5.966, "B", COL_SPUR, "isolated marker")
pC <- zoom(12.750, "C", COL_REAL, "supported cluster")

## ---- D: the separating statistic ---------------------------------------
hi[, kind := ifelse(support == 0, "isolated", "supported")]
pD <- ggplot(hi, aes(support + 1, lp)) +
  geom_hline(yintercept = BONF, linetype = "dashed", linewidth = 0.3,
             colour = "grey35") +
  geom_point(aes(fill = af), shape = 21, size = 2, stroke = 0.25,
             colour = "grey25",
             position = position_jitter(width = 0.06, height = 0, seed = 1)) +
  scale_x_log10(breaks = c(1, 2, 3, 6, 11, 21, 51),
                labels = c(0, 1, 2, 5, 10, 20, 50)) +
  scale_fill_viridis_c(option = "mako", direction = -1, end = 0.9,
                       name = "Allele\nfrequency") +
  labs(x = paste0("Other markers above the eigen line within ",
                  WIN / 1e3, " kb"),
       y = "−log10 p", title = panel_title("D")) +
  theme_pub() +
  theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 7),
        legend.background = element_rect(fill = alpha("white", 0.85),
                                         colour = NA))

fig <- (pA | pD) / (pB | pC)
ggsave(file.path(OUT, "gwas_interval_diagnostic.pdf"), fig,
       width = 10.4, height = 7, device = cairo_pdf)
ggsave(file.path(OUT, "gwas_interval_diagnostic.png"), fig,
       width = 10.4, height = 7, dpi = 300, bg = "white")
cat("\nwrote gwas_interval_diagnostic.{pdf,png} and TABLE_gwas_local_support.tsv",
    "to", OUT, "\n")
