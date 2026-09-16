## The two LOD conventions, on the pos-1 contrast of both crosses -------------
##
##   Rscript scripts/DIAG_cross_lod_conventions.R
##     -> plots/diagnostics/DIAG_cross_lod_conventions.{pdf,png}
##     -> plots/diagnostics/cross_lod_conventions.tsv
##
## The cross exports ship a LOD column, and xqtl_pipeline.R's recompute_lod()
## can regenerate it under either of two conventions that differ by exactly
## log10(2) = 0.30103:
##
##   package  what the exports use. PvalToLOD() expects a ONE-tailed p and
##            doubles its argument, but calcContrastStats hands it a TWO-tailed
##            p, so the result sits ~0.3 LOD below the textbook value.
##   chisq    the textbook 1-df LOD, z^2 / (2 ln 10).
##
## THE SHIPPED COLUMN IS ALREADY THE LOG-SPACE VALUE. This is the thing to be
## clear about, because "shipped versus recomputed" invites the wrong reading.
## The recomputation happened inside the export pipeline BEFORE the files were
## written, so comparing the shipped LOD to recompute_lod() is a reproducibility
## check, not a before-and-after. It agrees to 1.8e-12, which is worth knowing
## but is not the interesting comparison.
##
## The real before-and-after is against what the ORIGINAL pipeline produced.
## xqtl_stats.R used calcContrastStats(), which routes through PvalToLOD() on
## the LINEAR-scale p. That p underflows to exactly 0 above |z| ~ 38.5, and
## PvalToLOD(0) is qchisq(0, df = 1, lower.tail = FALSE) = Inf. So the original
## LOD is INFINITE for every marker above |z| = 38.5 -- 4,765 of them in the
## N2 x XZ1516 pos-1 contrast, whose log-space LOD runs from 321 to 710. Any
## shipped LOD above about 322 is therefore proof in itself that the log-space
## recomputation was applied: it could not have come from the p column.
##
## The third panel is where the two CONVENTIONS stop being a constant apart.
## PvalToLOD() maps |z| <= qnorm(0.25, lower = FALSE) = 0.6745 to 0 because
## qchisq's doubled argument exceeds 1 there, so below that the package LOD is
## pinned at zero while the chisq LOD is not.
##
## Reads the cross exports under data/, so it stays a DIAG script.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

CE   <- "data/cross_experiments"
DIAG <- "plots/diagnostics"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

SRC <- list(
  `N2 x XZ1516` = file.path(CE, "N2-XZ_export/plot_data",
    "N2_XZ1516_F2-2_contrast_HT115-pos1_10000_plot_DF.tsv.gz"),
  `JU1793 x JU2466` = file.path(CE, "JU1793-JU2466_export/plot_data",
    "JU1793_JU2466_F2-2_contrast_HT115g-POS1g_10000_plot_DF.tsv.gz"))
for (f in SRC) if (!file.exists(f)) stop("missing export: ", f, call. = FALSE)

## recompute_lod(), verbatim in substance from
## data/cross_experiments/N2-XZ_export/scripts/xqtl_pipeline.R
lod_from_z <- function(z, convention = c("package", "chisq")) {
  convention <- match.arg(convention)
  log_p <- log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
  if (convention == "chisq") return(z^2 / (2 * log(10)))
  l <- suppressWarnings(
    qchisq(log(2) + log_p, df = 1, lower.tail = FALSE, log.p = TRUE)) / (2 * log(10))
  l[is.nan(l)] <- 0
  l
}

D <- rbindlist(lapply(names(SRC), function(nm) {
  d <- fread(cmd = paste("gzcat", shQuote(SRC[[nm]])),
             select = c("chrom", "physical.position", "z", "p", "LOD"))
  d[, cross := nm]
  d[is.finite(z)]
}))
D[, `:=`(lod_package = lod_from_z(z, "package"),
         lod_chisq   = lod_from_z(z, "chisq"))]
## what PvalToLOD() does to the LINEAR-scale p the exports also ship: the
## original path, reproduced so the failure can be seen rather than described
D[, lod_naive := suppressWarnings(
    qchisq(pmin(2 * p, 1), df = 1, lower.tail = FALSE) / (2 * log(10)))]
msg(sprintf("markers: %s", paste(D[, .N, by = cross][, sprintf("%s %s", cross,
            format(N, big.mark = ","))], collapse = " | ")))

## --- does the shipped column reproduce? ------------------------------------
chk <- D[, .(n = .N,
             max_abs_diff = max(abs(LOD - lod_package)),
             max_rel_diff = max(abs(LOD - lod_package) / pmax(LOD, 1e-12)),
             cor = cor(LOD, lod_package)), by = cross]
cat("\n== shipped LOD against the package-convention recomputation ==\n")
cat("   (a reproducibility check -- the shipped column is ALREADY log space)\n")
print(chk)
cat("\n== the real before-and-after: PvalToLOD() on the linear p ==\n")
print(D[, .(n = .N,
            naive_infinite = sum(!is.finite(lod_naive)),
            pct = round(100 * mean(!is.finite(lod_naive)), 2),
            logspace_LOD_of_those = if (any(!is.finite(lod_naive)))
              sprintf("%.0f - %.0f", min(lod_package[!is.finite(lod_naive)]),
                      max(lod_package[!is.finite(lod_naive)])) else "-",
            max_finite_naive = round(max(lod_naive[is.finite(lod_naive)]), 1)), by = cross])

cat("\n== the two conventions ==\n")
print(D[, .(n = .N,
            median_gap = round(median(lod_chisq - lod_package), 5),
            gap_at_z_gt_2 = round(median((lod_chisq - lod_package)[abs(z) > 2]), 5),
            log10_2 = round(log10(2), 5),
            n_package_zero = sum(lod_package == 0),
            max_abs_z = round(max(abs(z)), 1)), by = cross])
cat("\n  p underflowed to exactly 0 (why LOD is recomputed in log space at all):\n")
print(D[, .(p_exactly_zero = sum(p == 0),
            pct = round(100 * mean(p == 0), 2),
            min_z_where_p_is_0 = if (any(p == 0)) round(min(abs(z)[p == 0]), 1) else NA_real_,
            max_abs_z = round(max(abs(z)), 1)), by = cross])
fwrite(D[, .(cross, chrom, physical.position, z, p, LOD, lod_package, lod_chisq)][
         order(cross, chrom, physical.position)],
       file.path(DIAG, "cross_lod_conventions.tsv.gz"), sep = "\t")

## --- how fast does the offset actually converge? ----------------------------
## It is log10(2) only in the limit. Worth pinning because the obvious reading of
## "the conventions differ by 0.30103" is that a threshold can be shifted by that
## constant, and near the threshold it cannot.
zz <- seq(0.5, 20, by = 0.001)
pk <- lod_from_z(zz, "package"); gp <- lod_from_z(zz, "chisq") - pk
gap_at_thr <- gp[which.min(abs(pk - 3.57))]
cat("\n== the offset is asymptotic, not constant ==\n")
for (g in c(0.29, 0.295, 0.30)) {
  i <- which(gp >= g)[1]
  cat(sprintf("  reaches %.3f at package LOD %6.2f (|z| = %.2f)\n", g, pk[i], zz[i]))
}
cat(sprintf("  at the genome-wide threshold LOD 3.57 the gap is %.4f, not %.5f\n",
            gap_at_thr, log10(2)))

## --- figure -----------------------------------------------------------------
theme_set(theme_bw(9) + theme(
  plot.title = element_text(face = "bold", size = 9.5),
  plot.subtitle = element_text(size = 7.2, colour = "grey30"),
  panel.grid.minor = element_blank(), legend.position = "bottom",
  strip.background = element_rect(fill = "grey93"),
  strip.text = element_text(size = 7.3, face = "bold")))
COL <- c(`N2 x XZ1516` = "#1B7837", `JU1793 x JU2466` = "#C4302B")
## 522k + 154k points; thin for the scatter, keep everything for the statistics
set.seed(1)
S <- D[, .SD[sample(.N, min(.N, 40000))], by = cross]

## Panel A: the ORIGINAL path against the log-space one. Infinite naive values
## cannot be drawn, so they are pinned to the top of the panel and marked.
ZCUT <- 38.5
CAP <- max(S$lod_package) * 1.06
A <- copy(S)[, `:=`(naive_plot = fifelse(is.finite(lod_naive), lod_naive, CAP),
                    blown = !is.finite(lod_naive))]
nblow <- A[, .(n = sum(blown)), by = cross][n > 0]
pA <- ggplot(A, aes(lod_package, naive_plot, colour = cross)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4, colour = "grey45") +
  geom_hline(yintercept = CAP, linewidth = 0.3, colour = "grey60", linetype = "dotted") +
  geom_point(aes(shape = blown), size = 0.5, alpha = 0.4) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), guide = "none") +
  ## count from the FULL data, not the thinned scatter -- S is a 40k sample per
  ## cross, so sum(A$blown) would report the sample's share and not the truth
  annotate("text", x = 0, y = CAP, hjust = 0, vjust = -0.5, size = 2.5, colour = "grey25",
           label = sprintf("Inf -- %s of %s markers, pinned here to be visible",
                           format(D[!is.finite(lod_naive), .N], big.mark = ","),
                           format(nrow(D), big.mark = ","))) +
  scale_colour_manual(values = COL, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "log-space LOD (what the exports ship)",
       y = "PvalToLOD() on the linear p\n(what the original pipeline gave)",
       title = "The real before-and-after: where the original LOD died",
       subtitle = paste(sprintf("HT115 vs pos-1. Above |z| = %.1f the shipped p column is exactly 0 and PvalToLOD(0) is Inf, so the", ZCUT),
                        "\noriginal LOD is infinite for every marker whose log-space LOD exceeds about 322. The N2 x XZ1516",
                        "\ncontrast has 4,765 of them, running to LOD 710; the JU cross has none, peaking at |z| = 25.4.",
                        "\nA shipped LOD above 322 is itself proof the log-space recomputation was applied."))

## Plotted as a DIFFERENCE. On a 0-700 axis a 0.3 offset is a fraction of a
## pixel, so a chisq-against-package scatter is the identity line and shows
## nothing; the gap has to be the y axis to be visible at all.
pB <- ggplot(S, aes(lod_package + 0.1, lod_chisq - lod_package, colour = cross)) +
  geom_hline(yintercept = log10(2), linewidth = 0.4, colour = "grey45",
             linetype = "dashed") +
  geom_point(size = 0.35, alpha = 0.35) +
  geom_vline(xintercept = 3.57, linewidth = 0.4, colour = "#C4302B") +
  annotate("text", x = min(S$lod_package) + 0.1, y = log10(2),
           label = "log10(2) = 0.30103 ", hjust = 0, vjust = -0.6, size = 2.5,
           colour = "grey30") +
  annotate("text", x = 3.57, y = 0.05, label = " genome-wide threshold, LOD 3.57",
           hjust = 0, size = 2.4, colour = "#C4302B") +
  scale_x_log10() +
  scale_colour_manual(values = COL, guide = "none") +
  labs(x = "package-convention LOD (log, +0.1)", y = "chisq minus package",
       title = "The conventions differ by log10(2) ASYMPTOTICALLY, and converge slowly",
       subtitle = paste(sprintf("The gap is %.4f at the genome-wide threshold of LOD 3.57, not 0.30103. It reaches 0.290 at LOD 5.2,",
                                gap_at_thr),
                        "\n0.295 at LOD 10.1 and 0.300 only at LOD 62.7. So for a marker near the threshold the two conventions",
                        "\nare 0.286 apart, and treating the offset as a flat 0.301 misstates it by 0.015 LOD."))

pC <- ggplot(S[abs(z) < 4], aes(abs(z), lod_chisq - lod_package, colour = cross)) +
  geom_hline(yintercept = log10(2), linewidth = 0.4, colour = "grey45", linetype = "dashed") +
  geom_vline(xintercept = qnorm(0.25, lower.tail = FALSE), linewidth = 0.4,
             colour = "grey45") +
  geom_point(size = 0.35, alpha = 0.35) +
  annotate("text", x = qnorm(0.25, lower.tail = FALSE), y = Inf,
           label = " |z| = 0.6745", hjust = 0, vjust = 1.6, size = 2.4, colour = "grey30") +
  scale_colour_manual(values = COL, guide = "none") +
  labs(x = "|z|", y = "chisq minus package",
       title = "Where the offset stops being constant, and why cross does not matter here",
       subtitle = paste("Below |z| = 0.6745 the doubled argument to qchisq exceeds 1, PvalToLOD() returns 0, and the",
                        "\npackage LOD is pinned there while the chisq LOD is not. Nothing significant lives here. The two",
                        "\ncrosses lie on ONE curve because the gap is a deterministic function of |z| -- colour separates",
                        "\nthem only in panel A, where the LOD ranges differ."))

fig <- pA / pB / pC
ggsave(file.path(DIAG, "DIAG_cross_lod_conventions.pdf"), fig, width = 7.5, height = 11)
ggsave(file.path(DIAG, "DIAG_cross_lod_conventions.png"), fig, width = 7.5, height = 11, dpi = 200)
msg("wrote DIAG_cross_lod_conventions.{pdf,png}")
