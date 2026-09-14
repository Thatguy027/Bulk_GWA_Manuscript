## What predicts per-strain NNLS error against MIP-seq? -----------------------
##
## Per-strain RMSD between the deconvolved and the published MIP-seq frequency
## has been reported throughout as a quality measure without ever being
## explained. This tests it against the obvious candidates: how abundant the
## strain is, how genetically close its nearest neighbour in the reference is,
## and how often the solver puts it at exactly zero.
##
## Both an absolute and a relative error are used, because they answer different
## questions. RMSD in frequency units is what propagates into a trait; RMSD
## divided by the strain's own mean frequency is what says whether the estimate
## is proportionally reliable, and the two point in OPPOSITE directions with
## abundance.
##
## The strongest predictor turns out to be the one that was not obvious: the
## number of markers at which a strain is the ONLY carrier of the alternate
## allele. That is what makes a strain identifiable to a non-negative least
## squares fit at all -- with no private marker, its column is a combination of
## others and the solver has nothing to anchor it with. Counts are deposited as
## baugh_strain_private_markers.tsv (computed from the deposited genotype matrix
## plus the grafted PB306, 1,237,106 markers, 103 strains).
##
## Reads from supplemental_data, so this one does stay inside the deposit-only
## rebuild.
##
## Writes plots/diagnostics/DIAG_baugh_rmsd_predictors.{pdf,png} and
## baugh_rmsd_predictors.tsv
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})

DEC <- "supplemental_data/deconvolution"; DIAG <- "plots/diagnostics"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)

nn <- fread(cmd = paste("gzcat", shQuote(file.path(DEC, "baugh_nnls_dep103_with_mipseq.tsv.gz"))))
d <- nn[!is.na(frq) & !is.na(published_frq) & strain != "N2",
        .(rmsd = sqrt(mean((frq - published_frq)^2)),
          bias = mean(frq - published_frq),
          mean_f = mean(published_frq),
          zero_frac = mean(frq == 0), n = .N), by = strain]
sim <- fread(file.path(DEC, "baugh_strain_similarity.tsv"))
d <- merge(d, sim, by = "strain", all.x = TRUE)
priv <- fread(file.path(DEC, "baugh_strain_private_markers.tsv"))
d <- merge(d, priv[, .(strain, n_private, n_alt)], by = "strain", all.x = TRUE)
d[, rel := rmsd / mean_f]
fwrite(d, file.path(DIAG, "baugh_rmsd_predictors.tsv"), sep = "\t")

VARS <- data.table(
  v    = c("n_private", "mean_f", "nn_ibs_wild", "zero_frac"),
  lab  = c("markers where the strain is the only alt carrier (log10)",
           "mean MIP-seq frequency", "identity by state to nearest wild neighbour",
           "fraction of samples the solver sets to zero"),
  ## private-marker counts span 0 to 154,000, so that one is plotted on a log
  ## axis; the rest are readable linear
  logx = c(TRUE, FALSE, FALSE, FALSE))

mk <- function(yv, ylab, logy) rbindlist(lapply(seq_len(nrow(VARS)), function(i) {
  x <- d[[VARS$v[i]]]
  data.table(strain = d$strain, x = if (VARS$logx[i]) log10(x + 1) else x,
             y = d[[yv]], panel = factor(VARS$lab[i], VARS$lab)) }))
L <- rbind(mk("rmsd", "", FALSE)[, metric := "RMSD (frequency units)"],
           mk("rel",  "", FALSE)[, metric := "RMSD / mean frequency"])
L[, metric := factor(metric, c("RMSD (frequency units)", "RMSD / mean frequency"))]
L <- L[is.finite(x) & is.finite(y)]

ANN <- L[, .(r = cor(x, y, method = "spearman"),
             p = cor.test(x, y, method = "spearman")$p.value), by = .(panel, metric)]
ANN[, lab := sprintf("rho = %+.2f%s", r, ifelse(p < 0.001, "  (p < 0.001)", sprintf("  (p = %.3f)", p)))]

theme_set(theme_bw(9) + theme(
  plot.title = element_text(face = "bold", size = 10),
  plot.subtitle = element_text(size = 7.6, colour = "grey30"),
  panel.grid.minor = element_blank(),
  strip.background = element_rect(fill = "grey93"),
  strip.text = element_text(size = 7.3, face = "bold")))

hl <- unique(L[, .SD[order(-y)][1:3], by = .(panel, metric)])
hl[, y := y * 0.94]   ## nudge labels clear of the annotation box
p <- ggplot(L, aes(x, y)) +
  geom_point(colour = "#2E4057", alpha = 0.7, size = 1.5) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              colour = "#C4302B", linewidth = 0.5) +
  geom_text_repel(data = hl, aes(label = strain), size = 2.2, colour = "grey30",
                  min.segment.length = 0, max.overlaps = Inf, seed = 1) +
  geom_label(data = ANN, aes(x = -Inf, y = Inf, label = lab), hjust = -0.04, vjust = 1.15,
             size = 2.4, colour = "grey15", inherit.aes = FALSE,
             fill = "white", alpha = 0.85, label.size = 0,
             label.padding = unit(1.6, "pt")) +
  scale_y_log10() +
  facet_grid(metric ~ panel, scales = "free", switch = "y") +
  labs(title = "What predicts a strain's deconvolution error?",
       subtitle = paste("Per-strain RMSD of the NNLS frequency against published MIP-seq, 99 strains, N2 excluded.",
                        "\nTop row absolute error, bottom row error relative to the strain's own abundance.",
                        "Red line is a loess fit. Both axes log where the spread demands it."),
       x = NULL, y = NULL)
ggsave(file.path(DIAG, "DIAG_baugh_rmsd_predictors.pdf"), p, width = 10, height = 6.4)
ggsave(file.path(DIAG, "DIAG_baugh_rmsd_predictors.png"), p, width = 10, height = 6.4, dpi = 200)

cat("Spearman of each predictor against the two error measures\n")
print(dcast(ANN, panel ~ metric, value.var = "r")[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x, 3) else x)])
ok <- is.finite(d$nn_ibs_wild) & is.finite(d$n_private)
cat("\nrank models for absolute RMSD\n")
m0 <- lm(rank(rmsd) ~ rank(mean_f) + rank(nn_ibs_wild) + rank(zero_frac), data = d[ok])
m1 <- lm(rank(rmsd) ~ rank(mean_f) + rank(nn_ibs_wild) + rank(zero_frac) +
                      rank(n_private), data = d[ok])
cat(sprintf("  without private-marker count: adj R2 = %.3f\n", summary(m0)$adj.r.squared))
cat(sprintf("  with    private-marker count: adj R2 = %.3f\n", summary(m1)$adj.r.squared))
print(round(summary(m1)$coefficients[, c(1, 4)], 4))
cat(sprintf("\nwrote %s/DIAG_baugh_rmsd_predictors.{pdf,png}\n", DIAG))
