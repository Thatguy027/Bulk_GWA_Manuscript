## How closely each reconstruction tracks the published traits, per interval --
##
## Spearman of -log10 p between a reconstruction and the published trait it
## reconstructs, computed on markers INSIDE each published interval only. Each
## reconstruction is compared only with its OWN published trait -- Slope
## reconstructions against published Slope, PC1 against published PC1. Neither
## the full trait-by-trait matrix nor the Slope-against-PC1 cross-comparisons
## are the question here.
##
## This measures agreement on the SHAPE of the local signal, which is a
## different thing from detection. A reconstruction can track the published
## trait closely inside an interval and still fail to clear a genome-wide
## threshold there -- the clearest case being PC1 V:1.35-2.76, where both pooled
## PC1 reconstructions sit near 0.8 against published PC1 while neither detects
## the locus.
##
## Reads from data/, so this is a diagnostic and outside the deposit-only build.
## Writes plots/diagnostics/DIAG_baugh_interval_correlations.{pdf,png} and
## baugh_interval_correlations.tsv
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(data.table); library(ggplot2)})

D <- "data/baugh/reanalysis_mappings"; DIAG <- "plots/diagnostics"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)

IVL <- c("Slope  IV:15.94-16.61", "Slope  V:15.66-17.62",
         "PC1  V:1.35-2.76", "PC1  V:15.78-18.07")
IV <- data.table(id = factor(IVL, levels = IVL), chr = c("IV", "V", "V", "V"),
                 lo = c(15939340, 15660911, 1345848, 15775895),
                 hi = c(16613710, 17615557, 2764788, 18065050))
CMPL <- c("delta Slope (MIP)", "log-ratio Slope (WGS)", "delta Slope (WGS)",
          "delta PC1 (MIP)",   "log-ratio PC1 (WGS)",   "delta PC1 (WGS)")
CMP <- data.table(trait = c("delta_slope_baugh","slope_nnls","delta_slope_nnls",
                            "delta_pc1_baugh","pc1_nnls","delta_pc1_nnls"),
                  cmp = factor(CMPL, levels = rev(CMPL)),
                  ref = rep(c("published_slope_baugh","published_pc1_baugh"), each = 3),
                  fam = rep(c("Slope traits","PC1 traits"), each = 3))

A <- rbindlist(lapply(list.files(D, full.names = TRUE), function(f) {
  d <- fread(cmd = paste("gzcat", shQuote(f)))
  d[, .(trait = trait[1], chr, ps, p = p_wald, lp = -log10(p_wald))] }))
A <- A[!is.na(p)]
BF <- -log10(0.05 / A[trait == A$trait[1], .N])

R <- rbindlist(lapply(seq_len(nrow(IV)), function(i) {
  W <- dcast(A[chr == IV$chr[i] & ps >= IV$lo[i] & ps <= IV$hi[i]],
             ps ~ trait, value.var = "lp")
  CMP[, .(interval = IV$id[i], cmp, fam,
          r = mapply(function(t, rf)
                cor(W[[t]], W[[rf]], method = "spearman", use = "pairwise.complete.obs"),
                trait, ref),
          detects = sapply(trait, function(t) sum(W[[t]] > BF) > 0))] }))
fwrite(R, file.path(DIAG, "baugh_interval_correlations.tsv"), sep = "\t")
cat("correlation to the published traits, inside each interval\n")
print(dcast(R, cmp ~ interval, value.var = "r")[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x, 2) else x)], nrows = 10)

theme_set(theme_bw(9) + theme(
  plot.title = element_text(face = "bold", size = 10),
  plot.subtitle = element_text(size = 7.6, colour = "grey30"),
  panel.grid.minor = element_blank(), legend.position = "top",
  legend.title = element_blank(), strip.background = element_rect(fill = "grey93"),
  strip.text = element_text(size = 7.4, face = "bold")))

p <- ggplot(R, aes(r, cmp, colour = fam, shape = detects)) +
  geom_vline(xintercept = c(0.5, 0.8), linetype = 3, colour = "grey72", linewidth = 0.3) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.2f", r)), vjust = -1.15, size = 2.4,
            show.legend = FALSE) +
  scale_colour_manual(values = c(`Slope traits` = "#2E4057",
                                 `PC1 traits`   = "#C4302B")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "clears Bonferroni in the interval",
                                `FALSE` = "does not")) +
  scale_x_continuous(limits = c(0, 1.04), breaks = seq(0, 1, 0.25)) +
  facet_wrap(~ interval, nrow = 1) +
  labs(title = "How closely each reconstruction tracks the published traits, inside each interval",
       subtitle = paste("Spearman of -log10 p across markers within the published interval, each reconstruction against its OWN published trait.",
                        "\nFilled points clear Bonferroni somewhere in that interval, open points do not -- tracking the published signal and detecting it are different things."),
       x = "Spearman correlation to the matching published trait", y = NULL)
ggsave(file.path(DIAG, "DIAG_baugh_interval_correlations.pdf"), p, width = 11.5, height = 4.2)
ggsave(file.path(DIAG, "DIAG_baugh_interval_correlations.png"), p, width = 11.5, height = 4.2, dpi = 200)
cat(sprintf("\nwrote %s/DIAG_baugh_interval_correlations.{pdf,png}\n", DIAG))
