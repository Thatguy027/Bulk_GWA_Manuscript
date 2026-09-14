## The four published Baugh QTL intervals, zoomed, with within-interval sharing
##
## Webster et al. report four significant intervals:
##   Slope  IV:15,939,340-16,613,710   and  V:15,660,911-17,615,557
##   PC1    V:1,345,848-2,764,788      and  V:15,775,895-18,065,050
##
## The two right-arm chromosome V intervals overlap almost entirely (15.78-17.62
## Mb is common to both), so they are drawn separately but should be read as one
## region seen through two traits rather than two independent findings.
##
## Each interval is drawn once per trait family -- the four Slope traits in one
## panel, the four PC1 traits in the other -- so that "published against
## reconstructed" is a within-panel comparison. Lines are the strongest marker
## per 10 kb, which keeps the shape without drawing 26,000 points.
##
## The correlation panel is Spearman of -log10 p between traits computed on
## markers INSIDE each interval only. It asks whether the traits agree on the
## shape of the local signal, which is a different question from whether they
## each clear a genome-wide threshold, and a trait can score high here while
## detecting nothing.
##
## Reads from data/, so this is a diagnostic and outside the deposit-only build.
## Writes plots/diagnostics/DIAG_baugh_published_intervals.{pdf,png} and
## baugh_published_interval_table.tsv
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

D <- "data/baugh/reanalysis_mappings"; DIAG <- "plots/diagnostics"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)

IV <- data.table(
  id  = factor(c("Slope  IV:15.94-16.61", "Slope  V:15.66-17.62",
                 "PC1  V:1.35-2.76", "PC1  V:15.78-18.07"),
               levels = c("Slope  IV:15.94-16.61", "Slope  V:15.66-17.62",
                          "PC1  V:1.35-2.76", "PC1  V:15.78-18.07")),
  chr = c("IV", "V", "V", "V"),
  lo  = c(15939340, 15660911, 1345848, 15775895),
  hi  = c(16613710, 17615557, 2764788, 18065050))

FAM <- data.table(
  trait = c("published_slope_baugh","delta_slope_baugh","slope_nnls","delta_slope_nnls",
            "published_pc1_baugh","delta_pc1_baugh","pc1_nnls","delta_pc1_nnls"),
  fam   = rep(c("Slope traits","PC1 traits"), each = 4),
  src   = rep(c("published (MIP)","delta (MIP)","log-ratio (WGS)","delta (WGS)"), 2))
COL <- c("published (MIP)" = "#1A1A1A", "delta (MIP)" = "#E08214",
         "log-ratio (WGS)" = "#2E7BB6", "delta (WGS)" = "#1A7F5A")

A <- rbindlist(lapply(list.files(D, full.names = TRUE), function(f) {
  d <- fread(cmd = paste("gzcat", shQuote(f)))
  d[, .(trait = trait[1], chr, ps, af, p = p_wald, lp = -log10(p_wald))] }))
A <- A[!is.na(p)]
BF <- -log10(0.05 / A[trait == A$trait[1], .N])
A <- merge(A, FAM, by = "trait")
A[, fam := factor(fam, c("Slope traits", "PC1 traits"))]

## --- per-interval summary --------------------------------------------------
TB <- rbindlist(lapply(seq_len(nrow(IV)), function(i) {
  R <- A[chr == IV$chr[i] & ps >= IV$lo[i] & ps <= IV$hi[i]]
  R[, .(interval = IV$id[i], peak = ps[which.max(lp)], max_lp = max(lp),
        af = af[which.max(lp)], n_bf = sum(lp > BF), n_mark = .N),
    by = .(trait, fam, src)] }))
fwrite(TB, file.path(DIAG, "baugh_published_interval_table.tsv"), sep = "\t")
cat("peak and count per published interval\n")
print(TB[, .(interval, src, fam, peak, max_lp = round(max_lp, 2), af, n_bf)], nrows = 40)

## --- zoom traces -----------------------------------------------------------
Z <- rbindlist(lapply(seq_len(nrow(IV)), function(i) {
  pad <- (IV$hi[i] - IV$lo[i]) * 0.25
  R <- A[chr == IV$chr[i] & ps >= IV$lo[i] - pad & ps <= IV$hi[i] + pad]
  R[, .(lp = max(lp)), by = .(trait, fam, src, bin = (ps %/% 1e4) * 1e4)
    ][, interval := IV$id[i]][] }))
BND <- IV[, .(interval = id, lo, hi)]

theme_set(theme_bw(8.5) + theme(
  plot.title = element_text(face = "bold", size = 9.5),
  plot.subtitle = element_text(size = 7.3, colour = "grey30"),
  panel.grid.minor = element_blank(), legend.position = "top",
  legend.title = element_blank(), legend.key.height = unit(8, "pt"),
  strip.background = element_rect(fill = "grey93"),
  strip.text = element_text(size = 7.4, face = "bold")))

pA <- ggplot(Z, aes(bin / 1e6, lp, colour = src)) +
  geom_rect(data = BND, aes(xmin = lo / 1e6, xmax = hi / 1e6, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey70", alpha = 0.22) +
  geom_hline(yintercept = BF, linetype = 2, colour = "grey40", linewidth = 0.3) +
  geom_line(linewidth = 0.42, alpha = 0.9) +
  scale_colour_manual(values = COL) +
  facet_grid(fam ~ interval, scales = "free_x") +
  labs(title = "A  The four published intervals, every trait",
       subtitle = "Strongest marker per 10 kb. Grey band is the published interval, drawn with 25% flanking. Dashed = Bonferroni.",
       x = "Position (Mb)", y = "-log10 p")

## --- within-interval correlation -------------------------------------------
ORD <- FAM$trait
CC <- rbindlist(lapply(seq_len(nrow(IV)), function(i) {
  R <- A[chr == IV$chr[i] & ps >= IV$lo[i] & ps <= IV$hi[i]]
  M <- dcast(R, ps ~ trait, value.var = "lp")
  C <- cor(as.matrix(M[, ..ORD]), method = "spearman", use = "pairwise.complete.obs")
  d <- as.data.table(as.table(C)); setnames(d, c("a", "b", "r"))
  d[, interval := IV$id[i]][] }))
lab <- setNames(FAM$src, FAM$trait)
lab2 <- setNames(paste0(FAM$src, "\n", sub(" traits", "", FAM$fam)), FAM$trait)
CC[, `:=`(a = factor(a, ORD), b = factor(b, rev(ORD)))]
pB <- ggplot(CC, aes(a, b, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 1.85,
            colour = ifelse(CC$r > 0.75, "white", "grey15")) +
  scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2E4057",
                       midpoint = 0.5, limits = c(0, 1), guide = "none") +
  scale_x_discrete(labels = lab2) + scale_y_discrete(labels = lab2) +
  facet_wrap(~ interval, nrow = 1) +
  labs(title = "B  Do the traits agree on the shape of the local signal?",
       subtitle = "Spearman of -log10 p across markers INSIDE each interval. High here does not imply either trait detects the locus.",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5, lineheight = 0.9),
        axis.text.y = element_text(size = 5, lineheight = 0.9), panel.grid = element_blank())

fig <- pA / pB + plot_layout(heights = c(1, 1.15))
ggsave(file.path(DIAG, "DIAG_baugh_published_intervals.pdf"), fig, width = 12, height = 9.5)
ggsave(file.path(DIAG, "DIAG_baugh_published_intervals.png"), fig, width = 12, height = 9.5, dpi = 200)
cat(sprintf("\nwrote %s/DIAG_baugh_published_intervals.{pdf,png}\n", DIAG))
