## Signal sharing across the eight Baugh trait scans -------------------------
##
## GEMMA LOCO scans of the eight columns of baugh_mapping_traits.csv, read from
## data/baugh/reanalysis_mappings/. Two questions: how much signal the traits
## share, and how close the pooled-WGS traits get to a scan on the published
## phenotypes.
##
## Reads from data/, not supplemental_data/, so this is a diagnostic and is not
## part of the deposit-only rebuild.
##
## Writes plots/diagnostics/DIAG_baugh_reanalysis_mappings.{pdf,png} and
## baugh_reanalysis_locus_table.tsv
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})

D    <- "data/baugh/reanalysis_mappings"
DIAG <- "plots/diagnostics"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)

ORD <- c("published_slope_baugh", "published_pc1_baugh",
         "delta_slope_baugh",     "delta_pc1_baugh",
         "slope_nnls",            "pc1_nnls",
         "delta_slope_nnls",      "delta_pc1_nnls")
PRETTY <- c(published_slope_baugh = "published Slope (MIP)",
            published_pc1_baugh   = "published PC1 (MIP)",
            delta_slope_baugh     = "delta Slope (MIP)",
            delta_pc1_baugh       = "delta PC1 (MIP)",
            slope_nnls            = "log-ratio Slope (WGS)",
            pc1_nnls              = "log-ratio PC1 (WGS)",
            delta_slope_nnls      = "delta Slope (WGS)",
            delta_pc1_nnls        = "delta PC1 (WGS)")
## the Manhattan grid has little vertical room per row, so its strips get a
## two-line form rather than the one-line labels used elsewhere
SHORT <- sub(" \\(", "\n(", PRETTY)

A <- rbindlist(lapply(list.files(D, full.names = TRUE), function(f) {
  d <- fread(cmd = paste("gzcat", shQuote(f)))
  d[, .(trait = trait[1], chr, ps, af, beta, p = p_wald, lp = -log10(p_wald))]
}))
A <- A[!is.na(p)]
NM <- A[trait == ORD[1], .N]
BF <- -log10(0.05 / NM)
A[, trait := factor(trait, ORD)]

infl <- A[, .(lambda = median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1),
              max_lp = max(lp), n_bf = sum(lp > BF)), by = trait][order(trait)]
cat(sprintf("markers %d   Bonferroni %.2f\n\n", NM, BF))
print(infl)
cat("\nGENOMIC INFLATION IS HIGH IN EVERY SCAN (lambda 1.16 to 1.56 on 99\n",
    "strains), so these thresholds are anti-conservative and the locus counts\n",
    "below should be read as a ranking rather than as calibrated significance.\n")

## --- loci ------------------------------------------------------------------
loci <- function(d, thr = BF, gap = 1e6) {
  x <- d[lp > thr][order(chr, ps)]
  if (!nrow(x)) return(NULL)
  x[, grp := cumsum(c(1, (diff(ps) > gap) | (head(chr, -1) != tail(chr, -1))))]
  x[, .(chr = chr[1], peak = ps[which.max(lp)], lp = max(lp), af = af[which.max(lp)],
        lo = min(ps), hi = max(ps), n = .N), by = grp][, grp := NULL][]
}
LT <- rbindlist(lapply(ORD, function(t) { L <- loci(A[trait == t])
  if (is.null(L)) NULL else cbind(trait = t, L) }))
fwrite(LT, file.path(DIAG, "baugh_reanalysis_locus_table.tsv"), sep = "\t")

## --- anchors: the two published peaks, and the delta-slope chrIV peak -------
anch <- data.table(
  label = c("published Slope peak\nV:15.92 Mb", "published PC1 peak\nV:15.93 Mb",
            "delta-Slope peak\nIV:16.22 Mb"),
  chr = c("V", "V", "IV"), ps = c(15917359, 15933722, 16218716))
AN <- A[anch, on = .(chr, ps)][, .(trait, label, lp)]

## --- genome-wide sharing, max per 100 kb bin -------------------------------
A[, bin := paste0(chr, "_", ps %/% 1e5)]
Bn <- A[, .(lp = max(lp)), by = .(trait, bin)]
M  <- dcast(Bn, bin ~ trait, value.var = "lp")
CM <- cor(as.matrix(M[, ..ORD]), method = "spearman", use = "pairwise.complete.obs")

theme_set(theme_bw(8.5) + theme(
  plot.title = element_text(face = "bold", size = 9.5),
  plot.subtitle = element_text(size = 7.3, colour = "grey30"),
  panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey93"),
  strip.text = element_text(size = 7.2, face = "bold")))

CHR <- c("I", "II", "III", "IV", "V", "X")
A[, chrf := factor(chr, CHR)]
pA <- ggplot(A[lp > 1], aes(ps / 1e6, lp)) +
  geom_point(aes(colour = lp > BF), size = 0.22, alpha = 0.55) +
  geom_hline(yintercept = BF, linetype = 2, colour = "grey45", linewidth = 0.3) +
  geom_vline(data = data.table(chrf = factor("V", CHR), x = 15.93),
             aes(xintercept = x), colour = "#C4302B", linewidth = 0.3, alpha = 0.7) +
  geom_vline(data = data.table(chrf = factor("IV", CHR), x = 16.22),
             aes(xintercept = x), colour = "#1A7F5A", linewidth = 0.3, alpha = 0.7) +
  facet_grid(trait ~ chrf, scales = "free_x", space = "free_x",
             labeller = labeller(trait = SHORT)) +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#2E4057"), guide = "none") +
  labs(title = "A  All eight scans, with the two anchor loci marked",
       subtitle = "Red line V:15.93 Mb (where the published traits peak); green line IV:16.22 Mb (where the delta-Slope traits peak). Dashed = Bonferroni.",
       x = "Position (Mb)", y = "-log10 p") +
  theme(panel.spacing.x = unit(1.5, "pt"), axis.text.x = element_text(size = 5.5),
        strip.text.y = element_text(size = 6, angle = 0, lineheight = 0.95))

CD <- as.data.table(as.table(CM)); setnames(CD, c("a", "b", "r"))
CD[, `:=`(a = factor(a, ORD), b = factor(b, rev(ORD)))]
pB <- ggplot(CD, aes(a, b, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.3,
            colour = ifelse(CD$r > 0.7, "white", "grey15")) +
  scale_fill_gradient(low = "#F2F2F2", high = "#2E4057", limits = c(0, 1), guide = "none") +
  scale_x_discrete(labels = PRETTY) + scale_y_discrete(labels = PRETTY) +
  labs(title = "B  Genome-wide signal sharing",
       subtitle = "Spearman of -log10 p, strongest marker per 100 kb bin",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6.4),
        axis.text.y = element_text(size = 6.4), panel.grid = element_blank())

AN[, trait := factor(trait, rev(ORD))]
pC <- ggplot(AN, aes(lp, trait)) +
  geom_vline(xintercept = BF, linetype = 2, colour = "grey45", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = lp, yend = trait), colour = "grey80", linewidth = 0.35) +
  geom_point(aes(colour = lp > BF), size = 2) +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#2E4057"), guide = "none") +
  scale_y_discrete(labels = PRETTY) +
  facet_wrap(~ label, nrow = 1) +
  labs(title = "C  Evidence at each anchor marker, every trait",
       subtitle = "Dashed = Bonferroni. A trait that recovers a locus clears it here.",
       x = "-log10 p", y = NULL) +
  theme(axis.text.y = element_text(size = 6.4))

fig <- pA / (pB | pC) + plot_layout(heights = c(1.35, 1))
ggsave(file.path(DIAG, "DIAG_baugh_reanalysis_mappings.pdf"), fig, width = 11.5, height = 10)
ggsave(file.path(DIAG, "DIAG_baugh_reanalysis_mappings.png"), fig, width = 11.5, height = 10, dpi = 200)
cat(sprintf("\nwrote %s/DIAG_baugh_reanalysis_mappings.{pdf,png}\n", DIAG))
