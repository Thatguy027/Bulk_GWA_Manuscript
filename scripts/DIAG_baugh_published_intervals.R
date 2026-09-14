## The four published Baugh QTL intervals, one comparison per panel ----------
##
## Webster et al. report four significant intervals:
##   Slope  IV:15,939,340-16,613,710   and  V:15,660,911-17,615,557
##   PC1    V:1,345,848-2,764,788      and  V:15,775,895-18,065,050
##
## The two right-arm chromosome V intervals overlap almost entirely (15.78-17.62
## Mb is common to both), so they are drawn separately but should be read as one
## region seen through two traits rather than two independent findings.
##
## Each panel is the published trait against ONE reconstruction, so the question
## "does this trait follow the published one here" is answered without reading
## four overlaid lines. Slope comparators are drawn against published Slope and
## PC1 comparators against published PC1. Lines are the strongest marker per
## 10 kb, which keeps the shape without drawing 26,000 points.
##
## The four irld genes the source study nominates are marked: irld-39 on
## chromosome IV, and irld-11, irld-57 and irld-52 on chromosome V. Each
## published interval contains exactly one, and irld-39 lies 16.5 kb from the
## chromosome IV peak marker.
##
## Within-interval correlations are a separate figure,
## scripts/DIAG_baugh_interval_correlations.R.
##
## Reads from data/, so this is a diagnostic and outside the deposit-only build.
## Writes plots/diagnostics/DIAG_baugh_published_intervals.{pdf,png} and
## baugh_published_interval_table.tsv
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
CMP <- data.table(
  trait = c("delta_slope_baugh","slope_nnls","delta_slope_nnls",
            "delta_pc1_baugh","pc1_nnls","delta_pc1_nnls"),
  cmp   = factor(CMPL, levels = CMPL),
  ref   = rep(c("published_slope_baugh", "published_pc1_baugh"), each = 3))

## --- the candidate irld genes named in the source study ---------------------
## Coordinates from WS283. Each published interval contains one, and irld-39
## sits 16.5 kb from the chromosome IV peak marker IV:16,218,716.
IRLD <- data.table(
  gene = c("irld-39", "irld-11", "irld-57", "irld-52"),
  chr  = c("IV", "V", "V", "V"),
  mid  = c((16233148 + 16236555) / 2, (1655145 + 1657403) / 2,
           (15724810 + 15726710) / 2, (15779679 + 15781247) / 2))
## a gene is drawn in every interval panel whose drawn window contains it
IRLD <- IRLD[IV, on = .(chr), allow.cartesian = TRUE
  ][, pad := (hi - lo) * 0.25
  ][mid >= lo - pad & mid <= hi + pad
  ][, .(gene, interval = id, mid)
  ][order(interval, mid)
  ## two of them are 55 kb apart, so labels alternate height within a panel
  ][, ylab := c(14.3, 12.3)[seq_len(.N) %% 2 + 1], by = interval][]

A <- rbindlist(lapply(list.files(D, full.names = TRUE), function(f) {
  d <- fread(cmd = paste("gzcat", shQuote(f)))
  d[, .(trait = trait[1], chr, ps, af, p = p_wald, lp = -log10(p_wald))] }))
A <- A[!is.na(p)]
BF <- -log10(0.05 / A[trait == A$trait[1], .N])

## --- per-interval summary, kept for the table ------------------------------
TB <- rbindlist(lapply(seq_len(nrow(IV)), function(i) {
  A[chr == IV$chr[i] & ps >= IV$lo[i] & ps <= IV$hi[i],
    .(peak = ps[which.max(lp)], max_lp = max(lp), af = af[which.max(lp)],
      n_bf = sum(lp > BF), n_mark = .N), by = trait][, interval := IV$id[i]][] }))
fwrite(TB, file.path(DIAG, "baugh_published_interval_table.tsv"), sep = "\t")

## --- one row per (comparison, interval) ------------------------------------
Z <- rbindlist(lapply(seq_len(nrow(IV)), function(i) {
  pad <- (IV$hi[i] - IV$lo[i]) * 0.25
  R <- A[chr == IV$chr[i] & ps >= IV$lo[i] - pad & ps <= IV$hi[i] + pad,
         .(lp = max(lp)), by = .(trait, bin = (ps %/% 1e4) * 1e4)]
  rbindlist(lapply(seq_len(nrow(CMP)), function(k) {
    rbind(R[trait == CMP$ref[k]][, role := "published"],
          R[trait == CMP$trait[k]][, role := "reconstruction"]
    )[, `:=`(cmp = CMP$cmp[k], interval = IV$id[i])][] })) }))
Z[, role := factor(role, c("published", "reconstruction"))]
BND <- IV[, .(interval = id, lo, hi)][rep(1:4, each = length(CMPL))][
  , cmp := factor(rep(CMPL, 4), levels = CMPL)][]

theme_set(theme_bw(8.5) + theme(
  plot.title = element_text(face = "bold", size = 10),
  plot.subtitle = element_text(size = 7.5, colour = "grey30"),
  panel.grid.minor = element_blank(), legend.position = "top",
  legend.title = element_blank(), legend.key.height = unit(8, "pt"),
  strip.background = element_rect(fill = "grey93"),
  strip.text = element_text(size = 7.2, face = "bold")))

p <- ggplot(Z, aes(bin / 1e6, lp, colour = role, linewidth = role)) +
  geom_rect(data = BND, aes(xmin = lo / 1e6, xmax = hi / 1e6, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey70", alpha = 0.22) +
  geom_hline(yintercept = BF, linetype = 2, colour = "grey40", linewidth = 0.3) +
  geom_vline(data = IRLD, aes(xintercept = mid / 1e6), inherit.aes = FALSE,
             colour = "#1A7F5A", linewidth = 0.4, linetype = 5) +
  geom_text(data = copy(IRLD)[, cmp := factor(CMPL[1], levels = CMPL)],
            aes(x = mid / 1e6, y = ylab, label = gene), inherit.aes = FALSE,
            colour = "#1A7F5A", size = 2.3, hjust = -0.12, fontface = "italic") +
  geom_line(alpha = 0.92) +
  scale_colour_manual(values = c(published = "#1A1A1A", reconstruction = "#C4302B")) +
  scale_linewidth_manual(values = c(published = 0.32, reconstruction = 0.5), guide = "none") +
  facet_grid(cmp ~ interval, scales = "free_x") +
  labs(title = "The four published Baugh intervals, one comparison per panel",
       subtitle = paste("Black is the matching published trait, red the reconstruction.",
                        "Strongest marker per 10 kb; grey band is the published interval with 25% flanking;",
                        "dashed line Bonferroni. Green dashed lines are the candidate irld genes."),
       x = "Position (Mb)", y = "-log10 p")
ggsave(file.path(DIAG, "DIAG_baugh_published_intervals.pdf"), p, width = 11, height = 10)
ggsave(file.path(DIAG, "DIAG_baugh_published_intervals.png"), p, width = 11, height = 10, dpi = 200)
cat(sprintf("wrote %s/DIAG_baugh_published_intervals.{pdf,png}\n", DIAG))
