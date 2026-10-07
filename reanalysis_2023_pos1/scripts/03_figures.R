## 2023 pos-1: the figures that depend on this deconvolution, redone ---------
##
##   Rscript reanalysis_2023_pos1/scripts/03_figures.R
##     -> reanalysis_2023_pos1/plots/*.png|pdf
##
## Two of the manuscript's figures are built from the 2023 pos-1 deconvolution
## and so change with it:
##
##   Figure 1B   the phenotype distribution across the panel
##   Figure S4   replicate reproducibility, all six pairs
##
## Figure 1C, the association scan, cannot be redone here -- it needs the
## cluster. mapping_traits_dp5.csv is what it should be run on.
##
## ON THE SCALE. Figure 1B is drawn on the deposited vst trait, which cannot be
## rebuilt (see 02_build_traits.R), so the panel here is on the delta scale --
## the RNAi-minus-control change in pool frequency, which is what the
## deconvolution gives directly. The shape of the distribution is comparable;
## the axis is not.
##
## Each figure is drawn twice, old beside new, because the point of this
## directory is the comparison and not the new figure on its own.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(patchwork); library(ggtext)
})

DAT <- "reanalysis_2023_pos1/data"
OUT <- "reanalysis_2023_pos1/plots"
PH  <- "supplemental_data/phenotypes"
CUT <- 5L
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source("scripts/figure_palette.R")
source("scripts/figure_theme.R")
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

panel_title <- function(letter, txt = NULL) {
  lt <- paste0("<span style='font-size:13pt;color:#111111'>**", letter, "**</span>")
  if (is.null(txt)) lt else paste0(lt, " <span style='color:#555555'>", txt, "</span>")
}

new_tr <- fread(file.path(DAT, "pool_reference_traits_dp5.csv"))
old_tr <- fread(cmd = paste("gzcat", shQuote(file.path(PH, "pos1_2023_association_traits.csv.gz"))))
setnames(old_tr, c("delta_ctrl_pos-1_T2", "vst_ctrl_pos-1_T2"), c("delta", "vst"))

## ===========================================================================
## A -- the phenotype distribution, old reference against pool reference
## ===========================================================================
dist_panel <- function(x, lab, sub) {
  d <- tibble(v = x[is.finite(x)])
  ggplot(d, aes(v)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey60") +
    geom_histogram(bins = 40, fill = COL_PT, colour = NA) +
    labs(x = "pos-1 response (delta frequency vs control)", y = "Wild isotypes",
         title = panel_title(lab, sprintf("%s, n = %d", sub, nrow(d)))) +
    theme_pub() + theme(plot.title = element_markdown())
}
pA <- dist_panel(old_tr$delta, "A", "deposited, 367-strain reference")
pB <- dist_panel(new_tr$delta_ctrl_pos1_T2, "B", "pool reference, 224 strains")

## the two phenotypes against each other, for the strains in both
cmp <- merge(new_tr[, .(strain, new = delta_ctrl_pos1_T2)],
             old_tr[, .(strain, old = delta)], by = "strain")[!is.na(new) & !is.na(old)]
rho <- cor(cmp$old, cmp$new, method = "spearman")
pC <- ggplot(cmp, aes(old, new)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey60") +
  geom_point(size = 1.8, alpha = 0.6, colour = COL_PT) +
  labs(x = "deposited (367-strain reference)", y = "pool reference (224 strains)",
       title = panel_title("C", sprintf("Spearman %.3f, n = %d", rho, nrow(cmp)))) +
  theme_pub() + theme(plot.title = element_markdown())

fig1 <- (pA | pB) / pC + plot_layout(heights = c(1, 1.2))
ggsave(file.path(OUT, "phenotype_old_vs_pool_reference.png"), fig1,
       width = 10, height = 8, dpi = 200, bg = "white")
ggsave(file.path(OUT, "phenotype_old_vs_pool_reference.pdf"), fig1,
       width = 10, height = 8, device = cairo_pdf)
msg("wrote phenotype_old_vs_pool_reference | Spearman ", round(rho, 3), " over ", nrow(cmp))

## ===========================================================================
## B -- replicate reproducibility on the corrected frequencies
## ===========================================================================
f <- fread(cmd = paste("gzcat", shQuote(file.path(DAT, "pool_reference_frequencies.csv.gz"))))[depth_cutoff == CUT]
ctrl <- f[grepl("^T2_ctrl", sample_info), .(ctrl_frq = mean(frq)), by = strain]
reps <- f[grepl("^T2_pos-1", sample_info)][ctrl, on = "strain"][
  , .(strain, replicate = sub("^T2_pos-1_", "rep", sample_info), delta = frq - ctrl_frq)]
## a strain with no control signal has no defined response, as in the traits
reps <- reps[strain %in% new_tr[!is.na(delta_ctrl_pos1_T2)]$strain]
wide <- dcast(reps, strain ~ replicate, value.var = "delta")
rn <- sort(setdiff(names(wide), "strain"))
msg("replicate panel: ", nrow(wide), " strains x ", length(rn), " replicates")

pairs_tbl <- rbindlist(lapply(combn(rn, 2, simplify = FALSE), function(p)
  wide[, .(strain, x = get(p[1]), y = get(p[2]),
           pair = paste0(p[1], "  vs  ", p[2]))]))
labs <- pairs_tbl[, .(rho = cor(x, y, method = "spearman"), n = .N), by = pair]
labs[, lab := sprintf("rho = %.2f   n = %d", rho, n)]
cat("\n== replicate pair correlations, pool reference ==\n")
print(as.data.frame(labs[order(-rho), .(pair, n, rho = round(rho, 3))]), row.names = FALSE)

fig2 <- ggplot(pairs_tbl, aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey60") +
  geom_point(size = 1.4, alpha = 0.55, colour = COL_PT) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, colour = "#D55E00", formula = y ~ x) +
  geom_text(data = labs, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.08, vjust = 1.6, size = 3.1, colour = "grey20") +
  facet_wrap(~pair, ncol = 3) +
  labs(x = "delta frequency vs control", y = "delta frequency vs control") +
  theme_pub() + theme(strip.text = element_text(face = "bold", size = 10))

ggsave(file.path(OUT, "replicate_reproducibility_pool_reference.png"), fig2,
       width = 10, height = 7, dpi = 200, bg = "white")
ggsave(file.path(OUT, "replicate_reproducibility_pool_reference.pdf"), fig2,
       width = 10, height = 7, device = cairo_pdf)
msg("wrote replicate_reproducibility_pool_reference")
