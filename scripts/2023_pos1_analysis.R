## 2023_pos1_analysis.R ------------------------------------------------------
## The 2023 pos-1 bulk competition experiment, re-plotted from the updated
## dataset. Two figures:
##
##   1  replicate reproducibility of the pos-1 delta-frequency phenotype
##      -> plots/SUPP_FIG_XX_original_pos1_dfreq_rep_correlation.{pdf,png}
##   2  the vst phenotype distribution and its association scan
##      -> plots/SUPP_FIG_XX_2023_pos1_vst_and_mapping.{pdf,png}
##
##   Rscript scripts/2023_pos1_analysis.R
##
## NOTE: the full analysis for this experiment lives outside this repository, at
##
##       /UCLA/Projects/bulkGWAS/lipid_RNAi/2023_pos1
##
## What is used here is only the exported end product of that analysis:
##   data/pos1_original/updated_analysis/final_dataset.csv
##       sample-level frequencies, 366 strains, T2 only, ctrl replicates A and B
##       and pos-1 replicates 1-4, each at three read-depth cutoffs (3, 5, 10).
##   data/pos1_original/updated_analysis/association_traits.csv
##       the per-strain association traits (delta_ctrl, vst, log2fc).
##   data/pos1_original/updated_analysis/vst_ctrl_pos-1_T2_loco_results.csv.gz
##       GEMMA LOCO results for vst_ctrl_pos-1_T2.
##
## Depth cutoff
## ------------
## final_dataset.csv carries three read-depth cutoffs. Cutoff 5 is the one the
## association traits were built from -- taking the mean pos-1 delta_ctrl per
## strain at cutoff 5 reproduces association_traits.csv's delta_ctrl_pos-1_T2 to
## 2.5e-16, against r = 0.9998 at cutoff 3 and 0.9988 at cutoff 10 -- so cutoff 5
## is used throughout, and the replicate figure is on the same footing as the map.
##
## The updated dataset has FOUR pos-1 replicates where the earlier version had
## two, so the replicate figure now shows all six pairwise comparisons rather
## than a single scatter.
##
## Duplicated strain entry -- JU1793
## ---------------------------------
## JU1793, and only JU1793, appears TWICE per sample in final_dataset.csv, at
## every depth cutoff and in both arms, with identical `sample` and
## `sample_info`. One row carries the real frequency and the other is ~1e-21,
## which is NNLS splitting one strain's abundance across two identical entries in
## the genotype reference rather than two distinct measurements.
##
## This matters because JU1793 is a cross parent. The shipped
## association_traits.csv AVERAGES the two rows, which halves the strain: its
## `ctrl_frq` of 0.006240 is the mean of four control rows, two of which are the
## spurious ~0, and its delta_ctrl_pos-1_T2 of 0.008570 follows from that.
##
## This script instead SUMS rows sharing a (strain, sample) key, recovering the
## strain's total assigned frequency, and recomputes the control baseline and
## delta from the collapsed table. For every other strain the collapse is a
## no-op and the recomputed delta reproduces the shipped value, which the script
## checks and reports. For JU1793 it roughly doubles the pos-1 delta, from
## 0.00857 to about 0.0171. The source files are left untouched; the correction
## is auditable here.
##
## The vst phenotype and the association scan are used AS SHIPPED -- the vst
## transform is not reproducible from this repository -- so JU1793's trait value
## in panel A and in the mapping still carries the averaging described above.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(ggtext)
})

DIR   <- "data/pos1_original/updated_analysis"
OUT   <- "plots"
DEPTH <- 5

COL_PT   <- "#2E4057"
COL_FIT  <- "#C4302B"
COL_PEAK <- "#C4302B"

CHROMS  <- c("I", "II", "III", "IV", "V", "X")
ALL_LEN <- c(I = 15072434, II = 15279421, III = 13783801,
             IV = 17493829, V = 20924180, X = 17718942)
fct_chr <- function(x) factor(as.character(x), levels = CHROMS)

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## ===========================================================================
## Figure 1 -- replicate reproducibility
## ===========================================================================
msg("replicate reproducibility")
dat <- read_csv(file.path(DIR, "final_dataset.csv"), show_col_types = FALSE)

## --- collapse rows that share a (strain, sample) key ------------------------
dup_keys <- dat %>% filter(depth_cutoff == DEPTH) %>%
  count(strain, sample) %>% filter(n > 1)
if (nrow(dup_keys)) {
  msg("  duplicated (strain, sample) keys: ", nrow(dup_keys),
      " rows across strain(s): ", paste(unique(dup_keys$strain), collapse = ", "))
}

collapsed <- dat %>% filter(depth_cutoff == DEPTH) %>%
  group_by(strain, sample, sample_info, time, rnai, replicate) %>%
  summarise(frq = sum(frq, na.rm = TRUE), n_rows = n(), .groups = "drop")

## control baseline per strain: mean over control replicates of the summed frq
base <- collapsed %>% filter(rnai == "ctrl") %>%
  group_by(strain) %>% summarise(ctrl_base = mean(frq, na.rm = TRUE), .groups = "drop")

reps <- collapsed %>% filter(rnai == "pos-1") %>%
  left_join(base, by = "strain") %>%
  transmute(strain, replicate = paste0("rep", replicate),
            delta_ctrl = frq - ctrl_base)

## does the recomputed delta reproduce the shipped trait for unduplicated strains?
shipped <- read_csv(file.path(DIR, "association_traits.csv"), show_col_types = FALSE) %>%
  transmute(strain, shipped = `delta_ctrl_pos-1_T2`)
chk <- reps %>% group_by(strain) %>%
  summarise(recomputed = mean(delta_ctrl, na.rm = TRUE), .groups = "drop") %>%
  inner_join(shipped, by = "strain") %>% filter(!is.na(shipped)) %>%
  mutate(dup = strain %in% dup_keys$strain)
cat("\n== recomputed vs shipped delta_ctrl_pos-1_T2 ==\n")
print(as.data.frame(chk %>% group_by(`duplicated strain` = dup) %>%
  summarise(n = n(), r = round(cor(recomputed, shipped), 6),
            `max abs diff` = signif(max(abs(recomputed - shipped)), 3),
            .groups = "drop")), row.names = FALSE)
if (any(chk$dup)) {
  cat("\n   the duplicated strain(s):\n")
  print(as.data.frame(chk %>% filter(dup) %>%
    transmute(strain, recomputed = round(recomputed, 5),
              shipped = round(shipped, 5),
              ratio = round(recomputed / shipped, 2))), row.names = FALSE)
}

rep_names <- sort(unique(reps$replicate))
wide <- reps %>% pivot_wider(names_from = replicate, values_from = delta_ctrl)

## every pair of replicates, long, so one facet per pair
pairs_tbl <- combn(rep_names, 2, simplify = FALSE) %>%
  map_dfr(function(p) {
    wide %>%
      transmute(strain, x = .data[[p[1]]], y = .data[[p[2]]],
                xr = p[1], yr = p[2]) %>%
      filter(!is.na(x), !is.na(y)) %>%
      mutate(pair = paste0(p[1], "  vs  ", p[2]))
  })

pair_stats <- pairs_tbl %>% group_by(pair) %>%
  summarise(rho = cor(x, y, method = "spearman"),
            n = n(), .groups = "drop") %>%
  mutate(lab = sprintf("rho = %.2f   n = %d", rho, n))

msg("  ", length(rep_names), " pos-1 replicates -> ", nrow(pair_stats), " pairs")

lim <- range(c(pairs_tbl$x, pairs_tbl$y), na.rm = TRUE)

p_rep <- ggplot(pairs_tbl, aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3,
              linetype = "dashed", colour = "grey65") +
  geom_point(size = 1.1, alpha = 0.45, colour = COL_PT) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = COL_FIT, linewidth = 0.5) +
  geom_text(data = pair_stats, aes(x = -Inf, y = Inf, label = lab),
            inherit.aes = FALSE, hjust = -0.07, vjust = 1.4, size = 3.1,
            colour = "grey20") +
  facet_wrap(~pair, nrow = 2) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "Δ frequency vs control", y = "Δ frequency vs control",
       title = "*pos-1* RNAi replicate reproducibility, 2023 pooled competition",
       subtitle = sprintf(paste0("Per-strain change in pool frequency relative to the ",
                                 "control, timepoint 2, read-depth cutoff %d. Every pair ",
                                 "of the four *pos-1* replicates; dashed line is y = x, ",
                                 "red is the fitted slope. Rows sharing a (strain, sample) ",
                                 "key are summed before the delta is taken, which affects ",
                                 "JU1793 only \u2014 see the script header."), DEPTH)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank(),
        plot.title = element_markdown(face = "bold", size = 12),
        plot.subtitle = element_markdown(size = 8.5, colour = "grey30"),
        plot.title.position = "plot")

for (ext in c("pdf", "png")) {
  f <- file.path(OUT, paste0("SUPP_FIG_XX_original_pos1_dfreq_rep_correlation.", ext))
  if (ext == "pdf") ggsave(f, p_rep, width = 10, height = 7, device = cairo_pdf)
  else              ggsave(f, p_rep, width = 10, height = 7, dpi = 300, bg = "white")
}
msg("  wrote SUPP_FIG_XX_original_pos1_dfreq_rep_correlation.{pdf,png}")

cat("\n== replicate pair correlations (Spearman) ==\n")
print(as.data.frame(pair_stats %>% transmute(pair, n, rho = round(rho, 3)) %>%
                      arrange(desc(rho))), row.names = FALSE)

## ===========================================================================
## Figure 2 -- vst phenotype distribution and the association scan
## ===========================================================================
msg("vst phenotype and association scan")
tr <- read_csv(file.path(DIR, "association_traits.csv"), show_col_types = FALSE) %>%
  transmute(strain, vst = `vst_ctrl_pos-1_T2`, delta = `delta_ctrl_pos-1_T2`)
## panel A and the scan use the shipped vst as-is; see the JU1793 note in the header

n_panel <- nrow(tr)
n_vst   <- sum(!is.na(tr$vst))
msg("  ", n_vst, " strains with a vst value, of ", n_panel, " in the trait file")

## (A) the phenotype distribution
p_dist <- tr %>% filter(!is.na(vst)) %>%
  ggplot(aes(vst)) +
  geom_histogram(bins = 32, fill = COL_PT, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey55") +
  labs(x = "*pos-1* response (vst)", y = "wild isotypes",
       title = "Phenotype distribution",
       subtitle = sprintf("n = %d of %d strains in the trait file", n_vst, n_panel)) +
  theme_classic(base_size = 11) +
  theme(axis.title.x = element_markdown(),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        plot.title.position = "plot")

## (B) the association scan
gw <- fread(file.path(DIR, "vst_ctrl_pos-1_T2_loco_results.csv.gz"),
            select = c("chr", "ps", "af", "beta", "p_wald"))
gw <- gw[chr %in% CHROMS]
gw[, neglog10p := -log10(p_wald)]

## Two thresholds. Bonferroni over every tested marker treats the markers as
## independent, which they are not -- this panel carries long haplotype blocks,
## so that line is far too strict. The eigen threshold divides alpha by the
## effective number of independent tests from the eigenvalues of the marker
## correlation matrix (Li & Ji 2005); see scripts/eigen_independent_tests.R.
source("scripts/gwas_thresholds.R")
TH <- gwas_thresholds("pos1_2023")
BF  <- TH$bonferroni
EIG <- TH$eigen
COL_THR <- "grey50"; COL_EIG <- "#1A7F5A"
msg("  ", nrow(gw), " markers | Bonferroni ", round(BF, 2),
    " | eigen ", round(EIG, 2), " (M_eff = ", round(TH$m_eff), ")")

## peak-preserving thinning: 464k markers is ~140 vertices per pixel
gwt <- gw %>% as_tibble() %>%
  mutate(chrom = fct_chr(chr), bin = floor(ps / 1e4)) %>%
  group_by(chrom, bin) %>% slice_max(neglog10p, n = 1, with_ties = FALSE) %>%
  ungroup()

top <- gw[order(-neglog10p)][, .SD[1], by = chr][order(-neglog10p)]
top[, `:=`(sig_bonferroni = neglog10p > BF, sig_eigen = neglog10p > EIG)]
cat("\n== per-chromosome peaks, vst_ctrl_pos-1_T2 ==\n")
print(as.data.frame(top[, .(chr, ps, af, beta = round(beta, 4),
                            neglog10p = round(neglog10p, 2),
                            sig_bonferroni, sig_eigen)]), row.names = FALSE)
cat("\nmarkers clearing Bonferroni: ", sum(gw$neglog10p > BF),
    " | clearing eigen: ", sum(gw$neglog10p > EIG), "\n", sep = "")

sig <- gwt %>% filter(neglog10p > BF)
sig_eig <- gwt %>% filter(neglog10p > EIG, neglog10p <= BF)
## one label per threshold, in the chromosome-I facet only
thr_lab <- tibble(
  chrom  = fct_chr(rep("I", 2)),
  pos.mb = 0.3,
  y      = c(BF, EIG),
  label  = c(sprintf("Bonferroni %.2f", BF),
             sprintf("eigen %.2f (M<sub>eff</sub> = %s)", EIG,
                     format(round(TH$m_eff), big.mark = ","))),
  col    = c(COL_THR, COL_EIG))

span <- bind_rows(
  tibble(chrom = fct_chr(names(ALL_LEN)), pos.mb = unname(ALL_LEN) / 1e6),
  tibble(chrom = fct_chr(names(ALL_LEN)), pos.mb = 0))

p_man <- ggplot(gwt, aes(ps / 1e6, neglog10p)) +
  geom_blank(data = span, aes(x = pos.mb, y = 0), inherit.aes = FALSE) +
  geom_hline(yintercept = BF, linetype = "dashed", linewidth = 0.3, colour = COL_THR) +
  geom_hline(yintercept = EIG, linetype = "dotted", linewidth = 0.45, colour = COL_EIG) +
  geom_point(size = 0.4, alpha = 0.5, colour = COL_PT) +
  { if (nrow(sig_eig)) geom_point(data = sig_eig, colour = COL_EIG, size = 0.8) } +
  { if (nrow(sig)) geom_point(data = sig, colour = COL_PEAK, size = 1.2) } +
  geom_richtext(data = thr_lab, aes(x = pos.mb, y = y, label = label),
                inherit.aes = FALSE, colour = thr_lab$col, size = 2.7,
                hjust = 0, vjust = -0.35, fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  facet_grid(. ~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, 25, 5), expand = expansion(mult = 0.02),
                     guide = guide_axis(check.overlap = TRUE)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = "Genomic Position (Mb)", y = "−log<sub>10</sub>*p*",
       title = "Association scan",
       subtitle = paste(strwrap(sprintf(paste0(
                  "GEMMA LOCO on vst_ctrl_pos-1_T2, %s markers, n = %d strains. Grey ",
                  "dashed is Bonferroni over every marker, which assumes they are ",
                  "independent; green dotted divides alpha by the %s effective ",
                  "independent tests estimated from the eigenvalues of the marker ",
                  "correlation matrix. %s markers clear Bonferroni, %s clear the eigen ",
                  "threshold."),
                  format(nrow(gw), big.mark = ","), TH$n_strain,
                  format(round(TH$m_eff), big.mark = ","),
                  format(sum(gw$neglog10p > BF), big.mark = ","),
                  format(sum(gw$neglog10p > EIG), big.mark = ",")),
                  width = 150), collapse = "<br>")) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        panel.spacing.x = grid::unit(4, "pt"),
        axis.title.y = element_markdown(),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.text.x = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_markdown(size = 8.5, colour = "grey30"),
        plot.title.position = "plot")

fig2 <- p_dist / p_man +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 14)))

ggsave(file.path(OUT, "SUPP_FIG_XX_2023_pos1_vst_and_mapping.pdf"), fig2,
       width = 11, height = 7, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_2023_pos1_vst_and_mapping.png"), fig2,
       width = 11, height = 7, dpi = 300, bg = "white")
msg("  wrote SUPP_FIG_XX_2023_pos1_vst_and_mapping.{pdf,png}")
