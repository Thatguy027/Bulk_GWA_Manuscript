## Figure 1 -- shared vocabulary and panel builders --------------------------
## Sourced by the Figure 1 variants. Not run on its own.
##
##   scripts/Figure1.R        the platform-validation figure (Baugh L1 data)
##   scripts/Figure1_pos1.R   validation panel + the 2023 pos-1 phenotype and map
##
## Every panel is a function returning a ggplot, and every one takes a `letter`
## and a `bare` flag, so panels can be recombined into any arrangement without
## editing them. `bare = TRUE` drops the descriptive title and subtitle and
## leaves the panel letter, which is what a main display wants; `bare = FALSE`
## keeps them, which is what a standalone or supplementary output wants.
##
## Replaces scripts/legacy/Figure1_legacy_{exploratory,clean}.R. Three things
## are different beyond the styling.
##
## 1  IT RUNS UNDER Rscript. The legacy scripts opened with
##      setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
##    which only resolves inside RStudio, so Figure 1 was the one figure in the
##    repository that could not be regenerated from the command line. These
##    scripts assume the working directory is the project root, as Figures 2-4
##    do, and take no action to change it.
##
## 2  THE NNLS IS CACHED. The legacy script recomputed the deconvolution on
##    every run -- crossprod on a 30 MB genotype matrix, then a per-sample
##    RcppML::nnls -- and then re-saved the cache it had just written. The
##    result was already committed as data/baugh/2024_processedBOOTs_with_MIP
##    .RData, so that file is now the default input and the expensive path runs
##    only on baugh_frequencies(refresh = TRUE).
##
## 3  THE PANEL LETTERS ARE REAL. The legacy script forced them with
##      plot_annotation(tag_levels = list(c("B")))
##    on a single plot, so the letter was a literal string that would not
##    follow a rearrangement. Here each panel carries its own letter.
##
## What the figure argues
## ----------------------
## That NNLS deconvolution of pooled whole-genome sequence recovers the same
## per-strain frequency dynamics as published targeted MIP-seq, so the pooled
## assay measures what it claims to. The test bed is the Baugh L1 starvation
## experiment, where both measurements exist for the same samples.
##
## Dropped from the legacy version: the PCA / PC1-vs-PC1 comparison (the
## exploratory script's Figure 1 was slope + PC1, and its PCA block appears
## twice, near-identically), the PC1-PC9 correlation panel, and a "delta rank
## vs sites shared with N2" diagnostic. None of it was in the last arrangement.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(ggtext)
})

OUT   <- "plots"
BAUGH <- "supplemental_data/deconvolution"
POS1  <- "supplemental_data/phenotypes"

CACHE <- file.path(BAUGH, "baugh_nnls_with_mipseq.RData")
BOOT  <- file.path(BAUGH, "2024bootstrapINPUT.Rdata")   # Dryad only
## The bootstrap array's first dimension is unnamed; its order is the column
## order of the genotype matrix. Rather than open that 31 MB matrix just to
## recover 102 strain names, the order is shipped as a one-column file.
SORDER <- file.path(BAUGH, "baugh_strain_order.txt")
MIP   <- file.path(BAUGH, "mipseq_frequencies.txt.gz")
DS    <- file.path(BAUGH, "baugh_downsampled_slopes.rda")
BPRED <- file.path(BAUGH, "baugh_bootstrap_array.rda")
BCACHE<- file.path(BAUGH, "cache_boot_slopes.rds")
FCACHE<- file.path(BAUGH, "cache_boot_freq.rds")

COL_PT   <- "#2E4057"   # the point colour used by the 2023 pos-1 figures
COL_FIT  <- "#C4302B"
COL_PEAK <- "#C4302B"
COL_THR  <- "grey50"
COL_EIG  <- "#1A7F5A"
COL_HIST <- "#69B3A2"   # the histogram fill the legacy figure used

CHROMS  <- c("I", "II", "III", "IV", "V", "X")
ALL_LEN <- c(I = 15072434, II = 15279421, III = 13783801,
             IV = 17493829, V = 20924180, X = 17718942)
fct_chr <- function(x) factor(as.character(x), levels = CHROMS)

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## panel letter, in the weight the other figures use
panel_title <- function(letter) {
  if (is.null(letter)) return(NULL)
  paste0("<span style='font-size:14pt;color:#111111'>**", letter, "**</span>")
}

## a title line with the panel letter in front of it, for bare = FALSE.
## The separator is two literal EM SPACE characters, not "&nbsp;": gridtext
## does not decode that entity, so the letter ran into the title text.
titled <- function(letter, txt) {
  if (is.null(letter)) return(txt)
  paste0(panel_title(letter), "\u2003\u2003", txt)
}

theme_pub <- function(base_size = 11.5) {
  theme_classic(base_size = base_size) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = base_size),
          panel.spacing.x = grid::unit(8, "pt"),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_markdown(size = base_size),
          plot.subtitle = element_markdown(size = base_size - 3,
                                           colour = "grey30"),
          plot.title.position = "plot",
          legend.key.size = grid::unit(9, "pt"))
}

## wrap a long subtitle: patchwork will not do it, and an unwrapped one is
## silently truncated at the panel edge
wrap_md <- function(txt, width = 120) {
  paste(strwrap(txt, width = width), collapse = "<br>")
}

## a Spearman label, computed here rather than by ggpubr::stat_cor, so the
## number in the caption and the number on the panel come from one place
rho_label <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  p  <- ct$p.value
  sprintf("rho = %.2f, %s, n = %d", unname(ct$estimate),
          if (p < 1e-4) "*p* < 1e-4" else sprintf("*p* = %.3g", p), sum(ok))
}

## ===========================================================================
## the Baugh frequencies: NNLS on pooled WGS, joined to published MIP-seq
##
## refresh = FALSE reads the committed result. refresh = TRUE redoes the
## deconvolution from the 30 MB genotype matrix, which is the only slow step in
## the whole figure, and rewrites the cache.
## ===========================================================================
baugh_frequencies <- function(refresh = FALSE) {
  if (!refresh && file.exists(CACHE)) {
    e <- new.env(); load(CACHE, e)
    msg("  frequencies from cache: ", CACHE, " (", nrow(e$wgs_mip_results),
        " rows)")
    return(as_tibble(e$wgs_mip_results))
  }

  msg("  recomputing NNLS from ", BOOT, " -- this is the slow path")
  stopifnot(requireNamespace("RcppML", quietly = TRUE))
  e <- new.env(); load(BOOT, e); load(file.path(BAUGH,
        "2024baugh_bootstrap_prediction.rda"), e)

  ## markers with any missing genotype are dropped rather than zero-filled: a
  ## zero here is a real homozygous-reference call, so imputing one would
  ## change the design matrix rather than just lose a marker
  keep <- which(apply(e$flipped_bootstrap_input[[1]], 1,
                      function(x) sum(is.na(x))) == 0)
  gt   <- e$flipped_bootstrap_input[[1]][keep, ]
  a_ct <- e$flipped_bootstrap_input[[2]][keep, ]

  fullGy  <- crossprod(gt, a_ct)
  fullGGp <- crossprod(gt)
  pred <- apply(fullGy, 2,
                function(x) as.vector(RcppML::nnls(fullGGp, matrix(x),
                                                   fast_nnls = TRUE)))
  pred <- apply(pred, 2, function(x) x / sum(x))   # frequencies, not weights

  predictions_df <- data.frame(pred) %>%
    mutate(strain = sub("_.*$", "", colnames(fullGGp))) %>%
    pivot_longer(-strain, names_to = "sample", values_to = "frq")

  strains <- sub("_.*$", "", colnames(e$flipped_bootstrap_input[[1]]))
  boot_df <- data.frame(strain = strains,
                        bootmean = apply(e$ab, c(1, 2), mean),
                        bootse = e$bootse) %>%
    pivot_longer(-strain, names_to = "stat_sample", values_to = "value") %>%
    separate(stat_sample, into = c("stat", "sample"), sep = "\\.") %>%
    pivot_wider(names_from = stat, values_from = value)

  pr_boot_df <- boot_df %>%
    left_join(predictions_df, by = c("sample", "strain")) %>%
    separate(sample, into = c("replicate", "day"), sep = "_",
             remove = FALSE, extra = "drop") %>%
    mutate(baseline = grepl("baseline", sample),
           day = as.numeric(gsub("d", "", day))) %>%
    select(sample, replicate, day, baseline, strain, bootmean, bootse, frq)

  elife_fq <- fread(MIP) %>% rename(strain = V1) %>%
    pivot_longer(-strain, names_to = "sample", values_to = "published_frq") %>%
    mutate(sample = gsub("_BL", "_d1_baseline", sample))

  wgs_mip_results <- left_join(pr_boot_df, elife_fq,
                               by = c("sample", "strain"))
  save(wgs_mip_results, file = CACHE)
  msg("  wrote ", CACHE)
  as_tibble(wgs_mip_results)
}

## ---------------------------------------------------------------------------
## per-strain, per-replicate slope of delta frequency against day, for each
## platform. Baseline samples and day 17 are excluded and deltas are taken
## against day 1, as in the published analysis.
## ---------------------------------------------------------------------------
## the slope of a simple regression in closed form. Identical to
## coef(lm(y ~ x))[2] but vectorisable and free of the cur_data() deprecation,
## which matters because the bootstrap path fits this 40,000 times.
ols_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2) return(NA_real_)
  x <- x[ok]; y <- y[ok]
  xc <- x - mean(x)
  d  <- sum(xc^2)
  if (d == 0) return(NA_real_)
  sum(xc * y) / d
}

platform_slopes <- function(freq) {
  base <- freq %>% filter(!baseline, day == 1) %>%
    select(replicate, strain, base_frq = frq, base_pubfrq = published_frq)

  delta <- freq %>% filter(!baseline, day != 17) %>%
    select(sample:strain, frq, published_frq) %>%
    left_join(base, by = c("replicate", "strain")) %>%
    mutate(delta_d1_wgs = frq - base_frq,
           delta_d1_mip = published_frq - base_pubfrq)

  fit <- function(d, col) {
    d %>% filter(is.finite(.data[[col]])) %>%
      group_by(replicate, strain) %>%
      summarise(slope = ols_slope(day, .data[[col]]), .groups = "drop")
  }

  wgs <- fit(delta, "delta_d1_wgs") %>% rename(wgs_slope = slope)
  mip <- fit(delta, "delta_d1_mip") %>% rename(mip_slope = slope)
  full_join(wgs, mip, by = c("replicate", "strain"))
}

## ===========================================================================
## panels
## ===========================================================================

## --- the platform comparison: MIP-seq slope against NNLS slope -------------
## N2 is excluded throughout: it is the reference strain, its genotype shares
## sites with every other strain in the pool, and its NNLS frequency is not
## identified on the same footing as the wild isolates'.
panel_slope <- function(slopes, letter = "A", bare = TRUE,
                        base_size = 11.5, boot_mat = NULL, level = 0.95) {
  d <- slopes %>% filter(strain != "N2") %>%
    group_by(strain) %>%
    summarise(mip = mean(mip_slope, na.rm = TRUE),
              wgs = mean(wgs_slope, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(mip), is.finite(wgs))

  lab <- rho_label(d$mip, d$wgs)
  bars <- NULL

  if (!is.null(boot_mat)) {
    bs <- boot_slope_summary(boot_mat, level = level)
    d  <- d %>% left_join(bs, by = "strain")
    br <- boot_rho(boot_mat, setNames(d$mip, d$strain), level = level)
    msg("  panel slope: bootstrap ", ncol(boot_mat), " replicates | rho ",
        sprintf("%.3f [%.3f, %.3f]", br$med, br$lo, br$hi),
        " over ", br$n, " strains | median interval width ",
        signif(median(d$hi - d$lo, na.rm = TRUE), 3))
    ## the interval on rho replaces the p value: a p value against rho = 0 is
    ## not the question here, how reproducible rho is under resampling is
    ## three decimals: the interval is narrow enough that two would print
    ## "0.97 [0.97, 0.97]" and say nothing
    lab <- sprintf("rho = %.3f [%.3f, %.3f], n = %d", br$med, br$lo, br$hi,
                   nrow(d))
    ## bars under the points, no caps: 98 strains cluster near zero and caps
    ## turn the cluster into a solid block
    bars <- geom_errorbar(aes(ymin = lo, ymax = hi), width = 0,
                          linewidth = 0.35, colour = COL_PT, alpha = 0.4)
  }

  msg("  panel slope: ", nrow(d), " strains | ", gsub("[*]", "", lab))

  ggplot(d, aes(mip, wgs)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.4, colour = COL_FIT) +
    bars +
    geom_point(size = 1.9, alpha = 0.65, colour = COL_PT) +
    geom_richtext(data = tibble(lab = lab), aes(x = -Inf, y = Inf, label = lab),
                  inherit.aes = FALSE, hjust = -0.06, vjust = 1.5,
                  size = 3.1, colour = COL_FIT, fill = NA, label.color = NA,
                  label.padding = grid::unit(rep(0, 4), "pt")) +
    coord_equal() +
    labs(x = "MIP-seq slope", y = "NNLS slope",
         title = if (bare) panel_title(letter)
                 else titled(letter, "**Platform agreement**"),
         subtitle = if (bare) NULL else wrap_md(paste0(
           "Per-strain rate of change in pool frequency, Baugh L1 starvation ",
           "time course, one point per wild isolate averaged over replicates. ",
           "Dashed line is y = x, not a fit. N2 excluded.",
           if (is.null(boot_mat)) "" else sprintf(paste0(
             " Vertical bars are %g%% percentile intervals over %d bootstrap ",
             "replicates of the deconvolution, with the slope recomputed ",
             "inside each."), level * 100, ncol(boot_mat)))))+
    theme_pub(base_size)
}

## --- the distribution of per-sample agreement -----------------------------
## The slope panel averages over samples, so it cannot show a sample that
## disagrees. This does: one Spearman correlation per sample, over strains.
panel_sample_rho <- function(freq, letter = "B", bare = TRUE,
                             base_size = 11.5) {
  d <- freq %>% group_by(sample) %>%
    summarise(rho = cor(frq, published_frq, method = "spearman",
                        use = "complete.obs"), .groups = "drop")
  med <- median(d$rho)
  msg("  panel sample rho: ", nrow(d), " samples | median ", round(med, 3),
      " | min ", round(min(d$rho), 3))

  ggplot(d, aes(rho)) +
    geom_histogram(binwidth = 0.02, fill = COL_HIST, colour = NA) +
    geom_vline(xintercept = med, linetype = "dashed", linewidth = 0.4,
               colour = COL_FIT) +
    labs(x = "Spearman's &rho;, per sample", y = "Samples",
         title = if (bare) panel_title(letter)
                 else titled(letter, "**Per-sample agreement**"),
         subtitle = if (bare) NULL else wrap_md(sprintf(paste0(
           "One correlation per sample across strains, %d samples. Dashed ",
           "line is the median, %.2f."), nrow(d), med))) +
    theme_pub(base_size) +
    theme(axis.title.x = element_markdown())
}

## --- how much sequencing depth the deconvolution needs --------------------
## Read counts are multinomially subsampled to fixed depths and the slopes
## recomputed; each point is the correlation of those slopes with the MIP
## slopes. See scripts/baugh_L1_DownSample_Counts.R for the subsampling.
panel_downsample <- function(slopes, letter = "C", bare = TRUE,
                             base_size = 11.5, drop_depth = 0.5) {
  e <- new.env(); load(DS, e)
  ds <- e$ds_predictions_df %>% filter(strain != "N2") %>%
    separate(sample, into = c("replicate", "day"), sep = "_",
             remove = FALSE, extra = "drop") %>%
    mutate(baseline = grepl("baseline", sample),
           day = as.numeric(gsub("d", "", day)))

  base <- ds %>% filter(!baseline, day == 1) %>%
    select(replicate, strain, ds_n, base_frq = ds_frq)

  ds_slope <- ds %>% filter(!baseline, day != 17) %>%
    left_join(base, by = c("replicate", "strain", "ds_n")) %>%
    mutate(delta = ds_frq - base_frq) %>%
    filter(is.finite(delta)) %>%
    group_by(replicate, strain, ds_n) %>%
    summarise(slope = ols_slope(day, delta), .groups = "drop") %>%
    group_by(strain, ds_n) %>%
    summarise(slope = mean(slope, na.rm = TRUE), .groups = "drop")

  mip <- slopes %>% filter(strain != "N2") %>% group_by(strain) %>%
    summarise(mip = mean(mip_slope, na.rm = TRUE), .groups = "drop")

  d <- ds_slope %>% left_join(mip, by = "strain") %>%
    filter(is.finite(slope), is.finite(mip)) %>%
    group_by(ds_n) %>%
    summarise(rho = cor(slope, mip, method = "spearman"), n = n(),
              .groups = "drop")

  ## the legacy figure dropped the 0.5x depth silently. It is dropped here too,
  ## but named, and the value is reported so the omission is auditable.
  if (!is.null(drop_depth)) {
    gone <- d %>% filter(ds_n %in% drop_depth)
    if (nrow(gone))
      msg("  panel downsample: dropping depth(s) ",
          paste(sprintf("%gx (rho %.3f)", gone$ds_n, gone$rho),
                collapse = ", "))
    d <- d %>% filter(!ds_n %in% drop_depth)
  }
  msg("  panel downsample: ",
      paste(sprintf("%gx %.3f", d$ds_n, d$rho), collapse = " | "))

  ## a real continuous axis on a log scale, with breaks read from the data.
  ## The legacy version used factor(ds_n) plus a hardcoded
  ## scale_x_discrete(labels = c("0.25","1","3","5","10")), which was a no-op
  ## relabel that would have silently mismatched had the depth set changed.
  ggplot(d, aes(ds_n, rho)) +
    geom_line(linewidth = 0.4, colour = "grey60") +
    geom_point(size = 2.6, colour = COL_PT) +
    scale_x_log10(breaks = d$ds_n, labels = function(x) paste0(x, "×")) +
    labs(x = "Sequencing depth", y = "Spearman's &rho; vs MIP-seq",
         title = if (bare) panel_title(letter)
                 else titled(letter, "**Depth requirement**"),
         subtitle = if (bare) NULL else wrap_md(paste0(
           "Reads subsampled to each depth and the slopes recomputed, then ",
           "correlated with the MIP-seq slopes. Agreement saturates by 3x."))) +
    theme_pub(base_size) +
    theme(axis.title.y = element_markdown(),
          panel.grid.minor.x = element_blank())
}

## ===========================================================================
## the 2023 pos-1 panels
##
## These are the two panels of SUPP_FIG_XX_2023_pos1_vst_and_mapping, built
## here so they can be embedded in Figure 1. The trait and the scan are used AS
## SHIPPED; see scripts/2023_pos1_analysis.R for the depth-cutoff choice, the
## JU1793 duplicate-entry correction and the caveat that the shipped vst still
## carries the averaged JU1793 value.
##
## NOTE: scripts/2023_pos1_analysis.R still builds its own copies of these two
## panels. The two are intentionally not yet merged -- that script's output is
## already approved -- so a change to the look of these panels has to be made
## in both places until they are converged.
## ===========================================================================
pos1_traits <- function() {
  read_csv(file.path(POS1, "pos1_2023_association_traits.csv.gz"), show_col_types = FALSE) %>%
    transmute(strain, vst = `vst_ctrl_pos-1_T2`)
}

## --- the pos-1 phenotype distribution over the mapping panel --------------
panel_pos1_dist <- function(letter = "B", bare = TRUE, base_size = 11.5,
                            mark = NULL) {
  tr <- pos1_traits()
  n_panel <- nrow(tr); n_vst <- sum(!is.na(tr$vst))
  msg("  panel pos-1 distribution: ", n_vst, " of ", n_panel,
      " strains have a vst value")

  d <- tr %>% filter(!is.na(vst))
  p <- ggplot(d, aes(vst)) +
    geom_histogram(bins = 32, fill = COL_PT, colour = "white",
                   linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
               colour = "grey55")

  ## optional: name particular strains in the distribution, for tying the
  ## mapping panel to strains used downstream
  if (!is.null(mark)) {
    mk <- d %>% filter(strain %in% mark)
    p <- p +
      geom_vline(data = mk, aes(xintercept = vst), linewidth = 0.5,
                 colour = COL_FIT) +
      ## horizontal, inside the panel. An angled label at y = Inf lands
      ## outside the panel and is clipped away without warning.
      geom_richtext(data = mk, aes(x = vst, y = Inf, label = strain),
                    inherit.aes = FALSE, colour = COL_FIT, size = 2.7,
                    hjust = -0.08, vjust = 1.6, fill = "white",
                    label.color = NA,
                    label.padding = grid::unit(c(1, 2, 1, 2), "pt"))
    msg("    marked: ", paste(mk$strain, collapse = ", "))
  }

  p +
    labs(x = "*pos-1* response (VST)", y = "Wild isotypes",
         title = if (bare) panel_title(letter)
                 else titled(letter, "**Phenotype distribution**"),
         subtitle = if (bare) NULL else
           sprintf("n = %d of %d strains in the trait file", n_vst, n_panel)) +
    theme_pub(base_size) +
    theme(axis.title.x = element_markdown())
}

## --- the pos-1 association scan ------------------------------------------
panel_pos1_manhattan <- function(letter = "C", bare = TRUE, base_size = 11.5) {
  gw <- fread("supplemental_data/mapping/pos1_2023_gemma_loco.csv.gz",
              select = c("chr", "ps", "af", "beta", "p_wald"))
  gw <- gw[chr %in% CHROMS]
  gw[, neglog10p := -log10(p_wald)]

  ## Two thresholds. Bonferroni over every tested marker treats the markers as
  ## independent, which they are not -- this panel carries long haplotype
  ## blocks -- so that line is far too strict. The eigen threshold divides
  ## alpha by the effective number of independent tests from the eigenvalues of
  ## the marker correlation matrix (Li & Ji 2005).
  source("scripts/gwas_thresholds.R", local = TRUE)
  TH  <- gwas_thresholds("pos1_2023")
  BF  <- TH$bonferroni
  EIG <- TH$eigen
  msg("  panel pos-1 scan: ", nrow(gw), " markers | Bonferroni ",
      round(BF, 2), " | eigen ", round(EIG, 2),
      " (M_eff = ", round(TH$m_eff), ")")

  ## peak-preserving thinning: bin, keep the maximum. 464k markers is ~140
  ## vertices per pixel, and thinning by row order would move the peaks.
  gwt <- gw %>% as_tibble() %>%
    mutate(chrom = fct_chr(chr), bin = floor(ps / 1e4)) %>%
    group_by(chrom, bin) %>% slice_max(neglog10p, n = 1, with_ties = FALSE) %>%
    ungroup()

  sig     <- gwt %>% filter(neglog10p > BF)
  sig_eig <- gwt %>% filter(neglog10p > EIG, neglog10p <= BF)

  ## one label per threshold, in the chromosome-I facet only: annotate() would
  ## draw them in all six
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

  ggplot(gwt, aes(ps / 1e6, neglog10p)) +
    geom_blank(data = span, aes(x = pos.mb, y = 0), inherit.aes = FALSE) +
    geom_hline(yintercept = BF, linetype = "dashed", linewidth = 0.3,
               colour = COL_THR) +
    geom_hline(yintercept = EIG, linetype = "dotted", linewidth = 0.45,
               colour = COL_EIG) +
    geom_point(size = 0.4, alpha = 0.5, colour = COL_PT) +
    { if (nrow(sig_eig)) geom_point(data = sig_eig, colour = COL_EIG,
                                    size = 0.8) } +
    { if (nrow(sig)) geom_point(data = sig, colour = COL_PEAK, size = 1.2) } +
    geom_richtext(data = thr_lab, aes(x = pos.mb, y = y, label = label),
                  inherit.aes = FALSE, colour = thr_lab$col, size = 2.7,
                  hjust = 0, vjust = -0.35, fill = NA, label.color = NA,
                  label.padding = grid::unit(rep(0, 4), "pt")) +
    facet_grid(. ~ chrom, scales = "free_x", space = "free_x") +
    scale_x_continuous(breaks = seq(0, 25, 5), expand = expansion(mult = 0.02),
                       guide = guide_axis(check.overlap = TRUE)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = "Genomic Position (Mb)", y = "&minus;log<sub>10</sub>*p*",
         title = if (bare) panel_title(letter)
                 else titled(letter, "**Association scan**"),
         subtitle = if (bare) NULL else wrap_md(sprintf(paste0(
           "GEMMA LOCO on vst_ctrl_pos-1_T2, %s markers, n = %d strains. Grey ",
           "dashed is Bonferroni over every marker; green dotted divides alpha ",
           "by the %s effective independent tests. %s markers clear ",
           "Bonferroni, %s clear the eigen threshold."),
           format(nrow(gw), big.mark = ","), TH$n_strain,
           format(round(TH$m_eff), big.mark = ","),
           format(sum(gw$neglog10p > BF), big.mark = ","),
           format(sum(gw$neglog10p > EIG), big.mark = ",")))) +
    theme_pub(base_size) +
    theme(axis.title.y = element_markdown(),
          axis.text.x = element_text(size = 8))
}

## ===========================================================================
## output
## ===========================================================================
write_fig <- function(fig, stem, width, height) {
  ggsave(file.path(OUT, paste0(stem, ".pdf")), fig, width = width,
         height = height, device = cairo_pdf)
  ggsave(file.path(OUT, paste0(stem, ".png")), fig, width = width,
         height = height, dpi = 300, bg = "white")
  msg("wrote ", stem, ".{pdf,png}")
}

## ===========================================================================
## bootstrap uncertainty on the NNLS slopes
##
## supplemental_data/deconvolution/baugh_bootstrap_array.rda carries `ab`, the full
## bootstrap array: 102 strains x 23 samples x 100 replicates of the NNLS
## deconvolution, already normalised to sum to 1 within each sample. The legacy
## scripts loaded it, averaged it to `bootmean`/`bootse`, carried both through
## every join, and then plotted neither.
##
## The uncertainty that belongs on the slope panel is uncertainty on a SLOPE,
## not on a frequency, so the whole slope calculation is redone inside each
## bootstrap replicate: delta against the day-1 sample of the same replicate
## arm, regress on day, average over arms. That gives 100 slopes per strain and
## the interval is their percentiles. Propagating `bootse` through the
## regression instead would need an error model for a quantity that is already
## available by resampling, so it is not done.
##
## The slope dimension of `ab` is unnamed; its order is the column order of the
## genotype matrix, which is why the 30 MB input has to be opened once. The
## result is cached, so that happens only on the first run.
## ===========================================================================
## Strain order for the bootstrap array. Reads the shipped one-column file;
## falls back to the genotype matrix if only that is available, so the
## function works from either the supplemental deposit or the full archive.
baugh_strain_order <- function() {
  ## the shipped file carries a "strain" header line
  if (file.exists(SORDER)) {
    x <- readLines(SORDER)
    return(if (identical(x[1], "strain")) x[-1] else x)
  }
  if (!file.exists(BOOT))
    stop("need ", SORDER, " or ", BOOT, " to recover the bootstrap array's\n",
         "  strain order; see DATA_AVAILABILITY.md")
  g <- new.env(); load(BOOT, g)
  sub("_.*$", "", colnames(g$flipped_bootstrap_input[[1]]))
}

boot_slope_matrix <- function(refresh = FALSE) {
  if (!refresh && file.exists(BCACHE)) {
    S <- readRDS(BCACHE)
    msg("  bootstrap slopes from cache: ", nrow(S), " strains x ", ncol(S),
        " replicates")
    return(S)
  }

  e <- new.env(); load(BPRED, e)
  ab <- e$ab
  samples <- dimnames(ab)[[2]]

  strains <- baugh_strain_order()
  stopifnot(length(strains) == dim(ab)[1])

  meta <- tibble(sample = samples,
                 replicate = sub("_.*$", "", samples),
                 day = as.numeric(sub("^.*_d([0-9]+).*$", "\\1", samples)),
                 baseline = grepl("baseline", samples)) %>%
    filter(!baseline, day != 17)          # as in platform_slopes()
  arms <- split(meta, meta$replicate)
  msg("  bootstrap: ", dim(ab)[3], " replicates x ", length(arms),
      " arms, days ", paste(sort(unique(meta$day)), collapse = "/"))

  S <- matrix(NA_real_, dim(ab)[1], dim(ab)[3],
              dimnames = list(strains, NULL))
  for (b in seq_len(dim(ab)[3])) {
    per_arm <- vapply(arms, function(a) {
      idx <- match(a$sample, samples)
      d1  <- idx[a$day == 1]
      ## column-wise recycling subtracts each strain's own day-1 frequency
      delta <- ab[, idx, b] - ab[, d1, b]
      xc <- a$day - mean(a$day)
      as.vector(delta %*% xc) / sum(xc^2)
    }, numeric(dim(ab)[1]))
    S[, b] <- rowMeans(per_arm)
  }

  saveRDS(S, BCACHE)
  msg("  wrote ", BCACHE)
  S
}

boot_slope_summary <- function(S, level = 0.95) {
  a <- (1 - level) / 2
  tibble(strain    = rownames(S),
         boot_mean = rowMeans(S, na.rm = TRUE),
         boot_sd   = apply(S, 1, sd, na.rm = TRUE),
         lo        = apply(S, 1, quantile, probs = a, na.rm = TRUE),
         hi        = apply(S, 1, quantile, probs = 1 - a, na.rm = TRUE))
}

## the correlation itself has a bootstrap distribution: recompute rho against
## the MIP slopes inside each replicate. This is the number a reader should be
## given, rather than a single rho with no indication of how stable it is.
boot_rho <- function(S, mip, level = 0.95) {
  ok  <- intersect(rownames(S), names(mip))
  Sk  <- S[ok, , drop = FALSE]
  mk  <- mip[ok]
  rho <- apply(Sk, 2, function(y) cor(mk, y, method = "spearman",
                                      use = "complete.obs"))
  a <- (1 - level) / 2
  list(rho = rho, med = median(rho),
       lo = unname(quantile(rho, a)), hi = unname(quantile(rho, 1 - a)),
       n = length(ok))
}

## ---------------------------------------------------------------------------
## bootstrap intervals on the individual sample-level frequencies
##
## The same array, summarised one level earlier: a percentile interval per
## (strain, sample) rather than per strain. These are the estimates that go
## into the deltas, which go into the slopes, so this is what the slope panel
## aggregates over. Cached because the reduction is slow, not because the
## inputs are large.
## ---------------------------------------------------------------------------
boot_freq_intervals <- function(refresh = FALSE, level = 0.95) {
  if (!refresh && file.exists(FCACHE)) {
    d <- readRDS(FCACHE)
    msg("  bootstrap frequency intervals from cache: ", nrow(d), " rows")
    return(d)
  }
  e <- new.env(); load(BPRED, e)
  ab <- e$ab
  samples <- dimnames(ab)[[2]]
  strains <- baugh_strain_order()
  stopifnot(length(strains) == dim(ab)[1])

  a  <- (1 - level) / 2
  q  <- function(pr) apply(ab, c(1, 2), quantile, probs = pr, na.rm = TRUE)
  lo <- q(a); hi <- q(1 - a)
  sdm <- apply(ab, c(1, 2), sd, na.rm = TRUE)
  mn  <- apply(ab, c(1, 2), mean, na.rm = TRUE)

  flat <- function(m, nm) {
    dimnames(m) <- list(strain = strains, sample = samples)
    as_tibble(as.table(m), .name_repair = "minimal") %>%
      setNames(c("strain", "sample", nm))
  }
  d <- reduce(list(flat(lo, "lo"), flat(hi, "hi"), flat(sdm, "boot_sd"),
                   flat(mn, "boot_mean")),
              left_join, by = c("strain", "sample"))
  saveRDS(d, FCACHE)
  msg("  wrote ", FCACHE)
  d
}

## ===========================================================================
## a note on aligning a panel that spans columns
##
## For a layout with two panels on top and one full-width panel below, USE
## NESTING -- (pA | pB) / pC -- not a flat plot_layout(design = ...) with an
## area() spanning both columns. Verified by rendering with coloured panel
## backgrounds and measuring the pixels:
##
##   flat design, area(2, 1, 2, 2) for C : row 1 starts at x = 220, C at 117
##   nested, (pA | pB) / pC              : both rows start at x = 117
##
## Under the flat design the spanning panel keeps its own axis gutter and its
## panel begins inside the column-1 gutter track, so C is wider than A + B by
## exactly A's gutter -- here 103 px, "NNLS slope" plus 0.000/0.002/0.004
## against "-log10 p" plus 0.0/2.5/5.0. Nesting puts both gutters in the same
## outer track, sized to the wider of the two, and the panel edges coincide.
##
## Padding the narrower gutter does NOT fix the flat case: widening C's gutter
## widens the outer track and shifts row 1 right by the same amount, so the gap
## survives. That was tried and removed; the layout is the fix.
##
## Note this is the opposite conclusion from Figure3_common.R, where the
## spanning panel is on TOP and the panels below carry no y-axis furniture at
## all. There the flat design is what aligns. The rule is not "always nest" or
## "always flatten" -- it is to measure.
## ===========================================================================

## ===========================================================================
## the depth curve at the level of individual samples
##
## panel_downsample() correlates SLOPES, which are fitted across three
## timepoints and averaged over five replicate arms, so a sample cannot be seen
## in it. This is the same question one level down: at each depth, how well
## does each individual sample's frequency estimate agree with the published
## MIP-seq frequencies for that sample, across strains?
##
## All six depths are kept here, including the 0.5x that panel_downsample()
## drops. That drop is inherited from the legacy script and the reason is
## visible in the aggregate: at the slope level 0.5x scores 0.870, ABOVE the
## 0.851 of the 1x that follows it, so the aggregate curve is not monotone. Per
## sample there is no such inversion, which suggests the aggregate wobble is
## noise in one downsampling draw rather than anything about 0.5x.
##
## No bootstrap intervals: `ab` holds replicates of the full-depth
## deconvolution only, so there is nothing to resample at a reduced depth.
## ===========================================================================
downsample_per_sample <- function(freq) {
  e <- new.env(); load(DS, e)
  e$ds_predictions_df %>% filter(strain != "N2") %>%
    left_join(freq %>% select(strain, sample, published_frq),
              by = c("strain", "sample")) %>%
    filter(is.finite(ds_frq), is.finite(published_frq)) %>%
    group_by(sample, ds_n) %>%
    summarise(rho = cor(ds_frq, published_frq, method = "spearman"),
              n = n(), .groups = "drop")
}

## full-depth per-sample agreement, the value the curves should approach
full_depth_per_sample <- function(freq) {
  freq %>% filter(strain != "N2") %>% group_by(sample) %>%
    summarise(rho = cor(frq, published_frq, method = "spearman",
                        use = "complete.obs"), .groups = "drop")
}

panel_downsample_samples <- function(freq, letter = "A", bare = TRUE,
                                     base_size = 11.5) {
  d  <- downsample_per_sample(freq)
  fd <- full_depth_per_sample(freq)

  med <- d %>% group_by(ds_n) %>%
    summarise(rho = median(rho), .groups = "drop")
  ## the sample that agrees worst, averaged over depths -- named rather than
  ## left as an anonymous low line
  worst <- d %>% group_by(sample) %>%
    summarise(m = median(rho), .groups = "drop") %>% slice_min(m, n = 1)
  hl <- d %>% filter(sample == worst$sample)

  msg("  panel per-sample depth: ", n_distinct(d$sample), " samples x ",
      n_distinct(d$ds_n), " depths")
  msg("    median rho: ",
      paste(sprintf("%gx %.3f", med$ds_n, med$rho), collapse = " | "))
  msg("    worst sample: ", worst$sample, " median ", round(worst$m, 3),
      " (range ", round(min(hl$rho), 3), "-", round(max(hl$rho), 3), ")")
  msg("    full depth per-sample rho: median ", round(median(fd$rho), 3),
      " | min ", round(min(fd$rho), 3))

  ggplot(d, aes(ds_n, rho)) +
    ## full-depth median, the ceiling the curves are climbing towards
    geom_hline(yintercept = median(fd$rho), linetype = "dashed",
               linewidth = 0.35, colour = "grey55") +
    geom_line(aes(group = sample), linewidth = 0.3, colour = "grey70",
              alpha = 0.8) +
    geom_line(data = hl, aes(group = sample), linewidth = 0.7,
              colour = COL_FIT) +
    geom_point(data = hl, size = 1.3, colour = COL_FIT) +
    geom_line(data = med, linewidth = 1, colour = COL_PT) +
    geom_point(data = med, size = 2.4, colour = COL_PT) +
    annotate("richtext", x = max(d$ds_n), y = median(fd$rho),
             label = sprintf("full depth, median %.2f", median(fd$rho)),
             hjust = 1, vjust = -0.5, size = 2.6, colour = "grey35",
             fill = NA, label.color = NA,
             label.padding = grid::unit(rep(0, 4), "pt")) +
    geom_richtext(data = hl %>% slice_max(ds_n, n = 1),
                  aes(label = sample), hjust = 1.1, vjust = 1.5, size = 2.6,
                  colour = COL_FIT, fill = NA, label.color = NA,
                  label.padding = grid::unit(rep(0, 4), "pt")) +
    scale_x_log10(breaks = sort(unique(d$ds_n)),
                  labels = function(x) paste0(x, "×")) +
    labs(x = "Sequencing depth", y = "Spearman's &rho; vs MIP-seq",
         title = if (bare) panel_title(letter)
                 else titled(letter, "**Per-sample agreement against depth**"),
         subtitle = if (bare) NULL else wrap_md(paste0(
           "One grey line per sample, ", n_distinct(d$sample),
           " samples; the dark line is the median across samples. Dashed line ",
           "is the median full-depth agreement. The lowest line is named."))) +
    theme_pub(base_size) +
    theme(axis.title.y = element_markdown(),
          panel.grid.minor.x = element_blank())
}
