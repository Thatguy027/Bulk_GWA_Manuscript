library(tidyverse)
library(ggpubr)
library(glue)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

dir.create("plot/cor_plots", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# load processed strain frequencies (strain x sample NNLS deconvolution output,
# see scripts/20240328_calculate_frequencies.R)
# ---------------------------------------------------------------------------

load("data/processed_strain_fq.rda")

fq <- processed_strain_frq %>%
  dplyr::filter(!is.na(strain))

# ---------------------------------------------------------------------------
# metadata correction: sample "43" (rde-3, supernatant, t2, bacteria pellet
# tube 10) and sample "44" (ctrl, supernatant, t2, tube 11) were both recorded
# with SAMPLE NUMBER 29 in meta_files/experiment_meta.tsv, so both draw from
# the same "S29" prediction column in processed_strain_fq.rda. Sample "S28" is
# present in the frequency table but was orphaned (no metadata match) -- this
# is the true rde-3 sample. Confirmed by comparing delta-frequency spearman
# correlations of S28/S29 against the independent rde-3 and ctrl replicates
# (S28 correlates more strongly with the rde-3 rep2/pellet-tube10 samples than
# S29 does), consistent with a SAMPLE NUMBER typo (should read 28, not 29) in
# the rde-3 supernatant row. The fix below re-labels the orphaned S28 rows as
# rde-3/supernatant/t2 and drops the mislabeled rde-3 copy of S29.
# ---------------------------------------------------------------------------

rde3_s28_fix <- fq %>%
  dplyr::filter(sample == "S28") %>%
  dplyr::mutate(
    `SAMPLE NAME`    = "43",
    `SAMPLE NUMBER`  = 28,
    BACTERIA_PELLET  = 10,
    sample_type      = "supernatant",
    Gene             = "rde-3",
    Timepoint        = 2
  )

fq <- fq %>%
  dplyr::filter(!(sample == "S29" & Gene == "rde-3")) %>%
  dplyr::bind_rows(rde3_s28_fix)

# ---------------------------------------------------------------------------
# t0 baseline: mean of the two t=0 input samples (0A, 0B)
# ---------------------------------------------------------------------------

t0 <- fq %>%
  dplyr::filter(Timepoint == 0) %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(t0_frq = mean(frq), .groups = "drop")

# ---------------------------------------------------------------------------
# delta frequency (t2 - t0) for every t2 sample, then assign replicate 1/2 per
# Gene x sample_type condition, ordered by SAMPLE NUMBER (block A = rep 1,
# block B = rep 2). Any extra samples beyond the first two per condition
# (e.g. the "33_b"/"34_b" repeat-library samples) are dropped, and conditions
# with fewer than two usable replicates (e.g. mig-6 supernatant, S48, which
# failed frequency deconvolution) are excluded and noted in the summary table.
# ---------------------------------------------------------------------------

# S19 (ctrl, pellet, tube 11, block A) and S20 (ric-3, supernatant, tube 1,
# block A) are excluded here, matching the precedent in
# scripts/20240328_calculate_frequencies.R. Both tubes have a "_b" re-prep
# (33_b = S54, 34_b = S55) -- the exclusion of S19/S20 in that script implies
# the originals were flagged as failed/low-quality and the "_b" samples are
# their replacements, not redundant extra replicates. Excluding S19/S20 here
# leaves exactly S40+S54 for ctrl-pellet and S41+S55 for ric-3-supernatant,
# which then get assigned rep1/rep2 by the usual SAMPLE NUMBER ordering below.
t2 <- fq %>%
  dplyr::filter(Timepoint == 2) %>%
  dplyr::filter(!sample %in% c("S19", "S20")) %>%
  dplyr::left_join(t0, by = "strain") %>%
  dplyr::mutate(delta_fq = frq - t0_frq) %>%
  dplyr::filter(!is.na(delta_fq))

rep_lookup <- t2 %>%
  dplyr::distinct(Gene, sample_type, sample, `SAMPLE NUMBER`) %>%
  dplyr::arrange(Gene, sample_type, `SAMPLE NUMBER`) %>%
  dplyr::group_by(Gene, sample_type) %>%
  dplyr::mutate(rep = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(rep <= 2)

t2 <- t2 %>%
  dplyr::inner_join(
    rep_lookup %>% dplyr::select(Gene, sample_type, sample, rep),
    by = c("Gene", "sample_type", "sample")
  )

all_genes <- sort(unique(t2$Gene))
all_preps <- sort(unique(t2$sample_type))

spearman_stat <- function(x, y) {
  suppressWarnings(cor.test(x, y, method = "spearman"))
}

cor_summary <- tibble::tibble()

# ---------------------------------------------------------------------------
# Comparison 1: replicate 1 vs replicate 2 delta frequency (t2 - t0), for
# every RNAi condition and the control, within each sample_type (pellet /
# supernatant)
#
# Not used downstream -- commented out so it doesn't run or write to
# plot/cor_plots or the summary table.
# ---------------------------------------------------------------------------

# for (gene in all_genes) {
#   for (prep in all_preps) {
#
#     cond_df <- t2 %>%
#       dplyr::filter(Gene == gene, sample_type == prep) %>%
#       dplyr::select(strain, rep, delta_fq)
#
#     n_reps <- length(unique(cond_df$rep))
#
#     if (n_reps < 2) {
#       cor_summary <- dplyr::bind_rows(cor_summary, tibble::tibble(
#         comparison = "replicate_vs_replicate", Gene = gene, sample_type = prep,
#         rep = NA_character_, n_strains = NA_integer_, spearman_rho = NA_real_,
#         p_value = NA_real_, note = glue::glue("skipped - only {n_reps} usable replicate(s)")
#       ))
#       next
#     }
#
#     wide_df <- cond_df %>%
#       tidyr::pivot_wider(names_from = rep, values_from = delta_fq, names_prefix = "rep") %>%
#       tidyr::drop_na(rep1, rep2)
#
#     ct <- spearman_stat(wide_df$rep1, wide_df$rep2)
#     rho <- round(unname(ct$estimate), 3)
#
#     cor_summary <- dplyr::bind_rows(cor_summary, tibble::tibble(
#       comparison = "replicate_vs_replicate", Gene = gene, sample_type = prep,
#       rep = "rep1_vs_rep2", n_strains = nrow(wide_df), spearman_rho = rho,
#       p_value = ct$p.value, note = NA_character_
#     ))
#
#     p <- ggplot(wide_df, aes(x = rep1, y = rep2)) +
#       geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
#       geom_point(alpha = 0.7) +
#       ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", color = "firebrick3") +
#       labs(
#         title = glue::glue("{gene} ({prep}), t2 - t0: replicate 1 vs replicate 2"),
#         x = "Delta frequency (t2 - t0), replicate 1",
#         y = "Delta frequency (t2 - t0), replicate 2"
#       ) +
#       theme_bw(14)
#
#     ggsave(p, filename = glue::glue("plot/cor_plots/{gene}_{prep}_t2-t0_rep1_vs_rep2.pdf"), height = 6, width = 6)
#   }
# }

# ---------------------------------------------------------------------------
# Comparison 2: (t2 RNAi condition - t2 ctrl) raw frequency [the RNAi
# condition's t2 frq minus ctrl's t2 frq -- no t0 baseline correction],
# replicate 1 vs replicate 2 -- checks whether the control-relative RNAi
# effect reproduces across the two replicate blocks, within each sample_type.
#
# ctrl_frq is the MEAN of ctrl's two replicates per strain x sample_type (not
# matched rep-by-rep to the RNAi replicate) -- this matches the "t2 - t2c"
# convention in scripts/20240328_calculate_frequencies.R (`mt2c`), where a
# single pooled ctrl baseline is subtracted from every individual RNAi
# replicate.
# ---------------------------------------------------------------------------

ctrl_mean_df <- t2 %>%
  dplyr::filter(Gene == "ctrl") %>%
  dplyr::group_by(strain, sample_type) %>%
  dplyr::summarise(ctrl_frq = mean(frq), .groups = "drop")

rnai_genes <- setdiff(all_genes, "ctrl")

for (gene in rnai_genes) {
  for (prep in all_preps) {

    rnai_df <- t2 %>%
      dplyr::filter(Gene == gene, sample_type == prep) %>%
      dplyr::select(strain, rep, rnai_frq = frq)

    diff_df <- dplyr::inner_join(
        rnai_df,
        ctrl_mean_df %>% dplyr::filter(sample_type == prep) %>% dplyr::select(strain, ctrl_frq),
        by = "strain"
      ) %>%
      dplyr::mutate(rnai_minus_ctrl = rnai_frq - ctrl_frq) %>%
      dplyr::select(strain, rep, rnai_minus_ctrl)

    n_reps <- length(unique(diff_df$rep))

    if (n_reps < 2) {
      cor_summary <- dplyr::bind_rows(cor_summary, tibble::tibble(
        comparison = "rnai_minus_ctrl_rep1_vs_rep2", Gene = gene, sample_type = prep,
        rep = NA_character_, n_strains = NA_integer_, spearman_rho = NA_real_,
        p_value = NA_real_, note = glue::glue("skipped - only {n_reps} usable matched replicate(s)")
      ))
      next
    }

    wide_df <- diff_df %>%
      tidyr::pivot_wider(names_from = rep, values_from = rnai_minus_ctrl, names_prefix = "rep") %>%
      tidyr::drop_na(rep1, rep2)

    ct <- spearman_stat(wide_df$rep1, wide_df$rep2)
    rho <- round(unname(ct$estimate), 3)

    cor_summary <- dplyr::bind_rows(cor_summary, tibble::tibble(
      comparison = "rnai_minus_ctrl_rep1_vs_rep2", Gene = gene, sample_type = prep,
      rep = "rep1_vs_rep2", n_strains = nrow(wide_df), spearman_rho = rho,
      p_value = ct$p.value, note = NA_character_
    ))

    p <- ggplot(wide_df, aes(x = rep1, y = rep2)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
      geom_point(alpha = 0.7) +
      ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", color = "firebrick3") +
      labs(
        title = glue::glue("{gene} ({prep}): ({gene} - ctrl) t2 frequency, replicate 1 vs replicate 2"),
        x = "(RNAi - ctrl) t2 frequency, replicate 1",
        y = "(RNAi - ctrl) t2 frequency, replicate 2"
      ) +
      theme_bw(14)

    ggsave(p, filename = glue::glue("plot/cor_plots/{gene}_{prep}_rnai-minus-ctrl_rep1_vs_rep2.pdf"), height = 6, width = 6)
  }
}

# ---------------------------------------------------------------------------
# summary correlation table
# ---------------------------------------------------------------------------

readr::write_tsv(cor_summary, "plot/cor_plots/correlation_summary_table.tsv")

# ---------------------------------------------------------------------------
# Flagged strains: for each RNAi gene x sample_type, rank all strains by their
# mean (RNAi - ctrl) t2 frequency (averaged across rep1/rep2) and flag the 5
# largest increases and 5 largest decreases. rank 1 = largest increase, rank
# n_strains = largest decrease.
# ---------------------------------------------------------------------------

rnai_vs_ctrl_mean <- t2 %>%
  dplyr::filter(Gene != "ctrl") %>%
  dplyr::select(strain, Gene, sample_type, rep, frq) %>%
  dplyr::inner_join(
    ctrl_mean_df %>% dplyr::select(strain, sample_type, ctrl_frq),
    by = c("strain", "sample_type")
  ) %>%
  dplyr::mutate(rnai_minus_ctrl = frq - ctrl_frq) %>%
  dplyr::group_by(strain, Gene, sample_type) %>%
  dplyr::summarise(
    mean_rnai_minus_ctrl = mean(rnai_minus_ctrl),
    n_reps = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_reps == 2) %>%
  dplyr::group_by(Gene, sample_type) %>%
  dplyr::mutate(
    rank = dplyr::min_rank(dplyr::desc(mean_rnai_minus_ctrl)),
    n_strains = dplyr::n()
  ) %>%
  dplyr::ungroup()

flagged_strains <- rnai_vs_ctrl_mean %>%
  dplyr::group_by(Gene, sample_type) %>%
  dplyr::group_modify(~ dplyr::bind_rows(
    dplyr::slice_max(.x, mean_rnai_minus_ctrl, n = 5, with_ties = FALSE) %>%
      dplyr::mutate(direction = "largest_increase"),
    dplyr::slice_min(.x, mean_rnai_minus_ctrl, n = 5, with_ties = FALSE) %>%
      dplyr::mutate(direction = "largest_decrease")
  )) %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_type, Gene, direction, strain, mean_rnai_minus_ctrl, rank, n_strains) %>%
  dplyr::arrange(sample_type, Gene, direction, rank)

readr::write_tsv(flagged_strains, "plot/cor_plots/top_bottom5_strains_by_gene.tsv")

# ---------------------------------------------------------------------------
# Heatmap (one per sample_type): rows = strain, columns = RNAi gene, tile
# color = that gene's rank (1..n_strains) for that strain, red (largest
# increase) - white - blue (largest decrease). Strains and genes are ordered
# by hierarchical clustering (Euclidean distance, complete linkage) on the
# rank matrix, rather than alphanumerically, to surface strains/genes with
# similar response patterns.
#
# Within each gene's own column, the top 3 / bottom 3 strains are marked with
# stacked black arrows: 3 arrows for rank 1 (or n_strains), 2 for rank 2 (or
# n_strains-1), 1 for rank 3 (or n_strains-2) -- up arrows for the best
# (largest increase), down arrows for the worst (largest decrease).
# ---------------------------------------------------------------------------

for (prep in all_preps) {

  flagged_ids <- flagged_strains %>%
    dplyr::filter(sample_type == prep) %>%
    dplyr::distinct(strain) %>%
    dplyr::pull(strain)

  if (length(flagged_ids) == 0) next

  heat_df <- rnai_vs_ctrl_mean %>%
    dplyr::filter(sample_type == prep, strain %in% flagged_ids)

  n_strains_prep <- max(heat_df$n_strains, na.rm = TRUE)

  # cluster strains (rows) and genes (columns) on the rank matrix
  mat <- heat_df %>%
    dplyr::select(strain, Gene, rank) %>%
    tidyr::pivot_wider(names_from = Gene, values_from = rank) %>%
    tibble::column_to_rownames("strain") %>%
    as.matrix()

  mat_imputed <- apply(mat, 2, function(col) {
    col[is.na(col)] <- mean(col, na.rm = TRUE)
    col
  })

  strain_order <- rownames(mat)[hclust(dist(mat_imputed))$order]
  gene_order <- colnames(mat)[hclust(dist(t(mat_imputed)))$order]

  heat_df <- heat_df %>%
    dplyr::mutate(
      strain = factor(strain, levels = strain_order),
      Gene = factor(Gene, levels = gene_order)
    )

  arrow_df <- heat_df %>%
    dplyr::filter(rank <= 3 | rank > n_strains - 3) %>%
    dplyr::mutate(
      n_arrows = dplyr::if_else(rank <= 3, 4 - rank, rank - (n_strains - 3)),
      symbol = dplyr::if_else(rank <= 3, "▲", "▼"),
      arrow_label = purrr::map2_chr(symbol, n_arrows, ~ paste(rep(.x, .y), collapse = "\n"))
    )

  hp <- ggplot(heat_df, aes(x = Gene, y = strain, fill = rank)) +
    geom_tile(color = "grey90") +
    geom_text(
      data = arrow_df, aes(x = Gene, y = strain, label = arrow_label),
      color = "black", size = 1.8, lineheight = 0.7, inherit.aes = FALSE
    ) +
    scale_fill_gradient2(
      low = "red", mid = "white", high = "blue",
      midpoint = (1 + n_strains_prep) / 2,
      name = "Rank\n(1 = largest\nincrease)"
    ) +
    theme_bw(12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(
      title = glue::glue("Flagged strains ({prep}): rank of (RNAi - ctrl) t2 frequency"),
      x = "Gene", y = "Strain"
    )

  ggsave(hp, filename = glue::glue("plot/cor_plots/flagged_strains_heatmap_{prep}.pdf"), height = 10, width = 7, device = cairo_pdf)
}
