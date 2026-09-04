library(tidyverse)
library(ggpubr)
library(ggrepel)
library(glue)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

dir.create("plot/paaby_comparison", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Our data: (RNAi - ctrl) t2 frequency, mean of replicate 1/2, for the genes
# shared with Paaby et al. 2015 -- pos-1 AND rpn-12 (both are shared RNAi
# clones between the two experiments, not just pos-1). Same conventions as
# scripts/20260813_delta_freq_correlations.R: S19/S20 excluded (bad
# libraries, replaced by the "33_b"/"34_b" re-preps), ctrl baseline is the
# mean of ctrl's two replicates per strain x sample_type, no t0 correction.
# ---------------------------------------------------------------------------

load("data/processed_strain_fq.rda")

fq <- processed_strain_frq %>%
  dplyr::filter(!is.na(strain))

t2 <- fq %>%
  dplyr::filter(Timepoint == 2) %>%
  dplyr::filter(!sample %in% c("S19", "S20"))

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

ctrl_mean_df <- t2 %>%
  dplyr::filter(Gene == "ctrl") %>%
  dplyr::group_by(strain, sample_type) %>%
  dplyr::summarise(ctrl_frq = mean(frq), .groups = "drop")

shared_genes <- c("pos-1", "rpn-12")

our_metric <- t2 %>%
  dplyr::filter(Gene %in% shared_genes) %>%
  dplyr::select(strain, Gene, sample_type, rep, frq) %>%
  dplyr::inner_join(ctrl_mean_df, by = c("strain", "sample_type")) %>%
  dplyr::mutate(rnai_minus_ctrl = frq - ctrl_frq) %>%
  dplyr::group_by(strain, Gene, sample_type) %>%
  dplyr::summarise(
    mean_rnai_minus_ctrl = mean(rnai_minus_ctrl),
    n_reps = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_reps == 2)

# ---------------------------------------------------------------------------
# Paaby et al. 2015 embryonic lethality data (Dryad doi:10.5061/dryad.d5j06),
# well-level manual/automated counts of adults/larvae/eggs per strain x RNAi
# vector. Embryonic lethality per well = eggs / (eggs + larvae) -- the
# fraction of laid eggs that failed to hatch into larvae. Wells with NA or
# zero (eggs + larvae) are dropped. Mean taken across all wells per strain x
# vector (pools plate/bacterial replicates and dates).
# ---------------------------------------------------------------------------

paaby_raw <- data.table::fread("paaby2015/emb_leth_data.txt")

paaby_metric <- paaby_raw %>%
  dplyr::filter(vector %in% shared_genes) %>%
  dplyr::filter(!is.na(eggs), !is.na(larvae), (eggs + larvae) > 0) %>%
  dplyr::mutate(leth = eggs / (eggs + larvae)) %>%
  dplyr::group_by(strain, vector) %>%
  dplyr::summarise(mean_leth = mean(leth), n_wells = dplyr::n(), .groups = "drop") %>%
  dplyr::rename(Gene = vector)

# ---------------------------------------------------------------------------
# overlap strains and comparison, per Gene x sample_type. Spearman rank
# correlation is used deliberately: only 9 strains are shared between the two
# strain panels, so a rank-based measure is more robust to the very different
# scales/units of the two phenotypes (ctrl-corrected pool frequency change vs
# embryonic lethality fraction) and to the small sample size -- but note that
# with n=9, correlation estimates and p-values are inherently low-powered and
# should be read as suggestive, not confirmatory.
# ---------------------------------------------------------------------------

overlap_strains <- intersect(unique(our_metric$strain), unique(paaby_metric$strain))

message(glue::glue("Overlapping strains (n={length(overlap_strains)}): {paste(sort(overlap_strains), collapse = ', ')}"))

comparison_df <- our_metric %>%
  dplyr::filter(strain %in% overlap_strains) %>%
  dplyr::inner_join(
    paaby_metric %>% dplyr::filter(strain %in% overlap_strains),
    by = c("strain", "Gene")
  )

readr::write_tsv(comparison_df, "plot/paaby_comparison/paaby_comparison_data.tsv")

cor_summary <- tibble::tibble()

for (gene in shared_genes) {
  for (prep in sort(unique(comparison_df$sample_type))) {

    pair_df <- comparison_df %>%
      dplyr::filter(Gene == gene, sample_type == prep)

    if (nrow(pair_df) < 3) {
      cor_summary <- dplyr::bind_rows(cor_summary, tibble::tibble(
        Gene = gene, sample_type = prep, n_strains = nrow(pair_df),
        spearman_rho = NA_real_, p_value = NA_real_,
        note = "skipped - fewer than 3 overlapping strains"
      ))
      next
    }

    ct <- suppressWarnings(cor.test(pair_df$mean_rnai_minus_ctrl, pair_df$mean_leth, method = "spearman"))
    rho <- round(unname(ct$estimate), 3)

    cor_summary <- dplyr::bind_rows(cor_summary, tibble::tibble(
      Gene = gene, sample_type = prep, n_strains = nrow(pair_df),
      spearman_rho = rho, p_value = ct$p.value, note = NA_character_
    ))

    p <- ggplot(pair_df, aes(x = mean_rnai_minus_ctrl, y = mean_leth)) +
      geom_point(size = 2.5) +
      ggrepel::geom_text_repel(aes(label = strain), size = 3, min.segment.length = 0, seed = 42) +
      ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", color = "firebrick3") +
      theme_bw(14) +
      theme(plot.title = element_text(size = 11)) +
      labs(
        title = glue::glue("{gene} ({prep}): this study vs Paaby et al. 2015 (n = {nrow(pair_df)} shared strains)"),
        x = "This study: mean (RNAi - ctrl) t2 frequency",
        y = "Paaby et al. 2015: mean embryonic lethality [eggs / (eggs+larvae)]"
      )

    ggsave(p, filename = glue::glue("plot/paaby_comparison/{gene}_{prep}_vs_paaby2015.pdf"), height = 6, width = 7.5)
  }
}

readr::write_tsv(cor_summary, "plot/paaby_comparison/paaby_comparison_summary.tsv")
