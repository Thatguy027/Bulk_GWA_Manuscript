library(tidyverse)
library(glue)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))

# ---------------------------------------------------------------------------
# Export the per-strain pos-1 RNAi "change in frequency" phenotype used in
# plots/original_pos1_strain_freq_change.pdf (built in pos1_script1.R,
# ~line 178-261). Reconstructs t2_ctrl/t2_pos1 from the saved bootstrap
# prediction table rather than rerunning the full genotype/bootstrap
# pipeline (which needs large intermediate files not otherwise touched
# here).
# ---------------------------------------------------------------------------

experiment_folder <- "../data/pos1_original/"

load(glue::glue("{experiment_folder}2024pub_prediction_bootstrap_combined_df.Rdata"))

t2_ctrl <- prediction_boot_df %>%
  dplyr::filter(time == "t2", condition %in% c("ctrlA", "ctrlB")) %>%
  dplyr::mutate(rep = ifelse(grepl("B", condition), "B", "A")) %>%
  dplyr::select(strain, btm_ctrl = bootmean, bts_ctrl = bootse, ctrl_frq = new_frq, control_cond = condition, time, rep)

t2_pos1 <- prediction_boot_df %>%
  dplyr::filter(time == "t2", condition %in% c("pos1B_pos1B", "pos1A_pos1A")) %>%
  dplyr::mutate(rep = ifelse(grepl("B", condition), "B", "A")) %>%
  dplyr::left_join(., t2_ctrl, by = c("strain", "time", "rep"))

# per strain x replicate (A/B): matches the two lines-per-strain shown in the
# figure (Control -> pos-1 RNAi, faceted by replicate)
pos1_strain_freq_change <- t2_pos1 %>%
  dplyr::transmute(
    strain,
    rep,
    ctrl_frq_boot = btm_ctrl,
    ctrl_frq_boot_se = bts_ctrl,
    pos1_frq_boot = bootmean,
    pos1_frq_boot_se = bootse,
    delta_freq = bootmean - btm_ctrl,
    direction = dplyr::if_else(delta_freq < 0, "decrease", "increase")
  ) %>%
  dplyr::arrange(strain, rep)

readr::write_tsv(pos1_strain_freq_change, glue::glue("{experiment_folder}pos1_strain_freq_change.tsv"))

# per strain, mean of the two replicates (convenience summary, not what's
# directly plotted in the figure but useful for downstream comparisons)
pos1_strain_freq_change_mean <- pos1_strain_freq_change %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(
    mean_delta_freq = mean(delta_freq),
    n_reps = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(mean_delta_freq))

readr::write_tsv(pos1_strain_freq_change_mean, glue::glue("{experiment_folder}pos1_strain_freq_change_mean.tsv"))
