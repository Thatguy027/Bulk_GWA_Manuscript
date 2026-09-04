library(tidyverse)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# ---------------------------------------------------------------------------
# Builds a single clean, analysis-ready sample sheet for the pooled_RNAi_expt
# bulk competition experiment, programmatically applying every known
# metadata correction on top of meta_files/experiment_meta.tsv, rather than
# hand-editing that source file (keeps every fix auditable here and
# re-derivable if experiment_meta.tsv is ever regenerated).
#
# Output: meta_files/experiment_meta_clean.tsv
#
# Key columns for downstream use:
#   sample_number  - corrected SAMPLE NUMBER (see rde3_s29_collision fix below)
#   sample_id      - "S<sample_number>", the join key used in
#                    data/processed_strain_fq.rda ("sample" column)
#   aser_file_a/b  - expected raw ASE-reader filenames in data/aser/, of the
#                    form S<sample_number>_A_1.table / S<sample_number>_B_1.table
#                    (see scripts/20240320_process_genotypes_n_counts.R, which
#                    parses "<sample_id>_<lane A/B>_<junk>.table")
#   exclude        - TRUE for samples that should be dropped before any
#                    replicate/delta-frequency analysis
#   exclude_reason - why
#   rep            - replicate index (1/2) within Gene x sample_type x
#                    Timepoint, assigned by ascending sample_number among
#                    non-excluded rows only
# ---------------------------------------------------------------------------

meta_raw <- data.table::fread("meta_files/experiment_meta.tsv") %>%
  dplyr::rename(
    sample_name     = `SAMPLE NAME`,
    sample_number   = `SAMPLE NUMBER`,
    bacteria_pellet = BACTERIA_PELLET,
    gene            = Gene,
    timepoint       = Timepoint
  )

# ---------------------------------------------------------------------------
# Fix 1: rde-3/ctrl SAMPLE NUMBER 29 collision (t2 supernatant).
#
# experiment_meta.tsv records BOTH sample_name "43" (rde-3, supernatant, t2,
# bacteria_pellet tube 10) and sample_name "44" (ctrl, supernatant, t2, tube
# 11) with sample_number 29, so both would join to the same "S29" prediction
# column in processed_strain_fq.rda. sample_number 28 is otherwise unused
# (the table jumps 27 -> 29) but IS present as an orphaned column in the
# frequency table, and correlates more strongly with the independent rde-3
# replicates than S29 does -- i.e. this is a typo (should read 28, not 29)
# on the rde-3 row, not a real collision. The true pairing is:
#   sample_name "43" (rde-3) -> sample_number 28 -> S28
#   sample_name "44" (ctrl)  -> sample_number 29 -> S29  (unchanged)
# ---------------------------------------------------------------------------

meta_fixed <- meta_raw %>%
  dplyr::mutate(
    sample_number = dplyr::if_else(sample_name == "43" & gene == "rde-3", 28L, as.integer(sample_number))
  )

# ---------------------------------------------------------------------------
# Fix 2 & 3: exclusions.
#
# - sample_number 19 (ctrl, pellet, tube 11, "33") and sample_number 20
#   (ric-3, supernatant, tube 1, "34") were flagged bad/low-quality; the
#   "_b" re-preps of the SAME tubes (sample_name "33_b" -> sample_number 54,
#   "34_b" -> sample_number 55) are their replacements, not extra redundant
#   replicates. Precedent: scripts/20240328_calculate_frequencies.R
#   explicitly filters `!sample %in% c("S19","S20")` wherever replicate
#   correlations/deltas are computed.
# - sample_number 48 (mig-6, supernatant, tube 19, "63") failed NNLS
#   frequency deconvolution (all-NA frq in processed_strain_fq.rda) and is
#   dropped everywhere downstream (e.g. `select(-supernatant_mig-6_S48)` in
#   the same legacy script).
# ---------------------------------------------------------------------------

exclusions <- tibble::tribble(
  ~sample_number, ~exclude_reason,
  19L, "flagged bad/low-quality original prep; replaced by '33_b' (sample_number 54)",
  20L, "flagged bad/low-quality original prep; replaced by '34_b' (sample_number 55)",
  48L, "failed NNLS frequency deconvolution (all-NA frq in processed_strain_fq.rda)"
)

meta_clean <- meta_fixed %>%
  dplyr::mutate(sample_id = paste0("S", sample_number)) %>%
  dplyr::left_join(exclusions, by = "sample_number") %>%
  dplyr::mutate(
    exclude       = !is.na(exclude_reason),
    aser_file_a   = paste0(sample_id, "_A_1.table"),
    aser_file_b   = paste0(sample_id, "_B_1.table")
  ) %>%
  dplyr::group_by(gene, sample_type, timepoint) %>%
  dplyr::arrange(sample_number, .by_group = TRUE) %>%
  dplyr::mutate(rep = dplyr::if_else(exclude, NA_integer_, cumsum(!exclude))) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sample_number) %>%
  dplyr::select(
    sample_number, sample_id, aser_file_a, aser_file_b, sample_name,
    bacteria_pellet, sample_type, gene, timepoint, rep, exclude, exclude_reason
  )

# sanity checks -------------------------------------------------------------

stopifnot(
  "sample_id must be unique after the sample_number fix" =
    !anyDuplicated(meta_clean$sample_id)
)

n_excluded <- sum(meta_clean$exclude)
message(glue::glue("{nrow(meta_clean)} samples total, {n_excluded} excluded, {nrow(meta_clean) - n_excluded} usable"))

readr::write_tsv(meta_clean, "meta_files/experiment_meta_clean.tsv")
