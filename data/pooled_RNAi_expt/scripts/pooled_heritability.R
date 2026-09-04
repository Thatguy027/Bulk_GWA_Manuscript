library(tidyverse)
library(sommer)
library(glue)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# ---------------------------------------------------------------------------
# Narrow-sense (additive) and broad-sense (additive+epistatic) heritability
# of the (t2 RNAi condition - t2 ctrl) phenotype for every RNAi gene x
# sample_type in the pooled_RNAi_expt dataset, using the same methodology as
# scripts/pos1_heritability.R (top-level scripts/): one value per strain,
# sommer::A.mat()/E.mat() relationship matrices (rescaled to mean(diag)==1,
# since C. elegans genotypes here are haploid-coded 0/1 and sommer's
# defaults assume diploid dosage), delta-method SEs (sommer::pin() isn't
# available in the installed sommer version; vpredict() only supports mmes
# objects, not mmer).
#   value ~ 1 + vs(strain1, Gu=A) + vs(strain2, Gu=E) + residual
#   narrow-sense h2 = V1 / (V1+V2+V3)        (additive only)
#   broad-sense  H2 = (V1+V2) / (V1+V2+V3)   (additive + epistatic; no
#                                              dominance term -- homozygous)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# phenotype: mean (RNAi - ctrl) t2 frequency per strain x Gene x sample_type,
# reconstructed with the same corrections established in
# scripts/20260813_delta_freq_correlations.R (rde-3/S28 mislabeling fix,
# S19/S20 exclusion, pooled ctrl baseline, no t0 correction) -- see that
# script for the full rationale behind each step.
# ---------------------------------------------------------------------------

load("data/processed_strain_fq.rda")

fq <- processed_strain_frq %>%
  dplyr::filter(!is.na(strain))

rde3_s28_fix <- fq %>%
  dplyr::filter(sample == "S28") %>%
  dplyr::mutate(
    `SAMPLE NAME`   = "43",
    `SAMPLE NUMBER` = 28,
    BACTERIA_PELLET = 10,
    sample_type     = "supernatant",
    Gene            = "rde-3",
    Timepoint       = 2
  )

fq <- fq %>%
  dplyr::filter(!(sample == "S29" & Gene == "rde-3")) %>%
  dplyr::bind_rows(rde3_s28_fix)

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
  dplyr::filter(n_reps == 2) # drops mig-6/supernatant (S48, failed deconvolution)

pheno_strains <- sort(unique(rnai_vs_ctrl_mean$strain))

# ---------------------------------------------------------------------------
# genotype matrix: LD-pruned (plink --indep-pairwise 50 10 0.8), MAF > 5%
# CeNDR SNP set, built for this experiment's 93-strain panel specifically
# (data/genotypes/ld_pruned_pooled/, --keep restricted to these 93 strains --
# a fresh prune, not reusing the pos1 panel's marker list, since which SNPs
# are informative/pass MAF depends on which strains are in the panel).
# ---------------------------------------------------------------------------

load("../genotypes/processed_genotype_matrix.Rda") # g, genotype_calls_vcf

pruned_ids <- readLines("../genotypes/ld_pruned_pooled/pooled_ld_pruned_markers.txt")
g_chrom_pos <- sub("_[^_]+$", "", rownames(g))

f_strains <- paste(pheno_strains, pheno_strains, sep = "_")
g_pruned <- g[g_chrom_pos %in% pruned_ids, f_strains]

M <- t(g_pruned)
rownames(M) <- pheno_strains
n_markers <- ncol(M)
n_strains <- nrow(M)

A <- sommer::A.mat(M)
A <- A / mean(diag(A))

E <- sommer::E.mat(M)
E <- E / mean(diag(E))

# ---------------------------------------------------------------------------
# fit the additive+epistatic model separately for each Gene x sample_type
# condition (A and E are built once above and reused for every fit)
# ---------------------------------------------------------------------------

delta_se_ratio3 <- function(grad, Sigma) {
  as.numeric(sqrt(t(grad) %*% Sigma %*% grad))
}

fit_heritability <- function(gene, prep) {
  cond_df <- rnai_vs_ctrl_mean %>%
    dplyr::filter(Gene == gene, sample_type == prep) %>%
    dplyr::transmute(strain1 = strain, strain2 = strain, value = mean_rnai_minus_ctrl)

  if (nrow(cond_df) < 10) {
    return(tibble::tibble(
      Gene = gene, sample_type = prep,
      metric = c("narrow_sense_h2_additive", "broad_sense_H2_additive_plus_epistatic"),
      estimate = NA_real_, se = NA_real_,
      V1_additive = NA_real_, V2_epistatic = NA_real_, V3_residual = NA_real_,
      n_strains = nrow(cond_df), n_markers = n_markers,
      note = "skipped - too few strains with usable data"
    ))
  }

  keep_strains <- cond_df$strain1
  m <- sommer::mmer(
    value ~ 1,
    random = ~ sommer::vs(strain1, Gu = A[keep_strains, keep_strains]) +
      sommer::vs(strain2, Gu = E[keep_strains, keep_strains]),
    data = cond_df, verbose = FALSE
  )

  vc <- summary(m)$varcomp
  V1 <- vc$VarComp[1]
  V2 <- vc$VarComp[2]
  V3 <- vc$VarComp[3]
  Sigma <- m$sigmaSE

  h2_narrow <- V1 / (V1 + V2 + V3)
  grad_narrow <- c((V2 + V3) / (V1 + V2 + V3)^2, -V1 / (V1 + V2 + V3)^2, -V1 / (V1 + V2 + V3)^2)
  h2_se <- delta_se_ratio3(grad_narrow, Sigma)

  H2_broad <- (V1 + V2) / (V1 + V2 + V3)
  grad_broad <- c(V3 / (V1 + V2 + V3)^2, V3 / (V1 + V2 + V3)^2, -(V1 + V2) / (V1 + V2 + V3)^2)
  H2_se <- delta_se_ratio3(grad_broad, Sigma)

  tibble::tibble(
    Gene = gene, sample_type = prep,
    metric = c("narrow_sense_h2_additive", "broad_sense_H2_additive_plus_epistatic"),
    estimate = c(h2_narrow, H2_broad), se = c(h2_se, H2_se),
    V1_additive = V1, V2_epistatic = V2, V3_residual = V3,
    n_strains = nrow(cond_df), n_markers = n_markers, note = NA_character_
  )
}

conditions <- rnai_vs_ctrl_mean %>% dplyr::distinct(Gene, sample_type)

heritability_summary <- purrr::map2_dfr(conditions$Gene, conditions$sample_type, fit_heritability)

print(heritability_summary, n = Inf)

dir.create("plot/cor_plots", recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(heritability_summary, "plot/cor_plots/pooled_heritability_summary.tsv")
