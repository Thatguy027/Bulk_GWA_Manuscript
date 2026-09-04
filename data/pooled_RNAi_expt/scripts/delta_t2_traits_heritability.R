library(tidyverse)
library(sommer)
library(glue)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# ---------------------------------------------------------------------------
# Narrow-sense (additive) and broad-sense (additive+epistatic) heritability
# of the traits in data/mapped_traits/delta_t2_traits_fix.tsv (output of
# scripts/20240328_calculate_frequencies.R's "t2 - t2c" section: each RNAi
# gene condition's t2 frequency minus the mean CTRL t2 frequency, per strain
# -- ctrl-baselined, unlike delta_t2pos_traits_fix.tsv which is pos-1-
# baselined and excludes pos-1 itself. This file includes pos-1 as a trait).
#
# Same approach as scripts/mapped_traits_heritability.R: replicate columns
# (e.g. "pellet_fog-2_S10" and "pellet_fog-2_S31") are averaged per strain
# per condition before fitting, one value per strain via sommer's A.mat/
# E.mat (rescaled to mean(diag)==1) -- see scripts/pos1_heritability.R for
# the full methodology rationale.
#
# NOTE: supernatant_ric-3 again has THREE replicate columns (S20, S41, S55)
# -- same S20-inclusion caveat as delta_t2pos_traits_fix.tsv (S20 is
# excluded elsewhere in this repo as a likely-bad original library; not
# excluded here, averaged as provided).
#
# NOTE: this trait file's 86-strain panel is NOT the same 86 strains as
# delta_t2pos_traits_fix.tsv's (the upstream nzero<30 QC filter is applied
# per-trait-set, and this file additionally includes pos-1 as a trait) -- a
# separate LD-pruned marker set was built for this panel
# (data/genotypes/ld_pruned_delta_t2/), not reused from the mapped_traits run.
#
#   value ~ 1 + vs(strain1, Gu=A) + vs(strain2, Gu=E) + residual
#   narrow-sense h2 = V1 / (V1+V2+V3)        (additive only)
#   broad-sense  H2 = (V1+V2) / (V1+V2+V3)   (additive + epistatic; no
#                                              dominance term -- homozygous)
# ---------------------------------------------------------------------------

traits_raw <- readr::read_tsv("data/mapped_traits/delta_t2_traits_fix.tsv", show_col_types = FALSE)
pheno_strains <- sort(traits_raw$strain)
traits_raw <- traits_raw %>% dplyr::arrange(strain)

trait_cols <- setdiff(colnames(traits_raw), "strain")

col_groups <- tibble::tibble(col = trait_cols) %>%
  tidyr::separate(col, into = c("sample_type", "Gene", "sample"), sep = "_", remove = FALSE, extra = "merge") %>%
  dplyr::mutate(condition = paste(sample_type, Gene, sep = "_"))

# mean across replicate columns per condition -> one column per condition
condition_means <- purrr::map_dfc(
  split(col_groups$col, col_groups$condition),
  ~ rowMeans(traits_raw[, .x, drop = FALSE])
)
colnames(condition_means) <- names(split(col_groups$col, col_groups$condition))

pheno_wide <- dplyr::bind_cols(
  tibble::tibble(strain = traits_raw$strain),
  condition_means
)

pheno_long <- pheno_wide %>%
  tidyr::pivot_longer(-strain, names_to = "trait", values_to = "value")

# ---------------------------------------------------------------------------
# genotype matrix: LD-pruned (plink --indep-pairwise 50 10 0.8), MAF > 5%
# CeNDR SNP set built for this trait file's 86-strain panel specifically
# (data/genotypes/ld_pruned_delta_t2/, --keep restricted to these 86 strains).
# ---------------------------------------------------------------------------

load("../genotypes/processed_genotype_matrix.Rda") # g, genotype_calls_vcf

pruned_ids <- readLines("../genotypes/ld_pruned_delta_t2/delta_t2_ld_pruned_markers.txt")
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
# fit the additive+epistatic model separately for each trait (condition
# means); A and E are built once above and reused for every fit
# ---------------------------------------------------------------------------

delta_se_ratio3 <- function(grad, Sigma) {
  as.numeric(sqrt(t(grad) %*% Sigma %*% grad))
}

fit_heritability <- function(trait_name) {
  cond_df <- pheno_long %>%
    dplyr::filter(trait == trait_name, !is.na(value)) %>%
    dplyr::transmute(strain1 = strain, strain2 = strain, value)

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
    trait = trait_name,
    metric = c("narrow_sense_h2_additive", "broad_sense_H2_additive_plus_epistatic"),
    estimate = c(h2_narrow, H2_broad), se = c(h2_se, H2_se),
    V1_additive = V1, V2_epistatic = V2, V3_residual = V3,
    n_strains = nrow(cond_df), n_markers = n_markers
  )
}

heritability_summary <- purrr::map_dfr(unique(pheno_long$trait), fit_heritability)

print(heritability_summary, n = Inf)

readr::write_tsv(heritability_summary, "data/mapped_traits/delta_t2_traits_heritability_summary.tsv")
