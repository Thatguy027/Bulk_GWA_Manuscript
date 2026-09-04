## Stage supplemental_data/ from the working tree ----------------------------
##
##   Rscript scripts/make_supplemental_data.R
##
## Copies every file the sixteen figure scripts read into supplemental_data/,
## under descriptive names, compressing what compresses. The result is
## self-contained: with only supplemental_data/ and scripts/, all sixteen
## figures rebuild.
##
## BUILD ORDER for a deposit from scratch
##   1  Rscript scripts/make_supplemental_data.R     (this script)
##   2  Rscript scripts/make_thinned_bundle.R        (needs the full bundle)
##   3  python3  scripts/sid2_ribbon_render.py       (writes sid2_per_residue.tsv)
##   4  Rscript scripts/sid2_variant_table.R         (writes the variant tables)
##   5  python3  scripts/sid2_zoom_render.py         (structure renders)
##
## Steps 3-5 write into supplemental_data/ themselves, so they are not staged
## here: there is one copy of each derived file and it lives in the deposit.
##
## THREE STAGED FILES HAVE PRODUCERS THAT NEED THE DRYAD ARCHIVE, not the
## deposit, so they are copied rather than rebuilt, and re-running their
## producer means re-running this script:
##   mapping/eigen_independent_tests.tsv    <- eigen_independent_tests.R
##   deconvolution/baugh_downsampled_slopes.rda <- baugh_L1_DownSample_Counts.R
##
## The two validation experiments are NOT staged here. Their deposit files are
## derived rather than copied, so they have their own builders:
##   deconvolution/dilution_*        <- make_experiments_deposit.R
##   deconvolution/simulation_nnls_frequencies.tsv.gz,
##   deconvolution/simulation_gwas_traits.tsv.gz
##                                   <- make_experiments_deposit.R
##   deconvolution/simulation_reported_r2.tsv
##                                   <- extract_sim_reported_r2.py  (the only
##                                      route to those values; see that script)
##   structure/sid2_parental_variants.tsv   <- pooled_cross_candidate_variation.R
##
## NOT copied here, because they have their own builders:
##   supplemental_data/mapping/pooled_cross_bundle_thinned.rds
##     <- scripts/make_thinned_bundle.R
##   supplemental_data/genotypes/sid2_region.{bed,bim,fam}
##     <- built below with plink2, subset to chrIII:13,679,000-13,682,000
##
## COMPRESSION POLICY
##   gzip: derived tables read by readr/data.table, which open .gz directly
##   leave alone: files already compressed (.rds/.rda/.RData/.csv.gz/.tsv.gz),
##     PDB files (Biopython cannot read gzip), the PLINK trio (plink needs it
##     uncompressed), tiny JSON, and the hand-scored primary data -- those stay
##     plain text so a reader can open them without tooling
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(dplyr) })

SD  <- "supplemental_data"
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## src, dest, gzip?
M <- tibble::tribble(
  ~src, ~dest, ~gz,
  # --- phenotypes -----------------------------------------------------------
  "data/pooled_RNAi_expt/reanalysis/vst_association_traits.csv",
    "phenotypes/pooled_vst_traits.csv", TRUE,
  "data/pos1_original/updated_analysis/association_traits.csv",
    "phenotypes/pos1_2023_association_traits.csv", TRUE,
  "data/pos1_original/updated_analysis/final_dataset.csv",
    "phenotypes/pos1_2023_sample_frequencies.csv", TRUE,
  "data/pos1_plate_phenotyping/pos1_phenotypes_first2rounds.tsv",
    "phenotypes/plate_scores_pos1.tsv", FALSE,
  "data/pooled_RNAi_expt/paaby2015/emb_leth_data.txt",
    "phenotypes/paaby2015_embryonic_lethality.txt", TRUE,
  # --- hatching assays (primary, hand-scored: left as plain text) -----------
  "data/plate_rnai_phenotyping/RNAi Sensitivity - JU1793 - November2025 gSZ177 - gSZ179 NILs.tsv",
    "hatching_assays/nil_series_hatching.tsv", FALSE,
  "data/plate_rnai_phenotyping/20260409_ju2466swap_plus_N2A_swap.csv",
    "hatching_assays/ju_allele_swaps_hatching.csv", FALSE,
  "data/plate_rnai_phenotyping/n2_swap.tsv",
    "hatching_assays/n2_allele_swaps_hatching.tsv", FALSE,
  "data/nil_ranges.bed",
    "hatching_assays/nil_introgression_ranges.bed", FALSE,
  # --- mapping --------------------------------------------------------------
  "data/pos1_original/updated_analysis/vst_ctrl_pos-1_T2_loco_results.csv.gz",
    "mapping/pos1_2023_gemma_loco.csv.gz", FALSE,
  "data/cross_experiments/JU1793-JU2466_export/plot_data/JU1793_JU2466_F2-2_contrast_HT115g-POS1g_10000_plot_DF.tsv.gz",
    "mapping/ju_cross_ht115_vs_pos1_scan.tsv.gz", FALSE,
  "data/eigen_independent_tests.tsv",
    "mapping/eigen_independent_tests.tsv", FALSE,
  # --- deconvolution validation --------------------------------------------
  "data/baugh/2024_processedBOOTs_with_MIP.RData",
    "deconvolution/baugh_nnls_with_mipseq.RData", FALSE,
  "data/baugh/2024baugh_bootstrap_prediction.rda",
    "deconvolution/baugh_bootstrap_array.rda", FALSE,
  "data/baugh/2026downsampled_dataset.rda",
    "deconvolution/baugh_downsampled_slopes.rda", FALSE,
  "data/baugh/MIPseq_frequencies.txt",
    "deconvolution/mipseq_frequencies.txt", TRUE,
  "data/baugh/cache_boot_slopes.rds",
    "deconvolution/cache_boot_slopes.rds", FALSE,
  "data/baugh/cache_boot_freq.rds",
    "deconvolution/cache_boot_freq.rds", FALSE,
  # --- structure ------------------------------------------------------------
  "data/structure/AF_sid2_wt_dimer/sid2wt.pdb",
    "structure/sid2_alphafold_dimer.pdb", FALSE,
  "data/structure_modeling/claude_docking/stage4_electrostatics/sid2_wt_elec.pdb",
    "structure/sid2_membrane_oriented.pdb", FALSE,
  "data/structure/SID2deeptmhmm/predicted_topologies.3line",
    "structure/sid2_deeptmhmm_topology.3line", FALSE,
  "plots/pooled_cross_intersection/TABLE_sid2_parental_variants.tsv",
    "structure/sid2_parental_variants.tsv", FALSE)

## the five AlphaFold confidence files
conf <- list.files("data/structure/AF_sid2_wt_dimer",
                   pattern = "summary_confidences_.*json$", full.names = TRUE)
M <- bind_rows(M, tibble::tibble(
  src = conf,
  dest = file.path("structure/sid2_alphafold_confidences", basename(conf)),
  gz = FALSE))

missing <- M$src[!file.exists(M$src)]
if (length(missing)) { cat("MISSING:\n"); print(missing); stop("inputs absent") }

raw <- out <- 0
for (i in seq_len(nrow(M))) {
  s <- M$src[i]; d <- file.path(SD, M$dest[i])
  dir.create(dirname(d), showWarnings = FALSE, recursive = TRUE)
  if (M$gz[i]) {
    d <- paste0(d, ".gz")
    con_in <- file(s, "rb"); con_out <- gzfile(d, "wb", compression = 9)
    repeat { x <- readBin(con_in, "raw", 1e6); if (!length(x)) break
             writeBin(x, con_out) }
    close(con_in); close(con_out)
  } else file.copy(s, d, overwrite = TRUE)
  raw <- raw + file.size(s); out <- out + file.size(d)
}
msg(nrow(M), " files copied: ", round(raw / 1e6, 1), " MB -> ",
    round(out / 1e6, 1), " MB")

## --- strain order for the bootstrap array ---------------------------------
## The array's first dimension is unnamed and its order is the column order of
## the 31 MB genotype matrix. Extract the 102 names once so the deposit does
## not have to carry that matrix.
BOOTIN <- "data/baugh/2024bootstrapINPUT.Rdata"
SORD   <- file.path(SD, "deconvolution", "baugh_strain_order.txt")
if (file.exists(BOOTIN)) {
  g <- new.env(); load(BOOTIN, g)
  strains <- sub("_.*$", "", colnames(g$flipped_bootstrap_input[[1]]))
  ## a header line, so the file is self-describing in the deposit
  writeLines(c("strain", strains), SORD)
  msg("strain order: ", length(strains), " names, ",
      round(file.size(SORD) / 1e3, 1), " kB (from ",
      round(file.size(BOOTIN) / 1e6, 1), " MB)")
  rm(g)
} else if (file.exists(SORD)) {
  msg("strain order: kept existing ", SORD)
} else msg("NOTE: ", BOOTIN, " absent -- strain order not written")

## --- genotype subset ------------------------------------------------------
## The two scripts that read genotypes want the sid-2 region only. Subsetting
## turns 56.5 MB of chromosome III into a few kB.
PL  <- "data/genotypes/CeNDR20210121_Plink/III"
GT  <- file.path(SD, "genotypes")
dir.create(GT, showWarnings = FALSE, recursive = TRUE)
if (file.exists(paste0(PL, ".bed"))) {
  st <- system2("plink2", c("--bfile", PL, "--chr", "III",
                            "--from-bp", 13679000, "--to-bp", 13682000,
                            "--make-bed", "--out", file.path(GT, "sid2_region"),
                            "--allow-extra-chr", "--silent"),
                stdout = TRUE, stderr = TRUE)
  if (!file.exists(file.path(GT, "sid2_region.bed"))) {
    cat(st, sep = "\n"); stop("plink2 subset failed")
  }
  n <- length(readLines(file.path(GT, "sid2_region.bim")))
  sz <- sum(file.size(list.files(GT, "^sid2_region", full.names = TRUE)))
  msg("genotype subset: ", n, " variants, ", round(sz / 1e3, 1), " kB (from ",
      round(sum(file.size(paste0(PL, c(".bed", ".bim", ".fam")))) / 1e6, 1),
      " MB)")
} else {
  msg("NOTE: ", PL, ".bed absent -- genotype subset not rebuilt")
}

tot <- sum(file.size(list.files(SD, recursive = TRUE, full.names = TRUE)))
msg("supplemental_data total: ", round(tot / 1e6, 2), " MB in ",
    length(list.files(SD, recursive = TRUE)), " files")
