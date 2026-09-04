## Build the deposit files for the two validation experiments ----------------
##
##   Rscript scripts/make_experiments_deposit.R
##
## Reads the raw experiment folders under data/experiments/ (Dryad-hosted, not
## in the repository) and writes the compact derived files the two earliest
## supplemental figures read. Companion to make_supplemental_data.R, which
## stages straight copies; these need transformation, so they live here.
##
## The reported-r-squared table is NOT built here -- it needs PDF text
## extraction, which R cannot do conveniently. Run
##   python3 scripts/extract_sim_reported_r2.py
## which writes deconvolution/simulation_reported_r2.tsv from the archived
## per-trait PDFs and explains why that is the only route to those numbers.
##
## Outputs, all under supplemental_data/deconvolution:
##   dilution_strain_sets.tsv                 strain -> isotype -> set
##   dilution_predictions_poolref.tsv.gz      170-strain reference, 19 samples
##   dilution_predictions_fullref.tsv.gz      540-strain CeNDR reference
##   dilution_predictions_regenotype.tsv.gz   540-strain, regenotyped VCF
##   dilution_predictions_bcref.tsv.gz        84-strain B+C reference (original)
##   simulation_nnls_frequencies.tsv.gz       7 traits x 8 depths x 327 strains
##   dilution_strain_similarity.tsv           per-strain relatedness, 170 rows
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse); library(data.table)
})

BC  <- "data/experiments/20211111_BulkGWA"
SIM <- "data/experiments/initial_sims"
OUT <- "supplemental_data/deconvolution"
ISO <- "data/pooled_RNAi_expt/meta_files/strain_isotype_lookup.tsv"

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
have <- function(p) if (!file.exists(p)) { msg("NOTE: ", p, " absent -- skipped"); FALSE } else TRUE

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## --- 1. strain -> set map --------------------------------------------------
## Kept at STRAIN level on purpose. JU1580 and JU1793 are the same isotype in
## different sets, and flattening to isotype here would hide that; the figure
## script resolves it explicitly and asserts that nothing else is ambiguous.
if (have(file.path(BC, "BC_strains.txt")) && have(ISO)) {
  iso <- fread(ISO)
  sets <- fread(file.path(BC, "BC_strains.txt"),
                col.names = c("strain", "set"), header = FALSE) %>%
    left_join(iso %>% select(strain, isotype), by = "strain") %>%
    transmute(strain, isotype, set) %>% arrange(set, strain)
  write_tsv(sets, file.path(OUT, "dilution_strain_sets.tsv"))
  amb <- sets %>% filter(!is.na(isotype)) %>% distinct(isotype, set) %>%
    count(isotype) %>% filter(n > 1) %>% pull(isotype)
  msg("dilution_strain_sets.tsv: ", nrow(sets), " strains, ",
      sum(is.na(sets$isotype)), " without an isotype (",
      paste(sets$strain[is.na(sets$isotype)], collapse = ", "), ")")
  msg("  isotypes spanning two sets: ",
      if (length(amb)) paste(amb, collapse = ", ") else "none")
}

## --- 2. the four prediction tables, wide -> long ---------------------------
specs <- tribble(
  ~src,                                                   ~dest,
  "Strain_predictions_genoSubset_phenotyped_strains.tsv", "dilution_predictions_poolref.tsv.gz",
  "Strain_predictions.tsv",                               "dilution_predictions_fullref.tsv.gz",
  "Strain_predictions.regenotype.tsv",                    "dilution_predictions_regenotype.tsv.gz")
for (i in seq_len(nrow(specs))) {
  p <- file.path(BC, specs$src[i])
  if (!have(p)) next
  d <- read_tsv(p, show_col_types = FALSE) %>%
    pivot_longer(-c(strain, set), names_to = "sample", values_to = "frequency") %>%
    arrange(sample, set, strain)
  write_tsv(d, file.path(OUT, specs$dest[i]))
  msg(specs$dest[i], ": ", nrow(d), " rows, ", n_distinct(d$strain),
      " strains x ", n_distinct(d$sample), " samples")
}

## --- 3. the original B+C-reference predictions -----------------------------
p <- file.path(BC, "BC_flipCommon_NArm_predictions.rda")
if (have(p)) {
  e <- new.env(); load(p, envir = e)
  d <- e$predictions_df %>% transmute(strain, sample, frequency = new_frq) %>%
    arrange(sample, strain)
  write_tsv(d, file.path(OUT, "dilution_predictions_bcref.tsv.gz"))
  msg("dilution_predictions_bcref.tsv.gz: ", nrow(d), " rows, ",
      n_distinct(d$strain), " strains x ", n_distinct(d$sample), " samples")
}

## --- 4. simulation coefficients -> frequencies -----------------------------
## The stored coefficients are on an expected-count scale: each depth's row
## sums to the depth itself. Dividing by the row sum gives frequencies, and the
## assertion below is what establishes that reading of the scale.
p <- file.path(SIM, "plots/nnls_estFreqs.RDS")
if (have(p)) {
  x <- readRDS(p)
  long <- imap_dfr(x, function(m, trait) {
    as.data.frame(m) %>% rownames_to_column("depth") %>%
      pivot_longer(-depth, names_to = "strain", values_to = "coefficient") %>%
      mutate(trait = trait, depth = as.integer(depth), strain = sub("_.*", "", strain))
  })
  chk <- long %>% group_by(trait, depth) %>%
    summarise(s = sum(coefficient), .groups = "drop") %>%
    mutate(rel = abs(s - depth) / depth)
  stopifnot(max(chk$rel) < 0.01)
  msg("simulation coefficient rows sum to their depth (max relative error ",
      signif(max(chk$rel), 2), ")")
  long <- long %>% group_by(trait, depth) %>%
    mutate(frequency = coefficient / sum(coefficient)) %>% ungroup() %>%
    select(trait, depth, strain, coefficient, frequency) %>%
    arrange(trait, desc(depth), strain)
  write_tsv(long, file.path(OUT, "simulation_nnls_frequencies.tsv.gz"))
  msg("simulation_nnls_frequencies.tsv.gz: ", nrow(long), " rows, ",
      n_distinct(long$trait), " traits x ", n_distinct(long$depth),
      " depths x ", n_distinct(long$strain), " strains")
  msg("  negative coefficients: ", sum(long$coefficient < 0), " (",
      sprintf("%.1f%%", 100 * mean(long$coefficient < 0)), ") at depths ",
      paste(sort(unique(long$depth[long$coefficient < 0])), collapse = ", "))
}

## --- 5. the GWAS input built from the simulated frequencies ----------------
p <- file.path(SIM, "processed_simFreq_traits.tsv")
if (have(p)) {
  d <- read_tsv(p, show_col_types = FALSE)
  write_tsv(d, file.path(OUT, "simulation_gwas_traits.tsv.gz"))
  msg("simulation_gwas_traits.tsv.gz: ", nrow(d), " strains x ",
      ncol(d) - 1, " trait-depth columns")
}

## --- 6. genetic similarity among the pooled strains ------------------------
## Panel D of the dilution figure asks whether the strains that resolve badly
## are the genetically similar ones, which needs a relatedness measure. The
## genotype matrix is 12 GB in memory and Dryad-hosted, so the answer is
## reduced to a 170-row table here and that table is what the figure reads --
## the deposit-only rule holds.
##
## IBS is the fraction of markers at which two strains carry the same allele.
## With 0/1 coding, mean(x == y) = (both-1 + both-0) / n, which is what the
## crossprod expression below computes. 500,000 markers are sampled under a
## fixed seed: relatedness at this scale is estimated to far better than the
## third decimal from that many sites, and using all 2.9M costs minutes for no
## change in the ranking that panel D depends on.
GM <- "data/genotypes/processed_genotype_matrix.Rda"
SETF <- file.path(OUT, "dilution_strain_sets.tsv")
if (have(GM) && file.exists(SETF)) {
  strains_in <- read_tsv(file.path(OUT, "dilution_predictions_poolref.tsv.gz"),
                         show_col_types = FALSE) %>% pull(strain) %>% unique()
  e <- new.env(); load(GM, envir = e)
  keep <- colnames(e$g) %in% paste(strains_in, strains_in, sep = "_")
  msg("similarity: ", sum(keep), " of ", length(strains_in),
      " pooled strains found in the genotype matrix")
  set.seed(1)
  idx <- sort(sample(nrow(e$g), 5e5))
  gs <- e$g[idx, keep, drop = FALSE]
  rm(e); invisible(gc())
  colnames(gs) <- sub("_.*", "", colnames(gs))
  ok <- complete.cases(gs) & apply(gs, 1, function(x) length(unique(x)) > 1)
  gs <- gs[ok, , drop = FALSE]
  n <- nrow(gs)
  cp <- crossprod(gs); cs <- colSums(gs)
  same <- (cp + (n - outer(cs, cs, "+") + cp)) / n
  diag(same) <- NA
  sim <- tibble(strain = colnames(gs),
                nn_ibs = apply(same, 1, max, na.rm = TRUE),
                nn_partner = colnames(gs)[apply(same, 1, which.max)],
                mean_ibs = rowMeans(same, na.rm = TRUE))
  write_tsv(sim, file.path(OUT, "dilution_strain_similarity.tsv"))
  msg("dilution_strain_similarity.tsv: ", nrow(sim), " strains from ", n,
      " variable markers | nearest-neighbour IBS ",
      sprintf("%.4f-%.4f", min(sim$nn_ibs), max(sim$nn_ibs)))
}

f <- list.files(OUT, pattern = "^(dilution|simulation)", full.names = TRUE)
msg("total: ", round(sum(file.size(f)) / 1024, 1), " kB in ", length(f), " files")
