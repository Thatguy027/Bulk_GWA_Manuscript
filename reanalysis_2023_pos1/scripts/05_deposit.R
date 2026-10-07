## 2023 pos-1: write the corrected analysis into supplemental_data/ ----------
##
##   Rscript reanalysis_2023_pos1/scripts/05_deposit.R
##
## Replaces the four deposited files the 2023 pos-1 experiment feeds, keeping
## every schema and column name byte-compatible so the sixteen scripts that
## read them do not need editing:
##
##   phenotypes/pos1_2023_sample_frequencies.csv.gz
##   phenotypes/pos1_2023_association_traits.csv.gz
##   mapping/pos1_2023_gemma_loco.csv.gz
##   mapping/eigen_independent_tests.tsv   (the pos1_2023 rows only)
##
## The previous contents are in git history; reanalysis_2023_pos1/ keeps the
## provenance and the old-versus-new comparison.
##
## WHAT CHANGES. The deconvolution reference goes from 367 columns -- the whole
## strain panel -- to the 224 isotypes actually in the pool, so 142 strains that
## were never in the tube stop taking mass from the ones that were. The trait
## formulae are unchanged and reproduce the old ones to machine precision when
## fed the old frequencies; see 02_build_traits.R.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(data.table); library(dplyr)})

DAT <- "reanalysis_2023_pos1/data"
MAP <- "reanalysis_2023_pos1/mapping"
PH  <- "supplemental_data/phenotypes"
MP  <- "supplemental_data/mapping"
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

freq <- fread(cmd = paste("gzcat", shQuote(file.path(DAT, "pool_reference_frequencies.csv.gz"))))
tr   <- fread(file.path(DAT, "pool_reference_traits_dp5.csv"))

## --- 1. sample-level frequencies, in the deposited schema ------------------
## The deposited file carries T2 rows only, with the T0 frequency as a column.
c_ps <- freq[, min(frq[frq > 0]), by = depth_cutoff]
setnames(c_ps, "V1", "ps"); c_ps[, ps := ps / 2]
t0   <- freq[grepl("^T0", sample_info), .(strain, depth_cutoff, t0_frq = frq)]
ctrl <- freq[grepl("^T2_ctrl", sample_info),
             .(ctrl_frq = mean(frq)), by = .(strain, depth_cutoff)]
out <- freq[grepl("^T2_", sample_info)][t0, on = .(strain, depth_cutoff)][
  ctrl, on = .(strain, depth_cutoff)][c_ps, on = "depth_cutoff"]
out[, c("time", "rnai", "replicate") := tstrsplit(sample_info, "_", fixed = TRUE)]
out[, `:=`(delta_t0 = frq - t0_frq, delta_ctrl = frq - ctrl_frq,
           log2fc_t0 = log2((frq + ps) / (t0_frq + ps)),
           log2fc_ctrl = log2((frq + ps) / (ctrl_frq + ps)))]
out <- out[, .(sample, sample_info, time, rnai, replicate, depth_cutoff, strain,
               frq, t0_frq, ctrl_frq, delta_t0, delta_ctrl, log2fc_t0, log2fc_ctrl)]
setorder(out, depth_cutoff, sample, strain)
fwrite(out, file.path(PH, "pos1_2023_sample_frequencies.csv.gz"), na = "NA", quote = FALSE)
msg("frequencies: ", format(nrow(out), big.mark = ","), " rows | ",
    uniqueN(out$strain), " strains x ", uniqueN(out$sample), " samples x ",
    uniqueN(out$depth_cutoff), " cutoffs")

## --- 2. association traits, under the deposited column names ---------------
tr_out <- tr[, .(strain,
                 `delta_ctrl_pos-1_T2`  = delta_ctrl_pos1_T2,
                 `vst_ctrl_pos-1_T2`    = vst_ctrl_pos1_T2,
                 `log2fc_ctrl_pos-1_T2` = log2fc_ctrl_pos1_T2,
                 negctrl_growth_HT115_delta_t0,
                 negctrl_abundance_log10 = abundance_log10)]
setorder(tr_out, strain)
fwrite(tr_out, file.path(PH, "pos1_2023_association_traits.csv.gz"), na = "NA", quote = FALSE)
msg("traits: ", nrow(tr_out), " strains | with a phenotype: ",
    sum(!is.na(tr_out$`vst_ctrl_pos-1_T2`)))

## --- 3. the scan ------------------------------------------------------------
gw <- fread(cmd = paste("gzcat", shQuote(file.path(MAP, "vst_ctrl_pos1_T2_loco_results.csv.gz"))))
gw[, trait := "vst_ctrl_pos-1_T2"]          # the deposited spelling
fwrite(gw, file.path(MP, "pos1_2023_gemma_loco.csv.gz"), na = "NA", quote = FALSE)
msg("scan: ", format(nrow(gw), big.mark = ","), " markers")

## --- 4. the eigenvalue table, pos1_2023 rows only --------------------------
## gwas_thresholds() reads the scope == "total" row, so both the per-chromosome
## rows and the total are rewritten for this panel. The threshold columns are
## filled -- SUPP_FIG_XX_gwas_peak_genotype_splits.R reads thr_bonferroni and
## thr_eigen_liji straight out of this table rather than deriving them. Columns
## the recomputation does not produce (the imputation fraction, the rank, the
## trace error, and the var995 variant) are left NA rather than carried over
## from the old panel, which would attach 231-strain numbers to a 184-strain
## row.
eig <- fread(file.path(DAT, "eigen_independent_tests_pool.tsv"))
tested <- gw[, .(n_tested = .N), by = .(chrom = chr)]
eig <- eig[tested, on = "chrom"]
thr <- function(m) -log10(0.05 / m)
new_rows <- rbind(
  eig[, .(panel = "pos1_2023", chrom, n_strain, n_tested, n_marker,
          frac_calls_imputed = NA_real_, rank_nonzero = NA_integer_,
          M_eff_liji, M_eff_var995 = NA_integer_, trace_err = NA_real_,
          scope = "chromosome", thr_bonferroni = NA_real_,
          thr_eigen_liji = NA_real_, thr_eigen_var995 = NA_real_)],
  data.table(panel = "pos1_2023", chrom = "all",
             n_strain = unique(eig$n_strain), n_tested = sum(eig$n_tested),
             n_marker = sum(eig$n_marker), frac_calls_imputed = NA_real_,
             rank_nonzero = NA_integer_, M_eff_liji = sum(eig$M_eff_liji),
             M_eff_var995 = NA_integer_, trace_err = NA_real_, scope = "total",
             thr_bonferroni = thr(sum(eig$n_tested)),
             thr_eigen_liji = thr(sum(eig$M_eff_liji)),
             thr_eigen_var995 = NA_real_))
old <- fread(file.path(MP, "eigen_independent_tests.tsv"))
stopifnot(identical(names(new_rows), names(old)))
fwrite(rbind(old[panel != "pos1_2023"], new_rows),
       file.path(MP, "eigen_independent_tests.tsv"), sep = "\t", na = "NA", quote = FALSE)
msg("eigen table: pos1_2023 rows replaced | n_strain ", unique(eig$n_strain),
    " | M_eff ", round(sum(eig$M_eff_liji)), " | markers ", sum(eig$n_tested))

source("scripts/gwas_thresholds.R")
th <- gwas_thresholds("pos1_2023")
msg("thresholds now: Bonferroni ", round(th$bonferroni, 2),
    " | eigen ", round(th$eigen, 2), " (M_eff ", round(th$m_eff), ")")
