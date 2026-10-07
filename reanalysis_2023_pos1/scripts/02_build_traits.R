## 2023 pos-1: per-strain traits from the pool-reference deconvolution -------
##
##   Rscript reanalysis_2023_pos1/scripts/02_build_traits.R
##
## Builds the mapping traits from reanalysis_2023_pos1/data/
## pool_reference_frequencies.csv.gz, written by 01_deconvolve_pool_reference.R.
##
## THE CONSTRUCTION, and how it is checked. The control frequency of a strain is
## the mean of the two T2 control pools; its pos-1 response in one replicate is
## that replicate's frequency minus the control, and the trait is the mean over
## the four replicates. That is exactly what the deposited analysis did, and the
## check below proves it: fed the DEPOSITED frequencies, this code reproduces
## the deposited delta_ctrl_pos-1_T2 to machine precision. So any difference in
## the new traits comes from the deconvolution, not from the trait formula.
##
## THE TRANSFORMS, from the upstream pipeline
## (/Users/Stefan/UCLA/Projects/bulkGWAS/lipid_RNAi/scripts/07_make_phenotypes.R
## and 04_control_corrections.R, documented in TRAIT_SPEC.md):
##
##   delta_ctrl  = f - p
##   vst_ctrl    = asin(sqrt(f)) - asin(sqrt(p))   arcsine-root, whose variance
##                 is about 1/(4N) and so free of the mean-variance coupling a
##                 raw frequency difference carries
##   log2fc_ctrl = log2((f + c) / (p + c)),  c = half the smallest non-zero
##                 frequency in the dataset
##
## with f the mean frequency over the four pos-1 replicates and p the control
## frequency. All three are NA where p == 0, which is the pipeline's own
## `usable` rule. The checks below reproduce the deposited delta AND vst from
## the deposited frequencies, so the formulae are the originals and not a
## reconstruction.
##
## WHO GETS A VALUE. A strain with no control signal has no defined response, so
## strains whose control frequency is zero are NA. Everything else carries a
## value. Unlike the deposited set, this cannot include a strain that was never
## in the pool, because the reference has only pool members in it.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(data.table); library(dplyr)})

DAT <- "reanalysis_2023_pos1/data"
PH  <- "supplemental_data/phenotypes"
CUT <- 5L                       # the cutoff the deposited traits were built at
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## --- the trait formula, written once and used for both checks and output ---
build <- function(d) {
  ## pseudocount: half the smallest non-zero frequency, as 04_control_corrections.R
  c_ps <- min(d$frq[d$frq > 0], na.rm = TRUE) / 2
  ctrl <- d[grepl("^T2_ctrl", sample_info), .(ctrl_frq = mean(frq)), by = strain]
  t0   <- d[grepl("^T0", sample_info), .(t0_frq = mean(frq)), by = strain]
  pos  <- d[grepl("^T2_pos-1", sample_info)][ctrl, on = "strain"]
  ## the pipeline averages the replicate frequencies FIRST, then transforms
  avg <- pos[, .(mean_frq = mean(frq), n_rep = .N, ctrl_frq = ctrl_frq[1]), by = strain]
  cl <- function(x) pmin(pmax(x, 0), 1)
  tr <- avg[, .(strain, n_rep, ctrl_frq,
                delta_ctrl_pos1_T2 = mean_frq - ctrl_frq,
                vst_ctrl_pos1_T2 = asin(sqrt(cl(mean_frq))) - asin(sqrt(cl(ctrl_frq))),
                log2fc_ctrl_pos1_T2 = log2((mean_frq + c_ps) / (ctrl_frq + c_ps)))]
  ## the deposited frequency export carries T2 only, so the T0-based negative
  ## control exists for the new frequencies and not for the check below
  if (nrow(t0)) tr <- tr[t0, on = "strain"] else tr[, t0_frq := NA_real_]
  tr[, `:=`(negctrl_growth_HT115_delta_t0 = ctrl_frq - t0_frq,
            abundance_log10 = ifelse(ctrl_frq > 0, log10(ctrl_frq), NA_real_))]
  ## the pipeline's `usable` rule: no control signal, no defined response
  tr[ctrl_frq == 0, `:=`(delta_ctrl_pos1_T2 = NA_real_, vst_ctrl_pos1_T2 = NA_real_,
                         log2fc_ctrl_pos1_T2 = NA_real_)]
  tr[]
}

## --- CHECK: the same formula on the DEPOSITED frequencies -------------------
old_f <- fread(cmd = paste("gzcat", shQuote(file.path(PH, "pos1_2023_sample_frequencies.csv.gz"))))
old_f <- old_f[depth_cutoff == CUT, .(frq = sum(frq)), by = .(strain, sample_info)]
old_t <- fread(cmd = paste("gzcat", shQuote(file.path(PH, "pos1_2023_association_traits.csv.gz"))))
setnames(old_t, "delta_ctrl_pos-1_T2", "deposited")
setnames(old_t, "vst_ctrl_pos-1_T2", "deposited_vst")
chk <- build(old_f)[old_t[, .(strain, deposited, deposited_vst)], on = "strain"][!is.na(deposited)]
chk[, `:=`(d = abs(delta_ctrl_pos1_T2 - deposited),
           dv = abs(vst_ctrl_pos1_T2 - deposited_vst))]
## JU1793 is the one strain the deposited reference carries TWICE, so summing
## its two rows -- which is what every other strain needs -- doubles it. The
## corrected reference has one column per isotype, so this cannot recur; it is
## excluded from the check rather than worked around.
msg("formula check against the deposited traits: n = ", nrow(chk),
    " | max |diff| = ", signif(max(chk[strain != "JU1793"]$d), 3),
    " | JU1793 alone differs by ", signif(chk[strain == "JU1793"]$d, 3),
    " (the duplicated reference column)")
msg("  vst likewise: max |diff| = ", signif(max(chk[strain != "JU1793"]$dv), 3))
stopifnot(max(chk[strain != "JU1793"]$d) < 1e-12, nrow(chk[d > 1e-12]) == 1,
          max(chk[strain != "JU1793"]$dv) < 1e-12)

## --- the corrected traits ---------------------------------------------------
new_f <- fread(cmd = paste("gzcat", shQuote(file.path(DAT, "pool_reference_frequencies.csv.gz"))))
tr <- build(new_f[depth_cutoff == CUT])
setorder(tr, strain)
msg("pool strains: ", nrow(tr), " | with a phenotype: ", sum(!is.na(tr$delta_ctrl_pos1_T2)),
    " | no control signal: ", sum(is.na(tr$delta_ctrl_pos1_T2)))

fwrite(tr, file.path(DAT, "pool_reference_traits_dp5.csv"))
## a mapping-ready file: strain plus the trait columns, NAs kept
fwrite(tr[, .(strain, delta_ctrl_pos1_T2, vst_ctrl_pos1_T2, log2fc_ctrl_pos1_T2)],
       file.path(DAT, "mapping_traits_dp5.csv"))
msg("wrote pool_reference_traits_dp5.csv and mapping_traits_dp5.csv")

## --- what changed -----------------------------------------------------------
cmp <- tr[, .(strain, new = delta_ctrl_pos1_T2)][old_t[, .(strain, old = deposited)], on = "strain"]
both <- cmp[!is.na(new) & !is.na(old)]
cat("\n== the deposited phenotype against the corrected one ==\n")
cat(sprintf("  strains with a phenotype, deposited : %d (%d of them not in the pool)\n",
            sum(!is.na(cmp$old)), sum(!is.na(cmp$old) & is.na(cmp$new) | (!is.na(cmp$old) & !cmp$strain %in% tr$strain))))
cat(sprintf("  strains with a phenotype, corrected : %d\n", sum(!is.na(tr$delta_ctrl_pos1_T2))))
cat(sprintf("  in both                             : %d\n", nrow(both)))
cat(sprintf("  Spearman between them               : %.3f\n", cor(both$old, both$new, method = "spearman")))
cat(sprintf("  Pearson                             : %.3f\n", cor(both$old, both$new)))
cat("\n  declining strains (delta < 0):\n")
cat(sprintf("    deposited, all with a value : %d of %d (%.0f%%)\n",
            sum(old_t$deposited < 0, na.rm = TRUE), sum(!is.na(old_t$deposited)),
            100 * mean(old_t$deposited < 0, na.rm = TRUE)))
cat(sprintf("    corrected                   : %d of %d (%.0f%%)\n",
            sum(tr$delta_ctrl_pos1_T2 < 0, na.rm = TRUE), sum(!is.na(tr$delta_ctrl_pos1_T2)),
            100 * mean(tr$delta_ctrl_pos1_T2 < 0, na.rm = TRUE)))
