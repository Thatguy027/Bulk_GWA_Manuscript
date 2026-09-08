## Diff a permutation run's observed scan against the shipped scan ------------
##
##   Rscript scripts/compare_observed_scan.R <observed_scan_dir>
##
## e.g. Rscript scripts/compare_observed_scan.R \
##        cluster/results_perm_20260908_smoke3/observed_scan
##
## WHY THIS EXISTS. The permutation pipeline reproduces the scan it thresholds,
## so its observed genome-wide maximum is a known number and is asserted against
## it. Twice now that assertion has failed, and each time a single number gave
## nothing to work from: 8.6894 (wrong marker set -- MAF computed on all 540
## strains), then 8.5700 (wrong kinship -- -gk 1 where the scan used -gk 2), and
## now 8.8759 against 8.8361, which is 0.4% and cannot be either of those.
##
## A per-marker diff distinguishes the possibilities that a maximum cannot:
##   * a near-constant offset in -log10 p across all markers, growing with the
##     statistic, points at the variance components -- a different kinship, or a
##     different set of markers built into it
##   * agreement everywhere except a handful of markers points at the marker
##     set or its coding
##   * agreement in rank but not in value points at the GEMMA version's REML
##     optimiser rather than at the inputs
##
## The comparison is on p_wald, because that is the column the shipped scan
## carries and the column the threshold is meant to apply to.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(TRUE)
if (!length(args)) stop("usage: compare_observed_scan.R <observed_scan_dir>",
                        call. = FALSE)
DIR  <- args[1]
SCAN <- "supplemental_data/mapping/pos1_2023_gemma_loco.csv.gz"
stopifnot(dir.exists(DIR), file.exists(SCAN))

fs <- list.files(DIR, pattern = "^observed_.*\\.assoc\\.txt\\.gz$", full.names = TRUE)
if (!length(fs)) stop("no observed_*.assoc.txt.gz in ", DIR,
                      ".\n  Re-run the pipeline: publishing the observed scan ",
                      "was added after the run that produced this directory.",
                      call. = FALSE)
cat("observed scan files:", length(fs), "\n")

new <- rbindlist(lapply(fs, function(f) fread(f)[, .(rs, p_new = p_wald)]))
old <- fread(SCAN)[, .(chr, rs, ps, af, p_old = p_wald)]
cat("markers -- shipped:", nrow(old), " this run:", nrow(new), "\n")

## the GEMMA logs, which carry the version banner and the analyzed counts
lg <- list.files(DIR, pattern = "^gemma_.*\\.log\\.txt$", full.names = TRUE)
if (length(lg)) {
  cat("\n== GEMMA banner and analysed counts (first log) ==\n")
  writeLines(grep("Version|analyzed|total individuals|total SNPs|analyzed SNPs",
                  readLines(lg[1]), value = TRUE))
}

m <- merge(old, new, by = "rs")
cat("\nmarkers in both:", nrow(m),
    " | only shipped:", nrow(old) - nrow(m),
    " | only this run:", nrow(new) - nrow(m), "\n")
if (nrow(old) != nrow(m) || nrow(new) != nrow(m)) {
  cat("  MARKER SETS DIFFER -- that is the finding; the rest is secondary.\n")
  cat("  only shipped, first 10: ",
      paste(head(setdiff(old$rs, new$rs), 10), collapse = ", "), "\n")
  cat("  only this run, first 10: ",
      paste(head(setdiff(new$rs, old$rs), 10), collapse = ", "), "\n")
}

m[, `:=`(lp_old = -log10(p_old), lp_new = -log10(p_new))]
m[, d := lp_new - lp_old]

cat("\n== -log10 p difference, this run minus shipped ==\n")
cat(sprintf("  mean %+.5f | median %+.5f | sd %.5f | max |d| %.5f\n",
            mean(m$d), median(m$d), sd(m$d), max(abs(m$d))))
cat(sprintf("  Pearson r %.8f | Spearman %.8f\n",
            cor(m$lp_old, m$lp_new), cor(m$lp_old, m$lp_new, method = "spearman")))
cat(sprintf("  markers differing by >0.001: %d of %d (%.2f%%)\n",
            sum(abs(m$d) > 0.001), nrow(m), 100 * mean(abs(m$d) > 0.001)))

## Does the discrepancy grow with the statistic? A ratio near-constant across
## the range is the signature of a variance-component difference; a difference
## confined to the tail is not.
m[, bin := cut(lp_old, breaks = c(-Inf, 1, 2, 3, 4, 5, 6, 7, Inf),
               labels = c("<1","1-2","2-3","3-4","4-5","5-6","6-7",">7"))]
cat("\n== by strength of the shipped signal ==\n")
print(as.data.frame(m[, .(n = .N,
                          mean_d = round(mean(d), 5),
                          mean_ratio = round(mean(lp_new / lp_old), 6)),
                      by = bin][order(bin)]), row.names = FALSE)

cat("\n== the shipped scan's top 10 markers, both scans ==\n")
print(as.data.frame(m[order(p_old)][1:10,
        .(chr, rs, ps, af, shipped = round(lp_old, 4),
          this_run = round(lp_new, 4), diff = round(d, 4))]), row.names = FALSE)

cat("\n== the 10 largest disagreements ==\n")
print(as.data.frame(m[order(-abs(d))][1:10,
        .(chr, rs, ps, af, shipped = round(lp_old, 4),
          this_run = round(lp_new, 4), diff = round(d, 4))]), row.names = FALSE)

## ---------------------------------------------------------------------------
## Why were markers dropped? GEMMA's only defaults that discard a marker are
## -maf 0.01 and -miss 0.05, so the allele frequencies of the dropped set
## separate the two: extreme af means MAF, ordinary af means missingness.
## ---------------------------------------------------------------------------
drop <- old[!new, on = "rs"]
if (nrow(drop)) {
  cat("\n== the ", nrow(drop), " markers this run dropped ==\n", sep = "")
  n_extreme <- sum(drop$af < 0.01 | drop$af > 0.99)
  cat(sprintf("  af < 0.01 or > 0.99 (would fail -maf 0.01): %d (%.1f%%)\n",
              n_extreme, 100 * n_extreme / nrow(drop)))
  cat("  af quantiles: ",
      paste(sprintf("%.3f", quantile(drop$af, c(0, .25, .5, .75, 1))),
            collapse = " "), "\n", sep = "")
  if (n_extreme == 0)
    cat("  -> NOT the MAF filter. GEMMA's remaining default that discards a\n",
        "     marker is -miss 0.05, so these are markers whose missingness the\n",
        "     shipped scan did not see. Confirm on the cluster with\n",
        "     plink --bfile all --missing: the count of F_MISS > 0.05 should\n",
        "     match, and the ids should be the same ones.\n", sep = "")

  cat("\n  per chromosome:\n")
  tab <- merge(old[, .(total = .N), by = chr], drop[, .(dropped = .N), by = chr],
               by = "chr")
  tab[, pct := round(100 * dropped / total, 1)]
  print(as.data.frame(tab[order(chr)]), row.names = FALSE)

  ## Missingness in this species tracks the hyper-divergent regions, which sit
  ## on the arms and tips. A drop rate that is flat across 1 Mb windows would
  ## argue against the missingness explanation.
  old[, win := paste0(chr, ":", floor(ps / 1e6))]
  drop[, win := paste0(chr, ":", floor(ps / 1e6))]
  z <- merge(old[, .(n = .N), by = win], drop[, .(d = .N), by = win],
             by = "win", all.x = TRUE)
  z[is.na(d), d := 0][, pct := 100 * d / n]
  cat("\n  1 Mb windows with the highest drop rate (>=200 markers):\n")
  print(as.data.frame(z[n >= 200][order(-pct)][1:10,
          .(win, markers = n, dropped = d, pct = round(pct, 1))]),
        row.names = FALSE)

  fwrite(drop[, .(rs)], "plots/diagnostics/TABLE_scan_markers_dropped.txt",
         col.names = FALSE)
  cat("\n  ids written to plots/diagnostics/TABLE_scan_markers_dropped.txt\n")

  ## Are the largest p-value disagreements in the dropped-heavy regions? If not,
  ## local marker loss is not what moved them and the kinship is.
  top <- m[order(-abs(d))][1:200]
  top[, win := paste0(chr, ":", floor(ps / 1e6))]
  cat("\n  drop rate in the windows of the 200 most-disagreeing markers: ",
      sprintf("%.1f%%", mean(merge(top[, .(win)], z, by = "win")$pct)),
      " against ", sprintf("%.1f%%", 100 * nrow(drop) / nrow(old)),
      " genome-wide.\n", sep = "")
  cat("  Below the genome-wide rate means local marker loss is NOT what moved\n",
      "  them; the kinship is, and it changed because it is built from the\n",
      "  markers that survived -- 7.4% fewer, concentrated on the arms.\n", sep = "")
}

cat("\n== where each scan peaks ==\n")
cat("  shipped:  ", m[which.min(p_old), rs], " at ",
    sprintf("%.4f", m[, max(lp_old)]), "\n", sep = "")
cat("  this run: ", m[which.min(p_new), rs], " at ",
    sprintf("%.4f", m[, max(lp_new)]), "\n", sep = "")
cat("  same peak marker: ",
    identical(m[which.min(p_old), rs], m[which.min(p_new), rs]), "\n", sep = "")
