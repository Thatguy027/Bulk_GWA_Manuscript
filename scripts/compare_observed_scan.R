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

cat("\n== where each scan peaks ==\n")
cat("  shipped:  ", m[which.min(p_old), rs], " at ",
    sprintf("%.4f", m[, max(lp_old)]), "\n", sep = "")
cat("  this run: ", m[which.min(p_new), rs], " at ",
    sprintf("%.4f", m[, max(lp_new)]), "\n", sep = "")
cat("  same peak marker: ",
    identical(m[which.min(p_old), rs], m[which.min(p_new), rs]), "\n", sep = "")
