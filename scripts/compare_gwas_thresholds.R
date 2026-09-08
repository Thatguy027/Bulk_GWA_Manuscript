## Two thresholds, side by side ----------------------------------------------
##
##   Rscript scripts/gwas_qtl_intervals.R                                # eigen
##   Rscript scripts/gwas_qtl_intervals.R --threshold <perm> --label perm
##   Rscript scripts/compare_gwas_thresholds.R eigen perm
##
##     -> plots/diagnostics/TABLE_threshold_comparison.tsv
##
## Answers the only question that matters when a threshold changes: which loci
## survive, which appear, which vanish, and does the answer for chromosome III
## change. It does NOT re-derive intervals -- it reads the two tables the
## interval script wrote, so the comparison cannot drift from what was computed.
##
## MATCHING. Loci from the two runs are the same locus if their peak markers
## are within MATCH_KB. A lower threshold admits more markers, which can merge
## two loci into one or shift a peak, so an exact position match would report
## spurious appearances and disappearances.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })

OUT <- "plots/diagnostics"
MATCH_KB <- 500

a_lab <- commandArgs(TRUE)[1]; b_lab <- commandArgs(TRUE)[2]
if (is.na(a_lab) || is.na(b_lab))
  stop("usage: compare_gwas_thresholds.R <label_a> <label_b>", call. = FALSE)
f <- function(l) file.path(OUT, sprintf("TABLE_gwas_qtl_intervals_%s.tsv", l))
for (l in c(a_lab, b_lab))
  if (!file.exists(f(l)))
    stop("missing ", f(l), "\n  run: Rscript scripts/gwas_qtl_intervals.R --label ", l,
         call. = FALSE)

A <- fread(f(a_lab)); B <- fread(f(b_lab))
cat(sprintf("%-6s threshold %.3f | %d loci | %d localising\n",
            a_lab, A$threshold[1], nrow(A), sum(A$localises)))
cat(sprintf("%-6s threshold %.3f | %d loci | %d localising\n",
            b_lab, B$threshold[1], nrow(B), sum(B$localises)))
if (A$threshold[1] == B$threshold[1])
  cat("  NOTE: the two runs used the SAME threshold; this is a no-op comparison\n")

## nearest-peak matching within MATCH_KB, one-to-one, closest first
A[, aid := .I]; B[, bid := .I]
pairs <- CJ(aid = A$aid, bid = B$bid)[
  A[, .(aid, chr, ap = peak_ps)], on = "aid"][
  B[, .(bid, bchr = chr, bp = peak_ps)], on = "bid"][
  chr == bchr][, d := abs(ap - bp)][d <= MATCH_KB * 1e3][order(d)]
used_a <- integer(0); used_b <- integer(0); keep <- integer(0)
for (i in seq_len(nrow(pairs))) {
  if (pairs$aid[i] %in% used_a || pairs$bid[i] %in% used_b) next
  used_a <- c(used_a, pairs$aid[i]); used_b <- c(used_b, pairs$bid[i]); keep <- c(keep, i)
}
m <- pairs[keep]

cat(sprintf("\nmatched within %.0f kb: %d | only in %s: %d | only in %s: %d\n",
            MATCH_KB, nrow(m), a_lab, nrow(A) - nrow(m), b_lab, nrow(B) - nrow(m)))

res <- rbindlist(list(
  if (nrow(m)) A[m$aid][, .(chr, peak_Mb = round(peak_ps/1e6, 3),
      status = "both", lp = round(peak_lp, 2),
      a_ld80_kb = ld80_kb, b_ld80_kb = B[m$bid]$ld80_kb,
      a_localises = localises, b_localises = B[m$bid]$localises)],
  A[!aid %in% m$aid][, .(chr, peak_Mb = round(peak_ps/1e6, 3),
      status = paste0("only_", a_lab), lp = round(peak_lp, 2),
      a_ld80_kb = ld80_kb, b_ld80_kb = NA_integer_,
      a_localises = localises, b_localises = NA)],
  B[!bid %in% m$bid][, .(chr, peak_Mb = round(peak_ps/1e6, 3),
      status = paste0("only_", b_lab), lp = round(peak_lp, 2),
      a_ld80_kb = NA_integer_, b_ld80_kb = ld80_kb,
      a_localises = NA, b_localises = localises)]), use.names = TRUE)
setorder(res, chr, peak_Mb)
fwrite(res, file.path(OUT, "TABLE_threshold_comparison.tsv"), sep = "\t")
cat("\n== locus by locus ==\n"); print(as.data.frame(res), row.names = FALSE)

## the question the manuscript actually asks
NIL <- c(13.658, 13.695)
cat("\n== chromosome III, under each threshold ==\n")
for (nm in c(a_lab, b_lab)) {
  d <- if (nm == a_lab) A else B
  c3 <- d[chr == "III"]
  cat(sprintf("  %-6s %d locus/loci on III", nm, nrow(c3)))
  if (!nrow(c3)) { cat(" -- nothing admitted\n"); next }
  cat("\n")
  for (i in seq_len(nrow(c3))) {
    ov80 <- !(c3$ld80_hi[i] < NIL[1] | c3$ld80_lo[i] > NIL[2])
    cat(sprintf("           peak %.3f Mb (lp %.2f), r2>=0.8 interval %.3f-%.3f (%d kb) -> %s\n",
        c3$peak_ps[i]/1e6, c3$peak_lp[i], c3$ld80_lo[i], c3$ld80_hi[i], c3$ld80_kb[i],
        ifelse(ov80, "overlaps the NIL interval",
               sprintf("%.2f Mb from the NIL interval",
                       min(abs(c3$ld80_lo[i] - NIL[2]), abs(NIL[1] - c3$ld80_hi[i]))))))
  }
}
cat("\nwrote ", file.path(OUT, "TABLE_threshold_comparison.tsv"), "\n", sep = "")
