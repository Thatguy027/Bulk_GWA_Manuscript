## The per-strain responsiveness test behind the 2023 pos-1 pilot ------------
##
##   Rscript scripts/DIAG_pos1_responsiveness_test.R
##     -> plots/diagnostics/pos1_responsiveness_variants.tsv
##
## A reviewer asked for a formal test of responsiveness in the 231-strain pos-1
## experiment. Two claims were drafted in reply: that 183 of 231 isolates
## declined, and that a per-strain test identifies 74 with individually
## significant declines. This recomputes both from the deposit.
##
## THE POPULATION CLAIM REPRODUCES EXACTLY. 183 of 231 strains carrying a
## variance-stabilised value have delta_ctrl_pos-1_T2 < 0, which is 79.2%, and
## the binomial sign test against 0.5 gives p < 2.22e-16. The same 183 comes out
## of the vst and log2fc columns, so it does not depend on the parameterisation.
##
## AND A FLOOR THAT OMITS THE PURGED STRAINS IS NOT CONSERVATIVE, IT IS WRONG.
## 81 of the 231 are absent from every pos-1 pool while present in control --
## the strongest response the experiment can show -- and a t-test on constant
## data returns nothing, so they fall out of every t-test construction. Counting
## them alongside the 60 that clear BH among the 150 testable strains gives a
## floor of 141 of 231, which is both defensible and higher than the 74.
##
## THE PER-STRAIN CLAIM OF 74 DOES NOT SURVIVE. It is reproducible, but only by
## pooling the three depth cutoffs, which turns 4 pos-1 libraries into 12
## observations and 2 control pools into 6. Those are not 12 libraries. They are
## the same 4 sequenced once and filtered at three read-depth thresholds, and
## the cutoff-3 and cutoff-10 values of a single library correlate at r = 0.998.
## The degrees of freedom are inflated about threefold and the test is
## anticonservative. The drafted sentence also describes a four-versus-two test,
## which is not what produces 74.
##
## What the design does support is printed below. The honest reading is that the
## experiment is well powered for a statement about the POPULATION of strains --
## two control pools and four treatment pools make the binomial decisive -- and
## poorly powered for a statement about any INDIVIDUAL strain, because a Welch
## test on 4 against 2 has almost no degrees of freedom.
##
## Reads only supplemental_data. Exploratory; nothing in the manuscript reads
## this file, and it is here so the number that does go in can be checked.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
PH   <- "supplemental_data/phenotypes"
DIAG <- "plots/diagnostics"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)

tr <- fread(cmd = paste("gzcat", shQuote(file.path(PH, "pos1_2023_association_traits.csv.gz"))))
keep <- tr[is.finite(`vst_ctrl_pos-1_T2`)]$strain
cat(sprintf("strains with a variance-stabilised value: %d of %d\n", length(keep), nrow(tr)))

## --- the population claim ---------------------------------------------------
v <- tr[strain %in% keep]
nd <- sum(v$`delta_ctrl_pos-1_T2` < 0)
bt <- binom.test(nd, length(keep), 0.5)
cat(sprintf("\nPOPULATION: %d of %d decline (%.1f%%), binomial sign test p = %s\n",
            nd, length(keep), 100 * nd / length(keep), format.pval(bt$p.value)))
cat(sprintf("  identical count from vst (%d) and log2fc (%d)\n",
            sum(v$`vst_ctrl_pos-1_T2` < 0), sum(v$`log2fc_ctrl_pos-1_T2` < 0)))

## --- the per-strain claim, every construction -------------------------------
raw <- fread(cmd = paste("gzcat", shQuote(file.path(PH, "pos1_2023_sample_frequencies.csv.gz"))))
## rows sharing a (strain, sample) key are summed: JU1793 appears twice per
## sample in this dataset, per METHODS.txt
agg <- function(x, ...) x[, .(frq = sum(frq), delta_ctrl = sum(delta_ctrl)), by = c(...)]

## returns a one-row data.table, not a list: data.table() recycles a bare list
## argument into ROWS rather than splicing it into columns
count <- function(construction, cutoff, independent, p, down) {
  q <- p.adjust(p, "BH")
  data.table(construction = construction, cutoff = cutoff, independent = independent,
             q05_down = sum(q < 0.05 & down, na.rm = TRUE),
             q05_up   = sum(q < 0.05 & !down, na.rm = TRUE),
             raw05_down = sum(p < 0.05 & down, na.rm = TRUE))
}
out <- list()
for (dc in c(3, 5, 10)) {
  d <- agg(raw[depth_cutoff == dc & strain %in% keep], "strain", "sample_info", "rnai")
  a <- d[, { t <- tryCatch(t.test(frq[rnai == "pos-1"], frq[rnai != "pos-1"]),
                           error = function(e) NULL)
             .(pv = if (is.null(t)) NA_real_ else t$p.value,
               down = mean(frq[rnai == "pos-1"]) < mean(frq[rnai != "pos-1"])) }, by = strain]
  out[[length(out) + 1]] <- count("Welch, 4 pos-1 vs 2 ctrl", as.character(dc), TRUE,
                                  a$pv, a$down)
  b <- d[rnai == "pos-1", { t <- tryCatch(t.test(delta_ctrl), error = function(e) NULL)
             .(pv = if (is.null(t)) NA_real_ else t$p.value, down = mean(delta_ctrl) < 0) }, by = strain]
  out[[length(out) + 1]] <- count("one-sample t on 4 delta_ctrl vs 0", as.character(dc),
                                  TRUE, b$pv, b$down)
}
## the construction that produces 74, reproduced so it can be seen to be wrong
d <- agg(raw[strain %in% keep], "strain", "sample_info", "rnai", "depth_cutoff")
a <- d[, { t <- tryCatch(t.test(frq[rnai == "pos-1"], frq[rnai != "pos-1"]),
                         error = function(e) NULL)
           .(pv = if (is.null(t)) NA_real_ else t$p.value,
             down = mean(frq[rnai == "pos-1"]) < mean(frq[rnai != "pos-1"])) }, by = strain]
out[[length(out) + 1]] <- count("Welch, 12 vs 6 (cutoffs pooled)", "3+5+10", FALSE,
                                a$pv, a$down)

R <- rbindlist(out)
fwrite(R, file.path(DIAG, "pos1_responsiveness_variants.tsv"), sep = "\t")
cat("\nPER STRAIN, every construction (BH across the 231 tests):\n")
print(R)

## --- the strains a t-test cannot see ---------------------------------------
## 81 of the 231 are absent from EVERY pos-1 pool while present in control. That
## is the strongest evidence of response the experiment can produce, and a
## t-test on constant data returns nothing at all, so they drop silently out of
## every construction above. Any floor that omits them is not conservative, it
## is wrong.
d5 <- agg(raw[depth_cutoff == 5 & strain %in% keep], "strain", "sample_info", "rnai")
zz <- merge(
  d5[rnai == "pos-1", .(sd_delta = sd(delta_ctrl), all_zero = all(frq == 0)), by = strain],
  d5[rnai != "pos-1", .(ctrl_frq = mean(frq)), by = strain], by = "strain")
purged <- zz[all_zero & ctrl_frq > 0, .N]
untest <- zz[sd_delta == 0 | is.na(sd_delta), .N]
sig <- R[construction == "one-sample t on 4 delta_ctrl vs 0" & cutoff == "5"]$q05_down
cat(sprintf("\nTHE STRAINS THE TEST CANNOT SEE.\n"))
cat(sprintf("  zero variance across the four pos-1 pools: %d of %d\n", untest, length(keep)))
cat(sprintf("  of those, absent from every pos-1 pool but present in control: %d\n", purged))
cat(sprintf("\n  defensible floor = %d purged + %d significant among the %d testable = %d of %d (%.0f%%)\n",
            purged, sig, length(keep) - untest, purged + sig, length(keep),
            100 * (purged + sig) / length(keep)))
fwrite(data.table(n_strains = length(keep), n_decline = nd,
                  binom_p = bt$p.value, n_purged = purged, n_untestable = untest,
                  n_sig_onesample_c5 = sig, floor_total = purged + sig),
       file.path(DIAG, "pos1_responsiveness_headline.tsv"), sep = "\t")

## --- why the last row is not a test ----------------------------------------
cc <- raw[strain %in% keep & rnai == "pos-1",
          .(r = cor(frq[depth_cutoff == 3], frq[depth_cutoff == 10])), by = sample_info]
cat("\nWHY THE 12-vs-6 ROW IS PSEUDO-REPLICATION.\n")
cat("The three cutoffs are one sequencing library filtered three ways, not three\n")
cat("libraries. Correlation between the cutoff-3 and cutoff-10 values of the SAME\n")
cat("library, per pos-1 pool:\n")
print(cc)
cat(sprintf("\nmedian r = %.4f. Treating these as independent inflates the degrees of\n",
            median(cc$r)))
cat("freedom about threefold, which is the whole difference between 6 and 74.\n")
