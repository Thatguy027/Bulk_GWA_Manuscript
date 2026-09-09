## Repeatability of the pooled pos-1 response ---------------------------------
##
##   Rscript scripts/pos1_repeatability.R
##     -> plots/diagnostics/TABLE_pos1_repeatability.tsv
##
## THE QUESTION. The draft asserts a broad-sense heritability of 0.32 for the
## pooled pos-1 response. No heritability estimate exists anywhere in this
## repository, and 0.32 passes the number checker only because it rounds onto an
## unrelated leakage rho of 0.326. The 2023 experiment has four pos-1 replicate
## pools against a control, so what these data CAN support is a repeatability --
## the intraclass correlation across those replicates.
##
## REPEATABILITY IS NOT HERITABILITY, and the difference matters here rather
## than being a formality. Two reasons, both structural:
##
##   The four replicates are replicate pools within one experiment, on one
##   batch of animals, at one timepoint. They span technical and within-batch
##   variation, not the environmental range a broad-sense estimate needs. R is
##   an upper bound on H2, not an estimate of it.
##
##   All four share ONE control baseline. ctrl_frq is the mean over both
##   control replicates and is verified below to be identical across all four
##   pos-1 replicates for every strain, so any error in the baseline enters all
##   four values as a common term and is counted as among-strain signal.
##
## THE ANSWER DEPENDS ENTIRELY ON SCALE, which is the substantive finding here
## and the reason this script reports three numbers instead of one. On the raw
## delta the strains separate almost perfectly, but they separate by ABUNDANCE:
## a strain's absolute frequency change tracks how much of the pool it occupies
## (Spearman ~0.90), so the between-strain variance is mostly "how common is
## this strain", not "how does it respond". Normalising that out by taking the
## log ratio to the control roughly halves the estimate. The mapped trait is
## variance-stabilised for exactly this reason, but the VST is defined per
## strain rather than per replicate, so it cannot be used here.
##
## WHICH NUMBER TO QUOTE. The log2 fold-change one. It is the only one of the
## three on a scale where a large and a small strain contribute comparably, and
## it is therefore the only one that answers the question the sentence in the
## draft is asking. Quote it as repeatability, with the replicate structure
## stated, and do not call it heritability.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})

FREQ  <- "supplemental_data/phenotypes/pos1_2023_sample_frequencies.csv.gz"
TRAIT <- "supplemental_data/phenotypes/pos1_2023_association_traits.csv.gz"
OUT   <- "plots/diagnostics"
CUTOFF <- 3L        # reported for all three; 3 is the headline
NBOOT  <- 500L
SEED   <- 100L

stopifnot(file.exists(FREQ), file.exists(TRAIT))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

f  <- fread(cmd = paste("gzcat", shQuote(FREQ)))
tr <- fread(cmd = paste("gzcat", shQuote(TRAIT)))

## The mapped set: strains with a pos-1 VST value. This is also the answer to
## the pool-size question -- the draft's "224" is not this number.
mapped <- tr[!is.na(`vst_ctrl_pos-1_T2`), strain]
cat("== the mapped set ==\n")
cat(sprintf("  trait file rows                : %d\n", nrow(tr)))
cat(sprintf("  with a pos-1 VST measurement   : %d\n", length(mapped)))
cat(sprintf("  without                        : %d\n\n", nrow(tr) - length(mapped)))
stopifnot(length(mapped) == 231L)

## ---------------------------------------------------------------------------
## the shared baseline, asserted rather than assumed
## ---------------------------------------------------------------------------
p <- f[rnai == "pos-1" & depth_cutoff == CUTOFF & strain %in% mapped]
shared <- p[, .(n_distinct = uniqueN(round(ctrl_frq, 12))), by = strain]
cat("== control baseline ==\n")
cat(sprintf("  strains whose ctrl_frq is identical across all four pos-1 replicates: %d of %d\n\n",
            sum(shared$n_distinct == 1L), nrow(shared)))
stopifnot(all(shared$n_distinct == 1L))

## ---------------------------------------------------------------------------
## the abundance confound, quantified
## ---------------------------------------------------------------------------
ab <- p[, .(abs_delta = mean(abs(delta_ctrl), na.rm = TRUE),
            ctrl      = mean(ctrl_frq,        na.rm = TRUE)), by = strain]
rho_ab <- cor(ab$abs_delta, ab$ctrl, method = "spearman")
cat("== abundance confound on the raw delta scale ==\n")
cat(sprintf("  Spearman |delta_ctrl| vs control frequency: %.3f\n", rho_ab))
ab[, q := cut(ctrl, quantile(ctrl, 0:4 / 4), include.lowest = TRUE,
              labels = paste0("Q", 1:4))]
print(ab[, .(n = .N, median_ctrl = median(ctrl),
             median_abs_delta = median(abs_delta)), by = q][order(q)])
cat("\n")

## ---------------------------------------------------------------------------
## repeatability, three scales
## ---------------------------------------------------------------------------
icc <- function(d) {
  m  <- lmer(y ~ 1 + (1 | strain), data = d, REML = TRUE)
  vc <- as.data.table(VarCorr(m))
  vg <- vc[grp == "strain",   vcov]
  ve <- vc[grp == "Residual", vcov]
  list(R = vg / (vg + ve), vg = vg, ve = ve, n = uniqueN(d$strain))
}

scales <- function(dc) {
  q <- f[rnai == "pos-1" & depth_cutoff == dc & strain %in% mapped]

  raw  <- q[, .(y = sum(delta_ctrl, na.rm = TRUE)), by = .(strain, replicate)]
  lfc  <- q[is.finite(log2fc_ctrl), .(y = log2fc_ctrl[1]), by = .(strain, replicate)]
  lfc  <- lfc[strain %in% lfc[, .N, by = strain][N == 4L, strain]]
  rank <- copy(raw)[, y := frank(y) / .N, by = replicate]

  rbindlist(lapply(list(delta_ctrl = raw, log2fc_ctrl = lfc, rank_within_rep = rank),
    function(d) { r <- icc(d)
      data.table(depth_cutoff = dc, n_strain = r$n, n_rep = 4L,
                 var_strain = r$vg, var_resid = r$ve, R = r$R) }),
    idcol = "scale")
}

res <- rbindlist(lapply(c(3L, 5L, 10L), scales))

## bootstrap over strains, on the scale actually being quoted
d <- f[rnai == "pos-1" & depth_cutoff == CUTOFF & strain %in% mapped &
       is.finite(log2fc_ctrl), .(y = log2fc_ctrl[1]), by = .(strain, replicate)]
d <- d[strain %in% d[, .N, by = strain][N == 4L, strain]]
set.seed(SEED)
bs <- vapply(seq_len(NBOOT), function(i) {
  ss <- sample(unique(d$strain), replace = TRUE)
  dd <- rbindlist(lapply(seq_along(ss), function(j) {
    z <- d[strain == ss[j]]; z$strain <- paste0("s", j); z }))
  m <- try(icc(dd)$R, silent = TRUE)
  if (inherits(m, "try-error")) NA_real_ else m
}, numeric(1))
ci <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)

cat("== repeatability by scale ==\n")
print(as.data.frame(res), row.names = FALSE, digits = 4)
cat(sprintf(
  "\n  QUOTE THIS: R = %.2f on log2 fold-change, %d strains x 4 replicates,\n",
  res[scale == "log2fc_ctrl" & depth_cutoff == CUTOFF, R],
  res[scale == "log2fc_ctrl" & depth_cutoff == CUTOFF, n_strain]))
cat(sprintf("              95%% bootstrap interval [%.2f, %.2f] over %d replicates.\n",
            ci[1], ci[2], sum(!is.na(bs))))
cat("              Repeatability, an upper bound on H2 -- not heritability.\n\n")
cat(sprintf("  The raw-delta %.2f is inflated: |delta| tracks abundance at rho %.2f.\n",
            res[scale == "delta_ctrl" & depth_cutoff == CUTOFF, R], rho_ab))

res[, `:=`(boot_lo = NA_real_, boot_hi = NA_real_)]
res[scale == "log2fc_ctrl" & depth_cutoff == CUTOFF,
    `:=`(boot_lo = ci[1], boot_hi = ci[2])]
res[, rho_abs_delta_vs_abundance := rho_ab]

f_out <- file.path(OUT, "TABLE_pos1_repeatability.tsv")
fwrite(res, f_out, sep = "\t")
cat("\nwrote ", f_out, "\n", sep = "")
