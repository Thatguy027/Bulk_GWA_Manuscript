#!/usr/bin/env Rscript
## Genome-wide maxima -> permutation threshold ----------------------------
##
## Input is one row per (trait, chromosome, permutation, max -log10 p). The
## genome-wide maximum for a permutation is the largest of its six chromosome
## maxima; the (1 - alpha) quantile of those is the threshold.
##
## Permutation id 0 is the OBSERVED phenotype, not a permutation. It is excluded
## from the null distribution and used for the empirical p-value:
##     p = (1 + #{permuted maxima >= observed}) / (1 + n_perm)
## The +1 in both places is deliberate. Without it an observed maximum larger
## than every permutation reports p = 0, which no finite number of permutations
## can support; the smallest honest value is 1/(1 + n_perm).

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})
args <- commandArgs(TRUE)
get <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) return(d); args[i + 1] }
maxf   <- get("--maxima"); alphas <- as.numeric(strsplit(get("--alpha", "0.05"), ",")[[1]])
n_perm <- as.integer(get("--n_perm", NA))
## The observed genome-wide maximum this pipeline reproduces is a KNOWN number:
## it is the maximum of the scan being thresholded. Asserting it is the only
## check that the model, kinship and phenotype all match that scan -- the marker
## count assertion in PLINK_CONVERT catches the panel and the marker set, and
## passed while the kinship was still wrong (-gk 1 gave 8.5700 against 8.8361).
## Pass --expect_observed_max 0 to disable, e.g. for a trait with no shipped scan.
expect_obs <- as.numeric(get("--expect_observed_max", "0"))
obs_tol    <- as.numeric(get("--observed_tol", "0.01"))

d <- fread(maxf, header = FALSE,
           col.names = c("trait", "chrom", "perm", "max_lp"))
d[, max_lp := suppressWarnings(as.numeric(max_lp))]
cat("rows read:", nrow(d), " | traits:", uniqueN(d$trait),
    " | chromosomes:", uniqueN(d$chrom), "\n")

na_rows <- d[is.na(max_lp)]
if (nrow(na_rows)) {
  cat("WARNING:", nrow(na_rows), "chromosome scans returned no usable p-value.\n")
  cat("  These are dropped. If the count is large the run is not trustworthy:\n")
  print(na_rows[, .N, by = .(trait, chrom)])
  d <- d[!is.na(max_lp)]
}

## a permutation must have a maximum from EVERY chromosome or its genome-wide
## maximum is not genome-wide
nchr <- uniqueN(d$chrom)
cover <- d[, .(n = .N), by = .(trait, perm)]
short <- cover[n < nchr]
if (nrow(short)) {
  cat("WARNING: dropping", nrow(short), "permutations missing a chromosome\n")
  d <- d[!short, on = .(trait, perm)]
}

gw <- d[, .(gw_max = max(max_lp)), by = .(trait, perm)]
fwrite(gw[order(trait, perm)], "permutation_maxima.tsv", sep = "\t")

res <- rbindlist(lapply(split(gw, gw$trait), function(x) {
  obs  <- x[perm == 0L, gw_max]
  null <- x[perm != 0L, gw_max]
  if (!length(obs)) obs <- NA_real_
  q <- quantile(null, 1 - alphas, names = FALSE, type = 7)
  data.table(trait = x$trait[1], n_permutations = length(null),
             alpha = alphas, threshold = round(q, 4),
             observed_max = round(obs, 4),
             empirical_p = if (is.na(obs)) NA_real_ else
               round((1 + sum(null >= obs)) / (1 + length(null)), 5),
             null_median = round(median(null), 4),
             null_max = round(max(null), 4))
}))
fwrite(res, "permutation_thresholds.tsv", sep = "\t")
cat("\n== thresholds ==\n"); print(as.data.frame(res), row.names = FALSE)

if (is.finite(expect_obs) && expect_obs > 0) {
  got <- unique(res$observed_max)
  if (length(got) != 1L || !is.finite(got) || abs(got - expect_obs) > obs_tol) {
    cat("\n")
    stop("observed genome-wide maximum is ", paste(got, collapse = ", "),
         ", expected ", expect_obs, " (tolerance ", obs_tol, ").\n",
         "  The permutation scan is not reproducing the scan it thresholds, so\n",
         "  its threshold does not apply to that scan. Check the GEMMA model\n",
         "  (-gk 1 vs -gk 2, -lmm), the kinship, and the phenotype column\n",
         "  before using any number in this run. --expect_observed_max 0 skips\n",
         "  this check for a trait with no shipped scan to compare against.",
         call. = FALSE)
  }
  cat("\nobserved maximum ", got, " matches the shipped scan (", expect_obs,
      ") within ", obs_tol, "\n", sep = "")
}

## Reference lines: the analytic thresholds this is meant to replace. Values are
## those the manuscript quotes for the 231-strain pos-1 panel; they are drawn
## for orientation only and are not recomputed here.
REF <- data.table(label = c("eigen (Li & Ji)", "Bonferroni"), value = c(4.60, 6.97))

pl <- ggplot(gw[perm != 0L], aes(gw_max)) +
  geom_histogram(bins = 40, fill = "grey78", colour = "grey35", linewidth = 0.25) +
  geom_vline(data = res[alpha == min(alphas)],
             aes(xintercept = threshold), colour = "#0E6B62", linewidth = 0.7) +
  geom_vline(data = REF, aes(xintercept = value, linetype = label),
             colour = "grey30", linewidth = 0.4) +
  geom_vline(data = unique(res[, .(trait, observed_max)]),
             aes(xintercept = observed_max), colour = "#B23A48",
             linewidth = 0.7, linetype = "longdash") +
  facet_wrap(~ trait, scales = "free") +
  scale_linetype_manual(values = c("eigen (Li & Ji)" = "dotted",
                                   "Bonferroni" = "dashed"), name = NULL) +
  labs(x = expression(paste("genome-wide maximum  ", -log[10], " p")),
       y = "permutations",
       title = "Permutation null for the genome-wide maximum",
       subtitle = paste0("teal = permutation threshold at alpha ", min(alphas),
                         "; red dashed = observed scan")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", plot.title.position = "plot")

ggsave("permutation_threshold.pdf", pl, width = 8, height = 5)
ggsave("permutation_threshold.png", pl, width = 8, height = 5, dpi = 300, bg = "white")
cat("\nwrote permutation_maxima.tsv, permutation_thresholds.tsv and the plot\n")
