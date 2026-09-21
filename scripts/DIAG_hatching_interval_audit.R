## Recompute every hatching percentage and interval the draft quotes -----------
##
##   Rscript scripts/DIAG_hatching_interval_audit.R
##
## The draft quotes binomial intervals on the hatching assays, and METHODS.txt
## declares Wilson score intervals. This prints Wilson AND Clopper-Pearson side
## by side for every strain the draft mentions, so a quoted interval can be
## traced to a convention rather than argued about. It is what established that
## the draft's allele-swap intervals are Clopper-Pearson while its NIL intervals
## are Wilson, and it is the check to rerun before accepting any edit to those
## numbers -- claude_science/percentage_audit/ carries the resulting edit list
## and its applied/outstanding status.
##
## Reads only supplemental_data, so it runs from a clone. Nothing in the
## manuscript reads its output; it prints.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
wilson <- function(x, n, z = qnorm(0.975)) {
  p <- x / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  sprintf("%.1f%% (%.1f-%.1f)", 100 * p, 100 * (c - h), 100 * (c + h))
}
cp <- function(x, n) sprintf("%.1f-%.1f",
  100 * qbeta(0.025, x, n - x + 1), 100 * qbeta(0.975, x + 1, n - x))

ju <- fread("supplemental_data/hatching_assays/ju_allele_swaps_hatching.csv")
setnames(ju, make.names(names(ju)))
ju[, hatched := n_plated - n_unhatched]
cat("=== JU allele swaps, pos-1 food ===\n")
for (s in c("JU1793","wSZ200","JU2466_A","JU2466_B","wSZ206","wSZ209","wSZ208")) {
  r <- ju[strain == s & condition == "pos"]
  if (!nrow(r)) next
  cat(sprintf("  %-9s %-16s Wilson %-22s CP %s\n", s, r$genotype[1],
              wilson(r$hatched, r$n_plated), cp(r$hatched, r$n_plated)))
}
p <- ju[strain %in% c("JU2466_A","JU2466_B") & condition == "pos",
        .(h = sum(hatched), n = sum(n_plated))]
cat(sprintf("  %-9s %-16s Wilson %-22s CP %s\n", "A+B", "pooled",
            wilson(p$h, p$n), cp(p$h, p$n)))

n2 <- fread("supplemental_data/hatching_assays/n2_allele_swaps_hatching.tsv")
n2[, hatched := n_plated - n_unhatched]
cat("\n=== N2 dose series ===\n")
for (d in c(0, 25, 50, 75, 100)) {
  a <- n2[strain == "N2" & condition == d]
  b <- n2[strain %in% c("wSZ203","wSZ204") & condition == d,
          .(h = sum(hatched), n = sum(n_plated))]
  cat(sprintf("  %3s%%  N2 %-22s | 96K pooled %-22s\n", d,
              wilson(a$hatched, a$n_plated), wilson(b$h, b$n)))
}
cat("\n=== NIL series, pos-1 ===\n")
nl <- fread("supplemental_data/hatching_assays/nil_series_hatching.tsv")
setnames(nl, make.names(names(nl)))
nl[, hatched := plated.embryo - unhatched]
for (s in c("JU2466","JU1793","wSZ153")) {
  r <- nl[Strain == s & condition == "pos1"]
  if (nrow(r)) cat(sprintf("  %-8s Wilson %s\n", s, wilson(r$hatched, r$plated.embryo)))
}
