## Build the thinned pooled/cross bundle shipped as supplemental data --------
##
##   Rscript scripts/make_thinned_bundle.R
##     data/pooled_cross_intersection/bundle.rds        (full, 258.7 MB)
##     -> supplemental_data/mapping/pooled_cross_bundle_thinned.rds  (~9.6 MB)
##
## The full bundle carries every marker: 8,297,244 rows across 18 cross
## contrast scans and 644,020 rows of pooled GWAS. Every figure that reads it
## bins to 5 kb (cross scans) or 10 kb (GWAS) before drawing, so the full
## resolution is never displayed. This keeps one row per bin -- the row with
## the highest LOD, so peaks are preserved exactly -- and cuts the file 27x.
##
## WHAT THIS COSTS, MEASURED
## scripts/compare_full_vs_thinned_bundle.R rebuilds Figure 2 from both and
## reports the difference. As of the 5 kb / 10 kb settings below:
##
##   every peak position and every peak LOD: IDENTICAL
##   the GWAS-peak-to-cross-interval distances: IDENTICAL
##   SUPP_FIG_XX_cross_contrast_panels: pixel-for-pixel identical
##   Figure 2: 95.7% of pixels identical, 83 of 6.3M differing by >32/255
##   QTL support interval WIDTHS: contract by 0-9 kb
##
## The width contraction is the one real difference and it has a mechanism:
## interval bounds are the outermost marker within a 5% LOD drop of the peak,
## and thinning can discard that marker, pulling each edge in by up to one bin.
## It is always a contraction, never an expansion, and it is under two bins.
## Interval widths quoted in the manuscript should come from whichever bundle
## the figures were built from -- they are not interchangeable to the kilobase.
##
## Peak-preserving thinning, not subsampling: taking every nth marker or the
## first per bin would move peak positions and lower peak LODs.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
})

FULL <- "data/pooled_cross_intersection/bundle.rds"
OUT  <- "supplemental_data/mapping/pooled_cross_bundle_thinned.rds"
BIN_SCAN <- 5e3     # cross contrast scans, matching the figures' display bin
BIN_GWAS <- 1e4     # pooled GWAS, matching the Manhattan panels

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

if (!file.exists(FULL)) {
  stop("full bundle not found at ", FULL, ".\n",
       "  It is not in the git repository -- see DATA_AVAILABILITY.md.\n",
       "  Rebuild it with scripts/pooled_cross_intersection_prep.R, or take it\n",
       "  from the Dryad archive.")
}
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

b <- readRDS(FULL)
msg("full bundle: ", round(file.size(FULL) / 1e6, 1), " MB, ",
    length(b$scans), " scans")

## bin, keep the maximum-LOD row per bin. slice_max with n = 1 and
## with_ties = FALSE keeps exactly one row and is stable.
thin_scan <- function(d, bin = BIN_SCAN) {
  d %>% mutate(.b = floor(physical.position / bin)) %>%
    group_by(chrom, .b) %>% slice_max(LOD, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(-.b) %>% arrange(chrom, physical.position)
}
thin_gwas <- function(d, bin = BIN_GWAS) {
  d %>% mutate(.b = floor(ps / bin)) %>%
    group_by(gene, chr, .b) %>% slice_max(neglog10p, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(-.b) %>% arrange(gene, chr, ps)
}

n_scan_before <- sum(vapply(b$scans, nrow, integer(1)))
n_gwas_before <- nrow(b$gwas)

bt <- b
bt$scans <- lapply(b$scans, thin_scan)
bt$gwas  <- thin_gwas(b$gwas)
bt$thinned <- list(bin_scan = BIN_SCAN, bin_gwas = BIN_GWAS,
                   rows_scan_before = n_scan_before,
                   rows_scan_after = sum(vapply(bt$scans, nrow, integer(1))),
                   rows_gwas_before = n_gwas_before,
                   rows_gwas_after = nrow(bt$gwas),
                   source = FULL, built = Sys.time())

## every peak must survive: assert it rather than trust it
for (k in names(b$scans)) {
  stopifnot(isTRUE(all.equal(max(b$scans[[k]]$LOD, na.rm = TRUE),
                             max(bt$scans[[k]]$LOD, na.rm = TRUE))))
}
msg("  peak LOD preserved in all ", length(b$scans), " scans")

saveRDS(bt, OUT, compress = "xz")
msg("cross scans: ", format(n_scan_before, big.mark = ","), " -> ",
    format(bt$thinned$rows_scan_after, big.mark = ","), " rows (",
    sprintf("%.1f%%", 100 * bt$thinned$rows_scan_after / n_scan_before), ")")
msg("pooled GWAS: ", format(n_gwas_before, big.mark = ","), " -> ",
    format(bt$thinned$rows_gwas_after, big.mark = ","), " rows (",
    sprintf("%.1f%%", 100 * bt$thinned$rows_gwas_after / n_gwas_before), ")")
msg("wrote ", OUT, "  ", round(file.size(OUT) / 1e6, 2), " MB (",
    sprintf("%.0fx", file.size(FULL) / file.size(OUT)), " smaller)")
