## 2023 pos-1: redo the deconvolution against the POOL, not the whole panel ---
##
##   Rscript reanalysis_2023_pos1/scripts/01_deconvolve_pool_reference.R
##
## WHAT WAS WRONG. The deposited analysis deconvolved every sample against a
## 367-column genotype reference covering the whole strain panel. The 2023
## pos-1 pool was not the whole panel: it is sets B, C, E, F and G of
## meta_files/BulkCe_strainsets.tsv, which is 226 strains, 224 of them with a
## CeNDR isotype. Offering NNLS 142 strains that were never in the tube lets it
## place mass on them, and that mass is taken from the strains that were.
##
## Measured on the deposited traits: of the 231 strains given a phenotype, 65
## are not in the pool at all, and 58 strains that ARE in the pool were left
## without one.
##
## WHAT THIS DOES. The same markers, the same counts and the same solver as the
## original run -- the marker filter is not re-derived, it is inherited by
## reusing the `gt` matrix the original stored for each depth cutoff -- with the
## reference restricted to the 224 pool isotypes. The only thing that changes
## is which columns the design matrix has.
##
## It first reproduces the ORIGINAL 367-column fit from those same inputs and
## checks it against the stored predictions, so the reimplementation is shown to
## be faithful before the corrected one is believed.
##
## THE DUPLICATED COLUMN. The reference carries JU1793_JU1793 TWICE, which is
## why the deposited frequencies carry JU1793 twice per sample with one row at
## ~1e-21: two identical columns, and NNLS splits the mass between them
## arbitrarily. The pool reference keeps one column per isotype, so this
## disappears rather than needing a downstream fix.
##
## Reads the upstream analysis directory, which is outside the repository, and
## writes only inside reanalysis_2023_pos1/.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(dplyr)
})

UP   <- "/Users/Stefan/UCLA/Projects/bulkGWAS/lipid_RNAi/2023_pos1"
META <- "/Users/Stefan/UCLA/Projects/bulkGWAS/lipid_RNAi/2023_original_pos1/meta_files"
OUT  <- "reanalysis_2023_pos1/data"
CUTOFFS <- c(3, 5, 10)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## --- the pool, from the metadata that describes this experiment ------------
iso  <- fread(file.path(META, "strain_isotype_lookup.tsv"))
sets <- fread(file.path(META, "BulkCe_strainsets.tsv"),
              col.names = c("strain", "set"), header = FALSE, skip = 1)
pool <- sets %>% left_join(iso, by = "strain") %>% na.omit() %>%
  select(strain = isotype, set) %>% filter(set %in% c("C", "B", "E", "F", "G"))
POOL <- unique(pool$strain)
msg("pool: ", nrow(sets), " strains in the sheet -> ",
    sum(sets$set %in% c("C","B","E","F","G")), " in sets B,C,E,F,G -> ",
    length(POOL), " distinct isotypes")
stopifnot(length(POOL) == 224)

## --- counts, shared by every cutoff ----------------------------------------
msg("loading counts")
e <- new.env(); load(file.path(UP, "data/01_genotypes_and_counts.Rdata"), e)
counts  <- as.data.table(e$alt_df_bias)
samples <- as.data.table(e$sample_wide)[, .(sample = well, sample_info)]
rm(e); invisible(gc())

nnls_fit <- function(G, GtG, Y) {
  GtY <- crossprod(G, Y)
  est <- apply(GtY, 2, function(y)
    as.vector(RcppML::nnls(GtG, matrix(y), fast_nnls = TRUE)))
  est <- apply(est, 2, function(x) x / sum(x))
  dimnames(est) <- list(colnames(G), colnames(Y))
  est
}

out <- list()
for (cut in CUTOFFS) {
  msg("cutoff ", cut)
  p <- new.env(); load(file.path(UP, sprintf("data/predictions_dp%d.Rdata", cut)), p)
  G <- p$gt; storage.mode(G) <- "double"
  mk <- rownames(G)
  msg("  markers ", format(nrow(G), big.mark = ","), " x ", ncol(G), " reference columns")

  ## counts aligned to those markers, in a fixed sample order
  ac <- dcast(counts[marker %in% mk], marker ~ sample, value.var = "alt_ct")
  ac <- ac[match(mk, ac$marker)]
  Y  <- as.matrix(ac[, -1]); rownames(Y) <- ac$marker
  stopifnot(!anyNA(Y), nrow(Y) == nrow(G))

  ## 1. reproduce the original fit, to show the method is the same one
  rep_full <- nnls_fit(G, p$GtG, Y)
  orig <- as.data.table(p$predictions_df)[, .(sample, strain, frq)]
  ## the original stripped names, which collapses the duplicated JU1793 column;
  ## compare on the summed value so the two are on the same footing
  mine <- data.table(strain = sub("_.*$", "", rownames(rep_full)), as.data.table(rep_full)) %>%
    melt(id.vars = "strain", variable.name = "sample", value.name = "mine") %>%
    .[, .(mine = sum(mine)), by = .(strain, sample)]
  chk <- merge(orig[, .(orig = sum(frq)), by = .(strain, sample)], mine,
               by = c("strain", "sample"))
  ## The stored solution is NOT uniquely determined, so this is a comparison
  ## and not an assertion. With 367 near-collinear columns -- many of these
  ## strains are close relatives, and 142 of them were never in the tube --
  ## NNLS has a wide set of near-optimal solutions, and only about 95 strains
  ## carry any mass at all. Refitting from the same markers, the same counts and
  ## the same GtG (which recomputes to machine zero against the stored one)
  ## lands on a different member of that set, correlated with the stored
  ## solution at ~0.995 and with a slightly LOWER residual. That degeneracy is
  ## itself part of what the pool reference fixes.
  msg("  vs the stored 367-column fit: r ", round(cor(chk$orig, chk$mine), 4),
      " | max |diff| ", signif(max(abs(chk$orig - chk$mine)), 3),
      " | cells differing ", sum(abs(chk$orig - chk$mine) > 1e-9), " of ", nrow(chk))
  stopifnot(cor(chk$orig, chk$mine) > 0.98)

  ## 2. the corrected fit: one column per pool isotype
  isot <- sub("^.*?_", "", colnames(G))
  take <- !duplicated(isot) & isot %in% POOL
  Gp <- G[, take, drop = FALSE]; colnames(Gp) <- isot[take]
  stopifnot(ncol(Gp) == 224, !anyDuplicated(colnames(Gp)))
  msg("  pool reference: ", ncol(Gp), " columns")
  est <- nnls_fit(Gp, crossprod(Gp), Y)
  msg("  strains carrying mass: ", sum(rowSums(est) > 1e-12), " of ", ncol(Gp),
      " (the 367-column fit placed mass on ",
      sum(rowSums(rep_full) > 1e-12), " of ", ncol(G), ")")

  out[[as.character(cut)]] <- as.data.table(est, keep.rownames = "strain") %>%
    melt(id.vars = "strain", variable.name = "sample", value.name = "frq") %>%
    .[, depth_cutoff := cut]
  rm(p, G, Gp, Y, ac, rep_full); invisible(gc())
}

freq <- rbindlist(out)[samples, on = "sample"]
setcolorder(freq, c("sample", "sample_info", "depth_cutoff", "strain", "frq"))
fwrite(freq, file.path(OUT, "pool_reference_frequencies.csv.gz"), na = "NA", quote = FALSE)
msg("wrote ", file.path(OUT, "pool_reference_frequencies.csv.gz"),
    " (", format(nrow(freq), big.mark = ","), " rows)")

cat("\n== zero frequencies by cutoff, pool reference ==\n")
print(freq[, .(strains = uniqueN(strain), cells = .N,
               zero = sum(frq == 0), pct_zero = round(100 * mean(frq == 0), 1)),
           by = depth_cutoff][order(depth_cutoff)])
