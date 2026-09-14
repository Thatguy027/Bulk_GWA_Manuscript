## Baugh L1-starvation traits for association mapping -------------------------
##
## Writes one row per strain with four traits, two per measurement platform:
##
##   slope_baugh   rate of change in pool frequency across days of L1
##   slope_nnls    starvation, from published targeted MIP-seq (baugh) and from
##                 NNLS deconvolution of pooled whole-genome sequence (nnls)
##   PC1_baugh     first principal component of the samples x strains frequency
##   PC1_nnls      matrix, per platform, as a per-strain loading
##
## SLOPES are platform_slopes() from Figure1_common.R -- the same function
## behind Figure 1A -- averaged over the five replicate arms. Baseline samples
## and day 17 are excluded and deltas are taken against day 1, as in the
## published analysis.
##
## PC1 restores the comparison the legacy exploratory script drew and the
## current Figure 1 dropped (scripts/legacy/Figure1_legacy_exploratory.R:146-171):
## spread the long frequencies to a samples x strains matrix, prcomp(center =
## TRUE), take rotation[, 1]. Because the columns are strains, the PC1 loading
## is a per-strain value. N2 is dropped and rows carrying any NA are dropped
## before spreading, which is what restricts the set to the 98 strains that have
## both measurements. All 23 samples are used, baselines included, as in the
## legacy block. The legacy script also had a second, near-identical block using
## log2 frequency for the NNLS side only; the raw-frequency version is used here
## so the two platforms are treated the same way.
##
## SIGN. A principal component's sign is arbitrary, so each PC1 is oriented to
## correlate positively with its own platform's slope. Without this the two PC1
## columns could point in opposite directions and the sign of a mapped effect
## would not be comparable between them. The orientation applied is reported.
##
## TWO REFERENCES. The deconvolution reference is not restricted to the strains
## the MIP panel measured: it carries three extra strains (CX11262, ECA348,
## NIC260). Both versions are produced so a scan can be run either way.
##
##   --source=deposited  the shipped cache. 98 strains, NO PB306, because the
##                       deposited input dropped it (see below). Kept as the
##                       record of what Figure 1 was built on.
##   --source=dep103     the DEPOSITED matrix plus PB306 grafted back, i.e. the
##                       full 103-strain pool. 99 traits.
##   --source=pool100    the same matrix cut to the 100 strains MIP measures.
##                       99 traits.
##
## dep103 and pool100 are the controlled comparison: one matrix, one marker set,
## PB306 in both, differing only in whether the three unmeasured strains are in
## the reference. Both are built by scripts/make_baugh_deposited_fits.R, which
## also verifies that the deposited matrix reproduces the shipped cache exactly
## (max absolute difference 0).
##
## TERMINOLOGY, because two different files get called "the reference".
##
##   the CACHE   supplemental_data/deconvolution/baugh_nnls_with_mipseq.RData,
##               68 KB, in this repository. Already-computed NNLS frequencies
##               joined to MIP. This is what --source=deposited reads and what
##               Figure 1 draws.
##   the INPUT   2024bootstrapINPUT.Rdata, the genotype matrix and counts the
##               cache was computed FROM. Dryad only (Figure1_common.R:60) --
##               not in this repository and not anywhere on this machine.
##
## The genotype matrix inside that input has never been inspected here. Its
## strain list is known only indirectly, from baugh_strain_order.txt (shipped
## precisely because the bootstrap array's dimension is unnamed) and from the
## cache's own strain set. Both say 102 strains.
##
## WHERE PB306 WENT. All seven 2022 input variants in the source project
## (/Users/Stefan/UCLA/Projects/bulkGWAS/baugh_wgs/cluster_data/) carry 103
## strains and every one of them includes PB306. The deposited 2024 input
## carries 102, and the difference is exactly PB306 and nothing else. Why it
## was dropped cannot be determined from here; the file that would say is the
## one that is missing.
##
## That matters because PB306 IS in the pool -- the MIP panel measures it. So
## the deposited reference disagrees with the pool definition in BOTH
## directions: it is missing PB306, a real pool strain it therefore cannot
## estimate at all, and it carries CX11262, ECA348 and NIC260, which the pool
## definition does not contain.
##
## CONSEQUENCE FOR THESE TWO FILES. The restricted frequencies come from
## scripts/DIAG_baugh_refit_restricted_reference.R, which starts from the 2022
## matrix and so HAS PB306, staged into the repository here so this script
## needs no external file. The restricted output therefore gains PB306 because
## of which source matrix it was built from, NOT because of the restriction.
## The two files are not a controlled A/B on the reference alone; that would
## need the deposited input back, or one source matrix cut both ways. The MIP
## columns are unaffected by any of this and should agree exactly between the
## two; the script checks that.
##
## Usage:  Rscript scripts/make_baugh_association_traits.R [--source=deposited|restricted]
## Writes: supplemental_data/phenotypes/baugh_association_traits.csv
##     or  supplemental_data/phenotypes/baugh_association_traits_restricted.csv
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

source("scripts/Figure1_common.R")

args <- commandArgs(trailingOnly = TRUE)
src <- sub("^--source=", "", grep("^--source=", args, value = TRUE))
if (!length(src)) src <- "deposited"
stopifnot(src %in% c("deposited", "dep103", "pool100"))

STAGED <- c(dep103  = file.path(BAUGH, "baugh_nnls_dep103_with_mipseq.tsv.gz"),
            pool100 = file.path(BAUGH, "baugh_nnls_pool100_with_mipseq.tsv.gz"))
OUTFILE <- switch(src,
  deposited = "supplemental_data/phenotypes/baugh_association_traits.csv",
  dep103    = "supplemental_data/phenotypes/baugh_association_traits_dep103.csv",
  pool100   = "supplemental_data/phenotypes/baugh_association_traits_pool100.csv")

message("source: ", src)
freq <- if (src == "deposited") baugh_frequencies() else {
  stopifnot(file.exists(STAGED[[src]]))
  readr::read_tsv(STAGED[[src]], show_col_types = FALSE)
}

## --- slopes, averaged over replicate arms ----------------------------------
slopes <- platform_slopes(freq) %>%
  filter(strain != "N2") %>%
  group_by(strain) %>%
  summarise(slope_nnls  = mean(wgs_slope, na.rm = TRUE),
            slope_baugh = mean(mip_slope, na.rm = TRUE),
            n_arms      = sum(is.finite(wgs_slope)),
            .groups = "drop") %>%
  mutate(across(c(slope_nnls, slope_baugh), ~ ifelse(is.nan(.x), NA_real_, .x)))

## --- PC1 per platform ------------------------------------------------------
pc1 <- function(value_col) {
  m <- freq %>%
    filter(strain != "N2") %>%
    select(sample, strain, frq, published_frq) %>%
    na.omit() %>%
    select(sample, strain, all_of(value_col)) %>%
    pivot_wider(names_from = strain, values_from = all_of(value_col)) %>%
    column_to_rownames("sample") %>%
    as.matrix()
  stopifnot(!anyNA(m))
  message(sprintf("  %-14s matrix %d samples x %d strains",
                  value_col, nrow(m), ncol(m)))
  r <- prcomp(m, center = TRUE)
  ve <- (r$sdev^2 / sum(r$sdev^2))[1]
  message(sprintf("  %-14s PC1 explains %.1f%% of variance", value_col, 100 * ve))
  tibble(strain = rownames(r$rotation), pc1 = r$rotation[, 1])
}

nnls_pc <- pc1("frq")          %>% rename(PC1_nnls  = pc1)
mip_pc  <- pc1("published_frq") %>% rename(PC1_baugh = pc1)

out <- slopes %>%
  inner_join(mip_pc,  by = "strain") %>%
  inner_join(nnls_pc, by = "strain")

## --- orient each PC1 against its own slope ---------------------------------
for (p in list(c("PC1_baugh", "slope_baugh"), c("PC1_nnls", "slope_nnls"))) {
  rho <- cor(out[[p[1]]], out[[p[2]]], method = "spearman", use = "complete.obs")
  if (is.finite(rho) && rho < 0) {
    out[[p[1]]] <- -out[[p[1]]]
    message(sprintf("  flipped %s (was rho = %+.3f against %s)", p[1], rho, p[2]))
  } else {
    message(sprintf("  kept    %s (rho = %+.3f against %s)", p[1], rho, p[2]))
  }
}

out <- out %>%
  select(strain, slope_baugh, PC1_baugh, slope_nnls, PC1_nnls) %>%
  arrange(strain)

stopifnot(!anyDuplicated(out$strain), !("N2" %in% out$strain))

dir.create(dirname(OUTFILE), recursive = TRUE, showWarnings = FALSE)
write_csv(out, OUTFILE)

message("\ncross-platform agreement (Spearman):")
message(sprintf("  slope_baugh vs slope_nnls  rho = %+.3f",
                cor(out$slope_baugh, out$slope_nnls, method = "spearman")))
message(sprintf("  PC1_baugh   vs PC1_nnls    rho = %+.3f",
                cor(out$PC1_baugh, out$PC1_nnls, method = "spearman")))
## if both exist, report how far the two references move the traits
other <- if (src == "pool100")
  "supplemental_data/phenotypes/baugh_association_traits_dep103.csv" else
  "supplemental_data/phenotypes/baugh_association_traits_pool100.csv"
if (file.exists(other)) {
  o <- readr::read_csv(other, show_col_types = FALSE)
  j <- inner_join(out, o, by = "strain", suffix = c("_this", "_other"))
  message(sprintf("\nagainst %s (%d shared strains):", basename(other), nrow(j)))
  for (tr in c("slope_baugh", "PC1_baugh", "slope_nnls", "PC1_nnls"))
    message(sprintf("  %-12s Spearman %+.4f", tr,
      cor(j[[paste0(tr, "_this")]], j[[paste0(tr, "_other")]], method = "spearman")))
  message("  (the two baugh columns are the same measurement and should be ~1)")
}

message(sprintf("\nwrote %s  (%d strains, %d traits)",
                OUTFILE, nrow(out), ncol(out) - 1))
