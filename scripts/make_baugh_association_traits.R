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
##   --source=deposited   the shipped cache, full reference, 98 strains
##   --source=restricted  refit with the reference cut to the MIP pool,
##                        99 strains (it gains PB306, which the deposited
##                        reference lacks)
##
## The restricted frequencies come from scripts/DIAG_baugh_refit_restricted_
## reference.R, staged into the repository so this script needs no external
## file. That refit starts from a near relative of the deposited input -- 103
## strains including PB306, where the deposit has 102 without it -- so the two
## outputs are not a controlled A/B on the reference alone. The MIP columns are
## unaffected and should agree between them; the script checks that.
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
stopifnot(src %in% c("deposited", "restricted"))

RESTRICTED_FREQ <- file.path(BAUGH, "baugh_nnls_restricted_with_mipseq.tsv.gz")
OUTFILE <- if (src == "restricted")
  "supplemental_data/phenotypes/baugh_association_traits_restricted.csv" else
  "supplemental_data/phenotypes/baugh_association_traits.csv"

message("source: ", src)
freq <- if (src == "restricted") {
  stopifnot(file.exists(RESTRICTED_FREQ))
  readr::read_tsv(RESTRICTED_FREQ, show_col_types = FALSE)
} else baugh_frequencies()

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
other <- if (src == "restricted")
  "supplemental_data/phenotypes/baugh_association_traits.csv" else
  "supplemental_data/phenotypes/baugh_association_traits_restricted.csv"
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
