## eigen_independent_tests.R -------------------------------------------------
## Effective number of independent tests for each mapping panel, from the
## eigenvalues of the marker correlation matrix.
##
##   Rscript scripts/eigen_independent_tests.R
##
## Why
## ---
## A Bonferroni threshold over every tested marker assumes the markers are
## independent. They are not: the C. elegans wild-isotype panel carries long
## haplotype blocks, so a strict Bonferroni over ~10^5-10^6 markers is far too
## conservative. The standard correction is to replace the marker count with an
## effective number of independent tests estimated from the LD structure.
##
## Method
## ------
## Per chromosome, take the panel's genotypes at the markers actually tested,
## standardise each marker, and obtain the eigenvalues of the marker-by-marker
## correlation matrix R. R is M x M with M up to ~10^5, which is not formable
## directly, but its non-zero eigenvalues are exactly those of the n x n matrix
##
##     A = Zc Zc' / (n - 1)
##
## where Zc is the n x M column-standardised genotype matrix. The remaining
## M - n eigenvalues are zero. trace(A) must equal M, which the script asserts.
##
## Two summaries are reported from those eigenvalues:
##
##   Li & Ji (2005)   M_eff = sum over all M eigenvalues of
##                    [ 1(lambda >= 1) + (lambda - floor(lambda)) ]
##                    Heredity 95:221-227. This is the citable default.
##
##   99.5% variance   the number of leading eigenvalues whose cumulative sum
##                    reaches 99.5% of M. Reported alongside because it is the
##                    convention in some C. elegans GWA pipelines.
##
## Values are summed across chromosomes, since LD does not span chromosomes.
##
## IMPORTANT CAVEAT
## ----------------
## With far more markers than strains, R is rank-deficient: it has at most n
## non-zero eigenvalues, so any eigenvalue-based M_eff is bounded above by
## roughly n per chromosome, i.e. by the sample size rather than by the LD
## structure alone. The estimate is therefore a ceiling on resolvable
## dimensions, not a marker count, and it should be read as "this panel can
## support about this many independent tests", which is the quantity the
## threshold needs.
##
## Output
##   data/eigen_independent_tests.tsv
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

PLINK  <- Sys.which("plink"); if (!nzchar(PLINK)) PLINK <- "~/bin/plink"
GENO   <- "data/genotypes/CeNDR20210121_Plink"
OUTTSV <- "data/eigen_independent_tests.tsv"
TMP    <- file.path(tempdir(), "eigen")
dir.create(TMP, showWarnings = FALSE, recursive = TRUE)

CHROMS  <- c("I", "II", "III", "IV", "V", "X")
VAR_FRAC <- 0.995

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## ---------------------------------------------------------------------------
## the panels to profile: which strains were mapped, and at which markers
## ---------------------------------------------------------------------------
pooled_tr <- read_csv("data/pooled_RNAi_expt/reanalysis/vst_association_traits.csv",
                      show_col_types = FALSE)
pos1_tr   <- read_csv("data/pos1_original/updated_analysis/association_traits.csv",
                      show_col_types = FALSE)

PANELS <- list(
  pooled_RNAi_expt = list(
    strains = pooled_tr$strain[!is.na(pooled_tr$`vst_ctrl_mig-6_T2`)],
    gwas    = "data/pooled_RNAi_expt/reanalysis/mapping/vst_ctrl_mig-6_T2_loco_results.csv.gz"),
  pos1_2023 = list(
    strains = pos1_tr$strain[!is.na(pos1_tr$`vst_ctrl_pos-1_T2`)],
    gwas    = "data/pos1_original/updated_analysis/vst_ctrl_pos-1_T2_loco_results.csv.gz")
)

## ---------------------------------------------------------------------------
## eigenvalues of the marker correlation matrix for one chromosome
## ---------------------------------------------------------------------------
chrom_eigen <- function(chrom, strains, markers, tag) {
  n_tested <- length(markers)
  keep <- file.path(TMP, "keep.txt"); ext <- file.path(TMP, "ext.txt")
  writeLines(paste(strains, strains), keep)       # .fam has FID == IID == strain
  writeLines(markers, ext)
  out <- file.path(TMP, paste0(tag, "_", chrom))

  cmd <- sprintf(paste("%s --bfile %s --keep %s --extract %s --recode A",
                       "--allow-extra-chr --out %s --silent"),
                 PLINK, file.path(GENO, chrom), keep, ext, out)
  st <- system(cmd)
  raw <- paste0(out, ".raw")
  if (st != 0 || !file.exists(raw)) {
    warning("plink failed for ", chrom, " (", tag, ")"); return(NULL)
  }

  d <- fread(raw)
  ## .raw: FID IID PAT MAT SEX PHENOTYPE then one column per marker
  Z <- as.matrix(d[, -(1:6)])
  n <- nrow(Z)

  ## Missing calls are mean-imputed rather than dropped. Dropping any marker
  ## with a single missing call discarded ~40-50% of the tested markers at this
  ## genotyping rate, which both understates the LD structure and silently
  ## changed the marker count the Bonferroni was computed from.
  cm <- colMeans(Z, na.rm = TRUE)
  allna <- is.nan(cm)
  if (any(allna)) { Z <- Z[, !allna, drop = FALSE]; cm <- cm[!allna] }
  idx <- which(is.na(Z), arr.ind = TRUE)
  n_imputed <- nrow(idx)
  if (n_imputed) Z[idx] <- cm[idx[, 2]]

  ## a marker with no variance in this panel carries no information
  sdv <- apply(Z, 2, sd)
  Z <- Z[, sdv > 0, drop = FALSE]
  M <- ncol(Z)
  if (M < 2 || n < 3) return(NULL)

  Zc <- scale(Z)                                  # centre and scale each marker
  A  <- tcrossprod(Zc) / (n - 1)                  # n x n, same non-zero eigenvalues as R
  lam <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  lam[lam < 0] <- 0                               # tiny negatives from rounding

  ## trace(A) must equal M; a large mismatch means the algebra or input is wrong
  tr_err <- abs(sum(lam) - M) / M
  stopifnot(tr_err < 1e-6)

  ## Li & Ji (2005) over all M eigenvalues; the M - n zeros contribute nothing
  lj <- sum(as.numeric(lam >= 1) + (lam - floor(lam)))

  ## leading eigenvalues reaching VAR_FRAC of the total
  cs <- cumsum(sort(lam, decreasing = TRUE)) / M
  n_var <- which(cs >= VAR_FRAC)[1]

  file.remove(list.files(TMP, pattern = paste0(tag, "_", chrom), full.names = TRUE))
  tibble(chrom = chrom, n_strain = n, n_tested = n_tested, n_marker = M,
         frac_calls_imputed = n_imputed / (n * max(M, 1)),
         rank_nonzero = sum(lam > 1e-8),
         M_eff_liji = lj, M_eff_var995 = n_var, trace_err = tr_err)
}

## ---------------------------------------------------------------------------
res <- imap_dfr(PANELS, function(p, nm) {
  msg(nm, ": ", length(p$strains), " strains")
  gw <- fread(p$gwas, select = c("chr", "rs"))
  gw <- gw[chr %in% CHROMS]
  msg("  ", nrow(gw), " tested markers")

  per <- map_dfr(CHROMS, function(ch) {
    mk <- gw[chr == ch, rs]
    if (!length(mk)) return(NULL)
    r <- chrom_eigen(ch, p$strains, mk, tag = nm)
    if (!is.null(r))
      msg("  ", ch, ": ", r$n_marker, " markers, rank ", r$rank_nonzero,
          ", Li&Ji ", round(r$M_eff_liji, 1), ", 99.5% ", r$M_eff_var995)
    r
  })
  per %>% mutate(panel = nm, .before = 1)
})

tot <- res %>% group_by(panel) %>%
  summarise(n_strain = max(n_strain), n_tested = sum(n_tested),
            n_marker = sum(n_marker),
            M_eff_liji = sum(M_eff_liji), M_eff_var995 = sum(M_eff_var995),
            .groups = "drop") %>%
  ## Bonferroni is quoted over every TESTED marker, matching the figures
  mutate(thr_bonferroni = -log10(0.05 / n_tested),
         thr_eigen_liji = -log10(0.05 / M_eff_liji),
         thr_eigen_var995 = -log10(0.05 / M_eff_var995))

write_tsv(bind_rows(res %>% mutate(scope = "chromosome"),
                    tot %>% mutate(scope = "total", chrom = "all")),
          OUTTSV)

cat("\n== per chromosome ==\n")
print(as.data.frame(res %>% transmute(panel, chrom, n_tested, n_used = n_marker,
  pct_imputed = round(100 * frac_calls_imputed, 2),
  rank = rank_nonzero, M_eff_liji = round(M_eff_liji, 1), M_eff_var995)),
  row.names = FALSE)

cat("\n== totals and thresholds ==\n")
print(as.data.frame(tot %>% transmute(panel, n_strain, n_tested, n_marker,
  M_eff_liji = round(M_eff_liji), M_eff_var995,
  `Bonferroni -log10p` = round(thr_bonferroni, 2),
  `eigen (Li&Ji)` = round(thr_eigen_liji, 2),
  `eigen (99.5%)` = round(thr_eigen_var995, 2))), row.names = FALSE)

msg("wrote ", OUTTSV)
