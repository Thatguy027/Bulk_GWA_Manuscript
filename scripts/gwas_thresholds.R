## gwas_thresholds.R ---------------------------------------------------------
## Significance thresholds for the pooled association scans, sourced by the
## figure scripts so the numbers live in one place.
##
##   source("scripts/gwas_thresholds.R")
##   th <- gwas_thresholds("pooled_RNAi_expt")   # or "pos1_2023"
##   th$bonferroni ; th$eigen ; th$m_eff ; th$n_marker
##
## Two thresholds are carried deliberately:
##
##   bonferroni   0.05 / (every tested marker). Assumes marker independence,
##                which is false in this panel, so it is conservative. Kept
##                because it is the number readers expect to see.
##   eigen        0.05 / (effective number of independent tests), where the
##                effective count comes from the eigenvalues of the marker
##                correlation matrix by Li & Ji (2005). This is the threshold
##                appropriate to the LD in the panel.
##
## Both are computed by scripts/eigen_independent_tests.R, which writes
## data/eigen_independent_tests.tsv. Run that first if the file is missing or
## the panels change.
## ---------------------------------------------------------------------------

gwas_thresholds <- function(panel,
                            file = "data/eigen_independent_tests.tsv",
                            method = c("liji", "var995")) {
  method <- match.arg(method)
  if (!file.exists(file))
    stop("missing ", file, " -- run scripts/eigen_independent_tests.R first",
         call. = FALSE)
  d <- utils::read.delim(file, check.names = FALSE)
  d <- d[d$scope == "total" & d$panel == panel, , drop = FALSE]
  if (!nrow(d))
    stop("panel '", panel, "' not in ", file,
         " (have: ", paste(unique(d$panel), collapse = ", "), ")", call. = FALSE)
  m_eff <- if (method == "liji") d$M_eff_liji else d$M_eff_var995
  list(panel      = panel,
       n_strain   = d$n_strain,
       n_marker   = d$n_tested,
       m_eff      = m_eff,
       method     = method,
       bonferroni = -log10(0.05 / d$n_tested),
       eigen      = -log10(0.05 / m_eff))
}
