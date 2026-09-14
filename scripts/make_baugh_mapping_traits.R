## The Baugh L1 mapping table ------------------------------------------------
##
## One file to run association scans from, assembled from the three trait
## builders rather than recomputing anything:
##
##   published_slope_baugh  eLife Slope, verbatim               MIP-seq
##   published_pc1_baugh    eLife PC1, verbatim                 MIP-seq
##   delta_slope_baugh      difference-based slope              MIP-seq
##   delta_pc1_baugh        difference-based PC1                MIP-seq
##   slope_nnls             published recipe, floored           pooled WGS
##   pc1_nnls               published recipe, floored           pooled WGS
##
## TWO PARAMETERISATIONS OF THE SAME MIP DATA. The published pair are log
## ratios against the baseline sample: log2(f_day / f_baseline), with PC1 from
## prcomp(scale, center) over the 20 rep-by-day columns and Slope the
## regression of that ratio on day excluding day 17. The delta pair are the
## difference-based traits this repository already used for Figure 1 --
## frequency minus its day-1 value, regressed on day. They are not
## interchangeable: they agree at Spearman 0.89 on slope but only 0.38 on PC1,
## so a scan on delta_pc1_baugh is not a scan on the published phenotype.
## Both are shipped because the delta parameterisation is the one that survives
## the deconvolution better (see below) while the published one is what the
## eLife mapping used.
##
## THE WGS COLUMNS USE THE PUBLISHED RECIPE, not the delta one, because after
## the low-frequency floor was tuned it agrees with the published traits far
## better: PC1 0.814 and slope 0.884, against 0.376 and 0.886 for the delta
## version. The delta slope is as good; the delta PC1 is not close.
##
## Agreement with the published traits, dep103 reference, n = 99:
##
##                        vs published_slope   vs published_pc1
##   delta_*  (MIP)             0.889                0.331
##   slope/pc1_nnls (WGS)       0.884                0.814
##   ceiling (MIP, 15 cols)     0.983                0.961
##
## The ceiling is MIP scored against itself on the 15 rep-by-day columns the
## WGS samples also cover, so it is the most any WGS trait could reach given
## the two missing samples. The WGS slope is at that ceiling for practical
## purposes; PC1 sits below it by about 0.15, which is deconvolution error at
## low frequency rather than anything the transform can fix.
##
## Usage:  Rscript scripts/make_baugh_mapping_traits.R [--nnls=dep103|pool100]
## Writes: supplemental_data/phenotypes/baugh_mapping_traits.csv
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse)})

args <- commandArgs(trailingOnly = TRUE)
src <- sub("^--nnls=", "", grep("^--nnls=", args, value = TRUE))
if (!length(src)) src <- "dep103"
stopifnot(src %in% c("dep103", "pool100"))

PH  <- "supplemental_data/phenotypes"
OUT <- file.path(PH, "baugh_mapping_traits.csv")

recipe <- read_csv(file.path(PH, sprintf("baugh_recipe_traits_%s.csv", src)),
                   show_col_types = FALSE)
delta  <- read_csv(file.path(PH, sprintf("baugh_association_traits_%s.csv", src)),
                   show_col_types = FALSE)

tab <- recipe %>%
  select(strain,
         published_slope_baugh = Slope_pub,
         published_pc1_baugh   = PC1_pub,
         slope_nnls, pc1_nnls  = PC1_nnls) %>%
  inner_join(delta %>% select(strain,
                              delta_slope_baugh = slope_baugh,
                              delta_pc1_baugh   = PC1_baugh),
             by = "strain") %>%
  select(strain, published_slope_baugh, published_pc1_baugh,
         delta_slope_baugh, delta_pc1_baugh, slope_nnls, pc1_nnls) %>%
  arrange(strain)

stopifnot(!anyDuplicated(tab$strain), !("N2" %in% tab$strain),
          !anyNA(tab %>% select(-strain)))
write_csv(tab, OUT)

cc <- function(a, b) cor(tab[[a]], tab[[b]], method = "spearman")
message(sprintf("reference: %s   strains: %d\n", src, nrow(tab)))
message("agreement with the published traits (Spearman):")
message(sprintf("  delta_slope_baugh  vs published_slope_baugh  %+.3f",
                cc("delta_slope_baugh", "published_slope_baugh")))
message(sprintf("  delta_pc1_baugh    vs published_pc1_baugh    %+.3f",
                cc("delta_pc1_baugh", "published_pc1_baugh")))
message(sprintf("  slope_nnls         vs published_slope_baugh  %+.3f",
                cc("slope_nnls", "published_slope_baugh")))
message(sprintf("  pc1_nnls           vs published_pc1_baugh    %+.3f",
                cc("pc1_nnls", "published_pc1_baugh")))
message(sprintf("\nwrote %s  (%d strains, 6 traits)", OUT, nrow(tab)))
