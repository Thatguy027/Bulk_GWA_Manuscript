## The Baugh L1 mapping table ------------------------------------------------
##
## One file to run association scans from: two platforms x two constructions x
## two traits, plus the published values themselves.
##
##   published_slope_baugh  eLife Slope, verbatim             MIP-seq
##   published_pc1_baugh    eLife PC1, verbatim               MIP-seq
##   delta_slope_baugh      delta construction                MIP-seq
##   delta_pc1_baugh        delta construction                MIP-seq
##   slope_nnls             log-ratio construction, floored   pooled WGS
##   pc1_nnls               log-ratio construction, floored   pooled WGS
##   delta_slope_nnls       delta construction                pooled WGS
##   delta_pc1_nnls         delta construction                pooled WGS
##
## THE TWO CONSTRUCTIONS DIFFER ONLY IN THE TRANSFORM. Both build a strains x
## replicate-day matrix, run prcomp(scale, center) and take PC1 = x[,1], and
## both fit the slope of that same value on day excluding day 17.
##
##   log-ratio   v = log2(f_day / f_baseline)   the published recipe. Requires
##               a floor because non-negative least squares returns exact zeros;
##               1/(4n) is used. Reproduces the published PC1 at Spearman 1.0000
##               and the published Slope at 0.9959 when run on MIP-seq.
##   delta       v = f_day - f_day1             no logarithm, so no floor is
##               needed and none is applied.
##
## THE REFERENCE SAMPLE MATTERS FAR MORE THAN THE TRANSFORM, which is why both
## constructions here reference the baseline. Holding the transform at identity
## and switching only the reference:
##
##                        ref = baseline        ref = day 1
##   MIP-seq PC1               0.903               0.686
##   WGS PC1                   0.815               0.697
##
## against a spread of 0.814 to 0.826 across every transform tried at a fixed
## baseline reference. Slope is nearly indifferent (0.889 either way for MIP).
## An earlier version of this table referenced day 1 for the delta construction
## and so reported a delta PC1 of 0.686 that was handicapped by the reference
## rather than by the missing logarithm.
##
## Agreement with the published traits, dep103 reference, n = 99:
##
##                              vs published Slope   vs published PC1
##   delta, MIP-seq                   0.889                0.903
##   delta, pooled WGS                0.869                0.815
##   log-ratio, pooled WGS            0.884                0.814
##   ceiling (MIP, matched cols)      0.983                0.961
##
## So with the reference fixed the two constructions are equivalent on pooled
## WGS -- 0.815 against 0.814 for PC1 -- and the delta one gets there without a
## floor or any tuning parameter, which is the reason to prefer it. The MIP
## delta PC1 reaching 0.903 where the WGS one reaches 0.815 locates the
## remaining shortfall in the deconvolution rather than in the trait definition.
##
## VARIANCE STABILISATION WAS TESTED AND IS NOT ADOPTED. asinh(f/c), which is
## linear near zero and logarithmic above it and so needs no floor, peaks at
## c = 0.01 -- about one equal share -- giving Slope 0.901 against 0.884 and
## PC1 0.826 against 0.814. A paired bootstrap over strains puts those gains at
## +0.017 (95% CI -0.002 to +0.038) and +0.012 (-0.020 to +0.044), so neither is
## distinguishable at n = 99, and c was chosen against the target it is being
## scored on. sqrt and arcsine-sqrt are slightly worse than identity. The
## finding is recorded rather than used.
##
## NOTE ON A CHANGED DEFINITION. An earlier version of this table built
## delta_pc1_baugh from the legacy exploratory PCA -- a samples x strains
## rotation on raw frequencies -- which agrees with the published PC1 at only
## 0.331. It is replaced here by the delta construction above, which is the same
## orientation and structure as the published recipe and differs from it only in
## the transform, so the four PC1 columns are now comparable to one another.
##
## Usage:  Rscript scripts/make_baugh_mapping_traits.R [--nnls=dep103|pool100]
## Writes: supplemental_data/phenotypes/baugh_mapping_traits.csv
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse)})

args <- commandArgs(trailingOnly = TRUE)
src <- sub("^--nnls=", "", grep("^--nnls=", args, value = TRUE))
if (!length(src)) src <- "dep103"
stopifnot(src %in% c("dep103", "pool100"))

DEC <- "supplemental_data/deconvolution"
PH  <- "supplemental_data/phenotypes"
OUT <- file.path(PH, "baugh_mapping_traits.csv")

ols <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); if (sum(ok) < 2) return(NA_real_)
  x <- x[ok]; y <- y[ok]; xc <- x - mean(x); s <- sum(xc^2)
  if (s == 0) NA_real_ else sum(xc * y) / s
}

nn <- read_tsv(file.path(DEC, sprintf("baugh_nnls_%s_with_mipseq.tsv.gz", src)),
               show_col_types = FALSE)
FLOOR <- 1 / (4 * n_distinct(nn$strain[!is.na(nn$frq)]))

## one construction, one frequency column
build <- function(col, mode, do_floor = FALSE) {
  l <- nn %>% rename(f = all_of(col)) %>% filter(!is.na(f)) %>%
    mutate(f = if (do_floor) pmax(f, FLOOR) else f,
           day2 = ifelse(baseline, "BL", paste0("d", day))) %>%
    select(strain, rep = replicate, day2, dnum = day, f)
  ## Both constructions reference the BASELINE sample. That is what makes them
  ## a clean A/B on the transform alone, and it matters far more than the
  ## transform does -- see the header.
  ref <- l %>% filter(day2 == "BL") %>% select(strain, rep, r = f)
  w <- l %>% filter(day2 != "BL") %>% left_join(ref, by = c("strain", "rep")) %>%
    mutate(v = if (mode == "logratio") log2(f / r) else f - r,
           colk = paste0(rep, "_", day2)) %>%
    filter(is.finite(v))
  m <- w %>% select(strain, colk, v) %>%
    pivot_wider(names_from = colk, values_from = v) %>%
    column_to_rownames("strain") %>% as.matrix()
  m <- m[, colSums(is.na(m)) < 0.1 * nrow(m), drop = FALSE]
  m <- m[complete.cases(m), , drop = FALSE]
  m <- m[, apply(m, 2, sd) > 0, drop = FALSE]   # the delta day-1 column is all zero
  p <- prcomp(m, scale. = TRUE, center = TRUE)
  message(sprintf("  %-14s %-9s %2d columns, %d strains", col, mode, ncol(m), nrow(m)))
  tibble(strain = rownames(m), pc1 = p$x[, 1]) %>%
    left_join(w %>% filter(dnum != 17) %>% group_by(strain) %>%
                summarise(sl = ols(dnum, v), .groups = "drop"), by = "strain")
}

message("building:")
d_mip  <- build("published_frq", "delta")
d_nnls <- build("frq",           "delta")
l_nnls <- build("frq",           "logratio", do_floor = TRUE)

pub <- read_tsv(file.path(PH, "baugh_published_traits.txt"), show_col_types = FALSE)

## every PC1 oriented to the published one so effect signs are comparable
orient <- function(t, ref) {
  j <- inner_join(t, ref, by = "strain")
  if (cor(j$pc1, j$PC1, method = "spearman") < 0) t$pc1 <- -t$pc1
  t
}
d_mip <- orient(d_mip, pub); d_nnls <- orient(d_nnls, pub); l_nnls <- orient(l_nnls, pub)

tab <- pub %>%
  transmute(strain, published_slope_baugh = Slope, published_pc1_baugh = PC1) %>%
  inner_join(d_mip  %>% transmute(strain, delta_slope_baugh = sl, delta_pc1_baugh = pc1), by = "strain") %>%
  inner_join(l_nnls %>% transmute(strain, slope_nnls        = sl, pc1_nnls        = pc1), by = "strain") %>%
  inner_join(d_nnls %>% transmute(strain, delta_slope_nnls  = sl, delta_pc1_nnls  = pc1), by = "strain") %>%
  filter(strain != "N2") %>% arrange(strain)

stopifnot(!anyDuplicated(tab$strain), !anyNA(tab %>% select(-strain)))
write_csv(tab, OUT)

cc <- function(a, b) cor(tab[[a]], tab[[b]], method = "spearman")
message(sprintf("\nreference %s, n = %d\n", src, nrow(tab)))
message("                     vs published Slope   vs published PC1")
for (p in list(c("delta_slope_baugh", "delta_pc1_baugh",  "delta, MIP-seq   "),
               c("delta_slope_nnls",  "delta_pc1_nnls",   "delta, WGS       "),
               c("slope_nnls",        "pc1_nnls",         "log-ratio, WGS   ")))
  message(sprintf("  %s        %+.3f             %+.3f", p[3],
                  cc(p[1], "published_slope_baugh"), cc(p[2], "published_pc1_baugh")))
message(sprintf("\nwrote %s  (%d strains, 8 traits)", OUT, nrow(tab)))
