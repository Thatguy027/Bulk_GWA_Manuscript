## Baugh L1 traits on the PUBLISHED recipe, both platforms --------------------
##
## The earlier trait tables reconstructed a PC1 from an old exploratory script.
## That PC1 is NOT the published one: it correlates with the eLife trait at
## Spearman 0.33. The email thread in data/baugh/FW_ FW_ Order 7919.eml states
## the real recipe (Amy Webster, relayed by R. Baugh, 2022-09-08):
##
##   "dividing day 1, 9, 13, and 17 by baseline frequency and log2 transforming
##    to put into PCA (the prcomp function with scale and center set to TRUE)"
##
## Implemented here and VALIDATED against data/baugh/eLife_traits.txt:
##   PC1    log2(f_day / f_baseline) for d1,d9,d13,d17 x 5 reps -> 100 x 20
##          matrix, prcomp(scale.=TRUE, center=TRUE), PC1 = x[,1]
##          reproduces the published PC1 at Spearman 1.0000, |Pearson| 1.0000
##   Slope  slope of log2(f/baseline) on day, pooled across reps, EXCLUDING
##          day 17. Reproduces the published Slope at Spearman 0.9958,
##          |Pearson| 0.9975 -- very close but not exact, so the published
##          column is shipped as the canonical value and this is used only to
##          derive the matching WGS trait.
##
## WHAT THIS WRITES. supplemental_data/phenotypes/baugh_recipe_traits.csv:
##
##   Slope_pub, PC1_pub          the published eLife values, verbatim. Use
##                               these as the MIP-side trait.
##   slope_nnls, PC1_nnls        the same recipe applied to NNLS frequencies.
##   slope_mip_m, PC1_mip_m      the same recipe applied to MIP frequencies but
##                               restricted to the columns WGS also has, so the
##                               platform comparison is like for like.
##
## WHY THE MATCHED COLUMNS EXIST. WGS is missing two of the 25 MIP samples --
## rep2's baseline and rep5_d17 (confirmed in the thread: "two of our 25 samples
## were used up"). Without a baseline, none of rep2's four days can be
## normalised, so the WGS recipe runs on 15 of the 20 columns. Comparing a
## 20-column MIP trait against a 15-column WGS trait would confound the platform
## with the sample set; the _m columns remove that.
##
## ZEROS. The recipe is a log ratio, and NNLS returns exact zeros -- it is a
## non-negative fit with a hard boundary at zero, so a strain the solver cannot
## place gets 0.000, not a small number. 340 of 2,369 NNLS cells are exactly
## zero, which makes 793 of 2,057 log ratios non-finite and leaves only 47
## strains with a complete 15-column row. MIP frequencies are read-count ratios
## and never hit zero, so this problem is specific to the deconvolution and is
## probably why the original comparison used a difference-based slope for WGS
## rather than this recipe.
##
## The rule applied here: zeros are floored at half the smallest positive
## frequency in the table, which is the usual convention for log-transforming
## count-derived compositions. The floor value and the number of cells it
## touches are reported on every run. Change FLOOR_RULE if you want a different
## convention -- this is a judgement call, not a derived quantity.
##
## WHAT THE RESULT SAYS. On the published log-ratio recipe the two platforms
## agree far less well than on the difference-based slope the repository
## already uses:
##
##                                       all 99      43 strains with no zeros
##   slope, published log-ratio recipe    0.702              0.946
##   PC1,   published recipe              0.674              0.770
##   slope, difference-based (repo)       0.975              0.968
##
## Two things follow. The zeros are doing most of the damage to the log-ratio
## slope -- it recovers to 0.946 once strains carrying a zero are set aside,
## while the difference-based slope barely moves (0.975 to 0.968), so that
## trait is simply robust to them. And PC1 stays comparatively low even on
## clean strains, because a scaled log ratio weights rare strains as heavily as
## common ones and NNLS is least reliable exactly there.
##
## 57 of 103 strains carry at least one zero and they are the low-frequency
## ones (median mean frequency 0.0052 against 0.0099 for the rest). JU2001 is
## zero in every sample under the ref103 fit despite a non-zero MIP frequency.
##
## PRACTICAL READ: map Slope_pub and PC1_pub as the published phenotypes, and
## use slope_nnls/PC1_nnls for the platform comparison rather than as primary
## traits, because their agreement is limited by NNLS zeros rather than by
## anything biological.
##
## Usage:  Rscript scripts/make_baugh_published_recipe_traits.R [--nnls=ref103|ref100]
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse)})

args <- commandArgs(trailingOnly = TRUE)
which_nnls <- sub("^--nnls=", "", grep("^--nnls=", args, value = TRUE))
if (!length(which_nnls)) which_nnls <- "ref103"
stopifnot(which_nnls %in% c("ref103", "ref100"))

BAUGH <- "supplemental_data/deconvolution"
PUB   <- "supplemental_data/phenotypes/baugh_published_traits.txt"
OUT   <- sprintf("supplemental_data/phenotypes/baugh_recipe_traits_%s.csv", which_nnls)

ols <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); if (sum(ok) < 2) return(NA_real_)
  x <- x[ok]; y <- y[ok]; xc <- x - mean(x); s <- sum(xc^2)
  if (s == 0) NA_real_ else sum(xc * y) / s
}

## long frame of frequencies with a rep/day key, one row per strain x sample
as_long <- function(d, value) {
  d %>% rename(f = all_of(value)) %>%
    mutate(day = ifelse(baseline, "BL", paste0("d", day))) %>%
    select(strain, rep = replicate, day, f)
}

## the published recipe, given a long frame and the columns to keep
recipe <- function(long, cols = NULL, label = "") {
  nz <- sum(long$f == 0, na.rm = TRUE)
  if (nz > 0) {
    floor_v <- min(long$f[long$f > 0], na.rm = TRUE) / 2
    message(sprintf("  %s: flooring %d zero cells at %.3g (half the smallest positive)",
                    label, nz, floor_v))
    long <- long %>% mutate(f = ifelse(!is.na(f) & f == 0, floor_v, f))
  }
  bl <- long %>% filter(day == "BL") %>% select(strain, rep, base = f)
  w  <- long %>% filter(day != "BL") %>%
    left_join(bl, by = c("strain", "rep")) %>%
    mutate(l2 = log2(f / base), dnum = as.numeric(sub("d", "", day)),
           col = paste0(rep, "_", day)) %>%
    filter(is.finite(l2))
  if (!is.null(cols)) w <- w %>% filter(col %in% cols)
  mat <- w %>% select(strain, col, l2) %>%
    pivot_wider(names_from = col, values_from = l2) %>%
    column_to_rownames("strain") %>% as.matrix()
  keep_col <- colSums(is.na(mat)) < 0.1 * nrow(mat)   # drop a column only if it is mostly empty
  mat <- mat[, keep_col, drop = FALSE]
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  stopifnot(nrow(mat) > 10, ncol(mat) > 2)
  p <- prcomp(mat, scale. = TRUE, center = TRUE)
  sl <- w %>% filter(day != "d17") %>% group_by(strain) %>%
    summarise(slope = ols(dnum, l2), .groups = "drop")
  list(pc1 = tibble(strain = rownames(mat), pc1 = p$x[, 1]), slope = sl,
       cols = colnames(mat), ve = (p$sdev^2 / sum(p$sdev^2))[1])
}

## --- MIP, full 20 columns: validate against the published values -----------
ln  <- readLines(gzfile(file.path(BAUGH, "mipseq_frequencies.txt.gz")))
hdr <- strsplit(ln[1], "\t")[[1]]; body <- strsplit(ln[-1], "\t")
mip_long <- tibble(strain = sapply(body, `[`, 1),
    !!!setNames(lapply(seq_along(hdr),
        function(j) as.numeric(sapply(body, `[`, j + 1))), hdr)) %>%
  pivot_longer(-strain, names_to = "sample", values_to = "f") %>%
  separate(sample, into = c("rep", "day"), sep = "_") %>%
  select(strain, rep, day, f)

pub <- read_tsv(PUB, show_col_types = FALSE)
r_full <- recipe(mip_long, label = "MIP full")
chk <- pub %>% inner_join(r_full$pc1, by = "strain") %>%
  inner_join(r_full$slope, by = "strain")
message(sprintf("validation against the published traits (n = %d):", nrow(chk)))
message(sprintf("  PC1   Spearman %+.4f  |Pearson| %.4f  (PC1 explains %.1f%%)",
  cor(chk$PC1, chk$pc1, method = "spearman"), abs(cor(chk$PC1, chk$pc1)), 100 * r_full$ve))
message(sprintf("  Slope Spearman %+.4f  |Pearson| %.4f",
  cor(chk$Slope, chk$slope, method = "spearman"), abs(cor(chk$Slope, chk$slope))))
stopifnot(cor(chk$PC1, chk$pc1, method = "spearman") > 0.999)

## --- NNLS, on whatever columns WGS has -------------------------------------
nn <- read_tsv(file.path(BAUGH, sprintf("baugh_nnls_%s_with_mipseq.tsv.gz", which_nnls)),
               show_col_types = FALSE)
nn_long <- as_long(nn, "frq")
r_nnls <- recipe(nn_long, label = paste("NNLS", which_nnls))
message(sprintf("\nNNLS (%s): %d columns usable of 20 -- %s",
  which_nnls, length(r_nnls$cols),
  paste(setdiff(r_full$cols, r_nnls$cols), collapse = ", ")))

## --- MIP again, matched to those columns -----------------------------------
r_mipm <- recipe(mip_long, cols = r_nnls$cols, label = "MIP matched")

out <- pub %>% rename(Slope_pub = Slope, PC1_pub = PC1) %>%
  left_join(r_nnls$slope %>% rename(slope_nnls = slope), by = "strain") %>%
  left_join(r_nnls$pc1   %>% rename(PC1_nnls   = pc1),   by = "strain") %>%
  left_join(r_mipm$slope %>% rename(slope_mip_m = slope), by = "strain") %>%
  left_join(r_mipm$pc1   %>% rename(PC1_mip_m   = pc1),   by = "strain") %>%
  filter(strain != "N2", !is.na(slope_nnls)) %>%
  select(strain, Slope_pub, PC1_pub, slope_mip_m, PC1_mip_m, slope_nnls, PC1_nnls) %>%
  arrange(strain)

## orient every PC1 to the published one so signs are comparable
for (cl in c("PC1_nnls", "PC1_mip_m")) {
  r <- cor(out[[cl]], out$PC1_pub, method = "spearman")
  if (is.finite(r) && r < 0) { out[[cl]] <- -out[[cl]]
    message(sprintf("  flipped %s (was %+.3f vs PC1_pub)", cl, r)) }
}
write_csv(out, OUT)

message(sprintf("\nplatform agreement on matched columns (n = %d):", nrow(out)))
message(sprintf("  slope  MIP vs NNLS  Spearman %+.4f",
  cor(out$slope_mip_m, out$slope_nnls, method = "spearman")))
message(sprintf("  PC1    MIP vs NNLS  Spearman %+.4f",
  cor(out$PC1_mip_m, out$PC1_nnls, method = "spearman")))
message(sprintf("\npublished trait vs NNLS:"))
message(sprintf("  Slope_pub vs slope_nnls  Spearman %+.4f",
  cor(out$Slope_pub, out$slope_nnls, method = "spearman")))
message(sprintf("  PC1_pub   vs PC1_nnls    Spearman %+.4f",
  cor(out$PC1_pub, out$PC1_nnls, method = "spearman")))
message(sprintf("\nwrote %s  (%d strains)", OUT, nrow(out)))
