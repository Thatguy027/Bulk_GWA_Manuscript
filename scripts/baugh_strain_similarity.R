## Identity-by-state within the Baugh NNLS reference -------------------------
##
##   Rscript scripts/baugh_strain_similarity.R
##     -> supplemental_data/deconvolution/baugh_strain_similarity.tsv
##
## NEEDS THE ARCHIVE. Reads data/baugh/2024bootstrapINPUT.Rdata (31 MB), which
## is Dryad-hosted and not in the deposit. Everything else in this repository
## runs from a clone; this does not. Restore the archive over the repository
## root and run it once; the 102-row table it writes is small enough to commit,
## after which baugh_leakage_vs_similarity.R runs from a clone again.
##
## WHY IT EXISTS. baugh_leakage_vs_similarity.R asks whether the strains the
## deconvolution resolves badly are the genetically similar ones, and needs a
## relatedness measure to ask it. The only such table in the deposit is
## dilution_strain_similarity.tsv, whose nearest-neighbour distances are to the
## closest strain among the 170 DILUTION strains. That is the wrong reference:
## NNLS confusability is a property of the design matrix actually solved, so
## the nearest confounder must be sought in the panel actually used. Only 66 of
## the Baugh panel's 98 strains appear in that table at all, and for only 24 is
## the listed neighbour also in the Baugh panel.
##
## THE RIGHT MATRIX IS THE DESIGN MATRIX. flipped_bootstrap_input[[1]] is the
## genotype matrix the Baugh NNLS was actually solved against -- markers x
## strains, column names "STRAIN_STRAIN". Using it rather than the CeNDR panel
## means the IBS computed here is similarity between the very columns whose
## near-collinearity is the thing under test, on the very markers the solver
## saw. Markers with any missing call are dropped first, exactly as
## Figure1_common.R does before building fullGGp, so the marker set matches the
## one the deconvolution used.
##
## COMPARABILITY WITH THE DILUTION TABLE. The IBS definition, the 500,000-marker
## sample and the seed are the same as in make_experiments_deposit.R, so the two
## tables can be read on the same scale. The columns are the same four, so
## baugh_leakage_vs_similarity.R reads either without changing anything but the
## path.
##
## NOT YET RUN. This script has never been executed against the real input --
## the archive is not on the machine it was written on. It therefore asserts
## its assumptions rather than trusting them: genotype coding, panel membership
## and the range of the result are all checked, and it stops rather than
## writing a table it cannot vouch for.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

BOOT <- "data/baugh/2024bootstrapINPUT.Rdata"
ORD  <- "supplemental_data/deconvolution/baugh_strain_order.txt"
OUT  <- "supplemental_data/deconvolution/baugh_strain_similarity.tsv"
N_MARKERS <- 5e5      # as in make_experiments_deposit.R
SEED      <- 1

if (!file.exists(BOOT))
  stop("need ", BOOT, " -- restore the Dryad archive over the repository root.\n",
       "  See DATA_AVAILABILITY.md. Nothing else in scripts/ requires it.",
       call. = FALSE)
stopifnot(file.exists(ORD))

msg <- function(...) cat("  ", ..., "\n", sep = "")

## ---------------------------------------------------------------------------
## the design matrix, exactly as the deconvolution saw it
## ---------------------------------------------------------------------------
e <- new.env(); load(BOOT, e)
if (!"flipped_bootstrap_input" %in% ls(e))
  stop("flipped_bootstrap_input not found in ", BOOT, call. = FALSE)

gt_all <- e$flipped_bootstrap_input[[1]]
msg("design matrix: ", nrow(gt_all), " markers x ", ncol(gt_all), " strains")

## drop markers with any missing call -- Figure1_common.R does this before
## crossprod, so the retained set is the one the NNLS actually used
keep <- which(apply(gt_all, 1, function(x) sum(is.na(x))) == 0)
gt   <- gt_all[keep, , drop = FALSE]
rm(gt_all); invisible(gc())
msg("complete markers: ", nrow(gt), " of ", length(keep) + 0, " retained")

strains <- sub("_.*$", "", colnames(gt))
if (anyDuplicated(strains))
  msg("NOTE: duplicated strain names in the reference: ",
      paste(unique(strains[duplicated(strains)]), collapse = ", "))

panel <- readLines(ORD)[-1]
if (!setequal(strains, panel))
  stop("the design matrix strains do not match baugh_strain_order.txt.\n",
       "  only in matrix: ", paste(setdiff(strains, panel), collapse = ", "), "\n",
       "  only in order file: ", paste(setdiff(panel, strains), collapse = ", "),
       call. = FALSE)

## ---------------------------------------------------------------------------
## IBS = fraction of markers at which two strains carry the same allele.
##
## With 0/1 coding, mean(x == y) = (both-one + both-zero) / n. both-one is
## crossprod(g); both-zero is n - (colsum_x + colsum_y) + both-one. That is the
## expression below, and it is the same one make_experiments_deposit.R uses, so
## the two similarity tables are on one scale.
##
## The coding is asserted rather than assumed: if this matrix turns out to be
## coded -1/0/1, or to carry dosages, the crossprod identity is wrong and would
## return a plausible-looking but meaningless number.
## ---------------------------------------------------------------------------
set.seed(SEED)
if (nrow(gt) > N_MARKERS) {
  idx <- sort(sample(nrow(gt), N_MARKERS))
  gt  <- gt[idx, , drop = FALSE]
  msg("sampled ", nrow(gt), " markers under seed ", SEED)
} else {
  msg("using all ", nrow(gt), " markers (fewer than the ", N_MARKERS, " cap)")
}

vals <- unique(as.vector(gt[seq_len(min(5000, nrow(gt))), , drop = FALSE]))
if (!all(vals %in% c(0, 1)))
  stop("genotypes are not 0/1 coded -- found ",
       paste(sort(vals), collapse = ", "),
       ".\n  The crossprod IBS identity below assumes 0/1; fix it before ",
       "trusting the output.", call. = FALSE)

## keep only markers that actually vary, since invariant sites inflate IBS
## uniformly and carry no information about who confounds whom
vary <- apply(gt, 1, function(x) length(unique(x)) > 1)
gt   <- gt[vary, , drop = FALSE]
n    <- nrow(gt)
msg("variable markers used: ", n)
if (n < 1e4)
  stop("only ", n, " variable markers -- too few for a stable IBS.", call. = FALSE)

colnames(gt) <- strains
cp   <- crossprod(gt)
cs   <- colSums(gt)
same <- (cp + (n - outer(cs, cs, "+") + cp)) / n
diag(same) <- NA

sim <- tibble(strain     = colnames(gt),
              nn_ibs     = apply(same, 1, max, na.rm = TRUE),
              nn_partner = colnames(same)[apply(same, 1, which.max)],
              mean_ibs   = rowMeans(same, na.rm = TRUE))

## a similarity outside these bounds means the identity above misfired
stopifnot(all(sim$nn_ibs > 0.5), all(sim$nn_ibs <= 1),
          all(sim$mean_ibs > 0.3), all(sim$mean_ibs < sim$nn_ibs))

write_tsv(sim, OUT)
msg("wrote ", OUT, ": ", nrow(sim), " strains")
msg("nearest-neighbour IBS ", sprintf("%.4f-%.4f", min(sim$nn_ibs), max(sim$nn_ibs)))

## ---------------------------------------------------------------------------
## how much did the reference matter? The whole reason for this script.
## ---------------------------------------------------------------------------
old <- "supplemental_data/deconvolution/dilution_strain_similarity.tsv"
if (file.exists(old)) {
  cmp <- sim %>% rename(new_ibs = nn_ibs, new_partner = nn_partner) %>%
    inner_join(read_tsv(old, show_col_types = FALSE) %>%
                 select(strain, old_ibs = nn_ibs, old_partner = nn_partner),
               by = "strain")
  cat("\n== the dilution table as a proxy: how wrong was it? ==\n")
  msg("strains in both tables: ", nrow(cmp))
  msg("same nearest neighbour named: ", sum(cmp$new_partner == cmp$old_partner),
      " of ", nrow(cmp))
  msg("Spearman of the two nn_ibs columns: ",
      sprintf("%+.3f", cor(cmp$new_ibs, cmp$old_ibs, method = "spearman")))
  msg("median change in nn_ibs: ",
      sprintf("%+.4f", median(cmp$new_ibs - cmp$old_ibs)))
  cat("\n  the ten strains whose nearest neighbour changed most:\n")
  print(as.data.frame(cmp %>% mutate(d = new_ibs - old_ibs) %>%
          arrange(desc(abs(d))) %>% head(10) %>%
          transmute(strain, old_partner, old_ibs = round(old_ibs, 4),
                    new_partner, new_ibs = round(new_ibs, 4),
                    change = round(d, 4))), row.names = FALSE)
  cat("\n  Now rerun: Rscript scripts/baugh_leakage_vs_similarity.R\n",
      "  It prefers this table when it exists, reports all 98 strains rather\n",
      "  than 66, and drops the well-posed/proxy split, which only existed\n",
      "  because the predictor was borrowed. Its pinned numbers WILL fail --\n",
      "  that is the point of pinning them. Update them to what it prints.\n",
      sep = "")
}
