## Re-fit the Baugh NNLS with and without a pool-restricted reference --------
##
## The genotype reference is not restricted to the strains the MIP panel
## measured, so it solves for 103 strains when the pool definition has 100.
## This fits both and saves the two frequency tables so they can be scored
## against MIP. The finding is written up in the header of
## scripts/DIAG_baugh_platform_discrepancy.R.
##
## REQUIRES a file outside this repository, from the source project:
##   /Users/Stefan/UCLA/Projects/bulkGWAS/baugh_wgs/cluster_data/
##     20220908_Baugh_BulkL1_Bootstrap_Input_flippedCommon_NAfix.RData
## That is a near relative of the deposited 2024bootstrapINPUT.Rdata rather
## than the same file -- it carries 103 strains including PB306, where the
## shipped cache has 102 and no PB306 -- so the absolute numbers here will not
## match the deposited fit exactly. The comparison between the two fits is
## internal to this script and unaffected.
##
## Writes refit.rds (full and restricted frequency tables) to the scratch path
## set below; nothing in the repository is modified.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse)})
S <- "/private/tmp/claude-501/-Users-Stefan-UCLA-Projects-RNAi-Manuscript-Bulk-GWA-Manuscript/914ad822-1756-45db-a32c-68329e5f52a4/scratchpad"
setwd("/Users/Stefan/UCLA/Projects/bulkGWAS/baugh_wgs/cluster_data")
e <- new.env(); load("20220908_Baugh_BulkL1_Bootstrap_Input_flippedCommon_NAfix.RData", e)
gt0 <- e$flipped_bootstrap_input[[1]]; ct <- e$flipped_bootstrap_input[[2]]
keep <- which(rowSums(is.na(gt0)) == 0)
gt0 <- gt0[keep, ]; ct <- ct[keep, ]
strains <- sub("_.*$", "", colnames(gt0))
cat("markers", nrow(gt0), " strains", ncol(gt0), "\n")

mip <- sub("\t.*$", "", readLines(gzfile("/Users/Stefan/UCLA/Projects/RNAi/Manuscript/Bulk_GWA_Manuscript/supplemental_data/deconvolution/mipseq_frequencies.txt.gz"))[-1])
inpool <- strains %in% mip
cat("reference strains also measured by MIP:", sum(inpool), " extra:", sum(!inpool),
    " ->", paste(strains[!inpool], collapse=" "), "\n\n")

fit <- function(cols, tag) {
  g <- gt0[, cols, drop = FALSE]
  Gy <- crossprod(g, ct); GG <- crossprod(g)
  p <- apply(Gy, 2, function(x) as.vector(RcppML::nnls(GG, matrix(x), fast_nnls = TRUE)))
  p <- apply(p, 2, function(x) x / sum(x))
  rownames(p) <- strains[cols]
  cat(tag, ": fitted", nrow(p), "strains\n")
  as_tibble(p, rownames = "strain") %>% pivot_longer(-strain, names_to="sample", values_to="frq")
}
full <- fit(seq_along(strains), "FULL     ")
rest <- fit(which(inpool),      "RESTRICTED")
saveRDS(list(full=full, rest=rest, strains=strains, inpool=inpool), file.path(S,"refit.rds"))
cat("done\n")
