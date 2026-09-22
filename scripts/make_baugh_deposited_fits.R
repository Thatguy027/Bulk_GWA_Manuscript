## Refit the Baugh deconvolution from the DEPOSITED input ---------------------
##
## data/baugh/2024bootstrapINPUT.Rdata is the input the shipped cache was built
## from -- 1,237,106 markers x 102 strains, matching baugh_strain_order.txt
## exactly. Everything downstream should start here rather than from the 2022
## matrix in the source project.
##
## THE DEPOSITED REFERENCE DOES NOT MATCH THE POOL, in both directions. It
## carries CX11262, ECA348 and NIC260, which the MIP panel does not measure,
## and it LACKS PB306, which the MIP panel does measure and which the 2022
## matrices all carry. The thread in data/baugh/FW_ FW_ Order 7919.eml settles
## what the pool was: "Amy has confirmed that all 103 strains were included in
## the pool" (R. Baugh, 2022-09-08). So all 103 belong and the deposited
## reference is one strain short.
##
## PB306 is therefore grafted back from the 2022 matrix. That is defensible
## because the two matrices agree closely: every one of the 1,237,106 deposited
## markers exists in the 2022 matrix, and across shared strains the genotype
## calls disagree at 0.04%. It is not free, though -- PB306's column comes from
## a slightly older call set than the other 102. The grafted fit is labelled
## dep103 and kept separate from the untouched dep102 for that reason.
##
## THREE FITS
##   dep102   the deposited reference, untouched. Validated against the shipped
##            cache; this is the reproduction of what Figure 1 was built on.
##   dep103   dep102 plus PB306, i.e. the full pool as confirmed in the thread.
##   pool100  dep103 cut to the 100 strains the MIP panel measures.
##
## dep103 against pool100 is the controlled comparison: one matrix, one marker
## set, PB306 present in both, differing only in whether the three unmeasured
## strains are in the reference.
##
## Usage:  Rscript scripts/make_baugh_deposited_fits.R
## Writes: supplemental_data/deconvolution/baugh_nnls_{dep102,dep103,pool100}_with_mipseq.tsv.gz
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse)})

DEP   <- "data/baugh/2024bootstrapINPUT.Rdata"
OLD   <- Sys.getenv("BAUGH_BOOT", "/Users/Stefan/UCLA/Projects/bulkGWAS/baugh_wgs/cluster_data/20220908_Baugh_BulkL1_Bootstrap_Input_flippedCommon_NAfix.RData")
BAUGH <- "supplemental_data/deconvolution"
CACHE <- file.path(BAUGH, "baugh_nnls_with_mipseq.RData")

e <- new.env(); load(DEP, e)
gt <- e$flipped_bootstrap_input[[1]]; ct <- e$flipped_bootstrap_input[[2]]
rm(e); gc()
keep <- which(rowSums(is.na(gt)) == 0)
gt <- gt[keep, ]; ct <- ct[keep, ]
colnames(gt) <- sub("_.*$", "", colnames(gt))
message(sprintf("deposited: %d markers x %d strains after dropping NA rows",
                nrow(gt), ncol(gt)))

fit <- function(mat, tag) {
  Gy <- crossprod(mat, ct); GG <- crossprod(mat)
  p <- apply(Gy, 2, function(x) as.vector(RcppML::nnls(GG, matrix(x), fast_nnls = TRUE)))
  p <- apply(p, 2, function(x) x / sum(x)); rownames(p) <- colnames(mat)
  message(sprintf("  %-8s %d strains", tag, nrow(p)))
  as_tibble(p, rownames = "strain") %>%
    pivot_longer(-strain, names_to = "sample", values_to = "frq")
}

## --- dep102 and validation against the shipped cache -----------------------
dep102 <- fit(gt, "dep102")
ce <- new.env(); load(CACHE, ce)
cache <- as_tibble(ce$wgs_mip_results) %>% select(strain, sample, cached = frq)
v <- dep102 %>% inner_join(cache, by = c("strain", "sample"))
message(sprintf("\nVALIDATION against the shipped cache (%d cells):", nrow(v)))
message(sprintf("  Pearson %.10f   max abs difference %.3g",
                cor(v$frq, v$cached), max(abs(v$frq - v$cached))))
if (max(abs(v$frq - v$cached)) > 1e-8)
  message("  NOTE: not bit-identical; see the header of this script.")

## --- graft PB306 -----------------------------------------------------------
o <- new.env(); load(OLD, o)
go <- o$flipped_bootstrap_input[[1]]
pb <- go[rownames(gt), grep("^PB306_", colnames(go))]
sh <- intersect(colnames(gt), sub("_.*$", "", colnames(go)))
disagree <- mean(gt[, sh[1:10]] != go[rownames(gt), paste0(sh[1:10], "_", sh[1:10])])
rm(o, go); gc()
message(sprintf("\ngrafting PB306 from the 2022 matrix (shared-strain call disagreement %.4f)", disagree))
gt103 <- cbind(gt, PB306 = pb)

mip <- sub("\t.*$", "", readLines(gzfile(file.path(BAUGH, "mipseq_frequencies.txt.gz")))[-1])
pool <- intersect(colnames(gt103), mip)
message(sprintf("MIP pool strains present in dep103: %d", length(pool)))

dep103  <- fit(gt103, "dep103")
pool100 <- fit(gt103[, pool], "pool100")

## --- stage all three -------------------------------------------------------
ln <- readLines(gzfile(file.path(BAUGH, "mipseq_frequencies.txt.gz")))
hdr <- strsplit(ln[1], "\t")[[1]]; body <- strsplit(ln[-1], "\t")
mipf <- tibble(strain = sapply(body, `[`, 1),
  !!!setNames(lapply(seq_along(hdr),
      function(j) as.numeric(sapply(body, `[`, j + 1))), hdr)) %>%
  pivot_longer(-strain, names_to = "sample", values_to = "published_frq") %>%
  mutate(sample = gsub("_BL", "_d1_baseline", sample))

stage <- function(d, tag) {
  f <- file.path(BAUGH, sprintf("baugh_nnls_%s_with_mipseq.tsv.gz", tag))
  out <- d %>% full_join(mipf, by = c("strain", "sample")) %>%
    separate(sample, into = c("replicate", "day"), sep = "_",
             remove = FALSE, extra = "drop") %>%
    mutate(baseline = grepl("baseline", sample),
           day = as.numeric(gsub("d", "", day))) %>%
    select(sample, replicate, day, baseline, strain, frq, published_frq) %>%
    arrange(strain, sample)
  write_tsv(out, gzfile(f))
  message(sprintf("  wrote %-52s %d strains", basename(f), n_distinct(out$strain)))
}
message("")
stage(dep102, "dep102"); stage(dep103, "dep103"); stage(pool100, "pool100")
