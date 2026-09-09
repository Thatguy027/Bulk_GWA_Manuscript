## Stage the genome-wide parental allele frequencies of both crosses ---------
##
##   Rscript scripts/make_cross_af_tables.R
##     data/cross_experiments/N2-XZ_export/afd/*.afd.tsv.gz
##     data/cross_experiments/JU1793-JU2466_export/afd/*.afd.tsv.gz
##     -> supplemental_data/mapping/cross_af_N2xXZ1516.tsv.gz      (~10.0 MB)
##        supplemental_data/mapping/cross_af_JU1793xJU2466.tsv.gz  (~1.2 MB)
##        supplemental_data/mapping/cross_af_samples.tsv
##
## WHY THIS EXISTS
## The thinned bundle carries LOD and contrast.beta for all 18 contrasts but
## not the allele counts they were computed from, so nothing in a clone could
## state a parental frequency. The exports do carry them, in two forms:
##
##   plot_data/*_plot_DF.tsv.gz   15 pairwise contrast files, 780 MB, in which
##                                each sample's 11 columns are repeated once
##                                per contrast it appears in -- five times over
##   afd/*.afd.tsv.gz             one file per sample, 158 MB, no duplication
##
## Both reduce to the same 20 columns of counts. Building from afd/ and keeping
## only chrom, position and each sample's p1/p2 gives 11.2 MB for both crosses
## at full marker resolution -- small enough to commit, so every figure and
## table can read a parental frequency from a clone.
##
## p1 IS ALWAYS PARENT 1, AND PARENT 1 IS THE FIRST NAME IN THE CROSS
## N2 x XZ1516: p1 = N2.  JU1793 x JU2466: p1 = JU1793, the resistant parent.
## The assertion below is the one from load_parent_freq() in Figure3_common.R:
## pos-1 selection must raise the resistant parent's frequency on the right arm
## of chromosome III relative to the HT115 control. If p1/p2 were ever swapped
## upstream, that comparison inverts and this script stops.
##
## Marker sets are identical across the samples of a cross -- asserted, not
## assumed -- so the join is a cbind and the row order is the exports' own.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

SD  <- "supplemental_data/mapping"
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

CROSS <- tibble::tribble(
  ~cross,          ~dir,
  "N2xXZ1516",     "data/cross_experiments/N2-XZ_export",
  "JU1793xJU2466", "data/cross_experiments/JU1793-JU2466_export")

missing <- CROSS$dir[!dir.exists(file.path(CROSS$dir, "afd"))]
if (length(missing))
  stop("no afd/ directory under: ", paste(missing, collapse = ", "),
       "\n  these live in the Dryad archive, not the repository", call. = FALSE)

dir.create(SD, recursive = TRUE, showWarnings = FALSE)

## --- the sample sheet, both crosses in one table ---------------------------
sheet <- purrr::pmap_dfr(CROSS, function(cross, dir) {
  fread(file.path(dir, "sample_sheet.tsv"), colClasses = "character") %>%
    as_tibble() %>%
    mutate(cross = cross,
           afd_id = NA_character_, .before = 1)
})

## --- one table per cross ---------------------------------------------------
build_cross <- function(cross, dir) {
  fs <- sort(Sys.glob(file.path(dir, "afd", "*.afd.tsv.gz")))
  msg(cross, ": ", length(fs), " samples")

  ids <- NULL
  key <- NULL
  cols <- list()
  for (f in fs) {
    samp <- sub("\\.afd\\.tsv\\.gz$", "", basename(f))
    d <- fread(f, select = c("chrom", "physical.position", "ID",
                             paste0("p1_", samp), paste0("p2_", samp)),
               colClasses = list(character = c("chrom", "ID")))
    if (is.null(ids)) {
      ids <- d$ID
      key <- d[, .(chrom, pos = physical.position)]
    } else if (!identical(ids, d$ID)) {
      stop("marker set of ", samp, " differs from ", basename(fs[1]),
           " -- the cbind join is not valid for this export", call. = FALSE)
    }
    cols[[paste0("p1_", samp)]] <- d[[paste0("p1_", samp)]]
    cols[[paste0("p2_", samp)]] <- d[[paste0("p2_", samp)]]
  }
  out <- cbind(key, as.data.table(cols))
  msg("    ", nrow(out), " markers x ", length(fs), " samples")
  out
}

## the resistant-parent assertion, in the one cross that has a resistant parent
## named in the manuscript: JU1793 is resistant, so the pos-1 pool must carry
## more of it than the HT115 control on the right arm of chromosome III
assert_parent_orientation <- function(d, cross) {
  if (cross != "JU1793xJU2466") return(invisible(NULL))
  nm <- names(d)
  p1_ht <- grep("^p1_.*HT115g$", nm, value = TRUE)
  p2_ht <- grep("^p2_.*HT115g$", nm, value = TRUE)
  p1_ps <- grep("^p1_.*POS1g$",  nm, value = TRUE)
  p2_ps <- grep("^p2_.*POS1g$",  nm, value = TRUE)
  if (length(c(p1_ht, p2_ht, p1_ps, p2_ps)) != 4)
    stop("expected one HT115g and one POS1g sample in ", cross, call. = FALSE)
  r <- d[chrom == "III" & pos >= 13e6]
  f_ps <- sum(r[[p1_ps]]) / sum(r[[p1_ps]] + r[[p2_ps]])
  f_ht <- sum(r[[p1_ht]]) / sum(r[[p1_ht]] + r[[p2_ht]])
  msg("    JU1793 frequency on III right arm (>= 13 Mb): pos-1 ",
      sprintf("%.3f", f_ps), " vs HT115 ", sprintf("%.3f", f_ht))
  if (!(f_ps > f_ht))
    stop("the pos-1 pool is not enriched for JU1793 on the right arm of\n",
         "  chromosome III (", sprintf("%.3f vs %.3f", f_ps, f_ht), ").\n",
         "  JU1793 is the RESISTANT parent, so pos-1 selection must raise its\n",
         "  frequency relative to the HT115 control. p1/p2 have most likely\n",
         "  been swapped upstream -- see the export README.", call. = FALSE)
  invisible(NULL)
}

written <- purrr::pmap_dfr(CROSS, function(cross, dir) {
  d <- build_cross(cross, dir)
  assert_parent_orientation(d, cross)
  dest <- file.path(SD, paste0("cross_af_", cross, ".tsv.gz"))
  fwrite(d, dest, sep = "\t", compress = "gzip")
  samps <- grep("^p1_", names(d), value = TRUE) %>% sub("^p1_", "", .)
  tibble::tibble(cross = cross, dest = dest,
                 markers = nrow(d), samples = length(samps),
                 MB = round(file.size(dest) / 1048576, 1))
})

## --- the sample key, so a p1_/p2_ column pair can be named ----------------
key_tbl <- purrr::pmap_dfr(CROSS, function(cross, dir) {
  fs <- sort(Sys.glob(file.path(dir, "afd", "*.afd.tsv.gz")))
  sh <- fread(file.path(dir, "sample_sheet.tsv"), colClasses = "character")
  tibble::tibble(afd_id = sub("\\.afd\\.tsv\\.gz$", "", basename(fs))) %>%
    mutate(cross = cross,
           sample_name = sub("_.*$", "", afd_id)) %>%
    left_join(as_tibble(sh), by = "sample_name") %>%
    transmute(cross, sample_name, afd_id,
              p1_column = paste0("p1_", afd_id),
              p2_column = paste0("p2_", afd_id),
              parent1, parent2, timepoint, condition)
})
write.table(key_tbl, file.path(SD, "cross_af_samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n== written ==\n"); print(as.data.frame(written), row.names = FALSE)
cat("\n== sample key ==\n")
print(as.data.frame(key_tbl %>% select(cross, sample_name, timepoint,
                                       condition, parent1, parent2)),
      row.names = FALSE)
cat("\np1 is parent1 throughout; frequency of parent 1 is p1 / (p1 + p2).\n")
