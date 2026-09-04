## Full vs thinned bundle: what the supplemental deposit costs ---------------
##
##   Rscript scripts/compare_full_vs_thinned_bundle.R
##     -> plots/legacy/compare_bundle/Figure2_FULL.{pdf,png}
##        plots/legacy/compare_bundle/Figure2_THINNED.{pdf,png}
##        plots/legacy/compare_bundle/Figure2_stacked.png
##        plots/legacy/compare_bundle/interval_comparison.tsv
##
## Figure 2 is built twice, once from each bundle, and the two are compared on
## the numbers that matter and on pixels. This is the audit behind the decision
## to ship the thinned bundle as supplemental data.
##
## REQUIRES THE FULL BUNDLE, which is not in the git repository and not in the
## deposit: data/pooled_cross_intersection/bundle.rds, 258.7 MB, from the Dryad
## archive or rebuilt with scripts/pooled_cross_intersection_prep.R. The script
## exits with a clear message rather than an error if it is absent.
##
## What it reports:
##   peak position and peak LOD per interval, both bundles
##   interval width per interval, and the difference
##   the fraction of identical pixels in the rendered PNG
##
## Expected result, as of the 5 kb / 10 kb thinning: peaks and peak LODs
## identical, widths contracted by 0-9 kb, 95.7% of pixels identical.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr)
})

FULL <- "data/pooled_cross_intersection/bundle.rds"
THIN <- "supplemental_data/mapping/pooled_cross_bundle_thinned.rds"
OUT  <- "plots/legacy/compare_bundle"
FIG2 <- "scripts/Figure2.R"

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

if (!file.exists(FULL)) {
  msg("full bundle not present at ", FULL)
  msg("  It is not in the repository or the supplemental deposit -- take it")
  msg("  from Dryad, or rebuild with scripts/pooled_cross_intersection_prep.R.")
  msg("  Nothing to compare; exiting cleanly.")
  quit(status = 0)
}
stopifnot(file.exists(THIN), file.exists(FIG2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## Build Figure 2 from a given bundle by rewriting only the bundle path and the
## output name, so the comparison uses the SHIPPED figure code rather than a
## copy of it that could drift.
build_from <- function(bundle, stem) {
  src <- read_file(FIG2)
  src <- str_replace_all(src, fixed(THIN), bundle)
  src <- str_replace_all(src, fixed('OUT <- "plots"'), sprintf('OUT <- "%s"', OUT))
  src <- str_replace_all(src, fixed('build(LOD_MIN, "Figure2"'),
                         sprintf('build(LOD_MIN, "%s"', stem))
  tmp <- tempfile(fileext = ".R"); write_file(src, tmp)
  log <- tempfile(fileext = ".log")
  st <- system2("Rscript", tmp, stdout = log, stderr = log)
  if (st != 0) { cat(read_file(log)); stop("build failed for ", stem) }
  read_file(log)
}

msg("building Figure 2 from the FULL bundle (", round(file.size(FULL)/1e6), " MB)")
log_full <- build_from(FULL, "Figure2_FULL")
msg("building Figure 2 from the THINNED bundle (",
    round(file.size(THIN)/1e6, 1), " MB)")
log_thin <- build_from(THIN, "Figure2_THINNED")

## --- the interval table each run prints ----------------------------------
parse_intervals <- function(log) {
  ln <- str_split(log, "\n")[[1]]
  i <- which(str_detect(ln, "== intervals on Figure 2"))
  if (!length(i)) stop("interval table not found in the log")
  rows <- ln[(i + 2):length(ln)]
  rows <- rows[str_detect(rows, "^\\s*(DN|UP)\\s")]
  read_table(paste(rows, collapse = "\n"),
             col_names = c("dir", "cond1", "vs", "cond2", "cross", "chrom",
                           "peak_Mb", "interval", "width_kb", "peak_LOD"),
             show_col_types = FALSE)
}
a <- parse_intervals(log_full); b <- parse_intervals(log_thin)
stopifnot(nrow(a) == nrow(b))

cmp <- a %>%
  transmute(track = paste(dir, cond1, vs, cond2), cross, chrom,
            peak_Mb_full = peak_Mb, peak_Mb_thin = b$peak_Mb,
            peak_LOD_full = peak_LOD, peak_LOD_thin = b$peak_LOD,
            interval_full = interval, interval_thin = b$interval,
            width_full = width_kb, width_thin = b$width_kb,
            width_delta = b$width_kb - width_kb)
write_tsv(cmp, file.path(OUT, "interval_comparison.tsv"))

cat("\n== interval comparison ==\n")
print(as.data.frame(cmp %>% select(cross, chrom, peak_Mb_full, peak_Mb_thin,
                                   peak_LOD_full, peak_LOD_thin,
                                   width_full, width_thin, width_delta)),
      row.names = FALSE)
cat("\npeak positions identical: ", all(cmp$peak_Mb_full == cmp$peak_Mb_thin),
    "\npeak LODs identical     : ", all(cmp$peak_LOD_full == cmp$peak_LOD_thin),
    "\nwidth change (kb)       : ", min(cmp$width_delta), " to ",
    max(cmp$width_delta), "\n", sep = "")

## --- pixels ---------------------------------------------------------------
if (requireNamespace("png", quietly = TRUE)) {
  pa <- png::readPNG(file.path(OUT, "Figure2_FULL.png"))
  pb <- png::readPNG(file.path(OUT, "Figure2_THINNED.png"))
  if (identical(dim(pa), dim(pb))) {
    d <- apply(abs(pa - pb), c(1, 2), max) * 255
    cat("identical pixels        : ", sprintf("%.2f%%", 100 * mean(d == 0)),
        "\npixels differing > 32   : ", sum(d > 32), " of ", length(d),
        "\nmax channel difference  : ", round(max(d)), "/255\n", sep = "")
    ## stacked image, full over thinned, for eyeballing
    st <- array(1, dim = c(dim(pa)[1] * 2 + 20, dim(pa)[2], 3))
    st[seq_len(dim(pa)[1]), , ] <- pa[, , 1:3]
    st[dim(pa)[1] + 20 + seq_len(dim(pb)[1]), , ] <- pb[, , 1:3]
    png::writePNG(st, file.path(OUT, "Figure2_stacked.png"))
    msg("wrote Figure2_stacked.png (full over thinned)")
  } else msg("PNG dimensions differ; skipping the pixel comparison")
}
msg("outputs in ", OUT)
