#!/usr/bin/env Rscript
## Reduce GEMMA scans to one row per trait -- RUN THIS ON THE CLUSTER ---------
##
##   Rscript reduce_scans.R <scan_dir> [out.tsv] [key.tsv] [panel]
##
## <panel> names which panel these scans are, and matters: the key holds one row
## per (panel, trait), so joining on trait alone multiplies every row by the
## number of panels. If it is omitted the script tries to read it from the scan
## directory name and stops if that is ambiguous, rather than silently emitting
## a table three times too long.
##
## 672 scans at 464,045 markers each is about 30 GB raw and 7 GB gzipped, and
## scoring this simulation needs roughly fifteen numbers per trait. So the scans
## are reduced where they are produced and only the summary comes back -- a few
## hundred kilobytes.
##
## <scan_dir> is searched recursively for GEMMA association output: *.assoc.txt,
## or the csv/gz form the repo's own exports use. The trait is taken from a
## `trait` column when one exists and from the file name otherwise, so it works
## whether the pipeline writes one file per trait or one file with a trait
## column.
##
## WHAT IS KEPT, and why each one.
##   lambda_gc            median chi-squared over the whole scan, divided by
##                        qchisq(0.5, 1). The panels differ in structure, so a
##                        shared threshold is not safe; this is the correction.
##                        It needs every marker, which is exactly why it has to
##                        be computed here rather than after transfer.
##   max_logp + position  the genome-wide peak, for a false-positive rate
##   best_logp_*          the peak within 100 kb, 500 kb and 1 Mb of the causal
##                        site. Three windows rather than one so the detection
##                        rule can be chosen after the fact instead of being
##                        baked in now.
##   peak_dist_bp         how far the genome-wide peak sits from the truth --
##                        localisation, which no marker count measures
##   max_logp_offtarget   the peak more than 1 Mb from the causal site, or on
##                        another chromosome. A detection that is really a
##                        false positive somewhere else shows up here.
##
## Null traits have no causal site, so the window columns are NA and the
## genome-wide max is the quantity of interest.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE)
if (!length(a)) stop("usage: Rscript reduce_scans.R <scan_dir> [out.tsv] [key.tsv]")
SCAN <- a[1]
OUT  <- if (length(a) >= 2) a[2] else "scan_summary.tsv"
KEY  <- if (length(a) >= 3) a[3] else "simulation_key.tsv"
PANEL <- if (length(a) >= 4) a[4] else NA_character_
WIN  <- c(1e5, 5e5, 1e6)

f <- list.files(SCAN, pattern = "\\.(assoc\\.txt|csv|csv\\.gz|tsv|tsv\\.gz)$",
                full.names = TRUE, recursive = TRUE)
if (!length(f)) stop("no scan files under ", SCAN)
cat(sprintf("[%s] %d scan files under %s\n", format(Sys.time(), "%H:%M:%S"), length(f), SCAN))

key <- if (file.exists(KEY)) fread(KEY) else NULL
if (is.null(key)) {
  cat("NOTE: no key found, window columns will be NA\n")
} else if (uniqueN(key$panel) > 1) {
  ## resolve the panel, or stop -- see the header
  if (is.na(PANEL)) {
    slug <- function(x) gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(x)))
    hit <- unique(key$panel)[slug(unique(key$panel)) %in% slug(basename(normalizePath(SCAN)))]
    if (length(hit) == 1L) PANEL <- hit
  }
  if (is.na(PANEL) || !PANEL %in% key$panel)
    stop("the key holds ", uniqueN(key$panel), " panels and the panel could not be\n",
         "  resolved from the directory name. Pass it as the 4th argument, one of:\n    ",
         paste(unique(key$panel), collapse = "\n    "), call. = FALSE)
  cat(sprintf("panel: %s\n", PANEL))
  key <- key[panel == PANEL]
}

pcol <- function(nm) {
  for (p in c("p_wald", "p_lrt", "p_score", "p")) if (p %in% nm) return(p)
  stop("no p-value column in: ", paste(nm, collapse = ", "))
}

one <- function(d, tr) {
  k <- if (!is.null(key)) key[trait == tr][1] else NULL
  chi <- qchisq(d$P, df = 1, lower.tail = FALSE)
  lam <- median(chi, na.rm = TRUE) / qchisq(0.5, 1)
  d[, lp := -log10(pmax(P, .Machine$double.xmin))]
  i <- which.max(d$lp)
  win <- setNames(rep(NA_real_, length(WIN)), sprintf("best_logp_%gkb", WIN / 1e3))
  bp1 <- NA_real_; pdist <- NA_real_; off <- NA_real_
  if (!is.null(k) && nrow(k) && !is.na(k$chrom) && nzchar(as.character(k$chrom))) {
    dist <- ifelse(d$CHR == k$chrom, abs(d$BP - k$pos), Inf)
    for (j in seq_along(WIN)) {
      inw <- dist <= WIN[j]
      if (any(inw)) win[j] <- max(d$lp[inw])
    }
    inw <- dist <= 1e6
    if (any(inw))  bp1 <- d$BP[inw][which.max(d$lp[inw])]
    if (any(!inw)) off <- max(d$lp[!inw])
    if (d$CHR[i] == k$chrom) pdist <- abs(d$BP[i] - k$pos)
  }
  ## do.call, because data.table() does NOT splice a bare list argument into
  ## columns -- it recycles it into as many ROWS as the list is long, which
  ## silently tripled every row here before the join guard caught it
  do.call(data.table, c(
    list(trait = tr, n_markers = nrow(d), lambda_gc = lam,
         max_logp = d$lp[i], peak_chrom = as.character(d$CHR[i]),
         peak_pos = d$BP[i]),
    as.list(win),
    list(best_pos_1mb = bp1, peak_dist_bp = pdist, max_logp_offtarget = off)))
}

res <- rbindlist(lapply(seq_along(f), function(j) {
  d <- fread(f[j], showProgress = FALSE)
  nm <- names(d)
  P <- pcol(nm)
  ## column names differ between GEMMA's own output and the repo's exports
  ch <- intersect(c("chr", "CHR", "chrom"), nm)[1]
  bp <- intersect(c("ps", "BP", "pos", "physical.position"), nm)[1]
  setnames(d, c(ch, bp, P), c("CHR", "BP", "P"), skip_absent = TRUE)
  out <- if ("trait" %in% nm) {
    rbindlist(lapply(split(d, d$trait), function(dd) one(dd, dd$trait[1])))
  } else {
    one(d, sub("\\.(assoc\\.txt|csv|csv\\.gz|tsv|tsv\\.gz)$", "", basename(f[j])))
  }
  if (j %% 50 == 0) cat(sprintf("[%s]   %d / %d\n", format(Sys.time(), "%H:%M:%S"), j, length(f)))
  out
}), fill = TRUE)

stopifnot(!anyDuplicated(res$trait))     # one row per trait out of the reducer
if (!is.null(key)) {
  res <- merge(key, res, by = "trait", all.y = TRUE)
  stopifnot(!anyDuplicated(res$trait))   # and the join must not multiply them
}
fwrite(res, OUT, sep = "\t")
cat(sprintf("[%s] wrote %s -- %d rows, %.1f KB\n", format(Sys.time(), "%H:%M:%S"),
            OUT, nrow(res), file.size(OUT) / 1024))
cat("\nthis is the only file that needs to come back.\n")
