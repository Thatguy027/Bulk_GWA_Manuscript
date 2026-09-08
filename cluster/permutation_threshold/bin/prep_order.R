#!/usr/bin/env Rscript
## Check the strain order against the phenotype file -------------------------
##
## The order comes from the SCAN's own oxford .sample file, not from anything
## recomputed here, because BIMBAM has no sample ids -- its dosage columns are
## positional. This script does not choose the order; it verifies that the order
## the scan used is usable with the phenotype file it is being given, and says
## how many strains carry each trait.
##
## Refuses two things:
##   * a strain in the scan's order that the phenotype file does not contain,
##     which would mean the two came from different experiments
##   * traits whose non-missing strain sets differ from each other, when more
##     than one is requested. GEMMA's -gk uses phenotype column 1 to choose the
##     individuals for the kinship, so a second trait with a different set would
##     be mapped against a kinship built for the first one's.

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(TRUE)
get <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) return(d); args[i + 1] }
pf     <- get("--pheno"); of <- get("--order")
traits <- trimws(strsplit(get("--traits", ""), ",")[[1]])
expect <- as.integer(get("--expect_individuals", "0"))
stopifnot(!is.null(pf), !is.null(of), length(traits) > 0, nzchar(traits[1]))

ord <- readLines(of)
ph  <- fread(pf)
idc <- names(ph)[1]

missing_trait <- setdiff(traits, names(ph))
if (length(missing_trait))
  stop("trait(s) not in ", pf, ": ", paste(missing_trait, collapse = ", "),
       "\n  available: ", paste(setdiff(names(ph), idc), collapse = ", "),
       call. = FALSE)

if (expect > 0 && length(ord) != expect)
  stop("the strain order has ", length(ord), " strains, expected ", expect,
       ".\n  This file defines which BIMBAM column is which strain and must be",
       "\n  the one the scan used.", call. = FALSE)

absent <- setdiff(ord, ph[[idc]])
if (length(absent))
  stop(length(absent), " strain(s) in the scan's order are not in the ",
       "phenotype file:\n  ", paste(head(absent, 10), collapse = ", "),
       if (length(absent) > 10) ", ..." else "",
       "\n  The order file and the phenotype file describe different panels.",
       call. = FALSE)

sets <- lapply(traits, function(t) ord[!is.na(ph[[t]][match(ord, ph[[idc]])])])
names(sets) <- traits
if (length(traits) > 1) {
  ref <- sets[[1]]
  bad <- traits[vapply(sets, function(s) !setequal(s, ref), logical(1))]
  if (length(bad))
    stop("these traits have a different set of phenotyped strains from '",
         traits[1], "': ", paste(bad, collapse = ", "),
         "\n  GEMMA's -gk takes its individuals from phenotype column 1, so one",
         "\n  kinship cannot serve both. Run them as separate jobs.",
         call. = FALSE)
}

lines <- c(paste("traits:", paste(traits, collapse = ", ")),
           paste("strains in the scan's order:", length(ord)),
           paste("strains in the phenotype file:", nrow(ph)),
           "",
           "per trait, strains with a value:",
           vapply(traits, function(t) sprintf("  %-34s %d", t, length(sets[[t]])), ""))
writeLines(lines, "panel_summary.txt")
writeLines(lines)
