#!/usr/bin/env Rscript
## The strain panel the scan actually used ----------------------------------
##
## Writes panel.txt (FID IID, plink --keep format) listing every strain with a
## non-missing value for at least one requested trait.
##
## WHY THIS EXISTS. Minor-allele frequency has to be computed among the
## PHENOTYPED strains. A marker at 5% across all 540 CeNDR isotypes can be below
## 5% among the 231 that carry a pos-1 value, and the reverse. Filtering on the
## full collection produced 519,341 markers where the scan had 464,045, and the
## observed genome-wide maximum came out at 8.69 instead of 8.84 -- a different
## marker set, so a threshold that did not apply to the scan it was for.
##
## plink --keep matches on FID and IID. The CeNDR fam uses the strain name for
## both, so both columns are the strain name here.

suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(TRUE)
get <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) d else args[i + 1] }
pf <- get("--pheno"); tr_arg <- get("--traits")
stopifnot(!is.null(pf), !is.null(tr_arg))

ph <- fread(pf)
idc <- names(ph)[1]
traits <- trimws(strsplit(tr_arg, ",")[[1]])
missing <- setdiff(traits, names(ph))
if (length(missing))
  stop("trait(s) not in ", pf, ": ", paste(missing, collapse = ", "),
       "\n  available: ", paste(setdiff(names(ph), idc), collapse = ", "),
       call. = FALSE)

## ONE TRAIT PER RUN, ENFORCED. The panel defines the strain set MAF is
## computed on, and different traits have different panels: in this file the
## three pos-1 traits have 231 strains each while
## negctrl_growth_HT115_delta_t0 has all 366. Taking the union across traits
## would compute MAF on a superset for at least one of them -- exactly the bug
## this script exists to prevent, reintroduced by the back door. Nothing here
## can express a per-trait marker set, so more than one trait is refused rather
## than silently averaged over.
n_by_trait <- vapply(traits, function(t) sum(!is.na(ph[[t]])), integer(1))
if (length(traits) > 1L && length(unique(n_by_trait)) > 1L)
  stop("traits requested together have different panels:\n",
       paste(sprintf("    %-34s %d strains", traits, n_by_trait), collapse = "\n"),
       "\n  MAF must be computed per panel, so run them separately:\n",
       paste(sprintf("    --traits '%s'", traits), collapse = "\n"), call. = FALSE)

keep <- ph[[idc]][Reduce(`|`, lapply(traits, function(t) !is.na(ph[[t]])))]
keep <- unique(keep[!is.na(keep) & nzchar(keep)])
if (!length(keep)) stop("no strain has a value for any requested trait", call. = FALSE)
stopifnot(length(keep) == max(n_by_trait))

fwrite(data.table(fid = keep, iid = keep), "panel.txt", sep = " ", col.names = FALSE)
writeLines(c(
  sprintf("traits: %s", paste(traits, collapse = ", ")),
  sprintf("strains in the phenotype file: %d", nrow(ph)),
  sprintf("strains with a value for at least one trait: %d", length(keep)),
  "",
  "per trait:",
  vapply(traits, function(t) sprintf("  %-34s %d", t, sum(!is.na(ph[[t]]))),
         character(1))), "panel_summary.txt")
cat("  panel:", length(keep), "strains from", nrow(ph), "in the phenotype file\n")
for (t in traits) cat(sprintf("    %-34s %d with a value\n", t, sum(!is.na(ph[[t]]))))
