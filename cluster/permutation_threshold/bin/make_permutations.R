#!/usr/bin/env Rscript
## Permuted phenotype columns for GEMMA, in .fam order --------------------
##
## Writes perm_<trait>_b<offset>.txt files, each a headerless matrix of
## phenotype columns with one row per individual in the .fam file. GEMMA's -n
## selects a column, so a batch file is a batch of permutations.
##
## COLUMN 1 OF THE FIRST BATCH IS THE OBSERVED PHENOTYPE, not a permutation.
## Its genome-wide maximum therefore comes out of exactly the same code path as
## the permutations, which is what makes the empirical p-value trustworthy: a
## separately-computed observed maximum could differ for reasons that have
## nothing to do with the null.
##
## ORDER IS EVERYTHING. BIMBAM has no sample IDs -- the dosage columns are
## positional, in .fam order. A phenotype file in a different order silently
## maps the wrong values to the wrong strains and produces a plausible,
## meaningless scan. The join below is by strain ID and asserts completeness.
##
## Strains with a missing phenotype are written as NA; GEMMA drops them, which
## is the same behaviour as the real scan.

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(TRUE)
get <- function(f, d = NULL) {
  i <- match(f, args); if (is.na(i)) return(d); args[i + 1]
}
fam_f  <- get("--fam");   pheno_f <- get("--pheno"); trait <- get("--trait")
n_perm <- as.integer(get("--n_perm", 1000))
batch  <- as.integer(get("--batch", 25))
seed   <- as.integer(get("--seed", 1))
stopifnot(!is.null(fam_f), !is.null(pheno_f), !is.null(trait))

fam <- fread(fam_f, header = FALSE)
setnames(fam, 1:2, c("fid", "iid"))
ph  <- fread(pheno_f)
idc <- names(ph)[1]
if (!trait %in% names(ph))
  stop("trait '", trait, "' not in ", pheno_f, ". Available: ",
       paste(setdiff(names(ph), idc), collapse = ", "), call. = FALSE)

m <- match(fam$iid, ph[[idc]])
cat("  .fam individuals:", nrow(fam), "\n")
cat("  matched in the phenotype file:", sum(!is.na(m)), "\n")
if (all(is.na(m)))
  stop("no .fam individual matches the phenotype ID column '", idc,
       "' -- check that the trait file uses the same strain names as the VCF",
       call. = FALSE)
obs <- as.numeric(ph[[trait]][m])
n_ok <- sum(!is.na(obs))
cat("  with a non-missing", trait, "value:", n_ok, "\n")
if (n_ok < 20) stop("only ", n_ok, " phenotyped individuals; refusing", call. = FALSE)

fwrite(data.table(fid = fam$fid, iid = fam$iid, value = obs),
       paste0("strains_", trait, ".tsv"), sep = "\t")

## Permute only among the individuals that HAVE a value, leaving the missing
## ones missing. Shuffling NAs into phenotyped positions would change the
## sample size from one permutation to the next and make the maxima
## incomparable.
set.seed(seed)
have <- which(!is.na(obs))
cols <- vector("list", n_perm + 1L)
cols[[1]] <- obs                       # observed, as permutation 0
for (i in seq_len(n_perm)) {
  v <- obs
  v[have] <- obs[sample(have)]
  cols[[i + 1L]] <- v
}

## split into batch files; the file name carries the id of its first column so
## GEMMA_PERM can recover a global permutation id without another lookup
starts <- seq(1, length(cols), by = batch)
for (s in starts) {
  e   <- min(s + batch - 1L, length(cols))
  dt  <- as.data.table(cols[s:e])
  ## permutation ids are 0-based: column 1 overall is the observed data
  fwrite(dt, sprintf("perm_%s_b%d.txt", trait, s - 1L),
         sep = "\t", col.names = FALSE, na = "NA")
}
cat("  wrote", length(starts), "batch files of up to", batch, "columns",
    "(", n_perm, "permutations plus the observed )\n")
