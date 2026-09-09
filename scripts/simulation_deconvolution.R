## The NNLS deconvolution simulation -------------------------------------------
##
##   Rscript scripts/simulation_deconvolution.R [experiment]
##
##   experiment = depth      (default) recovery against known input, by depth
##                panel      recovery as a function of reference panel size
##                dropout    recovery when a fraction of the panel is absent
##
##   -> plots/diagnostics/TABLE_simulation_<experiment>.tsv
##
## PROVENANCE. This is a port of the simulation described in METHODS.txt, which
## until now was documented but not archived. The original is kept verbatim at
## scripts/legacy/haploReg_original.R -- a console log covering two organisms
## and a year of exploratory work, which does not parse (an unclosed cor.test(
## at its line 178) and cannot be sourced. Lines 1-114, 142-186 and 201-242 of
## that file are the C. elegans simulation and are what this script rebuilds;
## its yeast sections, its scratch comparison of four NNLS solvers, and its
## vendored copy of the NMF package's internal .fcnnls are deliberately not
## carried over. The solver here is RcppML::nnls, which is what the surviving
## working code uses.
##
## WHAT THE ORIGINAL PINS THAT METHODS.txt COULD NOT STATE. Two [TO FILL] blocks
## are answerable from it:
##
##   the inverse chi-squared draw is extraDistr::rinvchisq(n, df = 12, scale = 1)
##
##   fitness maps to expected pooled frequency as eCount = rowSums(G %*% w) and
##   eFreq = eCount / sum(w), i.e. an allele's expected frequency is the
##   fitness-weighted mean of the strains carrying it
##
##   the bootstrap resamples MARKERS with replacement -- not reads and not
##   strains -- at set.seed(100), nboot = 100
##
## The marker reading is not a guess. supplemental_data/deconvolution/
## baugh_bootstrap_array.rda carries objects named `ab`, `blist` and `bootse`,
## which are the original's own variable names at its lines 94, 98 and 101, and
## `ab` is 102 x 23 x 100. Resampling rows of a markers x strains matrix is the
## only one of the three candidates that yields a strains x samples estimate per
## replicate; resampling strains would change the first dimension and reads
## would not index nrow(G) at all.
##
## ONE DISCREPANCY, CARRIED FORWARD RATHER THAN SILENTLY RESOLVED. The original
## fits its point estimate on the COMPLEMENT-STACKED matrix rbind(G, 1 - G) (its
## line 17), so that reference-allele counts enter as their own equations, but
## its bootstrap (lines 78-79) fits on the plain G. The published intervals are
## therefore intervals on a different estimator from the published point
## estimates. `stacked` below selects which, it is reported in the output table,
## and the two are not silently reconciled here because which one the archived
## Baugh run used is a question about that run, not about this code.
##
## REQUIRES THE GENOTYPE PANEL, so this does not run from a clone. See
## SYNC_MANIFEST.md Tier 1; set CENDR_PLINK to override the default location.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

PANEL <- Sys.getenv("CENDR_PLINK", "data/genotypes/CeNDR20210121_Plink")
OUT   <- "plots/diagnostics"
CHROMS <- c("I", "II", "III", "IV", "V", "X")
DEPTHS <- c(500, 100, 50, 30, 10, 5, 3, 1)
READ_LENGTH <- 300
GENOME_SIZE <- 120e6

## the original's constants, kept as literals so they are greppable
INVCHISQ_DF    <- 12
INVCHISQ_SCALE <- 1
NBOOT          <- 100
SEED           <- 100

need <- function(pkg) if (!requireNamespace(pkg, quietly = TRUE))
  stop("package '", pkg, "' is required; install it or run a different experiment",
       call. = FALSE)

## ---------------------------------------------------------------------------
## genotypes: markers x strains, biallelic, missing treated as reference
## ---------------------------------------------------------------------------
## The original recodes 2 -> 1 (its line 10). These are selfing wild isolates
## called as homozygous, so the alternate-allele dosage is 0 or 2 and the recode
## makes the design matrix an indicator of "carries the alternate allele". That
## is what makes the column sums interpretable as pool frequencies below.
load_panel <- function(panel = PANEL, chroms = CHROMS) {
  need("BEDMatrix")
  if (!dir.exists(panel))
    stop("genotype panel not found at '", panel, "'.\n",
         "  It is not in the repository -- see SYNC_MANIFEST.md, Tier 1.\n",
         "  Set CENDR_PLINK to point at your copy.", call. = FALSE)
  mats <- lapply(chroms, function(chr) {
    bed <- file.path(panel, paste0(chr, ".bed"))
    if (!file.exists(bed)) stop("missing ", bed, call. = FALSE)
    message("  reading ", chr)
    g <- t(as.matrix(BEDMatrix::BEDMatrix(bed)))
    g[is.na(g)] <- 0
    g[g == 2] <- 1
    g
  })
  do.call(rbind, mats)
}

## rbind(G, 1 - G): the reference allele gets its own equations, so both alleles
## constrain the fit rather than only the alternate.
stack_complement <- function(g) rbind(g, abs(g - 1))

## ---------------------------------------------------------------------------
## the forward model, then NNLS back
## ---------------------------------------------------------------------------
expected_freq <- function(g, w) as.vector(rowSums(g %*% w)) / sum(w)

simulate_counts <- function(efreq, depths = DEPTHS) {
  out <- vapply(depths, function(d) rbinom(length(efreq), size = d, prob = efreq),
                numeric(length(efreq)))
  colnames(out) <- as.character(depths)
  out
}

## Solved on the normal equations, which is what makes this tractable at
## ~1e6 markers: crossprod collapses the marker dimension once, and the solver
## then works on a strains x strains system.
deconvolve <- function(g, counts) {
  need("RcppML")
  gtg <- crossprod(g)
  gty <- crossprod(g, counts)
  est <- apply(gty, 2, function(y)
    as.vector(RcppML::nnls(gtg, matrix(y), fast_nnls = TRUE)))
  apply(est, 2, function(x) x / sum(x))
}

## Resamples MARKERS with replacement, nboot times, refitting from scratch each
## time. Returns the per-strain standard deviation across replicates.
bootstrap_se <- function(g, counts, nboot = NBOOT, seed = SEED) {
  need("RcppML"); need("abind")
  set.seed(seed)
  reps <- lapply(seq_len(nboot), function(i) {
    if (i %% 10 == 0) message("  bootstrap ", i, "/", nboot)
    idx <- sample.int(nrow(g), replace = TRUE)
    deconvolve(g[idx, , drop = FALSE], counts[idx, , drop = FALSE])
  })
  ab <- abind::abind(reps, along = 3)
  list(array = ab, se = apply(ab, c(1, 2), sd))
}

## ---------------------------------------------------------------------------
## experiments
## ---------------------------------------------------------------------------

## Recovery against known input across the depth series, with bootstrap SEs.
exp_depth <- function(g, stacked = TRUE, boot = TRUE) {
  need("extraDistr")
  set.seed(SEED)
  w <- extraDistr::rinvchisq(ncol(g), INVCHISQ_DF, INVCHISQ_SCALE)
  design <- if (stacked) stack_complement(g) else g
  efreq  <- expected_freq(design, w)
  counts <- simulate_counts(efreq)
  est    <- deconvolve(design, counts)
  truth  <- w / sum(w)

  se <- if (boot) bootstrap_se(design, counts)$se else NULL
  data.table(
    depth      = DEPTHS,
    yield      = (DEPTHS * GENOME_SIZE) / READ_LENGTH,
    r2         = apply(est, 2, function(x) cor(x, truth)^2),
    spearman   = apply(est, 2, function(x) cor(x, truth, method = "spearman")),
    mean_abs_err = apply(est, 2, function(x) mean(abs(x - truth))),
    median_se  = if (is.null(se)) NA_real_ else apply(se, 2, median),
    n_strain   = ncol(g),
    stacked    = stacked)
}

## Does a bigger reference panel help or hurt? The original samples a subset of
## strains, rebuilds the reference from that subset alone, and deconvolves a
## pool drawn from the same subset.
exp_panel <- function(g, sizes = c(50, 100, 200, 300, 400, 540), depth = 50,
                      n_rep = 20, stacked = TRUE) {
  need("extraDistr")
  set.seed(SEED)
  sizes <- sizes[sizes <= ncol(g)]
  rbindlist(lapply(sizes, function(k) {
    message("  panel size ", k)
    r <- vapply(seq_len(n_rep), function(i) {
      pick   <- sort(sample.int(ncol(g), k))
      gs     <- g[, pick, drop = FALSE]
      design <- if (stacked) stack_complement(gs) else gs
      w      <- extraDistr::rinvchisq(k, INVCHISQ_DF, INVCHISQ_SCALE)
      efreq  <- expected_freq(design, w)
      est    <- deconvolve(design, simulate_counts(efreq, depth))
      cor(as.vector(est), w / sum(w))^2
    }, numeric(1))
    data.table(panel_size = k, depth = depth, n_rep = n_rep,
               r2_mean = mean(r), r2_sd = sd(r),
               r2_min = min(r), r2_max = max(r), stacked = stacked)
  }))
}

## What happens to strains that are in the reference but absent from the pool?
## This is the leakage the dilution experiment measures, in simulation: a
## fraction of the panel is given zero fitness and the pool is built without it.
##
## The original gives every present strain fitness exactly 1 (its line 213,
## `(runif(ncol(g)) > .25) + 0`), which makes the pool uniform and leaves no
## variance among present strains to recover -- so recovery among them cannot be
## scored at all. Present strains are drawn from the same inverse chi-squared as
## everywhere else here, which keeps the absent set the only thing that changes.
exp_dropout <- function(g, absent_frac = c(0, 0.1, 0.25, 0.5), depth = 50,
                        stacked = TRUE) {
  need("extraDistr")
  set.seed(SEED)
  design <- if (stacked) stack_complement(g) else g
  rbindlist(lapply(absent_frac, function(f) {
    message("  absent fraction ", f)
    present <- runif(ncol(g)) >= f
    if (sum(present) < 2) return(NULL)
    w <- extraDistr::rinvchisq(ncol(g), INVCHISQ_DF, INVCHISQ_SCALE)
    w[!present] <- 0
    efreq  <- expected_freq(design, w)
    est    <- deconvolve(design, simulate_counts(efreq, depth))[, 1]
    data.table(absent_frac = f,
               n_absent    = sum(!present),
               n_present   = sum(present),
               depth       = depth,
               ## the share of inferred frequency that landed on strains which
               ## were not in the pool at all
               leaked      = sum(est[!present]),
               max_leak    = if (any(!present)) max(est[!present]) else 0,
               r2_present  = cor(est[present], w[present] / sum(w))^2,
               stacked     = stacked)
  }))
}

## ---------------------------------------------------------------------------

main <- function() {
  which <- commandArgs(trailingOnly = TRUE)
  which <- if (length(which)) which[1] else "depth"
  if (!which %in% c("depth", "panel", "dropout"))
    stop("experiment must be one of: depth, panel, dropout", call. = FALSE)

  dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
  message("loading genotype panel from ", PANEL)
  g <- load_panel()
  message("  ", nrow(g), " markers x ", ncol(g), " strains")

  res <- switch(which,
                depth   = exp_depth(g),
                panel   = exp_panel(g),
                dropout = exp_dropout(g))

  f <- file.path(OUT, paste0("TABLE_simulation_", which, ".tsv"))
  fwrite(res, f, sep = "\t")
  print(as.data.frame(res), row.names = FALSE)
  cat("\nwrote ", f, "\n", sep = "")
}

if (sys.nframe() == 0L) main()
