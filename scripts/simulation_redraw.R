## Redraw the seven-trait simulation and compare with the archived estimates --
##
##   NREP=3 RAW_REPS=1 SEED=1 Rscript scripts/simulation_redraw.R
##     -> $SP/redraw3_r2.tsv   (SP defaults to the working directory)
##
## The archived run's binomial draw was never seeded, so the question is not
## whether the estimates come back identical -- they cannot -- but how much the
## redraw costs, and whether the estimator the port uses is the same estimator
## the archive used. Three arms on IDENTICAL drawn counts separate the two:
##
##   fcnnls.raw   mcrals.fcnnls(counts, G)  -- the archived run's exact call
##   fcnnls.ne    mcrals.fcnnls(Gy, GGp)    -- same solver, normal equations
##   rcppml.ne    RcppML::nnls(GGp, Gy)     -- what simulation_deconvolution.R does
##
## fcnnls.raw vs the archive  = redraw noise (same solver, same formulation)
## fcnnls.ne  vs fcnnls.raw   = the cost of squaring the condition number
## rcppml.ne  vs fcnnls.ne    = solver difference at equal conditioning
##
## WHAT IT FOUND (3 replicates, 2,917,997 markers x 327 strains, 2026-09-09).
## The estimator is immaterial: fcnnls.raw and rcppml.ne agree to a median 4e-06
## and a max 1.1e-03 in r-squared. The formulation is not: fcnnls.ne diverges by
## up to 6.2e-02 and warns that it hit its iteration cap. The unseeded draw
## costs about 0.01 of r-squared at 10x and less deeper; archived values deviate
## from the redraw mean by a median -0.0002 (range -0.0245 to +0.0808) and sit
## inside the three-replicate range in 32 of 56 cells. Negative coefficients:
## fcnnls.raw emits 2 of 18,312 (0.011%, min -0.0021) against the archive's 139
## (0.76%), so the solver is implicated in those but does not account for them.
## Results are archived at plots/diagnostics/TABLE_simulation_redraw_r2.tsv.
##
## REQUIRES THE GENOTYPE PANEL and ~10 GB of memory: it holds G as a 7.1 GB
## double matrix so that the raw-design fcnnls can be called at all. Does not
## run from a clone. Roughly an hour with RAW_REPS=1.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(mdatools)
})
NREP   <- as.integer(Sys.getenv("NREP", "1"))
CHROMS <- c("I","II","III","IV","V","X")[seq_len(as.integer(Sys.getenv("MAXCHR","6")))]
SPD    <- Sys.getenv("SP", ".")
PANEL  <- "data/genotypes/CeNDR20210121_Plink"
DEC    <- "supplemental_data/deconvolution"
DEPTHS <- c(500,100,50,30,10,5,3,1)
BLK    <- 100000L
set.seed(as.integer(Sys.getenv("SEED","1")))
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep="")

est <- read_tsv(file.path(DEC,"simulation_nnls_frequencies.tsv.gz"), show_col_types=FALSE)
strains <- sort(unique(est$strain)); stopifnot(length(strains)==327)
tr <- read_tsv(file.path(DEC,"simulation_fitness_traits.tsv"), show_col_types=FALSE) %>%
  filter(strain %in% strains) %>% arrange(strain) %>%
  pivot_longer(-strain, names_to="trait", values_to="published") %>%
  group_by(trait) %>%
  mutate(fitness = replace_na(published - min(published, na.rm=TRUE), 0),
         input   = fitness/sum(fitness)) %>% ungroup()
TRAITS <- sort(unique(tr$trait))
W  <- sapply(TRAITS, function(t) tr$fitness[tr$trait==t])   # strain order = sorted
TV <- sapply(TRAITS, function(t) tr$input[tr$trait==t])
rownames(W) <- rownames(TV) <- strains

## ---- count markers, then fill one double matrix ----
nmk <- vapply(CHROMS, function(c) nrow(data.table::fread(
  file.path(PANEL, paste0(c,".bim")), select=1L, showProgress=FALSE)), 1L)
NMARK <- sum(nmk)
msg("allocating G: ", NMARK, " x 327 = ", round(NMARK*327*8/1024^3, 2), " GB")
G <- matrix(0, NMARK, 327)
off <- 0L
for (ci in seq_along(CHROMS)) {
  bm <- BEDMatrix::BEDMatrix(file.path(PANEL, paste0(CHROMS[ci],".bed")), simple_names=FALSE)
  keep <- match(paste0(strains,"_",strains), rownames(bm)); stopifnot(!anyNA(keep))
  nm <- ncol(bm)
  for (s in seq(1L, nm, by=BLK)) {
    e <- min(s+BLK-1L, nm)
    g <- t(bm[keep, s:e, drop=FALSE])
    g[is.na(g)] <- 0; g[g==2] <- 1
    G[(off+s):(off+e), ] <- g
  }
  off <- off + nm
  msg("  ", CHROMS[ci], ": ", nm, " markers")
}
msg("computing GGp and expected frequencies")
GGp   <- crossprod(G)
eFreq <- (G %*% W) / rep(colSums(W), each = NMARK)
stopifnot(all(eFreq >= 0), all(eFreq <= 1))

norm1 <- function(x) x / sum(x)
out <- list()
for (rep in seq_len(NREP)) for (ti in seq_along(TRAITS)) {
  cts <- vapply(DEPTHS, function(d) rbinom(NMARK, d, eFreq[, ti]), numeric(NMARK))
  Gy  <- crossprod(G, cts)

  ## the raw-matrix fcnnls recomputes crossprod(G) on every call, so it runs
  ## only on the replicates named by RAW_REPS -- it agreed with rcppml.ne to
  ## 1.5e-3 on a chromosome I pilot, and one full-scale replicate is enough to
  ## confirm that and to settle whether it emits the archive's negatives
  do_raw <- rep <= as.integer(Sys.getenv("RAW_REPS","1"))
  Kraw <- if (do_raw) mcrals.fcnnls(cts, G) else NULL
  Kne  <- mcrals.fcnnls(Gy,  GGp)                                 # 8 x 327
  Krcp <- t(apply(Gy, 2, function(y) as.vector(RcppML::nnls(GGp, matrix(y), fast_nnls=TRUE))))

  arms <- if (do_raw) c("fcnnls.raw","fcnnls.ne","rcppml.ne") else c("fcnnls.ne","rcppml.ne")
  for (arm in arms) {
    K <- switch(arm, fcnnls.raw=Kraw, fcnnls.ne=Kne, rcppml.ne=Krcp)
    out[[length(out)+1L]] <- tibble(
      rep = rep, trait = TRAITS[ti], arm = arm, depth = DEPTHS,
      r2  = vapply(seq_along(DEPTHS), function(j) cor(norm1(K[j,]), TV[,ti])^2, 0),
      n.neg = vapply(seq_along(DEPTHS), function(j) sum(K[j,] < 0), 0L),
      min.coef = vapply(seq_along(DEPTHS), function(j) min(K[j,]), 0))
  }
  msg("  rep ", rep, " ", TRAITS[ti], " done")
}
write_tsv(bind_rows(out), file.path(SPD, "redraw3_r2.tsv"))
msg("wrote redraw3_r2.tsv")
