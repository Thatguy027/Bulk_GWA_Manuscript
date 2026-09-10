## Generate the seven-trait simulation from source, seeded and replicated -----
##
##   NREP=10 SEED=20260909 Rscript scripts/make_simulation_seeded.R
##     -> supplemental_data/deconvolution/simulation_seeded_frequencies.tsv.gz
##        supplemental_data/deconvolution/simulation_seeded_r2.tsv
##
## WHY THIS EXISTS
## The 2021 run's binomial draw was never seeded, so its estimates could be
## recomputed against but never regenerated, and its single draw reported each
## trait's recovery as one exact number. Now that the fitness input is recovered
## (simulation_fitness_traits.tsv) the whole simulation runs from source, so it
## can be seeded AND replicated: the spread across replicates is the sampling
## variability the single draw could only hide.
##
## WHAT IS DELIBERATELY THE SAME AS 2021
##   plain G, not the complement-stacked design. METHODS.txt records the stacked
##   design as worth about 0.14 of r-squared at 1x, so using it would change the
##   conclusions rather than only the numbers.
##   fitness = the published trait value shifted by the trait's own minimum,
##   with NA sent to 0 -- haploReg_original.R lines 327-328.
##   the same eight depths and the same 327 strains.
##
## WHAT IS DELIBERATELY DIFFERENT
##   RcppML::nnls on the normal equations rather than mdatools::mcrals.fcnnls on
##   the raw design. scripts/simulation_redraw.R measured these as the same
##   estimator to a median of 4e-06 in r-squared, and RcppML emits no negative
##   coefficients where the archive carries 139 of 18,312.
##   a seed, and NREP replicates instead of one draw.
##
## REQUIRES THE GENOTYPE PANEL (~1 GB resident, blockwise) and does not run from
## a clone; its two outputs are small and are what the figure reads.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr)})

NREP  <- as.integer(Sys.getenv("NREP", "10"))
SEED  <- as.integer(Sys.getenv("SEED", "20260909"))
PANEL <- Sys.getenv("CENDR_PLINK", "data/genotypes/CeNDR20210121_Plink")
DEC   <- "supplemental_data/deconvolution"
CHROMS <- c("I","II","III","IV","V","X")
DEPTHS <- c(500,100,50,30,10,5,3,1)
BLK    <- 100000L
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
if (!dir.exists(PANEL)) stop("genotype panel not found at ", PANEL, call. = FALSE)
set.seed(SEED)

## the 327 strains and their fitness, from the deposit
est <- read_tsv(file.path(DEC, "simulation_nnls_frequencies.tsv.gz"), show_col_types = FALSE)
strains <- sort(unique(est$strain)); stopifnot(length(strains) == 327)
tr <- read_tsv(file.path(DEC, "simulation_fitness_traits.tsv"), show_col_types = FALSE) %>%
  filter(strain %in% strains) %>% arrange(strain) %>%
  pivot_longer(-strain, names_to = "trait", values_to = "published") %>%
  group_by(trait) %>%
  mutate(fitness = replace_na(published - min(published, na.rm = TRUE), 0),
         input   = fitness / sum(fitness)) %>% ungroup()
TRAITS <- sort(unique(tr$trait))
W  <- sapply(TRAITS, function(t) tr$fitness[tr$trait == t])
TV <- sapply(TRAITS, function(t) tr$input[tr$trait == t])
rownames(W) <- rownames(TV) <- strains

## genotypes streamed as raw; GGp and the expected frequencies accumulated
msg("reading genotypes")
gb <- list(); GGp <- matrix(0, 327, 327); eC <- matrix(0, 0, length(TRAITS))
for (chr in CHROMS) {
  bm <- BEDMatrix::BEDMatrix(file.path(PANEL, paste0(chr, ".bed")), simple_names = FALSE)
  keep <- match(paste0(strains, "_", strains), rownames(bm)); stopifnot(!anyNA(keep))
  for (s in seq(1L, ncol(bm), by = BLK)) {
    e <- min(s + BLK - 1L, ncol(bm))
    g <- t(bm[keep, s:e, drop = FALSE]); g[is.na(g)] <- 0; g[g == 2] <- 1
    storage.mode(g) <- "double"
    GGp <- GGp + crossprod(g); eC <- rbind(eC, g %*% W)
    gb[[length(gb) + 1L]] <- structure(as.raw(g), dim = dim(g))
  }
  msg("  ", chr)
}
NM <- sum(vapply(gb, nrow, 1L))
eFreq <- sweep(eC, 2, colSums(W), "/")
stopifnot(all(eFreq >= 0), all(eFreq <= 1))
msg("markers ", NM)

out <- list()
for (rep in seq_len(NREP)) {
  for (ti in seq_along(TRAITS)) {
    Gy <- matrix(0, 327, length(DEPTHS)); off <- 0L
    for (b in gb) {
      nb <- nrow(b)
      cts <- vapply(DEPTHS, function(d) rbinom(nb, d, eFreq[(off + 1L):(off + nb), ti]),
                    numeric(nb))
      Gy <- Gy + crossprod(matrix(as.numeric(b), nb, 327), cts)
      off <- off + nb
    }
    K <- apply(Gy, 2, function(y) as.vector(RcppML::nnls(GGp, matrix(y), fast_nnls = TRUE)))
    K <- sweep(K, 2, colSums(K), "/")
    out[[length(out) + 1L]] <- tibble(
      replicate = rep, trait = TRAITS[ti],
      depth = rep(DEPTHS, each = 327), strain = rep(strains, length(DEPTHS)),
      frequency = as.vector(K), input = rep(TV[, ti], length(DEPTHS)))
  }
  msg("replicate ", rep, " of ", NREP)
}
freq <- bind_rows(out)
stopifnot(nrow(freq) == NREP * 7 * 8 * 327, all(freq$frequency >= 0))

write_tsv(freq, file.path(DEC, "simulation_seeded_frequencies.tsv.gz"))
r2 <- freq %>% group_by(replicate, trait, depth) %>%
  summarise(r2 = cor(frequency, input)^2, .groups = "drop")
write_tsv(r2, file.path(DEC, "simulation_seeded_r2.tsv"))

cat(sprintf("\nseed %d, %d replicates, %d negative coefficients\n",
            SEED, NREP, sum(freq$frequency < 0)))
cat("\n== mean r-squared across replicates, and the replicate range ==\n")
print(as.data.frame(r2 %>% group_by(trait, depth) %>%
  summarise(mean = round(mean(r2), 4), lo = round(min(r2), 4), hi = round(max(r2), 4),
            .groups = "drop") %>%
  filter(depth %in% c(1, 10, 50, 500)) %>% arrange(trait, depth)), row.names = FALSE)
msg("wrote simulation_seeded_{frequencies.tsv.gz,r2.tsv}")
