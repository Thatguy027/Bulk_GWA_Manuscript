## Recompute the simulation's reported r-squared from the deposit -------------
##
##   Rscript scripts/simulation_recompute_r2.R
##     supplemental_data/deconvolution/simulation_fitness_traits.tsv
##     supplemental_data/deconvolution/simulation_nnls_frequencies.tsv.gz
##     supplemental_data/deconvolution/simulation_reported_r2.tsv
##     -> plots/diagnostics/TABLE_simulation_recomputed_r2.tsv
##
##   Exits non-zero if any of the 56 recomputed values disagrees with the
##   reported one at the two decimals the figures carry.
##
## WHAT THIS SETTLES
## The 56 r-squared in the manuscript were reported, not recomputed: they were
## read back out of text embedded in the 2021 per-trait PDFs by
## scripts/extract_sim_reported_r2.py, because METHODS.txt recorded that the
## simulation's fitness input had never been archived. It had not, but it was
## recoverable -- the seven-trait arm of scripts/legacy/haploReg_original.R used
## published trait values as fitness rather than a draw, and those files are now
## staged by scripts/make_simulation_fitness_table.R. Given them, every reported
## value reproduces exactly at two decimals, so the paragraph rests on a
## recomputation rather than on a transcription.
##
## THE ORIGINAL'S TWO TRANSFORMS (haploReg_original.R lines 327-328)
##
##   phenop <- apply(phenos, 2, function(x) x - min(x, na.rm = TRUE))
##   phenop[is.na(phenop)] <- 0
##
## The shift makes fitness non-negative, which the forward model requires. The
## NA rule is the substantive one: a strain with no published value for a trait
## is assigned fitness exactly 0, i.e. it is absent from that simulated pool.
## That is the origin of the per-trait strain subsets in
## simulation_gwas_traits.tsv.gz, and it means the r-squared are computed over
## all 327 strains with the unmeasured ones pinned at zero -- stated here
## because it flatters the correlation. The shift contributes one more zero of
## its own: the lowest-scoring measured strain lands at exactly 0 and is
## therefore absent from the pool too, indistinguishable from an unmeasured
## one. Both counts are reported below rather than collapsed.
##
## r-squared is cor(estimate, fitness)^2 per trait per depth, on the
## coefficient scale and against the shifted trait values, which is what the
## original plots at its line 385. Correlation is scale-invariant, so
## normalising fitness would not change it.
##
## RUNS FROM A CLONE. Deposit only; no genotype panel and no re-simulation. The
## binomial draw behind the archived estimates was never seeded, so those
## estimates cannot be regenerated -- but they do not need to be, because they
## are what is deposited.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
})

SD  <- "supplemental_data/deconvolution"
OUT <- "plots/diagnostics"
DP  <- 2   # the decimals the figure text carries

ph  <- read_tsv(file.path(SD, "simulation_fitness_traits.tsv"), show_col_types = FALSE)
est <- read_tsv(file.path(SD, "simulation_nnls_frequencies.tsv.gz"), show_col_types = FALSE)
rep <- read_tsv(file.path(SD, "simulation_reported_r2.tsv"), show_col_types = FALSE)

## the archive fixes the strain set: the union of trait strains intersected with
## the genotype panel, which the original does by matching against colnames(g)
panel <- sort(unique(est$strain))
if (!all(panel %in% ph$strain))
  stop("archived output carries strains absent from the trait table: ",
       paste(setdiff(panel, ph$strain), collapse = ", "), call. = FALSE)

meas <- ph %>% filter(strain %in% panel) %>% arrange(strain) %>%
  pivot_longer(-strain, names_to = "trait", values_to = "published")

fit <- ph %>%
  filter(strain %in% panel) %>%
  arrange(strain) %>%
  mutate(across(-strain, ~ .x - min(.x, na.rm = TRUE)),
         across(-strain, ~ tidyr::replace_na(.x, 0))) %>%
  pivot_longer(-strain, names_to = "trait", values_to = "fitness") %>%
  left_join(meas, by = c("strain", "trait"))

r2 <- est %>%
  inner_join(fit, by = c("strain", "trait")) %>%
  group_by(trait, depth) %>%
  summarise(n.strains       = n(),
            n.measured      = sum(!is.na(published)),
            n.zero.fitness  = sum(fitness == 0),
            r2.recomputed   = cor(coefficient, fitness)^2,
            .groups = "drop") %>%
  inner_join(rep, by = c("trait", "depth")) %>%
  rename(r2.reported = r2) %>%
  mutate(agrees = round(r2.recomputed, DP) == r2.reported) %>%
  arrange(trait, desc(depth))

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write_tsv(r2, file.path(OUT, "TABLE_simulation_recomputed_r2.tsv"))

cat(sprintf("strains in the archived run : %d\n", length(panel)))
cat(sprintf("trait x depth pairs compared: %d\n", nrow(r2)))
cat(sprintf("agree at %d dp               : %d\n", DP, sum(r2$agrees)))
cat(sprintf("max |recomputed - reported| : %.4f\n\n",
            max(abs(round(r2$r2.recomputed, DP) - r2$r2.reported))))

cat("per trait: measured strains, and strains left at fitness 0\n")
cat("(zero = the unmeasured strains PLUS the single lowest-scoring one, which\n")
cat(" the shift by the trait minimum sends to exactly 0)\n\n")
pres <- r2 %>% group_by(trait) %>%
  summarise(strains = first(n.strains), measured = first(n.measured),
            unmeasured = first(n.strains) - first(n.measured),
            zero.fitness = first(n.zero.fitness), .groups = "drop") %>%
  arrange(desc(measured))
print(as.data.frame(pres), row.names = FALSE)

if (!all(r2$agrees)) {
  cat("\nDISAGREEMENTS:\n")
  print(as.data.frame(r2 %>% filter(!agrees)), row.names = FALSE)
  stop("recomputed r-squared do not reproduce the reported values", call. = FALSE)
}

cat(sprintf("\nall %d reported r-squared reproduce at %d dp\n", nrow(r2), DP))
