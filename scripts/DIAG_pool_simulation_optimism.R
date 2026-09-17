## Why the deconvolution simulation looks better than the real data -----------
##
##   Rscript scripts/DIAG_pool_simulation_optimism.R
##     -> plots/diagnostics/DIAG_pool_simulation_optimism.{pdf,png}
##     -> plots/diagnostics/simulation_optimism_r2.tsv
##     -> plots/diagnostics/simulation_optimism_perstrain.tsv
##
## The seven-trait simulation reports r-squared rising above 0.95 from 10x. The
## panel optimiser was built on the opposite premise -- that strains are hard to
## tell apart and identifiability is the binding constraint. Both cannot be the
## whole story, so this asks what the simulation's number is actually measuring.
##
## PART A. THE ZERO BLOCK. The original transform sends a strain with no
## published value for a trait to fitness exactly 0, i.e. absent from the pool,
## and non-negative least squares returns exactly 0 for an absent strain because
## zero is a hard boundary of the feasible set. Those points sit at (0, 0) and
## are free. This recomputes every r-squared with them in and with them out.
##
## The internal control is mtDNA_ratio. It is 1.5% zeros where PC1 is 74.6%, so
## if the zero block is what inflates the number, mtDNA_ratio should barely move
## and PC1 should collapse. That is a prediction the data can refuse.
##
## PART B. WHAT PRIVATENESS BUYS IN A CLEAN WORLD. The simulation generates its
## reads from the same genotype matrix it then solves with, so the only error is
## binomial sampling -- no reference mismatch, no genotype-call error, no strain
## in the pool missing from the reference. If privateness predicts per-strain
## error STRONGLY here and weakly against real MIP-seq (rho = -0.42, adj R2 =
## 0.165 in DIAG_baugh_rmsd_predictors.R), then identifiability is doing its job
## and the real-data residual is reference quality -- which no panel choice
## fixes, and which means K* buys less than the optimiser's framing implies.
##
## Private markers are counted within the 327-strain REFERENCE, not within a
## trait's pool: the reference is what the solve inverts, so that is what sets
## identifiability. Genotypes are handled exactly as make_simulation_seeded.R
## does -- NA to 0, 2 to 1, the same six chromosomes.
##
## REQUIRES THE 2021 GENOTYPE PANEL. Exploratory, on the pool_optimization
## branch.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(BEDMatrix); library(ggplot2); library(patchwork)
})

DEC   <- "supplemental_data/deconvolution"
PANEL <- Sys.getenv("CENDR_PLINK_2021", "data/genotypes/CeNDR20210121_Plink")
CHROMS <- c("I", "II", "III", "IV", "V", "X")
BLK   <- 100000L
DIAG  <- "plots/diagnostics"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
if (!dir.exists(PANEL)) stop("2021 genotype panel not found at ", PANEL, call. = FALSE)

## --- fitness, with the original's two transforms ---------------------------
sd_ <- fread(cmd = paste("gzcat", shQuote(file.path(DEC, "simulation_seeded_frequencies.tsv.gz"))))
strains <- sort(unique(sd_$strain))
msg(sprintf("seeded simulation: %d strains, %d traits, %d depths, %d replicates",
            length(strains), uniqueN(sd_$trait), uniqueN(sd_$depth), uniqueN(sd_$replicate)))

## --- PART A: the zero block ------------------------------------------------
zf <- sd_[replicate == 1 & depth == max(depth),
          .(frac_zero = mean(input == 0)), by = trait]
r2 <- sd_[, {
  ok <- input > 0
  .(r2_all = cor(frequency, input)^2,
    r2_present = cor(frequency[ok], input[ok])^2,
    n_present = sum(ok))
}, by = .(trait, depth, replicate)]
R2 <- r2[, .(r2_all = mean(r2_all), r2_present = mean(r2_present),
             n_present = n_present[1]), by = .(trait, depth)]
R2 <- merge(R2, zf, by = "trait")
R2[, drop := r2_all - r2_present]
fwrite(R2, file.path(DIAG, "simulation_optimism_r2.tsv"), sep = "\t")

cat("\n== PART A: r-squared with the absent strains in, and out ==\n")
print(R2[, .(traits = .N, r2_all = round(mean(r2_all), 3),
             r2_present = round(mean(r2_present), 3),
             drop = round(mean(drop), 3)), by = depth][order(depth)])
## Dispersion of the PRESENT strains, which turns out to be what matters
disp <- sd_[replicate == 1 & depth == 10 & input > 0,
            .(cv = sd(input) / mean(input), log10_span = log10(max(input) / min(input))),
            by = trait]
R2 <- merge(R2, disp, by = "trait")
cat("\n== the control: does the drop track the size of the zero block? ==\n")
byt <- R2[depth == 10, .(trait, frac_zero = round(frac_zero, 3), n_present,
                         cv = round(cv, 2), log10_span = round(log10_span, 2),
                         r2_all = round(r2_all, 3), r2_present = round(r2_present, 3),
                         drop = round(drop, 3))][order(-drop)]
print(byt)
cat("\n  THE PREDICTION FAILED. The drop does not track how many strains are absent;\n")
cat("  it tracks how spread out the PRESENT ones are. A tight cluster of present\n")
cat("  strains has little variance of its own, so the distant zeros supply most of\n")
cat("  the range and removing them collapses the correlation. A well-spread trait\n")
cat("  barely notices. PC1 is 74.6%% absent and loses 0.026; etoposide is 58.4%%\n")
cat("  absent and loses 0.223.\n\n")
for (v in c("frac_zero", "n_present", "cv", "log10_span")) {
  ct <- suppressWarnings(cor.test(R2[depth == 10][[v]], R2[depth == 10]$drop, method = "spearman"))
  cat(sprintf("    Spearman(%-11s, drop) = %+.3f  p = %.3f\n", v, ct$estimate, ct$p.value))
}

## --- private markers in the 327-strain reference ---------------------------
PCACHE <- ".pool_opt_cache/sim2021_private.rds"
if (file.exists(PCACHE)) { PRIV <- readRDS(PCACHE); msg("private counts from cache") } else {
msg("counting private markers in the 327-strain reference")
npriv <- integer(length(strains))
for (chr in CHROMS) {
  bm <- BEDMatrix(file.path(PANEL, paste0(chr, ".bed")), simple_names = FALSE)
  keep <- match(paste0(strains, "_", strains), rownames(bm)); stopifnot(!anyNA(keep))
  for (s in seq(1L, ncol(bm), by = BLK)) {
    e <- min(s + BLK - 1L, ncol(bm))
    g <- bm[keep, s:e, drop = FALSE]; g[is.na(g)] <- 0; g[g == 2] <- 1
    cnt <- colSums(g)
    w <- which(cnt == 1)
    if (length(w)) npriv <- npriv + tabulate(max.col(t(g[, w, drop = FALSE]),
                                             ties.method = "first"),
                                             nbins = length(strains))
  }
  msg("  ", chr)
}
PRIV <- data.table(strain = strains, n_private = npriv)
saveRDS(PRIV, PCACHE) }
msg(sprintf("private markers per strain: median %.0f, range %d-%d, zero-private strains %d",
            median(PRIV$n_private), min(PRIV$n_private), max(PRIV$n_private),
            sum(PRIV$n_private == 0)))

## --- PART B: per-strain error against privateness --------------------------
## Only where the strain is actually IN the pool. An absent strain has error
## near zero by construction and would dilute the comparison to nothing.
E <- sd_[input > 0, .(abs_err = mean(abs(frequency - input)),
                      rel_err = mean(abs(frequency - input) / input),
                      mean_input = mean(input), n_obs = .N),
         by = .(strain, depth)]
E <- merge(E, PRIV, by = "strain")
fwrite(E, file.path(DIAG, "simulation_optimism_perstrain.tsv"), sep = "\t")

cat("\n== PART B: does privateness predict per-strain error in the CLEAN case? ==\n")
sp <- E[, {
  a <- suppressWarnings(cor.test(n_private, abs_err, method = "spearman"))
  r <- suppressWarnings(cor.test(n_private, rel_err, method = "spearman"))
  .(rho_abs = a$estimate, p_abs = a$p.value, rho_rel = r$estimate, p_rel = r$p.value,
    n_strains = .N)
}, by = depth][order(depth)]
print(sp[, .(depth, n_strains, rho_abs = round(rho_abs, 3), p_abs = signif(p_abs, 2),
             rho_rel = round(rho_rel, 3), p_rel = signif(p_rel, 2))])
cat("\n  adj R2 of rank(rel_err) ~ rank(n_private), by depth:\n")
for (d in sort(unique(E$depth)))
  cat(sprintf("    %3dx  %.3f\n", d,
      summary(lm(rank(rel_err) ~ rank(n_private), data = E[depth == d]))$adj.r.squared))
cat("\n  for comparison, real MIP-seq (DIAG_baugh_rmsd_predictors.R):\n")
cat("    rho = -0.42, adj R2 = 0.165 against per-strain RMSD\n")
cat("\n  THE TEST RETURNS NEGATIVE. Privateness explains about the same modest share\n")
cat("  of per-strain error in the perfectly specified simulation as it does against\n")
cat("  real MIP-seq. So the variance it leaves is NOT reference quality -- it is\n")
cat("  there even when the reference is exact. Most of it is abundance:\n\n")
E10 <- E[depth == 10]
cat(sprintf("    adj R2  privateness only  %.3f\n",
            summary(lm(rank(rel_err) ~ rank(n_private), E10))$adj.r.squared))
cat(sprintf("    adj R2  abundance only    %.3f\n",
            summary(lm(rank(rel_err) ~ rank(mean_input), E10))$adj.r.squared))
cat(sprintf("    adj R2  both              %.3f\n",
            summary(lm(rank(rel_err) ~ rank(n_private) + rank(mean_input), E10))$adj.r.squared))
cat("\n  Private-marker count is an incomplete measure of identifiability, not a\n")
cat("  measure of a broken reference. K* therefore constrains a real quantity but\n")
cat("  a partial one.\n")

## --- figure -----------------------------------------------------------------
theme_set(theme_bw(9) + theme(
  plot.title = element_text(face = "bold", size = 9.5),
  plot.subtitle = element_text(size = 7.2, colour = "grey30"),
  panel.grid.minor = element_blank(),
  strip.background = element_rect(fill = "grey93"),
  strip.text = element_text(size = 7, face = "bold")))

L <- melt(R2, id.vars = c("trait", "depth", "frac_zero"),
          measure.vars = c("r2_all", "r2_present"),
          variable.name = "set", value.name = "r2")
L[, set := factor(set, c("r2_all", "r2_present"),
                  c("all 327 strains (as reported)", "only strains present in the pool"))]
pA <- ggplot(L, aes(depth, r2, colour = set, group = interaction(trait, set))) +
  geom_line(linewidth = 0.4, alpha = 0.75) + geom_point(size = 0.8) +
  scale_x_log10(breaks = c(1, 3, 10, 30, 100, 500)) +
  scale_colour_manual(values = c("#C4302B", "#1B7837")) +
  facet_wrap(~ set) +
  labs(x = "simulated depth (x)", y = "r-squared",
       title = "Part A: the absent strains do inflate it, by 0.07 at 10x and 0.19 at 1x",
       subtitle = paste("One line per trait. Absent strains have fitness exactly 0 and NNLS returns exactly 0 for them, zero being",
                        "\na hard boundary of the feasible set. Dropping them takes the mean from 0.953 to 0.879 at 10x and from",
                        "\n0.741 to 0.547 at 1x -- real, but not the whole story, and not 'most of the recovery'."))

## the prediction was that the drop tracks the COUNT of zeros. It does not.
B <- melt(R2[depth == 10], id.vars = c("trait", "drop"),
          measure.vars = c("frac_zero", "cv"), variable.name = "x", value.name = "v")
B[, x := factor(x, c("frac_zero", "cv"),
                c("fraction of strains absent  (rho = -0.04, p = 0.96)",
                  "spread of the strains PRESENT, CV  (rho = -0.86, p = 0.02)"))]
pB <- ggplot(B, aes(v, drop)) +
  geom_point(colour = "#2E4057", size = 2) +
  geom_text(aes(label = trait), size = 1.9, vjust = -0.9, colour = "grey30") +
  facet_wrap(~ x, scales = "free_x") + expand_limits(y = c(0, 0.26)) +
  labs(x = NULL, y = "r-squared lost when absent strains are dropped",
       title = "The prediction failed, and the real driver is dispersion",
       subtitle = paste("I predicted the drop would track how MANY strains are absent. It does not (left). It tracks how spread out",
                        "\nthe present strains are (right): a tight cluster has little variance of its own, so the distant zeros supply",
                        "\nthe range and removing them collapses the fit. PC1 is 74.6% absent and loses 0.026."))

pC <- ggplot(E[depth %in% c(1, 10, 100)], aes(n_private + 1, rel_err)) +
  geom_point(colour = "#2E4057", alpha = 0.5, size = 1) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "#C4302B",
              linewidth = 0.5) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ paste0(depth, "x")) +
  labs(x = "private markers in the 327-strain reference (log10, +1)",
       y = "mean relative error, pool strains only",
       title = "Part B: privateness is no more predictive in a perfect world",
       subtitle = paste("The simulation solves with the same matrix it generated from, so the only error is binomial sampling.",
                        "\nPrivateness still explains only adj R2 0.21 here against 0.165 on real MIP-seq. The variance it leaves is",
                        "\nnot reference quality -- adding abundance takes the clean-case model from 0.21 to 0.47."))

fig <- (pA / (pB | pC)) + plot_layout(heights = c(1, 1.15))
ggsave(file.path(DIAG, "DIAG_pool_simulation_optimism.pdf"), fig, width = 11, height = 8)
ggsave(file.path(DIAG, "DIAG_pool_simulation_optimism.png"), fig, width = 11, height = 8, dpi = 200)
msg("wrote DIAG_pool_simulation_optimism.{pdf,png}")
