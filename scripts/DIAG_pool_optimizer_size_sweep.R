## How large should the panel be? ---------------------------------------------
##
##   Rscript scripts/DIAG_pool_optimizer_size_sweep.R [iters_per_size]
##     -> plots/diagnostics/DIAG_pool_optimizer_size_sweep.{pdf,png}
##     -> plots/diagnostics/pool_optimizer_size_sweep.tsv
##
## DIAG_pool_optimizer.R fixes the panel at 96 strains and optimises membership.
## Size is not given, though -- it is the other half of the design, and the two
## halves pull against each other:
##
##   more strains  -> more testable variation, so more mappable markers, and a
##                    larger association sample
##   more strains  -> every strain is harder to identify, because privateness is
##                    destroyed by adding a near neighbour, and at fixed
##                    sequencing budget each strain is also read less deeply
##
## So there should be a size beyond which the identifiability floor cannot be
## met at all, and, before that, a size beyond which each added strain buys
## little. This finds both.
##
## THREE QUANTITIES, and they answer different questions.
##   n_mappable    absolute, and it rises with size almost by construction --
##                 a bigger panel segregates more markers. On its own this
##                 would say "as large as possible", which is not useful.
##   vs_null       n_mappable over the mean of random panels OF THE SAME SIZE.
##                 This is what optimisation is worth, with the size effect
##                 divided out, and it is the honest measure of the search.
##   min_nondiv    the binding constraint. Where this falls below K* the panel
##                 is infeasible however good the other two look.
##
## Each size is warm-started from the previous size's solution, which both
## converges faster than a cold start and makes the series comparable -- the
## panels are nested rather than independently found.
##
## Every point is a LOWER BOUND: the per-size annealing budget is a fraction of
## what DIAG_pool_optimizer.R spends on 96 alone, and that run was still
## climbing when it stopped.
##
## Exploratory, on the pool_optimization branch.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

a <- commandArgs(TRUE)
ITERS   <- if (length(a) >= 1) as.integer(a[1]) else 1200L
SIZES   <- c(32L, 48L, 64L, 80L, 96L, 128L, 160L, 200L, 250L)
KSTAR   <- 1000L
MAC_MIN <- 10L
R2_MAX  <- 0.5
Q       <- 5L
NRAND   <- 150L
DIAG    <- "plots/diagnostics"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)
set.seed(1)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

P <- readRDS(".pool_opt_cache/precomputed.rds")
U <- P$strains
source("scripts/pool_optimizer_core.R")
pool_opt_attach(P)
uni <- priv_counts(seq_len(N_UNIV))
DIVBED <- Sys.getenv("CENDR_DIVERGENT",
  "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213/20231213_c_elegans_divergent_regions_strain.bed")
db <- fread(DIVBED, col.names = c("chrom", "start", "end", "strain"))[
  , .(div_mb = sum(end - start) / 1e6), by = strain]

out <- list(); prev <- integer(0); panels <- list()
for (k in SIZES) {
  msg(sprintf("=== n = %d ===", k))
  nullk <- replicate(NRAND, score(sample.int(N_UNIV, k))$n_mappable)
  S <- pool_opt_seed(k, fixed = prev)          # nested in the previous solution
  S <- pool_opt_repair(S, uni$nondiv)
  reached <- score(S)$min_nondiv
  ok <- reached >= KSTAR
  if (!ok) msg(sprintf("  INFEASIBLE at K* = %d: best achievable floor is %d", KSTAR, reached))
  ## anneal only inside the feasible set; if the floor cannot be met there is
  ## nothing to search over, so the seed is reported as-is
  z <- if (ok) pool_opt_anneal(S, ITERS) else list(S = S, best = score(S))
  S <- z$S; s <- z$best
  prev <- S; panels[[as.character(k)]] <- U[S]
  out[[length(out) + 1]] <- data.table(
    n = k, feasible = ok, n_mappable = s$n_mappable, n_mac = s$n_mac,
    frac_mappable = s$n_mappable / s$n_mac,
    null_mean = mean(nullk), vs_null = s$n_mappable / mean(nullk),
    min_nondiv = s$min_nondiv, med_nondiv = s$med_nondiv,
    pc1_share = s$pc1_share, eff_dim = s$eff_dim,
    med_div_mb = median(db[match(U[S], strain), div_mb], na.rm = TRUE))
  msg(sprintf("  %d mappable (%.3fx null), floor %d, PC1 %.3f",
              s$n_mappable, s$n_mappable / mean(nullk), s$min_nondiv, s$pc1_share))
}
R <- rbindlist(out)
R[, per_strain := n_mappable / n]
R[, marginal := c(NA_real_, diff(n_mappable) / diff(n))]
fwrite(R, file.path(DIAG, "pool_optimizer_size_sweep.tsv"), sep = "\t")
saveRDS(list(R = R, panels = panels), file.path(DIAG, "pool_optimizer_size_sweep.rds"))
cat("\n== size sweep ==\n")
print(R[, .(n, feasible, n_mappable, vs_null = round(vs_null, 3),
            min_nondiv, med_nondiv, pc1_share = round(pc1_share, 3),
            eff_dim = round(eff_dim, 1), med_div_mb = round(med_div_mb, 2),
            per_strain = round(per_strain, 1), marginal = round(marginal, 1))])

## --- figure -----------------------------------------------------------------
theme_set(theme_bw(9) + theme(
  plot.title = element_text(face = "bold", size = 9.5),
  plot.subtitle = element_text(size = 7.2, colour = "grey30"),
  panel.grid.minor = element_blank(), legend.position = "none"))
GRN <- "#1B7837"; RED <- "#C4302B"
mark <- function(p) p + geom_vline(xintercept = 96, linetype = "dotted",
                                   linewidth = 0.4, colour = "grey55")

p1 <- mark(ggplot(R, aes(n, n_mappable / 1e3)) +
  geom_line(aes(y = null_mean / 1e3), colour = "grey65", linewidth = 0.5) +
  geom_line(colour = GRN, linewidth = 0.6) +
  geom_point(aes(shape = feasible), colour = GRN, size = 1.8) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1))) +
  labs(x = "panel size", y = "mappable markers (thousands)",
       title = "Absolute yield rises with size",
       subtitle = "Grey is the random-panel mean at each size. On its own this axis just says 'as large as possible'.")

p2 <- mark(ggplot(R, aes(n, vs_null)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
  geom_line(colour = GRN, linewidth = 0.6) +
  geom_point(aes(shape = feasible), colour = GRN, size = 1.8) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1))) +
  labs(x = "panel size", y = "mappable markers / random panel of same size",
       title = "What the optimisation is worth, size divided out",
       subtitle = "Open circles are sizes where the identifiability floor could not be met, so no search was run.")

p3 <- mark(ggplot(R, aes(n, min_nondiv)) +
  geom_hline(yintercept = KSTAR, linetype = "dashed", linewidth = 0.4, colour = RED) +
  geom_line(colour = GRN, linewidth = 0.6) +
  geom_point(aes(shape = feasible), colour = GRN, size = 1.8) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  scale_y_log10()) +
  labs(x = "panel size", y = "worst strain's private markers\noutside divergent regions",
       title = "The constraint, and where it breaks",
       subtitle = sprintf("Dashed line is K* = %d. Adding strains destroys privateness, so this can only fall.", KSTAR))

p4 <- mark(ggplot(R[-1], aes(n, marginal)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_line(colour = GRN, linewidth = 0.6) + geom_point(colour = GRN, size = 1.8)) +
  labs(x = "panel size", y = "extra mappable markers\nper strain added",
       title = "Diminishing returns on size",
       subtitle = "Slope between consecutive sizes. Dotted line marks the 96 used in DIAG_pool_optimizer.R.")

fig <- (p1 | p2) / (p3 | p4)
ggsave(file.path(DIAG, "DIAG_pool_optimizer_size_sweep.pdf"), fig, width = 10.5, height = 7)
ggsave(file.path(DIAG, "DIAG_pool_optimizer_size_sweep.png"), fig, width = 10.5, height = 7, dpi = 200)
msg("wrote DIAG_pool_optimizer_size_sweep.{pdf,png}")
