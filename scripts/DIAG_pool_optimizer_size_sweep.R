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

CACHE <- ".pool_opt_cache/sweep_panels.rds"
cached <- if (file.exists(CACHE) && !nzchar(Sys.getenv("POOL_OPT_REFIT"))) readRDS(CACHE) else list()
out <- list(); prev <- integer(0); panels <- list()
for (k in SIZES) {
  msg(sprintf("=== n = %d ===", k))
  ## The null records the floor as well as the objective, because it turns out
  ## NO random panel meets the floor at any size -- at n = 96 the median random
  ## panel's worst strain carries 180 private markers against K* = 1000, and at
  ## n = 250 it carries 32. So vs_null is not "optimised against unoptimised",
  ## it is "constrained optimum against UNCONSTRAINED random". The null panels
  ## are buying their markers with an identifiability they do not have, which is
  ## why the ratio can fall below 1 at sizes where the floor binds hard.
  nulls <- replicate(NRAND, { Sr <- sample.int(N_UNIV, k)
    c(score(Sr)$n_mappable, min(priv_counts(Sr)$nondiv[Sr])) })
  nullk <- nulls[1, ]; null_floor <- nulls[2, ]
  S <- pool_opt_seed(k, fixed = prev)          # nested in the previous solution
  S <- pool_opt_repair(S, uni$nondiv)
  reached <- score(S)$min_nondiv
  ok <- reached >= KSTAR
  if (!ok) msg(sprintf("  INFEASIBLE at K* = %d: best achievable floor is %d", KSTAR, reached))
  ## Anneal at EVERY size. An earlier version skipped the search wherever K*
  ## could not be met, which made those sizes report a bare farthest-point seed
  ## and put them below the random mean -- an artefact of the code path, not a
  ## statement about optimisation. Where the floor is unreachable the search
  ## runs against the best floor that IS reachable, so vs_null stays a measure
  ## of the search at every size; `floor_used` records which floor applied.
  floor_used <- if (ok) KSTAR else reached
  KSTAR_SAVE <- KSTAR; KSTAR <<- floor_used
  key <- sprintf("n%d_i%d", k, ITERS)
  z <- if (!is.null(cached[[key]])) {
         msg("  reusing cached panel; POOL_OPT_REFIT=1 to refit")
         list(S = cached[[key]], best = score(cached[[key]]))
       } else pool_opt_anneal(S, ITERS)
  KSTAR <<- KSTAR_SAVE
  cached[[key]] <- z$S
  S <- z$S; s <- z$best
  prev <- S; panels[[as.character(k)]] <- U[S]
  out[[length(out) + 1]] <- data.table(
    n = k, feasible = ok, floor_used = floor_used,
    n_mappable = s$n_mappable, n_mac = s$n_mac,
    frac_mappable = s$n_mappable / s$n_mac,
    null_mean = mean(nullk), vs_null = s$n_mappable / mean(nullk),
    null_med_floor = median(null_floor), null_frac_feasible = mean(null_floor >= KSTAR),
    min_nondiv = s$min_nondiv, med_nondiv = s$med_nondiv,
    pc1_share = s$pc1_share, eff_dim = s$eff_dim,
    med_div_mb = median(db[match(U[S], strain), div_mb], na.rm = TRUE))
  msg(sprintf("  %d mappable (%.3fx null), floor %d, PC1 %.3f",
              s$n_mappable, s$n_mappable / mean(nullk), s$min_nondiv, s$pc1_share))
}
saveRDS(cached, CACHE)
R <- rbindlist(out)
R[, per_strain := n_mappable / n]
R[, marginal := c(NA_real_, diff(n_mappable) / diff(n))]
fwrite(R, file.path(DIAG, "pool_optimizer_size_sweep.tsv"), sep = "\t")
saveRDS(list(R = R, panels = panels), file.path(DIAG, "pool_optimizer_size_sweep.rds"))
cat("\n== size sweep ==\n")
cat(sprintf("random panels meeting K* = %d: %s (out of %d draws at each size)\n", KSTAR,
    paste(unique(sprintf("%.0f%%", 100 * R$null_frac_feasible)), collapse = ", "), NRAND))
print(R[, .(n, feasible, n_mappable, vs_null = round(vs_null, 3), null_med_floor,
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
       title = "Constrained optimum against UNCONSTRAINED random",
       subtitle = paste(sprintf(
         "Not a like-for-like ratio: almost no random panel meets the identifiability floor -- %s of %d draws at n = %d and %s everywhere",
         scales::percent(R$null_frac_feasible[1], accuracy = 1), NRAND, R$n[1],
         if (all(R$null_frac_feasible[-1] == 0)) "none" else "almost none"),
         sprintf("\nlarger, where the median random panel's worst strain carries %.0f private markers at n = %d and %.0f at n = %d.",
                 R$null_med_floor[R$n == 96], 96,
                 R$null_med_floor[nrow(R)], R$n[nrow(R)]),
         "\nThe null buys its markers with an identifiability it does not have, so a ratio below 1 means the constraint costs more",
         "\nthan the search gains -- not that the search failed. Open circles are sizes where even the constrained optimum misses K*."))

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
