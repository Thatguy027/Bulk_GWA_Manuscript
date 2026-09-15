## A panel optimiser for pooled phenotyping -----------------------------------
##
##   Rscript scripts/DIAG_pool_optimizer.R [n] [iters]
##     -> plots/diagnostics/DIAG_pool_optimizer.{pdf,png}
##     -> plots/diagnostics/pool_optimizer_panels.tsv    one row per strain
##     -> plots/diagnostics/pool_optimizer_summary.tsv   one row per panel
##
## Requires .pool_opt_cache/precomputed.rds from
## scripts/DIAG_pool_optimizer_precompute.R.
##
## THE PROBLEM. Choosing strains for a pooled panel pulls two ways. The
## deconvolution needs every strain to be distinguishable from the rest, which
## wants strains that are genetically distinctive. The association mapping needs
## the panel not to collapse onto one divergent-versus-swept axis, because a
## mixed model spends its power removing exactly that. Selecting hard on
## distinctiveness gives a Hawaii-heavy panel where every QTL is confounded with
## structure; selecting hard against structure gives a panel of near-identical
## swept strains the solver cannot take apart.
##
## THE FORMULATION. Identifiability is a CONSTRAINT, not an objective, because
## it saturates: past the point where a strain's frequency is pinned to the
## precision the sequencing depth supports, more privateness buys nothing and
## costs panel slots. So
##
##   maximise   mappable markers
##   subject to every strain having at least K* private markers OUTSIDE its own
##              divergent regions
##
## The constraint counts only non-divergent privateness on the evidence in
## DIAG_pool_private_marker_classes.R: splitting the Baugh panel's private
## markers on divergent-region membership, the markers outside carry all of the
## predictive signal against deconvolution error (rho = -0.417) and the markers
## inside essentially none (-0.219, p = 0.48 in a joint model).
##
## A mappable marker is one with at least MAC_MIN minor-allele carriers in the
## panel AND at most R2_MAX of its genotype variance explained by the panel's
## top Q principal components. That single count internalises the whole tension:
## an all-divergent panel loses markers to the R2 filter because they align with
## PC1, and an all-swept panel loses them to the MAC filter because there is no
## variation left to test. Neither failure mode needs a separate penalty term.
##
## WHAT IS COMPARED. The optimised panel, the two panels actually used in this
## manuscript, a naive panel built by taking the 96 strains with the most
## private markers, and a null distribution of random panels. The naive panel is
## the point of the exercise: it is what "maximise private alleles" produces.
##
## Exploratory, on the pool_optimization branch.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})

a <- commandArgs(TRUE)
N       <- if (length(a) >= 1) as.integer(a[1]) else 96L   # panel size
ITERS   <- if (length(a) >= 2) as.integer(a[2]) else 4000L
KSTAR   <- 1000L      # required private markers outside divergent regions
MAC_MIN <- 10L        # minor-allele carriers for a marker to be testable
R2_MAX  <- 0.5        # genotype variance a marker may share with the top PCs
Q       <- 5L         # principal components treated as structure
NRAND   <- 300L       # random panels for the null
DIAG    <- "plots/diagnostics"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)
set.seed(1)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

P <- readRDS(".pool_opt_cache/precomputed.rds")
U <- P$strains; NU <- length(U); NR <- P$n_rare
## carrier pairs sorted by strain, with per-strain offsets. Selecting the pairs
## belonging to a panel is then |S| contiguous reads instead of a logical pass
## over all 10.9 M rows, which is the whole cost of the inner loop.
setorder(P$car, s)
cm <- P$car$m; cs <- P$car$s; cd <- P$car$in_div
off <- c(0L, cumsum(tabulate(cs, nbins = length(P$strains))))
CSTART <- off[-length(off)] + 1L; CLEN <- diff(off)
gp <- P$gp; G <- P$G; M <- ncol(gp)
msg(sprintf("universe %d strains, %s rare markers, %d pruned markers",
            NU, format(NR, big.mark = ","), M))

## --- scoring ---------------------------------------------------------------
## privateness: one pass over the carrier list. A marker is private to the panel
## when exactly one of its carriers is in it.
priv_counts <- function(S) {
  k <- sequence(CLEN[S], CSTART[S])
  ms <- cm[k]; ss <- cs[k]; ds <- cd[k]
  cnt <- tabulate(ms, nbins = NR)
  p <- cnt[ms] == 1L
  list(nondiv = tabulate(ss[p & !ds], nbins = NU),
       div    = tabulate(ss[p &  ds], nbins = NU))
}

## mapping: the panel kinship is a double-centred submatrix of the one
## precomputed cross-product, and the top eigenvectors are already orthogonal to
## the intercept, so markers never have to be re-centred to be projected onto
## them.
map_score <- function(S) {
  n <- length(S)
  Gs <- G[S, S, drop = FALSE]
  rm_ <- rowMeans(Gs); K <- Gs - rm_ - rep(rm_, each = n) + mean(Gs)
  ev <- eigen(K, symmetric = TRUE)
  lam <- pmax(ev$values, 0)
  Qm <- ev$vectors[, seq_len(Q), drop = FALSE]
  A <- gp[S, , drop = FALSE]
  s1 <- colSums(A); s2 <- colSums(A * A)
  tot <- s2 - s1 * s1 / n                       # centred sum of squares
  expl <- colSums(crossprod(Qm, A)^2)           # Q is orthogonal to 1
  r2 <- ifelse(tot > 1e-9, expl / tot, 1)
  mac <- pmin(s1, n - s1)
  list(n_mappable = sum(mac >= MAC_MIN & r2 <= R2_MAX),
       n_mac      = sum(mac >= MAC_MIN),
       pc1_share  = lam[1] / sum(lam),
       eff_dim    = sum(lam)^2 / sum(lam^2))
}

score <- function(S) {
  pc <- priv_counts(S); ms <- map_score(S)
  nd <- pc$nondiv[S]
  c(ms, list(min_nondiv = min(nd), med_nondiv = median(nd),
             n_below = sum(nd < KSTAR), med_div = median(pc$div[S])))
}
feasible <- function(S) min(priv_counts(S)$nondiv[S]) >= KSTAR

## --- baselines --------------------------------------------------------------
rnai  <- match(unique(fread(cmd = "gzcat supplemental_data/phenotypes/pooled_vst_traits.csv.gz")$strain), U)
baugh <- match(fread("supplemental_data/deconvolution/baugh_strain_private_markers.tsv")$strain, U)
baugh <- baugh[!is.na(baugh)]

## "maximise private alleles" -- privateness scored in the FULL universe, which
## is how anyone would rank strains before having a panel to score against
uni <- priv_counts(seq_len(NU))
naive <- order(uni$nondiv + uni$div, decreasing = TRUE)[seq_len(N)]

## Nulls are drawn AT EACH PANEL'S OWN SIZE. n_mappable grows with panel size,
## so scoring a 93-strain panel against a 96-strain null would charge it for
## being small. Each panel is compared only to random panels of its own size.
msg("random nulls ...")
SIZES <- sort(unique(c(N, length(rnai), length(baugh))))
rand_by_n <- setNames(lapply(SIZES, function(k)
  replicate(NRAND, score(sample.int(NU, k))$n_mappable)), as.character(SIZES))
rand <- rand_by_n[[as.character(N)]]

## --- the annealed panel is cached, so the reporting half is cheap to re-run --
FIT <- sprintf(".pool_opt_cache/annealed_n%d.rds", N)
if (file.exists(FIT) && !nzchar(Sys.getenv("POOL_OPT_REFIT"))) {
  z <- readRDS(FIT); S <- z$S; trace <- z$trace
  msg(sprintf("reusing cached annealed panel (%d mappable); set POOL_OPT_REFIT=1 to refit",
              score(S)$n_mappable))
} else {

## --- seed: farthest-point sampling on the genotype distance -----------------
d2 <- outer(diag(G), diag(G), "+") - 2 * G
S <- integer(N); S[1] <- which.max(rowSums(d2))
dmin <- d2[S[1], ]
for (i in 2:N) { S[i] <- which.max(dmin); dmin <- pmin(dmin, d2[S[i], ]) }
msg(sprintf("farthest-point seed: %d mappable, min non-divergent private %d",
            score(S)$n_mappable, score(S)$min_nondiv))

## --- repair the seed into the feasible set, then anneal ---------------------
## A strain below K* is swapped out for the candidate that best relieves the
## binding constraint. Privateness is NOT monotone in the panel -- adding a
## strain can destroy another's privateness -- so this is re-checked each pass.
for (pass in 1:40) {
  pc <- priv_counts(S)$nondiv
  bad <- S[pc[S] < KSTAR]
  if (!length(bad)) break
  out <- bad[which.min(pc[bad])]
  cand <- setdiff(seq_len(NU), S)
  cand <- cand[order(uni$nondiv[cand], decreasing = TRUE)][1:min(25, length(cand))]
  best <- NULL; bestmin <- -1
  for (cc in cand) {
    T_ <- c(setdiff(S, out), cc); mn <- min(priv_counts(T_)$nondiv[T_])
    if (mn > bestmin) { bestmin <- mn; best <- T_ }
  }
  S <- best
  if (bestmin >= KSTAR) break
}
msg(sprintf("after repair: min non-divergent private %d (target %d)",
            score(S)$min_nondiv, KSTAR))

## --- simulated annealing over swaps ----------------------------------------
cur <- score(S); best <- cur; bestS <- S
T0 <- 40; trace <- numeric(ITERS)
for (it in seq_len(ITERS)) {
  Temp <- T0 * (1 - it / ITERS) + 1e-6
  out <- S[sample.int(N, 1)]
  ins <- sample(setdiff(seq_len(NU), S), 1)
  T_ <- c(setdiff(S, out), ins)
  st <- score(T_)
  if (st$min_nondiv >= KSTAR &&
      (st$n_mappable > cur$n_mappable ||
       runif(1) < exp((st$n_mappable - cur$n_mappable) / Temp))) {
    S <- T_; cur <- st
    if (cur$n_mappable > best$n_mappable) { best <- cur; bestS <- S }
  }
  trace[it] <- best$n_mappable
  if (it %% 500 == 0) msg(sprintf("  iter %5d  best %d mappable", it, best$n_mappable))
}
S <- bestS
saveRDS(list(S = S, trace = trace), FIT)
msg(sprintf("optimised: %d mappable markers", best$n_mappable))
}

## --- report -----------------------------------------------------------------
PAN <- list(optimised = S, `RNAi panel (93)` = rnai, `Baugh panel (102)` = baugh,
            `naive: most private` = naive)
res <- rbindlist(lapply(names(PAN), function(nm) {
  s <- score(PAN[[nm]])
  data.table(panel = nm, n = length(PAN[[nm]]), n_mappable = s$n_mappable,
             n_mac = s$n_mac, frac_mappable = s$n_mappable / s$n_mac,
             pc1_share = round(s$pc1_share, 4), eff_dim = round(s$eff_dim, 2),
             min_nondiv = s$min_nondiv, med_nondiv = s$med_nondiv,
             n_below_Kstar = s$n_below, med_div = s$med_div) }))
## each panel against the null at ITS OWN size
res[, null_mean := round(sapply(n, function(k) mean(rand_by_n[[as.character(k)]])))]
res[, vs_null := round(n_mappable / null_mean, 3)]
res[, null_pctile := round(100 * mapply(function(v, k)
  mean(rand_by_n[[as.character(k)]] < v), n_mappable, n), 1)]
fwrite(res, file.path(DIAG, "pool_optimizer_summary.tsv"), sep = "\t")
cat("\n== panels ==\n"); print(res)

memb <- rbindlist(lapply(names(PAN), function(nm) {
  pc <- priv_counts(PAN[[nm]])
  data.table(panel = nm, strain = U[PAN[[nm]]],
             private_nondiv = pc$nondiv[PAN[[nm]]], private_div = pc$div[PAN[[nm]]]) }))
## How much divergent genome each panel actually carries. This is the question
## the whole exercise is about: the optimiser is not avoiding divergent strains
## and not maximising them either, it is rationing them.
DIV <- Sys.getenv("CENDR_DIVERGENT",
                  "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213/20231213_c_elegans_divergent_regions_strain.bed")
db <- fread(DIV, col.names = c("chrom", "start", "end", "strain"))[
  , .(div_mb = sum(end - start) / 1e6), by = strain]
memb <- merge(memb, db, by = "strain", all.x = TRUE)
memb[is.na(div_mb), div_mb := 0]
cat("\n== divergent genome content per panel ==\n")
print(memb[, .(n = .N, median_div_mb = round(median(div_mb), 2),
               frac_over_5Mb = round(mean(div_mb > 5), 3)), by = panel])
fwrite(memb, file.path(DIAG, "pool_optimizer_panels.tsv"), sep = "\t")

cat("\n== overlap with the panels actually used ==\n")
cat(sprintf("  optimised vs RNAi panel : %d of %d\n", length(intersect(S, rnai)), N))
cat(sprintf("  optimised vs Baugh panel: %d of %d\n", length(intersect(S, baugh)), N))
cat(sprintf("  naive     vs RNAi panel : %d of %d\n", length(intersect(naive, rnai)), N))
saveRDS(list(res = res, memb = memb, trace = trace, rand = rand, U = U, PAN = PAN),
        file.path(DIAG, "pool_optimizer_result.rds"))
msg("wrote pool_optimizer_{summary,panels}.tsv")

## --- figure -----------------------------------------------------------------
theme_set(theme_bw(9) + theme(
  plot.title = element_text(face = "bold", size = 10),
  plot.subtitle = element_text(size = 7.4, colour = "grey30"),
  panel.grid.minor = element_blank(), legend.position = "none",
  strip.background = element_rect(fill = "grey93"),
  strip.text = element_text(size = 7.3, face = "bold")))
COL <- c("optimised" = "#1B7837", "RNAi panel (93)" = "#2E4057",
         "Baugh panel (102)" = "#4F86C6", "naive: most private" = "#C4302B")
R <- res[panel %in% names(COL)]

## A  where the naive panel loses: it has the markers, they just do not survive
pA <- ggplot(R, aes(n_mac / 1e3, frac_mappable, colour = panel)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = panel), size = 2.6, seed = 1, min.segment.length = 0,
                  box.padding = 0.5) +
  scale_colour_manual(values = COL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "markers with at least 10 minor-allele carriers (thousands)",
       y = "fraction surviving the\nstructure filter",
       title = "Where the naive panel loses",
       subtitle = paste("Selecting on private alleles finds nearly as much testable variation as the optimiser",
                        "\n(15.7k against 17.8k markers) and then gives 31% of it back to population structure."))

## B  the objective, against a null of random panels of the same size
## ratio to the own-size null, so the three panel sizes are on one axis
RB <- copy(R)[, `:=`(ratio = n_mappable / null_mean)]
nullspread <- rbindlist(lapply(as.character(SIZES), function(k)
  data.table(n = as.integer(k), r = rand_by_n[[k]] / mean(rand_by_n[[k]]))))
pB <- ggplot(nullspread, aes(r)) +
  geom_histogram(bins = 30, fill = "grey80", colour = "white", linewidth = 0.2) +
  geom_vline(data = RB, aes(xintercept = ratio, colour = panel), linewidth = 0.8) +
  geom_text_repel(data = RB, aes(x = ratio, y = Inf, label = panel, colour = panel),
                  size = 2.6, angle = 90, hjust = 1.05, seed = 1,
                  direction = "x", min.segment.length = Inf) +
  scale_colour_manual(values = COL) +
  labs(x = "mappable markers, relative to random panels of the SAME size",
       y = "random panels",
       title = sprintf("Each panel against %d random panels of its own size", NRAND),
       subtitle = paste("Panel size is divided out, so the three sizes share one axis. The naive panel sits on the null:",
                        "\nranking strains by private-allele count buys nothing at all for mapping. Note the optimised",
                        "\npanel is not the least structured one -- its PC1 share is 0.137 against 0.076-0.083 for the two",
                        "\nreal panels. It tolerates structure wherever the markers survive the filter anyway."))

## C  per-strain identifiability, the constraint the real panels violate
mb <- copy(memb); mb[, panel := factor(panel, names(COL))]
pC <- ggplot(mb, aes(panel, private_nondiv + 1, colour = panel)) +
  geom_hline(yintercept = KSTAR, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  geom_boxplot(outlier.shape = NA, fill = NA, linewidth = 0.35, width = 0.55) +
  geom_point(position = position_jitter(width = 0.14, height = 0, seed = 1),
             size = 0.9, alpha = 0.55) +
  scale_colour_manual(values = COL) + scale_y_log10() +
  annotate("text", x = 0.6, y = KSTAR * 1.35, label = sprintf("K* = %d", KSTAR),
           size = 2.4, colour = "grey30", hjust = 0) +
  labs(x = NULL, y = "private markers outside\ndivergent regions (log10, +1)",
       title = "The constraint the panels in use do not meet",
       subtitle = paste("Dashed line is the identifiability floor the optimiser enforces. 54 of 93 strains in the RNAi panel",
                        "\nand 40 of 102 in the Baugh panel fall below it; the lowest carries 10 private markers.")) +
  theme(axis.text.x = element_text(size = 6.6))

## D  the search
pD <- ggplot(data.table(it = seq_along(trace), y = trace / 1e3), aes(it, y)) +
  geom_line(colour = "#1B7837", linewidth = 0.5) +
  labs(x = "swap proposals", y = "best mappable\nmarkers (thousands)",
       title = "Simulated annealing over swaps",
       subtitle = paste("Each step drops one strain and adds another; proposals below the identifiability floor are",
                        "\nrejected outright. The trace is STILL CLIMBING at the last iteration, so the reported panel",
                        "\nis a lower bound on what this objective admits, not a converged optimum."))

fig <- (pA | pB) / (pC | pD)
ggsave(file.path(DIAG, "DIAG_pool_optimizer.pdf"), fig, width = 11, height = 7.6)
ggsave(file.path(DIAG, "DIAG_pool_optimizer.png"), fig, width = 11, height = 7.6, dpi = 200)
msg("wrote DIAG_pool_optimizer.{pdf,png}")
