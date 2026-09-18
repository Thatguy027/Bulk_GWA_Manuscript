## A simulated null for the cross scans, and what it says about interval width --
##
##   Rscript scripts/DIAG_cross_null_simulation.R
##     -> plots/diagnostics/cross_null_simulation.tsv
##     -> plots/diagnostics/cross_interval_calibration.tsv
##     -> plots/diagnostics/DIAG_cross_null_simulation.{pdf,png}
##
##   CROSS_SIM_NULL=200   null replicates       (default 100)
##   CROSS_SIM_QTL=60     replicates per effect (default 40)
##   CROSS_SIM_CORES=8    parallel workers      (default detectCores()-1, max 8)
##   CROSS_SIM_REFIT=1    ignore the cache and rerun
##
## WHY THIS EXISTS. The shipped support intervals are a fractional LOD drop --
## 10% of each chromosome's own peak in the interval tables, 5% in Figure 2 --
## and neither has a coverage property. The config comment says where the 10%
## came from: it was tuned on the Nov2024 chromosome III QTL, in two contrasts,
## so that sid-2 would fall inside (it needed 5.9% and 9.9%). That is a
## containment rule fitted with the answer known, on one locus. Its coverage is
## 100% there by construction and undefined everywhere else. The genome-wide
## threshold has the same problem from the other end: LOD 3.57 comes from
## effective.n.tests: 2000, a declared constant over 522,226 correlated markers.
##
## Neither can be calibrated from the data, because THE DESIGN HAS NO REPLICATE.
## T1 and T2 are not parallel pools: gen-11 worms seeded the RNAi and control
## conditions (T1), gen-12 progeny of those seeded T2, and gen-13 progeny were
## sequenced. T2 descends from T1. A T1-vs-T2 contrast of one condition carries
## a generation of real selection, so it is not a null, and the only genuine
## replication anywhere in the design is cross against cross.
##
## So the null has to be simulated, and it can be: AlphaSimR is already a
## pipeline dependency and every parameter of the cross is known.
##
## WHAT IS SIMULATED. Two inbred founders, opposite at every tracked locus, then
##   10 generations of random intercrossing at N = 12,500
##   split into an RNAi arm and a control arm at generation 11
##   2 further generations in each arm, the last expanded to 30,000
##   sequencing at the observed per-bin depth, then the pipeline's own estimator
## Loci are a 10 kb grid, which is the analysis resolution -- calcAFD bins at
## bin.width: 10000 -- so there is nothing to gain from simulating all 522,226
## markers. The map is the shipped one rescaled by map.expansion.factor: 0.2.
## THAT RESCALING IS CORRECT AND IS NOT A BUG: the shipped map spans 1,584 cM,
## which is the F10 AIL map; x0.2 gives 317 cM, the per-meiosis C. elegans map,
## and a per-meiosis map is exactly what a simulator that then runs the ten
## generations itself must be given. The AIL expansion is an emergent property
## of the repeated meioses, not an input.
##
## WHAT THE ESTIMATOR ACTUALLY DOES, which had to be established before the
## sequencing layer could be written. afd is NOT a per-bin quantity: within a
## 10 kb bin its SD is 0.0001 against a reported standard error of 0.0046, so
## calcAFD is a moving kernel. Reads per 10 kb bin are 1,316 (control) and 1,004
## (pos-1), and the shipped ndepth is 13,341 -- about ten bins' worth -- so the
## kernel spans roughly 100 kb and ndepth is the read count inside it, not a
## genome-wide total. min(nindv, ndepth) = min(9500, 13341) = 9500 is therefore
## right, and the individuals really are limiting. Two consequences:
##   the reported standard error prices SAMPLING correctly, and
##   no interval narrower than the ~100 kb kernel means anything as shipped,
## which on its own disqualifies the 0.05 Mb chromosome III interval the JU
## cross reports.
##
## WHAT THE STANDARD ERROR DOES NOT PRICE is drift after the arms separate. That
## is two generations, independent between arms; the ten AIL generations before
## the split are shared and cancel in the contrast. The simulation puts the
## post-split contribution at SD 0.0053 against a reported contrast standard
## error of 0.0066, so the honest figure is about 1.3x too small at census Ne.
## Ne below census -- RNAi killing, the transfer bottleneck, variance in
## reproductive success -- makes it worse, and Ne is the one parameter here that
## is not known.
##
## Reads the cross exports under data/ and the external map, so it is a DIAG
## script. Nothing in the manuscript reads its output.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(AlphaSimR); library(ggplot2); library(patchwork)
})

DIAG <- "plots/diagnostics"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)
CACHE <- file.path(DIAG, "cross_null_simulation_cache.rds")
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## --- the design, as reported by the experimenter -----------------------------
N_AIL     <- 12500L   # worms seeding each generation after generation 3
AIL_GEN   <- 10L      # generations of intercrossing before the conditions split
POST_GEN  <- 2L       # gen 11 -> gen 13, the sequenced generation
N_SEQ     <- 30000L   # gen-13 progeny sequenced
GRID      <- 1e4      # locus spacing = calcAFD's bin.width
EXPANSION <- 0.2      # map.expansion.factor: F10 AIL map -> per-meiosis map
NINDV     <- 9500L    # sample.size 10000 x sel.strength 0.95, the declared n
READS_BIN <- c(ctrl = 1316, rnai = 1004)   # measured medians, 10 kb bins
KERNEL_KB <- 100      # moving-kernel width implied by ndepth / reads-per-bin
THR_SHIP  <- 3.57     # the shipped genome-wide threshold

GMAP <- Sys.getenv("XQTL_GMAP",
  "/Users/Stefan/github_repos/xQTLSims/data/geneticMapXQTLsnplist.rds")
if (!file.exists(GMAP)) stop("genetic map not found: ", GMAP,
                             "\nset XQTL_GMAP to point at it", call. = FALSE)

NREP_NULL <- as.integer(Sys.getenv("CROSS_SIM_NULL", "100"))
NREP_QTL  <- as.integer(Sys.getenv("CROSS_SIM_QTL",  "40"))
CORES     <- as.integer(Sys.getenv("CROSS_SIM_CORES",
               as.character(min(8L, max(1L, parallel::detectCores() - 1L)))))
REFIT     <- nzchar(Sys.getenv("CROSS_SIM_REFIT"))
## selection coefficients against the susceptible homozygote, additive.
## pos-1 RNAi is embryonic-lethal, so the strong end of this range is the
## realistic one and the weak end is there to show the rules breaking.
S_SWEEP   <- c(0.15, 0.30, 0.60, 1.00)

## --- the locus grid and the per-meiosis map ---------------------------------
GM <- as.data.table(readRDS(GMAP))
grid <- GM[, {
  p <- seq(min(pos), max(pos), by = GRID)
  .(pos = p, cM = approx(pos, map, xout = p, rule = 2)$y * EXPANSION)
}, by = chrom]
grid[, cM := cM - min(cM), by = chrom]
CHRS   <- levels(grid$chrom)
genMap <- setNames(lapply(CHRS, function(c) grid[chrom == c]$cM / 100), CHRS)
NLOC   <- sapply(genMap, length)
msg(sprintf("grid: %s loci (%s), map %.0f cM per meiosis",
            format(sum(NLOC), big.mark = ","),
            paste(NLOC, collapse = "/"), sum(sapply(genMap, max)) * 100))

## locus index of the simulated QTL: the sid-2 region on chromosome III
QTL_CHR <- which(CHRS == "III")
QTL_POS <- 13680248L
QTL_IDX <- which.min(abs(grid[chrom == "III"]$pos - QTL_POS))
msg(sprintf("QTL placed at III:%s (grid locus %d of %d)",
            format(grid[chrom == "III"]$pos[QTL_IDX], big.mark = ","),
            QTL_IDX, NLOC[QTL_CHR]))

## --- one experiment ---------------------------------------------------------
## Returns per-locus true frequencies in both arms. Selection acts only in the
## RNAi arm and only over the two post-split generations, which is where it acts
## in the real experiment: the ten AIL generations are shared.
simulate_one <- function(seed, s = 0) {
  set.seed(seed)
  haplo <- setNames(lapply(NLOC, function(k)
    rbind(rep(0L, k), rep(0L, k), rep(1L, k), rep(1L, k))), CHRS)
  mp  <- newMapPop(genMap = genMap, haplotypes = haplo)
  SPx <- SimParam$new(mp); SPx$setTrackRec(FALSE)
  pop <- randCross(newPop(mp, simParam = SPx), nCrosses = 500, nProgeny = 20,
                   simParam = SPx)
  for (g in seq_len(AIL_GEN))
    pop <- randCross(pop, nCrosses = N_AIL, nProgeny = 1, simParam = SPx)

  ## one arm: POST_GEN generations, the last expanded to N_SEQ, optional
  ## fitness-weighted parent sampling on the QTL genotype
  one_arm <- function(p, sel) {
    for (g in seq_len(POST_GEN)) {
      n_out <- if (g == POST_GEN) N_SEQ else N_AIL
      if (sel > 0) {
        g_qtl <- pullSegSiteGeno(p, chr = QTL_CHR, simParam = SPx)[, QTL_IDX]
        w     <- 1 - sel * (1 - g_qtl / 2)        # additive, resistant = allele 1
        par   <- matrix(sample.int(p@nInd, 2L * n_out, replace = TRUE,
                                   prob = w), ncol = 2L)
        p <- makeCross(p, par, simParam = SPx)
      } else {
        p <- randCross(p, nCrosses = n_out, nProgeny = 1, simParam = SPx)
      }
    }
    unlist(lapply(seq_along(CHRS), function(i)
      colMeans(pullSegSiteGeno(p, chr = i, simParam = SPx)) / 2), use.names = FALSE)
  }
  list(rnai = one_arm(pop, s), ctrl = one_arm(pop, 0))
}

## --- the pipeline's estimator, applied to simulated truth -------------------
## Reads are drawn per 10 kb bin, then a moving mean of KERNEL_KB/10 bins stands
## in for calcAFD's kernel, then the standard error is formed the way the
## pipeline forms it: min(declared individuals, reads inside the kernel).
KW <- max(1L, round(KERNEL_KB * 1e3 / GRID))
chrom_vec <- rep(CHRS, NLOC)
smooth_by_chrom <- function(x) {
  unlist(lapply(CHRS, function(c) {
    v <- x[chrom_vec == c]
    frollmean(v, KW, align = "center", na.rm = TRUE, fill = NA)
  }), use.names = FALSE)
}
observe <- function(truth, reads_per_bin) {
  k <- rbinom(length(truth), reads_per_bin, truth)
  smooth_by_chrom(k / reads_per_bin)
}
lod_package <- function(z) {
  lp <- log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
  l  <- suppressWarnings(qchisq(log(2) + lp, df = 1, lower.tail = FALSE,
                                log.p = TRUE)) / (2 * log(10))
  l[is.nan(l)] <- 0; l
}
scan_of <- function(sim) {
  f_r <- observe(sim$rnai, READS_BIN[["rnai"]])
  f_c <- observe(sim$ctrl, READS_BIN[["ctrl"]])
  n_r <- min(NINDV, READS_BIN[["rnai"]] * KW)
  n_c <- min(NINDV, READS_BIN[["ctrl"]] * KW)
  se  <- sqrt(0.25 / n_r + 0.25 / n_c)
  z   <- (f_r - f_c) / se
  data.table(chrom = chrom_vec, pos = unlist(lapply(CHRS, function(c) grid[chrom == c]$pos)),
             beta = f_r - f_c, z = z, LOD = lod_package(z))[is.finite(LOD)]
}

## --- interval rules to be compared ------------------------------------------
## Each takes the scan of one chromosome and returns the interval around its
## peak. frac rules are the shipped ones; abs rules are the classical LOD drop;
## thr is the contiguous run clearing the genome-wide threshold.
interval <- function(d, rule, value, threshold) {
  i <- which.max(d$LOD)
  cut <- switch(rule,
                frac = d$LOD[i] * (1 - value),
                abs  = d$LOD[i] - value,
                thr  = threshold)
  lo <- i; while (lo > 1L      && d$LOD[lo - 1L] >= cut) lo <- lo - 1L
  hi <- i; while (hi < nrow(d) && d$LOD[hi + 1L] >= cut) hi <- hi + 1L
  list(peak = d$pos[i], lo = d$pos[lo], hi = d$pos[hi], lod = d$LOD[i])
}
RULES <- rbind(
  data.table(rule = "frac", value = 0.10, label = "10% drop (shipped tables)"),
  data.table(rule = "frac", value = 0.05, label = "5% drop (Figure 2)"),
  data.table(rule = "abs",  value = 1.50, label = "1.5 LOD (classical)"),
  data.table(rule = "abs",  value = 2.00, label = "2 LOD (config LOD.drop)"),
  data.table(rule = "thr",  value = NA_real_, label = "clears genome threshold"))

## --- run ---------------------------------------------------------------------
run_null <- function(i) {
  s <- scan_of(simulate_one(1000L + i, s = 0))
  data.table(rep = i, max_LOD = max(s$LOD),
             max_LOD_chrIII = max(s[chrom == "III"]$LOD),
             sd_beta = sd(s$beta), sd_z = sd(s$z))
}
run_qtl <- function(i, s_sel, threshold) {
  s <- scan_of(simulate_one(50000L + i * 17L, s = s_sel))
  d <- s[chrom == "III"][order(pos)]
  out <- lapply(seq_len(nrow(RULES)), function(k) {
    iv <- interval(d, RULES$rule[k], RULES$value[k], threshold)
    data.table(rep = i, s = s_sel, label = RULES$label[k],
               peak = iv$peak, lo = iv$lo, hi = iv$hi, lod = iv$lod,
               width_kb = (iv$hi - iv$lo) / 1e3,
               displacement_kb = abs(iv$peak - QTL_POS) / 1e3,
               covers = QTL_POS >= iv$lo & QTL_POS <= iv$hi,
               genome_max = max(s$LOD))
  })
  rbindlist(out)
}

if (!REFIT && file.exists(CACHE)) {
  R <- readRDS(CACHE); msg("loaded cache (CROSS_SIM_REFIT=1 to rerun)")
} else {
  msg(sprintf("null: %d replicates on %d cores", NREP_NULL, CORES))
  NUL <- rbindlist(parallel::mclapply(seq_len(NREP_NULL), run_null, mc.cores = CORES))
  thr_sim <- as.numeric(quantile(NUL$max_LOD, 0.95))
  msg(sprintf("simulated genome-wide 5%% threshold: LOD %.2f (shipped %.2f)",
              thr_sim, THR_SHIP))
  msg(sprintf("qtl: %d replicates x %d effects", NREP_QTL, length(S_SWEEP)))
  QTL <- rbindlist(lapply(S_SWEEP, function(ss)
    rbindlist(parallel::mclapply(seq_len(NREP_QTL), run_qtl, s_sel = ss,
                                 threshold = thr_sim, mc.cores = CORES))))
  R <- list(null = NUL, qtl = QTL, thr_sim = thr_sim)
  saveRDS(R, CACHE)
}
NUL <- R$null; QTL <- R$qtl; thr_sim <- R$thr_sim

## --- what it says ------------------------------------------------------------
cat("\n== the genome-wide threshold ==\n")
cat(sprintf("  shipped (effective.n.tests = 2000)      : LOD %.2f\n", THR_SHIP))
cat(sprintf("  simulated 95th pct of genome-wide max   : LOD %.2f\n", thr_sim))
cat(sprintf("  null replicates exceeding the shipped   : %d of %d (%.0f%%)\n",
            sum(NUL$max_LOD > THR_SHIP), nrow(NUL),
            100 * mean(NUL$max_LOD > THR_SHIP)))
cat(sprintf("  null sd(z) under the pipeline's own SE  : %.2f (should be 1)\n",
            median(NUL$sd_z)))
fwrite(NUL, file.path(DIAG, "cross_null_simulation.tsv"), sep = "\t")

cat("\n== interval rules, coverage of the true QTL position ==\n")
COV <- QTL[, .(n = .N, coverage = mean(covers),
               median_width_kb = median(width_kb),
               median_displacement_kb = median(displacement_kb),
               median_peak_LOD = median(lod)), by = .(s, label)]
setorder(COV, s, -coverage)
print(COV)
fwrite(COV, file.path(DIAG, "cross_interval_calibration.tsv"), sep = "\t")

cat("\n== displacement of the called peak from the true QTL ==\n")
print(QTL[label == "10% drop (shipped tables)",
          .(median_kb = round(median(displacement_kb), 1),
            p90_kb = round(quantile(displacement_kb, 0.9), 1),
            max_kb = round(max(displacement_kb), 1)), by = s])

## --- figure ------------------------------------------------------------------
theme_set(theme_bw(base_size = 9))
pA <- ggplot(NUL, aes(max_LOD)) +
  geom_histogram(bins = 30, fill = "grey75", colour = "grey30", linewidth = 0.2) +
  geom_vline(xintercept = THR_SHIP, colour = "firebrick", linewidth = 0.5) +
  geom_vline(xintercept = thr_sim, colour = "steelblue4", linewidth = 0.5,
             linetype = "dashed") +
  labs(title = "A  genome-wide maximum LOD under the simulated null",
       subtitle = sprintf("red: shipped threshold %.2f   blue: simulated 95%% %.2f",
                          THR_SHIP, thr_sim),
       x = "genome-wide max LOD (no QTL anywhere)", y = "null replicates")

pB <- ggplot(COV, aes(factor(s), coverage, fill = label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.4) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "B  does the interval contain the QTL?",
       subtitle = "dashed line is the 95% a confidence interval would have to reach",
       x = "selection coefficient against the susceptible homozygote",
       y = "coverage", fill = NULL) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6)) +
  guides(fill = guide_legend(nrow = 2))

pC <- ggplot(QTL, aes(factor(s), width_kb, colour = label)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.8, seed = 1),
             size = 0.5, alpha = 0.5) +
  scale_y_log10() +
  labs(title = "C  how wide each rule makes the interval",
       x = "selection coefficient", y = "interval width (kb, log)", colour = NULL) +
  theme(legend.position = "none")

pD <- ggplot(QTL[label == "10% drop (shipped tables)"],
             aes(factor(s), displacement_kb)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3, fill = "grey90") +
  geom_hline(yintercept = KERNEL_KB, colour = "firebrick", linetype = "dotted") +
  labs(title = "D  how far the called peak sits from the true QTL",
       subtitle = sprintf("dotted: the ~%d kb smoothing kernel", KERNEL_KB),
       x = "selection coefficient", y = "|called peak - true QTL| (kb)")

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "A simulated null for the AIL cross scans",
    subtitle = sprintf(paste("10 generations of intercross at N = %s, split at gen 11,",
                             "%d generations per arm, %s sequenced;\n%d null and %d x %d",
                             "QTL replicates"),
                       format(N_AIL, big.mark = ","), POST_GEN,
                       format(N_SEQ, big.mark = ","), nrow(NUL),
                       length(S_SWEEP), NREP_QTL))
ggsave(file.path(DIAG, "DIAG_cross_null_simulation.pdf"), fig,
       width = 11, height = 8.5, device = cairo_pdf)
ggsave(file.path(DIAG, "DIAG_cross_null_simulation.png"), fig,
       width = 11, height = 8.5, dpi = 200)
msg("wrote DIAG_cross_null_simulation.{pdf,png}")
