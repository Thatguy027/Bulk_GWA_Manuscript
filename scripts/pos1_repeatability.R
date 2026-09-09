## Repeatability of the pooled pos-1 response ---------------------------------
##
##   Rscript scripts/pos1_repeatability.R
##     -> plots/diagnostics/TABLE_pos1_repeatability.tsv
##
## THE QUESTION. The draft asserts a broad-sense heritability of 0.32 for the
## pooled pos-1 response. No heritability estimate exists anywhere in this
## repository, and 0.32 passes the number checker only because it rounds onto an
## unrelated leakage rho of 0.326. The 2023 experiment has four pos-1 replicate
## pools against a control, so what these data CAN support is a repeatability --
## the intraclass correlation across those replicates.
##
## IT MUST BE ON THE TRAIT THAT WAS MAPPED, which is vst_ctrl_pos-1_T2. That is
## the awkward part, because the variance-stabilising transform is not
## reproducible from this repository (METHODS.txt records this as [TO FILL]) and
## the deposit ships the VST per STRAIN, not per replicate. There is no direct
## way to compute an ICC on it.
##
## THE PROXY, AND WHY IT IS A FAIR ONE. The shipped VST is very nearly a
## monotone function of the shipped delta: Spearman 0.957 over the 231 mapped
## strains. So an isotonic regression of vst on delta, fitted on those 231
## shipped pairs, recovers the transform to R-squared 0.969 -- and applying that
## fitted map to each replicate's delta puts every replicate on the VST scale.
## That is the headline number below. It is a reconstruction rather than the
## transform itself, and if the real definition is recovered from the 2023
## analysis outside this repository, this should be recomputed and the [TO FILL]
## at METHODS.txt closed.
##
## THE ANSWER DEPENDS ON SCALE, and the spread is worth reporting:
##
##   raw delta_ctrl        ~0.95   variance still scales with abundance
##   VST (isotonic proxy)  ~0.91   the mapped trait
##   within-replicate rank ~0.85   ordering only; a lower bound on the above
##   log2fc_ctrl           ~0.44   abundance divided out entirely
##
## IS THE 0.91 JUST THE REPEATABILITY OF ABUNDANCE? No, and the control-vs-
## control null below is what rules that out. Both control replicates exist for
## the same strains, so their difference is a contrast with no response in it:
##
##   |null delta| vs abundance   Spearman +0.59  -- noise DOES scale with
##                                                  abundance, which is what a
##                                                  variance-stabilising
##                                                  transform exists to fix
##   signed null delta vs abundance      -0.01  -- but there is no systematic
##                                                  bias, only scale
##   sd(pos-1 delta) / sd(null delta)     ~18x  -- response dominates noise
##
## So the shipped VST's -0.56 correlation with abundance is NOT noise: the null
## contrast has no signed abundance trend at all. It is a real relationship in
## which rarer strains show weaker apparent response, which is what a floor on
## dynamic range looks like -- a strain at 1e-4 cannot drop far. Worth stating
## in the text as a limit on what the rare end of the panel can show, but it
## does not undermine the repeatability estimate, because the signal is roughly
## eighteen times the technical noise on the same scale.
##
## REPEATABILITY IS NOT HERITABILITY. R is an upper bound on H2, for two
## structural reasons asserted below rather than described: the four replicates
## are replicate pools within one experiment at one timepoint, spanning
## technical and within-batch variation rather than an environmental range; and
## all four share ONE control baseline, verified below to be identical across
## all four replicates for every strain, so baseline error is counted as
## among-strain signal.
##
## Depth cutoff 5 is the headline, because METHODS.txt records that cutoff 5 was
## used throughout; all three are reported.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})

FREQ   <- "supplemental_data/phenotypes/pos1_2023_sample_frequencies.csv.gz"
TRAIT  <- "supplemental_data/phenotypes/pos1_2023_association_traits.csv.gz"
OUT    <- "plots/diagnostics"
CUTOFF <- 5L
NBOOT  <- 500L
SEED   <- 100L

stopifnot(file.exists(FREQ), file.exists(TRAIT))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

f  <- fread(cmd = paste("gzcat", shQuote(FREQ)))
tr <- fread(cmd = paste("gzcat", shQuote(TRAIT)))
setnames(tr, c("delta_ctrl_pos-1_T2", "vst_ctrl_pos-1_T2",
               "negctrl_abundance_log10"), c("delta", "vst", "abund"))

## The mapped set. This is also the answer to the pool-size question: the
## draft's "224" is not this number and is not reproducible from the deposit.
mp <- tr[!is.na(vst)]
cat("== the mapped set ==\n")
cat(sprintf("  trait file rows              : %d\n", nrow(tr)))
cat(sprintf("  with a pos-1 VST measurement : %d\n", nrow(mp)))
cat(sprintf("  without                      : %d\n\n", nrow(tr) - nrow(mp)))
stopifnot(nrow(mp) == 231L)

## ---------------------------------------------------------------------------
## the shared baseline, asserted rather than assumed
## ---------------------------------------------------------------------------
p0 <- f[rnai == "pos-1" & depth_cutoff == CUTOFF & strain %in% mp$strain]
shared <- p0[, .(n = uniqueN(round(ctrl_frq, 12))), by = strain]
cat("== control baseline ==\n")
cat(sprintf("  ctrl_frq identical across all four pos-1 replicates: %d of %d strains\n\n",
            sum(shared$n == 1L), nrow(shared)))
stopifnot(all(shared$n == 1L))

## ---------------------------------------------------------------------------
## what the VST does and does not remove
## ---------------------------------------------------------------------------
ab <- p0[, .(abs_delta = mean(abs(delta_ctrl), na.rm = TRUE)), by = strain]
ab <- merge(ab, mp[, .(strain, vst, abund)], by = "strain")
rho_delta_ab <- cor(ab$abs_delta, 10^ab$abund, method = "spearman")
rho_vst_ab   <- cor(ab$vst,       ab$abund,    method = "spearman")
cat("== abundance ==\n")
cat(sprintf("  Spearman |delta_ctrl| vs control frequency : %+.3f\n", rho_delta_ab))
cat(sprintf("  Spearman shipped VST  vs log10 abundance   : %+.3f\n", rho_vst_ab))

## The control-vs-control null: a contrast with no response in it, so whatever
## abundance structure survives here is measurement, not biology.
cw <- dcast(f[rnai == "ctrl" & depth_cutoff == CUTOFF & strain %in% mp$strain,
              .(frq = sum(frq, na.rm = TRUE)), by = .(strain, replicate)],
            strain ~ replicate, value.var = "frq")
setnames(cw, 2:3, c("cA", "cB"))
cw[, `:=`(nd = cA - cB, cm = (cA + cB) / 2)]
cw <- cw[cm > 0]
sig <- p0[, .(d = sum(delta_ctrl, na.rm = TRUE)), by = strain]
snr <- sd(sig$d) / sd(cw$nd)
rho_null_abs    <- cor(abs(cw$nd), cw$cm,        method = "spearman")
rho_null_signed <- cor(cw$nd,      log10(cw$cm), method = "spearman")
cat(sprintf("  NULL |ctrlA - ctrlB| vs abundance          : %+.3f  (noise scales)\n",
            rho_null_abs))
cat(sprintf("  NULL  ctrlA - ctrlB  vs abundance          : %+.3f  (no signed bias)\n",
            rho_null_signed))
cat(sprintf("  sd(pos-1 delta) / sd(null delta)           : %.1fx\n\n", snr))
stopifnot(snr > 5)

## ---------------------------------------------------------------------------
## reconstruct the transform, then apply it per replicate
## ---------------------------------------------------------------------------
o     <- mp[order(delta)]
iso   <- isoreg(o$delta, o$vst)
to_vst <- approxfun(o$delta, iso$yf, rule = 2)
fit_r2 <- cor(to_vst(mp$delta), mp$vst)^2
cat("== reconstructed delta -> VST map ==\n")
cat(sprintf("  isotonic fit on the 231 shipped pairs: R2 = %.3f, Spearman = %.4f\n\n",
            fit_r2, cor(to_vst(mp$delta), mp$vst, method = "spearman")))
stopifnot(fit_r2 > 0.95)

icc <- function(d) {
  m  <- lmer(y ~ 1 + (1 | strain), data = d, REML = TRUE)
  vc <- as.data.table(VarCorr(m))
  vg <- vc[grp == "strain", vcov]; ve <- vc[grp == "Residual", vcov]
  list(R = vg / (vg + ve), vg = vg, ve = ve, n = uniqueN(d$strain))
}

by_scale <- function(dc) {
  q   <- f[rnai == "pos-1" & depth_cutoff == dc & strain %in% mp$strain]
  raw <- q[, .(y = sum(delta_ctrl, na.rm = TRUE)), by = .(strain, replicate)]
  lfc <- q[is.finite(log2fc_ctrl), .(y = log2fc_ctrl[1]), by = .(strain, replicate)]
  lfc <- lfc[strain %in% lfc[, .N, by = strain][N == 4L, strain]]

  sets <- list(vst_isotonic    = copy(raw)[, y := to_vst(y)],
               delta_ctrl      = raw,
               rank_within_rep = copy(raw)[, y := frank(y) / .N, by = replicate],
               log2fc_ctrl     = lfc)
  rbindlist(lapply(sets, function(d) { r <- icc(d)
    data.table(depth_cutoff = dc, n_strain = r$n, n_rep = 4L,
               var_strain = r$vg, var_resid = r$ve, R = r$R) }), idcol = "scale")
}

res <- rbindlist(lapply(c(3L, 5L, 10L), by_scale))

## bootstrap over strains, on the scale being quoted
d <- f[rnai == "pos-1" & depth_cutoff == CUTOFF & strain %in% mp$strain,
       .(y = sum(delta_ctrl, na.rm = TRUE)), by = .(strain, replicate)][, y := to_vst(y)]
set.seed(SEED)
bs <- vapply(seq_len(NBOOT), function(i) {
  ss <- sample(unique(d$strain), replace = TRUE)
  dd <- rbindlist(lapply(seq_along(ss), function(j) {
    z <- d[strain == ss[j]]; z$strain <- paste0("s", j); z }))
  r <- try(icc(dd)$R, silent = TRUE)
  if (inherits(r, "try-error")) NA_real_ else r
}, numeric(1))
ci <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)

cat("== repeatability by scale ==\n")
print(as.data.frame(res[order(depth_cutoff, -R)]), row.names = FALSE, digits = 4)

hd <- res[scale == "vst_isotonic" & depth_cutoff == CUTOFF]
cat(sprintf("\n  QUOTE THIS: R = %.2f on the mapped VST scale (isotonic reconstruction),\n", hd$R))
cat(sprintf("              95%% bootstrap interval [%.2f, %.2f] over %d replicates,\n",
            ci[1], ci[2], sum(!is.na(bs))))
cat(sprintf("              %d strains x 4 replicate pools, depth cutoff %d.\n", hd$n_strain, CUTOFF))
cat("              Repeatability, an upper bound on H2 -- not heritability.\n\n")
cat(sprintf("  Context: rank-only gives %.2f (a lower bound), and dividing abundance out\n",
            res[scale == "rank_within_rep" & depth_cutoff == CUTOFF, R]))
cat(sprintf("  entirely (log2 fold change) gives %.2f.\n",
            res[scale == "log2fc_ctrl" & depth_cutoff == CUTOFF, R]))
cat(sprintf("  The 0.91 is not abundance repeatability: the control-vs-control null carries\n"))
cat(sprintf("  no signed abundance trend (%+.2f) and the response is %.0fx its sd.\n",
            rho_null_signed, snr))

res[, `:=`(boot_lo = NA_real_, boot_hi = NA_real_)]
res[scale == "vst_isotonic" & depth_cutoff == CUTOFF,
    `:=`(boot_lo = ci[1], boot_hi = ci[2])]
res[, `:=`(rho_absdelta_vs_abundance = rho_delta_ab,
           rho_vst_vs_abundance      = rho_vst_ab,
           rho_null_abs_vs_abundance = rho_null_abs,
           rho_null_signed_vs_abund  = rho_null_signed,
           signal_to_null_sd         = snr,
           vst_map_r2                = fit_r2)]

fo <- file.path(OUT, "TABLE_pos1_repeatability.tsv")
fwrite(res, fo, sep = "\t")
cat("\nwrote ", fo, "\n", sep = "")
