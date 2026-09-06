## Leakage in the MIP-seq validation -----------------------------------------
##
##   Rscript scripts/baugh_leakage_vs_similarity.R
##     -> plots/diagnostics/baugh_leakage_vs_similarity.{pdf,png}
##     -> plots/diagnostics/TABLE_baugh_leakage.tsv
##
##   A  nearest-neighbour IBS against NNLS-MIP discrepancy, the direct port of
##      panel D of the dilution figure
##   B  the same, restricted to strains whose IBS nearest neighbour is actually
##      in this panel -- the only strains for which the predictor is well posed
##   C  per-strain agreement against IBS, the control panel
##   D  confusable pairs trade signal: residual correlation within pairs
##
## WHY THIS EXISTS. Panel D of SUPP_FIG_XX_dilution_validation asks whether the
## strains the deconvolution resolves badly are the genetically similar ones,
## and answers it with leakage into pools a strain is not in. That measure needs
## a designed zero. The Baugh L1 experiment has no designed zero -- every strain
## is in the pool -- but it has something the dilution experiment does not: an
## independent measurement of the same material, published MIP-seq. This script
## runs the same question against that second measurement.
##
## ---------------------------------------------------------------------------
## THE THREE THINGS THAT MAKE THIS WEAKER THAN THE DILUTION PANEL
##
## 1. MIP-SEQ IS NOT GROUND TRUTH. It is a second measurement with its own
##    error. A correlation between relatedness and NNLS-MIP disagreement shows
##    that the two platforms diverge on related strains; it does not by itself
##    say which is wrong. The asymmetry that makes NNLS the likely source is an
##    argument, not a measurement: MIP-seq genotypes strain-diagnostic sites
##    directly, so it has no confusability term, while NNLS solves a linear
##    system whose columns are near-collinear for related strains. Panel D is
##    the part of this figure that does not rest on that argument.
##
## 2. THE IBS TABLE IS FROM THE WRONG REFERENCE. nn_ibs in the deposit is
##    identity-by-state to the closest other strain among the 170 DILUTION
##    strains. NNLS confusability depends on the nearest confounder IN THE
##    REFERENCE ACTUALLY USED, and the Baugh reference is a different panel.
##    Only 66 of these 98 strains appear in that table at all, and for only 24
##    of them is the listed nearest neighbour also in this panel -- for the
##    other 42 the true nearest confounder is someone else and nn_ibs is an
##    overestimate of the distance to it. Panel B is the well-posed subset.
##    Recomputing IBS within the Baugh reference needs
##    data/genotypes/processed_genotype_matrix.Rda, which is Dryad-hosted; see
##    make_experiments_deposit.R for the recipe. That is the single change that
##    would most improve this analysis.
##
## 3. THE EFFECT IS IN ABSOLUTE UNITS, NOT FOLD UNITS. The correlation is with
##    the absolute discrepancy and does not survive on a log ratio (rho +0.18,
##    p 0.15). Leakage moves a QUANTITY of frequency mass between two columns,
##    and how large that quantity is depends on how abundant the pair is, not
##    on either strain's own share. So the honest statement is that relatedness
##    predicts how much frequency mass goes astray, not that related strains are
##    proportionally worse measured.
##
## 4. IT IS A THRESHOLD, NOT A GRADIENT -- and this is the main result, not a
##    caveat about it. Reporting a Spearman rho invites the reading that the
##    error climbs steadily with relatedness. It does not. Drop the six strains
##    above IBS 0.97 and the correlation falls to +0.26 (p 0.04) in panel A and
##    to +0.23 (p 0.33) in panel B. What is actually there is a step: those six
##    strains have a median gap of 10.8 per mille against 2.3 for the other 60,
##    a 4.6-fold difference at p 0.003. Below 0.97 IBS the deconvolution is not
##    measurably degraded by relatedness at this sample size; above it, it can
##    move a large share of one strain's signal onto its neighbour. Four of the
##    six are one clade -- NIC256, NIC271, JU782, NIC262 -- so the evidence is
##    thinner than n = 6 suggests, and a second such clade would be worth more
##    than any amount of extra sequencing on the rest of the panel.
##
## ---------------------------------------------------------------------------
## SCALES. The two platforms are not on a common scale: MIP-seq frequencies sum
## to about 1.6 per sample because N2 alone is assigned 0.573, against 0.032 by
## NNLS. N2 is dropped -- Figure1_common.R drops it everywhere for the same
## reason, that its genotype shares sites with every strain in the pool -- as
## are the three strains with no MIP-seq column (CX11262, ECA348, NIC260),
## leaving the 98 strains of Figure 1A. Both platforms are then renormalised to
## sum to one within each sample. The renormalisation is cosmetic: on raw
## frequencies the headline correlation is +0.414 rather than +0.391.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
})

## ggtext gives the markdown panel letters the other figure scripts use, but it
## is the only hard dependency here that a plain tidyverse install may lack, and
## nothing in this analysis needs rich text. Fall back to a bold plain title so
## the script runs wherever Figure1_common.R does and also where it does not.
HAS_GGTEXT <- requireNamespace("ggtext", quietly = TRUE)

DEC <- "supplemental_data/deconvolution"
OUT <- "plots/diagnostics"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

NNLS <- file.path(DEC, "baugh_nnls_with_mipseq.RData")
BOOT <- file.path(DEC, "baugh_bootstrap_array.rda")
ORD  <- file.path(DEC, "baugh_strain_order.txt")
CACHE<- file.path(DEC, "cache_boot_freq.rds")

## Prefer IBS computed inside the Baugh reference; fall back to the dilution
## table, which is a borrowed predictor. Which one is in use changes what this
## script can claim, so it is resolved once, here, and reported everywhere.
SIM_TRUE  <- file.path(DEC, "baugh_strain_similarity.tsv")
SIM_PROXY <- file.path(DEC, "dilution_strain_similarity.tsv")
TRUE_REF  <- file.exists(SIM_TRUE)
SIM       <- if (TRUE_REF) SIM_TRUE else SIM_PROXY
stopifnot(file.exists(NNLS), file.exists(BOOT), file.exists(ORD),
          file.exists(SIM), file.exists(CACHE))

panel_title <- function(letter) {
  if (HAS_GGTEXT)
    paste0("<span style='font-size:13pt;color:#111111'>**", letter, "**</span>")
  else letter
}

theme_pub <- function(base_size = 11) {
  title_el <- if (HAS_GGTEXT) ggtext::element_markdown(size = base_size + 0.5)
              else element_text(size = base_size + 2, face = "bold")
  theme_classic(base_size = base_size) +
    theme(axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = base_size - 1),
          plot.title = title_el,
          plot.title.position = "plot",
          legend.key.size = grid::unit(9, "pt"))
}

## Spearman with the estimate and p in one string, and a partial version that
## removes a third variable by ranking all three first. Abundance has to be
## removable because an ABSOLUTE discrepancy rises with abundance on its own --
## the console prints that correlation -- and relatedness could ride on it
## rather than cause anything. It does not: the partial is barely attenuated.
sp <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  list(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}
pspear <- function(x, y, z) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z)
  x <- rank(x[ok]); y <- rank(y[ok]); z <- rank(z[ok]); n <- length(x)
  rp <- (cor(x, y) - cor(x, z) * cor(y, z)) /
        sqrt((1 - cor(x, z)^2) * (1 - cor(y, z)^2))
  tt <- rp * sqrt((n - 3) / (1 - rp^2))
  list(rho = rp, p = 2 * pt(-abs(tt), n - 3), n = n)
}
fmt <- function(s) sprintf("rho %+.3f, p %.2g, n %d", s$rho, s$p, s$n)
lab <- function(s) sprintf("rho = %+.2f, p = %.1g", s$rho, s$p)

## ---------------------------------------------------------------------------
## THE PINNED NUMBERS
##
## As in the dilution script: the numbers this analysis reports are written out
## as literals and checked against what the deposit computes, so the figure and
## anything quoting it cannot drift apart silently. Update both together.
## ---------------------------------------------------------------------------
## These pin the PROXY run -- the one this script was written against. On the
## true-reference table every one of them is expected to change, so they are
## checked only when the proxy is in use. Replace them with the new values once
## baugh_strain_similarity.R has been run and the output reviewed.
PINNED <- list(n_strains = 98L, n_ibs = 66L, n_valid = 24L,
               rho_all = 0.391, rho_valid = 0.536, rho_ctrl_valid = -0.248,
               ibs_step = 0.97, n_above = 6L)

## ===========================================================================
## the data
## ===========================================================================
e <- new.env(); load(NNLS, e)
d0 <- as_tibble(e$wgs_mip_results)

no_mip <- d0 %>% group_by(strain) %>%
  summarise(k = sum(is.na(published_frq)), .groups = "drop") %>%
  filter(k > 0) %>% pull(strain)

freq <- d0 %>%
  filter(strain != "N2", !strain %in% no_mip) %>%
  group_by(sample) %>%
  mutate(wgs = frq / sum(frq), mip = published_frq / sum(published_frq)) %>%
  ungroup()
stopifnot(length(unique(freq$strain)) == PINNED$n_strains)

## per-strain summaries. d_abs is the leakage analogue; rho_s -- how well NNLS
## tracks MIP-seq across the 23 samples -- is the control, standing in for the
## dilution panel's own-pool recovery.
perstrain <- freq %>% group_by(strain) %>%
  summarise(mip_mean = mean(mip),
            wgs_mean = mean(wgs),
            d_signed = mean(wgs - mip),
            d_abs    = mean(abs(wgs - mip)),
            alr      = median(abs(log2((wgs + 1e-6) / (mip + 1e-6)))),
            ## JU2001 is flat on one platform across all 23 samples, so its
            ## rank correlation is undefined; it drops out of the control only
            ## and is why panel C reports n = 65 against panel A's n = 66.
            rho_s    = suppressWarnings(cor(wgs, mip, method = "spearman")),
            .groups  = "drop")

## ---------------------------------------------------------------------------
## bootstrap confusability, which needs no genotypes
##
## ab is 102 strains x 23 samples x 100 bootstrap replicates. Two strains whose
## columns the solver cannot separate must trade signal between them as the
## input is perturbed, so their estimates anti-correlate ACROSS REPLICATES.
## Correlating within a sample and averaging over samples on the Fisher z scale
## gives a confusability matrix in the RIGHT reference -- the Baugh one -- which
## is what nn_ibs is not.
##
## CAVEAT. METHODS.txt records that the resampling code for this bootstrap was
## never archived, so what was resampled is not known. That does not invalidate
## the reading -- whatever the perturbation, two strains that trade signal under
## it are trading signal -- but it does mean the matrix cannot be calibrated,
## which is why it is used to RANK pairs in panel D and never as an effect size.
## ---------------------------------------------------------------------------
b <- new.env(); load(BOOT, b)
ord <- readLines(ORD)[-1]
stopifnot(dim(b$ab)[1] == length(ord))

## the array carries no strain names, so the row order is asserted rather than
## assumed: bootse recomputed against the cached per-strain SDs must match
## exactly before anything is read out of ab by name
bse <- b$bootse; rownames(bse) <- ord
chk <- as_tibble(readRDS(CACHE))
stopifnot(max(abs(bse[cbind(chk$strain, chk$sample)] - chk$boot_sd)) == 0)

ab <- b$ab; dimnames(ab)[[1]] <- ord
keep <- ord %in% perstrain$strain
zsum <- 0
for (s in dimnames(ab)[[2]]) {
  cc <- suppressWarnings(cor(t(ab[keep, s, ])))
  cc[!is.finite(cc)] <- 0
  zsum <- zsum + atanh(pmin(pmax(cc, -0.999), 0.999))
}
conf_mat <- tanh(zsum / dim(ab)[2])
diag(conf_mat) <- NA
rownames(conf_mat) <- colnames(conf_mat) <- ord[keep]

similarity <- read_tsv(SIM, show_col_types = FALSE)
tab <- perstrain %>%
  left_join(tibble(strain   = rownames(conf_mat),
                   min_r    = apply(conf_mat, 1, min, na.rm = TRUE),
                   boot_partner = colnames(conf_mat)[apply(conf_mat, 1, which.min)]),
            by = "strain") %>%
  left_join(similarity, by = "strain") %>%
  mutate(well_posed = !is.na(nn_ibs) & nn_partner %in% perstrain$strain)

if (TRUE_REF) {
  cat("== IBS source: ", SIM_TRUE, " (the Baugh reference) ==\n", sep = "")
  ## every strain's nearest neighbour is by construction in this panel, so the
  ## proxy/well-posed split collapses and panels A and B become the same test
  stopifnot(all(tab$well_posed))
} else {
  cat("== IBS source: ", SIM_PROXY, " ==\n", sep = "")
  cat("   BORROWED PREDICTOR. Distances are to the closest of the 170 DILUTION\n")
  cat("   strains, not of this panel. Run scripts/baugh_strain_similarity.R\n")
  cat("   against the archive to replace it; see DATA_AVAILABILITY.md.\n")
  stopifnot(sum(!is.na(tab$nn_ibs)) == PINNED$n_ibs,
            sum(tab$well_posed)     == PINNED$n_valid)
}

has_ibs <- tab %>% filter(!is.na(nn_ibs))
valid   <- tab %>% filter(well_posed)

## ===========================================================================
## the tests
## ===========================================================================
s_all    <- with(has_ibs, sp(nn_ibs, d_abs))
s_allp   <- with(has_ibs, pspear(nn_ibs, d_abs, mip_mean))
s_valid  <- with(valid,   sp(nn_ibs, d_abs))
s_validp <- with(valid,   pspear(nn_ibs, d_abs, mip_mean))
s_ctrl   <- with(has_ibs, sp(nn_ibs, rho_s))
s_ctrlv  <- with(valid,   sp(nn_ibs, rho_s))
s_lr     <- with(has_ibs, sp(nn_ibs, alr))
s_abund  <- with(has_ibs, sp(mip_mean, d_abs))

if (!TRUE_REF)
  stopifnot(abs(s_all$rho   - PINNED$rho_all)        < 5e-3,
            abs(s_valid$rho - PINNED$rho_valid)      < 5e-3,
            abs(s_ctrlv$rho - PINNED$rho_ctrl_valid) < 5e-3)

cat("== the port of dilution panel D to the MIP-seq validation ==\n")
cat("  strains: ", PINNED$n_strains, " (98 of Figure 1A; N2 and 3 MIP-less dropped)\n", sep = "")
cat("  with an nn_ibs value          : ", nrow(has_ibs), "\n", sep = "")
cat("  with a WELL-POSED nn_ibs value: ", nrow(valid),
    "  <- nearest neighbour also in this panel\n\n", sep = "")

cat("  A  nn_ibs vs |NNLS-MIP|, all with IBS   : ", fmt(s_all),   "\n", sep = "")
cat("       abundance held constant            : ", fmt(s_allp),  "\n", sep = "")
cat("  B  nn_ibs vs |NNLS-MIP|, well-posed     : ", fmt(s_valid), "\n", sep = "")
cat("       abundance held constant            : ", fmt(s_validp),"\n", sep = "")
cat("  C  nn_ibs vs agreement,  all with IBS   : ", fmt(s_ctrl),  "  [control]\n", sep = "")
cat("     nn_ibs vs agreement,  well-posed     : ", fmt(s_ctrlv), "  [control]\n", sep = "")
cat("     nn_ibs vs |log2 ratio|               : ", fmt(s_lr),
    "  [scale-free: null]\n", sep = "")
cat("     abundance vs |NNLS-MIP|              : ", fmt(s_abund),
    "  [the confound]\n\n", sep = "")

## the IBS bins used by both the console table and the panel medians
BRK <- c(0.66, 0.90, 0.95, 0.97, 0.99, 1.0)

## the binned medians the panels draw, with their counts. The dilution script
## prints the same table for the same reason: the rise in panels A and B is
## carried by the top bin, and the top bin is small.
bin_report <- function(df, what) {
  b <- df %>% mutate(bin = cut(nn_ibs, BRK, include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarise(n = n(), median_gap_permille = signif(1000 * median(d_abs), 3),
              .groups = "drop") %>% filter(!is.na(bin))
  cat("  binned medians, ", what, ":\n", sep = "")
  print(as.data.frame(b), row.names = FALSE)
  cat("    NOTE: the top bin holds ", b$n[nrow(b)],
      " strains, so read the trend, not that bin.\n\n", sep = "")
}

## jackknife, because n = 24 invites one strain to carry the result
bin_report(has_ibs, "panel A")
bin_report(valid,   "panel B")

jk <- sapply(seq_len(nrow(valid)), function(i)
  cor(valid$nn_ibs[-i], valid$d_abs[-i], method = "spearman"))
cat(sprintf("  jackknife on the well-posed 24: rho %+.3f to %+.3f, no strain carries it\n",
            min(jk), max(jk)))
cat("    (most influential single strain: ", valid$strain[which.min(jk)], ")\n\n", sep = "")

## ---------------------------------------------------------------------------
## graded or a step? Drop everything above the threshold and refit. If the
## correlation survives, relatedness degrades the deconvolution continuously;
## if it collapses, the whole effect is a small number of very close pairs and
## should be reported that way.
## ---------------------------------------------------------------------------
STEP <- PINNED$ibs_step
below_a <- has_ibs %>% filter(nn_ibs <= STEP)
below_b <- valid   %>% filter(nn_ibs <= STEP)
above   <- has_ibs %>% filter(nn_ibs >  STEP)
if (!TRUE_REF) stopifnot(nrow(above) == PINNED$n_above)

s_below_a <- with(below_a, sp(nn_ibs, d_abs))
s_below_b <- with(below_b, sp(nn_ibs, d_abs))
wstep <- suppressWarnings(wilcox.test(above$d_abs, below_a$d_abs))

cat("== graded, or a step at IBS ", STEP, "? ==\n", sep = "")
cat("  panel A with the ", nrow(above), " strains above ", STEP, " dropped: ",
    fmt(s_below_a), "\n", sep = "")
cat("  panel B with them dropped                  : ", fmt(s_below_b), "\n", sep = "")
cat(sprintf("  above %.2f: n = %d, median gap %.2f per mille\n",
            STEP, nrow(above), 1000 * median(above$d_abs)))
cat(sprintf("  at or below : n = %d, median gap %.2f per mille\n",
            nrow(below_a), 1000 * median(below_a$d_abs)))
cat(sprintf("  Mann-Whitney p = %.2g, a %.1f-fold difference in median\n\n",
            wstep$p.value, median(above$d_abs) / median(below_a$d_abs)))
cat("  the strains above the threshold:\n")
print(as.data.frame(above %>% arrange(desc(d_abs)) %>%
        transmute(strain, nn_ibs = round(nn_ibs, 4), nn_partner, well_posed,
                  gap_permille = round(1000 * d_abs, 2),
                  agreement = round(rho_s, 2))), row.names = FALSE)
cat("  Four are one clade (NIC256, NIC271, JU782, NIC262), so this rests on\n",
    "  fewer independent observations than n = ", nrow(above), " suggests.\n\n", sep = "")

## ---------------------------------------------------------------------------
## panel D's test: do confusable pairs trade signal?
##
## This is the one test here that does not need MIP-seq to be right. If NNLS
## moves mass from one strain to a near-collinear neighbour, then the two
## strains' departures from MIP-seq must move in OPPOSITE directions sample by
## sample, whatever the true frequencies were. So: rank pairs by bootstrap
## confusability, and compare the residual correlation of the top 2% against
## every other pair.
## ---------------------------------------------------------------------------
resid_mat <- freq %>% mutate(diff = wgs - mip) %>%
  select(strain, sample, diff) %>%
  pivot_wider(names_from = strain, values_from = diff) %>%
  arrange(sample)
R <- cor(as.matrix(resid_mat[, -1]), method = "spearman")
diag(R) <- NA
common <- intersect(colnames(R), colnames(conf_mat))
R <- R[common, common]; Z <- conf_mat[common, common]

up  <- upper.tri(R)
rv  <- R[up]; zv <- Z[up]
cut <- quantile(zv, 0.02, na.rm = TRUE)
top <- zv <= cut
wt  <- suppressWarnings(wilcox.test(rv[top], rv[!top]))

cat("== do confusable pairs trade signal? ==\n")
cat(sprintf("  residual correlation, most-confusable 2%% of pairs: median %+.3f (n = %d)\n",
            median(rv[top], na.rm = TRUE), sum(top, na.rm = TRUE)))
cat(sprintf("  residual correlation, every other pair           : median %+.3f (n = %d)\n",
            median(rv[!top], na.rm = TRUE), sum(!top, na.rm = TRUE)))
cat(sprintf("  Wilcoxon p = %.2g\n\n", wt$p.value))

pairs_tbl <- {
  ix <- which(Z <= cut & up, arr.ind = TRUE)
  tibble(a = rownames(Z)[ix[, 1]], b = colnames(Z)[ix[, 2]],
         boot_r = Z[ix], resid_r = R[ix]) %>% arrange(boot_r)
}
cat("  the ten most confusable pairs:\n")
print(as.data.frame(pairs_tbl %>% head(10) %>%
        mutate(across(where(is.numeric), ~round(.x, 3)))), row.names = FALSE)

nic <- freq %>% filter(strain %in% c("NIC256", "NIC271")) %>% group_by(strain) %>%
  summarise(mip = mean(mip), nnls = mean(wgs), .groups = "drop") %>%
  mutate(ratio = nnls / mip)
cat(sprintf("\n  NIC256/NIC271, IBS %.4f: NNLS assigns %.2fx MIP-seq to NIC256 and\n",
            tab$nn_ibs[tab$strain == "NIC256"], nic$ratio[nic$strain == "NIC256"]))
cat(sprintf("  %.2fx to NIC271, and their residuals correlate %+.2f. Mass moved from\n",
            nic$ratio[nic$strain == "NIC271"],
            R["NIC256", "NIC271"]))
cat("  one to the other; the pair's total is close to right.\n\n")

## ===========================================================================
## panels
## ===========================================================================
binmed <- function(df) df %>%
  mutate(bin = cut(nn_ibs, BRK, include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(med = median(d_abs), mid = median(nn_ibs), .groups = "drop") %>%
  filter(!is.na(bin))

## The six strains above the threshold are named on the panels, because the
## rise in the binned median is theirs and a reader should be able to see that
## it is six points and which six, not infer a gradient from a rho.
ibs_panel <- function(df, letter, subtitle, stat) {
  bm <- binmed(df)
  ggplot(df, aes(nn_ibs, 1000 * d_abs)) +
    geom_vline(xintercept = STEP, linetype = "dotted", linewidth = 0.35,
               colour = "grey55") +
    geom_point(aes(fill = nn_ibs > STEP), shape = 21, size = 1.9, stroke = 0.25,
               colour = "grey30", alpha = 0.9, show.legend = FALSE) +
    scale_fill_manual(values = c(`TRUE` = "#B5446E", `FALSE` = "#3F7CAC")) +
    geom_text_repel(data = df %>% filter(nn_ibs > STEP), aes(label = strain),
                    size = 2.4, colour = "grey20", seed = 1,
                    min.segment.length = 0, segment.size = 0.25,
                    segment.colour = "grey60", box.padding = 0.4) +
    geom_line(data = bm, aes(mid, 1000 * med), inherit.aes = FALSE,
              colour = "grey25", linewidth = 0.5) +
    geom_point(data = bm, aes(mid, 1000 * med), inherit.aes = FALSE,
               colour = "grey15", size = 1.6, shape = 15) +
    annotate("text", x = min(df$nn_ibs), y = Inf, label = lab(stat),
             hjust = 0, vjust = 1.6, size = 2.7, colour = "grey20") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
    labs(x = "Identity-by-state to the closest other strain",
         y = "Absolute NNLS − MIP-seq gap (per mille)",
         title = panel_title(letter), subtitle = subtitle) +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 8, colour = "grey35"))
}

pA <- ibs_panel(has_ibs, "A",
                if (TRUE_REF)
                  sprintf("all %d strains, IBS within the Baugh reference",
                          nrow(has_ibs))
                else
                  sprintf("all %d strains with a (borrowed) IBS value",
                          nrow(has_ibs)),
                s_all)
pB <- ibs_panel(valid, "B",
                if (TRUE_REF)
                  "same data as A: with the true reference the split is empty"
                else
                  sprintf("the %d whose nearest neighbour is in this panel",
                          nrow(valid)),
                s_valid)

pC <- ggplot(has_ibs, aes(nn_ibs, rho_s)) +
  geom_point(aes(fill = well_posed), shape = 21, size = 1.9, stroke = 0.25,
             colour = "grey30", alpha = 0.9) +
  scale_fill_manual(values = c(`TRUE` = "#B5446E", `FALSE` = "grey75"),
                    labels = c(`TRUE` = "well posed", `FALSE` = "proxy IBS"),
                    name = NULL) +
  annotate("text", x = min(has_ibs$nn_ibs), y = -Inf,
           label = sprintf("all: %s\nwell posed: %s", lab(s_ctrl), lab(s_ctrlv)),
           hjust = 0, vjust = -0.4, size = 2.5, colour = "grey20") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(x = "Identity-by-state to the closest other strain",
       y = "Agreement with MIP-seq (Spearman)",
       title = panel_title("C"),
       subtitle = "the control: does similarity also degrade agreement?") +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 8, colour = "grey35"),
        legend.position = c(0.98, 0.02), legend.justification = c(1, 0))

pD <- ggplot(tibble(grp = ifelse(top, "most confusable 2%", "all other pairs"),
                    r = rv) %>% filter(is.finite(r)),
             aes(grp, r, fill = grp)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3,
             colour = "grey60") +
  geom_violin(colour = NA, alpha = 0.45, width = 0.85) +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.35,
               fill = "white") +
  scale_fill_manual(values = c(`most confusable 2%` = "#B5446E",
                               `all other pairs` = "grey70"), guide = "none") +
  annotate("text", x = 1.5, y = Inf,
           label = sprintf("Wilcoxon p = %.1g", wt$p.value),
           vjust = 1.6, size = 2.7, colour = "grey20") +
  labs(x = NULL, y = "Residual correlation within the pair",
       title = panel_title("D"),
       subtitle = "confusable pairs trade signal: one gains what the other loses") +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 8, colour = "grey35"))

fig <- (pA | pB) / (pC | pD)

## The other figure scripts write PDFs through cairo_pdf. Cairo is not always
## built into a mac R -- it needs XQuartz -- and ggsave reports success even
## when the device fails to load, leaving no file behind. Fall back to the
## default pdf device and say which one was used.
pdf_path <- file.path(OUT, "baugh_leakage_vs_similarity.pdf")
dev <- if (isTRUE(capabilities("cairo"))) cairo_pdf else pdf
ggsave(pdf_path, fig, width = 8.6, height = 7.4, device = dev)
if (!file.exists(pdf_path))
  ggsave(pdf_path, fig, width = 8.6, height = 7.4, device = pdf)
stopifnot(file.exists(pdf_path))
ggsave(file.path(OUT, "baugh_leakage_vs_similarity.png"), fig,
       width = 8.6, height = 7.4, dpi = 300, bg = "white")

write_tsv(tab %>% arrange(desc(d_abs)), file.path(OUT, "TABLE_baugh_leakage.tsv"))
cat("wrote baugh_leakage_vs_similarity.{pdf,png} and TABLE_baugh_leakage.tsv to ",
    OUT, "\n", sep = "")
