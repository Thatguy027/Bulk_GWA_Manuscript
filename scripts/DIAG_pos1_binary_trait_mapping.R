## What binarising the pos-1 response does to the 231-strain mapping ---------
##
##   Rscript scripts/DIAG_pos1_binary_trait_mapping.R
##     -> plots/diagnostics/DIAG_pos1_binary_trait_mapping.{pdf,png}
##
## A reader asked what happens if the 81 isolates that are absent from every
## pos-1 pool are called RNAi-responsive and the trait becomes binary. Two
## things change at once -- the SCALE (continuous to binary) and the STRAIN SET
## (the 81 censored strains become usable) -- so this crosses them into a 2x2
## and scans all four traits through ONE pipeline.
##
## THE ANSWER: binarising costs the mapping, and it is the binarisation that
## does it, not the censoring. Holding the strain set at the 150 testable
## strains, continuous -> binary drops chromosome III from LOD 5.63 to 0.21.
## Holding the scale continuous, dropping the 81 censored strains costs about
## one LOD. The 81 are INFORMATIVE; collapsing magnitude is what throws the
## signal away.
##
## AND THE THREE LOCI DO NOT BEHAVE ALIKE. Chromosome IV and X survive
## binarisation at reduced power -- they separate on WHETHER a strain responds
## (responsive rate 0.29 vs 0.69 by genotype at IV, 0.20 vs 0.65 at X).
## Chromosome III, the SID-2 locus, goes to nothing (0.48 vs 0.62) while
## separating cleanly on the continuous scale. It modifies DEGREE of response,
## not presence of it, which is what an uptake-efficiency variant should do.
##
## A RESULT THAT CONTRADICTS THE OBVIOUS READING. The censored strains are not
## the most extreme responders. Class means on the vst scale are
## non-responsive +0.022, censored -0.038, significant decliners -0.051. A
## strain that vanished has a LESS extreme value than one that fell sharply
## from a high control frequency, because the delta is bounded by where it
## started. Coding the censored strains as "most responsive" is mis-ordered,
## not merely lossy.
##
## CAVEAT ON ABSOLUTE LOD. There is no GEMMA binary here, so the scans are an
## EMMAX-style LOCO LMM: one REML delta per chromosome rather than GEMMA's
## per-marker lambda, kinship from thinned scan markers. Against the shipped
## GEMMA scan this reproduces Spearman 0.969 and the same top marker on five of
## six chromosomes, but peak LOD runs 1-2 lower (III 6.82 vs 8.68, IV 6.75 vs
## 8.84, X 7.23 vs 7.83). The four traits are comparable TO EACH OTHER, not to
## the published numbers. Run the traits in cluster/pos1_binary_mapping_traits.csv
## through the GEMMA pipeline for numbers that are.
##
## Reads the cached scans; set POS1_BIN_REFIT=1 to recompute them, which needs
## a local CeNDR PLINK panel (CENDR_PLINK_2021) and plink2 on PATH.
## Exploratory; nothing in the manuscript reads this.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
DIAG <- "plots/diagnostics"
PH   <- "supplemental_data/phenotypes"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

S   <- fread(file.path(DIAG, "CACHE_pos1_binary_scans.tsv.gz"))
SUM <- fread(file.path(DIAG, "CACHE_pos1_binary_summary.tsv"))
LAB <- c(vst = "continuous (vst), 231", binary = "binary, 231",
         vst_testable = "continuous, 150 testable", binary_testable = "binary, 150 testable")
S[, trait_lab := factor(LAB[trait], levels = LAB)]
BONF <- 6.97; EIG <- 4.60
CHRLEV <- c("I","II","III","IV","V","X")
S[, chr := factor(chr, levels = CHRLEV)]

## --- the response classes, recomputed so the panel is self-contained --------
tr   <- fread(cmd = paste("gzcat", shQuote(file.path(PH, "pos1_2023_association_traits.csv.gz"))))
keep <- tr[is.finite(`vst_ctrl_pos-1_T2`)]$strain
raw  <- fread(cmd = paste("gzcat", shQuote(file.path(PH, "pos1_2023_sample_frequencies.csv.gz"))))
d5   <- raw[depth_cutoff == 5 & strain %in% keep][
             , .(frq = sum(frq), delta_ctrl = sum(delta_ctrl)), by = .(strain, sample_info, rnai)]
zz <- merge(d5[rnai == "pos-1", .(all_zero = all(frq == 0)), by = strain],
            d5[rnai != "pos-1", .(ctrl = mean(frq)), by = strain], by = "strain")
purged <- zz[all_zero & ctrl > 0]$strain
bb <- d5[rnai == "pos-1", {
  t <- tryCatch(t.test(delta_ctrl), error = function(e) NULL)
  .(pv = if (is.null(t)) NA_real_ else t$p.value, down = mean(delta_ctrl) < 0) }, by = strain]
bb[, q := p.adjust(pv, "BH")]
sig <- bb[q < 0.05 & down]$strain
cls <- data.table(strain = keep)
cls[, class := fifelse(strain %in% purged, "censored\n(absent from every pos-1 pool)",
              fifelse(strain %in% sig, "significant decline", "not individually responsive"))]
cls <- merge(cls, tr[, .(strain, vst = `vst_ctrl_pos-1_T2`)], by = "strain")
cls[, class := factor(class, levels = c("not individually responsive", "censored\n(absent from every pos-1 pool)",
                                        "significant decline"))]
cmn <- cls[, .(m = mean(vst), n = .N), by = class][order(class)]

COL <- c("continuous (vst), 231" = "#1B3A6B", "binary, 231" = "#C4302B",
         "continuous, 150 testable" = "#7A9CC6", "binary, 150 testable" = "#E8A29C")

## --- A. mirrored Manhattan, continuous above the axis and binary below ------
M <- S[trait %in% c("vst", "binary")]
M[, lod := -log10(p)]
M[, y := fifelse(trait == "vst", lod, -lod)]
pA <- ggplot(M, aes(ps / 1e6, y, colour = trait_lab)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey30") +
  geom_hline(yintercept = c(BONF, -BONF), linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_hline(yintercept = c(EIG, -EIG), linetype = "dotted", linewidth = 0.3, colour = "grey60") +
  geom_point(size = 0.28, alpha = 0.5) +
  facet_grid(~ chr, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = COL, name = NULL) +
  scale_y_continuous(labels = abs) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "position (Mb)", y = "LOD",
       title = "The same 231 strains, scored continuously (above) and as responsive / not (below)",
       subtitle = sprintf(paste("Dashed = Bonferroni %.2f, dotted = eigen %.2f. Continuous clears the eigen threshold at %d markers and",
                                "\nBonferroni at %d; binary clears eigen at %d and Bonferroni at %d. Every chromosome's top marker changes."),
                          BONF, EIG, SUM[trait=="vst"]$n_eigen, SUM[trait=="vst"]$n_bonf,
                          SUM[trait=="binary"]$n_eigen, SUM[trait=="binary"]$n_bonf)) +
  theme_bw(9) + theme(legend.position = "top", panel.grid.minor = element_blank(),
                      panel.spacing = unit(1.5, "pt"))

## --- B. the 2x2 at the three reported loci ---------------------------------
FOC <- c("III:5965738" = "III:5.97 Mb", "IV:15323414" = "IV:15.32 Mb", "X:4875969" = "X:4.88 Mb")
L <- S[rs %in% names(FOC), .(rs, trait, lod = -log10(p))]
L[, locus := factor(FOC[rs], levels = FOC)]
L[, scale_ := factor(fifelse(grepl("^binary", trait), "binary", "continuous"),
                     levels = c("continuous", "binary"))]
L[, set_ := factor(fifelse(grepl("testable", trait), "150 testable", "all 231"),
                   levels = c("all 231", "150 testable"))]
pB <- ggplot(L, aes(scale_, lod, group = set_, colour = set_)) +
  geom_hline(yintercept = BONF, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_hline(yintercept = EIG, linetype = "dotted", linewidth = 0.3, colour = "grey60") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.8) +
  ## labels nudged apart by series -- the two lines converge at the binary end
  ## and the values collide if both sit directly above their point
  geom_text(aes(label = sprintf("%.2f", lod),
                vjust = fifelse(set_ == "all 231", -1.1, 1.9),
                hjust = fifelse(scale_ == "continuous", 1.15, -0.15)),
            size = 2.3, show.legend = FALSE) +
  facet_wrap(~ locus) +
  scale_colour_manual(values = c("all 231" = "#1B3A6B", "150 testable" = "#7A9CC6"), name = NULL) +
  expand_limits(y = c(-0.6, 8.6)) +
  scale_x_discrete(expand = expansion(mult = 0.35)) +
  labs(x = NULL, y = "LOD",
       title = "Crossing the scale against the strain set separates the two effects",
       subtitle = paste("Holding the strain set fixed, continuous -> binary is the steep drop. Holding the scale fixed,",
                        "\nremoving the 81 censored strains costs about one LOD -- so those strains are informative, and",
                        "\ncollapsing magnitude is what loses the signal. Chromosome III loses it completely; IV and X survive weakly.")) +
  theme_bw(9) + theme(legend.position = "top", panel.grid.minor = element_blank())

## --- C. QQ, and what the deflation says ------------------------------------
Q <- S[, {
  o <- sort(-log10(p))
  .(obs = o, exp = -log10(ppoints(length(o), a = 0)[length(o):1])) }, by = trait_lab]
lam <- SUM[, .(trait_lab = factor(LAB[trait], levels = LAB), lambda_gc)]
pC <- ggplot(Q, aes(exp, obs, colour = trait_lab)) +
  geom_abline(slope = 1, linewidth = 0.3, colour = "grey45") +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(values = COL, name = NULL) +
  labs(x = "expected LOD", y = "observed LOD",
       title = "The continuous scan is deflated; the binary scan is calibrated but flat",
       subtitle = paste(sprintf("lambda_gc: continuous %.3f and %.3f, binary %.3f and %.3f.",
                                lam[1]$lambda_gc, lam[3]$lambda_gc, lam[2]$lambda_gc, lam[4]$lambda_gc),
                        "The continuous scans sit BELOW the null line, so",
                        "\nthe mixed model is over-correcting and the published p values are conservative rather than",
                        "\nanti-conservative. Binarising fixes the calibration and removes the signal along with it.")) +
  theme_bw(9) + theme(legend.position = "top", panel.grid.minor = element_blank())

## --- D. why: the censored strains are not the extreme ones ------------------
pD <- ggplot(cls, aes(class, vst, colour = class)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey45") +
  ## seeded position, so the jittered cloud is identical on every rebuild
  geom_point(position = position_jitter(width = 0.18, height = 0, seed = 1),
             size = 0.7, alpha = 0.55) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45, linewidth = 0.4, colour = "black") +
  ## means printed under the cloud rather than on the crossbar they describe
  geom_text(data = cmn, aes(class, -0.175, label = sprintf("mean %.3f\nn = %d", m, n)),
            inherit.aes = FALSE, size = 2.5, lineheight = 0.95) +
  expand_limits(y = -0.20) +
  scale_colour_manual(values = c("#9A9A9A", "#C4302B", "#1B3A6B"), guide = "none") +
  labs(x = NULL, y = "vst_ctrl_pos-1_T2",
       title = "The censored strains are not the most extreme, so the binary coding is mis-ordered",
       subtitle = paste("A strain absent from every pos-1 pool has a LESS extreme continuous value than one that declined",
                        "\nsignificantly from a high control frequency, because the change is bounded by where it started.",
                        "\nCalling both of them simply 'responsive' puts them in one class that the genotypes do not predict.")) +
  theme_bw(9) + theme(panel.grid.minor = element_blank())

fig <- pA / pB / pC / pD + plot_annotation(tag_levels = "A")
ggsave(file.path(DIAG, "DIAG_pos1_binary_trait_mapping.pdf"), fig, width = 9, height = 15)
ggsave(file.path(DIAG, "DIAG_pos1_binary_trait_mapping.png"), fig, width = 9, height = 15, dpi = 180)
msg("wrote DIAG_pos1_binary_trait_mapping.{pdf,png}")
print(SUM)
