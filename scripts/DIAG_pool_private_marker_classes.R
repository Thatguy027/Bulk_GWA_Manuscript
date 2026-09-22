## Is a strain's private-marker count one quantity, or two? --------------------
##
##   Rscript scripts/DIAG_pool_private_marker_classes.R
##     -> plots/diagnostics/DIAG_pool_private_marker_classes.{pdf,png}
##     -> plots/diagnostics/pool_private_marker_classes.tsv
##
## DIAG_baugh_rmsd_predictors.R established that the number of markers at which
## a strain is the sole carrier of the alternate allele is the strongest single
## predictor of that strain's deconvolution error, and that it absorbs the IBS
## effect. It treated that count as one number. This asks whether it is really
## two numbers with different design consequences.
##
## THE DISTINCTION. A hyper-divergent haplotype generates private markers in
## bulk, but it does not generate them alone: the same region also carries
## markers that EVERY divergent strain shares against the reference. Those
## shared markers are what builds a population-structure axis, and an axis is
## what a mixed model spends its power removing. Rare variation scattered
## outside divergent regions carries no such companion: a singleton is close to
## orthogonal to the top eigenvectors by construction.
##
## So privateness inside a strain's own divergent regions is bought with a
## structure tax, and privateness outside them is not. If the two classes
## predict deconvolution error equally well, panel design should prefer the
## untaxed kind. If only the divergent kind predicts it, there is no free lunch
## and the tax has to be paid. That is the question here, and it is worth
## answering before any panel-optimisation machinery gets built on top of it.
##
## WHAT IS COMPUTED
##   n_private_div     private markers inside that strain's OWN divergent
##                     regions, as called by CaeNDR per strain
##   n_private_nondiv  private markers outside them
##   n_shared_div      markers inside its own divergent regions where it is alt
##                     but NOT alone -- the structure tax, counted directly
##   pc1_load, pc1_abs loading on the first eigenvector of the panel genotype
##                     correlation, the tax measured a second way
##
## Divergent calls are per strain, not a union across strains, so "inside a
## divergent region" means inside one of that strain's own, which is the right
## question for a marker it alone carries.
##
## SOURCES OUTSIDE THE REPOSITORY -- CENDR_DIVERGENT, and the deposited Baugh
## input under data/. Exploratory: this is on the pool_optimization branch and
## nothing in the manuscript depends on it.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})

DEP <- "data/baugh/2024bootstrapINPUT.Rdata"
OLD <- Sys.getenv("BAUGH_BOOT", "/Users/Stefan/UCLA/Projects/bulkGWAS/baugh_wgs/cluster_data/20220908_Baugh_BulkL1_Bootstrap_Input_flippedCommon_NAfix.RData")
DIV <- Sys.getenv("CENDR_DIVERGENT",
                  "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213/20231213_c_elegans_divergent_regions_strain.bed")
DEC  <- "supplemental_data/deconvolution"
DIAG <- "plots/diagnostics"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)
for (f in c(DEP, OLD, DIV)) if (!file.exists(f)) stop("missing source: ", f, call. = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## --- the dep103 genotype matrix, exactly as make_baugh_deposited_fits.R ------
e <- new.env(); load(DEP, e)
gt <- e$flipped_bootstrap_input[[1]]
rm(e); gc()
gt <- gt[rowSums(is.na(gt)) == 0, ]
colnames(gt) <- sub("_.*$", "", colnames(gt))
o <- new.env(); load(OLD, o)
go <- o$flipped_bootstrap_input[[1]]
gt <- cbind(gt, PB306 = go[rownames(gt), grep("^PB306_", colnames(go))])
rm(o, go); gc()
msg(sprintf("dep103 matrix: %d markers x %d strains", nrow(gt), ncol(gt)))

## --- privateness, with positions this time ---------------------------------
alt <- rowSums(gt)
pidx <- which(alt == 1)
carrier <- colnames(gt)[max.col(gt[pidx, ], ties.method = "first")]
mk <- tstrsplit(sub("_.*$", "", rownames(gt)[pidx]), ":", fixed = TRUE)
P <- data.table(strain = carrier, chrom = mk[[1]], pos = as.integer(mk[[2]]))
P[, `:=`(start = pos, end = pos)]
msg(sprintf("private markers: %d of %d (%.1f%%)", nrow(P), nrow(gt),
            100 * nrow(P) / nrow(gt)))

## --- CaeNDR per-strain divergent regions -----------------------------------
## BED is half-open and 0-based; +1 on start puts it on the VCF's 1-based
## coordinates, which is what the marker names carry.
bed <- fread(DIV, col.names = c("chrom", "start", "end", "strain"))
bed <- bed[strain %in% colnames(gt)]
bed[, start := start + 1L]
setkey(bed, strain, chrom, start, end)
msg(sprintf("divergent regions: %d over %d of %d panel strains (N2 has none, by definition)",
            nrow(bed), uniqueN(bed$strain), ncol(gt)))

hit <- foverlaps(P, bed, by.x = c("strain", "chrom", "start", "end"), nomatch = NULL)
P[, in_div := FALSE]
P[unique(hit[, .(strain, chrom, pos = i.start)]), in_div := TRUE,
  on = .(strain, chrom, pos)]

## 0-fill, because a strain with no private marker at all must show as 0 and
## not drop out of the table. PB306 is exactly that strain -- see the console.
priv <- P[, .(n_private = .N, n_private_div = sum(in_div),
              n_private_nondiv = sum(!in_div)), by = strain]
priv <- priv[data.table(strain = colnames(gt)), on = "strain"]
for (v in c("n_private", "n_private_div", "n_private_nondiv"))
  set(priv, which(is.na(priv[[v]])), v, 0L)
if (any(priv$n_private == 0L))
  msg("strains with NO private marker: ", paste(priv[n_private == 0L]$strain, collapse = ", "))

## --- the structure tax, counted directly -----------------------------------
## Markers inside a strain's own divergent regions where it carries alt but is
## NOT the only one to. These are the markers that give divergent strains a
## shared axis, and they are the price of the private markers beside them.
amk <- tstrsplit(sub("_.*$", "", rownames(gt)), ":", fixed = TRUE)
A <- data.table(chrom = amk[[1]], pos = as.integer(amk[[2]]), alt = alt,
                row = seq_len(nrow(gt)))
A[, `:=`(start = pos, end = pos)]
shared <- rbindlist(lapply(colnames(gt), function(s) {
  r <- which(gt[, s] == 1L & alt > 1L)
  if (!length(r)) return(data.table(strain = s, n_shared_div = 0L))
  x <- A[r][, strain := s]
  h <- foverlaps(x, bed, by.x = c("strain", "chrom", "start", "end"), nomatch = NULL)
  data.table(strain = s, n_shared_div = uniqueN(h$row))
}))
priv <- merge(priv, shared, by = "strain", all.x = TRUE)

## --- structure tax, measured a second way: loading on PC1 ------------------
pcs <- prcomp(t(gt[seq(1, nrow(gt), by = 20), ]), center = TRUE, scale. = FALSE)
pc <- data.table(strain = rownames(pcs$x), pc1 = pcs$x[, 1], pc2 = pcs$x[, 2])
pc[, `:=`(pc1_abs = abs(pc1 - median(pc1)))]
priv <- merge(priv, pc, by = "strain", all.x = TRUE)
msg(sprintf("PC1 explains %.1f%% of genotype variance (every 20th marker)",
            100 * pcs$sdev[1]^2 / sum(pcs$sdev^2)))

## --- validate against the deposited private-marker counts ------------------
dep <- fread(file.path(DEC, "baugh_strain_private_markers.tsv"))
v <- merge(priv, dep[, .(strain, dep_private = n_private)], by = "strain")
msg(sprintf("VALIDATION against baugh_strain_private_markers.tsv: %d strains, max abs diff %d",
            nrow(v), max(abs(v$n_private - v$dep_private))))
if (max(abs(v$n_private - v$dep_private)) != 0)
  stop("private-marker counts do not reproduce the deposited table", call. = FALSE)

## --- join the error measures ----------------------------------------------
nn <- fread(cmd = paste("gzcat", shQuote(file.path(DEC, "baugh_nnls_dep103_with_mipseq.tsv.gz"))))
d <- nn[!is.na(frq) & !is.na(published_frq) & strain != "N2",
        .(rmsd = sqrt(mean((frq - published_frq)^2)),
          mean_f = mean(published_frq), zero_frac = mean(frq == 0)), by = strain]
sim <- fread(file.path(DEC, "baugh_strain_similarity.tsv"))
d <- merge(merge(d, sim, by = "strain", all.x = TRUE), priv, by = "strain", all.x = TRUE)
d[, `:=`(rel = rmsd / mean_f, div_frac = n_private_div / n_private)]
fwrite(d, file.path(DIAG, "pool_private_marker_classes.tsv"), sep = "\t")
msg(sprintf("error table: %d strains", nrow(d)))

## --- the question ----------------------------------------------------------
sp <- function(x, y) { ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
                       sprintf("%+.3f (p = %s)", ct$estimate, signif(ct$p.value, 2)) }
cat("\n== Spearman against per-strain RMSD, 99 strains, N2 excluded ==\n")
for (v in c("n_private", "n_private_div", "n_private_nondiv", "n_shared_div",
            "div_frac", "pc1_abs", "nn_ibs_wild"))
  cat(sprintf("  %-18s %s\n", v, sp(d[[v]], d$rmsd)))

cat("\n== do the two classes carry the same information? ==\n")
cat(sprintf("  n_private_div vs n_private_nondiv   rho = %s\n",
            sp(d$n_private_div, d$n_private_nondiv)))
r <- function(f) summary(lm(f, data = d))$adj.r.squared
cat(sprintf("  rank RMSD ~ divergent private only      adj R2 = %.3f\n",
            r(rank(rmsd) ~ rank(n_private_div))))
cat(sprintf("  rank RMSD ~ non-divergent private only  adj R2 = %.3f\n",
            r(rank(rmsd) ~ rank(n_private_nondiv))))
cat(sprintf("  rank RMSD ~ both                        adj R2 = %.3f\n",
            r(rank(rmsd) ~ rank(n_private_div) + rank(n_private_nondiv))))
cat(sprintf("  rank RMSD ~ total private (the old fit) adj R2 = %.3f\n",
            r(rank(rmsd) ~ rank(n_private))))
cat("\n  both-classes model:\n")
print(round(summary(lm(rank(rmsd) ~ rank(n_private_div) + rank(n_private_nondiv),
                       data = d))$coefficients[, c(1, 4)], 4))

cat("\n== the tax, per unit benefit ==\n")
cat(sprintf("  shared-divergent markers per private-divergent marker: median %.1f\n",
            median(d$n_shared_div / pmax(d$n_private_div, 1), na.rm = TRUE)))
cat(sprintf("  n_shared_div vs n_private_div  rho = %s\n",
            sp(d$n_shared_div, d$n_private_div)))
cat(sprintf("  n_shared_div vs n_private_nondiv  rho = %s\n",
            sp(d$n_shared_div, d$n_private_nondiv)))

## Can you buy identifiability without buying structure? Only to the extent
## that benefit varies AT A FIXED tax. Regress log benefit on log tax and read
## the spread of the residual: that spread is the whole design margin.
cat("\n== how much benefit is available at fixed tax? ==\n")
fit <- lm(log10(n_private_nondiv + 1) ~ log10(n_shared_div + 1), data = d)
d[, eff := residuals(fit)]
cat(sprintf("  log-log slope %.2f, R2 %.3f -- tax explains %.0f%% of benefit\n",
            coef(fit)[2], summary(fit)$r.squared, 100 * summary(fit)$r.squared))
cat(sprintf("  residual spread: 10th-90th percentile is a %.0f-fold range in benefit at the same tax\n",
            10^diff(quantile(d$eff, c(0.1, 0.9)))))
cat("  most efficient (high benefit for their tax):\n")
print(d[order(-eff)][1:6, .(strain, n_private_nondiv, n_shared_div,
                            fold_vs_expected = round(10^eff, 2), rmsd = signif(rmsd, 3))])
cat("  least efficient:\n")
print(d[order(eff)][1:6, .(strain, n_private_nondiv, n_shared_div,
                           fold_vs_expected = round(10^eff, 2), rmsd = signif(rmsd, 3))])
cat(sprintf("\n  does efficiency predict error?  rho = %s\n", sp(d$eff, d$rmsd)))

## --- figure ----------------------------------------------------------------
theme_set(theme_bw(9) + theme(
  plot.title = element_text(face = "bold", size = 10),
  plot.subtitle = element_text(size = 7.6, colour = "grey30"),
  panel.grid.minor = element_blank(),
  strip.background = element_rect(fill = "grey93"),
  strip.text = element_text(size = 7.3, face = "bold"),
  legend.position = "none"))

## A/B  the two private-marker classes against RMSD, on matched axes
CL <- data.table(
  v   = c("n_private_nondiv", "n_private_div"),
  lab = c("private markers OUTSIDE the strain's divergent regions",
          "private markers INSIDE the strain's divergent regions"))
L <- rbindlist(lapply(seq_len(nrow(CL)), function(i)
  data.table(strain = d$strain, x = log10(d[[CL$v[i]]] + 1), y = d$rmsd,
             panel = factor(CL$lab[i], CL$lab))))
ANN <- L[, { ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
             .(lab = sprintf("rho = %+.2f   (p = %s)", ct$estimate,
                             signif(ct$p.value, 2))) }, by = panel]
hl <- L[, .SD[order(-y)][1:3], by = panel]
## the loess is fitted only where there are points. PB306 sits alone at x = 0
## with nothing between it and x = 2, and a smoother run over all of L draws a
## confident hook across that gap that no data supports.
pAB <- ggplot(L, aes(x, y)) +
  geom_point(colour = "#2E4057", alpha = 0.7, size = 1.5) +
  geom_smooth(data = L[x > 0], method = "loess", formula = y ~ x, se = FALSE,
              colour = "#C4302B", linewidth = 0.5) +
  geom_text_repel(data = hl, aes(label = strain), size = 2.2, colour = "grey30",
                  min.segment.length = 0, max.overlaps = Inf, seed = 1) +
  geom_label(data = ANN, aes(x = -Inf, y = Inf, label = lab), hjust = -0.04,
             vjust = 1.15, size = 2.5, colour = "grey15", inherit.aes = FALSE,
             fill = "white", alpha = 0.85, label.size = 0,
             label.padding = unit(1.6, "pt")) +
  scale_y_log10() + facet_wrap(~ panel) +
  labs(x = "private markers (log10, +1)", y = "per-strain RMSD vs MIP-seq",
       title = "Only the untaxed half of privateness predicts deconvolution error",
       subtitle = paste("99 strains, N2 excluded. Same total count as the earlier fit, split by whether each private",
                        "\nmarker falls inside that strain's own CaeNDR divergent regions. Red line is a loess fit."))

## C  benefit against tax -- the design plane
d[, tax_ratio := n_shared_div / pmax(n_private_nondiv, 1)]
lab3 <- d[order(-n_private_nondiv)][1:4]
lab3 <- rbind(lab3, d[order(-n_shared_div)][1:3], d[n_private == 0L])
## both axes span three orders of magnitude, so both are log; the linear
## version puts 90 of 99 strains in one corner
pC <- ggplot(d, aes(n_shared_div + 1, n_private_nondiv + 1)) +
  geom_point(aes(size = n_private_div + 1), colour = "#2E4057", alpha = 0.6) +
  geom_text_repel(data = unique(lab3), aes(label = strain), size = 2.2,
                  colour = "grey30", min.segment.length = 0, max.overlaps = Inf,
                  seed = 1) +
  scale_size_continuous(range = c(0.8, 4.5), trans = "log10") +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks(sides = "bl", size = 0.25,
                      short = unit(0.04, "cm"), mid = unit(0.07, "cm"),
                      long = unit(0.11, "cm")) +
  labs(x = "structure tax: shared alt markers in own divergent regions (log10, +1)",
       y = "benefit: private markers outside\ndivergent regions (log10, +1)",
       title = "The design plane: benefit and tax are correlated, but loosely",
       subtitle = paste("Point size is privateness INSIDE divergent regions -- the taxed kind, which the left panels show does not predict",
                        "\nerror once the untaxed count is known. Upper left is what an optimiser wants and is largely empty: the most",
                        "\nidentifiable strains are also the most divergent. The exploitable variation is VERTICAL, within a tax level."))

fig <- pAB / pC + plot_layout(heights = c(1, 1.05))
ggsave(file.path(DIAG, "DIAG_pool_private_marker_classes.pdf"), fig, width = 9.6, height = 8.2)
ggsave(file.path(DIAG, "DIAG_pool_private_marker_classes.png"), fig, width = 9.6, height = 8.2, dpi = 200)
msg("wrote DIAG_pool_private_marker_classes.{pdf,png}")
