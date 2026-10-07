## 2023 pos-1: the corrected association scan, against the deposited one -----
##
##   Rscript reanalysis_2023_pos1/scripts/04_mapping_figure.R
##
## Reads reanalysis_2023_pos1/mapping/vst_ctrl_pos1_T2_loco_results.csv.gz --
## GEMMA LOCO on the vst trait built from the pool-reference deconvolution --
## and draws it against supplemental_data/mapping/pos1_2023_gemma_loco.csv.gz,
## the deposited scan on the 367-reference trait.
##
## THRESHOLDS. Both the Bonferroni line and the eigenvalue line are recomputed
## for THIS panel. The effective number of independent tests is panel-specific:
## it falls with the number of strains, so reusing the deposited M_eff on a
## smaller panel would draw a line that belongs to a different experiment. The
## method is the one in scripts/eigen_independent_tests.R, reproduced here
## rather than imported because that script runs at top level against fixed
## paths: per chromosome, standardise the genotypes at the markers actually
## tested, take the eigenvalues of the marker correlation matrix through the
## n x n matrix that shares its non-zero eigenvalues, and sum the Li & Ji
## (2005) terms. trace(A) == M is asserted, as there.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(patchwork); library(ggtext)
})

DAT <- "reanalysis_2023_pos1/data"
MAP <- "reanalysis_2023_pos1/mapping"
OUT <- "reanalysis_2023_pos1/plots"
GENO <- "data/genotypes/CeNDR20210121_Plink"
PLINK <- path.expand("~/bin/plink")
TMP <- file.path(tempdir(), "eig"); dir.create(TMP, showWarnings = FALSE, recursive = TRUE)
VAR_FRAC <- 0.995
CHROMS <- c("I", "II", "III", "IV", "V", "X")
ALL_LEN <- c(I = 15072434, II = 15279421, III = 13783801,
             IV = 17493829, V = 20924180, X = 17718942)
source("scripts/figure_palette.R"); source("scripts/figure_theme.R")
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

new <- fread(cmd = paste("gzcat", shQuote(file.path(MAP, "vst_ctrl_pos1_T2_loco_results.csv.gz"))))
old <- fread(cmd = "gzcat supplemental_data/mapping/pos1_2023_gemma_loco.csv.gz")
tr  <- fread(file.path(DAT, "mapping_traits_dp5.csv"))
strains <- tr[!is.na(vst_ctrl_pos1_T2)]$strain

## The panel size is taken from GEMMA's own logs rather than assumed. GEMMA
## drops individuals with a missing phenotype, so the 224 strains plink kept
## become the 184 that carry a vst value -- and the eigenvalue threshold has to
## be computed on the individuals actually analysed, not on everyone genotyped.
logs <- list.files(MAP, pattern = "^gemma_LOCOmapping.*log\\.txt$", full.names = TRUE)
stopifnot(length(logs) == 6)
nind <- vapply(logs, function(f) {
  l <- grep("number of analyzed individuals", readLines(f), value = TRUE)
  as.integer(sub(".*=\\s*", "", l))
}, integer(1))
ntot <- vapply(logs, function(f) {
  l <- grep("number of total individuals", readLines(f), value = TRUE)
  as.integer(sub(".*=\\s*", "", l))
}, integer(1))
stopifnot(length(unique(nind)) == 1, length(unique(ntot)) == 1)
msg("GEMMA logs: ", unique(ntot), " individuals supplied, ", unique(nind), " analysed")
stopifnot(unique(nind) == length(strains))
msg("new scan: ", format(nrow(new), big.mark = ","), " markers | panel ", length(strains), " strains")

## --- M_eff for this panel, by the method of eigen_independent_tests.R -------
chrom_eigen <- function(chrom, strains, markers) {
  keep <- file.path(TMP, "keep.txt"); ext <- file.path(TMP, "ext.txt")
  writeLines(paste(strains, strains), keep); writeLines(markers, ext)
  out <- file.path(TMP, paste0("c_", chrom))
  st <- system(sprintf(paste("%s --bfile %s --keep %s --extract %s --recode A",
                             "--allow-extra-chr --out %s --silent"),
                       PLINK, file.path(GENO, chrom), keep, ext, out))
  raw <- paste0(out, ".raw")
  if (st != 0 || !file.exists(raw)) { warning("plink failed: ", chrom); return(NULL) }
  d <- fread(raw); Z <- as.matrix(d[, -(1:6)]); n <- nrow(Z)
  cm <- colMeans(Z, na.rm = TRUE); allna <- is.nan(cm)
  if (any(allna)) { Z <- Z[, !allna, drop = FALSE]; cm <- cm[!allna] }
  idx <- which(is.na(Z), arr.ind = TRUE); if (nrow(idx)) Z[idx] <- cm[idx[, 2]]
  Z <- Z[, apply(Z, 2, sd) > 0, drop = FALSE]; M <- ncol(Z)
  Zc <- scale(Z); A <- tcrossprod(Zc) / (n - 1)
  lam <- eigen(A, symmetric = TRUE, only.values = TRUE)$values; lam[lam < 0] <- 0
  stopifnot(abs(sum(lam) - M) / M < 1e-6)
  file.remove(list.files(TMP, pattern = paste0("c_", chrom), full.names = TRUE))
  data.table(chrom = chrom, n_strain = n, n_marker = M,
             M_eff_liji = sum(as.numeric(lam >= 1) + (lam - floor(lam))))
}

eig_file <- file.path(DAT, "eigen_independent_tests_pool.tsv")
if (file.exists(eig_file)) {
  eig <- fread(eig_file); msg("M_eff read from ", eig_file)
} else {
  msg("computing M_eff per chromosome (plink + eigen)")
  eig <- rbindlist(lapply(CHROMS, function(ch)
    chrom_eigen(ch, strains, new[chr == ch]$rs)))
  fwrite(eig, eig_file, sep = "\t", na = "NA", quote = FALSE)
}
M_EFF <- sum(eig$M_eff_liji)
BF  <- -log10(0.05 / nrow(new))
EIG <- -log10(0.05 / M_EFF)
msg("thresholds: Bonferroni ", round(BF, 2), " | eigen ", round(EIG, 2),
    " (M_eff = ", round(M_EFF), " over ", length(strains), " strains)")
print(eig)

## --- the two scans ----------------------------------------------------------
prep <- function(g, lab) g[, .(chrom = factor(chr, levels = CHROMS), ps,
                               lp = -log10(p_wald), panel = lab)]
old_bf <- -log10(0.05 / nrow(old)); old_eig <- -log10(0.05 / 1972)
d <- rbind(prep(old, "deposited\n367-strain reference\n231 strains"),
           prep(new, sprintf("pool reference\n224-strain reference\n%d strains", length(strains))))
d[, panel := factor(panel, levels = unique(panel))]
thr <- data.table(panel = factor(levels(d$panel), levels = levels(d$panel)),
                  bf = c(old_bf, BF), eig = c(old_eig, EIG))
span <- rbind(data.table(chrom = factor(CHROMS, levels = CHROMS), ps = unname(ALL_LEN), lp = 0),
              data.table(chrom = factor(CHROMS, levels = CHROMS), ps = 0, lp = 0))

fig <- ggplot(d, aes(ps / 1e6, lp)) +
  geom_blank(data = span[rep(1:nrow(span), each = 2)][
               , panel := rep(factor(levels(d$panel), levels = levels(d$panel)), nrow(span))],
             aes(ps / 1e6, lp)) +
  geom_hline(data = thr, aes(yintercept = bf), linetype = "dashed",
             linewidth = 0.3, colour = "grey45") +
  geom_hline(data = thr, aes(yintercept = eig), linetype = "dotted",
             linewidth = 0.45, colour = "#1B7837") +
  geom_point(size = 0.35, alpha = 0.5, colour = COL_PT) +
  facet_grid(panel ~ chrom, scales = "free_x", space = "free_x", switch = "x") +
  scale_x_continuous(breaks = seq(0, 25, 5), expand = expansion(mult = 0.02),
                     guide = guide_axis(check.overlap = TRUE)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = "Genomic position (Mb)", y = "&minus;log<sub>10</sub>*p*") +
  theme_pub() +
  theme(axis.title.y = element_markdown(), strip.placement = "outside",
        strip.text.y = element_text(size = 8, angle = 0, hjust = 0), panel.spacing.x = unit(1.5, "pt"))

ggsave(file.path(OUT, "mapping_old_vs_pool_reference.png"), fig,
       width = 13, height = 6, dpi = 200, bg = "white")
ggsave(file.path(OUT, "mapping_old_vs_pool_reference.pdf"), fig,
       width = 13, height = 6, device = cairo_pdf)
msg("wrote mapping_old_vs_pool_reference")

## --- what survived ----------------------------------------------------------
pk <- merge(old[, .(old = round(max(-log10(p_wald)), 2),
                    old_at = ps[which.max(-log10(p_wald))]), by = chr],
            new[, .(new = round(max(-log10(p_wald)), 2),
                    new_at = ps[which.max(-log10(p_wald))]), by = chr], by = "chr")
pk[, `:=`(old_over_bf = old > old_bf, new_over_bf = new > BF,
          new_over_eig = new > EIG)]
cat("\n== peak per chromosome ==\n"); print(pk)
cat(sprintf("\nthresholds  deposited: Bonferroni %.2f, eigen %.2f (M_eff 1972, 231 strains)\n",
            old_bf, old_eig))
cat(sprintf("            corrected: Bonferroni %.2f, eigen %.2f (M_eff %.0f, %d strains)\n",
            BF, EIG, M_EFF, length(strains)))
