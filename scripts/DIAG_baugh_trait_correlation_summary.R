## How far each reconstruction sits from the published trait, at three levels -
##
## One table, matched comparisons only: every Slope reconstruction is scored
## against the published Slope and every PC1 reconstruction against the
## published PC1. Cross-family comparisons are not the question.
##
## The three levels answer different things and they do not agree, which is the
## reason to keep them side by side:
##
##   trait_value   Spearman of the per-strain trait values against the published
##                 trait. How well the phenotype was reconstructed.
##   scan_genome   Spearman of -log10 p across the genome, strongest marker per
##                 100 kb bin. How well a scan on it reproduces the published
##                 scan.
##   <interval>    the same on markers inside each published interval, matched
##                 to the trait family: Slope reconstructions get the two Slope
##                 intervals, PC1 reconstructions the two PC1 intervals.
##
## Trait-value correlation runs 0.81 to 0.90 while genome-wide scan correlation
## runs 0.52 to 0.72, and the ORDERING reverses between them -- delta Slope
## (MIP) is second best on trait values and worst on the scan. Neither predicts
## detection; see scripts/DIAG_baugh_interval_correlations.R for that.
##
## Reads the GEMMA output under data/, so this is a diagnostic and sits outside
## the deposit-only rebuild.
##
## Writes plots/diagnostics/baugh_correlation_summary.tsv
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

D    <- "data/baugh/reanalysis_mappings"
PH   <- "supplemental_data/phenotypes"
DIAG <- "plots/diagnostics"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)

REC <- data.table(
  recon = c("delta Slope (MIP)", "log-ratio Slope (WGS)", "delta Slope (WGS)",
            "delta PC1 (MIP)",   "log-ratio PC1 (WGS)",   "delta PC1 (WGS)"),
  col   = c("delta_slope_baugh", "slope_nnls", "delta_slope_nnls",
            "delta_pc1_baugh",   "pc1_nnls",   "delta_pc1_nnls"),
  ref   = rep(c("published_slope_baugh", "published_pc1_baugh"), each = 3),
  fam   = rep(c("Slope", "PC1"), each = 3))

IV <- data.table(
  id  = c("Slope IV", "Slope V", "PC1 V-left", "PC1 V-right"),
  fam = c("Slope", "Slope", "PC1", "PC1"),
  chr = c("IV", "V", "V", "V"),
  lo  = c(15939340, 15660911, 1345848, 15775895),
  hi  = c(16613710, 17615557, 2764788, 18065050))

## --- level 1: the trait values --------------------------------------------
tr <- fread(file.path(PH, "baugh_mapping_traits.csv"))
REC[, trait_value := mapply(function(a, b)
  cor(tr[[a]], tr[[b]], method = "spearman"), col, ref)]

## --- the scans -------------------------------------------------------------
A <- rbindlist(lapply(list.files(D, full.names = TRUE), function(f) {
  d <- fread(cmd = paste("gzcat", shQuote(f)))
  d[, .(trait = trait[1], chr, ps, p = p_wald, lp = -log10(p_wald))] }))
A <- A[!is.na(p)]

## --- level 2: genome-wide --------------------------------------------------
A[, bin := paste0(chr, "_", ps %/% 1e5)]
M <- dcast(A[, .(lp = max(lp)), by = .(trait, bin)], bin ~ trait, value.var = "lp")
REC[, scan_genome := mapply(function(a, b)
  cor(M[[a]], M[[b]], method = "spearman", use = "pairwise.complete.obs"), col, ref)]

## --- level 3: inside each published interval, matched to the family --------
for (i in seq_len(nrow(IV))) {
  W <- dcast(A[chr == IV$chr[i] & ps >= IV$lo[i] & ps <= IV$hi[i]],
             ps ~ trait, value.var = "lp")
  REC[[IV$id[i]]] <- ifelse(REC$fam == IV$fam[i],
    mapply(function(a, b)
      cor(W[[a]], W[[b]], method = "spearman", use = "pairwise.complete.obs"),
      REC$col, REC$ref), NA_real_)
}

out <- REC[, .(recon, published_trait = ref, trait_value, scan_genome,
               `Slope IV`, `Slope V`, `PC1 V-left`, `PC1 V-right`)]
fwrite(out, file.path(DIAG, "baugh_correlation_summary.tsv"), sep = "\t")
print(out[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], nrows = 10)
cat(sprintf("\nwrote %s/baugh_correlation_summary.tsv\n", DIAG))
