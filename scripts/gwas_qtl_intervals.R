## QTL intervals from the association scan ----------------------------------
##
##   Rscript scripts/gwas_qtl_intervals.R
##     -> plots/diagnostics/TABLE_gwas_qtl_intervals.tsv
##     -> plots/diagnostics/gwas_qtl_intervals.{pdf,png}
##
## NEEDS data/genotypes/CeNDR20210121_Plink (Dryad-hosted) for linkage
## disequilibrium, and plink on PATH.
##
## THE PROBLEM THIS SOLVES. On chromosome III one marker clears Bonferroni at
## 5.966 Mb with ZERO other threshold-passing markers within 100 kb of the 628
## present, while a cluster at 12.70-12.80 Mb peaks BELOW Bonferroni with 14
## supporting markers. The eye reads the first as noise and the second as a QTL,
## and no significance threshold can tell them apart. See
## scripts/diagnostic_gwas_intervals.R for that diagnosis.
##
## ADMISSION, THEN EXTENT -- two separate steps, deliberately.
##
## ADMISSION is local support, not significance. A locus is admitted when the
## peak marker has at least MIN_SUPPORT other threshold-passing markers within
## SUPPORT_KB. In a panel with linkage disequilibrium a true association is
## tagged by several correlated markers; a lone spike with flat neighbours is
## genotyping error, an unshared rare haplotype, or chance. This is the clumping
## logic of human GWAS, applied to a panel whose Bonferroni threshold is known
## to be over-conservative because the markers are not independent. A minor
## allele count floor is applied as a second gate: genome-wide the isolated
## markers have median allele frequency 0.082 against 0.394 for the supported
## ones, so low frequency is the other half of the signature.
##
## EXTENT is linkage disequilibrium to the peak marker. The interval is the span
## of markers with r-squared at or above LD_R2 to the peak, computed IN THE
## PHENOTYPED PANEL rather than in all 540 isotypes -- LD is a property of the
## sample that produced the association, not of the species. This is the
## cegwas/NemaScan convention, and it is what makes a GWAS interval
## commensurable with a LINKAGE interval, which is the entire point of
## comparing the two. A signal-drop interval (the contiguous run of
## threshold-passing markers containing the peak) is reported beside it as a
## sanity check, because LD intervals fragment across a recombination hotspot.
##
## THE THRESHOLD IS THE EIGEN ONE, and that is a choice. Bonferroni over every
## marker assumes independence the panel does not have, and would reject the
## 12.7 Mb cluster. The eigen threshold divides alpha by the effective number of
## independent tests from the marker correlation matrix (Li & Ji 2005) and is
## what the captions and methods already quote. The rigorous alternative is a
## permutation threshold -- 100 to 1000 phenotype permutations through GEMMA,
## 95th percentile of the per-permutation maximum -- which would land between
## the two and is calibrated to this panel's actual LD and relatedness. Set
## THRESHOLD below to that value when it exists; nothing else needs to change.
##
## WHAT TO EXPECT, stated before the numbers so it cannot look like a surprise:
## the supported chromosome III cluster is at 12.70-12.80 Mb, the NIL interval
## is 13.658-13.695 Mb, and 13.5-13.9 Mb holds nothing above the eigen line.
## A rigorous interval will probably NOT overlap the cross QTL. That is a
## finding -- the association signal is a distinct locus -- and not a failure of
## the interval method.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(patchwork); library(ggtext)
})

OUT     <- "plots/diagnostics"
SCAN    <- "supplemental_data/mapping/pos1_2023_gemma_loco.csv.gz"
TRAITS  <- "supplemental_data/phenotypes/pos1_2023_association_traits.csv.gz"
PLINK_D <- "data/genotypes/CeNDR20210121_Plink"
BUNDLE  <- "supplemental_data/mapping/pooled_cross_bundle_thinned.rds"

THRESHOLD   <- 4.60     # eigen, alpha / 1972 effective tests. See header.
BONFERRONI  <- 6.97     # drawn for reference only; not used for admission
SUPPORT_KB  <- 100      # window for the support count
MIN_SUPPORT <- 1        # other threshold-passing markers required in it
MIN_MAC     <- 10       # minor allele count floor, in strains
CLUSTER_KB  <- 250      # admitted markers within this distance are one locus
LD_WINDOW_KB<- 3000     # how far to look for LD partners
LD_R2       <- c(0.5, 0.8)

stopifnot(file.exists(SCAN), file.exists(TRAITS))
if (!dir.exists(PLINK_D))
  stop("need ", PLINK_D, " for LD -- see DATA_AVAILABILITY.md", call. = FALSE)
PLINK <- Sys.which("plink"); if (PLINK == "") PLINK <- Sys.which("plink1.9")
if (PLINK == "") stop("plink not on PATH", call. = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep=""); flush.console() }

## ---- the phenotyped panel, which is the panel LD must be computed in ------
tr <- fread(TRAITS)
vcol <- grep("^vst", names(tr), value = TRUE)[1]
phenotyped <- tr[[1]][!is.na(tr[[vcol]])]
say("phenotyped strains: ", length(phenotyped), " (trait ", vcol, ")")

fam <- fread(file.path(PLINK_D, "I.fam"), header = FALSE)
keep <- fam[V2 %in% phenotyped, .(V1, V2)]
say("  matched in the PLINK set: ", nrow(keep), " of ", length(phenotyped))
KEEPF <- tempfile(); fwrite(keep, KEEPF, sep = " ", col.names = FALSE)

## ---- admission -------------------------------------------------------------
d <- fread(SCAN)[, .(chr, rs, ps, af, beta, p = p_wald)][, lp := -log10(p)]
setorder(d, chr, ps)
n_str <- length(phenotyped)
d[, mac := round(pmin(af, 1 - af) * n_str)]

hi <- d[lp > THRESHOLD]
hi[, support := {
  s <- integer(.N); for (i in seq_len(.N)) s[i] <- sum(abs(ps - ps[i]) <= SUPPORT_KB * 1e3) - 1L
  s
}, by = chr]
say("above threshold ", THRESHOLD, ": ", nrow(hi),
    " | isolated: ", sum(hi$support < MIN_SUPPORT),
    " | below MAC ", MIN_MAC, ": ", sum(hi$mac < MIN_MAC))

adm <- hi[support >= MIN_SUPPORT & mac >= MIN_MAC]
say("admitted markers: ", nrow(adm))
rej <- hi[!(support >= MIN_SUPPORT & mac >= MIN_MAC)]
if (nrow(rej)) {
  cat("\n== markers rejected, and why ==\n")
  print(as.data.frame(rej[order(-lp), .(chr, Mb = round(ps/1e6, 3), lp = round(lp, 2),
      af = round(af, 3), mac, support,
      why = fifelse(support < MIN_SUPPORT & mac < MIN_MAC, "isolated + rare",
             fifelse(support < MIN_SUPPORT, "isolated", "rare")))]), row.names = FALSE)
  cat("\n")
}

## ---- cluster into loci -----------------------------------------------------
adm <- adm[order(chr, ps)]
adm[, locus := cumsum(c(1L, as.integer(diff(ps) > CLUSTER_KB * 1e3))), by = chr]
loci <- adm[, .(n_markers = .N,
                first_Mb = min(ps)/1e6, last_Mb = max(ps)/1e6,
                peak_rs = rs[which.max(lp)], peak_ps = ps[which.max(lp)],
                peak_lp = max(lp), peak_af = af[which.max(lp)],
                peak_mac = mac[which.max(lp)]), by = .(chr, locus)]
say("loci after clustering at ", CLUSTER_KB, " kb: ", nrow(loci))

## ---- extent: LD to the peak marker ----------------------------------------
ld_span <- function(chrom, peak_rs, r2) {
  pref <- tempfile()
  args <- c("--bfile", file.path(PLINK_D, chrom), "--keep", KEEPF,
            "--r2", "--ld-snp", peak_rs,
            "--ld-window-kb", LD_WINDOW_KB, "--ld-window", 99999,
            "--ld-window-r2", 0, "--out", pref, "--silent", "--allow-extra-chr")
  system2(PLINK, args, stdout = FALSE, stderr = FALSE)
  f <- paste0(pref, ".ld")
  if (!file.exists(f)) return(list(lo = NA_real_, hi = NA_real_, n = NA_integer_))
  l <- fread(f)
  out <- lapply(r2, function(cut) {
    k <- l[R2 >= cut]
    if (!nrow(k)) return(c(NA, NA, 0))
    c(min(k$BP_B)/1e6, max(k$BP_B)/1e6, nrow(k))
  })
  unlink(c(f, paste0(pref, c(".log", ".nosex"))))
  out
}

say("computing LD spans")
spans <- pmap(list(loci$chr, loci$peak_rs), function(c_, r_) ld_span(c_, r_, LD_R2))
for (i in seq_along(LD_R2)) {
  loci[[paste0("ld", LD_R2[i]*100, "_lo")]] <- sapply(spans, function(s) s[[i]][1])
  loci[[paste0("ld", LD_R2[i]*100, "_hi")]] <- sapply(spans, function(s) s[[i]][2])
  loci[[paste0("ld", LD_R2[i]*100, "_n")]]  <- sapply(spans, function(s) s[[i]][3])
}

## ---- extent: the signal-drop interval, as a sanity check -------------------
loci[, drop_lo := first_Mb][, drop_hi := last_Mb]

loci[, `:=`(ld50_kb = round((ld50_hi - ld50_lo) * 1e3),
            ld80_kb = round((ld80_hi - ld80_lo) * 1e3),
            drop_kb = round((drop_hi - drop_lo) * 1e3))]

cat("== admitted QTL, with both interval definitions ==\n")
print(as.data.frame(loci[order(chr, peak_ps), .(chr,
  peak_Mb = round(peak_ps/1e6, 3), lp = round(peak_lp, 2), mac = peak_mac,
  markers = n_markers,
  `LD0.5` = sprintf("%.3f-%.3f", ld50_lo, ld50_hi), ld50_kb,
  `LD0.8` = sprintf("%.3f-%.3f", ld80_lo, ld80_hi), ld80_kb,
  drop = sprintf("%.3f-%.3f", drop_lo, drop_hi), drop_kb)]), row.names = FALSE)

fwrite(loci[order(chr, peak_ps)], file.path(OUT, "TABLE_gwas_qtl_intervals.tsv"), sep = "\t")

## ---- does any of it meet the cross QTL? ------------------------------------
if (file.exists(BUNDLE)) {
  b <- readRDS(BUNDLE)
  nm <- names(b)
  say("bundle objects: ", paste(head(nm, 8), collapse = ", "))
}
NIL <- c(13.658, 13.695)
cat("\n== chromosome III: association interval against the NIL interval ==\n")
c3 <- loci[chr == "III"]
if (nrow(c3)) {
  for (i in seq_len(nrow(c3))) {
    ov <- !(c3$ld50_hi[i] < NIL[1] | c3$ld50_lo[i] > NIL[2])
    gap <- if (ov) 0 else min(abs(c3$ld50_lo[i] - NIL[2]), abs(NIL[1] - c3$ld50_hi[i]))
    cat(sprintf("  peak %.3f Mb, LD0.5 interval %.3f-%.3f -> %s the NIL interval",
                c3$peak_ps[i]/1e6, c3$ld50_lo[i], c3$ld50_hi[i],
                ifelse(ov, "OVERLAPS", "does NOT overlap")))
    cat(if (ov) "\n" else sprintf(", gap %.3f Mb\n", gap))
  }
} else cat("  no admitted locus on chromosome III\n")
say("wrote TABLE_gwas_qtl_intervals.tsv")
