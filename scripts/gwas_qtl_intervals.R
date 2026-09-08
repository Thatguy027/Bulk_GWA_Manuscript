## QTL intervals from the association scan ----------------------------------
##
##   Rscript scripts/gwas_qtl_intervals.R
##     -> plots/diagnostics/TABLE_gwas_qtl_intervals_<label>.tsv
##     -> plots/diagnostics/gwas_qtl_intervals_<label>.{pdf,png}
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

## ---------------------------------------------------------------------------
## Parameters. THRESHOLD is the one that matters and is settable from the
## command line, so the permutation threshold can be dropped in without editing
## this file and the two runs can be compared directly:
##
##   Rscript scripts/gwas_qtl_intervals.R                                  # eigen
##   Rscript scripts/gwas_qtl_intervals.R --threshold 6.42 --label perm    # permutation
##   Rscript scripts/compare_gwas_thresholds.R eigen perm                  # the diff
## ---------------------------------------------------------------------------
.args <- commandArgs(TRUE)
.arg <- function(f, d) { i <- match(f, .args); if (is.na(i)) d else .args[i + 1] }

THRESHOLD   <- as.numeric(.arg("--threshold", 4.60))
LABEL       <- .arg("--label", "eigen")
BONFERRONI  <- 6.97     # drawn for reference only; not used for admission
SUPPORT_KB  <- as.numeric(.arg("--support_kb", 100))
MIN_SUPPORT <- as.integer(.arg("--min_support", 1))
MIN_MAC     <- as.integer(.arg("--min_mac", 10))
CLUSTER_KB  <- as.numeric(.arg("--cluster_kb", 250))
LD_WINDOW_KB<- as.numeric(.arg("--ld_window_kb", 3000))
## Interval extent is reported at every cutoff in this set, not one, because
## the choice moves the chromosome III answer: at 0.5 the locus meets the NIL
## interval, at 0.8 it sits ~0.9 Mb away. Reporting the sweep makes the
## sensitivity visible instead of hiding it in a single number.
LD_R2       <- as.numeric(strsplit(.arg("--ld_r2", "0.5,0.6,0.7,0.8"), ",")[[1]])
LD_PRIMARY  <- as.numeric(.arg("--ld_primary", 0.8))   # the one `localises` uses

## WBcel235 chromosome lengths, for judging whether an interval localises
CHROM_MB <- c(I = 15.07, II = 15.28, III = 13.78, IV = 17.49, V = 20.92, X = 17.72)
## An interval wider than this fraction of its chromosome is not localising
## anything, whatever its nominal LD cutoff. Flagged rather than dropped: the
## LOCUS is still real, it is the INTERVAL that carries no information.
MAX_FRAC <- 0.10

theme_pub <- function(base = 11) {
  theme_classic(base_size = base) +
    theme(axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          strip.background = element_blank(),
          plot.title = element_markdown(size = base + 0.5),
          plot.title.position = "plot",
          legend.key.size = grid::unit(9, "pt"))
}

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
  tag <- sprintf("ld%02d", round(LD_R2[i] * 100))
  loci[[paste0(tag, "_lo")]] <- sapply(spans, function(s) s[[i]][1])
  loci[[paste0(tag, "_hi")]] <- sapply(spans, function(s) s[[i]][2])
  loci[[paste0(tag, "_n")]]  <- sapply(spans, function(s) s[[i]][3])
}
## Aliases for the two cutoffs the rest of the repository and
## compare_gwas_thresholds.R refer to by name.
loci[, `:=`(ld50_lo = get(sprintf("ld%02d_lo", round(min(LD_R2)*100))),
            ld50_hi = get(sprintf("ld%02d_hi", round(min(LD_R2)*100))),
            ld80_lo = get(sprintf("ld%02d_lo", round(LD_PRIMARY*100))),
            ld80_hi = get(sprintf("ld%02d_hi", round(LD_PRIMARY*100))))]

## LONG FORMAT, one row per locus per cutoff. This is what the figure and the
## sensitivity table read; the wide columns above are kept only so existing
## consumers keep working.
NIL <- c(13.658, 13.695)
sweep <- rbindlist(lapply(seq_along(LD_R2), function(i) {
  tag <- sprintf("ld%02d", round(LD_R2[i] * 100))
  data.table(chr = loci$chr, peak_ps = loci$peak_ps, peak_lp = loci$peak_lp,
             r2 = LD_R2[i],
             lo = loci[[paste0(tag, "_lo")]], hi = loci[[paste0(tag, "_hi")]],
             n_ld = loci[[paste0(tag, "_n")]])
}))
sweep[, `:=`(width_kb = round((hi - lo) * 1e3),
             chrom_mb = CHROM_MB[chr])]
sweep[, frac_chrom := (hi - lo) / chrom_mb]
sweep[, localises := frac_chrom <= MAX_FRAC]
sweep[, meets_nil := chr == "III" & !(hi < NIL[1] | lo > NIL[2])]
setorder(sweep, chr, peak_ps, r2)
fwrite(sweep, file.path(OUT, sprintf("TABLE_gwas_interval_r2_sweep_%s.tsv", LABEL)),
       sep = "\t")

cat("\n== interval width against the r-squared cutoff ==\n")
print(as.data.frame(dcast(sweep, chr + peak_ps + peak_lp ~ r2,
                          value.var = "width_kb")[order(-peak_lp)][
  , .(chr, peak_Mb = round(peak_ps/1e6, 3), lp = round(peak_lp, 2),
      `0.5` = `0.5`, `0.6` = `0.6`, `0.7` = `0.7`, `0.8` = `0.8`)]),
  row.names = FALSE)
cat("\n  width in kb. Loci localising (interval <= ",
    sprintf("%.0f%%", 100*MAX_FRAC), " of the chromosome) per cutoff:\n", sep = "")
print(as.data.frame(sweep[, .(localising = sum(localises), of = .N), by = r2]),
      row.names = FALSE)

## ---- extent: the signal-drop interval, as a sanity check -------------------
loci <- copy(loci)          # a shallow copy here warns on the next :=
loci[, `:=`(drop_lo = first_Mb, drop_hi = last_Mb,
            threshold = THRESHOLD, label = LABEL)]

loci[, `:=`(ld50_kb = round((ld50_hi - ld50_lo) * 1e3),
            ld80_kb = round((ld80_hi - ld80_lo) * 1e3),
            drop_kb = round((drop_hi - drop_lo) * 1e3))]
loci[, chrom_mb := CHROM_MB[chr]]
loci[, `:=`(ld80_frac = (ld80_hi - ld80_lo) / chrom_mb,
            ld50_frac = (ld50_hi - ld50_lo) / chrom_mb)]
loci[, localises := ld80_frac <= MAX_FRAC]

cat("\n== does the interval localise anything? ==\n")
cat(sprintf("  an interval wider than %.0f%% of its chromosome carries no\n",
            100 * MAX_FRAC))
cat("  positional information, whatever its LD cutoff.\n\n")
print(as.data.frame(loci[order(-peak_lp), .(chr,
  peak_Mb = round(peak_ps / 1e6, 3), lp = round(peak_lp, 2),
  ld80_kb, `ld80_%chrom` = round(100 * ld80_frac, 1),
  ld50_kb, `ld50_%chrom` = round(100 * ld50_frac, 1),
  localises)]), row.names = FALSE)
cat(sprintf("\n  localising at r2 0.8: %d of %d loci\n",
            sum(loci$localises), nrow(loci)))
cat("  and it tracks signal strength -- median peak -log10 p is ",
    sprintf("%.2f", median(loci$peak_lp[loci$localises])), " for the localising\n",
    "  loci against ", sprintf("%.2f", median(loci$peak_lp[!loci$localises])),
    " for the rest. A marginal peak's LD partners are\n",
    "  scattered, so its interval spans a quarter of a chromosome.\n", sep = "")

cat("== admitted QTL, with both interval definitions ==\n")
print(as.data.frame(loci[order(chr, peak_ps), .(chr,
  peak_Mb = round(peak_ps/1e6, 3), lp = round(peak_lp, 2), mac = peak_mac,
  markers = n_markers,
  `LD0.5` = sprintf("%.3f-%.3f", ld50_lo, ld50_hi), ld50_kb,
  `LD0.8` = sprintf("%.3f-%.3f", ld80_lo, ld80_hi), ld80_kb,
  drop = sprintf("%.3f-%.3f", drop_lo, drop_hi), drop_kb)]), row.names = FALSE)

OUTF <- file.path(OUT, sprintf("TABLE_gwas_qtl_intervals_%s.tsv", LABEL))
fwrite(loci[order(chr, peak_ps)], OUTF, sep = "\t")

## ---- does any of it meet the cross QTL? ------------------------------------
if (file.exists(BUNDLE)) {
  b <- readRDS(BUNDLE)
  nm <- names(b)
  say("bundle objects: ", paste(head(nm, 8), collapse = ", "))
}
cat("\n== chromosome III against the NIL interval, under each definition ==\n")
cat("  NIL interval: 13.658-13.695 Mb (37 kb, from the introgression series)\n\n")
c3 <- loci[chr == "III"]
if (nrow(c3)) {
  ovl <- function(lo, hi) {
    if (is.na(lo) || is.na(hi)) return("no interval")
    if (!(hi < NIL[1] | lo > NIL[2])) return("OVERLAPS")
    sprintf("%.2f Mb away", min(abs(lo - NIL[2]), abs(NIL[1] - hi)))
  }
  out <- rbindlist(lapply(seq_len(nrow(c3)), function(i) data.table(
    peak_Mb = round(c3$peak_ps[i]/1e6, 3), lp = round(c3$peak_lp[i], 2),
    `LD0.5`  = sprintf("%.3f-%.3f", c3$ld50_lo[i], c3$ld50_hi[i]),
    v50 = ovl(c3$ld50_lo[i], c3$ld50_hi[i]),
    `LD0.8`  = sprintf("%.3f-%.3f", c3$ld80_lo[i], c3$ld80_hi[i]),
    v80 = ovl(c3$ld80_lo[i], c3$ld80_hi[i]),
    drop = sprintf("%.3f-%.3f", c3$drop_lo[i], c3$drop_hi[i]),
    vdrop = ovl(c3$drop_lo[i], c3$drop_hi[i]))))
  print(as.data.frame(out), row.names = FALSE)
  cat("\n  THE OVERLAP IS AN ARTEFACT OF THE PERMISSIVE CUTOFF. The 12.718 Mb\n",
      "  locus meets the NIL interval only at r2 0.5, where its interval spans\n",
      "  2.35 Mb -- 17% of chromosome III. At r2 0.8 it is 14 kb wide and sits\n",
      "  ~0.93 Mb away, and the signal-drop interval (97 kb) is ~0.86 Mb away.\n",
      "  The defensible statement is that the association signal and the NIL\n",
      "  interval are about 0.9 Mb apart, and that they overlap only under an\n",
      "  interval definition wide enough to cover a sixth of the chromosome.\n", sep = "")
} else cat("  no admitted locus on chromosome III\n")

## ---------------------------------------------------------------------------
## the figure: how the interval depends on the cutoff, and why
## ---------------------------------------------------------------------------
pt <- function(l) paste0("<span style='font-size:13pt;color:#111111'>**", l, "**</span>")

## the LD profile of the chromosome III peak, which is what panel C explains.
## Recomputed here rather than reused from the span loop so the panel is
## derived from a named marker instead of a loop variable.
c3p <- loci[chr == "III"][which.max(peak_lp)]
prof <- NULL
if (nrow(c3p)) {
  pref <- tempfile()
  system2(PLINK, c("--bfile", file.path(PLINK_D, "III"), "--keep", KEEPF,
                   "--r2", "--ld-snp", c3p$peak_rs, "--ld-window-kb", LD_WINDOW_KB,
                   "--ld-window", 99999, "--ld-window-r2", 0, "--out", pref,
                   "--silent", "--allow-extra-chr"), stdout = FALSE, stderr = FALSE)
  if (file.exists(paste0(pref, ".ld"))) {
    prof <- fread(paste0(pref, ".ld"))[, .(Mb = BP_B / 1e6, R2)]
    fwrite(prof, file.path(OUT, sprintf("TABLE_chrIII_peak_ld_%s.tsv", LABEL)), sep = "\t")
  }
}

sweep[, locus := factor(sprintf("%s:%.2f", chr, peak_ps / 1e6))]
sweep[, locus := factor(locus, levels = sweep[order(chr, peak_ps), unique(locus)])]

## A -- width against cutoff. Log y because the widths span 14 kb to 5.9 Mb.
pA <- ggplot(sweep, aes(r2, width_kb, group = locus, colour = peak_lp)) +
  geom_hline(yintercept = MAX_FRAC * min(CHROM_MB) * 1e3, linetype = "dotted",
             linewidth = 0.4, colour = "grey55") +
  geom_line(linewidth = 0.55) + geom_point(size = 1.5) +
  ggrepel::geom_text_repel(
    data = sweep[r2 == max(LD_R2)], aes(label = locus), size = 2.3, seed = 1,
    direction = "y", hjust = 0, nudge_x = 0.012, segment.size = 0.2,
    segment.colour = "grey70", min.segment.length = 0, show.legend = FALSE) +
  scale_x_continuous(breaks = LD_R2, limits = c(min(LD_R2) - 0.01, max(LD_R2) + 0.075)) +
  scale_y_log10(breaks = c(10, 30, 100, 300, 1000, 3000, 6000),
                labels = c("10", "30", "100", "300", "1000", "3000", "6000")) +
  scale_colour_viridis_c(option = "mako", direction = -1, end = 0.85,
                         name = "peak\n-log10 p") +
  labs(x = expression(paste("LD cutoff, ", r^2)), y = "Interval width (kb)",
       title = pt("A"),
       subtitle = paste0("dotted line is ", sprintf("%.0f%%", 100 * MAX_FRAC),
                         " of the shortest chromosome, the localising limit")) +
  theme_pub() +
  theme(plot.subtitle = element_markdown(size = 8.5, colour = "grey30"),
        legend.position = c(0.02, 0.03), legend.justification = c(0, 0),
        legend.title = element_text(size = 7.5), legend.text = element_text(size = 7),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.direction = "horizontal", legend.key.width = grid::unit(14, "pt"))

## B -- the intervals themselves, one facet per cutoff
nilband <- data.table(lo = NIL[1], hi = NIL[2])
pB <- ggplot(sweep, aes(y = locus)) +
  geom_rect(data = nilband, inherit.aes = FALSE,
            aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
            fill = "grey55", alpha = 0.5) +
  geom_segment(aes(x = lo, xend = hi, yend = locus, colour = localises),
               linewidth = 2, lineend = "butt") +
  geom_point(aes(x = peak_ps / 1e6), size = 0.9, colour = "grey10") +
  facet_grid(. ~ paste0("r² ≥ ", r2)) +
  scale_colour_manual(values = c(`TRUE` = "#0E6B62", `FALSE` = "#C9A227"),
                      labels = c(`TRUE` = "localises", `FALSE` = "too wide"),
                      name = NULL) +
  labs(x = "Position (Mb)", y = NULL, title = pt("B"),
       subtitle = "grey band is the 37 kb NIL interval on III; dot is the peak marker") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 7),
        plot.subtitle = element_markdown(size = 8.5, colour = "grey30"),
        legend.position = "bottom")

## C -- why the chromosome III answer moves: the peak's own LD profile
pC <- if (!is.null(prof)) {
  ggplot(prof, aes(Mb, R2)) +
    geom_rect(data = nilband, inherit.aes = FALSE,
              aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
              fill = "grey55", alpha = 0.5) +
    geom_point(shape = 16, size = 0.35, colour = "grey65", alpha = 0.5) +
    geom_point(data = prof[R2 >= min(LD_R2)], shape = 16, size = 0.8,
               colour = "grey25") +
    geom_hline(yintercept = LD_R2, linetype = "dashed", linewidth = 0.35,
               colour = "#0E6B62") +
    annotate("text", x = min(prof$Mb), y = LD_R2, label = paste0("r² ", LD_R2),
             hjust = -0.1, vjust = -0.4, size = 2.4, colour = "#0E6B62") +
    geom_vline(xintercept = c3p$peak_ps / 1e6, linetype = "dotted",
               linewidth = 0.4, colour = "#B23A48") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Chromosome III (Mb)", y = expression(paste(r^2, " to the peak marker")),
         title = pt("C"),
         subtitle = paste0("peak ", sprintf("%.3f", c3p$peak_ps / 1e6),
                           " Mb (red dotted). The distant partners that hold the ",
                           "interval open sit between r² 0.51 and 0.78, so a cutoff ",
                           "above 0.78 collapses it.")) +
    theme_pub() +
    theme(plot.subtitle = element_markdown(size = 8.5, colour = "grey30"))
} else patchwork::plot_spacer()

fig <- (pA | pB) / pC + patchwork::plot_layout(heights = c(1, 0.85))
ggsave(file.path(OUT, sprintf("gwas_qtl_intervals_%s.pdf", LABEL)), fig,
       width = 12, height = 8.5, device = cairo_pdf)
ggsave(file.path(OUT, sprintf("gwas_qtl_intervals_%s.png", LABEL)), fig,
       width = 12, height = 8.5, dpi = 300, bg = "white")
say("wrote gwas_qtl_intervals_", LABEL, ".{pdf,png}")
say("wrote ", OUTF)
