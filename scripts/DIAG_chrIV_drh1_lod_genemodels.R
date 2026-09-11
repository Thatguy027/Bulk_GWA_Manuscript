## Diagnostic -- the chrIV QTL over its gene models, with drh-1 marked --------
##
##   Rscript scripts/DIAG_chrIV_drh1_lod_genemodels.R
##     -> plots/diagnostics/DIAG_chrIV_drh1_lod_genemodels.{pdf,png}
##
## One panel. Above the zero line, the JU1793 allele frequency in each of the
## three sequenced pools -- HT115 control, mig-6 RNAi, pos-1 RNAi -- with 95%
## bands; below it, every protein-coding gene model in the window, packed into
## lanes so nothing overplots. drh-1 is drawn in red with a guide line running
## up through the traces.
##
## WHY FREQUENCY AND NOT LOD. The first version of this figure drew the two
## contrast LOD traces. Both are broad plateaus across the whole window,
## because a very large effect in a bulk-segregant scan carries almost no
## positional resolution, so the trace said nothing the interval bounds did not.
## The frequencies say more: they show the control baseline, how far each RNAi
## condition moves off it, and that the two conditions move by different
## amounts -- which is the actual result.
##
## WHY THIS EXISTS
## JU1793 carries niDf250, the 159 bp drh-1 deletion (Ashe et al. 2013), and
## JU2466 does not, so drh-1 segregates in this cross. The gene sits inside the
## HT115g-POS1g support interval and just outside the HT115g-MIG6g one. The
## figure is here to show that position honestly: against the LOD traces, and
## against the several hundred other genes in the same interval.
##
## WHAT IT DOES NOT SHOW. drh-1's only difference between these parents is the
## deletion; it carries no protein-altering SNV or indel, so it would score
## zero in a coding census. Gene models are drawn from the reference annotation
## and do not show the deletion.
##
## Source annotation and the cross plot data are outside the repository; the
## script names both and stops if either is missing.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(ggtext)
})

GFF <- "/Users/Stefan/UCLA/Genomics_Data/Annotations/c_elegans.PRJNA13758.WS283.csq.gff3.gz"
XP  <- "data/cross_experiments/JU1793-JU2466_export"
OUT <- "plots/diagnostics"
CHR <- "IV"
## 1.3 Mb, chosen so both chrIV peaks (5.75 Mb mig-6, 6.42 Mb pos-1) and drh-1
## (6.61 Mb) are in frame AND a 3 kb gene is still wide enough to draw. At the
## 3.3 Mb span of the support intervals a gene is four pixels and the models are
## meaningless.
WIN <- c(5.6e6, 6.9e6)
FOCUS <- "drh-1"
DRH1 <- c(6607376, 6613353)          # WS283 gene span
NIDF <- c(6607635, 6607793)          # niDf250, Ashe et al. 2013

COL_MIGR <- "#E08214"   # mig-6   (same hues as the cross-contrast figure)
COL_POS  <- "#7C6A9C"   # pos-1
COL_CTRL <- "#4D4D4D"   # HT115 control
COL_GENE <- "#5B6B78"
COL_FOC  <- "#C4302B"   # drh-1
COL_LOF  <- "#1F78A8"   # the other loss-of-function genes in the window

## ---------------------------------------------------------------------------
## Loss-of-function variation between the parents inside the window. PINNED,
## with the queries that produced it, so the figure does not re-scan a 10 GB
## VCF on every render.
##
##   SNV/indel, truncating only (stop gained/lost, frameshift, splice
##   acceptor/donor, start lost), homozygous-discordant between the parents:
##     bcftools query -f '%POS\t%INFO/BCSQ[\t%GT]' -s JU1793,JU2466 \
##       -r IV:5600000-6900000 CeNDR/20231213/bcsq.vcf.gz
##     -> 163 discordant homozygous sites in the window, 2 of them truncating
##
##   structural variants, PASS and discordant, overlapping a CDS:
##     bcftools query ... VCFs/WI.MANTAsv.soft-filter.vcf.gz
##     -> 8 PASS discordant SVs in the window, 1 homozygous CDS hit (drh-1)
##
## NOT included: a breakend call at IV:5,653,829 overlapping skn-1, genotyped
## 0/1 in JU1793. A heterozygous breakend in an inbred strain is not a credible
## loss-of-function call, and no deletion or duplication supports it.
LOF <- data.table(
  gene    = c("drh-1", "bec-1", "cpi-1"),
  kind    = c("159 bp deletion (*niDf250*)", "start lost", "stop gained"),
  carrier = c("JU1793", "JU1793", "JU2466"),
  af      = c(NA, 0.4182, 0.0050))

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
stopifnot(file.exists(GFF), dir.exists(XP))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------------------
## JU1793 allele frequency per pool, and the support intervals
## ---------------------------------------------------------------------------
## afd_<sample> is the JU1793 allele frequency; afd.se_ its standard error.
POOLS <- c("S32_JU1793_JU2466_2_HT115g" = "HT115 control",
           "S33_JU1793_JU2466_2_MIG6g"  = "*mig-6* RNAi",
           "S34_JU1793_JU2466_2_POS1g"  = "*pos-1* RNAi")
COLS <- setNames(c(COL_CTRL, COL_MIGR, COL_POS), POOLS)

af <- imap_dfr(POOLS, function(lab, smp) {
  f <- file.path(XP, "afd", paste0(smp, ".afd.tsv.gz"))
  stopifnot(file.exists(f))
  d <- fread(f, showProgress = FALSE)
  d[chrom == CHR, .(pool = lab, pos = physical.position,
                    f = get(paste0("afd_", smp)),
                    se = get(paste0("afd.se_", smp)))]
}) %>% filter(pos >= WIN[1], pos <= WIN[2], is.finite(f)) %>%
  mutate(pool = factor(pool, levels = POOLS),
         lo = pmax(0, f - 1.96 * se), hi = pmin(1, f + 1.96 * se))
msg("frequency: ", nrow(af), " rows across ", length(POOLS), " pools")

CONTRASTS <- c("HT115g-MIG6g" = "*mig-6* RNAi", "HT115g-POS1g" = "*pos-1* RNAi")
iv <- imap_dfr(CONTRASTS, function(lab, key) {
  f <- list.files(file.path(XP, "intervals"), pattern = key, full.names = TRUE)
  fread(f)[chrom == CHR][, .(pool = lab, lcon, rcon, peak.position, peak.LOD)]
}) %>% mutate(pool = factor(pool, levels = POOLS))
FMAX <- 1
at_drh1 <- af %>% group_by(pool) %>%
  slice_min(abs(pos - mean(DRH1)), n = 1, with_ties = FALSE) %>% ungroup()
msg("JU1793 frequency at drh-1: ",
    paste(sprintf("%s %.2f", gsub("[*]", "", at_drh1$pool), at_drh1$f),
          collapse = " | "))

## ---------------------------------------------------------------------------
## gene models: union exon structure per protein-coding gene
## ---------------------------------------------------------------------------
## length-preserving: str_match gives NA where the key is absent, which
## regmatches(regexpr(...)) does not
att <- function(x, k) stringr::str_match(x, paste0(k, "=([^;]*)"))[, 2]
gcmd <- sprintf("gzcat %s | awk -F'\\t' '$1==\"%s\" && ($3==\"gene\"||$3==\"mRNA\"||$3==\"exon\")'",
                shQuote(GFF), CHR)
g <- fread(cmd = gcmd, header = FALSE, sep = "\t", showProgress = FALSE,
           col.names = c("chrom","src","feat","start","end","score","strand","frame","attr"))

genes <- g[feat == "gene" & start <= WIN[2] & end >= WIN[1]]
genes[, `:=`(wb = att(attr, "Name"), locus = att(attr, "locus"),
             seqname = att(attr, "sequence_name"), biotype = att(attr, "biotype"))]
genes <- genes[biotype == "protein_coding"]
genes[, label := fifelse(!is.na(locus) & nzchar(locus), locus, seqname)]
msg("genes: ", nrow(genes), " protein-coding in IV:",
    format(WIN[1], big.mark = ","), "-", format(WIN[2], big.mark = ","))

mrna <- g[feat == "mRNA"]
mrna[, `:=`(tx = sub("^transcript:", "", att(attr, "ID")),
            gene = sub("^gene:", "", att(attr, "Parent")))]
tx2gene <- setNames(mrna$gene, mrna$tx)
ex <- g[feat == "exon"]
ex[, tx := sub("^transcript:", "", att(attr, "Parent"))]
ex[, wb := tx2gene[tx]]
ex <- ex[wb %in% genes$wb, .(wb, start, end)]
## union of exons across isoforms, so one model per gene
setorder(ex, wb, start)
ex[, grp := cumsum(c(TRUE, start[-1] > cummax(end)[-.N])), by = wb]
exu <- ex[, .(start = min(start), end = max(end)), by = .(wb, grp)]
msg("exons: ", nrow(ex), " -> ", nrow(exu), " merged blocks")

## ---- lane packing: first lane whose last gene ended far enough left --------
setorder(genes, start)
PAD <- diff(WIN) * 0.006          # keep a visible gap between neighbours
lane_end <- numeric(0)
genes[, lane := {
  out <- integer(.N)
  for (i in seq_len(.N)) {
    k <- which(lane_end < start[i] - PAD)[1]
    if (is.na(k)) { lane_end <<- c(lane_end, end[i]); k <- length(lane_end) }
    else lane_end[k] <<- end[i]
    out[i] <- k
  }
  out
}]
NLANE <- max(genes$lane)
msg("packed into ", NLANE, " lanes")

## map lanes into negative y, taking ~46% of the height the traces use
STEP <- (FMAX * 0.46) / NLANE
genes[, y := -lane * STEP]
exu <- merge(exu, genes[, .(wb, y, label)], by = "wb")
FOC <- genes[label == FOCUS]
stopifnot(nrow(FOC) == 1)
lof <- merge(LOF, genes[, .(gene = label, start, end, y, lane)], by = "gene")
stopifnot(nrow(lof) == nrow(LOF))
lof[, txt := paste0("**", gene, "** <span style='font-size:7pt'>", kind,
                    ", ", carrier,
                    fifelse(is.na(af), "", sprintf(" &middot; AF %.3f", af)),
                    "</span>")]
lof[, col := fifelse(gene == FOCUS, COL_FOC, COL_LOF)]
msg("loss-of-function genes in the window: ",
    paste(sprintf("%s (%s, %s)", lof$gene, lof$kind, lof$carrier), collapse = "; "))
msg(FOCUS, " in lane ", FOC$lane, " at IV:", format(FOC$start, big.mark = ","),
    "-", format(FOC$end, big.mark = ","))

EXH <- STEP * 0.36                # exon half-height
mb <- function(x) x / 1e6

## ---------------------------------------------------------------------------
p <- ggplot() +
  ## support intervals, clipped to the window
  geom_rect(data = iv, aes(xmin = mb(pmax(lcon, WIN[1])),
                           xmax = mb(pmin(rcon, WIN[2])), fill = pool),
            ymin = 0, ymax = FMAX * 1.04, alpha = 0.11, colour = NA) +
  ## the drh-1 guide, running the full height
  annotate("rect", xmin = mb(NIDF[1]), xmax = mb(NIDF[2]),
           ymin = -NLANE * STEP - STEP, ymax = FMAX * 1.04,
           fill = COL_FOC, alpha = 0.18, colour = NA) +
  annotate("segment", x = mb(mean(DRH1)), xend = mb(mean(DRH1)),
           y = -NLANE * STEP - STEP, yend = FMAX * 1.04,
           colour = COL_FOC, linewidth = 0.55, linetype = "21") +
  ## 0.5 = no distortion; the control sits well below it here
  annotate("segment", x = mb(WIN[1]), xend = mb(WIN[2]), y = 0.5, yend = 0.5,
           colour = "grey60", linewidth = 0.35, linetype = "31") +
  annotate("text", x = mb(WIN[2]), y = 0.515, hjust = 1, vjust = 0, size = 2.9,
           colour = "grey45", label = "0.5 = no allele-frequency distortion") +
  ## the three pools
  geom_ribbon(data = af, aes(mb(pos), ymin = lo, ymax = hi, fill = pool),
              alpha = 0.22, colour = NA, show.legend = FALSE) +
  geom_line(data = af, aes(mb(pos), f, colour = pool), linewidth = 0.6) +
  geom_point(data = iv, aes(mb(peak.position), 0.5, colour = pool),
             shape = 25, size = 2.1, show.legend = FALSE) +
  ## label for the lower half
  annotate("text", x = mb(WIN[1]) + mb(diff(WIN)) * 0.004,
           y = -NLANE * STEP - STEP * 0.45, hjust = 0, vjust = 0.5,
           size = 3.1, colour = "grey35",
           label = paste0(nrow(genes), " protein-coding gene models, ",
                          NLANE, " lanes")) +
  annotate("segment", x = mb(WIN[1]), xend = mb(WIN[2]), y = 0, yend = 0,
           colour = "grey35", linewidth = 0.4) +
  ## gene models
  geom_segment(data = genes, aes(x = mb(start), xend = mb(end), y = y, yend = y),
               colour = COL_GENE, linewidth = 0.3) +
  geom_rect(data = exu, aes(xmin = mb(start), xmax = mb(end),
                            ymin = y - EXH, ymax = y + EXH),
            fill = COL_GENE, colour = NA) +
  ## the loss-of-function genes on top, drh-1 in red and the rest in blue
  geom_segment(data = lof, aes(x = mb(start), xend = mb(end), y = y, yend = y,
                               colour = I(col)), linewidth = 0.5) +
  geom_rect(data = merge(exu, lof[, .(label = gene, col)], by = "label"),
            aes(xmin = mb(start), xmax = mb(end), ymin = y - EXH * 1.9,
                ymax = y + EXH * 1.9, fill = I(col)), colour = NA) +
  geom_richtext(data = lof, aes(mb(end) + mb(diff(WIN)) * 0.006, y, label = txt,
                                colour = I(col)),
                hjust = 0, size = 3.1, fill = alpha("white", 0.86),
                label.color = NA,
                label.padding = grid::unit(c(1.2, 2, 1.2, 2), "pt")) +
  scale_colour_manual(values = COLS, name = NULL) +
  scale_fill_manual(values = COLS, guide = "none") +
  scale_x_continuous(name = "Chromosome IV position (Mb)",
                     breaks = seq(5.6, 6.9, 0.1),
                     expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(
    name = "JU1793 allele frequency",
    breaks = seq(0, 1, 0.25),
    limits = c(-NLANE * STEP - STEP * 0.8, FMAX * 1.06),
    expand = expansion(mult = 0)) +
  coord_cartesian(xlim = mb(WIN), expand = TRUE) +
  labs(title = "**The chromosome IV locus, and the genes under it**",
       subtitle = paste0(
         "JU1793 &times; JU2466 F2 pools. The trace is the JU1793 allele ",
         "frequency with a 95% band, so a rise means selection for the ",
         "JU1793 allele. Shaded bands are the two contrasts' support ",
         "intervals, triangles their chrIV peaks.<br>",
         "The JU1793 allele sits near 0.25 in the untreated control across ",
         "this whole window, rises to about 0.78 under *mig-6* RNAi and about ",
         "0.61 under *pos-1* &mdash; enriched under **both** knockdowns. That ",
         "is the general-RNAi-response signature, and unlike the chrX locus.<br>",
         "Below the line, all ", nrow(genes), " protein-coding gene models in ",
         "the window. Every gene carrying loss-of-function variation between ",
         "the parents is marked: ",
         "<span style='color:", COL_FOC, "'>**drh-1**</span> in red (the red ",
         "stripe is *niDf250*), the others in ",
         "<span style='color:", COL_LOF, "'>**blue**</span>. Truncating SNVs ",
         "and indels plus PASS structural variants overlapping a CDS; ",
         "missense is not counted.<br>The traces are flat across megabases, ",
         "so position alone singles out none of them.")) +
  theme_classic(base_size = 11.5) +
  theme(axis.line.y = element_line(linewidth = 0.3),
        axis.line.x = element_blank(),
        axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_markdown(size = 13),
        plot.subtitle = element_markdown(size = 8.4, colour = "grey30",
                                         lineheight = 1.35),
        plot.title.position = "plot",
        legend.position = c(0.012, 0.62),
        legend.justification = c(0, 1),
        legend.text = element_markdown(size = 9.5),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key.size = grid::unit(11, "pt"))

ggsave(file.path(OUT, "DIAG_chrIV_drh1_lod_genemodels.pdf"), p,
       width = 18, height = 8.0, device = cairo_pdf)
ggsave(file.path(OUT, "DIAG_chrIV_drh1_lod_genemodels.png"), p,
       width = 18, height = 8.0, dpi = 300, bg = "white")
msg("wrote DIAG_chrIV_drh1_lod_genemodels.{pdf,png}")
