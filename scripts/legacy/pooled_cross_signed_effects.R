## pooled_cross_signed_effects.R -------------------------------------------
## Signed versions of the pooled GWAS and cross LOD scans, so the *direction*
## of each effect can be compared between the association and linkage
## experiments rather than only its strength.
##
##   Rscript scripts/pooled_cross_signed_effects.R
##
## ---------------------------------------------------------------------------
## Sign convention
## ---------------------------------------------------------------------------
## Each export phases allele-frequency deviation onto parent 1: N2 for the
## N2 x XZ1516 cross, JU1793 for the JU1793 x JU2466 cross. Its `contrast.beta`
## is (left condition - right condition) on that parent's allele frequency.
##
## The polarity was fixed empirically, not assumed. At the JU-cross chromosome
## III peak the JU1793 allele frequency rises from 0.33 in the HT115 control to
## 0.78 under pos-1 RNAi -- strongly enriched, i.e. the JU1793 allele is the
## resistant one -- and `contrast.beta` there is negative (-0.51). Hence
##
##     signed LOD = -sign(contrast.beta) * LOD
##
## and a POSITIVE signed LOD means the parent-1 allele (N2, or JU1793) is at
## higher frequency in the right-hand condition of the contrast name. For an
## "HT115 vs <gene>" contrast the right-hand condition is the RNAi treatment,
## so positive = the parent-1 allele confers resistance to that RNAi. This
## matches what the user asked for: N2 and JU1793 positive, XZ1516 and JU2466
## negative.
##
## Independently, this agrees with the pooled GWAS: JU1793 carries the
## alternate resistance allele at the chromosome III pooled peak, and the JU
## cross says its chromosome III allele is the resistant one.
##
## For the pooled GWAS, `beta` is the effect of allele1 = the ALT allele on the
## RNAi-minus-control trait, so positive beta = ALT confers resistance. Because
## the CeNDR reference allele IS the N2 allele, re-polarising to the reference
## gives a quantity directly comparable to the N2 x XZ1516 signed LOD:
##
##     signed GWAS = -sign(beta) * neglog10p     (positive = N2/REF resistant)
##
## FIG13 signs every cross by its own parent 1, as requested. FIG14 additionally
## re-polarises the JU cross onto the reference allele -- flipping the sign at
## every marker where JU1793 carries ALT -- so that all eight tracks share one
## axis meaning and the GWAS and both crosses can be read against each other.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(ggtext)
})

b   <- readRDS("data/pooled_cross_intersection/bundle.rds")
OUT <- "plots/pooled_cross_intersection"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

CHROMS <- c("I", "II", "III", "IV", "V", "X")
PAL <- c(pooled = "#2E4057", nxz = "#0B7A75", juju = "#D57A00", hi = "#C4302B")
CROSS_LAB <- c(N2xXZ1516 = "N2 × XZ1516", JU1793xJU2466 = "JU1793 × JU2466")
CROSS_P1  <- c(N2xXZ1516 = "N2", JU1793xJU2466 = "JU1793")
CROSS_P2  <- c(N2xXZ1516 = "XZ1516", JU1793xJU2466 = "JU2466")
COND_LAB  <- c(`ht115 vs mig6` = "HT115 vs *mig-6*",
               `ht115 vs pos1` = "HT115 vs *pos-1*",
               `pos1 vs mig6`  = "*pos-1* vs *mig-6*",
               `mig6 vs pos1`  = "*mig-6* vs *pos-1*")
KEY6 <- c("N2xXZ1516 | ht115 vs mig6", "N2xXZ1516 | ht115 vs pos1",
          "N2xXZ1516 | pos1 vs mig6",
          "JU1793xJU2466 | ht115 vs mig6", "JU1793xJU2466 | ht115 vs pos1",
          "JU1793xJU2466 | mig6 vs pos1")

LOCI <- tribble(
  ~locus,     ~chrom, ~pos,
  "I:1.9",    "I",     1876106,
  "III:12.4", "III",  12353680,
  "III:13.7", "III",  13679445,
  "IV:9.0",   "IV",    8975323,
  "V:13.8",   "V",    13805216,
  "V:14.6",   "V",    14647434,
  "X:6.9",    "X",     6904721,
  "X:9.7",    "X",     9705058) %>%
  mutate(chrom = factor(chrom, levels = CHROMS), pos.mb = pos / 1e6)

fct_chrom <- function(x) factor(as.character(x), levels = CHROMS)

## Peak-preserving thinning. An 11-inch genome panel has ~3300 usable pixels,
## so plotting 500k markers writes ~150 vertices per pixel: slow to render and
## a 13 MB vector file. Binning and keeping the most extreme value per bin
## preserves every peak, its position and its sign, while cutting vertices
## ~100x. `y` is the column whose absolute value defines "most extreme".
thin_extreme <- function(d, y, bp = 5e3, by = NULL) {
  g <- c(by, "chrom", ".bin")
  d %>% mutate(.bin = floor(pos / bp)) %>%
    group_by(across(all_of(g))) %>%
    slice_max(abs(.data[[y]]), n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(-.bin) %>% arrange(across(all_of(c(by, "chrom", "pos"))))
}
wrap_txt  <- function(x, n = 130) paste(strwrap(x, width = n), collapse = "\n")
wrap_md   <- function(x, n = 130) paste(strwrap(x, width = n), collapse = "<br>")

theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = base_size),
          strip.placement = "outside",
          panel.spacing.x = grid::unit(4, "pt"),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          plot.title.position = "plot",
          legend.key.size = grid::unit(9, "pt"))
}

## WBcel235/ce11 chromosome lengths. facet_grid(space = "free_x") sizes each
## panel by its own data range, so panels built from different marker sets (or
## from just the eight loci, as in loci_key) came out different widths and did
## not line up vertically. Forcing every genome panel to span the full
## chromosome via geom_blank fixes the alignment that this whole figure set
## depends on.
CHROM_LEN <- c(I = 15072434, II = 15279421, III = 13783801,
               IV = 17493829, V = 20924180, X = 17718942)
chrom_span <- tibble::tibble(chrom = factor(names(CHROM_LEN), levels = CHROMS),
                             pos.mb = unname(CHROM_LEN) / 1e6) %>%
  dplyr::bind_rows(tibble::tibble(chrom = factor(names(CHROM_LEN), levels = CHROMS),
                                  pos.mb = 0))

genome_facets <- function() list(
  geom_blank(data = chrom_span, aes(x = pos.mb, y = 0), inherit.aes = FALSE),
  facet_grid(. ~ chrom, scales = "free_x", space = "free_x", switch = "x"),
  scale_x_continuous(breaks = seq(0, 25, 5), expand = expansion(mult = 0.02),
                     guide = guide_axis(check.overlap = TRUE)))
axis_ctl <- function(strip = TRUE, x = TRUE) theme(
  strip.text = if (strip) NULL else element_blank(),
  axis.text.x = if (x) element_text(size = 6) else element_blank(),
  axis.ticks.x = if (x) element_line(linewidth = 0.3) else element_blank())
guide_lines <- function() geom_vline(data = LOCI, aes(xintercept = pos.mb),
  linetype = "dotted", linewidth = 0.25, colour = "grey55")

## slim labelled key row, so the dotted guides identify themselves
loci_key <- function() {
  ggplot(LOCI, aes(pos.mb, 0)) +
    geom_segment(aes(xend = pos.mb, y = 0, yend = 1), linetype = "dotted",
                 linewidth = 0.25, colour = "grey55") +
    geom_text(aes(y = 1.15, label = locus), angle = 90, hjust = 0, vjust = 0.5,
              size = 2.1, colour = "grey25") +
    genome_facets() +
    coord_cartesian(ylim = c(0, 4.4), clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_pub(9) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(), strip.text = element_blank(),
          plot.margin = margin(t = 2, r = 10, b = 0, l = 6))
}

save_fig <- function(p, name, width, height) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = width, height = height,
         device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = width, height = height,
         dpi = 300, bg = "white", limitsize = FALSE)
  cat("  wrote ", name, ".{pdf,png}\n", sep = "")
}

## ===========================================================================
## 1. signed cross scans, each on its own parent-1 polarity
## ===========================================================================
signed <- imap_dfr(b$scans[KEY6], function(s, k) {
  m <- b$scan_meta %>% filter(key == k)
  s %>% filter(is.finite(LOD), chrom %in% CHROMS) %>%
    transmute(key = k, cross = m$cross, label = m$label,
              chrom = fct_chrom(chrom), pos = physical.position,
              pos.mb = physical.position / 1e6,
              LOD, beta = contrast.beta,
              ## empirical polarity: negative beta = parent-1 allele enriched
              ## in the right-hand condition
              sLOD = -sign(contrast.beta) * LOD)
})

## ===========================================================================
## 2. re-polarise the JU cross onto the reference (N2) allele
## ===========================================================================
## The JU export is phased to JU1793. Where JU1793 carries ALT, "parent-1
## allele" and "reference allele" point opposite ways, so the sign must flip
## for the two crosses to be comparable on one axis.
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
pgt <- fread("data/pooled_cross_intersection/parents_gt.tsv.gz",
             col.names = c("chrom", "pos", "ref", "alt",
                           "N2", "XZ1516", "JU1793", "JU2466"))
pgt <- pgt[chrom %in% CHROMS & JU1793 %in% c("0/0", "1/1"),
           .(chrom, pos, ju_alt = as.integer(JU1793 == "1/1"))]

ju <- signed %>% filter(cross == "JU1793xJU2466") %>%
  mutate(chrom = as.character(chrom)) %>%
  left_join(pgt, by = c("chrom", "pos"))

hit <- mean(!is.na(ju$ju_alt))
msg("JU-cross markers matched to a called JU1793 genotype: ",
    sprintf("%.1f%%", 100 * hit), " (", sum(!is.na(ju$ju_alt)), " of ", nrow(ju), ")")

## flip where JU1793 carries ALT; drop markers we cannot orient
ju_ref <- ju %>% filter(!is.na(ju_alt)) %>%
  mutate(chrom = fct_chrom(chrom),
         sLOD_ref = ifelse(ju_alt == 1, -sLOD, sLOD))

## N2 is the reference by definition, so the N2 cross needs no flip
nxz_ref <- signed %>% filter(cross == "N2xXZ1516") %>% mutate(sLOD_ref = sLOD)
signed_ref <- bind_rows(nxz_ref, ju_ref %>% select(names(nxz_ref)))

## ===========================================================================
## 3. signed pooled GWAS
## ===========================================================================
gw <- b$gwas %>% filter(chr %in% CHROMS) %>%
  transmute(gene, chrom = fct_chrom(chr), pos.mb = ps / 1e6, ps,
            neglog10p, beta,
            ## positive = REF (= N2) allele raises resistance
            s_gwas = -sign(beta) * neglog10p)

## ===========================================================================
## tracks
## ===========================================================================
gw_track <- function(g, ycol, ylab, strip = FALSE, xax = FALSE) {
  d  <- gw %>% filter(gene == g) %>% mutate(pos = ps) %>% thin_extreme(ycol, bp = 1e4)
  pk <- b$gwas_sig %>% filter(gene == g) %>%
    transmute(chrom = fct_chrom(chr), pos.mb = ps / 1e6,
              y = -sign(beta) * neglog10p, beta)
  ggplot(d, aes(pos.mb, .data[[ycol]])) +
    guide_lines() +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey30") +
    geom_hline(yintercept = c(-1, 1) * b$GWAS_BF, linetype = "dashed",
               linewidth = 0.28, colour = "grey50") +
    geom_point(size = 0.18, colour = PAL[["pooled"]], alpha = 0.5, stroke = 0) +
    { if (nrow(pk)) geom_point(data = pk, aes(y = y), colour = PAL[["hi"]], size = 1.6) } +
    { if (nrow(pk)) geom_text(data = pk, aes(y = y, label = sprintf("%+.2f", beta)),
                              vjust = ifelse(pk$y > 0, -0.8, 1.7), size = 2.4,
                              colour = PAL[["hi"]], fontface = "bold") } +
    genome_facets() +
    scale_y_continuous(expand = expansion(mult = 0.14)) +
    labs(x = NULL, y = ylab,
         title = paste0("Pooled GWAS · *", g, "* RNAi")) +
    theme_pub(9) +
    theme(plot.title = element_markdown(face = "bold", size = 9,
                                        colour = PAL[["pooled"]]),
          axis.title.y = element_markdown(size = 7)) +
    axis_ctl(strip, xax)
}

cross_track <- function(dat, key, ycol, strip = FALSE, xax = FALSE) {
  d <- dat %>% filter(key == !!key) %>% thin_extreme(ycol)
  m <- b$scan_meta %>% filter(key == !!key)
  col <- if (m$cross == "N2xXZ1516") PAL[["nxz"]] else PAL[["juju"]]
  ggplot(d, aes(pos.mb, .data[[ycol]])) +
    guide_lines() +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey30") +
    geom_hline(yintercept = c(-1, 1) * b$CROSS_THR, linetype = "dashed",
               linewidth = 0.28, colour = "grey50") +
    ## one ribbon from zero, single group: a sign-keyed geom_area would build a
    ## polygon per sign flip and never finish on ~500k markers
    geom_ribbon(aes(ymin = pmin(.data[[ycol]], 0), ymax = pmax(.data[[ycol]], 0)),
                fill = col, alpha = 0.3, linewidth = 0) +
    geom_line(linewidth = 0.3, colour = "grey15") +
    genome_facets() +
    scale_y_continuous(expand = expansion(mult = 0.1)) +
    labs(x = NULL, y = "signed LOD",
         title = paste0(CROSS_LAB[[m$cross]], " · ", COND_LAB[[m$label]])) +
    theme_pub(9) +
    theme(plot.title = element_markdown(face = "bold", size = 9, colour = col),
          axis.title.y = element_text(size = 7)) +
    axis_ctl(strip, xax)
}

## ===========================================================================
## FIG 13  signed, each cross on its own parent-1 polarity
## ===========================================================================
fig13 <- if (nzchar(Sys.getenv("SKIP_BIG"))) NULL else
  loci_key() / gw_track("mig-6", "s_gwas", "signed<br>−log<sub>10</sub>*p*") /
  gw_track("pos-1", "s_gwas", "signed<br>−log<sub>10</sub>*p*") /
  cross_track(signed, KEY6[1], "sLOD") /
  cross_track(signed, KEY6[2], "sLOD") /
  cross_track(signed, KEY6[3], "sLOD") /
  cross_track(signed, KEY6[4], "sLOD") /
  cross_track(signed, KEY6[5], "sLOD") /
  cross_track(signed, KEY6[6], "sLOD", strip = TRUE, xax = TRUE) +
  plot_layout(heights = c(0.42, rep(1, 8))) +
  plot_annotation(
    title = "Signed scans: which parent's allele carries the effect",
    subtitle = wrap_md(paste0(
      "**Cross tracks** are signed −sign(contrast.beta) × LOD, so **up = the N2 or ",
      "JU1793 allele** is at higher frequency in the right-hand condition of the ",
      "contrast (for an HT115-vs-RNAi contrast, that means that parent's allele is the ",
      "resistant one) and **down = XZ1516 or JU2466**. **GWAS tracks** are signed ",
      "−sign(beta) × −log<sub>10</sub>*p*, so up = the reference (N2) allele raises ",
      "resistance; the number by each peak is its beta. Because N2 *is* the reference, the ",
      "three **N2 × XZ1516** tracks are already on the same polarity as the GWAS tracks and ",
      "can be read directly against them. The **JU1793 × JU2466** tracks are polarised on ",
      "JU1793, not the reference, so their signs are not comparable to the GWAS ",
      "genome-wide -- that comparison is only defined at a shared variant, which is what ",
      "FIG15 shows. Dashed lines are the genome-wide thresholds; the dotted verticals, ",
      "labelled in the key row at the top, are the eight loci that form the columns of ",
      "FIG15."), 132),
    caption = "Position (Mb)",
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_markdown(size = 8, colour = "grey25"),
                  plot.caption = element_text(hjust = 0.5, size = 9)))

SKIP_BIG <- nzchar(Sys.getenv("SKIP_BIG"))
if (!SKIP_BIG) save_fig(fig13, "FIG13_signed_scans_by_parent", 11, 14)

## ===========================================================================
## FIG 14 -- deliberately not produced
## ===========================================================================
## An earlier version re-polarised the JU cross onto the reference allele
## marker by marker, flipping signed LOD wherever JU1793 carried ALT, so that
## all eight tracks would share one axis. That figure is WRONG and has been
## removed rather than fixed.
##
## A cross LOD at a marker is a property of the surrounding haplotype block, not
## of that single SNP's REF/ALT identity. Adjacent markers alternate between
## JU1793-REF and JU1793-ALT constantly, so per-marker flipping makes the trace
## whipsaw between +LOD and -LOD thousands of times and the panel collapses into
## a solid block. The apparent "comparability" was an artefact.
##
## The comparison of GWAS effect direction against cross effect direction is
## only well defined AT A SHARED VARIANT: does the GWAS say the ALT allele is
## resistant, and does the cross say the ALT-carrying parent's allele is
## resistant? That is a per-locus question, and FIG15 answers it -- marking the
## cells where the cross actually segregates the variant, which is the only
## place the question has an answer.
##
## Note that no flipping is needed for N2 x XZ1516 at all: N2 IS the reference,
## so those three tracks in FIG13 are already on reference polarity and can be
## read directly against the GWAS tracks above them.
## ===========================================================================

## ===========================================================================
## FIG 15  effect direction at the named loci, GWAS vs cross, side by side
## ===========================================================================
at_locus <- function(chrom, pos) {
  cs <- signed_ref %>% filter(as.character(chrom) == !!chrom) %>%
    group_by(key, cross, label) %>%
    slice_min(abs(pos - !!pos), n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(track = paste0(CROSS_LAB[cross], " · ", label),
              value = sLOD_ref, kind = "cross", thr = b$CROSS_THR)
  gs <- gw %>% filter(as.character(chrom) == !!chrom, abs(ps - !!pos) <= 5e4) %>%
    group_by(gene) %>% slice_max(neglog10p, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(track = paste0("pooled GWAS · ", gene),
              value = s_gwas, kind = "gwas", thr = b$GWAS_BF)
  bind_rows(gs, cs)
}

dir_dat <- LOCI %>% rowwise() %>%
  mutate(d = list(at_locus(as.character(chrom), pos))) %>%
  ungroup() %>% select(locus, d) %>% unnest(d) %>%
  mutate(rel = value / thr,
         locus = factor(locus, levels = LOCI$locus),
         track = factor(track, levels = c(
           paste0("pooled GWAS · ", c("mig-6", "pos-1")),
           paste0(CROSS_LAB[["N2xXZ1516"]], " · ",
                  c("ht115 vs mig6", "ht115 vs pos1", "pos1 vs mig6")),
           paste0(CROSS_LAB[["JU1793xJU2466"]], " · ",
                  c("ht115 vs mig6", "ht115 vs pos1", "mig6 vs pos1")))))

## Cross LODs reach 200x their threshold while the pooled GWAS barely clears 1x,
## so a linear fill on that range renders every GWAS cell white and hides the
## direction this panel exists to show. Saturate at 6x on a sqrt scale: the sign
## stays legible in both assays and the printed number carries the magnitude.
cap <- 6
sat <- function(x) sign(x) * sqrt(pmin(abs(x), cap) / cap)

## Which cross rows can actually test a given pooled peak variant? Only where
## the two parents of that cross carry different alleles there. Anywhere else a
## sign agreement or disagreement with the GWAS is meaningless.
seg <- b$parent_gt_qtl %>%
  transmute(chrom, ps,
            N2xXZ1516     = N2     != XZ1516,
            JU1793xJU2466 = JU1793 != JU2466) %>%
  pivot_longer(c(N2xXZ1516, JU1793xJU2466), names_to = "cross",
               values_to = "segregating") %>%
  left_join(LOCI %>% transmute(locus, chrom = as.character(chrom), pos),
            by = c("chrom" = "chrom", "ps" = "pos")) %>%
  filter(!is.na(locus)) %>%
  transmute(locus, cross, segregating)

dir_dat <- dir_dat %>%
  mutate(cross = case_when(grepl("^N2", track) ~ "N2xXZ1516",
                           grepl("^JU", track) ~ "JU1793xJU2466",
                           TRUE ~ NA_character_)) %>%
  left_join(seg, by = c("locus", "cross"))

fig15 <- ggplot(dir_dat, aes(locus, fct_rev(track), fill = sat(rel))) +
  geom_tile(colour = "white", linewidth = 0.9) +
  ## a ring marks a cell where the cross segregates the pooled peak variant,
  ## i.e. where agreeing or disagreeing with the GWAS actually means something
  geom_point(data = dir_dat %>% filter(segregating),
             shape = 21, size = 7.5, stroke = 0.7,
             fill = NA, colour = "grey15") +
  geom_text(aes(label = ifelse(abs(rel) < 1, "·",
                               sprintf("%+.0f", pmin(abs(rel), 999) * sign(rel))),
                colour = abs(sat(rel)) > 0.72),
            size = 2.7, fontface = "bold") +
  annotate("segment", x = 0.5, xend = nrow(LOCI) + 0.5, y = 6.5, yend = 6.5,
           linewidth = 0.6, colour = "grey30") +
  scale_fill_gradient2(low = "#B2182B", mid = "grey96", high = "#2166AC",
                       midpoint = 0, limits = c(-1, 1),
                       breaks = sat(c(-cap, -1, 0, 1, cap)),
                       labels = c("≤−6×", "−1×", "0", "+1×", "≥+6×"),
                       name = "signed signal /\nthreshold") +
  scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "white"), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "Effect direction at each locus, on one polarity",
       subtitle = wrap_md(paste0(
         "Signal divided by that assay's own genome-wide threshold, signed so that ",
         "**blue = the reference (N2) allele raises resistance** and **red = the ",
         "alternate allele does**. Fill saturates at 6x so direction stays legible in ",
         "both assays; the printed number is the true multiple, and a dot marks a cell ",
         "below its own threshold where the sign is meaningless. A **ring** marks the ",
         "cells where that cross actually segregates the pooled peak variant -- only ",
         "there does agreeing or disagreeing with the GWAS mean anything. Elsewhere both ",
         "parents carry the same allele and the cross is necessarily reading a different, ",
         "linked variant. Pooled GWAS above the rule, crosses below."), 118)) +
  theme_pub(9) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(size = 8, angle = 30, hjust = 0),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_markdown(size = 8, colour = "grey25"),
        legend.key.height = grid::unit(15, "pt"))

save_fig(fig15, "FIG15_effect_direction_matrix", 10, 5.4)

## ===========================================================================
## numbers
## ===========================================================================
write_tsv(dir_dat %>% select(locus, track, kind, value, thr, rel),
          file.path(OUT, "TABLE_signed_effect_directions.tsv"))

cat("\n== signed effect at each named locus (value / own threshold) ==\n")
print(as.data.frame(dir_dat %>%
        transmute(locus, track, signed_value = round(value, 1),
                  x_threshold = round(rel, 1)) %>%
        pivot_wider(names_from = locus, values_from = c(signed_value, x_threshold)) %>%
        select(track, starts_with("x_threshold"))), row.names = FALSE)
cat("\nJU-cross marker orientation rate: ", sprintf("%.1f%%", 100 * hit), "\n", sep = "")
