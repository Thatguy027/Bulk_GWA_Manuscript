## pooled_cross_intersection_figures.R -------------------------------------
## Publication figures for the pooled-GWAS -> outlier -> cross story.
##
##   Rscript scripts/pooled_cross_intersection_prep.R      # once, builds the cache
##   Rscript scripts/pooled_cross_intersection_figures.R
##
## Every figure is written to plots/pooled_cross_intersection/ as both pdf
## (cairo, for the Greek/minus glyphs) and png at 300 dpi.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggtext)
  library(ggrepel)
})

b   <- readRDS("data/pooled_cross_intersection/bundle.rds")
OUT <- "plots/pooled_cross_intersection"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ===========================================================================
## Shared vocabulary
## ===========================================================================
CHROMS <- c("I", "II", "III", "IV", "V", "X")

PAL <- c(pooled  = "#2E4057",   # pooled GWAS
         nxz     = "#0B7A75",   # N2 x XZ1516
         juju    = "#D57A00",   # JU1793 x JU2466
         hi      = "#C4302B",   # peaks / support intervals
         gene    = "#4D4D4D",   # gene annotation (kept off the contrast hues)
         uptake  = "#762A83")

CROSS_COL <- c(`N2xXZ1516` = PAL[["nxz"]], `JU1793xJU2466` = PAL[["juju"]])
CROSS_LAB <- c(`N2xXZ1516`     = "N2 × XZ1516",
               `JU1793xJU2466` = "JU1793 × JU2466")

## The RNAi-vs-RNAi contrast is written in the opposite order in the two
## exports; both mean the same comparison, so they share a colour.
COND_LAB <- c(`ht115 vs mig6` = "HT115 vs *mig-6*",
              `ht115 vs pos1` = "HT115 vs *pos-1*",
              `pos1 vs mig6`  = "*pos-1* vs *mig-6*",
              `mig6 vs pos1`  = "*mig-6* vs *pos-1*")
COND_COL <- c(`ht115 vs mig6` = "#B2182B",
              `ht115 vs pos1` = "#2166AC",
              `pos1 vs mig6`  = "#1A9850",
              `mig6 vs pos1`  = "#1A9850")

## The six contrasts this story rests on.
KEY6 <- c("N2xXZ1516 | ht115 vs mig6", "N2xXZ1516 | ht115 vs pos1",
          "N2xXZ1516 | pos1 vs mig6",
          "JU1793xJU2466 | ht115 vs mig6", "JU1793xJU2466 | ht115 vs pos1",
          "JU1793xJU2466 | mig6 vs pos1")

## Loci the narrative names. `role` drives the colour coding in figs 8 and 10.
LOCI <- tribble(
  ~locus,     ~chrom, ~pos,      ~gene,     ~role,
  "I:1.9",    "I",     1876106,  NA,        "mig-6-specific",
  "III:12.4", "III",  12353680,  NA,        "pooled pos-1 QTL",
  "III:13.7", "III",  13679445,  "sid-2",   "shared uptake",
  "IV:9.0",   "IV",    8975323,  NA,        "shared uptake",
  "V:13.8",   "V",    13805216,  NA,        "mig-6-specific",
  "V:14.6",   "V",    14647434,  NA,        "pooled mig-6 QTL",
  "X:6.9",    "X",     6904721,  NA,        "mig-6-specific",
  "X:9.7",    "X",     9705058,  NA,        "shared uptake") %>%
  mutate(label = ifelse(is.na(gene), locus, paste0(locus, " (", gene, ")")))

ROLE_LEV <- c("pooled mig-6 QTL", "pooled pos-1 QTL", "shared uptake", "mig-6-specific")
ROLE_COL <- c(`pooled mig-6 QTL` = "#2E4057",
              `pooled pos-1 QTL` = "#7C6A9C",
              `shared uptake`    = "#8C8C8C",
              `mig-6-specific`   = "#C4302B")

## -- text helpers -----------------------------------------------------------
## Long subtitles overflow a patchwork annotation; wrap them explicitly.
wrap_txt <- function(x, n = 120) paste(strwrap(x, width = n), collapse = "\n")
wrap_md  <- function(x, n = 120) paste(strwrap(x, width = n), collapse = "<br>")

## -- theme ------------------------------------------------------------------
theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      strip.background    = element_blank(),
      strip.text          = element_text(face = "bold", size = base_size),
      strip.placement     = "outside",
      panel.spacing.x     = grid::unit(4, "pt"),
      axis.line           = element_line(linewidth = 0.3),
      axis.ticks          = element_line(linewidth = 0.3),
      plot.title          = element_text(face = "bold", size = base_size + 0.5),
      plot.subtitle       = element_text(size = base_size - 1, colour = "grey30"),
      plot.title.position = "plot",
      legend.key.size     = grid::unit(9, "pt"),
      legend.background   = element_blank())
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

genome_facets <- function() {
  list(
    geom_blank(data = chrom_span, aes(x = pos.mb, y = 0), inherit.aes = FALSE),
    facet_grid(. ~ chrom, scales = "free_x", space = "free_x", switch = "x"),
    scale_x_continuous(breaks = seq(0, 25, 5),
                       expand = expansion(mult = 0.02),
                       guide = guide_axis(check.overlap = TRUE)))
}

## Strip/axis suppression must be added AFTER theme_pub(), which resets strip.text.
axis_ctl <- function(strip = TRUE, x = TRUE) {
  theme(strip.text  = if (strip) NULL else element_blank(),
        axis.text.x = if (x) element_text(size = 6) else element_blank(),
        axis.ticks.x = if (x) element_line(linewidth = 0.3) else element_blank())
}

fct_chrom <- function(x) factor(as.character(x), levels = CHROMS)

save_fig <- function(p, name, width, height) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = width, height = height,
         device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = width, height = height,
         dpi = 300, bg = "white", limitsize = FALSE)
  cat("  wrote ", name, ".{pdf,png}  (", width, " x ", height, " in)\n", sep = "")
}

## -- tidy scan / interval frames -------------------------------------------
scan_long <- function(keys) {
  imap_dfr(b$scans[keys], function(s, k) {
    s %>% filter(is.finite(LOD)) %>%
      transmute(key = k, chrom = fct_chrom(chrom),
                pos.mb = physical.position / 1e6, LOD, beta = contrast.beta)
  }) %>%
    left_join(b$scan_meta %>% select(key, cross, label), by = "key") %>%
    mutate(cross = factor(cross, levels = names(CROSS_LAB)),
           contrast = factor(label, levels = names(COND_LAB)))
}

band_for <- function(sl, iv) {
  iv <- iv %>% filter(significant)
  pmap_dfr(iv %>% select(key, chrom, lcon, rcon), function(key, chrom, lcon, rcon) {
    sl %>% filter(key == !!key, as.character(chrom) == !!chrom,
                  pos.mb >= lcon / 1e6, pos.mb <= rcon / 1e6)
  })
}

gwas_long <- b$gwas %>%
  transmute(gene, chrom = fct_chrom(chr), pos.mb = ps / 1e6, ps, neglog10p, beta, af)

## vertical guides at the named loci, for whole-genome panels
loci_guides <- LOCI %>% mutate(chrom = fct_chrom(chrom), pos.mb = pos / 1e6)

guide_lines <- function() {
  geom_vline(data = loci_guides, aes(xintercept = pos.mb),
             linetype = "dotted", linewidth = 0.25, colour = "grey55")
}

## The dotted guides are useless unless the reader can see what they mark, so
## every stacked genome panel gets this slim key row above it, on the same
## chromosome facets. Labels are rotated because two of the loci sit 0.8 Mb
## apart on chromosome V and would collide horizontally.
loci_key <- function(base_size = 9) {
  ggplot(loci_guides, aes(pos.mb, 0)) +
    geom_segment(aes(xend = pos.mb, y = 0, yend = 1), linetype = "dotted",
                 linewidth = 0.25, colour = "grey55") +
    geom_text(aes(y = 1.15, label = label), angle = 90, hjust = 0, vjust = 0.5,
              size = 2.1, colour = "grey25") +
    genome_facets() +
    coord_cartesian(ylim = c(0, 5.6), clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_pub(base_size) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(), strip.text = element_blank(),
          plot.margin = margin(t = 2, r = 10, b = 0, l = 6))
}

cat("figures ->", OUT, "\n")

## ===========================================================================
## FIG 1  Master intersection: pooled GWAS above, both crosses below
## ===========================================================================
sl6 <- scan_long(KEY6)
bd6 <- band_for(sl6, b$ivs %>% filter(key %in% KEY6))

gw_track <- function(g, show_strip = FALSE, show_x = FALSE) {
  d  <- gwas_long %>% filter(gene == g)
  pk <- b$gwas_sig %>% filter(gene == g) %>%
    transmute(chrom = fct_chrom(chr), pos.mb = ps / 1e6, neglog10p)
  ggplot(d, aes(pos.mb, neglog10p)) +
    guide_lines() +
    geom_point(size = 0.18, colour = PAL[["pooled"]], alpha = 0.5, stroke = 0) +
    geom_hline(yintercept = b$GWAS_BF, linetype = "dashed",
               linewidth = 0.3, colour = "grey45") +
    { if (nrow(pk)) geom_point(data = pk, colour = PAL[["hi"]], size = 1.5) } +
    { if (nrow(pk)) geom_text(data = pk, aes(label = sprintf("%.2f", neglog10p)),
                              vjust = -0.7, size = 2.5,
                              colour = PAL[["hi"]], fontface = "bold") } +
    genome_facets() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(y = "−log<sub>10</sub>*p*", x = NULL,
         title = paste0("Pooled GWAS · *", g, "* RNAi (93 wild isotypes)")) +
    theme_pub(9) +
    theme(plot.title = element_markdown(face = "bold", size = 9.5),
          axis.title.y = element_markdown(size = 7.5)) +
    axis_ctl(show_strip, show_x)
}

cross_track <- function(key, show_strip = FALSE, show_x = FALSE) {
  d  <- sl6 %>% filter(key == !!key)
  bn <- bd6 %>% filter(key == !!key)
  iv <- b$ivs %>% filter(key == !!key, significant) %>%
    transmute(chrom = fct_chrom(chrom), pos.mb = peak.position / 1e6, peak.LOD)
  m  <- b$scan_meta %>% filter(key == !!key)
  ggplot(d, aes(pos.mb, LOD)) +
    guide_lines() +
    { if (nrow(bn)) geom_ribbon(data = bn, aes(ymin = 0, ymax = LOD),
                                fill = PAL[["hi"]], alpha = 0.5) } +
    geom_line(linewidth = 0.3, colour = "grey15") +
    geom_hline(yintercept = b$CROSS_THR, linetype = "dashed",
               linewidth = 0.3, colour = "grey45") +
    { if (nrow(iv)) geom_point(data = iv, aes(y = peak.LOD),
                               colour = PAL[["hi"]], size = 0.9) } +
    genome_facets() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(y = "LOD", x = NULL,
         title = paste0(CROSS_LAB[[m$cross]], " · ", COND_LAB[[m$label]])) +
    theme_pub(9) +
    theme(plot.title = element_markdown(face = "bold", size = 9,
                                        colour = CROSS_COL[[m$cross]]),
          axis.title.y = element_text(size = 7.5)) +
    axis_ctl(show_strip, show_x)
}

fig1 <- loci_key() / gw_track("mig-6") / gw_track("pos-1") /
  cross_track(KEY6[1]) / cross_track(KEY6[2]) / cross_track(KEY6[3]) /
  cross_track(KEY6[4]) / cross_track(KEY6[5]) /
  cross_track(KEY6[6], show_strip = TRUE, show_x = TRUE) +
  plot_layout(heights = c(0.42, rep(1, 8))) +
  plot_annotation(
    title = "Pooled association mapping and cross-based linkage mapping over the same genome",
    subtitle = wrap_txt(paste0(
      "Top two tracks: pooled RNAi competition GWAS (n = 93 isotypes), dashed line ",
      "Bonferroni. Lower six: F2 bulk-segregant LOD scans in the two crosses built from ",
      "pooled outliers, dashed line genome-wide threshold, shaded band the support ",
      "interval. The dotted verticals, labelled in the key row at the top, mark the eight ",
      "loci this figure set follows; they are the rows of FIG10 and the columns of FIG15. ",
      "Note the y scales differ between tracks; LOD and -log10 p are not on a common ",
      "scale."), 130),
    caption = "Position (Mb)",
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8.5, colour = "grey25"),
                  plot.caption = element_text(hjust = 0.5, size = 9)))

save_fig(fig1, "FIG1_master_intersection", 11, 14)

## -- condensed variant: the three tracks that carry the argument ------------
fig1b <- loci_key() / gw_track("mig-6") /
  cross_track("N2xXZ1516 | pos1 vs mig6") /
  cross_track("JU1793xJU2466 | mig6 vs pos1", show_strip = TRUE, show_x = TRUE) +
  plot_layout(heights = c(0.42, 1, 1, 1)) +
  plot_annotation(
    title = "The *mig-6* vs *pos-1* contrast recovers the pooled *mig-6* QTL",
    subtitle = wrap_md(paste0(
      "The pooled *mig-6* GWAS peak, V:14,647,434 (−log<sub>10</sub>*p* = 7.71), sits on the ",
      "chromosome-V ridge that the N2 × XZ1516 *pos-1*-vs-*mig-6* contrast maps at LOD 730. ",
      "The same contrast in JU1793 × JU2466 is flat there and instead peaks on ",
      "chromosome X, so the two crosses expose different alleles. Dotted verticals are ",
      "labelled in the key row at the top."), 125),
    caption = "Position (Mb)",
    theme = theme(plot.title = element_markdown(face = "bold", size = 12),
                  plot.subtitle = element_markdown(size = 8.5, colour = "grey25"),
                  plot.caption = element_text(hjust = 0.5, size = 9)))

save_fig(fig1b, "FIG1b_condensed_intersection", 11, 6.4)

## ===========================================================================
## FIG 2  Pooled phenotypes: ranked strains, cross parents highlighted
## ===========================================================================
PARENTS <- c("XZ1516", "JU1793", "JU2466")
ph <- b$pheno %>% filter(sample_type == "pellet", Gene %in% c("mig-6", "pos-1"))

rank_panel <- function(g) {
  d <- ph %>% filter(Gene == g) %>%
    arrange(pheno) %>% mutate(i = row_number(), is_par = strain %in% PARENTS)
  ggplot(d, aes(i, pheno)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_segment(aes(xend = i, yend = 0), linewidth = 0.28, colour = "grey72") +
    geom_point(aes(colour = is_par, size = is_par)) +
    geom_text_repel(data = d %>% filter(is_par), aes(label = strain),
                    size = 3, fontface = "bold", seed = 42,
                    min.segment.length = 0, box.padding = 0.5,
                    segment.size = 0.3, colour = "grey15") +
    scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = PAL[["hi"]]),
                        guide = "none") +
    scale_size_manual(values = c(`FALSE` = 1.1, `TRUE` = 2.4), guide = "none") +
    labs(x = "wild isotype, ranked", y = "RNAi − control Δ frequency",
         title = paste0("*", g, "* RNAi")) +
    theme_pub(10) + theme(plot.title = element_markdown(face = "bold"))
}

fig2 <- (rank_panel("mig-6") | rank_panel("pos-1")) +
  plot_annotation(
    title = "Pooled competition identifies the resistance and sensitivity extremes",
    subtitle = wrap_txt(paste0(
      "Per-strain change in pool frequency under RNAi relative to the control baseline, ",
      "timepoint 2, pellet fraction (n = 93). Highlighted strains became cross parents. ",
      "XZ1516 is the most sensitive isotype for both targets; JU1793 the most resistant ",
      "for mig-6 and second for pos-1."), 118),
    tag_levels = "a",
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8.5, colour = "grey25"),
                  plot.tag = element_text(face = "bold", size = 12)))

save_fig(fig2, "FIG2_pooled_phenotype_ranks", 10, 4.8)

## ===========================================================================
## FIG 3  pos-1 vs mig-6 phenotype, parents named
## ===========================================================================
sc <- ph %>% select(strain, Gene, pheno) %>%
  pivot_wider(names_from = Gene, values_from = pheno) %>%
  rename(mig6 = `mig-6`, pos1 = `pos-1`) %>%
  filter(!is.na(mig6), !is.na(pos1)) %>%
  mutate(is_par = strain %in% PARENTS)

rho <- cor(sc$mig6, sc$pos1, method = "spearman", use = "complete.obs")
ct  <- suppressWarnings(cor.test(sc$mig6, sc$pos1, method = "spearman"))

fig3 <- ggplot(sc, aes(mig6, pos1)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey40",
              fill = "grey88", linewidth = 0.4) +
  geom_point(aes(colour = is_par, size = is_par)) +
  geom_text_repel(data = sc %>% filter(is_par), aes(label = strain),
                  size = 3.4, fontface = "bold", seed = 7,
                  min.segment.length = 0, box.padding = 0.7, segment.size = 0.3) +
  scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = PAL[["hi"]]), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.5, `TRUE` = 3), guide = "none") +
  annotate("richtext", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.25,
           label = sprintf("Spearman ρ = %.2f<br>*p* = %.1g", rho, ct$p.value),
           size = 3.2, fill = NA, label.color = NA) +
  labs(x = "*mig-6* RNAi − control Δ frequency",
       y = "*pos-1* RNAi − control Δ frequency",
       title = "Sensitivity is correlated across targets, and the extremes are shared",
       subtitle = wrap_md(paste0(
         "Each point one wild isotype. XZ1516 is the most sensitive for both targets and ",
         "JU1793 among the most resistant, which is why they were chosen as cross parents. ",
         "JU2466, near the middle, is the JU-cross counter-parent."), 78)) +
  theme_pub(10) +
  theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(),
        plot.subtitle = element_markdown(size = 8.5, colour = "grey25"))

save_fig(fig3, "FIG3_phenotype_scatter", 6.6, 6)

## ===========================================================================
## FIG 4  Phenotype heatmap: every strain x every RNAi target
## ===========================================================================
hm <- b$pheno %>% filter(sample_type == "pellet") %>% select(strain, Gene, pheno)
ord  <- hm %>% group_by(strain) %>% summarise(m = mean(pheno, na.rm = TRUE)) %>%
  arrange(m) %>% pull(strain)
gord <- hm %>% group_by(Gene) %>% summarise(s = sd(pheno, na.rm = TRUE)) %>%
  arrange(desc(s)) %>% pull(Gene)

lim <- max(abs(hm$pheno), na.rm = TRUE)
fig4 <- hm %>%
  mutate(strain = factor(strain, levels = ord), Gene = factor(Gene, levels = gord)) %>%
  ggplot(aes(strain, Gene, fill = pheno)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#8C2D04", mid = "grey96", high = "#08519C",
                       midpoint = 0, limits = c(-lim, lim),
                       name = "RNAi − ctrl\nΔ frequency") +
  labs(x = NULL, y = NULL,
       title = "Pooled RNAi sensitivity across 93 wild isotypes and 10 targets",
       subtitle = wrap_txt(paste0(
         "Pellet fraction, timepoint 2. Strains ordered by mean response (most sensitive at ",
         "left); targets by response variance. Cross parents in bold red. XZ1516 sits at the ",
         "sensitive extreme for nearly every target, which is what makes it a general ",
         "rather than target-specific outlier."), 135)) +
  theme_pub(9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.4,
                                   face = ifelse(ord %in% PARENTS, "bold", "plain"),
                                   colour = ifelse(ord %in% PARENTS, PAL[["hi"]], "grey35")),
        axis.text.y = element_text(size = 8, face = "italic"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "right", legend.title = element_text(size = 7),
        legend.text = element_text(size = 6.5),
        legend.key.height = grid::unit(22, "pt"))

save_fig(fig4, "FIG4_phenotype_heatmap", 11, 3.9)

## ===========================================================================
## FIG 5  Mirrored pooled Manhattan, with cross peaks ticked underneath
## ===========================================================================
mir <- gwas_long %>% mutate(y = ifelse(gene == "mig-6", neglog10p, -neglog10p))

cross_ticks <- b$ivs %>% filter(key %in% KEY6, significant) %>%
  transmute(chrom = fct_chrom(chrom), pos.mb = peak.position / 1e6,
            cross = factor(cross, levels = names(CROSS_LAB)))

ymax <- max(mir$neglog10p, na.rm = TRUE)
fig5 <- ggplot(mir, aes(pos.mb, y)) +
  guide_lines() +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey25") +
  geom_hline(yintercept = c(-1, 1) * b$GWAS_BF, linetype = "dashed",
             linewidth = 0.3, colour = "grey45") +
  geom_point(aes(colour = gene), size = 0.2, alpha = 0.55, stroke = 0) +
  geom_point(data = b$gwas_sig %>%
               transmute(chrom = fct_chrom(chr), pos.mb = ps / 1e6,
                         y = ifelse(gene == "mig-6", neglog10p, -neglog10p), gene),
             aes(y = y), colour = PAL[["hi"]], size = 1.7) +
  geom_point(data = cross_ticks, aes(y = -ymax * 1.16, shape = cross, fill = cross),
             size = 1.5, colour = "grey20", stroke = 0.2) +
  scale_colour_manual(values = c(`mig-6` = PAL[["pooled"]], `pos-1` = "#7C6A9C"),
                      name = NULL, labels = c("mig-6 (above axis)", "pos-1 (below axis)")) +
  scale_shape_manual(values = c(24, 21), name = "cross QTL peak",
                     labels = unname(CROSS_LAB)) +
  scale_fill_manual(values = unname(CROSS_COL), name = "cross QTL peak",
                    labels = unname(CROSS_LAB)) +
  genome_facets() +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.18, 0.06))) +
  labs(x = "Position (Mb)", y = "−log<sub>10</sub> *p*",
       title = "Pooled GWAS for both RNAi targets, with cross QTL peaks marked below",
       subtitle = wrap_txt(paste0(
         "mig-6 above the axis, pos-1 below. Red points are the genome-wide significant ",
         "peaks; dashed lines Bonferroni. Symbols along the bottom mark every significant ",
         "support-interval peak from the six focal cross contrasts, showing how little of ",
         "the crosses' signal the panel recovers."), 130)) +
  theme_pub(10) +
  theme(axis.title.y = element_markdown(), plot.subtitle = element_text(size = 8.5),
        legend.position = "top", legend.box = "horizontal",
        legend.margin = margin(b = -4))

save_fig(fig5, "FIG5_mirrored_manhattan", 11, 5)

## ===========================================================================
## FIG 6  Cross design: which cross can test which pooled QTL
## ===========================================================================
pg <- b$parent_gt_qtl %>%
  mutate(locus = sprintf("%s:%s\n(pooled %s QTL)", chrom,
                         format(ps, big.mark = ",", trim = TRUE), gwas_gene)) %>%
  select(locus, ref, alt, N2, XZ1516, JU1793, JU2466) %>%
  pivot_longer(c(N2, XZ1516, JU1793, JU2466), names_to = "strain", values_to = "gt") %>%
  mutate(allele = case_when(gt == "0/0" ~ ref, gt == "1/1" ~ alt, TRUE ~ "?"),
         state  = case_when(gt == "0/0" ~ "reference allele",
                            gt == "1/1" ~ "alternate allele (resistance-associated)",
                            TRUE        ~ "no call"),
         strain = factor(strain, levels = c("N2", "XZ1516", "JU1793", "JU2466")))

p6a <- ggplot(pg, aes(strain, locus, fill = state)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = allele), size = 3.6, fontface = "bold",
            colour = ifelse(pg$state == "reference allele", "grey20", "white")) +
  scale_fill_manual(values = c(`reference allele` = "grey88",
                               `alternate allele (resistance-associated)` = PAL[["hi"]],
                               `no call` = "grey60"), name = NULL,
                    breaks = c("alternate allele (resistance-associated)",
                               "reference allele")) +
  labs(x = NULL, y = NULL, title = "Parental alleles at the pooled QTL",
       subtitle = wrap_txt(paste0(
         "Only JU1793 carries the alternate, resistance-associated allele. N2 x XZ1516 is ",
         "uninformative for both pooled QTL, so any QTL it maps must be a different locus."),
         58)) +
  theme_pub(9) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = 7.5, lineheight = 0.95),
        legend.position = "bottom", legend.direction = "vertical",
        plot.subtitle = element_text(size = 7.6, colour = "grey25"))

gtp <- b$panel_gt %>%
  transmute(gwas_gene, chrom, ps, strain, genotype, allele1, allele0) %>%
  left_join(b$pheno %>% filter(sample_type == "pellet") %>%
              select(strain, Gene, pheno), by = "strain") %>%
  filter(Gene == gwas_gene, !is.na(genotype), genotype != "het") %>%
  mutate(grp = factor(ifelse(genotype == allele1, "alternate", "reference"),
                      levels = c("reference", "alternate")),
         panel = sprintf("%s / %s:%s", gwas_gene, chrom,
                         format(ps, big.mark = ",", trim = TRUE)))

p6b <- ggplot(gtp, aes(grp, pheno)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_boxplot(aes(fill = grp), outlier.shape = NA, width = 0.55,
               linewidth = 0.35, alpha = 0.6) +
  geom_jitter(width = 0.14, height = 0, size = 1.2, alpha = 0.6, colour = "grey30") +
  geom_point(data = gtp %>% filter(strain %in% PARENTS),
             colour = PAL[["hi"]], size = 2.6) +
  geom_text_repel(data = gtp %>% filter(strain %in% PARENTS), aes(label = strain),
                  size = 2.8, fontface = "bold", seed = 3,
                  min.segment.length = 0, box.padding = 0.6, segment.size = 0.3) +
  scale_fill_manual(values = c(reference = "grey80", alternate = PAL[["hi"]]),
                    guide = "none") +
  facet_wrap(~panel, scales = "free_y", nrow = 1) +
  labs(x = "allele at the pooled QTL peak", y = "RNAi − control Δ frequency",
       title = "Pooled QTL effect, and where the parents sit",
       subtitle = wrap_txt(paste0(
         "Each point one isotype, phenotype for the matching RNAi target. The alternate ",
         "allele raises the phenotype at both loci. Note the free y scales."), 62)) +
  theme_pub(9) +
  theme(plot.subtitle = element_text(size = 7.6, colour = "grey25"))

inf <- b$inform %>%
  mutate(chrom = fct_chrom(chrom), pos.mb = bin / 1e6,
         cross = factor(cross, levels = names(CROSS_LAB)))

p6c <- ggplot(inf, aes(pos.mb, discordance, colour = cross)) +
  guide_lines() +
  geom_line(linewidth = 0.3, alpha = 0.85) +
  scale_colour_manual(values = unname(CROSS_COL), labels = unname(CROSS_LAB), name = NULL) +
  genome_facets() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.06))) +
  labs(x = "Position (Mb)", y = "parental\ndiscordance",
       title = "What each cross can see",
       subtitle = wrap_txt(paste0(
         "Fraction of CeNDR isotype variants per 100 kb at which the two parents carry ",
         "different alleles. A cross is blind where this approaches zero: JU1793 x JU2466 ",
         "is uninformative over most of chromosomes I and IV and the right of V, which is ",
         "where N2 x XZ1516 carries its strongest peaks. Dotted verticals are the eight ",
         "loci of FIG10/FIG15, labelled in FIG1's key row."),
         135)) +
  theme_pub(9) +
  theme(legend.position = "top", legend.margin = margin(b = -6),
        plot.subtitle = element_text(size = 8, colour = "grey25"))

fig6 <- (p6a | p6b) / p6c +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(
    title = "Crosses were built from pooled outliers, and their parents decide which loci are testable",
    tag_levels = "a",
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.tag = element_text(face = "bold", size = 12)))

save_fig(fig6, "FIG6_cross_design", 12, 8.6)

## ===========================================================================
## FIG 7  All six focal cross scans, one facet grid
## ===========================================================================
trk_order <- b$scan_meta %>% filter(key %in% KEY6) %>%
  mutate(cross = factor(cross, levels = names(CROSS_LAB))) %>%
  arrange(cross, label) %>%
  mutate(track = paste0(CROSS_LAB[as.character(cross)], "\n", label)) %>%
  pull(track)

sl6f <- sl6 %>%
  mutate(track = factor(paste0(CROSS_LAB[as.character(cross)], "\n", label),
                        levels = trk_order))
bd6f <- bd6 %>% left_join(sl6f %>% distinct(key, track), by = "key")

fig7 <- ggplot(sl6f, aes(pos.mb, LOD)) +
  guide_lines() +
  geom_ribbon(data = bd6f, aes(ymin = 0, ymax = LOD), fill = PAL[["hi"]], alpha = 0.5) +
  geom_line(linewidth = 0.28, colour = "grey15") +
  geom_hline(yintercept = b$CROSS_THR, linetype = "dashed",
             linewidth = 0.28, colour = "grey45") +
  facet_grid(track ~ chrom, scales = "free", space = "free_x", switch = "x") +
  scale_x_continuous(breaks = seq(0, 25, 5), expand = expansion(mult = 0.02),
                     guide = guide_axis(check.overlap = TRUE)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Position (Mb)", y = "LOD",
       title = "F2 bulk-segregant LOD scans, both crosses, harmonised",
       subtitle = wrap_txt(paste0(
         "Support intervals (shaded) use a drop of 10% of each chromosome's own peak. ",
         "Dashed line: genome-wide threshold (alpha = 0.05 over 2000 effective tests), LOD ",
         sprintf("%.2f", b$CROSS_THR),
         ". Every row has its own y scale. Dotted verticals are the eight loci of ",
         "FIG10/FIG15, labelled in FIG1's key row."), 130)) +
  theme_pub(9) +
  theme(strip.text.y = element_text(size = 7.5, angle = 0, hjust = 0, lineheight = 1),
        panel.spacing.y = grid::unit(5, "pt"),
        plot.subtitle = element_text(size = 8.5, colour = "grey25"))

save_fig(fig7, "FIG7_cross_scans_all", 11.5, 8.4)

## ===========================================================================
## FIG 8  Locus zooms: pooled GWAS + both crosses + gene track
## ===========================================================================
## colour = which contrast, linetype = which cross. Two small stacked legends,
## so the panel stays narrow enough to tile four across a page.
zoom_locus <- function(lc, half = 1.2e6, legend = TRUE) {
  lo <- lc$pos - half; hi <- lc$pos + half; ch <- lc$chrom

  gwz <- b$gwas %>% filter(chr == ch, ps >= lo, ps <= hi) %>% mutate(pos.mb = ps / 1e6)
  pg <- ggplot(gwz, aes(pos.mb, neglog10p, colour = gene)) +
    geom_hline(yintercept = b$GWAS_BF, linetype = "dashed",
               linewidth = 0.28, colour = "grey45") +
    geom_point(size = 0.35, alpha = 0.6, stroke = 0) +
    geom_vline(xintercept = lc$pos / 1e6, linetype = "dotted", linewidth = 0.35) +
    scale_colour_manual(values = c(`mig-6` = PAL[["pooled"]], `pos-1` = "#7C6A9C"),
                        name = "pooled GWAS") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    coord_cartesian(xlim = c(lo, hi) / 1e6) +
    labs(x = NULL, y = "pooled\n−log<sub>10</sub>*p*") +
    theme_pub(8.5) +
    theme(axis.title.y = element_markdown(size = 7.5, lineheight = 1),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = if (legend) "right" else "none",
          legend.title = element_text(size = 7), legend.text = element_text(size = 7))

  crz <- sl6 %>% filter(as.character(chrom) == ch,
                        pos.mb >= lo / 1e6, pos.mb <= hi / 1e6)
  ivz <- b$ivs %>% filter(key %in% KEY6, significant, chrom == ch, rcon >= lo, lcon <= hi)

  pc <- ggplot(crz, aes(pos.mb, LOD, colour = contrast, linetype = cross, group = key)) +
    { if (nrow(ivz)) geom_rect(data = ivz, inherit.aes = FALSE,
                               aes(xmin = lcon / 1e6, xmax = rcon / 1e6,
                                   ymin = -Inf, ymax = Inf),
                               fill = PAL[["hi"]], alpha = 0.06) } +
    geom_hline(yintercept = b$CROSS_THR, linetype = "dashed",
               linewidth = 0.28, colour = "grey45") +
    geom_line(linewidth = 0.42) +
    geom_vline(xintercept = lc$pos / 1e6, linetype = "dotted", linewidth = 0.35) +
    scale_colour_manual(values = COND_COL, labels = COND_LAB, name = "contrast",
                        breaks = names(COND_LAB)) +
    scale_linetype_manual(values = c(`N2xXZ1516` = "solid",
                                     `JU1793xJU2466` = "21"),
                          labels = unname(CROSS_LAB), name = "cross") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
    coord_cartesian(xlim = c(lo, hi) / 1e6) +
    labs(x = NULL, y = "cross LOD") +
    theme_pub(8.5) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 7.5),
          legend.position = if (legend) "right" else "none",
          legend.box = "vertical", legend.spacing.y = grid::unit(1, "pt"),
          legend.title = element_text(size = 7),
          legend.text = element_markdown(size = 7)) +
    guides(colour = guide_legend(order = 1, override.aes = list(linetype = "solid")),
           linetype = guide_legend(order = 2))

  gz <- b$rnai_genes %>% filter(chrom == ch, end >= lo, start <= hi)
  pgene <- ggplot(gz) +
    geom_segment(aes(x = start / 1e6, xend = end / 1e6, y = 1, yend = 1),
                 linewidth = 3.2, lineend = "butt",
                 colour = ifelse(gz$class == "dsRNA uptake", PAL[["uptake"]], PAL[["gene"]])) +
    { if (nrow(gz)) geom_text_repel(
        aes(x = (start + end) / 2e6, y = 1, label = gene),
        colour = ifelse(gz$class == "dsRNA uptake", PAL[["uptake"]], PAL[["gene"]]),
        size = 2.7, fontface = "italic", seed = 11, direction = "y",
        ylim = c(1.05, 1.6), min.segment.length = 0, segment.size = 0.25,
        box.padding = 0.25) } +
    geom_vline(xintercept = lc$pos / 1e6, linetype = "dotted", linewidth = 0.35) +
    coord_cartesian(xlim = c(lo, hi) / 1e6, ylim = c(0.85, 1.7)) +
    labs(x = sprintf("chromosome %s (Mb)", ch), y = NULL,
         caption = if (nrow(gz))
           paste0("RNAi-pathway gene positions (uptake genes purple). Positional only \u2014 ",
                  "see FIG12 for which carry parental coding variation.") else
           "no RNAi-pathway gene in window") +
    theme_pub(8.5) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.caption = element_text(hjust = 0, size = 6.5, colour = "grey40"))

  (pg / pc / pgene) +
    plot_layout(heights = c(1, 1.6, 0.7), guides = "collect") +
    plot_annotation(title = lc$label, subtitle = lc$role,
                    theme = theme(plot.title = element_text(face = "bold", size = 10),
                                  plot.subtitle = element_text(size = 8, colour = "grey35")))
}

for (i in seq_len(nrow(LOCI))) {
  lc <- LOCI[i, ]
  save_fig(wrap_elements(zoom_locus(lc)),
           sprintf("FIG8_zoom_%s", gsub("[^A-Za-z0-9]+", "_", lc$locus)), 7, 5.4)
}

## the four loci the argument leans on, in one sheet: legend once, on the last
fig8 <- wrap_plots(lapply(seq_along(c("III:13.7", "V:14.6", "I:1.9", "X:6.9")), function(j) {
  l <- c("III:13.7", "V:14.6", "I:1.9", "X:6.9")[j]
  wrap_elements(zoom_locus(LOCI %>% filter(locus == l), legend = (j == 4)))
}), nrow = 1, widths = c(1, 1, 1, 1.55)) +
  plot_annotation(
    title = "Four loci, three lines of evidence each",
    subtitle = wrap_txt(paste0(
      "Top: pooled GWAS in the same window. Middle: cross LOD, colour by contrast and ",
      "line style by cross, shaded where a support interval covers the window. Bottom: ",
      "RNAi-pathway genes. Dotted vertical is the named position. y scales are per panel."),
      150),
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8.5, colour = "grey25")))

save_fig(fig8, "FIG8_locus_zooms_combined", 20, 5.8)

## ===========================================================================
## FIG 9  Contrast decomposition: what cancels and what does not
## ===========================================================================
## A locus acting on dsRNA uptake shifts the pos-1 and mig-6 pools the same way,
## so it largely cancels in the RNAi-vs-RNAi contrast; a target-specific locus
## survives it.
dec_region <- tribble(
  ~chrom, ~lo,    ~hi,     ~keys,  ~title,                                 ~note,
  "III",  12.0e6, 14.2e6,  "nxz",  "chromosome III · *sid-2* region",
  "shared: LOD 815 and 589 collapse to 48 in the contrast",
  "V",    13.2e6, 14.9e6,  "nxz",  "chromosome V · pooled *mig-6* QTL",
  "*mig-6*-specific: *pos-1* is flat, the contrast is maximal",
  "I",    1.2e6,  2.6e6,   "nxz",  "chromosome I · 1.9 Mb",
  "*mig-6*-specific and cross-only: pooled −log<sub>10</sub>*p* = 1.9",
  "X",    5.4e6,  8.2e6,   "juju", "chromosome X · 6.9-7.2 Mb",
  "*mig-6*-specific, JU cross only: pooled −log<sub>10</sub>*p* = 2.7")

nxz3 <- KEY6[1:3]; juju3 <- KEY6[4:6]

dec_panel <- function(r) {
  keys <- if (r$keys == "nxz") nxz3 else juju3
  d <- scan_long(keys) %>%
    filter(as.character(chrom) == r$chrom, pos.mb >= r$lo / 1e6, pos.mb <= r$hi / 1e6)
  gz <- b$rnai_genes %>% filter(chrom == r$chrom, end >= r$lo, start <= r$hi)

  ggplot(d, aes(pos.mb, LOD, colour = contrast)) +
    geom_hline(yintercept = b$CROSS_THR, linetype = "dashed",
               linewidth = 0.28, colour = "grey45") +
    { if (nrow(gz)) geom_vline(data = gz, aes(xintercept = start / 1e6),
                               linetype = "dotted", linewidth = 0.32,
                               colour = PAL[["gene"]]) } +
    { if (nrow(gz)) geom_text(data = gz, aes(x = start / 1e6, y = Inf, label = gene),
                              inherit.aes = FALSE, angle = 90, hjust = 1.05,
                              vjust = -0.4, size = 2.5, fontface = "italic",
                              colour = PAL[["gene"]]) } +
    geom_line(linewidth = 0.5) +
    scale_colour_manual(values = COND_COL, labels = COND_LAB, name = NULL,
                        breaks = names(COND_LAB)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
    labs(x = sprintf("chromosome %s (Mb)", r$chrom), y = "LOD",
         title = r$title, subtitle = r$note) +
    theme_pub(9) +
    theme(plot.title = element_markdown(face = "bold", size = 9),
          plot.subtitle = element_markdown(size = 7.8, colour = "grey30"),
          legend.position = "top", legend.margin = margin(b = -7),
          legend.text = element_markdown(size = 7))
}

fig9 <- (dec_panel(dec_region[1, ]) | dec_panel(dec_region[2, ])) /
  (dec_panel(dec_region[3, ]) | dec_panel(dec_region[4, ])) +
  plot_annotation(
    title = "Contrasting two RNAi conditions separates shared dsRNA uptake from target-specific modifiers",
    subtitle = wrap_md(paste0(
      "Panels a-c: N2 × XZ1516. Panel d: JU1793 × JU2466. A locus that changes dsRNA uptake ",
      "moves the *pos-1* and *mig-6* pools together, so it largely cancels in the ",
      "RNAi-vs-RNAi contrast (green); a target-specific locus does not. Grey dotted lines ",
      "mark RNAi-pathway gene positions, which is annotation and not a candidate call \u2014 ",
      "of these only *sid-2* carries coding variation that distinguishes it within its ",
      "interval (FIG12). Every panel has its own y scale."), 140),
    tag_levels = "a",
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_markdown(size = 8.5, colour = "grey25"),
                  plot.tag = element_text(face = "bold", size = 12)))

save_fig(fig9, "FIG9_contrast_decomposition", 11.5, 8.2)

## ===========================================================================
## FIG 10  Evidence matrix: locus x experiment
## ===========================================================================
lod_at <- function(key, ch, ps) {
  d <- b$scans[[key]]; d <- d[d$chrom == ch & is.finite(d$LOD), ]
  if (!nrow(d)) return(NA_real_)
  d$LOD[which.min(abs(d$physical.position - ps))]
}
gw_at <- function(g, ch, ps, half = 5e4) {
  v <- b$gwas$neglog10p[b$gwas$gene == g & b$gwas$chr == ch & abs(b$gwas$ps - ps) <= half]
  if (!length(v)) NA_real_ else max(v)
}

## pretty row label, markdown-free (axis text is plain)
pretty_cond <- function(x) x %>%
  gsub("ht115", "HT115", .) %>% gsub("mig6", "mig-6", .) %>% gsub("pos1", "pos-1", .)

ev <- pmap_dfr(LOCI %>% select(locus, label, chrom, pos, role),
               function(locus, label, chrom, pos, role) {
  bind_rows(
    tibble(track = "pooled GWAS · mig-6", value = gw_at("mig-6", chrom, pos), scale = "gwas"),
    tibble(track = "pooled GWAS · pos-1", value = gw_at("pos-1", chrom, pos), scale = "gwas"),
    map_dfr(KEY6, function(k) {
      m <- b$scan_meta %>% filter(key == k)
      tibble(track = paste0(CROSS_LAB[[m$cross]], " · ", pretty_cond(m$label)),
             value = lod_at(k, chrom, pos), scale = "lod")
    })) %>%
    mutate(locus = label, role = role, .before = 1)
})

trk_levels <- c("pooled GWAS · mig-6", "pooled GWAS · pos-1",
                paste0(CROSS_LAB[["N2xXZ1516"]], " · ",
                       pretty_cond(c("ht115 vs mig6", "ht115 vs pos1", "pos1 vs mig6"))),
                paste0(CROSS_LAB[["JU1793xJU2466"]], " · ",
                       pretty_cond(c("ht115 vs mig6", "ht115 vs pos1", "mig6 vs pos1"))))

ev <- ev %>%
  mutate(track = factor(track, levels = rev(trk_levels)),
         locus = factor(locus, levels = LOCI$label),
         rel   = ifelse(scale == "gwas", value / b$GWAS_BF, value / b$CROSS_THR),
         shown = ifelse(scale == "gwas", sprintf("%.1f", value),
                        ifelse(value >= 100, sprintf("%.0f", value), sprintf("%.1f", value))))

fig10 <- ggplot(ev, aes(locus, track, fill = pmin(rel, 40))) +
  geom_tile(colour = "white", linewidth = 0.9) +
  geom_text(aes(label = shown, colour = pmin(rel, 40) > 12),
            size = 2.9, fontface = "bold") +
  annotate("segment", x = 0.5, xend = nrow(LOCI) + 0.5, y = 6.5, yend = 6.5,
           linewidth = 0.6, colour = "grey30") +
  scale_fill_gradientn(colours = c("grey96", "#D9E7F1", "#6BAED6", "#2171B5", "#08306B"),
                       trans = "sqrt", name = "signal /\nthreshold",
                       breaks = c(1, 5, 15, 40), labels = c("1×", "5×", "15×", "≥40×")) +
  scale_colour_manual(values = c(`FALSE` = "grey25", `TRUE` = "white"), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "Which experiment sees which locus",
       subtitle = wrap_txt(paste0(
         "Cells give the pooled GWAS -log10 p (top two rows, above the rule) or the cross ",
         "LOD at that position. Shading is the value divided by that assay's own ",
         "genome-wide threshold, so 1x is bare significance and the crosses run to >100x. ",
         "The strip below colours how we read each locus."), 120)) +
  theme_pub(9) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(size = 8, angle = 30, hjust = 0),
        axis.text.y = element_text(size = 8),
        plot.subtitle = element_text(size = 8, colour = "grey25"),
        legend.key.height = grid::unit(16, "pt"))

role_strip <- ggplot(LOCI %>% mutate(locus = factor(label, levels = LOCI$label),
                                     role = factor(role, levels = ROLE_LEV)),
                     aes(locus, 1, fill = role)) +
  geom_tile(colour = "white", linewidth = 0.9) +
  scale_fill_manual(values = ROLE_COL, name = NULL, breaks = ROLE_LEV) +
  labs(x = NULL, y = NULL) +
  theme_pub(9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), legend.position = "bottom",
        legend.margin = margin(t = -4))

save_fig(fig10 / role_strip + plot_layout(heights = c(1, 0.075)),
         "FIG10_evidence_matrix", 10, 5.6)

## ===========================================================================
## FIG 11  The workflow, drawn
## ===========================================================================
box <- function(x, y, w, h, fill, label) tibble(x, y, w, h, fill, label)
bx <- bind_rows(
  box(1.05, 3.0, 1.8, 1.2, "#E8EEF4",
      "**Pooled competition**<br>93 wild isotypes<br>10 RNAi targets"),
  box(3.25, 3.0, 1.8, 1.2, "#E8EEF4",
      "**Per-strain phenotype**<br>RNAi − control<br>Δ pool frequency"),
  box(5.45, 3.72, 1.8, 1.0, "#DCE6F1",
      "**GWAS**<br>2 genome-wide QTL"),
  box(5.45, 2.32, 1.8, 1.0, "#F6E3CE",
      "**Outlier strains**<br>XZ1516, JU1793"),
  box(7.75, 3.72, 1.9, 1.2, "#DBEEDD",
      "**Verified**<br>chrV *mig-6* QTL<br>chrIII uptake locus"),
  box(7.75, 2.18, 1.9, 1.45, "#F6E3CE",
      "**F2 crosses**<br>N2 × XZ1516<br>JU1793 × JU2466<br>bulk-segregant xQTL"),
  box(10.1, 2.18, 1.9, 1.45, "#FBE3E1",
      "**New QTL**<br>chrI 1.9 Mb<br>chrV 13.8 Mb<br>chrX 6.9 Mb"))

ar <- tribble(
  ~x,    ~y,   ~xend, ~yend,
  1.95,  3.00, 2.35,  3.00,
  4.15,  3.15, 4.55,  3.60,
  4.15,  2.85, 4.55,  2.45,
  6.35,  2.32, 6.80,  2.28,
  8.70,  2.91, 8.70,  3.12,
  8.70,  2.18, 9.15,  2.18)

fig11 <- ggplot() +
  geom_tile(data = bx, aes(x, y, width = w, height = h, fill = fill),
            colour = "grey55", linewidth = 0.4) +
  geom_richtext(data = bx, aes(x, y, label = label), size = 3.1,
                fill = NA, label.color = NA, lineheight = 1.35) +
  geom_segment(data = ar, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = grid::unit(6, "pt"), type = "closed"),
               linewidth = 0.45, colour = "grey30") +
  annotate("richtext", x = 8.70, y = 4.55, size = 3, hjust = 0.5,
           label = "*crosses verify what 93 isotypes could only hint at*",
           fill = NA, label.color = NA, colour = "grey25") +
  annotate("richtext", x = 10.1, y = 1.30, size = 3, hjust = 0.5,
           label = "*and reach loci the panel cannot power*",
           fill = NA, label.color = NA, colour = "grey25") +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0.1, 11.15), ylim = c(1.05, 4.8)) +
  labs(title = "Pooled mapping nominates; crosses resolve",
       subtitle = wrap_txt(paste0(
         "The pooled assay yields both a genome-wide scan and a ranked list of extreme ",
         "strains. The extremes become cross parents, and the crosses both confirm pooled ",
         "QTL and expose loci that 93 isotypes cannot power."), 120)) +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 9, colour = "grey30", hjust = 0,
                                     margin = margin(b = 6)),
        plot.margin = margin(10, 10, 10, 10))

save_fig(fig11, "FIG11_workflow_schematic", 11, 4.4)

## ===========================================================================
## Numbers behind the figures
## ===========================================================================
write_tsv(ev %>% select(locus, role, track, value, scale, rel),
          file.path(OUT, "TABLE_evidence_matrix.tsv"))
write_tsv(b$isect_B %>% arrange(desc(cross_LOD)),
          file.path(OUT, "TABLE_gwas_in_cross_intervals.tsv"))
write_tsv(b$isect_A, file.path(OUT, "TABLE_cross_lod_at_gwas_peaks.tsv"))
write_tsv(b$isect_C %>% arrange(desc(cross_LOD)),
          file.path(OUT, "TABLE_rnai_genes_in_intervals.tsv"))
write_tsv(b$parent_gt_qtl, file.path(OUT, "TABLE_parental_genotypes_at_qtl.tsv"))
write_tsv(b$ivs, file.path(OUT, "TABLE_cross_support_intervals.tsv"))

cat("\ndone. ", length(list.files(OUT, pattern = "\\.pdf$")), " figures, ",
    length(list.files(OUT, pattern = "\\.tsv$")), " tables in ", OUT, "\n", sep = "")
