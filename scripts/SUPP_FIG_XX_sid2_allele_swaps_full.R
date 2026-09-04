## Supplement -- the full sid-2 editing series -------------------------------
##
##   Rscript scripts/SUPP_FIG_XX_sid2_allele_swaps_full.R
##     -> plots/SUPP_FIG_XX_sid2_allele_swaps_full.{pdf,png}
##
##   A  every construct, both food conditions
##   B  why the N94A edits exist: T96 is the threonine of an N-x-T sequon
##   C  AlphaFold model confidence, which limits what panel B of Figure 4 can
##      claim
##
## Figure 4 shows the residue-96 swaps on pos-1 only. This carries everything
## held back from it: the HT115 controls, the N94A constructs, the sequon that
## motivated them, and the model confidence.
##
## THE GLYCOSYLATION HYPOTHESIS, AND WHY IT IS NOT IN FIGURE 4
## T96 is the threonine of an N94-C95-T96 sequon, so T96K removes a predicted
## N-glycosylation site. N94A was built to test that: it removes the same
## sequon while leaving residue 96 alone, and should therefore phenocopy T96K
## if loss of the glycan is what matters. It does not.
##
##   background   construct   motif   pos-1 hatching
##   JU1793       96T (wt)    NxT     0.948
##   JU1793       94A         AxT     0.993      <- no glycan, fully resistant
##   JU1793       96K         NxK     0.531
##   JU2466       96K (wt)    NxK     0.045
##   JU2466       96T         NxT     0.184
##   JU2466       94A         AxK     0.425      <- helps more than 96T does
##
## Two things follow. The AxT construct cannot be glycosylated at N94 and is
## fully resistant, so the glycan is not required for resistance. And removing
## Asn94 gains more in JU2466 (+0.38) than restoring Thr96 does (+0.14), while
## in JU1793 it gains almost nothing (+0.045) -- so the two positions interact
## rather than acting through one shared modification.
##
## A description that fits all six numbers is that Lys96 is deleterious and an
## unglycosylated Asn94 makes it worse, either of which can be relieved
## independently. That is epistasis between two neighbouring residues, not a
## mechanism, and each genotype is a single plate, so it is reported here and
## not built into a main figure.
##
## RNAi DOSE: 50% pos-1 RNAi bacteria, diluted with HT115. The dose is NOT
## recorded in the source data file -- the condition column carries only
## "pos"/"ht115" -- so it comes from the lab record, not from the data. Figure
## 4B uses 25% and its hatching percentages are therefore not comparable with
## these; see FIGURE_CAPTIONS.txt.
##
## CAVEAT: one plate per strain per condition throughout. Intervals are Wilson
## binomial intervals on that plate's embryo count and describe counting
## uncertainty, not between-plate variability. No construct is replicated,
## except that JU2466_A and JU2466_B are independent lines of the same edit.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggtext)
  library(png)
})

OUT    <- "plots"
PERRES <- "data/structure/sid2_per_residue.tsv"
SWAP   <- "data/plate_rnai_phenotyping/20260409_ju2466swap_plus_N2A_swap.csv"
RIBBON <- "plots/assets/sid2_ecd_ribbon_plddt.png"
AFDIR  <- "data/structure/AF_sid2_wt_dimer"

COL_PT   <- "#2E4057"
COL_MARK <- "#C4302B"
PLDDT <- c(`Very high (>90)` = "#0053D6", `Confident (70-90)` = "#65CBF3",
           `Low (50-70)` = "#FFDB13", `Very low (<50)` = "#FF7D45")
TOPO_COL <- c(`Signal peptide` = "grey72", `Extracellular` = "#9EC5DE",
              `TM helix` = "#37474F", `Cytoplasmic` = "grey88")

panel_title <- function(letter, txt = NULL) {
  lt <- paste0("<span style='font-size:14pt;color:#111111'>**", letter,
               "**</span>")
  if (is.null(txt)) lt else paste0(lt, " ", txt)
}
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_markdown(size = base_size + 0.5),
          plot.subtitle = element_markdown(size = base_size - 2.5,
                                           colour = "grey30"),
          plot.title.position = "plot",
          legend.key.size = grid::unit(9, "pt"))
}
wrap_md <- function(txt, width = 96)
  paste(strwrap(txt, width = width), collapse = "<br>")
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## ===========================================================================
## A -- every construct, both conditions
## ===========================================================================
msg("panel A: full editing series")
z <- qnorm(0.975)
sw <- read_csv(SWAP, show_col_types = FALSE) %>%
  transmute(strain, genotype, motif = `glycosylation motif`,
            cond = ifelse(condition == "ht115", "HT115", "*pos-1*"),
            n = n_plated, hatched = n_plated - n_unhatched) %>%
  mutate(p = hatched / n,
         lo = pmax(0, (p + z^2 / (2 * n) -
                         z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
                     (1 + z^2 / n)),
         hi = pmin(1, (p + z^2 / (2 * n) +
                         z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
                     (1 + z^2 / n)),
         cond = factor(cond, levels = c("HT115", "*pos-1*")),
         background = ifelse(grepl("^JU1793", genotype),
                             "JU1793 background", "JU2466 background"))

## every construct should have both conditions; say so if one does not
gaps <- sw %>% count(genotype) %>% filter(n < 2)
msg("  ", n_distinct(sw$genotype), " constructs | missing a condition: ",
    if (nrow(gaps)) paste(gaps$genotype, collapse = ", ") else "none")

ORD <- c("JU1793[96T]", "JU1793[94A]", "JU1793[96K]",
         "JU2466[96T]", "JU2466[94A]", "JU2466_A[96K]", "JU2466_B[96K]")
stopifnot(setequal(unique(sw$genotype), ORD))
sw <- sw %>% mutate(genotype = factor(genotype, levels = rev(ORD)))

lab <- sw %>% filter(cond == "*pos-1*") %>%
  transmute(genotype, background,
            txt = sprintf("%s, n = %d", sub(".*\\[", "", sub("\\]", "", motif)),
                          n))

cat("\n== all constructs, embryos hatched (Wilson 95% CI) ==\n")
print(as.data.frame(sw %>%
  transmute(strain, genotype, motif, condition = gsub("\\*", "", cond),
            embryos = n, hatched = round(p, 3),
            CI = sprintf("%.3f-%.3f", lo, hi)) %>%
  arrange(condition, genotype)), row.names = FALSE)

pd <- position_dodge(width = 0.72)
pA <- ggplot(sw, aes(p, genotype, fill = cond)) +
  geom_col(position = pd, width = 0.68, colour = "grey25", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), position = pd, orientation = "y",
                width = 0.26, linewidth = 0.4, colour = "grey15") +
  geom_text(data = lab, aes(x = 1.03, y = genotype, label = txt),
            inherit.aes = FALSE, hjust = 0, size = 2.7, colour = "grey30") +
  facet_grid(background ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(HT115 = "grey80", `*pos-1*` = COL_PT),
                    name = NULL) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.30), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = 0)) +
  labs(x = "Embryos hatched", y = NULL,
       title = panel_title("A", "**Every construct, both food conditions**"),
       subtitle = wrap_md(paste0(
         "Control hatching is 97.7&ndash;100% for all seven constructs, so ",
         "the *pos-1* differences are not a property of the edits. N94A ",
         "removes the same sequon as T96K but does not phenocopy it: it is ",
         "neutral in JU1793 and gains more in JU2466 than restoring 96T does. ",
         "One plate per strain per condition."))) +
  theme_pub() +
  theme(strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold", size = 8.5,
                                    hjust = 0),
        axis.text.y = element_text(size = 8.8),
        legend.position = "top", legend.text = element_markdown(),
        legend.margin = margin(b = -6))

## ===========================================================================
## B -- the sequon that motivated the N94A edits
## ===========================================================================
msg("panel B: sequon context")
pr <- read_tsv(PERRES, show_col_types = FALSE) %>%
  mutate(topology = factor(topology, levels = names(TOPO_COL)))
L <- max(pr$resid)

AA3 <- c(ALA="A",ARG="R",ASN="N",ASP="D",CYS="C",GLN="Q",GLU="E",GLY="G",
         HIS="H",ILE="I",LEU="L",LYS="K",MET="M",PHE="F",PRO="P",SER="S",
         THR="T",TRP="W",TYR="Y",VAL="V")
aa <- AA3[pr$resname]

## N-x-S/T with x != P, derived from the same table the figure uses
seqons <- tibble(pos = pr$resid[which(aa == "N")]) %>%
  filter(pos + 2 <= L) %>%
  mutate(x = aa[match(pos + 1, pr$resid)], t = aa[match(pos + 2, pr$resid)],
         motif = paste0("N", x, t),
         topo = pr$topology[match(pos, pr$resid)]) %>%
  filter(x != "P", t %in% c("S", "T")) %>%
  mutate(side = ifelse(topo == "Cytoplasmic",
                       "Cytoplasmic (inaccessible)", "Accessible"),
         focal = pos == 94)
msg("  sequons: ", nrow(seqons), " (", sum(seqons$side == "Accessible"),
    " accessible) | ", paste(sprintf("%s%d", seqons$motif, seqons$pos),
                             collapse = " "))

dom <- (function(v) { r <- rle(as.character(v))
  tibble(value = factor(r$values, levels = names(TOPO_COL)),
         end = cumsum(r$lengths),
         start = cumsum(r$lengths) - r$lengths + 1) })(pr$topology)

pB <- ggplot() +
  geom_rect(data = dom, aes(xmin = start - 0.5, xmax = end + 0.5,
                            ymin = 0, ymax = 1, fill = value),
            colour = "grey30", linewidth = 0.25) +
  geom_point(data = seqons, aes(pos, -0.45, shape = side), size = 2,
             colour = "grey25", fill = "white", stroke = 0.5) +
  geom_point(data = seqons %>% filter(focal), aes(pos, -0.45), size = 2.8,
             shape = 21, colour = COL_MARK, fill = COL_MARK) +
  geom_segment(data = tibble(pos = 96), aes(x = pos, xend = pos, y = 1,
                                            yend = 1.45),
               linewidth = 0.4, colour = COL_MARK) +
  geom_richtext(data = tibble(pos = 96), aes(pos, 1.5, label = "T96K"),
                colour = COL_MARK, size = 3.1, vjust = 0, fill = NA,
                label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  geom_richtext(data = tibble(x = 94),
                aes(x, -0.95, label = "N94&ndash;C95&ndash;**T96**"),
                colour = COL_MARK, size = 3, vjust = 1, fill = NA,
                label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  scale_fill_manual(values = TOPO_COL, name = NULL) +
  scale_shape_manual(values = c(Accessible = 21,
                                `Cytoplasmic (inaccessible)` = 4),
                     name = "N-x-S/T sequon:") +
  scale_x_continuous(breaks = c(1, seq(50, 300, 50), L),
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(-1.3, 1.95), expand = expansion(0)) +
  guides(fill = guide_legend(order = 1, nrow = 1),
         shape = guide_legend(order = 2, nrow = 1)) +
  labs(x = "SID-2 residue", y = NULL,
       title = panel_title("B", "**T96 is the threonine of a predicted N-glycosylation sequon**"),
       subtitle = wrap_md(paste0(
         "Nine N-x-S/T motifs; the four on the cytoplasmic side cannot be ",
         "glycosylated. T96K removes the N94 sequon, which is why N94A was ",
         "built as its test. Panel A shows that test failing."))) +
  theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom", legend.box = "vertical",
        legend.box.just = "left", legend.margin = margin(t = 2),
        legend.spacing.y = grid::unit(6, "pt"),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8.2))

## ===========================================================================
## C -- model confidence
## ===========================================================================
msg("panel C: model confidence")
conf <- list.files(AFDIR, pattern = "summary_confidences_.*json$",
                   full.names = TRUE) %>%
  map_dfr(function(f) {
    j <- jsonlite::fromJSON(f)
    tibble(model = sub(".*_(\\d+)\\.json$", "\\1", f),
           iptm = j$iptm, ptm = j$ptm,
           disordered = j$fraction_disordered)
  }) %>% arrange(model)
msg("  ipTM ", sprintf("%.2f-%.2f", min(conf$iptm), max(conf$iptm)),
    " | pTM ", sprintf("%.2f-%.2f", min(conf$ptm), max(conf$ptm)),
    " | disordered ", sprintf("%.0f-%.0f%%", 100 * min(conf$disordered),
                              100 * max(conf$disordered)))

dom_pl <- pr %>% group_by(topology) %>%
  summarise(mean_plddt = mean(plddt), n = n(), .groups = "drop") %>%
  mutate(band = cut(mean_plddt, c(-Inf, 50, 70, 90, Inf),
                    labels = rev(names(PLDDT))))
cat("\n== per-domain mean pLDDT ==\n")
print(as.data.frame(dom_pl %>% transmute(topology, residues = n,
        mean_plddt = round(mean_plddt, 1), band)), row.names = FALSE)

img <- readPNG(RIBBON)
pC <- ggplot() +
  annotation_raster(img, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                    interpolate = TRUE) +
  coord_fixed(ratio = dim(img)[1] / dim(img)[2], xlim = c(0, 1),
              ylim = c(0, 1), expand = FALSE) +
  labs(title = panel_title("C", "**Model confidence limits what the structure can be used for**"),
       subtitle = wrap_md(sprintf(paste0(
         "The same extracellular domain as Figure 4B, coloured by pLDDT ",
         "instead. Mean pLDDT by region: extracellular %.1f, signal peptide ",
         "%.1f, TM helix %.1f, cytoplasmic %.1f &mdash; only the ",
         "extracellular domain reaches the confident band, and the variant ",
         "neighbourhood is confident within it (N94 70.5, C95 78.0, T96 ",
         "79.3). The prediction is a dimer, but across all five models ipTM ",
         "is %.2f&ndash;%.2f and pTM %.2f&ndash;%.2f with %.0f&ndash;%.0f%% ",
         "of the model called disordered, so the interface is not evidence of ",
         "dimerisation and is not shown anywhere in this manuscript."),
         dom_pl$mean_plddt[dom_pl$topology == "Extracellular"],
         dom_pl$mean_plddt[dom_pl$topology == "Signal peptide"],
         dom_pl$mean_plddt[dom_pl$topology == "TM helix"],
         dom_pl$mean_plddt[dom_pl$topology == "Cytoplasmic"],
         min(conf$iptm), max(conf$iptm), min(conf$ptm), max(conf$ptm),
         100 * min(conf$disordered), 100 * max(conf$disordered)))) +
  theme_void(base_size = 11) +
  theme(plot.title = element_markdown(size = 11.5),
        plot.subtitle = element_markdown(size = 8.5, colour = "grey30"),
        plot.title.position = "plot")

key <- tibble(lab = factor(names(PLDDT), levels = names(PLDDT)))
p_key <- ggplot(key, aes(1, lab, fill = lab)) +
  geom_tile(width = 0.4, height = 0.6) +
  scale_fill_manual(values = PLDDT, guide = "none") +
  scale_y_discrete(limits = rev(levels(key$lab)), position = "right") +
  labs(x = NULL, y = NULL, title = "pLDDT") +
  theme_void(base_size = 8.5) +
  theme(axis.text.y.right = element_text(hjust = 0, size = 7.4,
                                         colour = "grey20",
                                         margin = margin(l = 3)),
        plot.title = element_text(size = 7.6, colour = "grey20", hjust = 0),
        plot.margin = margin(2, 2, 2, 2))

fig <- (pA | (pB / (pC + inset_element(p_key, left = 0, bottom = 0.6,
                                       right = 0.2, top = 1)) +
                plot_layout(heights = c(1, 1.15)))) +
  plot_layout(widths = c(1, 1.15))

ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_allele_swaps_full.pdf"), fig,
       width = 13.4, height = 6.6, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_allele_swaps_full.png"), fig,
       width = 13.4, height = 6.6, dpi = 300, bg = "white")
msg("wrote SUPP_FIG_XX_sid2_allele_swaps_full.{pdf,png}")
