## Supplement -- what the 37 kb NIL interval contains -------------------------
##
##   Rscript scripts/SUPP_FIG_XX_nil_interval_genes.R
##     -> plots/SUPP_FIG_XX_nil_interval_genes.{pdf,png}
##
## Gene models across 13.6577-13.6950 Mb of chromosome III, the interval the
## NIL series of Figure 3 resolves, with every site where the two parents of
## the cross differ drawn above them: protein-altering sites as lollipops,
## everything else as a rug so the denominator is visible.
##
## WHY THE RUG IS THERE. "Two missense changes in sid-2" is only an argument
## if the reader can see it is two OF something. The parents differ at 27 sites
## across these 37 kb; 2 are protein-altering, 3 synonymous, 22 intronic, UTR
## or unannotated. Drawing only the lollipops would let the panel be read as a
## search that found two hits rather than a census that found two.
##
## Colour is which parent carries the alternate allele, matching Figure 3's
## palette. Impact class is the lollipop head's outline, and the count of
## HIGH-impact differences -- zero -- is stated on the panel rather than left
## to be inferred from an absence.
##
## Runs from a clone: reads only the three staged tables. They are extracted
## from the BCSQ-annotated CeNDR VCF and a WormBase GFF3 by
## scripts/make_nil_interval_tables.R, which documents the consequence-to-
## impact mapping.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse); library(ggtext)
})

MAP <- "supplemental_data/mapping"
OUT <- "plots"
FROM <- 13657700; TO <- 13695000
COL_JU1793 <- "#F34C00"; COL_JU2466 <- "#40B4AB"; COL_REGION <- "#4D4D4D"
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

genes <- read_tsv(file.path(MAP, "nil_interval_genes.tsv"), show_col_types = FALSE)
exons <- read_tsv(file.path(MAP, "nil_interval_exons.tsv"), show_col_types = FALSE)
var   <- read_tsv(file.path(MAP, "nil_interval_parent_variants.tsv"), show_col_types = FALSE) %>%
  mutate(impact = factor(impact, levels = c("HIGH","MODERATE","LOW","MODIFIER")),
         ## 0/0 in JU1793 means JU2466 carries the alternate allele
         alt.parent = ifelse(gt.JU1793 == "0/0", "JU2466", "JU1793"))
stopifnot(nrow(genes) > 0, nrow(var) > 0)

## ---------------------------------------------------------------------------
## pack genes into rows so nothing overlaps -- sid-2 sits inside dyf-2 on the
## opposite strand, so a single row is not an option
## ---------------------------------------------------------------------------
## The six non-coding genes here are 52-260 bp -- at 37 kb across a page they
## are a pixel wide, and their names collide into illegibility. They get their
## own thin row of ticks, counted on the panel and named in the caption; the
## protein-coding genes get the labelled rows.
PAD <- 700
coding <- genes %>% filter(biotype == "protein_coding") %>% arrange(start) %>%
  mutate(row = NA_integer_)
ends <- numeric(0)
for (i in seq_len(nrow(coding))) {
  free <- which(ends < coding$start[i] - PAD)
  r <- if (length(free)) min(free) else length(ends) + 1L
  coding$row[i] <- r; ends[r] <- coding$end[i]
}
NROW_G <- max(coding$row)
g <- coding %>% mutate(y = -row)
nc <- genes %>% filter(biotype != "protein_coding") %>%
  mutate(y = -(NROW_G + 0.62))
exons <- exons %>% inner_join(g %>% select(wbgene, y), by = "wbgene")
msg("coding genes ", nrow(g), " in ", NROW_G, " rows; non-coding ", nrow(nc))

GH <- 0.26                                   # half-height of an exon box
RUG_Y  <- 0.30
LOLL_0 <- 0.75
LOLL_1 <- 1.75

lol <- var %>% filter(impact %in% c("HIGH","MODERATE")) %>%
  mutate(lab = sprintf("*%s* %s", gene, gsub(">", "&rarr;", aa.change)))
rug <- var %>% filter(!impact %in% c("HIGH","MODERATE"))
n_high <- sum(var$impact == "HIGH")

## stagger lollipop labels if two sit close together
lol <- lol %>% arrange(pos) %>%
  mutate(lab.y = LOLL_1 + 0.30 + 0.34 * (seq_len(n()) %% 2))

p <- ggplot() +
  ## the interval itself
  annotate("rect", xmin = FROM/1e6, xmax = TO/1e6,
           ymin = -NROW_G - 1.20, ymax = LOLL_1 + 1.15,
           fill = COL_REGION, alpha = 0.05) +
  ## the non-coding genes, as ticks on their own row
  geom_rect(data = nc, aes(xmin = start/1e6, xmax = end/1e6,
                           ymin = y - 0.13, ymax = y + 0.13),
            fill = "grey60", colour = "grey40", linewidth = 0.2) +
  annotate("richtext", x = FROM/1e6, y = -NROW_G - 1.00,
           label = sprintf("%d non-coding genes, 52&ndash;260 bp (%s)", nrow(nc),
                           paste(sort(nc$label), collapse = ", ")),
           size = 2.3, colour = "grey45", hjust = 0, vjust = 0.5,
           fill = NA, label.color = NA, label.padding = grid::unit(rep(0, 4), "pt")) +
  ## gene bodies: a thin line for the span, boxes for the exons
  geom_segment(data = g, aes(x = start/1e6, xend = end/1e6, y = y, yend = y),
               linewidth = 0.35, colour = "grey45") +
  geom_rect(data = exons, aes(xmin = start/1e6, xmax = end/1e6,
                              ymin = y - GH, ymax = y + GH),
            fill = "grey82", colour = "grey35", linewidth = 0.25) +
  geom_richtext(data = g, aes(x = (start + end)/2e6, y = y - GH - 0.20,
                              label = sprintf("*%s* %s", label, ifelse(strand == "+", "&rarr;", "&larr;"))),
                size = 2.5, colour = "grey20", vjust = 1, fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  ## the rug: every other site the parents differ at
  geom_segment(data = rug, aes(x = pos/1e6, xend = pos/1e6, y = RUG_Y - 0.13,
                               yend = RUG_Y + 0.13, colour = alt.parent),
               linewidth = 0.5) +
  annotate("richtext", x = FROM/1e6, y = RUG_Y + 0.20,
           label = sprintf("%d further differences &mdash; %d synonymous, %d intronic, UTR or unannotated",
                           nrow(rug), sum(var$impact == "LOW"), sum(var$impact == "MODIFIER")),
           size = 2.35, colour = "grey40", hjust = 0, vjust = 0,
           fill = NA, label.color = NA, label.padding = grid::unit(rep(0, 4), "pt")) +
  ## the lollipops
  geom_segment(data = lol, aes(x = pos/1e6, xend = pos/1e6, y = LOLL_0, yend = LOLL_1),
               linewidth = 0.4, colour = "grey35") +
  geom_point(data = lol, aes(x = pos/1e6, y = LOLL_1, fill = alt.parent, shape = impact),
             size = 3.1, stroke = 0.6, colour = "grey15") +
  geom_richtext(data = lol, aes(x = pos/1e6, y = lab.y, label = lab),
                size = 2.5, colour = "grey15", vjust = 0,
                fill = alpha("white", 0.8), label.color = NA,
                label.padding = grid::unit(rep(1, 4), "pt")) +
  annotate("richtext", x = TO/1e6, y = LOLL_1 + 1.05,
           label = sprintf("**%d HIGH-impact difference%s**", n_high, ifelse(n_high == 1, "", "s")),
           size = 2.6, colour = if (n_high == 0) "grey35" else "#9E4257",
           hjust = 1, vjust = 1, fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt")) +
  ## showing an empty HIGH key invites the reader to hunt for a class that is
  ## not there; the annotation states the zero instead
  scale_shape_manual(values = c(HIGH = 23, MODERATE = 21), drop = FALSE,
                     name = "Impact", limits = c("HIGH","MODERATE"),
                     guide = if (n_high > 0) "legend" else "none") +
  scale_fill_manual(values = c(JU1793 = COL_JU1793, JU2466 = COL_JU2466),
                    name = "Alternate allele carried by") +
  scale_colour_manual(values = c(JU1793 = COL_JU1793, JU2466 = COL_JU2466), guide = "none") +
  scale_x_continuous(labels = function(x) sprintf("%.3f", x),
                     breaks = scales::pretty_breaks(6),
                     expand = expansion(mult = 0.012)) +
  scale_y_continuous(breaks = NULL, expand = expansion(mult = 0.02)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3.1), order = 1)) +
  labs(x = "Chromosome III (Mb)", y = NULL,
       title = "The 37 kb interval the NIL series resolves",
       subtitle = paste0("Gene models, with every site where JU1793 and JU2466 differ drawn above them. ",
                         "Of ", nrow(var), " differences, ", nrow(lol), " alter a protein.")) +
  theme_classic(base_size = 11.5) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.3),
        axis.ticks.x = element_line(linewidth = 0.3),
        plot.title = element_markdown(size = 12),
        plot.subtitle = element_markdown(size = 8.6, colour = "grey35"),
        plot.title.position = "plot",
        legend.position = "bottom", legend.box = "horizontal",
        legend.key.size = grid::unit(9, "pt"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8.4))

ggsave(file.path(OUT, "SUPP_FIG_XX_nil_interval_genes.pdf"), p,
       width = 9.2, height = 4.4, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_nil_interval_genes.png"), p,
       width = 9.2, height = 4.4, dpi = 300, bg = "white")
msg("wrote SUPP_FIG_XX_nil_interval_genes.{pdf,png}")

cat("\n== protein-coding genes in the interval ==\n")
print(as.data.frame(g %>%
        transmute(gene = label, strand, start, end, kb = round((end - start)/1000, 1))),
      row.names = FALSE)
cat("\n== non-coding genes (drawn as ticks, named in the caption) ==\n")
print(as.data.frame(nc %>% transmute(gene = label, biotype, strand,
                                     bp = end - start)), row.names = FALSE)
cat("\n== the protein-altering differences ==\n")
print(as.data.frame(lol %>% select(pos, ref, alt, gene, aa.change, alt.parent)), row.names = FALSE)
