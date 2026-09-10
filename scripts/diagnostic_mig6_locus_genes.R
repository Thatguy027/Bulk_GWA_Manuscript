## Diagnostic -- gene content and parental variation at the top mig-6 QTL ----
##
##   Rscript scripts/diagnostic_mig6_locus_genes.R
##     -> plots/diagnostics/DIAG_mig6_locus_<chrom>_<Mb>.{pdf,png}
##        plots/diagnostics/TABLE_mig6_locus_census.tsv
##
## The census Figure S18 runs on the 37 kb NIL interval, applied to each
## independent HT115-vs-mig-6 cross QTL above LOD 500 in a 100 kb window on
## the peak. Tables come from scripts/make_mig6_locus_tables.R; this reads only
## those, so it builds from a clone.
##
## COLOUR IS THE PARENT CARRYING THE ALTERNATE ALLELE, and the two crosses share
## one pair of colours rather than taking four: parent 1 (N2, JU1793) pink,
## parent 2 (XZ1516, JU2466) green. p1/p2 is the same ordering the cross allele
## frequency tables use, so a reader who has seen those does not have to relearn
## it here.
##
## WHY ONLY SOME GENES ARE LABELLED. These windows hold 16-31 protein-coding
## genes. Labelling all of them at 100 kb across a page produces a band of
## overlapping italics, so only genes carrying a protein-altering difference are
## named -- those are the ones the panel exists to point at. The rest are drawn
## and counted.
##
## DIVERGENT REGIONS are shaded where they occur. In these four windows none do,
## which is itself worth seeing: the thing that actually limits the census here
## is per-parent missingness, so the no-call count is printed on each panel.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse); library(ggtext)})

MAP <- "supplemental_data/mapping"
OUT <- "plots/diagnostics"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
COL_P1 <- "#E15A97"   # N2, JU1793
COL_P2 <- "#2F8F5B"   # XZ1516, JU2466
COL_DIV <- "#A85B18"
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

genes <- read_tsv(file.path(MAP,"mig6_locus_genes.tsv"), show_col_types = FALSE)
exons <- read_tsv(file.path(MAP,"mig6_locus_exons.tsv"), show_col_types = FALSE)
var   <- read_tsv(file.path(MAP,"mig6_locus_variants.tsv"), show_col_types = FALSE)
div   <- read_tsv(file.path(MAP,"mig6_locus_divergent.tsv"), show_col_types = FALSE)
summ  <- read_tsv(file.path(MAP,"mig6_locus_summary.tsv"), show_col_types = FALSE)
P1 <- c(N2xXZ1516 = "N2", JU1793xJU2466 = "JU1793")
P2 <- c(N2xXZ1516 = "XZ1516", JU1793xJU2466 = "JU2466")

## an inframe indel has no single-residue field in BCSQ, so name the class
var <- var %>% mutate(
  change = ifelse(is.na(aa.change), consequence, gsub(">", "&rarr;", aa.change)),
  parent.rank = ifelse(alt.parent == P1[cross], "p1", "p2"))

one_panel <- function(L) {
  s  <- summ  %>% filter(locus == L)
  g0 <- genes %>% filter(locus == L)
  ex <- exons %>% filter(locus == L)
  vv <- var   %>% filter(locus == L)
  dd <- div   %>% filter(locus == L)
  FROM <- s$from; TO <- s$to
  p1 <- P1[[s$cross]]; p2 <- P2[[s$cross]]
  pal <- setNames(c(COL_P1, COL_P2), c(p1, p2))

  ## pack the coding genes; non-coding go on one tick row beneath
  coding <- g0 %>% filter(biotype == "protein_coding") %>% arrange(start) %>% mutate(row = NA_integer_)
  ends <- numeric(0)
  for (i in seq_len(nrow(coding))) {
    free <- which(ends < coding$start[i] - 1500)
    r <- if (length(free)) min(free) else length(ends) + 1L
    coding$row[i] <- r; ends[r] <- coding$end[i]
  }
  NR <- max(coding$row, 1)
  cg <- coding %>% mutate(y = -row)
  nc <- g0 %>% filter(biotype != "protein_coding") %>% mutate(y = -(NR + 0.62))
  ex <- ex %>% inner_join(cg %>% select(wbgene, y), by = "wbgene")

  hits <- vv %>% filter(impact %in% c("HIGH","MODERATE"))
  named <- cg %>% semi_join(hits %>% distinct(gene), by = c("label" = "gene"))
  rug   <- vv %>% filter(!impact %in% c("HIGH","MODERATE"))
  GH <- 0.26; RUG_Y <- 0.30; L0 <- 0.75; L1 <- 1.75
  ## One label per GENE, not per variant. Labelling each change individually
  ## overlaps into illegibility wherever a gene carries several -- trpp-10 has
  ## four -- and the per-variant detail is in TABLE_mig6_locus_census.tsv and
  ## the staged variant table anyway.
  glab <- hits %>% group_by(gene) %>%
    summarise(x = mean(pos), n = n(),
              one = first(change[!grepl("_", change)]), .groups = "drop") %>%
    arrange(x) %>%
    mutate(txt = ifelse(n > 1, sprintf("*%s* &times;%d", gene, n),
                        sprintf("*%s* %s", gene, coalesce(one, ""))),
           lab.y = L1 + 0.30 + 0.34 * (seq_len(n()) %% 3))

  p <- ggplot() +
    annotate("rect", xmin = FROM/1e6, xmax = TO/1e6,
             ymin = -NR - 1.25, ymax = L1 + 1.45, fill = "#4D4D4D", alpha = 0.05)
  if (nrow(dd))
    p <- p + geom_rect(data = dd, aes(xmin = start/1e6, xmax = end/1e6,
                                      ymin = -NR - 1.25, ymax = L1 + 1.45),
                       fill = COL_DIV, alpha = 0.13, inherit.aes = FALSE)
  p <- p +
    annotate("segment", x = s$peak.bp/1e6, xend = s$peak.bp/1e6,
             y = -NR - 1.25, yend = L1 + 1.15, linetype = "dashed",
             linewidth = 0.4, colour = "grey40") +
    annotate("richtext", x = s$peak.bp/1e6, y = L1 + 1.30, label = "peak",
             size = 2.3, colour = "grey35", hjust = 0.5, vjust = 0,
             fill = NA, label.color = NA, label.padding = grid::unit(rep(0,4),"pt")) +
    geom_segment(data = cg, aes(x = start/1e6, xend = end/1e6, y = y, yend = y),
                 linewidth = 0.3, colour = "grey55") +
    geom_rect(data = ex, aes(xmin = start/1e6, xmax = end/1e6,
                             ymin = y - GH, ymax = y + GH),
              fill = "grey84", colour = "grey40", linewidth = 0.2) +
    geom_rect(data = nc, aes(xmin = start/1e6, xmax = end/1e6,
                             ymin = y - 0.12, ymax = y + 0.12),
              fill = "grey65", colour = NA) +
    geom_richtext(data = named, aes(x = (start+end)/2e6, y = y - GH - 0.16,
                                    label = sprintf("*%s*", label)),
                  size = 2.25, colour = "grey15", vjust = 1,
                  fill = alpha("white", 0.75), label.color = NA,
                  label.padding = grid::unit(rep(0.5,4),"pt")) +
    geom_segment(data = rug, aes(x = pos/1e6, xend = pos/1e6,
                                 y = RUG_Y - 0.12, yend = RUG_Y + 0.12,
                                 colour = alt.parent), linewidth = 0.35) +
    geom_segment(data = hits, aes(x = pos/1e6, xend = pos/1e6, y = L0, yend = L1),
                 linewidth = 0.35, colour = "grey40") +
    geom_point(data = hits, aes(x = pos/1e6, y = L1, fill = alt.parent),
               shape = 21, size = 2.5, stroke = 0.5, colour = "grey20") +
    geom_richtext(data = glab, aes(x = x/1e6, y = lab.y, label = txt),
                  size = 2.15, colour = "grey20", vjust = 0, hjust = 0.5,
                  fill = alpha("white", 0.8), label.color = NA,
                  label.padding = grid::unit(rep(0.6,4),"pt")) +
    annotate("richtext", x = FROM/1e6, y = RUG_Y + 0.19,
             label = sprintf("%d further differences", nrow(rug)),
             size = 2.2, colour = "grey45", hjust = 0, vjust = 0,
             fill = NA, label.color = NA, label.padding = grid::unit(rep(0,4),"pt")) +
    annotate("richtext", x = FROM/1e6, y = -NR - 1.08,
             label = sprintf("%d non-coding genes%s", nrow(nc),
                             if (nrow(dd)) sprintf(" &middot; <span style='color:%s'>shaded: divergent region</span>", COL_DIV) else ""),
             size = 2.1, colour = "grey50", hjust = 0, vjust = 0.5,
             fill = NA, label.color = NA, label.padding = grid::unit(rep(0,4),"pt")) +
    scale_fill_manual(values = pal, name = "Alternate allele in", drop = FALSE) +
    scale_colour_manual(values = pal, guide = "none") +
    scale_x_continuous(labels = function(x) sprintf("%.2f", x),
                       breaks = scales::pretty_breaks(6),
                       expand = expansion(mult = 0.012)) +
    scale_y_continuous(breaks = NULL, expand = expansion(mult = 0.02)) +
    labs(x = sprintf("Chromosome %s (Mb)", s$chrom), y = NULL,
         title = sprintf("%s &mdash; LOD %.0f, &Delta;freq %+.2f",
                         gsub("_", " ", L), s$peak.LOD, s$dfreq),
         subtitle = sprintf(paste("%d protein-coding genes &middot; %d of %d sites differ &middot;",
                                  "**%d protein-altering, %d HIGH** &middot;",
                                  "no-call %s %d / %s %d of %d"),
                            s$coding, s$differ, s$sites, nrow(hits), s$HIGH,
                            p1, s$nocall.p1, p2, s$nocall.p2, s$sites)) +
    theme_classic(base_size = 10.5) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.x = element_line(linewidth = 0.3),
          plot.title = element_markdown(size = 11),
          plot.subtitle = element_markdown(size = 7.6, colour = "grey35"),
          plot.title.position = "plot",
          legend.position = "bottom", legend.key.size = grid::unit(8,"pt"),
          legend.text = element_text(size = 7.6), legend.title = element_text(size = 7.8))

  stem <- sprintf("DIAG_mig6_locus_%s_%s", s$chrom, sub("\\.", "-", sprintf("%.2f", s$peak.Mb)))
  h <- 2.6 + 0.42 * NR
  ggsave(file.path(OUT, paste0(stem, ".pdf")), p, width = 9.2, height = h, device = cairo_pdf)
  ggsave(file.path(OUT, paste0(stem, ".png")), p, width = 9.2, height = h, dpi = 300, bg = "white")
  msg("wrote ", stem, " (", NR, " gene rows)")
  tibble(locus = L, file = paste0(stem, ".png"), gene.rows = NR,
         labelled = nrow(named), protein.altering = nrow(hits))
}

## the summary table already carries parent names and per-parent no-call counts
stopifnot(all(c("nocall.p1","nocall.p2","sites","from","to","peak.bp") %in% names(summ)))

built <- bind_rows(lapply(summ$locus[order(-summ$peak.LOD)], one_panel))
write_tsv(summ %>% left_join(built, by = "locus"),
          file.path(OUT, "TABLE_mig6_locus_census.tsv"))
cat("\n== built ==\n"); print(as.data.frame(built), row.names = FALSE)
