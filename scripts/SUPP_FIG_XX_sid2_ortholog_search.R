## Supplement -- how far out SID-2 can be compared, and what that is worth ----
##
##   python3 scripts/make_sid2_ortholog_tables.py   # once, stages the tables
##   Rscript scripts/SUPP_FIG_XX_sid2_ortholog_search.R
##     -> plots/SUPP_FIG_XX_sid2_ortholog_search.{pdf,png}
##
##   A  reciprocal-best-hit blastp of C. elegans SID-2 against 20 nematode
##      proteomes spanning the phylum, banded by how far out the species sits
##   B  per-position conservation over the Elegans group, placing residues 94,
##      95 and 96 against the rest of the ectodomain
##
## This is the search and the calibration behind
## SUPP_FIG_XX_sid2_ortholog_conservation, which shows what residue 96 does
## across the species this search finds.
##
## PANEL A: THE CONSERVATION DOES NOT FADE, IT STOPS AT THE GENUS
## Eight of the 19 comparator proteomes clear E < 1e-5 and every one of them is
## a Caenorhabditis. Nothing outside the genus comes close, including
## Diploscapter pachys, the sister genus, at E = 4.1. Two independent resources
## agree: the UniRef50 cluster containing G5EEV9 has exactly one member, and
## NCBI's ortholog set for sid-2 within Nematoda is empty. So the comparison in
## the companion figure is bounded by the genus, and that bound is the reason
## this figure exists: a reader should not take "conserved" to mean more than
## the data can carry.
##
## Two Caenorhabditis proteomes also fail (C. bovis E = 2.2, C. auriculariae
## E = 0.15). Absence of a hit in one proteome is weak evidence about the gene
## and may only mean its annotation is incomplete, so no gene loss is claimed.
##
## PANEL B: WHY ONE CONSERVED COLUMN WOULD NOT HAVE BEEN ENOUGH
## Thr is 13.3% of this ectodomain, so a single conserved Thr in a 42-48%
## identity alignment is weak evidence. The three-residue window controls for
## that, because all three columns come from the same alignment:
##
##   N94  the constrained Asn of the sequon   6/6 conserved   92nd percentile
##   C95  the unconstrained X of the sequon   0/6 conserved   16th percentile
##   T96  the constrained Ser/Thr             5/6, 6/6 as ST  82nd percentile
##
## against an ectodomain background of 2.23/6 (37.1%), with only 15.0% of
## positions conserved in all six. A misaligned window does not produce
## 100% / 0% / 83% across three adjacent columns.
##
## AND THE CALIBRATION THAT CUTS THE OTHER WAY: three of the five ectodomain
## sequons in C. elegans SID-2 are intact in all six orthologs, so a fully
## conserved sequon is the norm in this protein and the one at 94 is a member
## of the conserved majority rather than a standout.
##
## The conservation statistics here are over the six FULL-LENGTH UniProt
## Elegans-group orthologs, which are the only ones with an alignment good
## enough to score every position. The companion figure's wider species set is
## scored only at the three-residue window, where an HSP can be checked
## directly.
##
## Runs from a clone: reads only the staged tables.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggtext)
})

ST  <- "supplemental_data/structure"
OUT <- "plots"
E_ORTH <- 1e-5
ECD <- c(21, 193)

COL_ORTH <- "#2E4057"    # cleared the ortholog threshold
COL_NO   <- "#B8C2CA"    # did not
COL_FOC  <- "#9E4257"    # the sequon positions

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
## proper ordinals: the same helper SUPP_FIG_XX_sid2_local_charge.R uses, after
## that figure shipped "82th"
ord <- function(n) vapply(n, function(k) {
  k <- round(k)
  suff <- if (k %% 100 %in% 11:13) "th"
          else switch(as.character(k %% 10), "1" = "st", "2" = "nd",
                      "3" = "rd", "th")
  paste0(k, suff)
}, character(1))
panel_title <- function(letter, txt = NULL) {
  lt <- paste0("<span style='font-size:13pt;color:#111111'>**", letter, "**</span>")
  if (is.null(txt)) lt else paste0(lt, " ", txt)
}
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_markdown(size = base_size + 0.5),
          plot.subtitle = element_markdown(size = base_size - 2.5, colour = "grey30"),
          plot.title.position = "plot",
          legend.key.size = grid::unit(8, "pt"))
}
## plain text only: strwrap splits HTML tags
wrap_md <- function(txt, width = 72)
  paste(strwrap(txt, width = width), collapse = "<br>")
ital <- function(txt) {
  for (g in c("Caenorhabditis", "elegans", "Diploscapter", "pachys", "sid-2"))
    txt <- gsub(g, paste0("*", g, "*"), txt, fixed = TRUE)
  gsub("\\*\\*", "", txt)
}

srch <- read_tsv(file.path(ST, "sid2_ortholog_search.tsv"), show_col_types = FALSE)
cons <- read_tsv(file.path(ST, "sid2_ortholog_conservation.tsv"), show_col_types = FALSE)
stopifnot(nrow(srch) == 20, nrow(cons) == 311)
N_ORTH <- unique(cons$n_orthologs)
stopifnot(length(N_ORTH) == 1)

## ===========================================================================
## A -- how far out a SID-2 ortholog is detectable
## ===========================================================================
## the query itself is not plotted: its E is 0 and it is the reference
BANDS <- c("Elegans group", "Japonica group", "basal *Caenorhabditis*",
           "outside *Caenorhabditis*")
sa <- srch %>%
  filter(species != "Caenorhabditis elegans") %>%
  mutate(nlp = -log10(pmax(evalue, 1e-300)),
         band = factor(case_when(depth_rank == 1 ~ BANDS[1],
                                 depth_rank == 2 ~ BANDS[2],
                                 depth_rank == 3 ~ BANDS[3],
                                 TRUE            ~ BANDS[4]), levels = BANDS),
         lab = paste0("*", species, "*  <span style='color:grey45'>",
                      sprintf("%.0f%% id, %.0f%% cov", percent_identity,
                              query_coverage), "</span>"),
         lab = fct_reorder(lab, depth_rank * 1000 - nlp, .desc = TRUE))
msg("panel A: ", sum(sa$is_ortholog), " of ", nrow(sa),
    " comparator proteomes yield an ortholog at E < ", E_ORTH)

pA <- ggplot(sa, aes(nlp, lab, colour = is_ortholog)) +
  geom_vline(xintercept = -log10(E_ORTH), linetype = "dashed",
             linewidth = 0.4, colour = "grey45") +
  geom_segment(aes(x = 0, xend = nlp, yend = lab), linewidth = 0.5) +
  geom_point(size = 2.1) +
  ## the threshold is named in the axis title, not annotated in the panel:
  ## annotate() draws once per facet, so it appeared four times and two of
  ## them landed on the lollipops
  scale_colour_manual(values = c(`TRUE` = COL_ORTH, `FALSE` = COL_NO),
                      guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.08))) +
  facet_grid(band ~ ., scales = "free_y", space = "free_y") +
  labs(x = paste0("&minus;log<sub>10</sub> *E* of the reciprocal best hit ",
                  "<span style='color:grey45'>(dashed: *E* = 10<sup>&minus;5</sup>)</span>"),
       y = NULL,
       title = panel_title("A", "**Outside *Caenorhabditis*, SID-2 has no detectable ortholog**"),
       subtitle = ital(wrap_md(paste0(
         "Reciprocal-best-hit blastp of C. elegans SID-2 against 20 UniProt ",
         "reference proteomes, 488,718 proteins. The query itself is not ",
         "drawn, so 19 comparators are: eight clear E < 1e-5 and every one of ",
         "them is a Caenorhabditis. Nothing beyond the genus clears it, ",
         "including Diploscapter pachys, the sister genus, at E = 4.1. Two ",
         "Caenorhabditis proteomes also fail, which is weak evidence about ",
         "the gene and may only say their annotation is incomplete."), 96))) +
  theme_pub(10) +
  theme(axis.text.y = element_markdown(size = 7.6),
        axis.title.x = element_markdown(size = 9),
        strip.text.y = element_markdown(size = 7, angle = 0, hjust = 0),
        strip.background = element_rect(fill = "grey96", colour = NA),
        panel.spacing.y = grid::unit(3, "pt"),
        plot.subtitle = element_markdown(size = 7.2, colour = "grey30",
                                         lineheight = 1.3))

## ===========================================================================
## B -- what conservation at 94, 95 and 96 is worth against the background
## ===========================================================================
ec <- cons %>% filter(in_ectodomain)
dist <- ec %>% count(n_conserved)
MARK <- ec %>% filter(ce_pos %in% c(94, 95, 96)) %>%
  mutate(lab = paste0(ce_aa, ce_pos, "<br><span style='font-size:6pt'>",
                      ord(ecd_percentile), " pct</span>"))
msg("panel B: ectodomain mean ",
    sprintf("%.2f/%d", mean(ec$n_conserved), N_ORTH), " (",
    sprintf("%.1f%%", 100 * mean(ec$n_conserved) / N_ORTH), "), ",
    sum(ec$n_conserved == N_ORTH), " of ", nrow(ec), " conserved in all")

pB <- ggplot(dist, aes(factor(n_conserved), n)) +
  geom_col(fill = COL_NO, width = 0.72) +
  geom_text(aes(label = n), vjust = -0.5, size = 2.6, colour = "grey35") +
  geom_richtext(data = MARK %>% left_join(dist, by = "n_conserved"),
                aes(x = factor(n_conserved), y = n, label = lab),
                inherit.aes = FALSE, vjust = -0.55, size = 2.5,
                colour = COL_FOC, fill = "white", label.color = NA,
                label.padding = grid::unit(c(1, 2, 1, 2), "pt")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.42))) +
  labs(x = paste0("Elegans-group orthologs sharing the *C. elegans* residue (of ",
                  N_ORTH, ")"),
       y = "Ectodomain positions",
       title = panel_title("B", "**Constrained positions high, free one lowest**"),
       subtitle = ital(wrap_md(sprintf(paste0(
         "All %d ectodomain positions (%d-%d), binned by how many of the %d ",
         "full-length Elegans-group orthologs share the C. elegans residue. ",
         "The background is %.2f of %d (%.0f%%) and only %.0f%% of positions ",
         "are conserved in all. Calibration that cuts the other way: three of ",
         "the five ectodomain sequons are intact in all %d orthologs, so a ",
         "fully conserved sequon is the norm in this protein and the one at 94 ",
         "is not a standout."),
         nrow(ec), ECD[1], ECD[2], N_ORTH, mean(ec$n_conserved), N_ORTH,
         100 * mean(ec$n_conserved) / N_ORTH,
         100 * sum(ec$n_conserved == N_ORTH) / nrow(ec), N_ORTH), 86))) +
  theme_pub(10) +
  ## axis.title = element_markdown() does NOT render markdown; the per-axis
  ## elements do. Verified: it printed the asterisks literally.
  theme(axis.title.x = element_markdown(size = 8.6),
        axis.title.y = element_markdown(size = 9),
        plot.subtitle = element_markdown(size = 7.2, colour = "grey30",
                                         lineheight = 1.3))

## ===========================================================================
fig <- pA / pB + plot_layout(heights = c(1.5, 1))

ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_ortholog_search.pdf"), fig,
       width = 9.6, height = 9.4, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_ortholog_search.png"), fig,
       width = 9.6, height = 9.4, dpi = 300, bg = "white")
msg("wrote SUPP_FIG_XX_sid2_ortholog_search.{pdf,png}")
