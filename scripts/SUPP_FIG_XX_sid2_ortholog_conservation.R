## Supplement -- how far out SID-2 can be compared, and what 96 does ---------
##
##   python3 scripts/make_sid2_ortholog_tables.py   # once, stages the tables
##   Rscript scripts/SUPP_FIG_XX_sid2_ortholog_conservation.R
##     -> plots/SUPP_FIG_XX_sid2_ortholog_conservation.{pdf,png}
##
##   A  reciprocal-best-hit blastp of C. elegans SID-2 against twenty nematode
##      proteomes, ordered by how far out the species sits, with the E < 1e-5
##      line that an ortholog call had to clear
##   B  the residues aligned to C. elegans 88-104 in every ortholog, with the
##      N94-x-T96 sequon columns marked
##   C  per-position conservation over the Elegans group, with 94, 95 and 96
##      placed against the rest of the ectodomain
##
## THE QUESTION. sid-2 96K is a polymorphism inside C. elegans, so calling 96T
## ancestral needs outgroups. Asking for more of them turns into the opposite
## question -- how far out does the comparison still work at all -- and that has
## a sharper answer than expected.
##
## THREE BREAKDOWN POINTS, NESTED
##   the sequon        breaks at the base of the genus. N94-x-[ST] is intact in
##                     all six Elegans-group orthologs; C. angaria, basal
##                     Angaria group, has Q-G-F there and no sequon.
##   the column        stops being readable one step earlier. C. japonica's
##                     ortholog is split across two proteome entries and the
##                     N-terminal fragment aligns into a TTDT repeat with a
##                     two-to-three residue offset, so that column is drawn as
##                     AMBIGUOUS, not as a substitution.
##   orthology itself  stops immediately outside Caenorhabditis. Nothing beyond
##                     the genus clears E < 1e-5, including Diploscapter
##                     pachys, the sister genus, at E = 4.1. UniRef50 puts
##                     G5EEV9 in a cluster of one and NCBI lists no nematode
##                     ortholog, so three resources agree.
##
## SO THE CLAIM HAS A DEPTH. 96T is ancestral at the depth of the Elegans group
## and 96K is derived within C. elegans. This figure is built so it cannot be
## read as "conserved across nematodes", because it is not.
##
## WHY PANEL B IS THE ARGUMENT AND PANEL C IS ITS CALIBRATION
## One conserved column in a 42-48% identity alignment of a Thr-rich region is
## weak on its own -- Thr is 13.3% of this ectodomain. The three-residue window
## controls for that, because all three columns come from the same alignment:
## the sequon's two constrained positions are conserved and its unconstrained
## middle is not, which a bad alignment does not produce.
##
## AND THE CALIBRATION THAT CUTS THE OTHER WAY, which panel C states: three of
## the five ectodomain sequons in C. elegans SID-2 are intact in all six
## orthologs. A fully conserved sequon is the norm in this protein, so the one
## at 94 is a member of the conserved majority rather than a standout.
##
## NO FUNCTIONAL CLAIM. The glycosylation hypothesis for this sequon was tested
## with N94A and failed -- see SUPP_FIG_XX_sid2_allele_swaps_full.R. What is
## shown here is history, not mechanism.
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
WIN <- 88:104            # the window panel B draws
ECD <- c(21, 193)

COL_ORTH <- "#2E4057"    # a hit that cleared the ortholog threshold
COL_NO   <- "#B8C2CA"    # one that did not
COL_KEEP <- "#1A7F5A"    # sequon intact
COL_LOST <- "#B03A2E"    # sequon broken
COL_AMB  <- "#C9C2B4"    # column not readable
COL_ID   <- "#7E9BB5"    # matches C. elegans
COL_DIFF <- "#EDE7DC"    # does not
COL_FOC  <- "#9E4257"    # the sequon columns

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
## plain text only through here: strwrap splits HTML tags
wrap_md <- function(txt, width = 72)
  paste(strwrap(txt, width = width), collapse = "<br>")
ital <- function(txt) {
  for (g in c("Caenorhabditis", "elegans", "briggsae", "nigoni", "tropicalis",
              "remanei", "latens", "brenneri", "japonica", "angaria", "bovis",
              "auriculariae", "Diploscapter", "pachys", "sid-2"))
    txt <- gsub(g, paste0("*", g, "*"), txt, fixed = TRUE)
  gsub("\\*\\*", "", txt)
}

srch <- read_tsv(file.path(ST, "sid2_ortholog_search.tsv"), show_col_types = FALSE)
aln  <- read_tsv(file.path(ST, "sid2_ortholog_alignment.tsv"), show_col_types = FALSE)
cons <- read_tsv(file.path(ST, "sid2_ortholog_conservation.tsv"), show_col_types = FALSE)
stopifnot(nrow(srch) == 20, nrow(cons) == 311, nrow(aln) == 311 * 8)

N_ORTH <- unique(cons$n_orthologs)
stopifnot(length(N_ORTH) == 1)
msg("conservation computed over ", N_ORTH, " Elegans-group orthologs")

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
    " comparator proteomes yield an ortholog")

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
         "Reciprocal-best-hit blastp of C. elegans SID-2 against twenty ",
         "UniProt reference proteomes, 488,718 proteins. The query itself is ",
         "not drawn, so 19 comparators are: eight clear E < 1e-5 and ",
         "every one of them is a Caenorhabditis. Nothing beyond the genus ",
         "clears it, including Diploscapter pachys, the sister genus, at ",
         "E = 4.1. Two Caenorhabditis proteomes also fail, which is weak ",
         "evidence about the gene and may only say their annotation is ",
         "incomplete."), 96))) +
  theme_pub(10) +
  theme(axis.text.y = element_markdown(size = 7.6),
        axis.title.x = element_markdown(size = 9),
        strip.text.y = element_markdown(size = 7, angle = 0, hjust = 0),
        strip.background = element_rect(fill = "grey96", colour = NA),
        panel.spacing.y = grid::unit(3, "pt"),
        plot.subtitle = element_markdown(size = 7.2, colour = "grey30",
                                         lineheight = 1.3))

## ===========================================================================
## B -- the residues aligned to the sequon window
## ===========================================================================
ORDER <- c("Caenorhabditis elegans", "Caenorhabditis briggsae",
           "Caenorhabditis nigoni", "Caenorhabditis tropicalis",
           "Caenorhabditis remanei", "Caenorhabditis latens",
           "Caenorhabditis brenneri", "Caenorhabditis japonica",
           "Caenorhabditis angaria")
ce_row <- cons %>% filter(ce_pos %in% WIN) %>%
  transmute(ce_pos, ce_aa, species = "Caenorhabditis elegans",
            aligned_aa = ce_aa, identical = TRUE, ambiguous = FALSE)
w <- aln %>% filter(ce_pos %in% WIN) %>%
  select(ce_pos, ce_aa, species, aligned_aa, identical, ambiguous) %>%
  bind_rows(ce_row) %>%
  mutate(species = factor(species, levels = rev(ORDER)),
         fill = case_when(ambiguous ~ "ambiguous",
                          identical ~ "identical",
                          TRUE ~ "different"),
         ## an unreadable row is drawn as unreadable: printing the residues a
         ## mis-anchored alignment happened to produce invites reading them
         shown = if_else(ambiguous, "?", aligned_aa))
stopifnot(!anyNA(w$species))

## is the sequon intact in each row, using that row's own aligned residues
sq <- w %>% filter(ce_pos %in% c(94, 95, 96)) %>%
  select(species, ce_pos, aligned_aa, ambiguous) %>%
  pivot_wider(names_from = ce_pos, values_from = aligned_aa,
              names_prefix = "p") %>%
  mutate(state = case_when(
    ambiguous ~ "ambiguous",
    p94 == "N" & p96 %in% c("S", "T") & p95 != "P" ~ "intact",
    TRUE ~ "broken"))
msg("panel B: sequon intact in ",
    sum(sq$state == "intact"), " of ", nrow(sq), " rows drawn (",
    sum(sq$state == "broken"), " broken, ", sum(sq$state == "ambiguous"),
    " ambiguous)")

pB <- ggplot(w, aes(factor(ce_pos), species)) +
  annotate("rect", xmin = which(WIN == 94) - 0.5, xmax = which(WIN == 96) + 0.5,
           ymin = -Inf, ymax = Inf, fill = COL_FOC, alpha = 0.10) +
  geom_tile(aes(fill = fill), colour = "white", linewidth = 0.9) +
  geom_text(aes(label = shown,
                colour = fill == "identical" & ce_pos %in% c(94, 96)),
            size = 2.9, fontface = "bold", show.legend = FALSE) +
  annotate("rect", xmin = which(WIN == 94) - 0.5, xmax = which(WIN == 96) + 0.5,
           ymin = 0.45, ymax = length(ORDER) + 0.55, fill = NA,
           colour = COL_FOC, linewidth = 0.7) +
  geom_point(data = sq, aes(x = length(WIN) + 1.1, y = species, shape = state),
             inherit.aes = FALSE, size = 2.2, colour = COL_ORTH) +
  scale_fill_manual(values = c(identical = COL_ID, different = COL_DIFF,
                               ambiguous = COL_AMB), name = NULL,
                    breaks = c("identical", "different", "ambiguous"),
                    labels = c("matches *C. elegans*", "differs",
                               "column not readable")) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey25")) +
  scale_shape_manual(values = c(intact = 16, broken = 4, ambiguous = 1),
                     name = "N94-x-[ST]",
                     breaks = c("intact", "broken", "ambiguous"),
                     labels = c("intact", "broken", "unreadable")) +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.8))) +
  labs(x = "*C. elegans* SID-2 residue", y = NULL,
       title = panel_title("B", "**The sequon holds across the Elegans group and breaks at the base of the genus**"),
       subtitle = ital(wrap_md(paste0(
         "Residues aligned to C. elegans 88-104. The shaded columns are the ",
         "N94-x-T96 sequon. Its two constrained positions hold while the ",
         "middle one varies freely; C. angaria has Q-G-F. C. japonica's ",
         "ortholog is split across two proteome entries and its N-terminal ",
         "fragment aligns into a TTDT repeat offset by two to three residues, ",
         "so that row is drawn as unreadable rather than as substitutions."), 96))) +
  theme_pub(10) +
  theme(axis.text.y = element_markdown(size = 7.4, face = "italic"),
        axis.text.x = element_text(size = 6.8),
        axis.title.x = element_markdown(size = 9),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_markdown(size = 6.8),
        legend.title = element_text(size = 7.2),
        plot.subtitle = element_markdown(size = 7.2, colour = "grey30",
                                         lineheight = 1.3))

## ===========================================================================
## C -- what conservation at 94, 95 and 96 is worth against the background
## ===========================================================================
ec <- cons %>% filter(in_ectodomain)
dist <- ec %>% count(n_conserved)
MARK <- ec %>% filter(ce_pos %in% c(94, 95, 96)) %>%
  mutate(lab = paste0(ce_aa, ce_pos, "<br><span style='font-size:6pt'>",
                      ord(ecd_percentile), " pct</span>"))
msg("panel C: ectodomain mean ",
    sprintf("%.2f/%d", mean(ec$n_conserved), N_ORTH), " (",
    sprintf("%.1f%%", 100 * mean(ec$n_conserved) / N_ORTH), "), ",
    sum(ec$n_conserved == N_ORTH), " of ", nrow(ec), " conserved in all")

pC <- ggplot(dist, aes(factor(n_conserved), n)) +
  geom_col(fill = COL_NO, width = 0.72) +
  geom_text(aes(label = n), vjust = -0.5, size = 2.6, colour = "grey35") +
  geom_richtext(data = MARK %>%
                  left_join(dist, by = "n_conserved"),
                aes(x = factor(n_conserved), y = n, label = lab),
                inherit.aes = FALSE, vjust = -0.55, size = 2.5,
                colour = COL_FOC, fill = "white", label.color = NA,
                label.padding = grid::unit(c(1, 2, 1, 2), "pt")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(x = paste0("Elegans-group orthologs sharing the *C. elegans* residue (of ",
                  N_ORTH, ")"),
       y = "Ectodomain positions",
       title = panel_title("C", "**Constrained positions high, free one lowest**"),
       subtitle = ital(wrap_md(sprintf(paste0(
         "All %d ectodomain positions (%d-%d), binned by how many of the %d ",
         "orthologs share the C. elegans residue. The background is %.2f of ",
         "%d (%.0f%%) and only %.0f%% of positions are conserved in all. ",
         "Calibration that cuts the other way: three of the five ectodomain ",
         "sequons are intact in all %d orthologs, so a fully conserved sequon ",
         "is the norm in this protein and the one at 94 is not a standout."),
         nrow(ec), ECD[1], ECD[2], N_ORTH, mean(ec$n_conserved), N_ORTH,
         100 * mean(ec$n_conserved) / N_ORTH,
         100 * sum(ec$n_conserved == N_ORTH) / nrow(ec), N_ORTH), 78))) +
  theme_pub(10) +
  ## axis.title = element_markdown() does NOT render markdown; the per-axis
  ## elements do. Verified: it printed the asterisks literally.
  theme(axis.title.x = element_markdown(size = 8.4),
        axis.title.y = element_markdown(size = 9),
        plot.subtitle = element_markdown(size = 7.2, colour = "grey30",
                                         lineheight = 1.3))

## ===========================================================================
fig <- pA / (pB | pC) + plot_layout(heights = c(1, 1.05))

ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_ortholog_conservation.pdf"), fig,
       width = 11.6, height = 9.2, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_ortholog_conservation.png"), fig,
       width = 11.6, height = 9.2, dpi = 300, bg = "white")
msg("wrote SUPP_FIG_XX_sid2_ortholog_conservation.{pdf,png}")
