## Supplement -- residue 96 across Caenorhabditis, and what removing it does --
##
##   python3 scripts/make_sid2_ortholog_tables.py   # once, stages the tables
##   Rscript scripts/SUPP_FIG_XX_sid2_ortholog_conservation.R
##     -> plots/SUPP_FIG_XX_sid2_ortholog_conservation.{pdf,png}
##
##   A  the residues aligned to C. elegans 88-104 in every Caenorhabditis
##      species where that window can be read, with the N94-x-T96 sequon marked
##   B  the editing series that asks whether the conserved sequon matters:
##      N94A removes it while leaving residue 96 alone
##
## The search behind panel A, and the calibration of what its conservation is
## worth, are in SUPP_FIG_XX_sid2_ortholog_search.
##
## THE QUESTION. sid-2 96K is a polymorphism inside C. elegans, so calling 96T
## ancestral needs outgroups. This is 14 of them.
##
## WHAT PANEL A SHOWS
##   residue 96 is Ser or Thr in 12 of the 14 species whose window can be
##   read -- 11 Thr, 1 Ser -- and LYSINE IN NONE OF THEM
##   the N94-x-[ST] sequon is intact in 11 of 14
##   of the three losses, only C. doughertyi changes residue 96 itself, to Ala;
##   C. afra keeps Thr96 and loses Asn94, and C. sp54 changes both
##
## So the constraint sits on residue 96 being small and hydroxylated rather than
## on the sequon as a unit. C. afra is a NATURAL AxT, which is the same
## construct panel B shows to be fully resistant.
##
## EVERY ROW CARRIES ITS OWN CONFIDENCE. The measure is the local identity of
## the +/-10 residue block around the window, and the floor is 37% because that
## is where the six UniProt Elegans-group orthologs already sat -- not a
## threshold picked to get an answer. Rows below it go in the lower facet and
## are counted neither way, which is where the two species that appear to carry
## Lys96 sit: at 32% and 21% block identity those calls are not trustworthy,
## and they are shown rather than hidden.
##
## WHY PANEL B IS IN A CONSERVATION FIGURE
## A conserved sequon invites the reading that the glycan matters. Panel B is
## the experiment that refuses it, on the same page rather than as a cross
## reference. N94A removes the sequon while leaving residue 96 alone, so it
## should phenocopy 96K if the glycan is what 96K destroys. It does not:
##
##   JU1793  NxT (wild type)  94.8%      JU2466  NxK (wild type)   5.4% / 3.7%
##   JU1793  AxT (N94A)       99.3%      JU2466  NxT (96T swap)   18.4%
##   JU1793  NxK (96K swap)   53.1%      JU2466  AxK (N94A)       42.5%
##
## The AxT construct has no sequon and is FULLY RESISTANT, so the glycan is not
## required. Losing the sequon by removing Asn94 is harmless; losing it by
## substituting Lys96 costs 42 points. It is the lysine, not the missing glycan
## -- which is what panel A's residue-96 column independently implies. In the
## JU2466 background removing Asn94 gains more (+37) than restoring Thr96 does
## (+13), so the two positions interact rather than acting through one shared
## modification.
##
## Each genotype is a SINGLE PLATE at 50% pos-1 RNAi, a different dose from
## Figure 4B, so these percentages are not comparable with that figure's.
## SUPP_FIG_XX_sid2_allele_swaps_full.R carries the same data with the HT115
## controls and the full reasoning.
##
## SO THE CLAIM HAS A DEPTH, AND IT IS THE GENUS. 96T is ancestral across
## Caenorhabditis and 96K is derived within C. elegans. It is NOT conserved
## beyond the genus: nothing outside Caenorhabditis is alignable at all, which
## SUPP_FIG_XX_sid2_ortholog_search shows.
##
## Runs from a clone: reads only the staged tables.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggtext)
})

ST   <- "supplemental_data/structure"
SWAP <- "supplemental_data/hatching_assays/ju_allele_swaps_hatching.csv"
OUT  <- "plots"
WIN  <- 88:104
FLOOR <- 37

COL_ID   <- "#7E9BB5"   # matches C. elegans
COL_DIFF <- "#EDE7DC"   # does not
COL_FOC  <- "#9E4257"   # the sequon columns
COL_KEEP <- "#1A7F5A"
COL_LOST <- "#B03A2E"
## the four motif states of the editing series
COL_MOTIF <- c(NxT = "#2E4057", AxT = "#1A7F5A",
               NxK = "#B03A2E", AxK = "#E08A3C")

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
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
  for (g in c("Caenorhabditis", "elegans", "doughertyi", "afra", "sid-2", "pos-1"))
    txt <- gsub(g, paste0("*", g, "*"), txt, fixed = TRUE)
  gsub("\\*\\*", "", txt)
}

sv   <- read_tsv(file.path(ST, "sid2_ortholog_window_survey.tsv"),
                 show_col_types = FALSE)
cons <- read_tsv(file.path(ST, "sid2_ortholog_conservation.tsv"),
                 show_col_types = FALSE)
stopifnot(nrow(sv) == 40, nrow(cons) == 311)

## ===========================================================================
## A -- the window across every species where it can be read
## ===========================================================================
CE_WIN <- cons %>% filter(ce_pos %in% WIN) %>% arrange(ce_pos)
stopifnot(nrow(CE_WIN) == length(WIN))
ce_chr <- CE_WIN$ce_aa

readable <- sv %>% filter(!is.na(window_88_104))
msg("panel A: ", nrow(readable), " species with a readable window; ",
    sum(readable$confident), " clear the ", FLOOR, "% floor, of which ",
    sum(readable$confident & readable$sequon), " keep the sequon")

long <- readable %>%
  mutate(rowlab = paste0("*C. ", species, "*"),
         aa = strsplit(window_88_104, "")) %>%
  select(rowlab, species, source, percent_identity, block_identity,
         confident, sequon, aa) %>%
  unnest_longer(aa, indices_to = "i") %>%
  mutate(ce_pos = WIN[i], ce_aa = ce_chr[i],
         fill = if_else(aa == ce_aa, "identical", "different"))
stopifnot(nrow(long) == nrow(readable) * length(WIN))

ce_row <- tibble(rowlab = "*C. elegans*", species = "elegans", source = "query",
                 percent_identity = 100, block_identity = NA_real_,
                 confident = TRUE, sequon = TRUE,
                 ce_pos = WIN, i = seq_along(WIN),
                 ce_aa = ce_chr, aa = ce_chr, fill = "identical")

BAND_OK <- paste0("readable, and at or above the ", FLOOR, "% block-identity floor")
BAND_NO <- "below the floor &mdash; not counted either way"
grid <- bind_rows(ce_row, long) %>%
  mutate(band = factor(if_else(confident, BAND_OK, BAND_NO),
                       levels = c(BAND_OK, BAND_NO)),
         rowlab = fct_reorder(rowlab, percent_identity))

state <- grid %>%
  filter(ce_pos %in% c(94, 95, 96)) %>%
  select(rowlab, band, confident, ce_pos, aa) %>%
  pivot_wider(names_from = ce_pos, values_from = aa, names_prefix = "p") %>%
  mutate(st = case_when(!confident ~ "not counted",
                        p94 == "N" & p96 %in% c("S", "T") & p95 != "P" ~ "intact",
                        TRUE ~ "broken"))
msg("  sequon markers: ",
    paste(sprintf("%s %d", names(table(state$st)), table(state$st)),
          collapse = ", "))

pA <- ggplot(grid, aes(factor(ce_pos), rowlab)) +
  geom_tile(aes(fill = fill), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = aa), size = 2.6, colour = "grey20") +
  annotate("rect", xmin = which(WIN == 94) - 0.5, xmax = which(WIN == 96) + 0.5,
           ymin = -Inf, ymax = Inf, fill = NA, colour = COL_FOC,
           linewidth = 0.7) +
  geom_point(data = state, aes(x = length(WIN) + 1.2, y = rowlab, shape = st,
                               colour = st), inherit.aes = FALSE, size = 2.1) +
  facet_grid(band ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(identical = COL_ID, different = COL_DIFF),
                    name = NULL, breaks = c("identical", "different"),
                    labels = c("matches *C. elegans*", "differs")) +
  scale_shape_manual(values = c(intact = 16, broken = 4, `not counted` = 1),
                     name = "N94-x-[ST]",
                     breaks = c("intact", "broken", "not counted")) +
  scale_colour_manual(values = c(intact = COL_KEEP, broken = COL_LOST,
                                 `not counted` = "grey55"), guide = "none") +
  scale_x_discrete(expand = expansion(add = c(0.5, 2))) +
  labs(x = "*C. elegans* SID-2 residue", y = NULL,
       title = panel_title("A", "**Residue 96 is Ser or Thr across *Caenorhabditis*, and never Lys**"),
       subtitle = ital(wrap_md(paste0(
         "Residues aligned to C. elegans 88-104 in every Caenorhabditis ",
         "species whose ortholog has an HSP spanning the window; rows ordered ",
         "by identity to C. elegans. The boxed columns are the N94-x-T96 ",
         "sequon. Among the 14 species at or above the block-identity floor, ",
         "residue 96 is Thr in 11, Ser in 1, and Lys in none, and the sequon ",
         "is intact in 11. Only C. doughertyi changes residue 96 itself, to ",
         "Ala; C. afra keeps Thr96 and loses Asn94, so it is a natural AxT, ",
         "the construct in panel B. Two species below the floor appear to ",
         "carry Lys96, but at 32% and 21% block identity those calls are not ",
         "trustworthy."), 108))) +
  theme_pub(10) +
  theme(axis.text.y = element_markdown(size = 7.2),
        axis.text.x = element_text(size = 6.8),
        axis.title.x = element_markdown(size = 9),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_markdown(size = 6.8, angle = 0, hjust = 0),
        strip.background = element_rect(fill = "grey96", colour = NA),
        panel.spacing.y = grid::unit(3, "pt"),
        legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_markdown(size = 7),
        legend.title = element_text(size = 7.4),
        plot.subtitle = element_markdown(size = 7.2, colour = "grey30",
                                         lineheight = 1.3))

## ===========================================================================
## B -- the experiment that tests whether the conserved sequon matters
## ===========================================================================
msg("panel B: the N94A editing series")
stopifnot(file.exists(SWAP))
z <- qnorm(0.975)
## Wilson interval, the same expression SUPP_FIG_XX_sid2_allele_swaps_full.R
## uses, so the two figures cannot disagree about an error bar
sw <- read_csv(SWAP, show_col_types = FALSE) %>%
  filter(condition == "pos") %>%
  transmute(genotype,
            motif = sub(".*\\[", "", sub("\\]", "", `glycosylation motif`)),
            n = n_plated, hatched = n_plated - n_unhatched) %>%
  mutate(p = hatched / n,
         lo = pmax(0, (p + z^2 / (2 * n) -
                         z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)),
         hi = pmin(1, (p + z^2 / (2 * n) +
                         z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)),
         background = if_else(grepl("^JU1793", genotype),
                              "*JU1793* background", "*JU2466* background"),
         sequon = if_else(startsWith(motif, "N"), "intact", "removed"),
         motif = factor(motif, levels = c("NxT", "AxT", "NxK", "AxK")))
ORD <- c("JU1793[96T]", "JU1793[94A]", "JU1793[96K]",
         "JU2466[96T]", "JU2466[94A]", "JU2466_A[96K]", "JU2466_B[96K]")
stopifnot(setequal(sw$genotype, ORD), nrow(sw) == 7)
sw <- sw %>% mutate(genotype = factor(genotype, levels = rev(ORD)),
                    lab = paste0(motif, "  "))

cat("\n== the editing series on pos-1 (Wilson 95% CI) ==\n")
print(as.data.frame(sw %>% arrange(match(genotype, ORD)) %>%
  transmute(genotype, motif, embryos = n, hatched,
            percent = sprintf("%.1f", 100 * p),
            CI = sprintf("%.1f-%.1f", 100 * lo, 100 * hi))), row.names = FALSE)

axt <- sw %>% filter(genotype == "JU1793[94A]")
nxk <- sw %>% filter(genotype == "JU1793[96K]")
stopifnot(axt$p > 0.95, nxk$p < 0.60)
msg("  removing the sequon (AxT) leaves ", sprintf("%.1f%%", 100 * axt$p),
    " hatching; substituting Lys96 (NxK) leaves ", sprintf("%.1f%%", 100 * nxk$p),
    " -- so the glycan is not what 96K destroys")

pB <- ggplot(sw, aes(p, genotype, colour = motif)) +
  ## geom_errorbar with orientation, not geom_errorbarh: the latter is
  ## deprecated in ggplot2 4.0 and warns on every build
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                linewidth = 0.5) +
  geom_point(aes(shape = sequon), size = 2.6, fill = "white", stroke = 0.7) +
  geom_text(aes(label = lab, x = 0), hjust = 1.05, size = 2.7,
            show.legend = FALSE) +
  facet_grid(background ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = COL_MOTIF, name = NULL) +
  scale_shape_manual(values = c(intact = 16, removed = 21), guide = "none") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(-0.15, 1.02), breaks = seq(0, 1, 0.25),
                     expand = expansion(0)) +
  labs(x = "Embryos hatched on *pos-1* RNAi", y = NULL,
       title = panel_title("B", "**Removing the conserved sequon does not phenocopy 96K**"),
       subtitle = ital(wrap_md(paste0(
         "N94A removes the sequon while leaving residue 96 alone, so it should ",
         "phenocopy 96K if the glycan is what 96K destroys. It does not: AxT ",
         "has no sequon and is fully resistant, so losing the sequon by ",
         "removing Asn94 is harmless while losing it by substituting Lys96 ",
         "costs 42 points. It is the lysine, not the missing glycan -- which ",
         "is what panel A's residue-96 column independently implies. Wilson ",
         "95% intervals; one plate per genotype at 50% pos-1 RNAi, a different ",
         "dose from Figure 4B, so these are not comparable with it."), 100))) +
  theme_pub(10) +
  theme(axis.text.y = element_text(size = 7.4, colour = "grey25"),
        axis.title.x = element_markdown(size = 9),
        strip.text.y = element_markdown(size = 7, angle = 0, hjust = 0),
        strip.background = element_rect(fill = "grey96", colour = NA),
        panel.spacing.y = grid::unit(3, "pt"),
        legend.position = "bottom", legend.margin = margin(t = -2),
        legend.text = element_text(size = 7.4),
        plot.subtitle = element_markdown(size = 7.2, colour = "grey30",
                                         lineheight = 1.3))

## ===========================================================================
fig <- pA / pB + plot_layout(heights = c(1.5, 1))

ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_ortholog_conservation.pdf"), fig,
       width = 10.6, height = 10.4, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_ortholog_conservation.png"), fig,
       width = 10.6, height = 10.4, dpi = 300, bg = "white")
msg("wrote SUPP_FIG_XX_sid2_ortholog_conservation.{pdf,png}")
