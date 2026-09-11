## Supplement -- residue 96 across Caenorhabditis, against the phenotype ------
##
##   python3 scripts/make_sid2_ortholog_tables.py   # once, stages the tables
##   Rscript scripts/SUPP_FIG_XX_sid2_ortholog_conservation.R
##     -> plots/SUPP_FIG_XX_sid2_ortholog_conservation.{pdf,png}
##
##   A  the residues aligned to C. elegans 88-104 in every Caenorhabditis
##      species where that window can be read, ordered on a reference
##      phylogeny, with the N94-x-T96 sequon marked and each species'
##      published response to ingested dsRNA beside it
##   B  the editing series that asks whether the conserved sequon matters:
##      N94A removes it while leaving residue 96 alone
##
## The search behind panel A, and the calibration of what its conservation is
## worth, are in SUPP_FIG_XX_sid2_ortholog_search.
##
## THE POINT OF PUTTING THE PHENOTYPE NEXT TO THE ALIGNMENT
## Residue 96 is Ser or Thr in 12 of the 14 species whose window can be read,
## and Lys in none, so 96K is derived within C. elegans. That is a statement
## about history. It is NOT a statement about mechanism, and the sensitivity
## column is what shows the difference: of the nine species with both a
## readable window and a published call, five respond to ingested dsRNA and
## four do not, and THE SEQUON DOES NOT SEPARATE THEM. C. briggsae,
## C. remanei, C. brenneri and C. tropicalis all carry the intact motif and are
## insensitive; C. elegans, C. kamaaina, C. portoensis and C. wallacei carry it
## and respond. And C. afra, which has lost Asn94 and is a natural AxT,
## RESPONDS -- the comparative mirror of panel B, where the AxT construct is
## fully resistant.
##
## Nuez & Felix 2012 reached the same conclusion from the phenotype alone: a
## minimum of four gains or losses of environmental RNAi within the genus.
##
## THE ISOLATE IS USUALLY NOT THE ONE THAT WAS SEQUENCED, and that matters more
## here than it normally would, because the same paper reports intraspecific
## variation in C. elegans. Rows where the two differ are marked. C. elegans
## (N2), C. tropicalis (JU1373) and C. drosophilae (DF5077) are the matches.
##
## PROVISIONAL NAMES WERE MAPPED, NOT ASSUMED. Nuez & Felix scored species as
## C. sp. 6, sp. 7, sp. 11 and so on; C. sp. 11 is now C. tropicalis and
## C. sp. 10 is C. doughertyi. Every call is mapped through Felix, Braendle &
## Cutter 2014, and the mapping is deposited as sid2_species_name_map.tsv so it
## can be checked. Species the paper tested but does not call in prose are
## absent rather than guessed at.
##
## THE PHYLOGENY IS BORROWED, NOT ESTIMATED. Rows are ordered by the Open Tree
## of Life induced subtree (opentree16.1), topology only, drawn as a cladogram.
## Nothing here re-estimates a tree: a SID-2 gene tree at 42-48% identity would
## be the wrong object to draw. Species with no position in that tree -- the
## Caenorhabditis Genomes Project's unnamed sp. NN isolates, and C. oiwi --
## cannot be placed and are listed below the tree instead of being slotted in
## somewhere convenient.
##
## WHY PANEL B IS IN A CONSERVATION FIGURE
## A conserved sequon invites the reading that the glycan matters. N94A removes
## it while leaving residue 96 alone, so it should phenocopy 96K if the glycan
## is what 96K destroys. It does not: AxT has no sequon and is fully resistant,
## so losing the sequon by removing Asn94 is harmless while losing it by
## substituting Lys96 costs 42 points. It is the lysine, not the missing glycan.
## Significance is Fisher's exact test against the wild type of the same
## background, on the hatched/unhatched counts.
##
## Each genotype is a SINGLE PLATE at 50% pos-1 RNAi, a different dose from
## Figure 4B, so these percentages are not comparable with that figure's, and
## the p values describe one plate against one plate.
## SUPP_FIG_XX_sid2_allele_swaps_full.R carries the same data with the HT115
## controls and the full reasoning.
##
## Runs from a clone: reads only the staged tables.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggtext)
  library(ape)
})

ST   <- "supplemental_data/structure"
SWAP <- "supplemental_data/hatching_assays/ju_allele_swaps_hatching.csv"
OUT  <- "plots"
WIN  <- 88:104
FLOOR <- 37

## panel A geometry, all on one continuous x so the tree and the grid share it
TREE_L <- -15.0  # left edge of the cladogram
TREE_R <- -7.4    # right edge; the band from here to LAB_X holds the names
LAB_X  <- -0.5    # names are RIGHT-aligned here, just left of the grid
X_SEQ  <- length(WIN) + 1.4    # the sequon marker column
X_RNAI <- length(WIN) + 2.8    # the sensitivity column

COL_ID   <- "#7E9BB5"   # matches C. elegans
COL_DIFF <- "#EDE7DC"   # does not
COL_FOC  <- "#9E4257"   # the sequon columns
COL_KEEP <- "#1A7F5A"
COL_LOST <- "#B03A2E"
COL_RNAI <- c(sensitive = "#1A7F5A", `weakly sensitive` = "#9CC5A1",
              insensitive = "#B03A2E")
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
  for (g in c("Caenorhabditis", "elegans", "briggsae", "remanei", "brenneri",
              "tropicalis", "kamaaina", "portoensis", "wallacei", "afra",
              "doughertyi", "oiwi", "sid-2", "pos-1"))
    txt <- gsub(g, paste0("*", g, "*"), txt, fixed = TRUE)
  gsub("\\*\\*", "", txt)
}
pfmt <- function(p) if_else(p < 1e-4, sprintf("%.0e", p), sprintf("%.3f", p))

sv   <- read_tsv(file.path(ST, "sid2_ortholog_window_survey.tsv"), show_col_types = FALSE)
cons <- read_tsv(file.path(ST, "sid2_ortholog_conservation.tsv"), show_col_types = FALSE)
rnai <- read_tsv(file.path(ST, "sid2_env_rnai_sensitivity.tsv"), show_col_types = FALSE)
stopifnot(nrow(sv) == 40, nrow(cons) == 311, nrow(rnai) == 12)

## ===========================================================================
## A -- the window, ordered on the reference topology, against the phenotype
## ===========================================================================
CE_WIN <- cons %>% filter(ce_pos %in% WIN) %>% arrange(ce_pos)
ce_chr <- CE_WIN$ce_aa
stopifnot(length(ce_chr) == length(WIN), ce_chr[WIN == 96] == "T")

readable <- sv %>% filter(!is.na(window_88_104)) %>%
  bind_rows(tibble(source = "query", species = "elegans",
                   percent_identity = 100, block_identity = NA_real_,
                   confident = TRUE, sequon = TRUE,
                   window_88_104 = paste(ce_chr, collapse = "")))
msg("panel A: ", nrow(readable), " rows (including C. elegans); ",
    sum(readable$confident), " at or above the ", FLOOR, "% floor")

## ---- the reference topology decides the row order -------------------------
tr <- read.tree(file.path(ST, "sid2_species_tree.nwk"))
tr$tip.label <- sub("^Caenorhabditis_", "", sub("_ott[0-9]+$", "", tr$tip.label))
placed <- intersect(tr$tip.label, readable$species)
tr <- ladderize(keep.tip(tr, placed))
n_tip <- Ntip(tr)
## ape's own plotting coordinates: reliable, and no branch lengths are invented
invisible(plot.phylo(tr, plot = FALSE))
env <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
xx <- env$xx; yy <- env$yy
## tips get their tree y; species with no position sit below, after a gap
tip_y <- setNames(yy[seq_len(n_tip)], tr$tip.label)
unplaced <- setdiff(readable$species, placed)
GAP <- 1.6
un_y <- setNames(seq(min(tip_y) - GAP, by = -1, length.out = length(unplaced)),
                 unplaced[order(match(unplaced, readable$species[order(-readable$percent_identity)]))])
row_y <- c(tip_y, un_y)
msg("  ", length(placed), " species placed on the topology, ",
    length(unplaced), " without a position: ",
    paste(sort(unplaced), collapse = ", "))

## rescale the tree's x into the negative margin
xr <- range(xx)
tx <- TREE_L + (xx - xr[1]) / diff(xr) * (TREE_R - TREE_L)
edges <- tibble(parent = tr$edge[, 1], child = tr$edge[, 2]) %>%
  mutate(x0 = tx[parent], x1 = tx[child], y0 = yy[parent], y1 = yy[child])

## ---- the grid -------------------------------------------------------------
long <- readable %>%
  mutate(aa = strsplit(window_88_104, "")) %>%
  select(species, source, percent_identity, block_identity, confident, aa) %>%
  unnest_longer(aa, indices_to = "i") %>%
  mutate(ce_pos = WIN[i], ce_aa = ce_chr[i],
         fill = if_else(aa == ce_aa, "identical", "different"),
         y = row_y[species])
stopifnot(!anyNA(long$y))

state <- long %>% filter(ce_pos %in% c(94, 95, 96)) %>%
  select(species, y, confident, ce_pos, aa) %>%
  pivot_wider(names_from = ce_pos, values_from = aa, names_prefix = "p") %>%
  mutate(st = case_when(!confident ~ "not counted",
                        p94 == "N" & p96 %in% c("S", "T") & p95 != "P" ~ "intact",
                        TRUE ~ "broken"))

## ---- the phenotype column -------------------------------------------------
sens <- readable %>% select(species, confident) %>%
  left_join(rnai %>% select(species, response, tested_strain, same_isolate),
            by = "species") %>%
  filter(!is.na(response)) %>%
  mutate(y = row_y[species],
         mismatch = !is.na(same_isolate) & !same_isolate)
msg("  phenotype column: ", nrow(sens), " species with a published call (",
    sum(sens$response %in% c("sensitive", "weakly sensitive")), " respond, ",
    sum(sens$response == "insensitive"), " do not); ",
    sum(sens$mismatch), " tested on a different isolate from the sequenced one")
ov <- sens %>% filter(confident)
msg("  of those, ", nrow(ov), " are confidently aligned: ",
    sum(ov$response %in% c("sensitive", "weakly sensitive")), " respond, ",
    sum(ov$response == "insensitive"), " do not")
stopifnot(nrow(ov) >= 8,
          sum(ov$response %in% c("sensitive", "weakly sensitive")) >= 3,
          sum(ov$response == "insensitive") >= 3)

## row labels: species name, dimmed when the window is below the floor
lab <- readable %>%
  transmute(species, y = row_y[species], confident,
            txt = if_else(confident,
                          paste0("*C. ", species, "*"),
                          paste0("<span style='color:grey60'>*C. ", species,
                                 "*</span>")))

pA <- ggplot() +
  ## the cladogram
  geom_segment(data = edges, aes(x = x0, xend = x0, y = y0, yend = y1),
               linewidth = 0.35, colour = "grey45") +
  geom_segment(data = edges, aes(x = x0, xend = x1, y = y1, yend = y1),
               linewidth = 0.35, colour = "grey45") +
  ## the residue grid
  geom_tile(data = long, aes(i, y, fill = fill), colour = "white",
            linewidth = 0.7, width = 0.95, height = 0.85) +
  geom_text(data = long, aes(i, y, label = aa), size = 2.4, colour = "grey20") +
  annotate("rect", xmin = which(WIN == 94) - 0.5, xmax = which(WIN == 96) + 0.5,
           ymin = min(row_y) - 0.6, ymax = max(row_y) + 0.6, fill = NA,
           colour = COL_FOC, linewidth = 0.7) +
  ## the sequon marker and the phenotype
  geom_point(data = state, aes(X_SEQ, y, shape = st, colour = st), size = 1.9) +
  geom_tile(data = sens, aes(X_RNAI, y, fill = response), colour = "white",
            linewidth = 0.7, width = 0.95, height = 0.85) +
  geom_point(data = sens %>% filter(mismatch), aes(X_RNAI, y),
             shape = 8, size = 1.0, colour = "white", stroke = 0.8) +
  ## row labels, in the gap between tree and grid
  geom_richtext(data = lab, aes(LAB_X, y, label = txt), hjust = 1,
                size = 2.35, fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  scale_fill_manual(values = c(identical = COL_ID, different = COL_DIFF,
                               COL_RNAI), name = NULL,
                    breaks = c("identical", "different", "sensitive",
                               "weakly sensitive", "insensitive"),
                    labels = c("matches *C. elegans*", "differs",
                               "responds to ingested dsRNA",
                               "weakly", "does not respond")) +
  scale_shape_manual(values = c(intact = 16, broken = 4, `not counted` = 1),
                     name = "N94-x-[ST]",
                     breaks = c("intact", "broken", "not counted"),
                     labels = c("intact", "broken", "below floor")) +
  scale_colour_manual(values = c(intact = COL_KEEP, broken = COL_LOST,
                                 `not counted` = "grey55"), guide = "none") +
  scale_x_continuous(breaks = seq_along(WIN), labels = WIN,
                     limits = c(TREE_L, X_RNAI + 0.8), expand = expansion(0)) +
  scale_y_continuous(expand = expansion(add = 0.9)) +
  guides(fill = guide_legend(order = 1, nrow = 2, byrow = TRUE),
         shape = guide_legend(order = 2, nrow = 2)) +
  labs(x = "*C. elegans* SID-2 residue", y = NULL,
       title = panel_title("A", "**Residue 96 is conserved, and does not predict the phenotype**"),
       subtitle = ital(wrap_md(paste0(
         "Residues aligned to C. elegans 88-104, ordered on the Open Tree of ",
         "Life topology (cladogram, left; no branch lengths). The boxed ",
         "columns are the N94-x-T96 sequon; the rightmost column is the ",
         "published response to ingested dsRNA. Residue 96 is Thr in 11 of ",
         "the 14 readable species, Ser in 1 and Lys in none -- yet of the 9 ",
         "with a published call, 5 respond and 4 do not, and the sequon does ",
         "not separate them: C. briggsae, C. remanei, C. brenneri and ",
         "C. tropicalis carry the intact motif and are insensitive. C. afra, ",
         "which has lost Asn94 and is a natural AxT, responds. Dimmed names ",
         "are below the block-identity floor and counted neither way; an ",
         "asterisk marks a phenotype scored on a different isolate from the ",
         "sequenced one; species with no position in the reference tree are ",
         "listed below it."), 116))) +
  theme_pub(10) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 6.4),
        axis.title.x = element_markdown(size = 9),
        legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_markdown(size = 6.8),
        legend.title = element_text(size = 7.2),
        legend.margin = margin(t = -2),
        plot.subtitle = element_markdown(size = 7, colour = "grey30",
                                         lineheight = 1.3))

## ===========================================================================
## B -- the experiment, with significance against the matching wild type
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
            n = n_plated, hatched = n_plated - n_unhatched,
            unhatched = n_unhatched) %>%
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

## Fisher's exact test against the wild type of the SAME background. The
## reference is the unedited strain: JU1793[96T] and JU2466_A[96K].
REF <- c(`*JU1793* background` = "JU1793[96T]",
         `*JU2466* background` = "JU2466_A[96K]")
refs <- sw %>% filter(genotype %in% REF) %>%
  transmute(background, ref_hatched = hatched, ref_unhatched = unhatched,
            ref_p = p, ref_genotype = genotype)
stopifnot(nrow(refs) == 2)
sw <- sw %>%
  left_join(refs, by = "background") %>%
  mutate(is_ref = genotype == ref_genotype,
         p_fisher = pmap_dbl(list(hatched, unhatched, ref_hatched,
                                  ref_unhatched, is_ref),
                             function(h, u, rh, ru, r)
                               if (r) NA_real_ else
                                 fisher.test(matrix(c(h, u, rh, ru), 2))$p.value),
         ## an arrow, so "p = 0.004" beside AxT cannot be read as impairment:
         ## AxT hatches HIGHER than its wild type, not lower
         dir = case_when(is_ref ~ "", p > ref_p ~ "\u2191", TRUE ~ "\u2193"),
         genotype = factor(genotype, levels = rev(ORD)),
         lab = paste0(motif, "  "),
         sig = if_else(is_ref, "reference",
                       paste0(dir, " *p* = ", pfmt(p_fisher))))

cat("\n== the editing series on pos-1 (Wilson 95% CI, Fisher vs the same background's wild type) ==\n")
print(as.data.frame(sw %>% arrange(match(genotype, ORD)) %>%
  transmute(genotype, motif, embryos = n, hatched,
            percent = sprintf("%.1f", 100 * p),
            CI = sprintf("%.1f-%.1f", 100 * lo, 100 * hi),
            `vs wild type` = if_else(is.na(p_fisher), "(reference)",
                                     sprintf("%.3g", p_fisher)))),
  row.names = FALSE)

axt <- sw %>% filter(genotype == "JU1793[94A]")
nxk <- sw %>% filter(genotype == "JU1793[96K]")
stopifnot(axt$p > 0.95, nxk$p < 0.60, nxk$p_fisher < 1e-10)
msg("  AxT ", sprintf("%.1f%%", 100 * axt$p), " (p = ",
    signif(axt$p_fisher, 3), " vs wild type) against NxK ",
    sprintf("%.1f%%", 100 * nxk$p), " (p = ", signif(nxk$p_fisher, 3),
    ") -- the glycan is not what 96K destroys")

pB <- ggplot(sw, aes(p, genotype, colour = motif)) +
  ## geom_errorbar with orientation, not geom_errorbarh: the latter is
  ## deprecated in ggplot2 4.0 and warns on every build
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                linewidth = 0.5) +
  geom_point(aes(shape = sequon), size = 2.6, fill = "white", stroke = 0.7) +
  geom_text(aes(label = lab, x = 0), hjust = 1.05, size = 2.7,
            show.legend = FALSE) +
  geom_richtext(aes(x = 1.06, label = sig), hjust = 0, size = 2.4,
                colour = "grey30", fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt"),
                show.legend = FALSE) +
  facet_grid(background ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = COL_MOTIF, name = NULL) +
  scale_shape_manual(values = c(intact = 16, removed = 21), guide = "none") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(-0.15, 1.42), breaks = seq(0, 1, 0.25),
                     expand = expansion(0)) +
  labs(x = "Embryos hatched on *pos-1* RNAi", y = NULL,
       title = panel_title("B", "**Removing the conserved sequon does not phenocopy 96K**"),
       subtitle = ital(wrap_md(paste0(
         "N94A removes the sequon while leaving residue 96 alone, so it should ",
         "phenocopy 96K if the glycan is what 96K destroys. It does not: AxT ",
         "has no sequon and is fully resistant, while NxK in the same ",
         "background costs 42 points. It is the lysine, not the missing glycan ",
         "-- which is what panel A's residue-96 column independently implies. ",
         "Wilson 95% intervals; p is Fisher's exact test against the wild type ",
         "of the same background, with an arrow for the direction: AxT hatches ",
         "HIGHER than wild type, not lower. One plate per genotype at 50% ",
         "pos-1 RNAi, a ",
         "different dose from Figure 4B, so these are not comparable with it ",
         "and each p describes one plate against one plate."), 104))) +
  theme_pub(10) +
  theme(axis.text.y = element_text(size = 7.4, colour = "grey25"),
        axis.title.x = element_markdown(size = 9),
        strip.text.y = element_markdown(size = 7, angle = 0, hjust = 0),
        strip.background = element_rect(fill = "grey96", colour = NA),
        panel.spacing.y = grid::unit(3, "pt"),
        legend.position = "bottom", legend.margin = margin(t = -2),
        legend.text = element_text(size = 7.4),
        plot.subtitle = element_markdown(size = 7, colour = "grey30",
                                         lineheight = 1.3))

## ===========================================================================
fig <- pA / pB + plot_layout(heights = c(1.7, 1))

ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_ortholog_conservation.pdf"), fig,
       width = 10.2, height = 10.8, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_ortholog_conservation.png"), fig,
       width = 10.2, height = 10.8, dpi = 300, bg = "white")
msg("wrote SUPP_FIG_XX_sid2_ortholog_conservation.{pdf,png}")
