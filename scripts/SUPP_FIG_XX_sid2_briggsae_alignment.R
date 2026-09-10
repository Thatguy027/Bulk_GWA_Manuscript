## Supplement -- SID-2 across two species, with elegans variation on top ------
##
##   Rscript scripts/SUPP_FIG_XX_sid2_briggsae_alignment.R
##     -> plots/SUPP_FIG_XX_sid2_briggsae_alignment.{pdf,png}
##
##   A  per-residue conservation between C. elegans and C. briggsae SID-2,
##      with the topology beneath and the annotated residues marked
##   B  the 17 missense variants segregating in the elegans population, by
##      allele frequency, with the XZ1516 and JU2466 states called out
##
## WHY THE COMPARISON IS THE RIGHT ONE. C. briggsae is insensitive to
## environmental RNAi and a C. elegans sid-2 transgene confers sensitivity on it
## (Winston et al. 2007), so the two species bracket a functional difference
## that sid-2 is sufficient to explain. Residues conserved between them are
## candidates for what SID-2 needs; residues that differ are candidates for what
## makes the elegans protein work where the briggsae one does not.
##
## THE HISTIDINES ARE LABELLED "IMPLICATED", NOT "CRITICAL". H32, H168 and H175
## are the residues McEwan et al. 2012 mutated, and their triple His->Arg mutant
## internalised MORE dsRNA than wild type, not less. All three differ in
## C. briggsae -- and two of them differ TO ARGININE, which is the substitution
## that paper made. That is a genuinely interesting coincidence and it cuts
## against the simple reading rather than for it: if His->Arg increases uptake,
## the briggsae state at those positions cannot by itself explain why briggsae
## does not take dsRNA up. The panel shows the observation and the caption says
## what it does not settle.
##
## Runs from a clone: reads only the staged tables written by
## scripts/make_sid2_alignment_tables.py.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse); library(ggtext); library(patchwork)})

ST  <- "supplemental_data/structure"
OUT <- "plots"
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

aln <- read_tsv(file.path(ST, "sid2_species_alignment.tsv"), show_col_types = FALSE)
mis <- read_tsv(file.path(ST, "sid2_population_missense.tsv"), show_col_types = FALSE)
stopifnot(nrow(aln) == 311, nrow(mis) > 0)

COL_ID   <- "#B8C4BE"   # identical between species
COL_DIFF <- "#A85B18"   # differs
COL_GAP  <- "#E4E0D8"
COL_HIS  <- "#1B6C7A"   # the implicated histidines
COL_FOCAL<- "#9E4257"   # T96 / N94
COL_XZ   <- "#6A5A8C"   # XZ1516 differs from N2

## ---- annotated residues ---------------------------------------------------
ANN <- tribble(
  ~resid, ~label,   ~cls,
  32,  "H32",  "implicated His",
  168, "H168", "implicated His",
  175, "H175", "implicated His",
  34,  "D34",  "qt13 allele",
  94,  "N94",  "sequon / focal",
  96,  "T96",  "sequon / focal") %>%
  left_join(aln %>% select(resid = ce_pos, ce_aa, cb_aa, state), by = "resid") %>%
  mutate(txt = sprintf("%s&rarr;%s", ce_aa, ifelse(cb_aa == "-", "gap", cb_aa))) %>%
  arrange(resid) %>%
  ## H32/D34 and N94/T96 are two residues apart, so one label of each pair has
  ## to sit higher or they overprint
  mutate(tier = c(0, 1, 0, 1, 0, 1)[rank(resid)],
         stem.top = 0.62 + 0.20 * tier,
         lab.y    = stem.top + 0.06)

## ---- panel A: conservation and topology -----------------------------------
topo_bands <- aln %>%
  mutate(grp = cumsum(topology != lag(topology, default = first(topology)))) %>%
  group_by(grp, topology) %>%
  summarise(from = min(ce_pos), to = max(ce_pos), .groups = "drop")
TOPO_COL <- c(`signal peptide` = "#C08A2E", extracellular = "#3E8E5A",
              transmembrane = "#B5623C", cytoplasmic = "#1B6C7A")

xz <- aln %>% filter(xz1516_differs_from_n2 == "TRUE")

pA <- ggplot() +
  geom_tile(data = aln, aes(x = ce_pos, y = 0, fill = state), width = 1, height = 0.62) +
  ## the topology strip
  geom_rect(data = topo_bands, aes(xmin = from - 0.5, xmax = to + 0.5,
                                   ymin = -0.72, ymax = -0.46, fill = topology),
            colour = NA) +
  ## every band gets a label; the two narrow ones are abbreviated rather than
  ## dropped, because "which part is the membrane" is the point of the strip
  geom_richtext(data = topo_bands %>%
                  mutate(lab = ifelse(to - from > 25, topology,
                                      recode(topology, `signal peptide` = "SP",
                                             transmembrane = "TM", .default = topology))),
                aes(x = (from + to)/2, y = -0.59, label = lab),
                size = 2.2, colour = "white", fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0,4),"pt")) +
  ## XZ1516 ticks below
  geom_segment(data = xz, aes(x = ce_pos, xend = ce_pos, y = -0.95, yend = -0.78),
               colour = COL_XZ, linewidth = 0.6) +
  annotate("richtext", x = 2, y = -1.02,
           label = sprintf("<span style='color:%s'>**XZ1516 differs from N2** at %d of the 8 curated sites</span>",
                           COL_XZ, nrow(xz)),
           size = 2.3, hjust = 0, vjust = 1, fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0,4),"pt")) +
  ## annotated residues above
  geom_segment(data = ANN, aes(x = resid, xend = resid, y = 0.33, yend = stem.top),
               colour = "grey35", linewidth = 0.4) +
  geom_point(data = ANN, aes(x = resid, y = stem.top, colour = cls), size = 2.4) +
  geom_richtext(data = ANN, aes(x = resid, y = lab.y,
                                label = sprintf("**%s** %s", label, txt), colour = cls),
                size = 2.15, vjust = 0, hjust = 0.5,
                fill = alpha("white", 0.85), label.color = NA,
                label.padding = grid::unit(rep(0.6,4),"pt")) +
  scale_fill_manual(values = c(identical = COL_ID, different = COL_DIFF, gap = COL_GAP,
                               TOPO_COL), name = NULL,
                    breaks = c("identical","different","gap")) +
  scale_colour_manual(values = c(`implicated His` = COL_HIS, `qt13 allele` = "grey25",
                                 `sequon / focal` = COL_FOCAL), guide = "none") +
  scale_x_continuous(limits = c(0, 320), breaks = seq(0, 300, 50), expand = expansion(0)) +
  scale_y_continuous(limits = c(-1.15, 1.35), breaks = NULL, expand = expansion(0)) +
  labs(x = NULL, y = NULL,
       title = "**A** &nbsp;SID-2 conservation, *C. elegans* against *C. briggsae*",
       subtitle = sprintf(paste("47.3%% identity over 296 aligned positions.",
                                "All three histidines implicated in dsRNA uptake differ;",
                                "N94 and T96 do not.")))+
  theme_classic(base_size = 10.5) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_markdown(size = 11.5),
        plot.subtitle = element_markdown(size = 8, colour = "grey35"),
        plot.title.position = "plot",
        legend.position = "bottom", legend.key.size = grid::unit(8,"pt"),
        legend.text = element_text(size = 7.6))

## ---- panel B: population missense variation -------------------------------
mis2 <- mis %>%
  mutate(cons = ifelse(cb_aa == from_aa, "conserved in *C. briggsae*",
                       "differs in *C. briggsae*"),
         lab = sprintf("%s%d%s", from_aa, resid, to_aa),
         common = af >= 0.05)

pB <- ggplot(mis2, aes(resid, af)) +
  geom_rect(data = topo_bands, inherit.aes = FALSE,
            aes(xmin = from - 0.5, xmax = to + 0.5, ymin = -Inf, ymax = Inf,
                fill = topology), alpha = 0.08) +
  geom_segment(aes(xend = resid, y = 0, yend = af, colour = cons), linewidth = 0.45) +
  geom_point(aes(colour = cons), size = 2.2) +
  geom_richtext(data = mis2 %>% filter(common),
                aes(label = lab), size = 2.2, vjust = -0.55, colour = "grey15",
                fill = alpha("white", 0.75), label.color = NA,
                label.padding = grid::unit(rep(0.5,4),"pt")) +
  geom_vline(data = ANN, aes(xintercept = resid), linetype = "dotted",
             linewidth = 0.3, colour = "grey55") +
  scale_fill_manual(values = TOPO_COL, guide = "none") +
  scale_colour_manual(values = c(`conserved in *C. briggsae*` = COL_ID,
                                 `differs in *C. briggsae*` = COL_DIFF), name = NULL) +
  scale_x_continuous(limits = c(0, 320), breaks = seq(0, 300, 50), expand = expansion(0)) +
  scale_y_continuous(limits = c(0, 0.58), breaks = seq(0, 0.5, 0.1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "SID-2 residue (*C. elegans* numbering)", y = "CeNDR allele frequency",
       title = "**B** &nbsp;Missense variation in the *C. elegans* population",
       subtitle = paste("All 17 missense variants in the sid-2 span across 540 CeNDR isotypes;",
                        "labelled where AF &ge; 5%. Dotted lines mark the panel A residues.")) +
  theme_classic(base_size = 10.5) +
  theme(axis.line = element_line(linewidth = 0.3),
        plot.title = element_markdown(size = 11.5),
        plot.subtitle = element_markdown(size = 8, colour = "grey35"),
        axis.title.x = element_markdown(), plot.title.position = "plot",
        legend.position = "bottom", legend.key.size = grid::unit(8,"pt"),
        legend.text = element_markdown(size = 7.6))

fig <- pA / pB + plot_layout(heights = c(1, 1.15))
ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_briggsae_alignment.pdf"), fig,
       width = 9.2, height = 6.4, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_briggsae_alignment.png"), fig,
       width = 9.2, height = 6.4, dpi = 300, bg = "white")
msg("wrote SUPP_FIG_XX_sid2_briggsae_alignment.{pdf,png}")

cat("\n== the annotated residues across species ==\n")
print(as.data.frame(ANN %>% select(resid, label, cls, ce_aa, cb_aa, state)), row.names = FALSE)
cat("\n== population missense, by conservation in C. briggsae ==\n")
print(as.data.frame(mis2 %>% count(cons, common)), row.names = FALSE)
cat("\n== XZ1516 sites ==\n")
print(as.data.frame(xz %>% select(ce_pos, ce_aa, n2_aa, xz1516_aa, cb_aa, state)), row.names = FALSE)
