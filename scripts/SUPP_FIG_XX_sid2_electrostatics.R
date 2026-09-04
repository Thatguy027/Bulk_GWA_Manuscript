## Supplement -- SID-2 surface charge and where T96 sits --------------------
##
##   python3 scripts/sid2_ribbon_render.py            # once, builds the assets
##   Rscript scripts/SUPP_FIG_XX_sid2_electrostatics.R
##     -> plots/SUPP_FIG_XX_sid2_electrostatics.{pdf,png}
##
##   A  the ectodomain coloured by residue class, with T96 and the published
##      uptake-critical residues marked
##   B  net charge of the ectodomain against pH, with and without T96K
##   C  how close the uptake-critical residues are to T96, against the null
##
## WHAT THIS FIGURE ARGUES, AND WHAT IT DOES NOT
## It argues that the compartment SID-2 works in, not the variant, dominates
## the charge on its dsRNA-facing surface -- and that T96K nudges that surface
## in the direction the phenotype implies. It does NOT estimate a binding
## affinity, a Delta-Delta-G, or a change in dsRNA occupancy. Panel B is
## arithmetic over the residue composition; panels A and C are geometry. None
## of it is a docking result.
##
## WHY THERE IS NO POTENTIAL MAP HERE
## data/structure_modeling/claude_docking/stage4_electrostatics contains one:
## a Coulomb potential with a distance-dependent dielectric, sampled on a
## plane 55 A above the protein, reported in arbitrary units, with no ionic
## screening. Its headline panel is the difference between the T96 and K96
## surfaces, which is necessarily a smooth positive blob centred on residue 96
## -- adding +1 e there cannot produce anything else -- so it carries no
## information beyond "lysine is positive". A plane slice also does not report
## what a dsRNA at the surface would feel. If a potential map is wanted later,
## it should be a screened (Debye-Huckel or Poisson-Boltzmann) potential on the
## solvent-accessible surface, in kT/e, on one fixed backbone.
##
## A CORRECTION TO THE EARLIER stage4 NUMBERS
## That analysis reported the ectodomain flipping to strongly net POSITIVE at
## pH 4.4. It reached that by assigning fixed fractional charges of -0.240 to
## Asp and -0.334 to Glu, which correspond to pKa values near 3.4-3.7. With
## the standard side-chain pKa values used here (Asp 3.90, Glu 4.07), Asp is
## still about 76% ionised at pH 4.4, and the ectodomain comes out close to
## ELECTRONEUTRAL rather than positive. The corrected numbers are in panel B
## and are the ones to quote.
##
## THE DOCKING IN THAT FOLDER IS NOT USABLE AND IS NOT USED
## The HADDOCK runs returned positive HADDOCK scores (+98 to +200, where
## favourable is negative), interaction energies of order 1e4 and
## restraint-violation energies of 480-1380, with 4-13 structures per cluster.
## The two arms were also run against different active-residue lists (91
## residues each), so they are not comparable even in principle.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggtext)
  library(png)
})

OUT    <- "plots"
PERRES <- "supplemental_data/structure/sid2_per_residue.tsv"
ASSET  <- "plots/assets/sid2_ecd_ribbon_charge.png"

ECD <- c(21, 193)          # DeepTMHMM extracellular span
PH_NEUTRAL <- 7.4          # body/cytosol reference
PH_LUMEN   <- 4.4          # C. elegans intestinal lumen, the pH used by the
                           # earlier stage4 analysis

COL_JU1793 <- "#F34C00"
COL_FUNC   <- "#16324A"
COL_WT     <- "#2E4057"
COL_MUT    <- "#F34C00"

CLASS_COL <- c(`Basic (Arg, Lys)` = "#3B6FB6", `His` = "#8E6BAF",
               `Acidic (Asp, Glu)` = "#C0392B", `Other` = "#D5DCE1")

## published alleles with a reported dsRNA-uptake effect (UniProt G5EEV9, via
## data/structure_modeling/claude_docking/stage3b_dimer/dimerization_report.txt
## -- VERIFY these residue numbers against UniProt before publication)
FUNC <- tribble(
  ~pos, ~effect,
  32,   "reduced uptake",
  34,   "qt13, complete RNAi resistance",
  168,  "no detectable uptake",
  175,  "reduced uptake")

## side-chain pKa values, Nozaki & Tanford / standard set
PKA <- tribble(~res, ~pka, ~sign,
               "ASP", 3.90, -1, "GLU", 4.07, -1, "HIS",  6.04, +1,
               "CYS", 8.18, -1, "TYR", 10.46, -1, "LYS", 10.54, +1,
               "ARG", 12.48, +1)

panel_title <- function(letter, txt = NULL) {
  lt <- paste0("<span style='font-size:13pt;color:#111111'>**", letter,
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
          legend.key.size = grid::unit(8, "pt"))
}
wrap_md <- function(txt, width = 72)
  paste(strwrap(txt, width = width), collapse = "<br>")
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

pr <- read_tsv(PERRES, show_col_types = FALSE)
ecd <- pr %>% filter(resid >= ECD[1], resid <= ECD[2])
msg("ectodomain ", ECD[1], "-", ECD[2], ": ", nrow(ecd), " residues")

## ===========================================================================
## A -- residue class on the structure
## ===========================================================================
msg("panel A: residue-class ribbon")
stopifnot(file.exists(ASSET))
img <- readPNG(ASSET)

key <- tibble(lab = factor(names(CLASS_COL), levels = names(CLASS_COL)))
p_key <- ggplot(key, aes(1, lab, fill = lab)) +
  geom_tile(width = 0.4, height = 0.62) +
  scale_fill_manual(values = CLASS_COL, guide = "none") +
  scale_y_discrete(limits = rev(levels(key$lab)), position = "right") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 8) +
  theme(axis.text.y.right = element_text(hjust = 0, size = 7.2,
                                         colour = "grey20",
                                         margin = margin(l = 3)),
        plot.margin = margin(1, 1, 1, 1))

pA <- (ggplot() +
  annotation_raster(img, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                    interpolate = TRUE) +
  coord_fixed(ratio = dim(img)[1] / dim(img)[2], xlim = c(0, 1),
              ylim = c(0, 1), expand = FALSE) +
  labs(title = panel_title("A", "**Charged residues on the ectodomain**"),
       subtitle = wrap_md(paste0(
         "AlphaFold3 model, residues ", ECD[1], "&ndash;", ECD[2],
         " of chain A, coloured by residue class. Histidine is separated ",
         "because it is the class that titrates between neutral pH and the ",
         "acidic gut lumen. T96 in orange; residues with a published uptake ",
         "defect in white."), 78)) +
  theme_void(base_size = 11) +
  theme(plot.title = element_markdown(size = 11.5),
        plot.subtitle = element_markdown(size = 8.2, colour = "grey30"),
        plot.title.position = "plot")) +
  inset_element(p_key, left = 0, bottom = 0.60, right = 0.26, top = 1)

## ===========================================================================
## B -- net charge against pH
## ===========================================================================
msg("panel B: net charge vs pH")
counts <- ecd %>% count(resname) %>% deframe()
get <- function(r) unname(ifelse(is.na(counts[r]), 0, counts[r]))

## Henderson-Hasselbalch over the side chains. The termini are omitted: the
## ectodomain is not a free polypeptide -- it runs into the TM helix -- so a
## free alpha-carboxylate would be an artefact of where we cut it.
net_charge <- function(ph, extra_lys = 0) {
  q <- 0
  for (i in seq_len(nrow(PKA))) {
    r <- PKA$res[i]; k <- PKA$pka[i]; sgn <- PKA$sign[i]
    n <- get(r) + if (r == "LYS") extra_lys else 0
    q <- q + if (sgn > 0) n / (1 + 10^(ph - k)) else -n / (1 + 10^(k - ph))
  }
  q
}

cat("\n== ectodomain titratable composition ==\n")
print(as.data.frame(PKA %>% mutate(n = map_dbl(res, get)) %>%
                      select(res, pka, sign, n)), row.names = FALSE)

curve <- tibble(ph = seq(2, 12, by = 0.02)) %>%
  mutate(`96T (JU1793)` = map_dbl(ph, net_charge),
         `96K (JU2466)` = map_dbl(ph, ~net_charge(.x, extra_lys = 1))) %>%
  pivot_longer(-ph, names_to = "allele", values_to = "q") %>%
  mutate(allele = factor(allele, levels = c("96T (JU1793)", "96K (JU2466)")))

pI <- function(extra) {
  uniroot(function(x) net_charge(x, extra), c(2, 12))$root
}
q74 <- net_charge(PH_NEUTRAL); q44 <- net_charge(PH_LUMEN)
msg("  net charge 96T: pH ", PH_NEUTRAL, " ", sprintf("%+.1f", q74),
    " e | pH ", PH_LUMEN, " ", sprintf("%+.1f", q44), " e | swing ",
    sprintf("%+.1f", q44 - q74), " e")
msg("  isoelectric point: 96T ", sprintf("%.2f", pI(0)), " | 96K ",
    sprintf("%.2f", pI(1)))
msg("  T96K adds +1 e = ", sprintf("%.0f%%", 100 / (q44 - q74)),
    " of the pH-driven swing")

q44_mut <- net_charge(PH_LUMEN, extra_lys = 1)
msg("  at pH ", PH_LUMEN, ": 96T ", sprintf("%+.1f", q44), " e -> 96K ",
    sprintf("%+.1f", q44_mut), " e (sign change: ",
    ifelse(sign(q44) != sign(q44_mut), "yes", "no"), ")")

marks <- tibble(ph = c(PH_NEUTRAL, PH_LUMEN),
                lab = c("pH 7.4", "pH 4.4<br>gut lumen"))
## the isoelectric points, which land either side of the lumenal pH
pis <- tibble(allele = factor(c("96T (JU1793)", "96K (JU2466)"),
                              levels = levels(curve$allele)),
              ph = c(pI(0), pI(1)), q = 0)

pB <- ggplot(curve, aes(ph, q, colour = allele)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey55") +
  geom_vline(data = marks, aes(xintercept = ph), inherit.aes = FALSE,
             linetype = "dotted", linewidth = 0.4, colour = "grey45") +
  geom_line(linewidth = 0.8) +
  geom_point(data = pis, aes(ph, q, colour = allele), size = 2.1,
             show.legend = FALSE) +
  geom_richtext(data = marks, aes(x = ph, y = Inf, label = lab),
                inherit.aes = FALSE, hjust = -0.06, vjust = 1.1, size = 2.7,
                colour = "grey30", fill = "white", label.color = NA,
                label.padding = grid::unit(c(1, 2, 1, 2), "pt")) +
  scale_colour_manual(values = c(`96T (JU1793)` = COL_WT,
                                 `96K (JU2466)` = COL_MUT), name = NULL) +
  scale_x_continuous(breaks = seq(2, 12, 2)) +
  labs(x = "pH", y = "Net charge of the ectodomain (e)",
       title = panel_title("B", "**At gut pH the variant changes the sign of the surface**"),
       subtitle = wrap_md(sprintf(paste0(
         "Henderson&ndash;Hasselbalch over the side chains of residues %d",
         "&ndash;%d; termini omitted. The ectodomain carries %+.1f e at pH ",
         "7.4 but only %+.1f e at the gut-lumen pH of 4.4 &mdash; its ",
         "isoelectric point, %.2f, is essentially the pH it works at (dots). ",
         "Against that near-neutral background the +1 e from T96K takes the ",
         "surface from %+.1f e to %+.1f e and the isoelectric point to %.2f. ",
         "The variant is a small change in absolute terms, but at lumenal pH ",
         "it is the difference between a slightly negative and a slightly ",
         "positive face &mdash; and 96K is the allele with efficient uptake."),
         ECD[1], ECD[2], q74, q44, pI(0), q44, q44_mut, pI(1)), 78)) +
  theme_pub() +
  theme(legend.position = c(0.02, 0.14), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "white", colour = NA))

## ===========================================================================
## C -- distance of the uptake-critical residues from T96
## ===========================================================================
msg("panel C: distance to T96")
xyz <- ecd %>% select(resid, x, y, z) %>% column_to_rownames("resid") %>%
  as.matrix()
p96 <- xyz[as.character(96), ]
d <- tibble(resid = as.integer(rownames(xyz)),
            dist = sqrt(rowSums((xyz - matrix(p96, nrow(xyz), 3,
                                              byrow = TRUE))^2))) %>%
  filter(resid != 96)

fd <- FUNC %>% left_join(d, by = c("pos" = "resid")) %>%
  mutate(lab = paste0(ecd$resname[match(pos, ecd$resid)], pos))
med <- median(d$dist)

## two null comparisons, both reported because both are unimpressive
set.seed(42)
obs_mean <- mean(fd$dist)
perm <- replicate(20000, mean(sample(d$dist, nrow(fd))))
p_perm <- mean(perm <= obs_mean)
frac20 <- mean(d$dist <= 20)
k20 <- sum(fd$dist <= 20)
p_bin <- sum(dbinom(k20:nrow(fd), nrow(fd), frac20))
msg("  observed mean ", round(obs_mean, 1), " A vs median ", round(med, 1),
    " A | permutation p = ", sprintf("%.3f", p_perm))
msg("  ", k20, " of ", nrow(fd), " within 20 A; ",
    sprintf("%.0f%%", 100 * frac20), " of the ectodomain is | binomial p = ",
    sprintf("%.3f", p_bin))

cat("\n== distance from T96 (Ca-Ca) ==\n")
print(as.data.frame(fd %>% transmute(residue = lab,
        `distance (A)` = round(dist, 1), effect)), row.names = FALSE)

pC <- ggplot(d, aes(dist)) +
  geom_histogram(binwidth = 2, fill = "grey85", colour = "white",
                 linewidth = 0.2) +
  geom_vline(xintercept = med, linetype = "dashed", linewidth = 0.4,
             colour = "grey35") +
  geom_vline(data = fd, aes(xintercept = dist), linewidth = 0.6,
             colour = COL_FUNC) +
  geom_richtext(data = fd %>% mutate(y = Inf),
                aes(x = dist, y = y, label = lab), inherit.aes = FALSE,
                angle = 90, hjust = 1.05, vjust = -0.2, size = 2.6,
                colour = COL_FUNC, fill = "white", label.color = NA,
                label.padding = grid::unit(c(1, 1, 1, 1), "pt")) +
  annotate("richtext", x = med, y = 0, label = "median", hjust = -0.1,
           vjust = -0.6, angle = 90, size = 2.6, colour = "grey35",
           fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt")) +
  scale_x_continuous(breaks = seq(0, 70, 10)) +
  labs(x = "C&alpha; distance from residue 96 (&Aring;)",
       y = "Ectodomain residues",
       title = panel_title("C", "**Close, but not more than chance would give**"),
       subtitle = wrap_md(sprintf(paste0(
         "Distances from residue 96 to every other ectodomain residue, with ",
         "the four published uptake-critical residues marked. Three of the ",
         "four are nearer than the median (%.1f &Aring;), but %.0f%% of the ",
         "ectodomain lies within 20 &Aring;, so %d of 4 landing there is not ",
         "surprising: binomial *p* = %.2f, and a permutation test on their ",
         "mean distance gives *p* = %.2f. Read this as spatial context for ",
         "T96, not as evidence of a shared site."),
         med, 100 * frac20, k20, p_bin, p_perm), 78)) +
  theme_pub() +
  theme(axis.title.x = element_markdown())

## ===========================================================================
fig <- pA / (pB | pC) + plot_layout(heights = c(0.85, 1))

ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_electrostatics.pdf"), fig,
       width = 11.2, height = 8.4, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_electrostatics.png"), fig,
       width = 11.2, height = 8.4, dpi = 300, bg = "white")
msg("wrote SUPP_FIG_XX_sid2_electrostatics.{pdf,png}")
