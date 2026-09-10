## Figure 4 -- sid-2, the candidate gene -------------------------------------
##
##   python3 scripts/sid2_ribbon_render.py    # once, secondary structure
##   python3 scripts/sid2_zoom_render.py      # once, the two panel-A images
##   Rscript scripts/Figure4_sid2.R           -> plots/Figure4_sid2.{pdf,png}
##
##   A  the residue-96 allele swaps in JU1793 and JU2466, at full-strength
##      pos-1 RNAi
##   B  the same swap in N2, at 25% pos-1 RNAi -- a third background, with two
##      independently edited lines
##   C  the SID-2 ectodomain with the lumenal face up, and a zoom on T96
##      against the residues already known to be required for dsRNA uptake
##   D  sid-2 coding variation in the wild population: the protein drawn
##      vertically with its topology, a local-net-charge strip on panel C's
##      scale to its left, and every annotated protein-altering variant listed
##      with its CeNDR allele frequency
##
## RNAI DOSE DIFFERS BETWEEN THE TWO EXPERIMENTS
## Panel A is 50% pos-1 RNAi bacteria and panel B is 25%. Neither dose is
## recorded in the source data files, so both come from the lab record. The
## 25% dose in B is the only dilution with dynamic range in the N2 background
## -- see scripts/n2_swap_panels.R -- and the consequence is that hatching
## percentages are NOT comparable between panels A and B. FIGURE_CAPTIONS.txt
## says so where a reader will see it.
##
## Panels carry letters only. Everything that was a panel title or subtitle is
## now in FIGURE_CAPTIONS.txt, which is where a reader of the manuscript will
## meet it.
##
## scripts/Figure4.R is untouched and still builds Figure4_clean.
##
## PANEL A ORIENTATION
## Both images come from the membrane-oriented model
## (data/structure_modeling/claude_docking/stage4_electrostatics/
## sid2_wt_elec.pdb), which the earlier stage3 work rotated onto the membrane
## normal: z = 0 is the bilayer centre and +z is the intestinal lumen. That
## frame is kept and +z is up the page, so the picture is anatomically right
## rather than an arbitrary view. The bilayer is marked on the overview and the
## dashed ring is the 20 A zoom region.
##
## In the zoom the ribbon is deliberately grey: it is scaffold. Side chains are
## drawn for T96 and for the three uptake-critical residues that share its
## face, and the dashed lines are Ca-Ca distances -- the same numbers quoted in
## the text. Nearest heavy-atom distances are shorter (D34 9.9, H32 14.6, H168
## 16.9 A) and are printed by the render script for the caption.
##
## Secondary structure is assigned from the backbone with the Kabsch & Sander
## hydrogen-bond criterion; there is no DSSP binary here. It is corroborated
## independently: the criterion calls residues 194-211 100% helix, exactly the
## span DeepTMHMM calls the transmembrane helix.
##
## WHAT IS NOT CLAIMED
##   the dimer          ipTM 0.45-0.46 across all five models, so the
##                      interface is not credible and is not drawn
##   any docking        the HADDOCK runs in data/structure_modeling returned
##                      positive scores with violated restraints and
##                      mismatched restraint sets between arms; unusable
##   glycosylation      T96 is the threonine of an N94-C95-T96 sequon, but the
##                      N94A test of that failed -- see
##                      SUPP_FIG_XX_sid2_allele_swaps_full.R
##   a panel-wide effect  the variant has none: Wilcoxon p = 0.99 across the
##                      231-strain panel -- see
##                      SUPP_FIG_XX_sid2_allele_in_panel.R
##   proximity as evidence  three of four uptake-critical residues are nearer
##                      than the median, but 41% of the ectodomain is within
##                      20 A (binomial p = 0.19). Spatial context only.
##
## Published alleles are from the UniProt annotation for G5EEV9 recorded in
## data/structure_modeling/claude_docking/stage3b_dimer/dimerization_report.txt
## -- VERIFY the residue numbers against UniProt before publication.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggtext)
  library(png)
})

## the N2 swap: data prep, colours and the dose reasoning
source("scripts/n2_swap_panels.R")

OUT  <- "plots"
SWAP <- "supplemental_data/hatching_assays/ju_allele_swaps_hatching.csv"
VAR  <- "supplemental_data/structure/sid2_variants_cendr.tsv"
PERRES2 <- "supplemental_data/structure/sid2_per_residue.tsv"
LOCALQ  <- "supplemental_data/structure/sid2_local_charge.tsv"
## Panel C is coloured by LOCAL NET CHARGE, not by secondary structure.
##
## The claim the panel makes is a charge claim: in a pathway where dsRNA
## recognition is known to be electrostatic and sequence-independent (McEwan et
## al. 2012, whose triple His->Arg mutant internalised MORE dsRNA than wild
## type), T96 sits in the most positive solvent-exposed pocket of an otherwise
## acidic lumenal domain, and T96K adds another positive charge there. Colouring
## by beta-strand and coil showed the fold instead, which is not the argument
## and is not in dispute.
##
## Both PNGs carry their own colourbar, so this script adds no colour key --
## the colour encodes a continuous quantity and cannot be named in a subtitle
## the way the two residue classes could.
##
## The secondary-structure renders are kept for the model-confidence supplement
## (scripts/sid2_zoom_render.py writes sid2_zoom_t96.png), which is where the
## annotated-residue classes are shown.
OVER <- "plots/assets/sid2_overview_charge.png"
ZOOM <- "plots/assets/sid2_zoom_charge.png"

## the strain colours as defined for this manuscript; do not re-map these
COL_JU1793 <- "#F34C00"
COL_JU2466 <- "#40B4AB"
COL_EDIT   <- "grey70"
SSE_COL <- c(Strand = "#4E6E8A", Coil = "#C4D0D8")
TOPO_COL2 <- c(`Signal peptide` = "grey72", `Extracellular` = "#9EC5DE",
               `TM helix` = "#37474F", `Cytoplasmic` = "grey88")

## letters only: the descriptive text lives in FIGURE_CAPTIONS.txt
panel_title <- function(letter) {
  paste0("<span style='font-size:13pt;color:#111111'>**", letter, "**</span>")
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
wrap_md <- function(txt, width = 78)
  paste(strwrap(txt, width = width), collapse = "<br>")
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## ===========================================================================
## A -- overview and zoom, as ONE panel
##
## Both images go into a single ggplot on an inch-scaled coordinate system with
## coord_fixed, rather than two plots in a patchwork. That way the panel gets
## one title and one subtitle spanning both images, and each image keeps its
## own aspect ratio -- annotation_raster stretches to whatever rectangle it is
## given, so the rectangles are computed from the pixel dimensions.
## ===========================================================================
msg("panel A: structure")
stopifnot(file.exists(OVER), file.exists(ZOOM))
im_ov <- readPNG(OVER); im_zm <- readPNG(ZOOM)
asp <- function(im) dim(im)[2] / dim(im)[1]          # width / height
msg("  overview ", dim(im_ov)[2], "x", dim(im_ov)[1],
    " (aspect ", round(asp(im_ov), 3), ") | zoom ",
    dim(im_zm)[2], "x", dim(im_zm)[1], " (aspect ", round(asp(im_zm), 3), ")")

H    <- 3.15                       # common display height, inches
GAP  <- 0.30
W_OV <- H * asp(im_ov)
W_ZM <- H * asp(im_zm)
X_ZM <- W_OV + GAP
TOTAL <- X_ZM + W_ZM
CAP_Y <- -0.20
msg("  composite: ", round(TOTAL, 2), " x ", round(H, 2), " in")

caps <- tibble(x = c(W_OV / 2, X_ZM + W_ZM / 2), y = CAP_Y,
               ## the left caption does NOT repeat the quantity: the key
               ## directly below it already names "net charge within 12 A (e),
               ## pH 4.4", and saying it twice in adjacent lines reads as an
               ## error rather than as emphasis
               lab = c("Ectodomain", "T96 pocket: K93, K132"))

## the diverging key, under the left image. QLIM must equal the renderer's.
QLIM  <- 2
KEY_W <- W_OV * 0.82
KEY_Y <- CAP_Y - 0.42
KEY_H <- 0.13
ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(64)
key <- tibble(i = seq_along(ramp),
              xmin = (i - 1) / length(ramp) * KEY_W,
              xmax = i / length(ramp) * KEY_W,
              col = ramp)

p_struct <- ggplot() +
  annotation_raster(im_ov, xmin = 0, xmax = W_OV, ymin = 0, ymax = H,
                    interpolate = TRUE) +
  annotation_raster(im_zm, xmin = X_ZM, xmax = TOTAL, ymin = 0, ymax = H,
                    interpolate = TRUE) +
  geom_richtext(data = caps, aes(x, y, label = lab), size = 2.8,
                colour = "grey30", vjust = 1, fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  ## The charge key, drawn here rather than in the render.
  ##
  ## Panel C's colour now encodes a CONTINUOUS quantity, so unlike the
  ## secondary-structure version it cannot be named in a subtitle -- it needs a
  ## scale. The renders do not carry one, so the first version of this panel
  ## shipped a diverging colour map with no key at all.
  ##
  ## The ramp is RdBu NOT reversed, matching Normalize(-QLIM, +QLIM) in
  ## scripts/sid2_charge_render.py: matplotlib's RdBu runs red at the low end to
  ## blue at the high end, so red is NEGATIVE charge and blue POSITIVE.
  ## Reversing it here would silently invert the key while leaving the structure
  ## colours untouched -- the failure would look like a result.
  geom_rect(data = key, aes(xmin = xmin, xmax = xmax, ymin = KEY_Y,
                            ymax = KEY_Y + KEY_H), fill = key$col) +
  annotate("text", x = 0, y = KEY_Y - 0.07, hjust = 0, size = 2.2,
           colour = "grey30", label = paste0("\u2212", QLIM)) +
  annotate("text", x = KEY_W, y = KEY_Y - 0.07, hjust = 1, size = 2.2,
           colour = "grey30", label = paste0("+", QLIM)) +
  annotate("text", x = 0, y = KEY_Y + KEY_H + 0.11, hjust = 0, size = 2.3,
           colour = "grey30",
           label = "net charge within 12 \u00c5 (e), pH 4.4") +
  ## a little x padding: the captions are centred under each image and the
  ## left one is flush at x = 0, so without it clip = "off" let it run off
  coord_fixed(ratio = 1, xlim = c(-0.45, TOTAL + 0.1),
              ylim = c(KEY_Y - 0.20, H), expand = FALSE, clip = "off") +
  labs(title = panel_title("C")) +
  theme_void(base_size = 11) +
  theme(plot.title = element_markdown(size = 11.5),
        plot.subtitle = element_markdown(size = 8.3, colour = "grey30"),
        plot.title.position = "plot",
        plot.margin = margin(2, 4, 2, 2))

## ===========================================================================
## B -- the residue-96 swaps, in the Figure4_clean style
## ===========================================================================
msg("panel B: residue-96 swaps")
z <- qnorm(0.975)
raw <- read_csv(SWAP, show_col_types = FALSE) %>%
  mutate(n_hatched = n_plated - n_unhatched,
         genotype = str_replace(genotype, "_A", ""))

LEV <- c("JU1793[96T]", "JU1793[96K]", "JU2466[96K]", "JU2466[96T]")
BAR <- c(`JU1793[96T]` = COL_JU1793, `JU1793[96K]` = COL_EDIT,
         `JU2466[96K]` = COL_JU2466, `JU2466[96T]` = COL_EDIT)

df_pos <- raw %>%
  filter(condition == "pos",
         strain %in% c("JU1793", "JU2466_A", "wSZ200", "wSZ206")) %>%
  mutate(p = n_hatched / n_plated,
         lo = pmax(0, (p + z^2 / (2 * n_plated) -
                z * sqrt(p * (1 - p) / n_plated + z^2 / (4 * n_plated^2))) /
                (1 + z^2 / n_plated)),
         hi = pmin(1, (p + z^2 / (2 * n_plated) +
                z * sqrt(p * (1 - p) / n_plated + z^2 / (4 * n_plated^2))) /
                (1 + z^2 / n_plated)),
         genotype = factor(genotype, levels = LEV))
stopifnot(nrow(df_pos) == 4, !any(is.na(df_pos$genotype)))

ctrl <- raw %>% filter(condition == "ht115")
msg("  HT115 controls (not plotted): ", nrow(ctrl), " constructs, ",
    sprintf("%.3f-%.3f", min(ctrl$n_hatched / ctrl$n_plated),
            max(ctrl$n_hatched / ctrl$n_plated)))

fisher_pair <- function(s1, s2) {
  d1 <- raw %>% filter(condition == "pos", strain == s1)
  d2 <- raw %>% filter(condition == "pos", strain == s2)
  fisher.test(matrix(c(d1$n_hatched, d1$n_unhatched,
                       d2$n_hatched, d2$n_unhatched), nrow = 2,
                     byrow = TRUE))$p.value
}
cmp <- tribble(~s1, ~s2, ~x1, ~x2,
               "JU1793", "wSZ200", 1, 2,
               "JU2466_A", "wSZ206", 3, 4) %>%
  rowwise() %>%
  mutate(pval = fisher_pair(s1, s2),
         label = if (pval < 0.001) "*p* < 0.001" else sprintf("*p* = %.3f", pval)) %>%
  ungroup()
msg("  Fisher: ", paste(sprintf("%s vs %s p=%.3g", cmp$s1, cmp$s2, cmp$pval),
                        collapse = " | "))

Y_BAR <- 1.04; TIP <- 0.035
segs <- cmp %>% rowwise() %>%
  reframe(x = c(x1, x1, x2), xend = c(x2, x1, x2),
          y = c(Y_BAR, Y_BAR, Y_BAR),
          yend = c(Y_BAR, Y_BAR + TIP, Y_BAR + TIP))
txt <- cmp %>% mutate(x = (x1 + x2) / 2, y = Y_BAR + TIP + 0.015)

cat("\n== residue-96 swaps on pos-1 (Wilson 95% CI) ==\n")
print(as.data.frame(df_pos %>%
  transmute(strain, genotype, embryos = n_plated, hatched = round(p, 3),
            CI = sprintf("%.3f-%.3f", lo, hi))), row.names = FALSE)

p_ju <- ggplot(df_pos, aes(genotype, p, fill = genotype)) +
  geom_col(colour = "black", width = 0.7, linewidth = 0.35) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, colour = "black",
                linewidth = 0.45) +
  geom_segment(data = segs, aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_richtext(data = txt, aes(x = x, y = y, label = label),
                inherit.aes = FALSE, size = 3, vjust = 0, fill = NA,
                label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  scale_fill_manual(values = BAR, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, 0.25), limits = c(0, 1.20),
                     expand = expansion(0)) +
  ## n above each bar's interval rather than in the tick label: the axis
  ## labels are already two lines when angled
  geom_text(aes(y = hi + 0.035, label = paste0("n = ", n_plated)),
            size = 2.7, colour = "grey25") +
  labs(x = NULL, y = "Embryos hatched",
       title = panel_title("A")) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.margin = margin(t = 4, r = 8, b = 4, l = 6))

## ===========================================================================
## C -- the same swap in N2, at 25% pos-1
##
## Shown at 25% because that is the only dose with dynamic range in this
## background: at 0% every genotype hatches and at 50% and above every
## genotype is dead. The full dose series is
## scripts/SUPP_FIG_XX_n2_swap_dose.R.
## ===========================================================================
msg("panel C: N2 swap at ", FOCAL_DOSE, "% pos-1")
n2 <- n2_swap_data()
n2_t <- n2_swap_tests(n2) %>% filter(dose == FOCAL_DOSE)
n25 <- n2 %>% filter(dose == FOCAL_DOSE)
msg("  ", nrow(n25), " genotypes | hatching ",
    paste(sprintf("%s %.3f", n25$strain, n25$p), collapse = " | "))
msg("  Fisher vs N2: ",
    paste(sprintf("%s p=%.3g", n2_t$strain, n2_t$fisher_p), collapse = " | "))

YC <- 0.44; TIPC <- 0.02
segC <- n2_t %>% mutate(x2 = match(line, levels(n2$line))) %>% rowwise() %>%
  reframe(x = c(1, 1, x2), xend = c(x2, 1, x2),
          y = c(YC, YC, YC) + (x2 - 2) * 0.07,
          yend = c(YC, YC + TIPC, YC + TIPC) + (x2 - 2) * 0.07)
lbC <- n2_t %>% mutate(x2 = match(line, levels(n2$line)),
                       x = (1 + x2) / 2,
                       y = YC + TIPC + 0.006 + (x2 - 2) * 0.07)

p_n2 <- ggplot(n25, aes(line, p, fill = line)) +
  geom_col(colour = "black", width = 0.66, linewidth = 0.35) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.22, colour = "black",
                linewidth = 0.45) +
  geom_segment(data = segC, aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_richtext(data = lbC, aes(x = x, y = y, label = fmt_p(fisher_p)),
                inherit.aes = FALSE, size = 2.9, vjust = 0, fill = NA,
                label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  scale_fill_manual(values = setNames(c(COL_N2, COL_EDIT, COL_EDIT),
                                      levels(n2$line)), guide = "none") +
  geom_text(aes(y = hi + 0.022, label = paste0("n = ", n_plated)),
            size = 2.7, colour = "grey25") +
  ## genotype only, as in panel A. The two edited lines share a genotype
  ## label because they are independent lines of the same edit.
  scale_x_discrete(labels = function(l) sub("  .*$", "", l)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 0.60), expand = expansion(0)) +
  labs(x = NULL, y = "Embryos hatched",
       title = panel_title("B")) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

## ===========================================================================
## D -- sid-2 coding variation in the wild population
##
## The protein runs vertically with residue 1 at the top, the topology as a
## filled bar, and each protein-altering variant labelled to the right with a
## frequency bar. Built from supplemental_data/structure/sid2_variants_cendr.tsv; run
## scripts/sid2_variant_table.R to rebuild it.
##
## COVERAGE. The amino-acid annotation covers variants segregating among the
## four cross parents, so this is not an exhaustive catalogue of sid-2 coding
## variation: 81 variants of any kind segregate in the 3 kb span in CeNDR and
## 8 protein-altering ones are annotated here. The panel says so rather than
## implying completeness.
## ===========================================================================
msg("panel D: wild variation")
stopifnot(file.exists(VAR))
vr <- read_tsv(VAR, comment = "#", show_col_types = FALSE) %>%
  mutate(focal = label == "T96K")
msg("  JU1793 vs JU2466 differ at ", sum(vr$parents_differ), " of ", nrow(vr),
    " sites: ", paste(vr$label[vr$parents_differ], collapse = ", "))
pr2 <- read_tsv(PERRES2, show_col_types = FALSE) %>%
  mutate(topology = factor(topology, levels = names(TOPO_COL2)))
LEN <- max(pr2$resid)
msg("  ", nrow(vr), " protein-altering variants | AF ",
    sprintf("%.3f-%.3f", min(vr$af), max(vr$af)))
msg("  ", paste(sprintf("%s %.0f%%", vr$label, 100 * vr$af), collapse = " | "))

dom2 <- (function(v) { r <- rle(as.character(v))
  tibble(value = factor(r$values, levels = names(TOPO_COL2)),
         end = cumsum(r$lengths),
         start = cumsum(r$lengths) - r$lengths + 1) })(pr2$topology)

## x geometry, in arbitrary units: the topology bar, then the labels, then the
## frequency bars. Residue is on y, reversed so residue 1 is at the top.
X_BAR <- c(0, 0.5); X_LAB <- 0.78; X_FRQ <- c(1.36, 2.55)
## A local-net-charge strip immediately left of the topology bar, on the SAME
## ramp and the SAME limits as panel C, so the two panels can be read against
## each other -- that is the entire point of putting it here. It therefore
## needs no key of its own; panel C's key serves both, and duplicating it would
## invite the two from drifting apart.
X_CHG <- c(-0.86, -0.16)
## Colours are precomputed to hex rather than mapped through a second fill
## scale: the panel already uses fill for the topology, and ggnewscale is not a
## dependency of this repository. The key in panel C does exactly the same.
chg <- read_tsv(LOCALQ, show_col_types = FALSE) %>%
  transmute(resid, q = q_local_pH44,
            col = ramp[pmax(1, pmin(length(ramp),
                     round((pmin(pmax(q, -QLIM), QLIM) + QLIM) /
                           (2 * QLIM) * (length(ramp) - 1)) + 1))])
## The charge is defined only where there is a model to measure it in: the
## AlphaFold ectodomain, residues 21-188. Panel D draws all 311 residues, so
## the strip covers a little over half the protein and the rest is left blank
## rather than filled with a zero that would read as "neutral here".
CHG_RANGE <- range(chg$resid)
stopifnot(nrow(chg) == 168, CHG_RANGE[1] == 21, CHG_RANGE[2] == 188)
X_P1 <- 3.16; X_P2 <- 3.62      # the two cross-parent columns
## Callout rows are evenly spaced down the panel and joined to the residue
## by a leader, rather than sitting at the residue's own height: four of the
## eight variants fall between residues 141 and 153 and their labels and
## frequency bars overlapped completely when placed at true scale.
vr <- vr %>% arrange(residue) %>%
  mutate(row = seq(14, LEN - 14, length.out = n()),
         xend = X_FRQ[1] + af * diff(X_FRQ))

pD <- ggplot() +
  ## the charge strip, and an outline showing how far the model reaches
  ## height slightly over 1 so adjacent residues abut with no hairline seam:
  ## at this scale a seam is a white line, and white is also the middle of a
  ## diverging ramp, so seams would read as neutral charge or as missing data
  geom_tile(data = chg, aes(x = mean(X_CHG), y = resid),
            fill = chg$col, width = diff(X_CHG), height = 1.02) +
  annotate("rect", xmin = X_CHG[1], xmax = X_CHG[2],
           ymin = CHG_RANGE[1] - 0.5, ymax = CHG_RANGE[2] + 0.5,
           fill = NA, colour = "grey45", linewidth = 0.25) +

  geom_rect(data = dom2,
            aes(xmin = X_BAR[1], xmax = X_BAR[2],
                ymin = start - 0.5, ymax = end + 0.5, fill = value),
            colour = "grey30", linewidth = 0.25) +
  ## tick across the topology bar at the true residue, then a leader out to
  ## the evenly spaced callout row
  geom_segment(data = vr,
               aes(x = X_BAR[1], xend = X_BAR[2], y = residue, yend = residue),
               linewidth = 0.5,
               colour = ifelse(vr$focal, COL_JU1793, "grey35")) +
  geom_segment(data = vr,
               aes(x = X_BAR[2], xend = X_LAB - 0.04, y = residue, yend = row),
               linewidth = 0.3,
               colour = ifelse(vr$focal, COL_JU1793, "grey55")) +
  geom_richtext(data = vr, aes(X_LAB, row, label = label),
                colour = ifelse(vr$focal, COL_JU1793, "grey15"),
                fontface = ifelse(vr$focal, "bold", "plain"),
                size = 2.9, hjust = 0, vjust = 0.5, fill = NA,
                label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  ## frequency bars
  geom_segment(data = vr, aes(x = X_FRQ[1], xend = X_FRQ[2], y = row,
                              yend = row),
               linewidth = 3.0, colour = "grey92", lineend = "butt") +
  geom_segment(data = vr, aes(x = X_FRQ[1], xend = xend, y = row,
                              yend = row),
               linewidth = 3.0, lineend = "butt",
               colour = ifelse(vr$focal, COL_JU1793, "grey55")) +
  geom_richtext(data = vr, aes(X_FRQ[2] + 0.06, row,
                               label = sprintf("%.0f%%", 100 * af)),
                colour = ifelse(vr$focal, COL_JU1793, "grey25"),
                size = 2.7, hjust = 0, vjust = 0.5, fill = NA,
                label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  ## the allele each cross parent carries. Rows where they differ are banded,
  ## which is the whole point of the two columns: only two of the eight
  ## protein-altering variants segregate in this cross.
  geom_rect(data = vr %>% filter(parents_differ),
            aes(xmin = X_P1 - 0.22, xmax = X_P2 + 0.22,
                ymin = row - 11, ymax = row + 11),
            fill = COL_JU1793, alpha = 0.10) +
  geom_richtext(data = vr, aes(X_P1, row, label = ju1793_aa),
                colour = ifelse(vr$parents_differ, COL_JU1793, "grey45"),
                fontface = ifelse(vr$parents_differ, "bold", "plain"),
                size = 2.9, hjust = 0.5, vjust = 0.5, fill = NA,
                label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  geom_richtext(data = vr, aes(X_P2, row, label = ju2466_aa),
                colour = ifelse(vr$parents_differ, COL_JU2466, "grey45"),
                fontface = ifelse(vr$parents_differ, "bold", "plain"),
                size = 2.9, hjust = 0.5, vjust = 0.5, fill = NA,
                label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  ## column headers, in the margin reserved above residue 1
  geom_richtext(data = tibble(
      x = c(mean(X_CHG), mean(X_FRQ), X_P1, X_P2), y = -7,
      lab = c("Net charge", "CeNDR frequency", "JU1793", "JU2466"),
      col = c("grey30", "grey30", COL_JU1793, COL_JU2466)),
      aes(x, y, label = lab),
      colour = c("grey30", "grey30", COL_JU1793, COL_JU2466),
      size = 2.6, hjust = 0.5, vjust = 0.5, fill = NA, label.color = NA,
      label.padding = grid::unit(rep(0, 4), "pt")) +
  ## The ends of the charge column, so it reads as a scaled quantity rather than
  ## as decoration; the ramp itself is keyed in panel C. geom_text, NOT
  ## geom_richtext: gridtext parses its label as markdown and a bare "+" becomes
  ## a list item, which fails at tag dispatch.
  geom_text(data = tibble(x = X_CHG, y = c(-2.5, -2.5),
                          lab = c("\u2212", "+")),
            aes(x, y, label = lab), size = 2.6, colour = "grey45",
            hjust = 0.5, vjust = 0.5) +
  scale_fill_manual(values = TOPO_COL2, name = NULL) +
  scale_x_continuous(limits = c(X_CHG[1] - 0.04, X_P2 + 0.30),
                     expand = expansion(0)) +
  ## 311 dropped from the breaks: it collided with the 300 tick
  scale_y_reverse(breaks = c(1, seq(50, 300, 50)),
                  limits = c(LEN + 6, -16), expand = expansion(0)) +
  labs(x = NULL, y = "SID-2 residue", title = panel_title("D"),
       subtitle = paste("Left strip: local net charge on panel C's scale,",
                        "over the modelled ectodomain (21&ndash;188) only")) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "bottom", legend.margin = margin(t = -4),
        legend.text = element_text(size = 7.6),
        plot.margin = margin(4, 4, 4, 4))

## ===========================================================================
## The two experiments on top, the structure across the bottom with room to
## breathe. The structure is a coord_fixed composite of aspect ~1.4, so across
## the full width it is height-limited and sits centred with white space either
## side; the generous row height is what stops it being a thin strip.
##
## Parentheses matter: `p_struct | p_ju / p_n2` would not group as intended, and
## `... + plot_layout(...)` attaches to the last operand unless the whole
## expression is wrapped.
## objects are named for what they contain, not for their letter: the letters
## were reordered once already and a layout written in terms of pA/pB/pC put
## the structure in the top-left without erroring
## The structure and the variant summary share the bottom row: the structure
## composite is wide and short of vertical content, the variant panel is tall
## and narrow, and pairing them uses the space the structure alone left empty.
fig <- ((p_ju | p_n2) / (p_struct | pD)) +
  plot_layout(heights = c(1, 1.5)) &
  theme(plot.tag = element_blank())

ggsave(file.path(OUT, "Figure4_sid2.pdf"), fig, width = 12.2, height = 9.0,
       device = cairo_pdf)
ggsave(file.path(OUT, "Figure4_sid2.png"), fig, width = 12.2, height = 9.0,
       dpi = 300, bg = "white")
msg("wrote Figure4_sid2.{pdf,png}")
