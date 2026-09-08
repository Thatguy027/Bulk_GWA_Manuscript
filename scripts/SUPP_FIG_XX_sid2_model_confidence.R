## Figure 4C, rebuilt --------------------------------------------------------
##
##   Rscript claude_science/scripts/alt_figure4C_plddt.R  -> claude_science/plots/Figure4C_rebuild.{pdf,png}
##
## Two sub-panels, both vector, replacing the two raster renders:
##
##   i   the membrane-oriented ectodomain as a Ca trace coloured by pLDDT, so
##       the reader can see that the scaffold around T96 is a low-confidence
##       model rather than being told so in the caption. Same frame as
##       sid2_membrane_oriented.pdb: +z is the intestinal lumen, up the page.
##
##   ii  the null the proximity claim is measured against: every ectodomain
##       residue's Ca distance from T96, with the four published
##       uptake-critical residues marked. 41.9% of the ectodomain is within
##       20 A, so 3 of 4 is what you expect by chance (binomial p = 0.20).
##       Panel C previously showed only the three near residues and their
##       distances, which reads as evidence; this states the limit in the
##       figure instead of the caption.
##
## Distances reproduce the existing render script exactly:
##   Ca-Ca      D34 13.5, H32 16.4, H168 19.2, H175 37.8 A
##   nearest heavy-atom  9.9, 14.6, 16.9, 36.7 A
##
## Colours and typography follow scripts/Figure4_sid2.R; the pLDDT ramp is the
## canonical AlphaFold four-bin palette, which readers already know.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggtext)
})

RES <- "supplemental_data/structure/sid2_panelC_residues.tsv"
PUB <- "supplemental_data/structure/sid2_panelC_published.tsv"
CARTOON <- "plots/assets/sid2_cartoon_plddt.png"
ZOOM <- "plots/assets/sid2_zoom_plddt.png"
OUT <- "plots"

COL_FOCAL <- "#F34C00"
PLDDT_COL <- c(`very low (<50)` = "#FF7D45", `low (50-70)` = "#FFDB13",
               `confident (70-90)` = "#65CBF3", `very high (>90)` = "#0053D6")

panel_title <- function(letter)
  paste0("<span style='font-size:13pt;color:#111111'>**", letter, "**</span>")

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

res <- read_tsv(RES, show_col_types = FALSE)
## Two classes, deliberately not merged. The uptake-critical set of McEwan,
## Weisman & Hunter 2012 (Mol Cell 47:746) is three extracellular HISTIDINES --
## H32, H168, H175 -- targeted because imidazole protonates only in the acidic
## conditions SID-2 requires. That paper never mentions residue 34; D34 is the
## qt13 loss-of-function allele, a different line of evidence.
pub <- read_tsv(PUB, show_col_types = FALSE) |>
  mutate(cls = factor(cls, levels = c("uptake histidine", "qt13 allele")))
CLS_COL <- c(`uptake histidine` = "#16324A", `qt13 allele` = "#7B5EA7")

bin_plddt <- function(v)
  cut(v, breaks = c(-Inf, 50, 70, 90, Inf), labels = names(PLDDT_COL))
res <- res |> mutate(band = bin_plddt(plddt))

## ---- i. the cartoon -------------------------------------------------------
## The cartoon is rendered by claude_science/scripts/alt_render_plddt_cartoon.py, which reuses the
## ribbon geometry and camera of scripts/sid2_zoom_render.py and colours by
## pLDDT instead of secondary structure. It is embedded the same way the
## released Figure 4 embeds its structure images: annotation_raster on an
## inch-scaled coord_fixed canvas, so the image keeps its own aspect ratio.
im <- png::readPNG(CARTOON)
asp <- dim(im)[2] / dim(im)[1]

H_IM <- 4.05                      # image height, inches
W_IM <- H_IM * asp
KEY_Y <- -0.30                    # the pLDDT key sits under the image
KEY_H <- 0.13

key <- tibble(band = factor(names(PLDDT_COL), levels = names(PLDDT_COL)),
              col = unname(PLDDT_COL),
              ## one column: the panel is narrow and a 2x2 key ran the
              ## "confident (70-90)" label into the swatch beside it
              row = 0:3,
              slot = 0) |>
  mutate(xmin = slot * (W_IM / 2) + 0.02,
         xmax = xmin + KEY_H,
         ymin = KEY_Y - row * (KEY_H + 0.10),
         ymax = ymin + KEY_H)

p_i <- ggplot() +
  annotation_raster(im, xmin = 0, xmax = W_IM, ymin = 0, ymax = H_IM,
                    interpolate = TRUE) +
  geom_rect(data = key, aes(xmin = xmin, xmax = xmax, ymin = ymin,
                            ymax = ymax, fill = band), colour = NA) +
  geom_text(data = key, aes(x = xmax + 0.07, y = (ymin + ymax) / 2,
                            label = band), hjust = 0, size = 2.35,
            colour = "grey25") +
  annotate("text", x = 0.02, y = KEY_Y + KEY_H + 0.14, hjust = 0,
           label = "model confidence (pLDDT)", size = 2.5, colour = "grey30") +
  scale_fill_manual(values = PLDDT_COL, guide = "none") +
  coord_fixed(ratio = 1, xlim = c(-0.05, W_IM + 0.05),
              ylim = c(KEY_Y - 3 * (KEY_H + 0.10) - 0.06, H_IM),
              expand = FALSE, clip = "off") +
  labs(title = panel_title("C"),
       subtitle = sprintf(
         "Ectodomain cartoon, lumenal face up<br>T96 sits at the edge of the modelled core (pLDDT %.0f)",
         res$plddt[res$resid == 96])) +
  theme_void(base_size = 11) +
  theme(plot.title = element_markdown(size = 11.5),
        plot.subtitle = element_markdown(size = 8.3, colour = "grey30"),
        plot.title.position = "plot",
        plot.margin = margin(2, 6, 2, 6))

## ---- ii. the T96 zoom -----------------------------------------------------
## Rendered by claude_science/scripts/alt_render_plddt_zoom.py: the released zoom's camera, scaffold
## radius, sticks and Ca-Ca annotations, with the scaffold carrying pLDDT in
## both hue and opacity so the disordered loop recedes instead of shouting.
## No colour key here -- the key on the cartoon panel serves both.
im_z <- png::readPNG(ZOOM)
asp_z <- dim(im_z)[2] / dim(im_z)[1]
W_Z <- H_IM * asp_z

p_zoom <- ggplot() +
  annotation_raster(im_z, xmin = 0, xmax = W_Z, ymin = 0, ymax = H_IM,
                    interpolate = TRUE) +
  coord_fixed(ratio = 1, xlim = c(-0.05, W_Z + 0.05), ylim = c(0, H_IM),
              expand = FALSE, clip = "off") +
  labs(subtitle = paste0("T96 environment: side chains and C&alpha;-C&alpha; distances for<br>",
                         "the two published residues in range; H175 is 37.8 &#197; away")) +
  theme_void(base_size = 11) +
  theme(plot.subtitle = element_markdown(size = 8.3, colour = "grey30"),
        plot.title.position = "plot",
        plot.margin = margin(2, 6, 2, 6))

## ---- iii. the null --------------------------------------------------------
ecd <- res |> filter(resid != 96)
med <- median(ecd$d_ca_t96)
frac20 <- mean(ecd$d_ca_t96 <= 20)
## the test is on the published uptake set only -- three histidines, of which
## two lie within 20 A. Folding D34 in made it 3 of 4 and p = 0.20.
his <- pub |> filter(cls == "uptake histidine")
bt <- binom.test(sum(his$d_ca <= 20), nrow(his), frac20, alternative = "greater")

## sorted by distance: D34 13.5, H32 16.4, H168 19.2, H175 37.8; heights are
## staggered so the two labels 2.9 A apart on x do not collide
pub_ii <- pub |> arrange(d_ca) |> mutate(y = c(15.6, 10.8, 13.2, 10.8))

p_ii <- ggplot(ecd, aes(d_ca_t96)) +
  annotate("rect", xmin = -Inf, xmax = 20, ymin = -Inf, ymax = Inf,
           fill = "#F34C00", alpha = 0.055) +
  geom_histogram(binwidth = 2, boundary = 0, fill = "grey78",
                 colour = "white", linewidth = 0.25) +
  geom_vline(xintercept = 20, colour = COL_FOCAL, linewidth = 0.4) +
  geom_vline(xintercept = med, colour = "grey35", linewidth = 0.4,
             linetype = "22") +
  geom_point(data = pub_ii, aes(d_ca, y, colour = cls), inherit.aes = FALSE,
             size = 2.4, shape = 21, fill = "white", stroke = 0.7) +
  geom_text(data = pub_ii, aes(d_ca, y + 1.4, label = label, colour = cls),
            inherit.aes = FALSE, size = 2.7, fontface = "bold",
            show.legend = FALSE) +
  scale_colour_manual(values = CLS_COL, name = NULL) +
  annotate("text", x = 19.2, y = 21.5, hjust = 1, size = 2.7,
           colour = COL_FOCAL, label = sprintf("%.0f%% of the ectodomain", 100 * frac20)) +
  annotate("text", x = 20.8, y = 21.5, hjust = 0, size = 2.7,
           colour = "grey35", label = sprintf("median %.1f \u00c5", med)) +
  scale_x_continuous("C&alpha; distance from T96 (\u00c5)",
                     breaks = seq(0, 60, 10), expand = expansion(c(0.01, 0.02))) +
  scale_y_continuous("Ectodomain residues", expand = expansion(c(0, 0.08))) +
  labs(subtitle = sprintf(paste0(
    "The uptake-critical histidines are not closer to T96 than chance:<br>",
    "%d of %d within 20 &#197; (binomial *p* = %.2f)"),
    sum(his$d_ca <= 20), nrow(his), bt$p.value)) +
  theme_pub() +
  theme(axis.title.x = element_markdown(size = 9),
        axis.title.y = element_text(size = 9),
        plot.subtitle = element_markdown(size = 8.3, colour = "grey30"),
        legend.position = "inside",
        legend.position.inside = c(0.99, 0.72),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 7.4),
        legend.background = element_rect(fill = "white", colour = NA))

fig <- (p_i | p_zoom | p_ii) + plot_layout(widths = c(0.60, 0.72, 1.30))

dir.create(OUT, showWarnings = FALSE)
ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_model_confidence.pdf"), fig, width = 12.2, height = 4.9,
       device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_model_confidence.png"), fig, width = 12.2, height = 4.9,
       dpi = 300, bg = "white")

cat(sprintf("ectodomain n=%d  median=%.1f  frac<=20A=%.3f  binom p=%.2f\n",
            nrow(ecd), med, frac20, bt$p.value))
cat("pLDDT bands:\n"); print(table(res$band))
