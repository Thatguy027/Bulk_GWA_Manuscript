## Supplement -- the simulation that started the project ---------------------
##
##   Rscript scripts/SUPP_FIG_XX_simulation_depth.R
##     -> plots/SUPP_FIG_XX_simulation_depth.{pdf,png}
##
##   A  reported r-squared of NNLS-estimated against known input frequency,
##      against simulated sequencing depth, one line per simulated trait
##   B  the same estimates plotted against the 500x estimate, faceted by depth
##
## THE SIMULATION. Each wild isolate was assigned a fitness value drawn from an
## inverse chi-squared distribution; the expected pooled allele frequencies such
## a population would produce were computed; observed alt-allele counts were
## simulated by binomial sampling at depths 1, 3, 5, 10, 30, 50, 100 and 500x;
## and those counts were deconvolved back to per-strain frequencies by NNLS.
## Seven published C. elegans traits with validated QTL supplied the fitness
## structure: Albendazole_q75.TOF, PC1, assay_norm, mtDNA_ratio, value,
## amsacrine_f.L1 and etoposide_median.TOF.
##
## WHAT IS ARCHIVED AND WHAT IS NOT -- read before quoting anything here
## ---------------------------------------------------------------------
## Only the NNLS OUTPUT survives. data/experiments/initial_sims holds the
## estimated coefficients (nnls_estFreqs.RDS) and the original per-trait plots.
## The simulation script itself, the drawn fitness values and the "expected"
## input frequencies are NOT in the archive, and the trait directory the
## original script reads (../traits_with_validated_qtl/) is not either.
##
## Consequently:
##   - Panel A CANNOT be recomputed. Its r-squared values are transcribed from
##     the text embedded in the original per-trait PDFs, which is the only
##     surviving record of the comparison against the known input. They are in
##     supplemental_data/deconvolution/simulation_reported_r2.tsv, and that file
##     names its own provenance.
##   - Panel B IS recomputed here, but against the 500x estimate rather than
##     against the truth, because the truth is not archived. It is a convergence
##     panel, not an accuracy panel, and is labelled as such. At 500x the
##     reported r-squared is 1.00 for all seven traits, which is what makes the
##     500x estimate a usable stand-in -- but it is a stand-in.
##
## To restore panel A as a recomputed figure, the expected input frequencies
## (327 strains x 7 traits) and the simulation script are what is needed.
##
## SCALE. The stored coefficients are not frequencies: each depth's row sums to
## the depth itself (500x sums to 500.03, 1x to 1.00), so they are on an
## expected-count scale. frequency = coefficient / sum(coefficient) within a
## trait and depth, which is how the deposit's frequency column was built.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggtext)
})

OUT <- "plots"
DEC <- "supplemental_data/deconvolution"
R2  <- file.path(DEC, "simulation_reported_r2.tsv")
FRQ <- file.path(DEC, "simulation_nnls_frequencies.tsv.gz")

stopifnot(file.exists(R2), file.exists(FRQ))

## a muted qualitative set; seven traits need seven distinguishable hues
TRAIT_COL <- c(
  `PC1`                  = "#1B6C7A",
  `value`                = "#3E8E5A",
  `amsacrine_f.L1`       = "#7A9A3B",
  `assay_norm`           = "#C08A2E",
  `etoposide_median.TOF` = "#B5623C",
  `Albendazole_q75.TOF`  = "#9E4257",
  `mtDNA_ratio`          = "#6A5A8C")

panel_title <- function(letter)
  paste0("<span style='font-size:13pt;color:#111111'>**", letter, "**</span>")

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = base_size - 1),
          plot.title = element_markdown(size = base_size + 0.5),
          plot.subtitle = element_markdown(size = base_size - 2.5,
                                           colour = "grey30"),
          plot.title.position = "plot",
          legend.key.size = grid::unit(9, "pt"))
}

DEPTHS <- c(1, 3, 5, 10, 30, 50, 100, 500)

## ===========================================================================
## A -- reported accuracy against the known input
## ===========================================================================
r2 <- read_tsv(R2, show_col_types = FALSE) %>%
  mutate(trait = factor(trait, levels = names(TRAIT_COL)))
stopifnot(nrow(r2) == 56, !anyNA(r2$trait))

at1 <- r2 %>% filter(depth == 1) %>% arrange(desc(r2))
cat("== reported r-squared at 1x, the depth the text claims ==\n")
print(as.data.frame(at1 %>% transmute(trait, r2)), row.names = FALSE)
cat(sprintf("  range %.2f-%.2f, median %.2f\n\n",
            min(at1$r2), max(at1$r2), median(at1$r2)))

## the depth at which every trait first reaches 0.95 and stays there
reach <- r2 %>% arrange(trait, depth) %>% group_by(trait) %>%
  summarise(first95 = {
    ok <- r2 >= 0.95
    idx <- which(ok & rev(cumprod(rev(ok))) == 1)[1]   # first depth from which it never drops
    if (is.na(idx)) NA_integer_ else depth[idx]
  }, .groups = "drop")
cat("== lowest depth from which r-squared stays >= 0.95 ==\n")
print(as.data.frame(reach), row.names = FALSE)
cat("  worst trait needs ", max(reach$first95, na.rm = TRUE), "x\n\n", sep = "")

lab <- r2 %>% filter(depth == 1)
pA <- ggplot(r2, aes(depth, r2, colour = trait)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.35,
             colour = "grey55") +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.5) +
  annotate("text", x = 1.05, y = 0.958, label = "r² = 0.95", hjust = 0,
           size = 2.7, colour = "grey40") +
  scale_x_log10(breaks = DEPTHS, labels = paste0(DEPTHS, "×")) +
  scale_y_continuous(limits = c(0.5, 1.005), breaks = seq(0.5, 1, 0.1)) +
  scale_colour_manual(values = TRAIT_COL, name = NULL) +
  labs(x = "Simulated sequencing depth", y = "r² vs known input frequency",
       title = panel_title("A")) +
  theme_pub() +
  theme(legend.position = c(0.985, 0.02), legend.justification = c(1, 0),
        legend.text = element_text(size = 7.6),
        legend.background = element_rect(fill = alpha("white", 0.85),
                                         colour = NA),
        panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey92"))

## ===========================================================================
## B -- convergence on the 500x estimate
## ===========================================================================
frq <- read_tsv(FRQ, show_col_types = FALSE)
stopifnot(nrow(frq) == 7 * 8 * 327)

ref <- frq %>% filter(depth == 500) %>% select(trait, strain, ref = frequency)
cmp <- frq %>% filter(depth != 500) %>%
  inner_join(ref, by = c("trait", "strain")) %>%
  mutate(depth_lab = factor(paste0(depth, "×"),
                            levels = paste0(setdiff(DEPTHS, 500), "×")))

conv <- cmp %>% group_by(depth) %>%
  summarise(r2 = cor(frequency, ref)^2, .groups = "drop") %>% arrange(depth)
cat("== convergence on the 500x estimate (derived, all traits pooled) ==\n")
print(as.data.frame(conv %>% mutate(r2 = round(r2, 3))), row.names = FALSE)
cat("\n")

ann <- conv %>% mutate(depth_lab = factor(paste0(depth, "×"),
                                          levels = levels(cmp$depth_lab)),
                       lab = sprintf("r² = %.2f", r2))

pB <- ggplot(cmp, aes(frequency, ref)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.35, colour = "grey60") +
  geom_point(shape = 16, size = 0.5, alpha = 0.25, colour = "grey20") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.12, vjust = 1.5, size = 2.7, colour = "grey20") +
  ## Both axes are the same quantity, so they must share limits: with
  ## scales = "free_x" the y = x line lands somewhere different in every facet
  ## and the comparison it is there to support becomes unreadable.
  facet_wrap(~ depth_lab, nrow = 1) +
  scale_x_continuous(breaks = c(0, 0.01, 0.02, 0.03),
                     labels = c("0", ".01", ".02", ".03")) +
  scale_y_continuous(breaks = scales::pretty_breaks(4)) +
  labs(x = "Estimated strain frequency at the stated depth",
       y = "Estimate at 500×",
       title = panel_title("B")) +
  theme_pub() +
  theme(panel.spacing.x = grid::unit(7, "pt"),
        axis.text.x = element_text(size = 7.2))

fig <- pA / pB + plot_layout(heights = c(1.3, 1))

ggsave(file.path(OUT, "SUPP_FIG_XX_simulation_depth.pdf"), fig,
       width = 9.2, height = 6.4, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_simulation_depth.png"), fig,
       width = 9.2, height = 6.4, dpi = 300, bg = "white")
cat("wrote SUPP_FIG_XX_simulation_depth.{pdf,png}\n")

## the caveat, printed so it cannot be missed by anyone re-running this
neg <- frq %>% filter(coefficient < 0)
cat(sprintf("\nNOTE: %d of %d coefficients (%.1f%%) are negative, at depths %s.\n",
            nrow(neg), nrow(frq), 100 * nrow(neg) / nrow(frq),
            paste(sort(unique(neg$depth)), collapse = ", ")))
cat("A strict non-negative solver cannot return these, so the archived values\n")
cat("carry solver or post-processing noise at the 1e-2 level. See the caption.\n")
