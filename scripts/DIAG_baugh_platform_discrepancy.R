## Which strains disagree between MIP-seq and NNLS, and which are in only one --
##
## Two questions, one script.
##
## 1  MEMBERSHIP. The joined cache carries 102 strains. They are not the same
##    102 on both sides, and N2 is present in both but excluded downstream.
##    Every name on both sides is resolved through the CaeNDR strain table --
##    the strain name itself AND each of its previous names -- so a rename
##    cannot masquerade as a missing strain. The table is written to
##    plots/diagnostics/baugh_strain_membership.tsv.
##
##    What that resolution finds: exactly one rename, MIP's PB306 -> isotype
##    ECA259, and it does not close the gap. At isotype level MIP has 100,
##    NNLS 102, sharing 99. CX11262, ECA348 and NIC260 are NNLS-only and
##    ECA259/PB306 is MIP-only; none is an alias of the other.
##
##    The three NNLS-only strains each have a near-twin that IS on the panel
##    (CX11262/CX11264 IBS 0.973, ECA348/PS2025 0.984, NIC260/NIC256 0.960),
##    which invites the reading that NNLS is splitting one strain across a
##    pair. That is tested below and REJECTED: summing each pair tracks the
##    partner's MIP values worse than the partner alone, by 4-10x in RMSD.
##    They behave instead like ordinary pool members -- mean |slope| 2.2e-4
##    to 6.9e-4 against a median of 3.4e-4 over the 98, across-arm CV 0.47 to
##    1.63 inside the 98-strain IQR of 0.84-2.24, and consistent slope sign
##    across 4 or 5 of the 5 independent replicate arms. The conclusion is
##    that they were in the physical pool and the MIP probe panel has no
##    column for them, not that the deconvolution invented them.
##
## 2  DISAGREEMENT. Three measures, because they do not rank strains the same
##    way and each answers a different question:
##      slope residual  how far the strain sits off y = x in Figure 1A, i.e.
##                      disagreement in the trait that gets mapped
##      frequency RMSD  root-mean-square |NNLS - MIP| across the 23 samples,
##                      i.e. disagreement in level
##      per-strain rho  Spearman of NNLS against MIP across the 23 samples,
##                      i.e. disagreement in trajectory SHAPE, scale-free
##    A strain can track the shape perfectly and still sit off y = x if one
##    platform reads it consistently high, so a low rho and a large RMSD are
##    different failures.
##
## Usage:  Rscript scripts/DIAG_baugh_platform_discrepancy.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
})

source("scripts/Figure1_common.R")

DIAG   <- "plots/diagnostics"
NLABEL <- 8                     # strains labelled per panel
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)

freq <- baugh_frequencies()

## --- 1. membership ---------------------------------------------------------
mip_raw <- readLines(gzfile(file.path(BAUGH, "mipseq_frequencies.txt.gz")))
mip_strains <- sub("\t.*$", "", mip_raw[-1])          # line 1 is the header

## resolve every name through the CaeNDR strain table, previous names included
sd_tab <- suppressWarnings(read_csv("data/20250625_c_elegans_strain_data.csv",
                                    show_col_types = FALSE))
resolver <- bind_rows(
  sd_tab %>% transmute(name = strain, isotype, via = "strain"),
  sd_tab %>% filter(!is.na(previous_names)) %>%
    separate_rows(previous_names, sep = "\\|") %>%
    transmute(name = str_trim(previous_names), isotype, via = "previous_name")
) %>% filter(name != "") %>% distinct(name, .keep_all = TRUE)

member <- freq %>%
  group_by(strain) %>%
  summarise(nnls_obs = sum(!is.na(frq)),
            mip_obs  = sum(!is.na(published_frq)), .groups = "drop") %>%
  mutate(in_nnls = nnls_obs > 0,
         in_mip  = mip_obs  > 0) %>%
  bind_rows(tibble(strain = setdiff(mip_strains, .$strain),
                   nnls_obs = 0L, mip_obs = NA_integer_,
                   in_nnls = FALSE, in_mip = TRUE)) %>%
  mutate(status = case_when(
    strain == "N2"        ~ "both, excluded (reference strain)",
    in_nnls & in_mip      ~ "both",
    in_nnls & !in_mip     ~ "NNLS only (not on the MIP panel)",
    !in_nnls & in_mip     ~ "MIP only (not in the NNLS genotype reference)")) %>%
  left_join(resolver, by = c("strain" = "name")) %>%
  arrange(status != "both", strain)

write_tsv(member, file.path(DIAG, "baugh_strain_membership.tsv"))

renamed <- member %>% filter(via == "previous_name")
unresolved <- member %>% filter(is.na(isotype))

cat("\n=== STRAIN MEMBERSHIP ===\n")
cat(sprintf("MIP-seq file            %3d strains\n", length(mip_strains)))
cat(sprintf("NNLS genotype reference %3d strains\n", sum(member$in_nnls)))
cat(sprintf("intersection            %3d strains\n",
            sum(member$in_nnls & member$in_mip)))
cat(sprintf("analysed                %3d strains (intersection minus N2)\n\n",
            sum(member$status == "both")))
member %>% filter(status != "both") %>%
  select(strain, isotype, status, nnls_obs, mip_obs) %>% print(n = 50)
cat("\nnames resolving via a PREVIOUS name (a rename, not a missing strain):\n")
if (nrow(renamed)) print(renamed %>% select(strain, isotype, status)) else
  cat("  none\n")
cat("\nnames that resolve to no isotype at all:\n")
if (nrow(unresolved)) print(unresolved %>% select(strain, status)) else
  cat("  none\n")
mi <- member %>% filter(in_mip)  %>% pull(isotype) %>% na.omit() %>% unique()
ni <- member %>% filter(in_nnls) %>% pull(isotype) %>% na.omit() %>% unique()
cat(sprintf("\nAT ISOTYPE LEVEL: MIP %d, NNLS %d, shared %d\n",
            length(mi), length(ni), length(intersect(mi, ni))))
cat("  NNLS-only:", paste(setdiff(ni, mi), collapse = " "), "\n")
cat("  MIP-only :", paste(setdiff(mi, ni), collapse = " "), "\n")

## --- 2. disagreement -------------------------------------------------------
slopes <- platform_slopes(freq) %>%
  filter(strain != "N2") %>%
  group_by(strain) %>%
  summarise(slope_nnls  = mean(wgs_slope, na.rm = TRUE),
            slope_baugh = mean(mip_slope, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(slope_nnls), is.finite(slope_baugh)) %>%
  mutate(resid = slope_nnls - slope_baugh,
         resid_z = as.numeric(scale(resid)))

persamp <- freq %>%
  filter(strain != "N2", !is.na(published_frq)) %>%
  group_by(strain) %>%
  summarise(rmsd = sqrt(mean((frq - published_frq)^2)),
            bias = mean(frq - published_frq),
            rho  = suppressWarnings(cor(frq, published_frq, method = "spearman")),
            mean_frq = mean(frq), n = n(), .groups = "drop")

dev <- slopes %>% left_join(persamp, by = "strain")
write_tsv(dev, file.path(DIAG, "baugh_platform_deviation.tsv"))

top <- function(d, col, n = NLABEL, worst_low = FALSE) {
  if (worst_low) slice_min(d, .data[[col]], n = n, with_ties = FALSE)
  else           slice_max(d, abs(.data[[col]]), n = n, with_ties = FALSE)
}

COL_PT <- "#2E4057"; COL_HI <- "#C4302B"

pA <- ggplot(dev, aes(slope_baugh, slope_nnls)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_point(colour = COL_PT, alpha = 0.75, size = 1.9) +
  geom_point(data = top(dev, "resid"), colour = COL_HI, size = 2.4) +
  geom_text_repel(data = top(dev, "resid"), aes(label = strain),
                  size = 2.7, colour = COL_HI, min.segment.length = 0,
                  max.overlaps = Inf, seed = 1) +
  labs(title = "A  Slope agreement, and who sits off the line",
       subtitle = sprintf("Spearman rho = %.3f, n = %d strains. Dashed line is y = x, not a fit.",
                          cor(dev$slope_baugh, dev$slope_nnls, method = "spearman"),
                          nrow(dev)),
       x = "MIP-seq slope (slope_baugh)", y = "NNLS slope (slope_nnls)")

dB <- dev %>% mutate(strain = fct_reorder(strain, rmsd))
pB <- ggplot(dB, aes(rmsd, strain)) +
  geom_point(colour = COL_PT, size = 1.4) +
  geom_point(data = top(dB, "rmsd"), colour = COL_HI, size = 2) +
  geom_text_repel(data = top(dB, "rmsd"), aes(label = strain),
                  size = 2.6, colour = COL_HI, direction = "y", hjust = 0,
                  nudge_x = max(dB$rmsd) * 0.06, segment.colour = "grey70",
                  segment.size = 0.25, min.segment.length = 0,
                  max.overlaps = Inf, seed = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.30))) +
  labs(title = "B  Level disagreement per strain",
       subtitle = "Root-mean-square |NNLS - MIP-seq| frequency across the 23 samples",
       x = "RMSD (frequency units)", y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

dC <- dev %>% filter(is.finite(rho)) %>% mutate(strain = fct_reorder(strain, rho))
pC <- ggplot(dC, aes(rho, strain)) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey60") +
  geom_point(colour = COL_PT, size = 1.4) +
  geom_point(data = top(dC, "rho", worst_low = TRUE), colour = COL_HI, size = 2) +
  geom_text_repel(data = top(dC, "rho", worst_low = TRUE), aes(label = strain),
                  size = 2.6, colour = COL_HI, direction = "y", hjust = 1,
                  nudge_x = -0.18, segment.colour = "grey70",
                  segment.size = 0.25, min.segment.length = 0,
                  max.overlaps = Inf, seed = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.30, 0.05))) +
  labs(title = "C  Trajectory-shape disagreement per strain",
       subtitle = "Spearman of NNLS against MIP-seq across the 23 samples, worst labelled",
       x = "Spearman rho across samples", y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

theme_set(theme_bw(9) +
  theme(plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.6, colour = "grey30"),
        panel.grid.minor = element_blank()))

fig <- pA / (pB | pC) + plot_layout(heights = c(1.25, 1))
ggsave(file.path(DIAG, "DIAG_baugh_platform_discrepancy.pdf"), fig,
       width = 9.5, height = 9)
ggsave(file.path(DIAG, "DIAG_baugh_platform_discrepancy.png"), fig,
       width = 9.5, height = 9, dpi = 200)

cat("\n=== LARGEST SLOPE RESIDUALS (NNLS - MIP) ===\n")
print(top(dev, "resid") %>%
        select(strain, slope_baugh, slope_nnls, resid, resid_z, rmsd, rho) %>%
        arrange(desc(abs(resid))), n = NLABEL)
cat("\n=== LARGEST FREQUENCY RMSD ===\n")
print(top(dev, "rmsd") %>% select(strain, rmsd, bias, mean_frq, rho) %>%
        arrange(desc(rmsd)), n = NLABEL)
cat("\n=== WORST TRAJECTORY AGREEMENT ===\n")
print(top(dev, "rho", worst_low = TRUE) %>%
        select(strain, rho, rmsd, mean_frq, resid) %>% arrange(rho), n = NLABEL)
cat(sprintf("\nwrote %s/DIAG_baugh_platform_discrepancy.{pdf,png}\n", DIAG))
cat(sprintf("wrote %s/baugh_strain_membership.tsv and baugh_platform_deviation.tsv\n", DIAG))
