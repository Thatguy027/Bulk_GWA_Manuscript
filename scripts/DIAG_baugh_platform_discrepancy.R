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
##    WHY ARE THEY IN THE REFERENCE AT ALL? They should not be, on the face of
##    it: the pool is defined by what MIP measured, so the deconvolution should
##    solve only for pool strains. The reference is indeed NOT pool-restricted.
##    The source genotype matrix
##      /Users/Stefan/UCLA/Projects/bulkGWAS/baugh_wgs/cluster_data/
##        20220908_Baugh_BulkL1_Bootstrap_Input_flippedCommon_NAfix.RData
##    is 1,240,618 markers x 103 strains, and exactly 100 of those 103 are in
##    the MIP pool definition. The 3 extras are CX11262, ECA348 and NIC260.
##    (That file also contains PB306, which the shipped cache does not, so it
##    is a near relative of the deposited input rather than the same file --
##    the repo loads 2024bootstrapINPUT.Rdata, which is Dryad-only and carries
##    102 strains. The source project holds several variants: flippedCommon,
##    dupNswitch and cleanGenotypes, each with and without an NA fix.)
##
##    WHETHER THE EXTRAS SHOULD BE DROPPED was settled by re-fitting NNLS twice
##    from that matrix, once with all 103 strains and once restricted to the
##    100 the MIP panel measured, and scoring both against MIP.
##
##    The 3 extras hold 3.21% of the pool mass. Freeing it does NOT spread it
##    evenly -- an even split over 100 strains would be +0.00032 each, while
##    PS2025 gains +0.0095 (30x), CX11264 +0.0033 and NIC256 +0.0023. All three
##    near-twins (ECA348/PS2025 IBS 0.984, CX11262/CX11264 0.973,
##    NIC260/NIC256 0.960) land in the top five gainers and together absorb 52%
##    of the freed mass. Taken alone that looks like the extras had been
##    splitting their twins abundance.
##
##    MIP says otherwise, and unanimously. Every one of the three twins gets
##    WORSE on every measure when its extra is removed:
##
##       twin     extra      RMSD full -> restricted   bias full -> restricted
##       PS2025   ECA348     0.00271 -> 0.01210        -0.00164 -> +0.00787
##       CX11264  CX11262    0.00138 -> 0.00288        -0.00096 -> +0.00237
##       NIC256   NIC260     0.01690 -> 0.01950        +0.01560 -> +0.01790
##
##    The bias sign flip is the point: under the full reference the twins read
##    slightly BELOW their independently measured MIP values, and under the
##    restricted reference they read above. Forcing the extras mass onto the
##    nearest genetic neighbour overshoots what MIP actually measured, so the
##    extras are not holding their twins signal -- they are holding material
##    the twin does not account for. And restricting the reference is worse on
##    every properly computed statistic, not neutral:
##
##                        slope rho     cell rho    cell RMSD (N2 excluded)
##       full (103)         0.9745       0.8243          0.00473
##       restricted (100)   0.9500       0.7945          0.00510
##
##    The full-reference slope rho of 0.9745 reproduces the shipped cache's
##    0.974, so the re-fit is scoring the same thing Figure 1A does.
##
##    TWO STATISTICS, DO NOT CONFUSE THEM. The slope rho of 0.974 is one value
##    per strain -- the slope of frequency change across days, averaged over
##    five replicate arms -- which is Figure 1A and the trait that gets mapped.
##    The cell rho of 0.824 is every strain x sample cell scored raw, and it
##    matches the median per-sample agreement of 0.84 that Figure 1B already
##    reports. The gap is averaging: a slope pools fifteen measurements, which
##    removes most of the per-sample noise, exactly as the depth supplement
##    states. Both numbers are real and they describe different things.
##
##    ALSO SUPERSEDED: an earlier version of this header quoted RMSD
##    0.05361 -> 0.05346 and called the change negligible. That RMSD was 99.2%
##    N2 and carried no information about the question. MIP puts N2 at mean
##    frequency 0.5726 against NNLS's 0.0389, so N2's 23 cells dominated the
##    squared error entirely. N2 is excluded from every downstream analysis for
##    the identifiability reason already documented, and excluding it here gives
##    the 0.00473 above.
##
##    WHY N2's MIP VALUE IS ON A DIFFERENT SCALE. The MIP table is not a simplex
##    over its 100 strains. Per-sample sums run 1.558 to 1.686 (mean 1.610),
##    with N2 alone at 0.540 to 0.631 in every sample. Drop N2 and the 99 wild
##    isolates sum to 0.982 to 1.069, mean 1.036 -- a near-simplex. So the N2
##    column is an extra ~0.57 layered on top rather than a pool share on the
##    same scale as the others, which is what produced its RMSD of 0.534 against
##    an NNLS value of 0.039. Excluding N2 is therefore right for a second,
##    concrete reason beyond the identifiability argument. The residual 3.6%
##    excess over 1 in the remaining 99 is unexplained and worth a look before
##    the MIP frequencies are used on their own scale anywhere.
##
##    CONCLUSION: the reference is not pool-restricted, which is worth fixing
##    upstream, but restricting it is not the fix -- by the independent MIP
##    yardstick it makes the three affected strains substantially worse. The
##    reading is that the three correspond to real material in the pool that
##    the MIP panel has no column for.
##
##    SUPERSEDED: an earlier version of this header rejected the splitting
##    hypothesis by summing each pair and comparing against the partner alone.
##    That test was mis-specified -- it scored against the shipped cache rather
##    than re-fitting, and assumed the mass was conserved within the pair when
##    NNLS renormalises across all strains. The re-fit above is the right test
##    and reaches the same conclusion on much better evidence.
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
