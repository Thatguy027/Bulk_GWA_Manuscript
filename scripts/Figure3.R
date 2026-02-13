library(tidyverse)
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))


# # Bottom panel
library(ggplot2)
library(patchwork)
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))

# Read BED file
bed <- read.table("../relevantNILs/nil_ranges.bed", header = FALSE,
                  col.names = c("chr", "start", "end", "szname", "geno")) %>%
  dplyr::mutate(background_start = 1,
                background_end = start,
                background_geno = ifelse(geno == "JU1793", "JU2466",
                                         ifelse(geno == "JU2466", "JU1793", NA)))

# nil phenotypes
nil_p <- data.table::fread("../relevantNILs/RNAi Sensitivity - JU1793 - November2025 gSZ177 - gSZ179 NILs.tsv") %>%
  dplyr::filter(Strain %in% c("JU2466","JU1793","","wSZ153","wSZ176","wSZ191","wSZ196")) %>%
  dplyr::mutate(condition = ifelse(condition == "ht115","HT115", "pos-1"))

# ---- 1) strain set + order from the RNAi plot ----
strain_levels <- nil_p %>%
  dplyr::filter(Strain %in% c("JU2466","JU1793","wSZ153","wSZ176","wSZ191","wSZ196")) %>%  # (dropped "" — add back if you truly have it)
  dplyr::distinct(Strain) %>%
  dplyr::pull(Strain)

# if you want the same visual "top-to-bottom" order as your barplot (rev),
# we’ll use this reversed order for BOTH plots:
strain_levels_rev <- rev(strain_levels)

nil_p2 <- nil_p %>%
  dplyr::filter(Strain %in% strain_levels) %>%
  dplyr::mutate(
    Strain = factor(Strain, levels = strain_levels_rev)
  ) %>%
  dplyr::filter(condition!="HT115")

# drop wSZ153 from the ordering + both datasets
strain_levels2     <- setdiff(strain_levels, "wSZ153")
strain_levels_rev2 <- rev(strain_levels2)

cutoff  <- 13.6e6
chr_end <- 13783801
base_size <- 18

# ---- 1) strain set + order from the RNAi plot ----
strain_levels <- nil_p %>%
  dplyr::filter(Strain %in% c("JU2466","JU1793","wSZ153","wSZ176","wSZ191","wSZ196")) %>%  # (dropped "" — add back if you truly have it)
  dplyr::distinct(Strain) %>%
  dplyr::pull(Strain)

# if you want the same visual "top-to-bottom" order as your barplot (rev),
# we’ll use this reversed order for BOTH plots:
strain_levels_rev <- rev(strain_levels)

# drop wSZ153 from the ordering + both datasets
strain_levels2     <- setdiff(strain_levels, "wSZ153")
strain_levels_rev2 <- rev(strain_levels2)

# --- BED plot (left) ---
bed_filt <- bed %>%
  dplyr::filter(szname %in% strain_levels2) %>%
  dplyr::mutate(
    szname = factor(szname, levels = strain_levels_rev2),
    y = as.numeric(szname),
    background_start_cutoff = cutoff,
    background_end_chr      = chr_end,
    start_cutoff = if_else(szname %in% c("JU1793","JU2466"), cutoff, start)
  )

p_geno <- ggplot(bed_filt) +
  geom_rect(
    aes(xmin = background_start_cutoff/1e6, xmax = background_end_chr/1e6,
        ymin = y - 0.4, ymax = y + 0.4, fill = background_geno),
    color = "black", alpha = 1
  ) +
  geom_rect(
    aes(xmin = start_cutoff/1e6, xmax = end/1e6,
        ymin = y - 0.4, ymax = y + 0.4, fill = geno),
    color = "black", alpha = 1
  ) +
  geom_vline(xintercept = c(13.657700, 13.695), linewidth = 0.7) +
  scale_fill_manual(values = c(JU1793 = "red", JU2466 = "blue")) +
  coord_cartesian(xlim = c(cutoff/1e6, chr_end/1e6)) +
  theme_bw(base_size = base_size) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(x = "Genomic Position (Mb)", fill = NULL)

# --- Bar plot (right) ---
pd <- position_dodge(width = 0.9)

nil_p2 <- nil_p %>%
  dplyr::filter(Strain %in% strain_levels2) %>%
  dplyr::mutate(Strain = factor(Strain, levels = strain_levels_rev2)) %>%
  dplyr::filter(condition != "HT115")

p_hatch <- nil_p2 %>%
  dplyr::mutate(
    fraction_hatched = (`plated embryo` - unhatched) / `plated embryo`,
    bi_error_hatched = sqrt((fraction_hatched * (1 - fraction_hatched)) / (`plated embryo`)) * 1.96
  ) %>%
  ggplot(aes(x = fraction_hatched, y = Strain)) +
  geom_col(position = pd, fill = "gray70") +
  geom_errorbar(
    aes(xmin = fraction_hatched - bi_error_hatched,
        xmax = fraction_hatched + bi_error_hatched),
    width = 0.25, colour = "black", alpha = 0.9, linewidth = 0.5,
    position = pd
  ) +
  scale_y_discrete(limits = strain_levels_rev2) +
  theme_bw(base_size = base_size) +   # <-- match left plot border
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  xlab("Fraction hatched")

result_file <- "../cross_experiments/Nov2024_JU_cross_pos_mig/plots/JU1793_JU2466_F2-2_contrast_HT115g-POS1g_10000_plot_DF.tsv"

result_df <- data.table::fread(result_file) %>%
  dplyr::select(chrom, physical.position, p)

effective.n.tests <- 2000
sig_line <- -log10(0.05 / effective.n.tests)

p_manhattan <- ggplot(result_df, aes(x = physical.position/1e6, y = -log10(p))) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_hline(yintercept = sig_line, linewidth = 0.6, linetype = "dashed") +
  facet_grid(. ~ chrom, scales = "free_x", space = "free_x") +
  labs(x = "Genomic Position (Mb)", y = expression(-log[10](p))) +
  theme_bw(base_size = base_size) +
  theme(legend.position = "none")

p_geno <- p_geno + theme(legend.position = "none")
p_hatch <- p_hatch + theme(legend.position = "none")
p_manhattan <- p_manhattan + theme(legend.position = "none")

# reassemble
# Add panel tags: A for top row, B for bottom row (upper-left of each row)
(p_manhattan / (p_geno | p_hatch)) +
  plot_layout(heights = c(1, 1), widths = c(1, 2)) +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(face = "bold"))
