library(tidyverse)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))


result_file <- "plots/JU1793_JU2466_F2-2_contrast_HT115g-POS1g_10000_plot_DF.tsv"

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
