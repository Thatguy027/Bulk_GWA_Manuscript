## Supplement -- the full NIL hatching experiment ---------------------------
## Every strain and both food conditions from the experiment behind Figure 3C,
## so the HT115 control arm is on the record: hatching is 97-100% for every
## strain on control food, which is what makes the pos-1 differences in
## Figure 3C attributable to the RNAi rather than to the introgressions.
##
##   Rscript scripts/SUPP_FIG_XX_nil_hatching_full.R
##
## Figure 3C shows five strains on pos-1 only. This adds the four further NILs
## measured in the same experiment (wSZ192-wSZ195) and wSZ153, and both
## conditions for all of them. Strains are labelled here, unlike Figure 3,
## because the point of a supplement is to be readable on its own.
##
## Strains are ordered by their pos-1 hatching, most resistant at the top, which
## lays out the allelic series rather than a nominal order.
##
## RNAi DOSE: 50% pos-1 RNAi bacteria, diluted with HT115. The dose is NOT
## recorded in the source data file -- the condition column carries only
## "pos"/"ht115" -- so it comes from the lab record, not from the data. Figure
## 4B uses 25% and its hatching percentages are therefore not comparable with
## these; see FIGURE_CAPTIONS.txt.
##
## CAVEAT: one plate per strain per condition. Intervals are Wilson binomial
## intervals on that plate's embryo count, so they describe counting
## uncertainty, not between-plate variability. No strain is replicated, so the
## ordering below should be read as a series, not as a set of tested contrasts.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggtext)
})

OUT  <- "plots"
PHEN <- paste0("data/plate_rnai_phenotyping/",
               "RNAi Sensitivity - JU1793 - November2025 gSZ177 - gSZ179 NILs.tsv")
NILB <- "data/nil_ranges.bed"

COL_CTRL <- "grey80"
COL_RNAI <- "grey40"

z <- qnorm(0.975)

phen <- fread(PHEN) %>% as_tibble() %>%
  transmute(strain = Strain,
            cond = ifelse(condition == "ht115", "HT115", "*pos-1*"),
            n = `plated embryo`, hatched = `plated embryo` - unhatched) %>%
  mutate(p = hatched / n,
         ## Wilson, which behaves near 0 and 1 where the Wald interval does not
         lo = pmax(0, (p + z^2 / (2 * n) -
                         z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)),
         hi = pmin(1, (p + z^2 / (2 * n) +
                         z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)))

## order by pos-1 hatching, most resistant at the top
ord <- phen %>% filter(cond == "*pos-1*") %>% arrange(p) %>% pull(strain)
phen <- phen %>%
  mutate(strain = factor(strain, levels = ord),
         cond = factor(cond, levels = c("HT115", "*pos-1*")))

## which strains have a mapped introgression, for the caption
bed <- fread(NILB, header = FALSE,
             col.names = c("chr", "start", "end", "strain", "geno")) %>% as_tibble()
mapped <- intersect(levels(phen$strain), bed$strain)
unmapped <- setdiff(levels(phen$strain), bed$strain)

cat("strains: ", nlevels(phen$strain),
    " | with an introgression in nil_ranges.bed: ", length(mapped),
    " | without: ", paste(unmapped, collapse = ", "), "\n", sep = "")

pd <- position_dodge(width = 0.68)
fig <- ggplot(phen, aes(p, strain, fill = cond)) +
  geom_col(position = pd, width = 0.66, colour = "grey25", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), position = pd, orientation = "y",
                width = 0.26, linewidth = 0.4, colour = "grey15") +
  geom_text(aes(x = 0.015, label = paste0("n=", n)), position = pd,
            hjust = 0, size = 2.5, colour = "grey25") +
  scale_fill_manual(values = c(HT115 = COL_CTRL, `*pos-1*` = COL_RNAI), name = NULL) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.04), expand = expansion(mult = 0)) +
  labs(x = "Embryos hatched", y = NULL,
       title = "NIL hatching, full experiment",
       subtitle = paste(strwrap(paste0(
         "Every strain and both food conditions. Control hatching is 97-100% ",
         "throughout, so the pos-1 differences are not a property of the ",
         "introgressions themselves. Bars are one plate per strain per ",
         "condition with Wilson 95% intervals on that plate's embryo count; ",
         "no strain is replicated. Strains ordered by pos-1 hatching."),
         width = 96), collapse = "\n")) +
  theme_classic(base_size = 11.5) +
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        plot.title.position = "plot",
        legend.position = "top",
        legend.text = element_markdown(size = 10),
        legend.margin = margin(b = -4),
        legend.key.size = grid::unit(10, "pt"))

ggsave(file.path(OUT, "SUPP_FIG_XX_nil_hatching_full.pdf"), fig,
       width = 8, height = 6, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_nil_hatching_full.png"), fig,
       width = 8, height = 6, dpi = 300, bg = "white")
cat("wrote SUPP_FIG_XX_nil_hatching_full.{pdf,png}\n")

cat("\n== hatching, fraction (Wilson 95% CI) ==\n")
print(as.data.frame(phen %>%
  transmute(strain, condition = gsub("\\*", "", cond), embryos = n,
            hatched = round(p, 3), CI = sprintf("%.3f-%.3f", lo, hi)) %>%
  arrange(condition, desc(strain))), row.names = FALSE)
