## Diagnostic -- the panel split by genotype at two loci ---------------------
##
##   Rscript scripts/diagnostic_genotype_splits.R
##     -> plots/diagnostics/genotype_splits.{pdf,png}
##     -> plots/diagnostics/TABLE_genotype_splits.tsv
##
## THE QUESTION. Two loci, asked separately and then together:
##
##   IV:15,323,414   the pooled GWAS peak, -log10 p 8.84, the strongest marker
##                   in the scan and the only region that survives a stringent
##                   threshold
##   III:13,680,248  sid-2 T96K, the allele the cross and the NILs implicate,
##                   which has no marginal effect in this panel (p = 0.24)
##
## Splitting the panel at each and comparing the pooled pos-1 response says how
## much of the phenotype each allele accounts for on its own, and the 2x2 says
## whether they are independent -- i.e. whether more than one locus is
## segregating for this trait, which the chromosome IV signal already implies.
##
## DIRECTION. The vst trait is positive for strains that GAINED pool frequency
## under pos-1 RNAi, i.e. RESISTANT. Higher = more resistant. Checked, not
## assumed: it correlates +0.410 (p = 7.8e-6) with the ordinal plate resistance
## score in SUPP_FIG_plate_vs_paaby_vs_pos1original.R.
##
## Genotypes come from the CeNDR plink set under data/, so this is a diagnostic
## rather than a figure script -- a figure script must run from a clone.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(patchwork); library(ggtext)
})

OUT     <- "plots/diagnostics"
PLINK_D <- "data/genotypes/CeNDR20210121_Plink"
TRAITS  <- "data/pos1_original/updated_analysis/association_traits.csv"
TRAIT   <- "vst_ctrl_pos-1_T2"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

SITES <- tribble(
  ~key,     ~chrom, ~pos,      ~label,
  "chrIV",  "IV",   15323414L, "IV:15.32 Mb (GWAS peak)",
  "sid2",   "III",  13680248L, "sid-2 T96K")

## --- genotypes -------------------------------------------------------------
dose_at <- function(chrom, pos) {
  scratch <- file.path(tempdir(), paste0("split_", chrom, "_", pos))
  dir.create(scratch, showWarnings = FALSE, recursive = TRUE)
  tmp <- file.path(scratch, "site")
  writeLines(paste0(chrom, ":", pos), paste0(tmp, ".snp"))
  st <- system2("plink2", c("--bfile", file.path(PLINK_D, chrom),
                            "--extract", paste0(tmp, ".snp"),
                            "--export", "A", "--out", tmp, "--allow-extra-chr"),
                stdout = TRUE, stderr = TRUE)
  if (!file.exists(paste0(tmp, ".raw"))) {
    cat(st, sep = "\n"); stop("plink2 produced no .raw for ", chrom, ":", pos)
  }
  raw <- read_table(paste0(tmp, ".raw"), show_col_types = FALSE)
  col <- grep(paste0("^", chrom, ":", pos), names(raw), value = TRUE)
  stopifnot(length(col) == 1)
  counted <- sub(".*_", "", col)
  msg("  ", chrom, ":", pos, " counted allele ", counted)
  tibble(strain = raw$IID, dose = raw[[col]], counted = counted)
}

ph <- fread(TRAITS) %>% as_tibble() %>%
  transmute(strain, vst = .data[[TRAIT]]) %>% filter(!is.na(vst))
msg("phenotyped strains: ", nrow(ph))

gt <- SITES %>% mutate(g = map2(chrom, pos, dose_at)) %>%
  select(key, label, g) %>% unnest(g)

## Homozygous only. These are inbred isotypes, so a heterozygous call is a
## genotyping artefact rather than a biological state, and averaging it into
## either group would blur the very contrast being drawn.
d <- gt %>% inner_join(ph, by = "strain") %>%
  mutate(geno = case_when(dose == 2 ~ paste0(counted, counted),
                          dose == 0 ~ "alt/alt",
                          TRUE      ~ NA_character_))
het <- d %>% filter(is.na(geno))
if (nrow(het)) msg("dropping ", nrow(het), " non-homozygous call(s)")
d <- d %>% filter(!is.na(geno))

## --- one locus at a time ---------------------------------------------------
per_site <- d %>% group_by(key, label) %>%
  group_modify(function(x, ...) {
    g  <- sort(unique(x$geno))
    a  <- x$vst[x$geno == g[1]]; b <- x$vst[x$geno == g[2]]
    wt <- suppressWarnings(wilcox.test(a, b))
    tibble(geno_a = g[1], n_a = length(a), mean_a = mean(a), med_a = median(a),
           geno_b = g[2], n_b = length(b), mean_b = mean(b), med_b = median(b),
           diff_mean = mean(b) - mean(a),
           pct_of_range = 100 * abs(mean(b) - mean(a)) / diff(range(x$vst)),
           wilcox_p = wt$p.value)
  }) %>% ungroup()

cat("\n== the panel split at each locus (pooled pos-1 response, VST) ==\n")
print(as.data.frame(per_site %>%
  transmute(locus = label,
            `group A` = sprintf("%s (n=%d)", geno_a, n_a),
            `mean A` = sprintf("%+.4f", mean_a),
            `group B` = sprintf("%s (n=%d)", geno_b, n_b),
            `mean B` = sprintf("%+.4f", mean_b),
            `B - A` = sprintf("%+.4f", diff_mean),
            `% of range` = sprintf("%.1f", pct_of_range),
            `Wilcoxon p` = signif(wilcox_p, 3))), row.names = FALSE)

## --- the two together ------------------------------------------------------
w <- d %>% select(key, strain, geno, vst) %>%
  pivot_wider(names_from = key, values_from = geno) %>%
  filter(!is.na(chrIV), !is.na(sid2))
msg("strains with a homozygous call at BOTH loci: ", nrow(w))

joint <- w %>% group_by(chrIV, sid2) %>%
  summarise(n = n(), mean_vst = mean(vst), .groups = "drop")
cat("\n== jointly, by genotype at both loci ==\n")
print(as.data.frame(joint %>%
  mutate(mean_vst = sprintf("%+.4f", mean_vst))), row.names = FALSE)

## Are the two loci correlated across the panel? If they were, a marginal
## effect at one could be borrowed from the other.
ct <- table(w$chrIV, w$sid2)
cat("\n== 2x2 genotype table ==\n"); print(ct)
fi <- fisher.test(ct)
cat(sprintf("  Fisher p = %.3g  (odds ratio %.2f) -- tests whether the two\n",
            fi$p.value, unname(fi$estimate)))
cat("  loci are associated in the panel, not whether either affects the trait\n")

## additive two-locus fit, purely descriptive
fit <- lm(vst ~ chrIV + sid2, data = w)
cat("\n== additive model, vst ~ chrIV + sid-2 ==\n")
print(summary(fit)$coefficients)
cat(sprintf("  adjusted R2 = %.3f\n", summary(fit)$adj.r.squared))

## --- the marginal statistic, from the scan itself -------------------------
##
## The split above is a two-group comparison the panel supports. This is the
## marginal association the SCAN reports at the same marker, which is the
## number a reader will look for and the one that has to be quoted honestly:
## the single-marker mixed-model test, its rank among all 464,045 markers, and
## the variance a single-locus fit accounts for.
SCAN <- "supplemental_data/mapping/pos1_2023_gemma_loco.csv.gz"
gw <- fread(SCAN)[, .(chr, ps, af, beta, se, p_wald)]
gw[, lp := -log10(p_wald)]
setorder(gw, p_wald)
gw[, rank := .I]

marg <- SITES %>% rowwise() %>%
  mutate(m = list(gw[chr == chrom & ps == pos])) %>% ungroup() %>%
  mutate(found = map_int(m, nrow))
if (any(marg$found != 1))
  msg("NOTE: ", sum(marg$found != 1), " site(s) not in the scan's marker set")

r2_of <- function(k) {
  x <- d %>% filter(key == k)
  summary(lm(vst ~ geno, data = x))$r.squared
}

marg_tab <- marg %>% filter(found == 1) %>%
  mutate(af = map_dbl(m, "af"), beta = map_dbl(m, "beta"),
         se = map_dbl(m, "se"), p_wald = map_dbl(m, "p_wald"),
         lp = map_dbl(m, "lp"), rank = map_int(m, ~as.integer(.x$rank)),
         r2_single = map_dbl(key, r2_of)) %>%
  select(key, label, af, beta, se, p_wald, lp, rank, r2_single)

cat("\n== the MARGINAL statistic at each locus, from the scan ==\n")
print(as.data.frame(marg_tab %>%
  transmute(locus = label,
            af = sprintf("%.3f", af),
            beta = sprintf("%+.4f", beta),
            se = sprintf("%.4f", se),
            `p (Wald)` = signif(p_wald, 3),
            `-log10 p` = sprintf("%.2f", lp),
            `rank of 464,045` = format(rank, big.mark = ","),
            `single-locus R2` = sprintf("%.3f", r2_single))),
  row.names = FALSE)
cat("\n  The mixed model corrects for relatedness; the Wilcoxon split above\n")
cat("  does not, so the two are not required to agree and the gap between\n")
cat("  them is itself informative.\n")

write_tsv(per_site, file.path(OUT, "TABLE_genotype_splits.tsv"))
write_tsv(marg_tab, file.path(OUT, "TABLE_genotype_splits_marginal.tsv"))

## --- figure ----------------------------------------------------------------
pt <- function(l) paste0("<span style='font-size:13pt;color:#111111'>**", l, "**</span>")
th <- theme_classic(base_size = 11) +
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_markdown(size = 11.5),
        plot.subtitle = element_markdown(size = 8.6, colour = "grey30"),
        plot.title.position = "plot", legend.position = "none")

one <- function(k, letter) {
  x <- d %>% filter(key == k)
  s <- per_site %>% filter(key == k)
  ggplot(x, aes(geno, vst)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed",
               colour = "grey70") +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = "grey92",
                 colour = "grey35", linewidth = 0.35) +
    geom_point(position = position_jitter(width = 0.13, height = 0, seed = 1),
               size = 1.1, alpha = 0.55, colour = "#2C3E50") +
    labs(x = NULL, y = "*pos-1* response (VST)", title = pt(letter),
         subtitle = paste0(s$label, "<br>&Delta; mean ",
                           sprintf("%+.4f", s$diff_mean), ", Wilcoxon *p* = ",
                           signif(s$wilcox_p, 2))) +
    th + theme(axis.title.y = element_markdown())
}

pj <- ggplot(w %>% mutate(cell = paste0("IV ", chrIV, "\nsid-2 ", sid2)),
             aes(cell, vst)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed",
             colour = "grey70") +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = "grey92",
               colour = "grey35", linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.13, height = 0, seed = 2),
             size = 1.1, alpha = 0.55, colour = "#2C3E50") +
  labs(x = NULL, y = NULL, title = pt("C"),
       subtitle = "both loci, homozygous calls only") +
  th

fig <- (one("chrIV", "A") | one("sid2", "B") | pj) + plot_layout(widths = c(1, 1, 1.5))
ggsave(file.path(OUT, "genotype_splits.pdf"), fig, width = 10.5, height = 3.9,
       device = cairo_pdf)
ggsave(file.path(OUT, "genotype_splits.png"), fig, width = 10.5, height = 3.9,
       dpi = 300, bg = "white")
msg("wrote genotype_splits.{pdf,png}")
