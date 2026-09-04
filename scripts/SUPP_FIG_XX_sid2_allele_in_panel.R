## Supplement -- the T96K allele has no marginal effect in the mapping panel -
##
##   Rscript scripts/SUPP_FIG_XX_sid2_allele_in_panel.R
##     -> plots/SUPP_FIG_XX_sid2_allele_in_panel.{pdf,png}
##
##   A  the 2023 pooled pos-1 phenotype, split by sid-2 residue 96
##   B  the association scan around sid-2, showing the T96K marker itself
##
## THE POINT OF THIS FIGURE IS A NEGATIVE RESULT
## sid-2 T96K was identified by the cross and the NILs, not by the GWAS, and
## this is the check of whether it also acts as a marginal-effect variant
## across the wild panel. It does not. The variant is common (36% in the
## panel), it is genotyped in 230 of the 231 phenotyped strains, and it splits
## the phenotype distribution essentially not at all: Wilcoxon p = 0.99, and
## the marker's own association p is 0.24, ranking about 101,000th of 464,045.
##
## That is not evidence against the allele. It is what a background-dependent
## effect looks like from a marginal test. The editing experiments put the
## T96K contribution at roughly half of the JU1793-JU2466 difference, in a
## specific pair of backgrounds; a 36%-frequency variant with an effect that
## size and that context-dependent would not be expected to surface in a
## marginal scan of 231 strains, and it does not.
##
## Note the two tests disagree in an informative way: Wilcoxon p = 0.99 but
## Welch t p = 0.16. The means differ in the direction the crosses predict
## (96K lower, i.e. more sensitive) while the medians do not, because the 96T
## group carries the resistant tail. The apparent mean difference is a tail
## effect, not a shift.
##
## ALLELE ORIENTATION, CHECKED
## The site is chrIII:13,680,248 C>A, from
## supplemental_data/structure/sid2_parental_variants.tsv: REF C =
## 96T, ALT A = 96K. Genotypes come from the CeNDR 20210121 PLINK set. The
## orientation was verified against the four strains whose allele is known
## independently: JU1793 and N2 are 96T, JU2466 and XZ1516 are 96K. Both cross
## parents that behave as sensitive carry 96K.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(ggtext)
})

OUT   <- "plots"
DIR   <- "supplemental_data"
PLINK <- "supplemental_data/genotypes/sid2_region"
SITE  <- 13680248L
GENE  <- c(13679000, 13682000)   # sid-2, for the panel B window

COL_96T <- "#F34C00"   # the JU1793 allele
COL_96K <- "#40B4AB"   # the JU2466 allele
COL_PT  <- "#2E4057"
COL_THR <- "grey50"; COL_EIG <- "#1A7F5A"

KNOWN <- tribble(~strain, ~allele,
                 "JU1793", "96T", "N2", "96T",
                 "JU2466", "96K", "XZ1516", "96K")

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
wrap_md <- function(txt, width = 84)
  paste(strwrap(txt, width = width), collapse = "<br>")
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## ===========================================================================
## genotypes at the T96K site
## ===========================================================================
msg("genotypes at III:", SITE)
## plink2 is fussy about the doubled separator R's tempfile() can produce, so
## build the scratch path explicitly and report plink's own message on failure
scratch <- file.path(tempdir(), "sid2_t96k")
dir.create(scratch, showWarnings = FALSE, recursive = TRUE)
tmp <- file.path(scratch, "t96k")
writeLines(paste0("III:", SITE), paste0(tmp, ".snp"))
st <- system2("plink2", c("--bfile", PLINK, "--extract", paste0(tmp, ".snp"),
                          "--export", "A", "--out", tmp,
                          ## roman-numeral chromosome names are "extra" to
                          ## plink2 and rejected without this
                          "--allow-extra-chr"),
              stdout = TRUE, stderr = TRUE)
if (!file.exists(paste0(tmp, ".raw"))) {
  cat("\n-- plink2 output --\n"); cat(st, sep = "\n"); cat("\n")
  stop("plink2 did not produce ", tmp, ".raw")
}

raw <- read_table(paste0(tmp, ".raw"), show_col_types = FALSE)
dose_col <- grep(paste0("^III:", SITE), names(raw), value = TRUE)
stopifnot(length(dose_col) == 1)
## plink2 --export A names the column with the allele it counted; here that is
## the REF C, so 2 = homozygous C = 96T and 0 = homozygous A = 96K
counted <- sub(".*_", "", dose_col)
msg("  counted allele: ", counted)
stopifnot(counted == "C")

gt <- raw %>% transmute(strain = IID, dose = .data[[dose_col]]) %>%
  mutate(allele = case_when(dose == 2 ~ "96T", dose == 0 ~ "96K",
                            TRUE ~ NA_character_))

## orientation check against strains whose allele is known independently
chk <- KNOWN %>% left_join(gt, by = "strain", suffix = c("_expected", "_plink"))
cat("\n== allele orientation check ==\n")
print(as.data.frame(chk %>% transmute(strain, expected = allele_expected,
                                      from_plink = allele_plink)),
      row.names = FALSE)
stopifnot(all(chk$allele_expected == chk$allele_plink))
msg("  orientation confirmed on all ", nrow(chk), " known strains")

## ===========================================================================
## A -- phenotype split by allele
## ===========================================================================
tr <- read_csv(file.path(DIR, "phenotypes/pos1_2023_association_traits.csv.gz"), show_col_types = FALSE) %>%
  transmute(strain, vst = `vst_ctrl_pos-1_T2`) %>% filter(is.finite(vst))
d <- tr %>% inner_join(gt, by = "strain") %>% filter(!is.na(allele)) %>%
  mutate(allele = factor(allele, levels = c("96T", "96K")))
msg("  phenotyped strains with a genotype: ", nrow(d), " of ", nrow(tr))

sm <- d %>% group_by(allele) %>%
  summarise(n = n(), mean = mean(vst), median = median(vst), sd = sd(vst),
            .groups = "drop")
wt <- wilcox.test(vst ~ allele, data = d)
tt <- t.test(vst ~ allele, data = d)
r2 <- summary(lm(vst ~ allele, data = d))$r.squared
af <- sm$n[sm$allele == "96K"] / sum(sm$n)
msg("  96K frequency in the phenotyped set: ", round(af, 3))
msg("  Wilcoxon p = ", signif(wt$p.value, 3), " | Welch t p = ",
    signif(tt$p.value, 3), " | r2 = ", signif(r2, 3))

cat("\n== phenotype by sid-2 residue 96 ==\n")
print(as.data.frame(sm %>% mutate(across(mean:sd, ~round(.x, 4)))),
      row.names = FALSE)

lab <- sm %>% transmute(allele, y = -Inf,
                        txt = sprintf("%s<br>n = %d", allele, n))
mk <- d %>% inner_join(KNOWN %>% select(strain), by = "strain")

pA <- ggplot(d, aes(allele, vst, fill = allele)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_violin(width = 0.85, colour = "grey35", linewidth = 0.3,
              alpha = 0.35, trim = FALSE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, colour = "grey20",
               linewidth = 0.35, fill = "white") +
  geom_jitter(width = 0.13, height = 0, size = 1, alpha = 0.45,
              colour = COL_PT) +
  ## the four strains whose allele is independently known
  geom_point(data = mk, size = 2.4, shape = 21, fill = "white",
             colour = "grey15", stroke = 0.6) +
  ggrepel::geom_text_repel(data = mk, aes(label = strain), size = 2.7,
                           colour = "grey15", box.padding = 0.4,
                           min.segment.length = 0, seed = 1,
                           segment.size = 0.25) +
  scale_fill_manual(values = c(`96T` = COL_96T, `96K` = COL_96K),
                    guide = "none") +
  scale_x_discrete(labels = function(a)
    sprintf("%s\nn = %d", a, sm$n[match(a, sm$allele)])) +
  labs(x = "*sid-2* residue 96", y = "*pos-1* response (VST)",
       title = panel_title("A", "**The allele does not split the panel**"),
       subtitle = wrap_md(sprintf(paste0(
         "Every phenotyped strain in the 2023 panel (%d of %d genotyped), ",
         "split by residue 96. 96K is at %.0f%% frequency. The distributions ",
         "are not distinguishable: Wilcoxon *p* = %.2f. The means differ in ",
         "the direction the crosses predict (96K more sensitive) but only ",
         "because the 96T group holds the resistant tail &mdash; Welch *p* = ",
         "%.2f, *r*&sup2; = %.3f. Open circles are the four strains whose ",
         "allele is known independently."),
         nrow(d), nrow(tr), 100 * af, wt$p.value, tt$p.value, r2))) +
  theme_pub() +
  theme(axis.title.x = element_markdown(), axis.title.y = element_markdown())

## ===========================================================================
## B -- the scan around sid-2
## ===========================================================================
msg("panel B: local scan")
gw <- fread(file.path(DIR, "mapping/pos1_2023_gemma_loco.csv.gz"),
            select = c("chr", "ps", "af", "beta", "p_wald"))
gw <- gw[chr == "III"][, nlp := -log10(p_wald)]
source("scripts/gwas_thresholds.R", local = TRUE)
TH <- gwas_thresholds("pos1_2023")

WIN <- 3e5
loc <- gw[ps > SITE - WIN & ps < SITE + WIN]
site <- gw[ps == SITE]
rank_site <- sum(gw$nlp > site$nlp) + 1
msg("  T96K marker: -log10p ", round(site$nlp, 3), " | rank ", rank_site,
    " of ", nrow(gw), " on chrIII | local max ", round(max(loc$nlp), 2))

thr_lab <- tibble(x = min(loc$ps) / 1e6,
                  y = c(TH$bonferroni, TH$eigen),
                  label = c(sprintf("Bonferroni %.2f", TH$bonferroni),
                            sprintf("eigen %.2f", TH$eigen)),
                  col = c(COL_THR, COL_EIG))

pB <- ggplot(loc, aes(ps / 1e6, nlp)) +
  annotate("rect", xmin = GENE[1] / 1e6, xmax = GENE[2] / 1e6,
           ymin = -Inf, ymax = Inf, fill = "grey88", alpha = 0.6) +
  annotate("richtext", x = mean(GENE) / 1e6, y = Inf, label = "*sid-2*",
           vjust = 1.4, size = 3, colour = "grey30", fill = NA,
           label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt")) +
  geom_hline(yintercept = TH$eigen, linetype = "dotted", linewidth = 0.45,
             colour = COL_EIG) +
  geom_hline(yintercept = TH$bonferroni, linetype = "dashed",
             linewidth = 0.3, colour = COL_THR) +
  geom_point(size = 0.7, alpha = 0.5, colour = COL_PT) +
  geom_point(data = site, size = 2.6, colour = COL_96K) +
  ggrepel::geom_text_repel(data = site, aes(label = "T96K"), size = 3,
                           colour = COL_96K, nudge_y = 1.2,
                           min.segment.length = 0, seed = 2,
                           segment.size = 0.3) +
  ## a data frame, not annotate(): annotate() treats the length-4
  ## label.padding unit as an aesthetic and rejects vectorised y/label
  geom_richtext(data = thr_lab, aes(x = x, y = y, label = label),
                inherit.aes = FALSE, colour = thr_lab$col, hjust = 0,
                vjust = -0.35, size = 2.6, fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  scale_y_continuous(limits = c(0, max(TH$bonferroni, max(loc$nlp)) * 1.18),
                     expand = expansion(mult = c(0.02, 0))) +
  labs(x = "Chromosome III position (Mb)", y = "&minus;log<sub>10</sub>*p*",
       title = panel_title("B", "**Nor is it a signal in the scan**"),
       subtitle = wrap_md(sprintf(paste0(
         "GEMMA LOCO on the same trait, %.0f kb either side of the variant. ",
         "The T96K marker itself reaches &minus;log<sub>10</sub>*p* = %.2f ",
         "(*p* = %.2f), ranking %s of %s markers on chromosome III; nothing ",
         "within 300 kb clears either threshold. The chromosome III peak for ",
         "this trait is at 5.97 Mb, %.1f Mb away. sid-2 was found by the ",
         "cross and the NILs, not here."),
         WIN / 1e3, site$nlp, site$p_wald,
         format(rank_site, big.mark = ","), format(nrow(gw), big.mark = ","),
         abs(gw[order(-nlp)][1]$ps - SITE) / 1e6))) +
  theme_pub() +
  theme(axis.title.y = element_markdown())

## ===========================================================================
fig <- pA | pB
ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_allele_in_panel.pdf"), fig,
       width = 11.4, height = 5.0, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_allele_in_panel.png"), fig,
       width = 11.4, height = 5.0, dpi = 300, bg = "white")
msg("wrote SUPP_FIG_XX_sid2_allele_in_panel.{pdf,png}")
