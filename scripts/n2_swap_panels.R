## N2 sid-2 residue-96 swap: shared data prep and panel builders -------------
## Sourced by scripts/SUPP_FIG_XX_n2_swap_dose.R and by scripts/Figure4_sid2.R.
## Not run on its own.
##
## THE EXPERIMENT
## N2 carries 96T, the same allele as JU1793. Two independent lines, wSZ203
## and wSZ204, were edited to 96K, the JU2466 allele. All three were assayed
## across a dilution series of pos-1 RNAi bacteria: 0, 25, 50, 75 and 100%.
##
## WHY THE DOSE SERIES MATTERS
## The effect is only visible at 25%. At 0% everything hatches and at 50% and
## above everything is dead in every genotype, so those doses are at the
## ceiling and the floor and can carry no difference:
##
##   % pos-1   96T (N2)   96K (pooled)   difference   Fisher p
##        0      0.985        0.992        -0.006      0.43
##       25      0.323        0.040        +0.282      8e-25
##       50      0.060        0.034        +0.026      0.11
##       75      0.034        0.027        +0.007      0.64
##      100      0.000        0.007        -0.007      0.56
##
## So "the biggest difference is at 25%" is not a matter of degree: it is the
## only dose at which the assay has any dynamic range in this background. A
## single-dose experiment at 100% would have found nothing. That is worth
## saying in the text, because it is also the likely reason a sub-maximal dose
## is needed to see sid-2 alleles at all in a resistant background.
##
## DIRECTION, AND WHY IT MATTERS
## 96K makes N2 MORE sensitive, which is the same direction as JU2466 (96K,
## sensitive) against JU1793 (96T, resistant). This is the third genetic
## background to show it and the only one with two independently derived
## edited lines, which agree closely (0.044 and 0.037 at 25%).
##
## CAVEAT: one plate per strain per dose. Intervals are Wilson binomial
## intervals on that plate's embryo count, so they describe counting
## uncertainty, not between-plate variability. wSZ203 and wSZ204 are
## independent lines of the same edit, so they are the closest thing here to a
## biological replicate.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggtext)
})

N2SWAP <- "data/plate_rnai_phenotyping/n2_swap.tsv"

## N2 has no strain colour in this manuscript; parental keeps the dark
## reference colour and edited lines go grey, as in Figure 4B
COL_N2   <- "#2E4057"
COL_EDIT <- "grey70"
FOCAL_DOSE <- 25L

n2_swap_data <- function() {
  z <- qnorm(0.975)
  read_tsv(N2SWAP, show_col_types = FALSE) %>%
    mutate(dose = as.integer(condition),
           hatched = n_plated - n_unhatched,
           p = hatched / n_plated,
           lo = pmax(0, (p + z^2 / (2 * n_plated) -
                  z * sqrt(p * (1 - p) / n_plated + z^2 / (4 * n_plated^2))) /
                  (1 + z^2 / n_plated)),
           hi = pmin(1, (p + z^2 / (2 * n_plated) +
                  z * sqrt(p * (1 - p) / n_plated + z^2 / (4 * n_plated^2))) /
                  (1 + z^2 / n_plated)),
           allele = ifelse(genotype == "N2[96T]", "96T", "96K"),
           allele = factor(allele, levels = c("96T", "96K")),
           line = factor(paste0(genotype, "  ", strain),
                         levels = c("N2[96T]  N2", "N2[96K]  wSZ203",
                                    "N2[96K]  wSZ204")))
}

## Fisher's exact test of each edited line against N2 at the same dose. The
## reference is JOINED, not filtered inside a rowwise mutate -- `dose == dose`
## there compares the column with itself and silently builds a 5-row table.
n2_swap_tests <- function(d) {
  ref <- d %>% filter(strain == "N2") %>%
    select(dose, r_h = hatched, r_n = n_plated)
  d %>% filter(strain != "N2") %>% left_join(ref, by = "dose") %>%
    rowwise() %>%
    mutate(ft = list(fisher.test(matrix(c(hatched, n_plated - hatched,
                                          r_h, r_n - r_h),
                                        nrow = 2, byrow = TRUE))),
           fisher_p = ft$p.value) %>%
    ungroup() %>% select(-ft)
}

fmt_p <- function(p) ifelse(p < 0.001, "*p* < 0.001", sprintf("*p* = %.3f", p))
