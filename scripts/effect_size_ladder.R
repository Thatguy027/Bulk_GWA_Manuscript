## Effect sizes across the three chromosome III experiments -------------------
##
##   Rscript scripts/effect_size_ladder.R
##     -> plots/diagnostics/TABLE_effect_size_ladder.tsv
##
## Does the NIL interval's effect match the residue-96 allele swap, and either
## match the JU1793 x JU2466 chromosome III QTL?
##
## THE PROBLEM THAT DECIDES THE METHOD. The two hatching experiments do not
## share a scale. JU2466 hatches 0.357 under pos-1 in the NIL series and 0.045
## in the allele-swap series -- the same strain, the same nominal 50% dose, an
## eightfold difference in surviving fraction. Raw percentage-point effects are
## therefore not comparable between the two, and any statement that they are
## would be an artefact of assay stringency. Every effect below is reported
## twice: as a raw difference in hatched fraction, and normalised to the
## JU1793-JU2466 span measured WITHIN THAT SAME EXPERIMENT, where 0 is JU2466
## and 1 is JU1793.
##
## THE QTL IS NOT ON EITHER SCALE. An allele-frequency shift in a selected pool
## is a selection response; converting it to a hatching difference needs the
## selection intensity and the number of generations, and METHODS.txt carries
## those as [TO FILL]. So the QTL is compared on direction and rank only, which
## is what the data support, and the arithmetic is not faked.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse)})

HA <- "supplemental_data/hatching_assays"
OUT <- "plots/diagnostics"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## Wilson interval, and a Newcombe interval for a difference of proportions
wilson <- function(x, n, z = qnorm(0.975)) {
  p <- x/n
  c(lo = max(0, (p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/(1 + z^2/n)),
    hi = min(1, (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/(1 + z^2/n)))
}
## Newcombe's hybrid-score interval for a difference of two proportions.
## unname() matters: wilson() returns named elements, and c(lo = <named "lo">)
## would produce "lo.lo" and silently break lookup by name.
diff_ci <- function(x1, n1, x2, n2) {
  a <- unname(wilson(x1, n1)); b <- unname(wilson(x2, n2))
  p1 <- x1/n1; p2 <- x2/n2; d <- p1 - p2
  c(diff = d,
    lo = d - sqrt((p1 - a[1])^2 + (b[2] - p2)^2),
    hi = d + sqrt((a[2] - p1)^2 + (p2 - b[1])^2))
}

## ---- the NIL series: fraction HATCHED under pos-1 -------------------------
nil <- read_tsv(file.path(HA, "nil_series_hatching.tsv"), show_col_types = FALSE) %>%
  filter(condition == "pos1") %>%
  transmute(experiment = "NIL series", strain = Strain,
            n = `plated embryo`, hatched = `plated embryo` - unhatched,
            p = hatched/n)
bed <- read_tsv(file.path(HA, "nil_introgression_ranges.bed"), col_names =
                  c("chrom","start","end","strain","introgression"), show_col_types = FALSE)

## ---- the N2 allele swaps, a dose series -----------------------------------
## N2 carries 96T and two independently edited lines carry 96K, across
## 0/25/50/75/100% pos-1 food. There is no second parental genotype here, so
## this experiment has no JU1793-JU2466-style span and its effects are reported
## raw only -- the normalised column is left empty rather than filled with a
## quantity that means something different from the one above it.
n2 <- read_tsv(file.path(HA, "n2_allele_swaps_hatching.tsv"), show_col_types = FALSE) %>%
  transmute(experiment = "N2 allele swap", strain, genotype, dose = condition,
            n = n_plated, hatched = n_plated - n_unhatched, p = hatched/n)

## ---- the JU allele swaps --------------------------------------------------
swap <- read_csv(file.path(HA, "ju_allele_swaps_hatching.csv"), show_col_types = FALSE) %>%
  filter(condition == "pos") %>%
  transmute(experiment = "JU allele swap", strain, genotype,
            n = n_plated, hatched = n_plated - n_unhatched, p = hatched/n)

## ---- one normalised scale per experiment ----------------------------------
## anchors: JU1793 = 1, JU2466 = 0, measured inside the same experiment.
## JU2466 appears as two independent lines in the swap series; pool them.
anch <- function(d, hi_strain, lo_strains) {
  hi <- d %>% filter(strain %in% hi_strain) %>% summarise(p = sum(hatched)/sum(n)) %>% pull(p)
  lo <- d %>% filter(strain %in% lo_strains) %>% summarise(p = sum(hatched)/sum(n)) %>% pull(p)
  c(lo = lo, hi = hi, span = hi - lo)
}
A_nil  <- anch(nil,  "JU1793", "JU2466")
A_swap <- anch(swap, "JU1793", c("JU2466_A","JU2466_B"))
cat(sprintf("NIL series      JU2466 %.3f -> JU1793 %.3f   span %.3f\n", A_nil["lo"],  A_nil["hi"],  A_nil["span"]))
cat(sprintf("JU allele swap  JU2466 %.3f -> JU1793 %.3f   span %.3f\n", A_swap["lo"], A_swap["hi"], A_swap["span"]))
cat(sprintf("\nSAME STRAINS, DIFFERENT ASSAYS: JU2466 %.3f vs %.3f (%.1f-fold), JU1793 %.3f vs %.3f\n",
            A_nil["lo"], A_swap["lo"], A_nil["lo"]/A_swap["lo"], A_nil["hi"], A_swap["hi"]))

norm <- function(p, A) (p - A["lo"]) / A["span"]
nil  <- nil  %>% mutate(rel = norm(p, A_nil))
swap <- swap %>% mutate(rel = norm(p, A_swap))

cat("\n== NIL series, ordered by hatching ==\n")
print(as.data.frame(nil %>% arrange(desc(p)) %>%
  left_join(bed %>% select(strain, start, end), by = "strain") %>%
  transmute(strain, n, hatched = round(p, 3), rel = round(rel, 3),
            introgression = ifelse(is.na(start), "-", sprintf("%.3f-%.3f Mb", start/1e6, end/1e6)))),
  row.names = FALSE)
cat("\n== JU allele swaps ==\n")
print(as.data.frame(swap %>% arrange(desc(p)) %>%
  transmute(genotype, n, hatched = round(p, 3), rel = round(rel, 3))), row.names = FALSE)

## ---- the contrasts that answer the question ------------------------------
get <- function(d, s) d %>% filter(strain %in% s) %>%
  summarise(x = sum(hatched), n = sum(n)) %>% as.list()
mk <- function(label, exper, a, b, A, note) {
  d <- diff_ci(a$x, a$n, b$x, b$n)
  span <- if (is.null(A)) NA_real_ else A[["span"]]
  tibble(experiment = exper, contrast = label,
         raw.diff = d[["diff"]], raw.lo = d[["lo"]], raw.hi = d[["hi"]],
         frac.of.span = d[["diff"]]/span,
         span.lo = d[["lo"]]/span, span.hi = d[["hi"]]/span, note = note)
}
## the N2 series: 96T minus 96K at each dose, the two 96K lines pooled
gn2 <- function(gt, dose) n2 %>% filter(genotype == gt, dose == !!dose) %>%
  summarise(x = sum(hatched), n = sum(n)) %>% as.list()
res <- bind_rows(
  mk("JU1793 -> wSZ191 (JU2466 alleles at the 37 kb interval only)", "NIL series",
     get(nil,"JU1793"), get(nil,"wSZ191"), A_nil,
     "the interval's own effect, on an otherwise JU1793 genome"),
  mk("JU1793 -> wSZ196 (JU2466 alleles distal to the interval only)", "NIL series",
     get(nil,"JU1793"), get(nil,"wSZ196"), A_nil,
     "the distal region alone, as the control on the contrast above"),
  mk("JU1793 -> wSZ176 (JU2466 alleles across interval + distal)", "NIL series",
     get(nil,"JU1793"), get(nil,"wSZ176"), A_nil,
     "the union of the two segments above"),
  mk("JU1793 -> JU2466 (whole genome)", "NIL series",
     get(nil,"JU1793"), get(nil,"JU2466"), A_nil, "the full parental difference"),
  ## Both swap rows are stated the same way the NIL rows are: the COST of
  ## carrying the JU2466 allele, so 96T minus 96K in each background. Computing
  ## one as wild-minus-swapped and the other as swapped-minus-wild would put
  ## opposite signs on the same biological statement and make the table read as
  ## a contradiction of the prose.
  mk("96T vs 96K in the JU1793 background (cost of the JU2466 residue)", "JU allele swap",
     get(swap,"JU1793"), get(swap,"wSZ200"), A_swap,
     "residue 96 alone, JU1793 background"),
  mk("96T vs 96K in the JU2466 background (cost of the JU2466 residue)", "JU allele swap",
     get(swap,"wSZ206"), get(swap,c("JU2466_A","JU2466_B")), A_swap,
     "residue 96 alone, JU2466 background"),
  mk("JU1793 -> JU2466 (whole genome)", "JU allele swap",
     get(swap,"JU1793"), get(swap,c("JU2466_A","JU2466_B")), A_swap,
     "the full parental difference"),
  mk("96T vs 96K in the N2 background, 25% pos-1 food", "N2 allele swap",
     gn2("N2[96T]", 25), gn2("N2[96K]", 25), NULL,
     "two independently edited 96K lines pooled; no second parent, so no span"),
  mk("96T vs 96K in the N2 background, 50% pos-1 food", "N2 allele swap",
     gn2("N2[96T]", 50), gn2("N2[96K]", 50), NULL,
     "N2 is already near-fully sensitive at this dose, so the window has closed"),
  mk("96T vs 96K in the N2 background, no pos-1 food", "N2 allele swap",
     gn2("N2[96T]", 0), gn2("N2[96K]", 0), NULL,
     "the negative control: no RNAi, so no effect expected"))
res <- res %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
write_tsv(res, file.path(OUT, "TABLE_effect_size_ladder.tsv"))

cat("\n== effect sizes, raw and as a fraction of that experiment's parental span ==\n")
print(as.data.frame(res %>% transmute(experiment, contrast,
        raw = sprintf("%+.3f [%+.3f, %+.3f]", raw.diff, raw.lo, raw.hi),
        `of span` = ifelse(is.na(frac.of.span), "-",
                           sprintf("%.0f%% [%.0f, %.0f]",
                                   100*frac.of.span, 100*span.lo, 100*span.hi)))),
      row.names = FALSE)

cat("\n== the N2 swap across the whole dose series ==\n")
print(as.data.frame(n2 %>% group_by(genotype, dose) %>%
  summarise(n = sum(n), hatched = round(sum(hatched)/sum(n), 3), .groups = "drop") %>%
  pivot_wider(names_from = genotype, values_from = c(n, hatched)) %>%
  mutate(effect = round(`hatched_N2[96T]` - `hatched_N2[96K]`, 3)) %>%
  arrange(dose)), row.names = FALSE)

## ---- additivity check on the NIL series ----------------------------------
pI <- nil$p[nil$strain=="wSZ191"]; pD <- nil$p[nil$strain=="wSZ196"]
pB <- nil$p[nil$strain=="wSZ176"]; p0 <- nil$p[nil$strain=="JU1793"]
cat(sprintf("\n== additivity of the two chromosome III segments ==\n"))
cat(sprintf("  interval alone      %.3f  (cost %.3f)\n", pI, p0-pI))
cat(sprintf("  distal alone        %.3f  (cost %.3f)\n", pD, p0-pD))
cat(sprintf("  both, additive      %.3f  predicted\n", p0-(p0-pI)-(p0-pD)))
cat(sprintf("  both, observed      %.3f  (wSZ176)\n", pB))
cat(sprintf("  excess over additive %+.3f\n", pB-(p0-(p0-pI)-(p0-pD))))

## ---- the QTL: direction and rank only ------------------------------------
cqf <- read_tsv("plots/TABLE_cross_qtl_full.tsv", show_col_types = FALSE) %>%
  filter(cross == "JU1793xJU2466", contrast == "ht115 vs pos1", separated)
cat("\n== JU1793 x JU2466, HT115 vs pos-1: every separated QTL by |dfreq| ==\n")
print(as.data.frame(cqf %>% arrange(desc(abs(dfreq))) %>%
  transmute(chrom, peak.Mb = round(peak.Mb,2), LOD = round(peak.LOD),
            dfreq = round(dfreq,3), rank = peak.rank) %>% head(8)), row.names = FALSE)
cat("\nSign convention: dfreq is parent 1 (JU1793) frequency in pool a minus pool b,\n")
cat("for 'ht115 vs pos1' that is HT115 minus pos-1. A negative value means the\n")
cat("JU1793 allele is ENRICHED in the pos-1-selected pool, i.e. JU1793 confers\n")
cat("resistance -- the same direction as both hatching experiments.\n")
