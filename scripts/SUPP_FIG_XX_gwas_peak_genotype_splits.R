## Supplement -- the panel split at each pos-1 pilot peak marker -------------
##
##   Rscript scripts/SUPP_FIG_XX_gwas_peak_genotype_splits.R
##     -> plots/SUPP_FIG_XX_gwas_peak_genotype_splits.{pdf,png}
##
## WHAT THIS SHOWS
## Figure 1C is a scan: it reports where the association is, not how large the
## underlying difference is or how many strains carry it. This is the same four
## peak markers drawn as the thing being tested -- the pooled pos-1 response of
## the strains carrying each homozygous genotype.
##
##   A  IV:15,323,414   the scan's strongest marker, -log10 p 8.84
##   B  X:4,875,969     the single chromosome X marker over Bonferroni
##   C  III:5,965,738   the centre-of-III marker, over Bonferroni but with no
##                      eigen-passing neighbour within 100 kb
##   D  III:12,718,465  the right-arm locus, under Bonferroni and over the
##                      eigen threshold, with 14 supporting markers
##   E  r-squared between all six pairs of those markers, lower triangle
##
## WHY C AND D ARE BOTH HERE
## They are the two chromosome III signals, and the repository admits the
## weaker one. That is the interval-admission rule in FIGURE_REPORT.md, not an
## oversight: a lone marker with no support around it is what genotyping error
## and unshared haplotypes produce, and 14 markers agreeing over 100 kb is what
## a locus produces. This figure lets a reader see what the rule discards.
##
## THE OBSERVATION THIS FIGURE MAKES
## For three of the four markers, twice the mixed-model effect size -- the
## difference the model predicts between the two homozygotes -- lands within a
## few percent of the difference actually observed. For III:5,965,738 it is
## 2.7 times larger than the observed split. An effect that only appears once
## the kinship correction is applied is the signature of a marker tracking
## relatedness rather than a locus, which is the same conclusion the admission
## rule reaches from the marker's isolation. The two arguments are independent.
##
## PANEL E: LD, AND WHY IT IS COMPUTED HERE RATHER THAN WITH PLINK
## r-squared is the squared Pearson correlation of the dosages. For these
## strains that is not an approximation to plink's --r2: there are no
## heterozygous calls at these four markers, so the correlation-based and
## haplotypic estimators coincide exactly. Verified against
## `plink --r2 inter-chr` on data/genotypes/CeNDR20210121_Plink, restricted to
## the same 231 strains -- agreement to 3e-08, which is plink's own printed
## precision. Those six values are pinned in R2_PLINK below and asserted, so
## the deposit cannot drift away from the panel it was checked against.
##
## Computing it here means the figure needs neither the archive nor a plink
## binary. LD is computed IN THE PHENOTYPED PANEL, not in all 540 isotypes,
## which is the convention scripts/gwas_qtl_intervals.R sets and states: LD is
## a property of the sample being analysed.
##
## DIRECTION, CHECKED NOT ASSUMED
## The vst trait is positive for strains that GAINED pool frequency under pos-1
## RNAi, i.e. RESISTANT: higher = more resistant. It correlates +0.410
## (p = 7.8e-6) with the ordinal plate resistance score in
## SUPP_FIG_plate_vs_paaby_vs_pos1original.R.
##
## ALLELE ORIENTATION, CHECKED NOT ASSUMED
## GEMMA's reported af matches the frequency of the allele plink2 did NOT count
## at all four markers, i.e. GEMMA counted the minor allele, and this script
## asserts that rather than trusting it. Every beta is positive, so at all four
## markers the MINOR allele is the resistant one.
##
## GENOTYPES
## supplemental_data/genotypes/gwas_peak_genotypes.tsv, staged by
## make_supplemental_data.R from the CeNDR 20210121 PLINK set. A long dosage
## table rather than a PLINK subset, because the four markers sit on three
## chromosomes and because a TSV needs no plink2 binary to read -- so this
## figure builds from a clone with nothing installed.
##
## HETEROZYGOUS CALLS
## These are inbred isotypes, so a heterozygous call is a genotyping artefact
## rather than a biological state. There are none at these four markers; the
## filter is kept, and reports itself if that ever changes.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(ggtext)
})

OUT   <- "plots"
DIR   <- "supplemental_data"
GENO  <- file.path(DIR, "genotypes/gwas_peak_genotypes.tsv")
PHENO <- file.path(DIR, "phenotypes/pos1_2023_association_traits.csv.gz")
SCAN  <- file.path(DIR, "mapping/pos1_2023_gemma_loco.csv.gz")
EIGEN <- file.path(DIR, "mapping/eigen_independent_tests.tsv")
TRAIT <- "vst_ctrl_pos-1_T2"

SUPPORT_WIN <- 1e5     # the admission rule's neighbourhood, +/- 100 kb

COL_MAJ <- "#B8C2CA"   # the common allele
COL_MIN <- "#2E4057"   # the minor allele, which is the resistant one here
COL_PT  <- "#2E4057"
COL_OK  <- "#1A7F5A"   # admitted by the interval rule
COL_NO  <- "#B03A2E"   # discarded by it
COL_LD_LO <- "#EEF2F5"
COL_LD_HI <- "#2E4057"

## the interval work's LD criterion, for scale: every pair here is far below it
LD_INTERVAL_R2 <- 0.7

## plink --r2 inter-chr on the archive, same 231 strains. Pinned so that a
## change in the deposited dosages cannot silently move panel E.
R2_PLINK <- tribble(
  ~a,             ~b,             ~r2,
  "X:4875969",    "III:5965738",  0.00400474,
  "X:4875969",    "III:12718465", 0.0320088,
  "X:4875969",    "IV:15323414",  0.0671796,
  "III:5965738",  "III:12718465", 0.0000947867,
  "III:5965738",  "IV:15323414",  0.0026374,
  "III:12718465", "IV:15323414",  0.112313)

## the four markers, in the order the panels run
SITES <- tribble(
  ~key,      ~chrom, ~pos,      ~short,         ~note,
  "IV",      "IV",   15323414L, "IV:15.32 Mb",  "The scan's strongest marker, on the right arm.",
  "X",       "X",     4875969L, "X:4.88 Mb",    "The only chromosome X marker over Bonferroni, 28% along the chromosome.",
  "III_mid", "III",   5965738L, "III:5.97 Mb",  "The centre of the chromosome.",
  "III_arm", "III",  12718465L, "III:12.72 Mb", "The right arm, 375 kb proximal to sid-2.")

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
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
## strwrap will break an HTML tag across lines and gridtext then renders the
## raw markup, so ONLY plain text goes through here. Subscripts and Greek are
## unicode for the same reason; markup that has to survive (the coloured
## admission status) goes in the title, which is not wrapped.
wrap_md <- function(txt, width = 72)
  paste(strwrap(txt, width = width), collapse = "<br>")
## emphasis is added AFTER wrapping, for the same reason: strwrap must never
## see the asterisks
emph <- function(txt) {
  txt <- gsub("\\bp\\b", "*p*", txt)
  gsub("sid-2", "*sid-2*", txt)
}

## ===========================================================================
## inputs
## ===========================================================================
stopifnot(file.exists(GENO), file.exists(PHENO), file.exists(SCAN),
          file.exists(EIGEN))

ph <- read_csv(PHENO, show_col_types = FALSE) %>%
  transmute(strain, vst = .data[[TRAIT]]) %>%
  filter(!is.na(vst))
msg("phenotyped strains: ", nrow(ph))
stopifnot(nrow(ph) == 231)

gt <- read_tsv(GENO, show_col_types = FALSE)
stopifnot(setequal(gt$marker, paste0(SITES$chrom, ":", SITES$pos)))

scan_all <- fread(SCAN, showProgress = FALSE) %>% as_tibble() %>%
  mutate(marker = paste0(chr, ":", ps), lp = -log10(p_wald))
n_marker <- nrow(scan_all)

eig <- read_tsv(EIGEN, show_col_types = FALSE) %>%
  filter(panel == "pos1_2023", chrom == "all")
stopifnot(nrow(eig) == 1, eig$n_marker == n_marker)
THR_BONF  <- eig$thr_bonferroni
THR_EIGEN <- eig$thr_eigen_liji
msg("thresholds: Bonferroni ", sprintf("%.3f", THR_BONF),
    " | eigen (Li & Ji, M_eff ", eig$M_eff_liji, ") ",
    sprintf("%.3f", THR_EIGEN))

## the scan's own numbers for these four markers
st <- SITES %>%
  mutate(marker = paste0(chrom, ":", pos)) %>%
  left_join(scan_all %>% select(marker, af, beta, p_wald, lp), by = "marker")
stopifnot(!anyNA(st$lp))

## how many eigen-passing markers sit within the admission window
st <- st %>%
  mutate(support = map2_int(chrom, pos, function(c, p)
    sum(scan_all$chr == c & abs(scan_all$ps - p) <= SUPPORT_WIN &
        scan_all$lp >= THR_EIGEN & scan_all$ps != p)),
    over_bonf = lp >= THR_BONF,
    admitted = support > 0 & (over_bonf | lp >= THR_EIGEN))

## ===========================================================================
## genotypes, and the orientation check
## ===========================================================================
d <- gt %>% inner_join(ph, by = "strain")

het <- d %>% filter(dose == 1)
if (nrow(het)) msg("dropping ", nrow(het), " heterozygous call(s)") else
  msg("no heterozygous calls at these markers, as expected for isotypes")
d <- d %>% filter(!is.na(dose), dose != 1)

## GEMMA's af should equal the frequency of ONE of the two alleles. Establish
## which, per marker, rather than assuming: the effect allele is the one whose
## frequency matches, and beta refers to it.
orient <- d %>%
  group_by(marker, allele_counted, allele_other) %>%
  summarise(n_called = n(), freq_counted = mean(dose) / 2, .groups = "drop") %>%
  left_join(st %>% select(marker, af, beta, lp), by = "marker") %>%
  mutate(effect_allele = case_when(
           abs(freq_counted - af)       < 0.005 ~ allele_counted,
           abs((1 - freq_counted) - af) < 0.005 ~ allele_other,
           TRUE ~ NA_character_))
if (anyNA(orient$effect_allele))
  stop("GEMMA's af matches neither allele's frequency at: ",
       paste(orient$marker[is.na(orient$effect_allele)], collapse = ", "),
       " -- the deposited genotypes and the scan disagree, so no panel here ",
       "can be labelled")
msg("allele orientation reproduces GEMMA's af at all ", nrow(orient),
    " markers; effect allele is the one plink2 did not count in ",
    sum(orient$effect_allele == orient$allele_other), " of ", nrow(orient))

st <- st %>% left_join(orient %>% select(marker, allele_counted, allele_other,
                                         effect_allele, n_called),
                       by = "marker")

## label each strain by its homozygous genotype, and mark which group carries
## the effect allele
d <- d %>%
  left_join(st %>% select(marker, effect_allele), by = "marker") %>%
  mutate(hom = if_else(dose == 2, allele_counted, allele_other),
         geno = paste0(hom, "/", hom),
         is_effect = hom == effect_allele)

## ===========================================================================
## the split at each marker
## ===========================================================================
stats <- d %>%
  group_by(marker) %>%
  group_modify(function(x, ...) {
    e <- x$vst[x$is_effect]; o <- x$vst[!x$is_effect]
    wt <- suppressWarnings(wilcox.test(e, o))
    tibble(n_eff = length(e), n_oth = length(o),
           mean_eff = mean(e), mean_oth = mean(o),
           med_eff = median(e), med_oth = median(o),
           obs_diff = mean(e) - mean(o),
           range_all = diff(range(x$vst)),
           wilcox_p = wt$p.value)
  }) %>% ungroup() %>%
  left_join(st, by = "marker") %>%
  mutate(pred_diff = 2 * beta,                 # the model's two-homozygote gap
         pct_of_range = 100 * abs(obs_diff) / range_all,
         infl = pred_diff / obs_diff)

cat("\n== the panel split at each peak marker (pooled pos-1 response, VST) ==\n")
print(as.data.frame(stats %>% arrange(match(marker, st$marker)) %>%
  transmute(marker,
            `effect allele` = effect_allele,
            af = sprintf("%.3f", af),
            `n eff` = n_eff, `n other` = n_oth,
            `mean eff` = sprintf("%+.4f", mean_eff),
            `mean other` = sprintf("%+.4f", mean_oth),
            observed = sprintf("%+.4f", obs_diff),
            `2*beta` = sprintf("%+.4f", pred_diff),
            `ratio` = sprintf("%.2f", infl),
            `% range` = sprintf("%.1f", pct_of_range),
            `Wilcoxon p` = signif(wilcox_p, 3),
            `-log10 p` = sprintf("%.2f", lp),
            support = support,
            admitted = admitted)), row.names = FALSE)

## The observation the header describes, asserted so the figure cannot quietly
## stop making it: three markers agree with the model, the centre-of-III one
## does not.
agree <- stats %>% filter(marker != "III:5965738")
stopifnot(all(abs(agree$infl - 1) < 0.15))
odd <- stats %>% filter(marker == "III:5965738")
stopifnot(odd$infl > 2)
msg("model-vs-observed: ", nrow(agree), " markers within 15% of 1.0; ",
    "III:5,965,738 inflated ", sprintf("%.1f", odd$infl), "x")

## ===========================================================================
## panels
## ===========================================================================
YL <- range(d$vst) + c(-0.04, 0.06) * diff(range(d$vst))

build <- function(i) {
  s  <- st[i, ]
  ss <- stats %>% filter(marker == s$marker)
  x  <- d %>% filter(marker == s$marker) %>%
    ## the common allele first, the effect allele second
    mutate(geno = fct_reorder(geno, as.integer(is_effect)))

  lab <- x %>% count(geno, is_effect) %>%
    mutate(txt = paste0(geno, "<br><span style='font-size:6.6pt'>n = ", n,
                        "</span>"))
  ## the coloured status lives in the TITLE, which is not wrapped, so its
  ## markup survives
  status <- sprintf(
    "<span style='font-size:8.2pt;color:%s'>%s</span>",
    if (s$admitted) COL_OK else COL_NO,
    if (s$admitted) "admitted" else "discarded")
  thr <- if (s$over_bonf) "over Bonferroni" else
    "under Bonferroni but over the eigen threshold"
  pfmt <- function(p) if (p < 1e-3) sprintf("%.0e", p) else sprintf("%.3f", p)

  sub <- emph(paste(
    wrap_md(s$note, 46),
    wrap_md(sprintf("\u2212log\u2081\u2080 p %.2f, %s, with %d eigen-passing neighbour%s within 100 kb.",
                    s$lp, thr, s$support, if (s$support == 1) "" else "s"), 46),
    wrap_md(sprintf("Effect allele %s, frequency %.3f. Wilcoxon p = %s.",
                    s$effect_allele, s$af, pfmt(ss$wilcox_p)), 46),
    wrap_md(sprintf("Observed gap %+.4f; the model's 2\u03b2 is %+.4f.",
                    ss$obs_diff, ss$pred_diff), 46),
    sep = "<br>"))

  ggplot(x, aes(geno, vst)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
    geom_boxplot(aes(fill = is_effect), width = 0.52, outlier.shape = NA,
                 linewidth = 0.32, colour = "grey35", alpha = 0.85) +
    ## seeded, or the figure does not reproduce byte-for-byte -- the same
    ## convention as every other jittered panel in scripts/
    geom_point(position = position_jitter(width = 0.14, height = 0, seed = 1),
               size = 0.85, alpha = 0.6, colour = COL_PT) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 1.9,
                 fill = "white", colour = "grey20", stroke = 0.5) +
    scale_fill_manual(values = c(`FALSE` = COL_MAJ, `TRUE` = COL_MIN),
                      guide = "none") +
    scale_x_discrete(labels = setNames(lab$txt, lab$geno)) +
    coord_cartesian(ylim = YL) +
    labs(x = NULL,
         y = if (i == 1) "Pooled *pos-1* response (VST)" else NULL,
         title = panel_title(LETTERS[i],
                             paste0("**", s$short, "** &nbsp;", status)),
         subtitle = sub) +
    theme_pub(10) +
    theme(axis.text.x = element_markdown(size = 8.4, lineheight = 1.25),
          axis.title.y = element_markdown(size = 9.2),
          plot.subtitle = element_markdown(size = 6.9, colour = "grey30",
                                           lineheight = 1.3),
          plot.title = element_markdown(size = 10.5))
}

## ===========================================================================
## E -- LD between the four markers
## ===========================================================================
msg("panel E: r-squared between the four markers")

## one row per strain, one column per marker, in panel order
W <- d %>% select(marker, strain, dose) %>%
  pivot_wider(names_from = marker, values_from = dose) %>%
  select(strain, all_of(st$marker))
M <- as.matrix(W[, -1])
R2 <- cor(M, use = "pairwise.complete.obs")^2

## the pinned plink values, asserted in whichever order they were recorded
chk <- R2_PLINK %>%
  mutate(mine = map2_dbl(a, b, ~ R2[.x, .y]), diff = abs(mine - r2))
cat("\n== r-squared: this script against plink --r2 inter-chr ==\n")
print(as.data.frame(chk %>% transmute(
  pair = paste(a, b, sep = " / "),
  plink = signif(r2, 6), computed = signif(mine, 6),
  difference = signif(diff, 3))), row.names = FALSE)
if (max(chk$diff) > 1e-6)
  stop("r-squared disagrees with the pinned plink values by ",
       signif(max(chk$diff), 3), " -- the deposited dosages have changed, so ",
       "panel E is no longer the checked quantity")
msg("  agrees with plink to ", signif(max(chk$diff), 2),
    "; pairwise n = ", paste(sort(unique(c(
      outer(seq_len(ncol(M)), seq_len(ncol(M)),
            Vectorize(function(i, j) sum(complete.cases(M[, c(i, j)]))))))),
      collapse = " and "))

## lower triangle only, in panel order, so a cell maps onto a pair of panels.
## Columns are markers 1..n-1 and rows are markers 2..n: carrying all n levels
## on both axes leaves an empty top row and an empty right column.
lv <- sub(" Mb$", "", st$short)
n <- nrow(st)
ld <- expand_grid(row = seq_len(n), col = seq_len(n)) %>%
  filter(row > col) %>%
  mutate(a = st$marker[row], b = st$marker[col],
         r2 = map2_dbl(a, b, ~ R2[.x, .y]),
         rl = factor(lv[row], levels = rev(lv[2:n])),
         cl = factor(lv[col], levels = lv[1:(n - 1)]),
         lab = ifelse(r2 < 0.001, format(r2, digits = 2, scientific = TRUE),
                      sprintf("%.3f", r2)))
LD_MAX <- max(ld$r2)
msg("  strongest pair ", ld$a[which.max(ld$r2)], " / ",
    ld$b[which.max(ld$r2)], " at r2 ", sprintf("%.3f", LD_MAX),
    " -- ", sprintf("%.0f", LD_INTERVAL_R2 / LD_MAX),
    "x below the interval criterion of ", LD_INTERVAL_R2)

pE <- ggplot(ld, aes(cl, rl, fill = r2)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_text(aes(label = lab, colour = r2 > 0.06), size = 3.1,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_gradient(low = COL_LD_LO, high = COL_LD_HI,
                      limits = c(0, LD_MAX),
                      breaks = c(0, LD_MAX / 2, LD_MAX),
                      labels = sprintf("%.3f", c(0, LD_MAX / 2, LD_MAX)),
                      name = expression(r^2)) +
  scale_colour_manual(values = c(`FALSE` = "grey25", `TRUE` = "white")) +
  scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
  coord_fixed() +
  labs(x = NULL, y = NULL,
       title = panel_title("E", "**LD between the four markers**"),
       subtitle = emph(wrap_md(sprintf(paste0(
         ## the flag is spelled out in the header and the caption; gridtext
         ## turns a double hyphen into an en dash, so it stays out of here
         "Squared correlation of the dosages, positions in Mb, same 231 ",
         "strains; agrees exactly with plink's own estimator. Note the ",
         "scale: the ",
         "strongest pair is %.3f and the interval work calls markers linked ",
         "at r-squared %.1f, so no pair here is close to being one signal. ",
         "The two chromosome III markers, 6.75 Mb apart, are the least ",
         "correlated pair of the six."),
         LD_MAX, LD_INTERVAL_R2), 88))) +
  theme_pub(10) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_text(size = 8.2, colour = "grey20"),
        plot.subtitle = element_markdown(size = 6.9, colour = "grey30",
                                         lineheight = 1.3),
        plot.title = element_markdown(size = 10.5),
        legend.position = "right",
        legend.title = element_text(size = 8.4),
        legend.text = element_text(size = 7.2),
        legend.key.height = grid::unit(20, "pt"))

## ===========================================================================
panels <- map(seq_len(nrow(st)), build)
## the triangle is a small square panel, so it is centred under the four
## splits rather than stretched across their full width
fig <- wrap_plots(panels, nrow = 1) /
  (plot_spacer() | pE | plot_spacer()) +
  plot_layout(heights = c(1, 0.62))

ggsave(file.path(OUT, "SUPP_FIG_XX_gwas_peak_genotype_splits.pdf"), fig,
       width = 11.4, height = 6.9, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_gwas_peak_genotype_splits.png"), fig,
       width = 11.4, height = 6.9, dpi = 300, bg = "white")
msg("wrote SUPP_FIG_XX_gwas_peak_genotype_splits.{pdf,png}")
