## Supplement -- the simulation that started the project ---------------------
##
##   Rscript scripts/SUPP_FIG_XX_simulation_depth.R
##     -> plots/SUPP_FIG_XX_simulation_depth.{pdf,png}
##
##   A  r-squared of NNLS-estimated against KNOWN INPUT frequency, against
##      simulated sequencing depth, one line per simulated trait
##   B  the same estimates plotted against the known input, faceted by depth
##
## THE SIMULATION. Each wild isolate was assigned a fitness value taken from one
## of seven published C. elegans traits with validated QTL -- the trait value
## itself, not a draw -- shifted by the trait's own minimum so fitness is
## non-negative, with strains carrying no value for a trait set to zero and so
## absent from that trait's pool. The expected pooled allele frequencies such a
## population would produce were computed; observed alt-allele counts were
## simulated by binomial sampling at depths 1, 3, 5, 10, 30, 50, 100 and 500x;
## and those counts were deconvolved back to per-strain frequencies by NNLS. The
## seven traits are Albendazole_q75.TOF, PC1, assay_norm, mtDNA_ratio, value,
## amsacrine_f.L1 and etoposide_median.TOF.
##
## BOTH PANELS ARE COMPUTED AGAINST THE KNOWN INPUT (2026-09-09)
## -------------------------------------------------------------
## This figure used to report panel A rather than compute it, and to substitute
## the 500x estimate for the truth in panel B, because the simulation's fitness
## input was thought lost. It was not: the seven-trait arm of
## scripts/legacy/haploReg_original.R uses published trait values as fitness,
## and those files are now deposited as
## supplemental_data/deconvolution/simulation_fitness_traits.tsv by
## scripts/make_simulation_fitness_table.R.
##
## So the known input is available and both panels use it. The old record is not
## discarded -- it is asserted against. Every one of the 56 r-squared computed
## here reproduces simulation_reported_r2.tsv, the values transcribed from text
## embedded in the original per-trait PDFs, at the two decimals those PDFs
## carry; the stopifnot below fails if that ever stops being true. See
## scripts/simulation_recompute_r2.R, which is the standalone check.
##
## Two consequences for reading the panels. Panel B now covers all eight depths
## rather than seven: 500x was previously excluded because it WAS the reference,
## and is now just another depth. And the truth has many exact zeros -- the
## strains carrying no published value for a trait -- so the pile of points on
## y = 0 is real, and is where NNLS assigns frequency to a strain that was not
## in the pool at all.
##
## SCALE. The stored coefficients are not frequencies: each depth's row sums to
## the depth itself (500x sums to 500.03, 1x to 1.00), so they are on an
## expected-count scale. frequency = coefficient / sum(coefficient) within a
## trait and depth, which is how the deposit's frequency column was built.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggtext)
})

OUT <- "plots"
DEC <- "supplemental_data/deconvolution"
R2  <- file.path(DEC, "simulation_reported_r2.tsv")            # the 2021 record
FRQ <- file.path(DEC, "simulation_nnls_frequencies.tsv.gz")   # the 2021 estimates
FIT <- file.path(DEC, "simulation_fitness_traits.tsv")
SR2 <- file.path(DEC, "simulation_seeded_r2.tsv")             # this figure's source
SFQ <- file.path(DEC, "simulation_seeded_frequencies.tsv.gz")

stopifnot(file.exists(R2), file.exists(FRQ), file.exists(FIT),
          file.exists(SR2), file.exists(SFQ))

## seven traits need seven distinguishable hues, and the previous set ran
## teal -> green -> olive -> gold -> rust -> wine -> purple, which collapses to
## about three under deuteranopia. Seven of the eight Okabe-Ito hues instead.
TRAIT_COL <- c(
  `PC1`                  = "#0072B2",
  `value`                = "#009E73",
  `amsacrine_f.L1`       = "#56B4E9",
  `assay_norm`           = "#E69F00",
  `etoposide_median.TOF` = "#D55E00",
  `Albendazole_q75.TOF`  = "#CC79A7",
  `mtDNA_ratio`          = "#999999")

panel_title <- function(letter)
  paste0("<span style='font-size:13pt;color:#111111'>**", letter, "**</span>")

## one theme and one type scale; see scripts/figure_theme.R
source("scripts/figure_theme.R")

DEPTHS <- c(1, 3, 5, 10, 30, 50, 100, 500)

## ===========================================================================
## the known input: fitness, then the pool frequency it implies
## ===========================================================================
## The original's two transforms (haploReg_original.R lines 327-328): shift each
## trait by its own minimum so fitness is non-negative, then send NA to 0 so a
## strain with no published value is absent from that trait's pool. Dividing by
## the column total puts it on the same frequency scale as the estimates, which
## is what lets panel B draw a y = x line.
## The seeded, replicated run is what both panels plot. It is generated from
## the genotype panel and the recovered fitness input by
## scripts/make_simulation_seeded.R, so it is reproducible from a seed rather
## than recovered from a 2021 RDS whose draw was never seeded.
sr2 <- read_tsv(SR2, show_col_types = FALSE)
sfq <- read_tsv(SFQ, show_col_types = FALSE)
NREP <- dplyr::n_distinct(sr2$replicate)
stopifnot(nrow(sr2) == NREP * 7 * 8, nrow(sfq) == NREP * 7 * 8 * 327)
## the estimator emits no negative coefficients; the old archive carried 139
stopifnot(all(sfq$frequency >= 0))
msg <- function(...) cat(...) # local, for the notes below

frq <- read_tsv(FRQ, show_col_types = FALSE)
stopifnot(nrow(frq) == 7 * 8 * 327)

panel_strains <- sort(unique(frq$strain))
truth <- read_tsv(FIT, show_col_types = FALSE) %>%
  filter(strain %in% panel_strains) %>%
  pivot_longer(-strain, names_to = "trait", values_to = "published") %>%
  group_by(trait) %>%
  mutate(fitness = published - min(published, na.rm = TRUE),
         fitness = replace_na(fitness, 0),
         input   = fitness / sum(fitness)) %>%
  ungroup() %>%
  select(trait, strain, published, fitness, input)
stopifnot(nrow(truth) == 7 * 327, !anyNA(truth$input))

acc <- frq %>%
  inner_join(truth, by = c("trait", "strain")) %>%
  group_by(trait, depth) %>%
  summarise(r2 = cor(frequency, input)^2, .groups = "drop")
stopifnot(nrow(acc) == 56)

## ===========================================================================
## A -- accuracy against the known input, computed
## ===========================================================================
## The transcribed record is the check on this, not the source of it: if the
## recomputation ever stops reproducing the 2-dp values read out of the 2021
## PDFs, that is a real regression and this figure must not build.
rep <- read_tsv(R2, show_col_types = FALSE)
chk <- acc %>% inner_join(rep, by = c("trait", "depth"),
                          suffix = c(".computed", ".reported"))
stopifnot(nrow(chk) == 56,
          all(round(chk$r2.computed, 2) == chk$r2.reported))
cat("== all 56 computed r-squared reproduce simulation_reported_r2.tsv at 2 dp ==\n\n")

## Panel A is the mean over replicates with a min-max ribbon. The single
## 2021 draw reported each recovery as one exact number; the ribbon is the
## sampling variability that number could only conceal, and it matters -- a
## trait can straddle the 0.95 line, so "the lowest depth at which all seven
## reach 0.95" is one depth in some draws and another in others and cannot
## honestly be quoted as a single depth.
##
## BASIS: panel A is computed over the strains that were actually IN each
## trait's pool (input > 0). The transform sends a strain with no published
## value for a trait to fitness exactly 0, so it is absent from that pool, and
## NNLS returns exactly 0 for an absent strain because zero is a hard boundary
## of the feasible set. Those strains sit on the origin and inflate r-squared
## without saying anything about measurement precision: they answer "did it
## notice the strain was missing", not "how precisely did it measure the
## strains that were there". Across the seven traits 49.4% of strain-trait
## cells are absent, so the inflation is large -- the all-strains mean at 10x
## is 0.95 against 0.88 here.
##
## The all-strains series is still computed below as `r2_all`, because the 2021
## record is on that basis and the check against it is the regression guard on
## this whole figure. It is printed, not plotted.
## lo and hi BEFORE r2: summarise() evaluates in order, so naming the mean `r2`
## first would shadow the replicate column and collapse the band to zero width.
sr2_present <- sfq %>% filter(input > 0) %>%
  group_by(trait, depth, replicate) %>%
  summarise(r2 = cor(frequency, input)^2, .groups = "drop")
stopifnot(nrow(sr2_present) == NREP * 7 * 8)

r2 <- sr2_present %>% group_by(trait, depth) %>%
  summarise(lo = min(r2), hi = max(r2), r2 = mean(r2), .groups = "drop") %>%
  mutate(trait = factor(trait, levels = names(TRAIT_COL)))
stopifnot(any(r2$hi > r2$lo))   # a zero-width band means the bug is back
stopifnot(!anyNA(r2$trait))

r2_all <- sr2 %>% group_by(trait, depth) %>%
  summarise(lo = min(r2), hi = max(r2), r2 = mean(r2), .groups = "drop")

## the archived run is on the all-strains basis, so it is checked against the
## all-strains band -- comparing it with the present-only band would be a
## category error, not a regression
arch <- acc %>% rename(archive = r2)
chk2 <- r2_all %>% inner_join(arch, by = c("trait","depth")) %>%
  mutate(inside = archive >= lo & archive <= hi)
cat(sprintf("== the 2021 draw falls inside the %d-replicate all-strains range in %d of %d cells ==\n",
            NREP, sum(chk2$inside), nrow(chk2)))
cat(sprintf("== basis: absent strains are %.1f%% of strain-trait cells; ",
            100 * mean(sfq$input == 0)))
cat(sprintf("mean r2 at 10x is %.3f all-strains vs %.3f present-only ==\n\n",
            r2_all$r2[r2_all$depth == 10] %>% mean(),
            r2$r2[r2$depth == 10] %>% mean()))

at1 <- r2 %>% filter(depth == 1) %>% arrange(desc(r2))
cat("== accuracy at 1x, the depth the text claims ==\n")
print(as.data.frame(at1 %>% transmute(trait, r2 = round(r2, 4))), row.names = FALSE)
cat(sprintf("  range %.2f-%.2f, median %.2f\n\n",
            min(at1$r2), max(at1$r2), median(at1$r2)))

## the 1.00s in the transcribed record are rounding, and the figure should not
## imply otherwise -- state the true ceiling
at500 <- r2 %>% filter(depth == 500)
cat(sprintf("== at 500x no trait reaches 1.000: %.4f (%s) to %.4f (%s) ==\n\n",
            min(at500$r2), at500$trait[which.min(at500$r2)],
            max(at500$r2), at500$trait[which.max(at500$r2)]))

## the depth at which every trait first reaches 0.95 and stays there
reach <- r2 %>% arrange(trait, depth) %>% group_by(trait) %>%
  summarise(first95 = {
    ok <- r2 >= 0.95
    idx <- which(ok & rev(cumprod(rev(ok))) == 1)[1]   # first depth from which it never drops
    if (is.na(idx)) NA_integer_ else depth[idx]
  }, .groups = "drop")
cat("== lowest depth from which the MEAN stays >= 0.95 ==\n")
print(as.data.frame(reach), row.names = FALSE)

## The mean crossing 0.95 is not the claim worth making, because a trait can sit
## on the line. Robustness -- clearing in EVERY replicate -- is what the
## caption states, so it is what gets asserted here. On the present-only basis
## this is a markedly harder bar than the all-strains basis made it look.
rob <- sr2_present %>% group_by(trait, depth) %>%
  summarise(all_reps = all(r2 >= 0.95), n_of = sum(r2 >= 0.95),
            reps = dplyr::n(), .groups = "drop")
by_depth <- rob %>% group_by(depth) %>%
  summarise(robust = sum(all_reps), .groups = "drop") %>% arrange(depth)
cat("\n== traits clearing 0.95 in EVERY replicate, by depth ==\n")
print(as.data.frame(by_depth), row.names = FALSE)
borderline <- rob %>% filter(depth == 30, !all_reps)
if (nrow(borderline))
  cat(sprintf("  at 30x %s clears 0.95 in only %d of %d replicates\n",
              borderline$trait[1], borderline$n_of[1], borderline$reps[1]))
cat(sprintf("  four of seven clear 0.95 in every replicate by 30x; %d by 100x; all seven only at 500x\n",
            by_depth$robust[by_depth$depth == 100]))
cat(sprintf("  lowest single-replicate r2: %.4f at 30x, %.4f at 50x, %.4f at 100x\n\n",
            min(sr2_present$r2[sr2_present$depth == 30]),
            min(sr2_present$r2[sr2_present$depth == 50]),
            min(sr2_present$r2[sr2_present$depth == 100])))
stopifnot(by_depth$robust[by_depth$depth == 10] == 2,
          by_depth$robust[by_depth$depth == 30] == 4,
          by_depth$robust[by_depth$depth == 50] == 5,
          by_depth$robust[by_depth$depth == 100] == 6,
          by_depth$robust[by_depth$depth == 500] == 7)

lab <- r2 %>% filter(depth == 1)
pA <- ggplot(r2, aes(depth, r2, colour = trait)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.35,
             colour = "grey55") +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = trait), alpha = 0.18,
              colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.3) +
  annotate("text", x = 1.05, y = 0.978, label = "r² = 0.95", hjust = 0,
           size = 2.7, colour = "grey40") +
  scale_x_log10(breaks = DEPTHS, labels = paste0(DEPTHS, "×")) +
  scale_y_continuous(limits = c(0.2, 1.005), breaks = seq(0.2, 1, 0.1)) +
  scale_colour_manual(values = TRAIT_COL, name = NULL) +
  scale_fill_manual(values = TRAIT_COL, guide = "none") +
  labs(x = "Simulated sequencing depth",
       y = "r² vs known input frequency,\nstrains present in the pool",
       title = panel_title("A"),
       subtitle = sprintf(paste("Accuracy against the simulated input, over the strains actually IN each",
                                "trait's pool.<br>Strains absent from a pool are returned as exactly zero and",
                                "are excluded: they inflate r²<br>without measuring anything. Mean of %d seeded",
                                "replicates, band spans them; dashed line r² = 0.95"), NREP)) +
  theme_pub() +
  theme(legend.position = c(0.985, 0.02), legend.justification = c(1, 0),
        legend.text = element_text(size = 7.6),
        legend.background = element_rect(fill = alpha("white", 0.85),
                                         colour = NA),
        panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey92"))

## ===========================================================================
## B -- the estimates against the known input, every depth
## ===========================================================================
## All eight depths appear here now. 500x used to be excluded because it was
## itself the reference; against the true input it is simply the deepest point.
## One replicate, not an average: averaging point clouds would narrow the
## scatter and misrepresent what a single experiment at that depth looks like.
cmp <- sfq %>% filter(replicate == 1) %>%
  mutate(depth_lab = factor(paste0(depth, "×"), levels = paste0(DEPTHS, "×")))
stopifnot(nrow(cmp) == 7 * 8 * 327)

pooled <- cmp %>% group_by(depth) %>%
  summarise(r2 = cor(frequency, input)^2, .groups = "drop") %>% arrange(depth)
cat("== accuracy against the known input, all seven traits pooled ==\n")
print(as.data.frame(pooled %>% mutate(r2 = round(r2, 3))), row.names = FALSE)
cat("\n")

## how much of the pooled fit is carried by strains that were not in the pool
zero <- cmp %>% group_by(depth) %>%
  summarise(absent = sum(input == 0),
            absent.nonzero.est = sum(input == 0 & frequency > 0),
            r2.present.only = cor(frequency[input > 0], input[input > 0])^2,
            .groups = "drop") %>% arrange(depth)
cat("== the same, restricted to strains actually in the pool ==\n")
print(as.data.frame(zero %>% mutate(r2.present.only = round(r2.present.only, 3))),
      row.names = FALSE)
cat("\n")

## Why panel B's cloud sits off the y = x line at low depth: NNLS puts weight on
## strains that were not in the pool, so the strains that WERE are
## underestimated by the same amount. Quantified as the share of each pool's
## mass landing on absent strains, and as the slope of input on estimate.
leak <- cmp %>% group_by(depth) %>%
  summarise(pct.mass.on.absent = 100 * sum(frequency[input == 0]) / 7,
            slope.input.on.est = coef(lm(input ~ frequency))[2],
            .groups = "drop") %>% arrange(depth)
cat("== frequency leaked onto strains that were not in the pool ==\n")
print(as.data.frame(leak %>% mutate(pct.mass.on.absent = round(pct.mass.on.absent, 1),
                                    slope.input.on.est = round(slope.input.on.est, 3))),
      row.names = FALSE)
cat("\n")

ann <- pooled %>% mutate(depth_lab = factor(paste0(depth, "×"),
                                            levels = levels(cmp$depth_lab)),
                         lab = sprintf("r² = %.2f", r2))

pB <- ggplot(cmp, aes(frequency, input)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.35, colour = "grey60") +
  geom_point(shape = 16, size = 0.5, alpha = 0.25, colour = "grey20") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.12, vjust = 1.5, size = 2.7, colour = "grey20") +
  ## Both axes are the same quantity -- an estimated and a true frequency on the
  ## same scale -- so they must share limits: with scales = "free_x" the y = x
  ## line lands somewhere different in every facet and the comparison it is
  ## there to support becomes unreadable.
  facet_wrap(~ depth_lab, nrow = 1) +
  scale_x_continuous(breaks = c(0, 0.01, 0.02, 0.03),
                     labels = c("0", ".01", ".02", ".03")) +
  scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03),
                     labels = c("0", ".01", ".02", ".03")) +
  labs(x = "Estimated strain frequency at the stated depth",
       y = "Known input frequency",
       title = panel_title("B"),
       subtitle = paste("Points on y = 0 are strains with no published value",
                        "for that trait, absent from the simulated pool")) +
  theme_pub() +
  theme(panel.spacing.x = grid::unit(7, "pt"),
        axis.text.x = element_text(size = 7.2))

fig <- pA / pB + plot_layout(heights = c(1.3, 1))

ggsave(file.path(OUT, "SUPP_FIG_XX_simulation_depth.pdf"), fig,
       width = 9.2, height = 6.4, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_simulation_depth.png"), fig,
       width = 9.2, height = 6.4, dpi = 300, bg = "white")
cat("wrote SUPP_FIG_XX_simulation_depth.{pdf,png}\n")

## the caveat, printed so it cannot be missed by anyone re-running this
neg <- frq %>% filter(coefficient < 0)
cat(sprintf("\nNOTE: this figure's own estimates carry %d negative coefficients.\n",
            sum(sfq$frequency < 0)))
cat(sprintf("The 2021 archive, no longer plotted here, carried %d of %d (%.1f%%) at depths %s;\n",
            nrow(neg), nrow(frq), 100 * nrow(neg) / nrow(frq),
            paste(sort(unique(neg$depth)), collapse = ", ")))
cat("a strict non-negative solver cannot return those, and they are one reason\n")
cat("the figure is built from the seeded run instead. See the caption.\n")
