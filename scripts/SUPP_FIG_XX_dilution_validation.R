## Supplement -- the DNA dilution validation ---------------------------------
##
##   Rscript scripts/SUPP_FIG_XX_dilution_validation.R
##     -> plots/SUPP_FIG_XX_dilution_validation.{pdf,png}
##
##   A  pure pools: the fraction of each single-set pool assigned to its own
##      set, and the fraction misassigned to each of the other three
##   B  the B-into-C dilution series, BC1 to BC7, on the pool-wide reference
##   C  the same series with the reference restricted to sets B and C and
##      renormalised, which is the analysis the original figure showed
##
## THE EXPERIMENT. 174 wild isolates were split into four sets of roughly equal
## size (A, B, C, D). Genomic DNA from each set was pooled, and the pools were
## sequenced both pure -- A0/A1/A2 and so on, three libraries per set -- and as
## a seven-step titration of set B against set C, BC1 to BC7. Alt-allele counts
## were taken with GATK ASEReadCounter and deconvolved to per-strain
## frequencies by NNLS against the CeNDR genotype matrix.
##
## This is the first of three validations of the deconvolution, and the only
## one on real sequence with a designed input:
##   simulation                  known input, simulated counts
##   THIS, the DNA dilution      known input, real counts
##   Figure 1A                   real input, real counts, published MIP-seq
##
## THE NOMINAL MIXING RATIOS ARE NOT RECORDED. Nothing in the transferred
## experiment folder states the intended B:C proportions for BC1-BC7, so this
## figure shows that recovery is MONOTONIC and COMPLEMENTARY, not that it is
## ACCURATE. Panel B is the ordering; it is not an accuracy plot, and no
## regression against a nominal series is drawn because there is no nominal
## series on disk. If the design is recovered from the lab record, a
## nominal-against-observed panel is the obvious addition and the data here
## support it directly.
##
## TWO GENOTYPE REFERENCES, DELIBERATELY. Panels A and B use the reference
## restricted to the 170 strains that are actually in the pools; panel C uses
## the reference restricted to the 84 strains of sets B and C alone, which is
## what the original analysis did. The two disagree by a few percent on the
## same samples -- see the console output -- and that disagreement is itself
## worth reporting, because it is the size of the reference-choice effect.
## supplemental_data/deconvolution also carries the 540-strain full-panel and
## regenotyped references for anyone who wants to extend that comparison.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggtext)
  library(ggrepel)
})

OUT  <- "plots"
DEC  <- "supplemental_data/deconvolution"
SETS <- file.path(DEC, "dilution_strain_sets.tsv")
POOL <- file.path(DEC, "dilution_predictions_poolref.tsv.gz")
BCR  <- file.path(DEC, "dilution_predictions_bcref.tsv.gz")

stopifnot(file.exists(SETS), file.exists(POOL), file.exists(BCR))

## sets B and C are the titrated pair and carry the colour; A and D are the
## untitrated sets and read as the off-target floor
SET_COL <- c(A = "grey72", B = "#B5446E", C = "#3F7CAC", D = "grey52")

## ---------------------------------------------------------------------------
## THE PLOTTED FRACTIONS, PINNED
##
## The set fractions the figure draws are written out here as literals and
## checked against what the deposit computes, so the numbers quoted in the
## caption and in the manuscript text cannot drift from the figure without this
## script failing. Update BOTH the literals and the caption if the deposit is
## ever rebuilt with a different reference or a corrected set assignment.
##
## Panel B -- 170-strain pool reference (sets A and D also present, not pinned)
## Panel C -- 84-strain B+C reference, renormalised; the original analysis
## ---------------------------------------------------------------------------
PINNED_B <- tibble::tribble(          # panel B, pool reference
  ~sample, ~B,    ~C,
  "BC1",   0.121, 0.701,
  "BC2",   0.157, 0.669,
  "BC3",   0.211, 0.614,
  "BC4",   0.395, 0.413,
  "BC5",   0.530, 0.281,
  "BC6",   0.678, 0.144,
  "BC7",   0.724, 0.100)

PINNED_C <- tibble::tribble(          # panel C, B+C reference, renormalised
  ~sample, ~B,     ~C,
  "BC1",   0.1713, 0.8287,
  "BC2",   0.2095, 0.7905,
  "BC3",   0.2847, 0.7153,
  "BC4",   0.4775, 0.5225,
  "BC5",   0.6200, 0.3800,
  "BC6",   0.7706, 0.2294,
  "BC7",   0.8433, 0.1567)

## NOMINAL DNA MIXING RATIOS ARE UNKNOWN. If the lab record turns up, put the
## intended B fraction for BC1-BC7 here and the figure gains an
## observed-against-nominal comparison; until then there is nothing to compare
## the recovery against and the panels show ordering, not accuracy.
NOMINAL_B <- NULL

check_pinned <- function(computed, pinned, what, tol) {
  cmp <- computed %>% filter(set %in% c("B", "C")) %>%
    select(sample, set, f) %>%
    pivot_wider(names_from = set, values_from = f) %>%
    inner_join(pinned, by = "sample", suffix = c("_got", "_pin"))
  stopifnot(nrow(cmp) == nrow(pinned))
  d <- max(abs(c(cmp$B_got - cmp$B_pin, cmp$C_got - cmp$C_pin)))
  if (d > tol)
    stop(what, ": pinned values are stale, max difference ", signif(d, 3),
         " (tolerance ", tol, "). Update the literals and the caption.")
  cat(what, ": pinned values agree, max difference ", signif(d, 2), "\n", sep = "")
  invisible(TRUE)
}

panel_title <- function(letter)
  paste0("<span style='font-size:13pt;color:#111111'>**", letter, "**</span>")

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = base_size - 1),
          plot.title = element_markdown(size = base_size + 0.5),
          plot.title.position = "plot",
          legend.key.size = grid::unit(9, "pt"))
}

## ---------------------------------------------------------------------------
## ONE ISOTYPE IS ASSIGNED TO TWO SETS, AND IT HAS TO BE RESOLVED FIRST
##
## Strain JU1580 and strain JU1793 are the same isotype (JU1793), and the
## experiment put them in different sets -- JU1580 in B, JU1793 in D. The
## deconvolution works at isotype level and returns ONE frequency for them, so
## the source prediction tables carry isotype JU1793 twice, once labelled B and
## once labelled D, and a naive set sum adds that single estimate to both.
##
## Uncorrected this is not cosmetic: JU1793's estimate is 5-19% of set D's
## apparent frequency in BC1-BC7 and it RISES across the titration, so set D
## appears to drift upward when in fact one set-B isotype is being counted in
## it. The rise is the tell -- D is not titrated and cannot respond to the
## series.
##
## Resolved to set B, which is what the original analysis did (its B+C
## reference contains JU1793 and not JU1580) and what the estimate itself says:
## JU1793 tracks set B across the titration. The assertion below fails if any
## other isotype ever becomes ambiguous.
## ---------------------------------------------------------------------------
AMBIGUOUS <- c(JU1793 = "B")

sets_raw <- read_tsv(SETS, show_col_types = FALSE)
iso2set <- sets_raw %>% filter(!is.na(isotype)) %>%
  distinct(isotype, set) %>%
  mutate(set = ifelse(isotype %in% names(AMBIGUOUS),
                      AMBIGUOUS[isotype], set)) %>%
  distinct(isotype, set)
stopifnot(!anyDuplicated(iso2set$isotype))

pool <- read_tsv(POOL, show_col_types = FALSE) %>%
  select(-set) %>% distinct() %>%
  left_join(iso2set, by = c("strain" = "isotype"))
stopifnot(!anyNA(pool$set),
          nrow(pool) == n_distinct(pool$strain) * n_distinct(pool$sample))
cat("pool reference: ", n_distinct(pool$strain), " isotypes x ",
    n_distinct(pool$sample), " samples", 
    " (JU1793 resolved to set ", AMBIGUOUS[["JU1793"]], ")\n\n", sep = "")

setf <- pool %>% group_by(sample, set) %>%
  summarise(f = sum(frequency), .groups = "drop")

## ===========================================================================
## A -- pure pools
## ===========================================================================
pure <- setf %>%
  filter(grepl("^[ABCD][012]$", sample)) %>%
  mutate(expected = substr(sample, 1, 1),
         kind = ifelse(set == expected, "Own set", "Other sets"))

own <- pure %>% filter(kind == "Own set")
cat("== pure pools: fraction assigned to the pool's own set ==\n")
print(as.data.frame(own %>% group_by(expected) %>%
  summarise(n = n(), mean = round(mean(f), 3),
            range = sprintf("%.3f-%.3f", min(f), max(f)), .groups = "drop")),
  row.names = FALSE)
off <- pure %>% filter(kind == "Other sets")
cat(sprintf("  off-target per wrong set: %.3f-%.3f (mean %.3f)\n",
            min(off$f), max(off$f), mean(off$f)))
cat(sprintf("  total off-target per pool: %.3f-%.3f\n\n",
            min(tapply(off$f, off$sample, sum)),
            max(tapply(off$f, off$sample, sum))))

pA <- ggplot(pure, aes(expected, f)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey85") +
  geom_point(data = ~ filter(.x, kind == "Other sets"),
             aes(colour = set), shape = 1, size = 1.5, stroke = 0.5,
             position = position_jitter(width = 0.13, height = 0, seed = 1)) +
  geom_point(data = ~ filter(.x, kind == "Own set"),
             aes(fill = set), shape = 21, size = 2.6, stroke = 0.35,
             colour = "grey20",
             position = position_jitter(width = 0.07, height = 0, seed = 1)) +
  scale_fill_manual(values = SET_COL, guide = "none") +
  scale_colour_manual(values = SET_COL, name = "Assigned to") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Pure pool (three libraries each)",
       y = "Pool assigned to set",
       title = panel_title("A")) +
  theme_pub() +
  theme(legend.position = c(0.5, 0.42), legend.direction = "horizontal",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8))

## ===========================================================================
## B -- the dilution series, pool-wide reference
## ===========================================================================
bc <- setf %>% filter(grepl("^BC[1-7]$", sample)) %>%
  mutate(step = as.integer(sub("BC", "", sample)))

cat("== dilution series, 170-strain pool reference ==\n")
print(as.data.frame(bc %>% select(sample, set, f) %>%
  pivot_wider(names_from = set, values_from = f) %>%
  mutate(across(-sample, ~ round(.x, 3)))), row.names = FALSE)

mono <- bc %>% filter(set %in% c("B", "C")) %>% arrange(set, step) %>%
  group_by(set) %>%
  summarise(monotonic = all(diff(f) > 0) || all(diff(f) < 0),
            spearman = cor(step, f, method = "spearman"), .groups = "drop")
cat("  monotonic across BC1-BC7: ",
    paste(sprintf("%s %s (rho %+.0f)", mono$set, mono$monotonic, mono$spearman),
          collapse = " | "), "\n\n", sep = "")

check_pinned(bc, PINNED_B, "panel B", tol = 0.001)

## value labels pointing at each B and C dot, as the original figure did.
## seed is fixed because ggrepel's placement is stochastic and the repo
## requires every figure to be byte-identical run to run.
lab_b <- bc %>% filter(set %in% c("B", "C")) %>%
  mutate(lab = sprintf("%.2f", f))

pB <- ggplot(bc, aes(step, f, colour = set)) +
  geom_line(data = ~ filter(.x, set %in% c("B", "C")), linewidth = 0.6) +
  geom_line(data = ~ filter(.x, set %in% c("A", "D")), linewidth = 0.4,
            linetype = "dotted") +
  geom_point(size = 2) +
  geom_text_repel(data = lab_b, aes(label = lab), size = 2.5,
                  seed = 1, box.padding = 0.28, point.padding = 0.18,
                  min.segment.length = 0, segment.size = 0.22,
                  segment.colour = "grey60", direction = "y",
                  nudge_x = ifelse(lab_b$set == "B", -0.28, 0.28),
                  show.legend = FALSE) +
  scale_x_continuous(breaks = 1:7, labels = paste0("BC", 1:7)) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_colour_manual(values = SET_COL, name = NULL) +
  labs(x = "Titration step", y = "Pool assigned to set",
       title = panel_title("B")) +
  theme_pub() +
  theme(legend.position = c(0.5, 0.99), legend.direction = "horizontal",
        legend.text = element_text(size = 8))

## ===========================================================================
## C -- the same series on the B+C reference, as originally analysed
## ===========================================================================
bcref <- read_tsv(BCR, show_col_types = FALSE) %>%
  left_join(iso2set, by = c("strain" = "isotype")) %>%
  filter(!is.na(set)) %>%
  group_by(sample, set) %>% summarise(f = sum(frequency), .groups = "drop") %>%
  mutate(step = as.integer(sub("BC", "", sample)))

cat("== dilution series, 84-strain B+C reference (the original analysis) ==\n")
print(as.data.frame(bcref %>% select(sample, set, f) %>%
  pivot_wider(names_from = set, values_from = f) %>%
  mutate(across(-sample, ~ round(.x, 4)))), row.names = FALSE)

## the reference-choice effect, on the same samples
side <- bc %>% filter(set %in% c("B", "C")) %>%
  select(sample, set, poolref = f) %>%
  left_join(bcref %>% select(sample, set, bcref = f), by = c("sample", "set")) %>%
  mutate(poolref_renorm = poolref / ave(poolref, sample, FUN = sum),
         diff = poolref_renorm - bcref)
cat("  reference-choice effect, after renormalising the pool reference",
    " within B+C:\n", sep = "")
cat(sprintf("    max difference %.3f, mean %.3f\n\n",
            max(abs(side$diff)), mean(abs(side$diff))))

check_pinned(bcref, PINNED_C, "panel C", tol = 0.001)

lab_c <- bcref %>% mutate(lab = sprintf("%.2f", f))

pC <- ggplot(bcref, aes(step, f, colour = set)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3,
             colour = "grey70") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  geom_text_repel(data = lab_c, aes(label = lab), size = 2.5,
                  seed = 1, box.padding = 0.28, point.padding = 0.18,
                  min.segment.length = 0, segment.size = 0.22,
                  segment.colour = "grey60", direction = "y",
                  nudge_x = ifelse(lab_c$set == "B", -0.3, 0.3),
                  show.legend = FALSE) +
  scale_x_continuous(breaks = 1:7, labels = paste0("BC", 1:7)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_colour_manual(values = SET_COL, name = NULL) +
  labs(x = "Titration step", y = "Set share, renormalised",
       title = panel_title("C")) +
  theme_pub() +
  theme(legend.position = c(0.5, 0.99), legend.direction = "horizontal",
        legend.text = element_text(size = 8))

## ===========================================================================
## HOW FAR OFF WAS THE INFERENCE?
##
## Three metrics, in decreasing order of how much they assume.
##
## 1. CONSERVATION OF THE TITRATED PAIR -- assumes nothing about the design.
##    Only B and C are titrated against each other, so whatever the intended
##    proportions were, the pair's combined share of the pool must stay constant
##    across BC1-BC7: DNA moved from C to B cannot change how much of the pool
##    is B-or-C. Any wobble in B+C is inference error, and it is measurable
##    without knowing a single nominal ratio. This is the metric to quote.
##
## 2. DEVIATION FROM A NOMINAL SERIES -- assumes the design. Computed the
##    moment NOMINAL_B is filled in; until then it is reported against evenly
##    spaced reference series, and labelled as conditional, so the machinery is
##    demonstrably working and the number is available the day the lab record
##    turns up. It is NOT evidence about accuracy while the design is unknown.
##
## 3. LINEARITY -- assumes the design was an even series, which is a guess.
##    Reported because it is the shape the recovered values take, not because
##    an even series is known to be what was mixed.
## ===========================================================================
cat("== metric 1: conservation of the titrated pair (assumes no design) ==\n")
pair <- bc %>% select(sample, set, f) %>%
  pivot_wider(names_from = set, values_from = f) %>%
  mutate(pair = B + C, other = A + D)
print(as.data.frame(pair %>% transmute(sample, B = round(B, 4), C = round(C, 4),
      `B+C` = round(pair, 4), `A+D` = round(other, 4))), row.names = FALSE)
cons <- with(pair, c(mean = mean(pair), sd = sd(pair),
                     max_dev = max(abs(pair - mean(pair)))))
cat(sprintf(paste0("  B+C mean %.4f, sd %.4f (%.2f%% of the mean), ",
                   "max deviation %.4f (%.2f%%)\n"),
            cons[["mean"]], cons[["sd"]], 100 * cons[["sd"]] / cons[["mean"]],
            cons[["max_dev"]], 100 * cons[["max_dev"]] / cons[["mean"]]))
cat("  -> as B is titrated into C the pair's total is conserved to under 1.5%,",
    "\n     which bounds the inference error without invoking the design.\n")

cat("\n== metric 2: deviation from a nominal series (assumes the design) ==\n")
dev_against <- function(observed, nominal, label) {
  d <- observed - nominal
  cat(sprintf("  %-22s RMSE %.4f | max |dev| %.4f | bias %+.4f\n",
              label, sqrt(mean(d^2)), max(abs(d)), mean(d)))
  invisible(c(rmse = sqrt(mean(d^2)), max = max(abs(d)), bias = mean(d)))
}
obs_B <- bcref %>% filter(set == "B") %>% arrange(step) %>% pull(f)
if (!is.null(NOMINAL_B)) {
  stopifnot(length(NOMINAL_B) == 7)
  dev_against(obs_B, NOMINAL_B, "NOMINAL_B (recorded)")
} else {
  cat("  NOMINAL_B is NULL -- the design is not on disk. Reported below",
      "against\n  evenly spaced reference series, which is an ASSUMPTION and",
      "not the design:\n")
  dev_against(obs_B, (1:7) / 8,                 "even 1/8..7/8")
  dev_against(obs_B, seq(0.10, 0.90, length.out = 7), "even 0.10..0.90")
  dev_against(obs_B, seq(0.15, 0.85, length.out = 7), "even 0.15..0.85")
  cat("  All three land at RMSE 0.044-0.048, so the metric is not sensitive to\n",
      "  which even series is assumed -- but none of them is known to be right.\n",
      sep = "")
}

cat("\n== metric 3: shape of the recovered series ==\n")
fit <- lm(f ~ step, data = filter(bcref, set == "B"))
cat(sprintf("  B share against step: R-squared %.4f, slope %+.4f per step, ",
            summary(fit)$r.squared, coef(fit)[["step"]]))
cat(sprintf("residual sd %.4f\n", sd(resid(fit))))
cat(sprintf("  Spearman rho against step: %+.0f\n\n",
            cor(filter(bcref, set == "B")$step,
                filter(bcref, set == "B")$f, method = "spearman")))

fig <- pA | pB | pC

ggsave(file.path(OUT, "SUPP_FIG_XX_dilution_validation.pdf"), fig,
       width = 11, height = 3.7, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_dilution_validation.png"), fig,
       width = 11, height = 3.7, dpi = 300, bg = "white")
cat("wrote SUPP_FIG_XX_dilution_validation.{pdf,png}\n")
