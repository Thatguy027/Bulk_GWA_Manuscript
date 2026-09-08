## Diagnostic -- off-pool leakage against coverage and reference size ---------
##
##   Rscript scripts/diagnostic_downsample_leakage.R
##     -> plots/diagnostics/downsample_leakage.{pdf,png}
##     -> plots/diagnostics/TABLE_downsample_leakage.tsv
##     -> plots/diagnostics/TABLE_downsample_leakage_slopes.tsv
##     -> plots/diagnostics/TABLE_absent_strain_leakage.tsv
##     -> plots/diagnostics/CACHE_absent_strain_ibs.tsv   (rebuilt with --refresh-ibs)
##
## WHAT THIS ADDS THAT diagnostic_reference_size.R DOES NOT. That script
## measures DISCREPANCY -- |NNLS - MIP| for strains that are in the pool. This
## one measures LEAKAGE: frequency assigned to strains that are NOT in the pool
## at all. They are different quantities and only the second is what goes wrong
## when a real experiment cannot name its own members. At reference size 102
## leakage is zero by construction, because every reference column is a pool
## member; at 540 there are 438 columns that should carry nothing.
##
## THE QUESTION BEHIND IT. If you hand the solver 438 candidates that are absent,
## how much of the pool does it hand back to them, does coverage make it worse,
## and -- the part that decides whether any of this matters for the manuscript --
## does it still matter once you restrict the reference to the true pool?
##
## NO NNLS IS RUN HERE. diagnostic_reference_size.R saved every fitted frequency
## for all 540 reference columns across the 5 x 7 design, including the columns
## for absent strains, which its own tables then dropped in an inner_join
## against the MIP-seq measurements. Those predictions are the input here, so
## this script is cheap and cannot drift from the numbers already reported.
##
## FREQUENCIES ARE NORMALISED TO 1 PER SAMPLE (checked below), so leakage reads
## directly as a fraction of the pool and the four reference sizes are on one
## scale.
##
## IBS IS COMPUTED WITHOUT THE FLIP MASK, DELIBERATELY. Leakage is compared
## against each absent strain's identity-by-state to its closest pool member,
## which needs genotypes for all 540. The archive's design matrix is recoded so
## the counted allele is the minor one in the 102-strain panel, and
## diagnostic_reference_size.R recovers that mask before building its enlarged
## reference. Here it is unnecessary: IBS counts markers at which two strains
## MATCH, and complementing a marker in both strains leaves the match unchanged.
## A per-marker flip is a relabelling, so IBS is invariant to it, and the raw
## CeNDR calls give the same answer as the oriented ones. The pool strains'
## values are correlated against the published similarity table as a check that
## this reasoning holds in practice, not just on paper.
##
## WHAT IT FOUND -- pinned here so the reading does not drift.
##
## 1. LEAKAGE IS LARGE IN AGGREGATE. At the full 540-strain reference and full
##    depth, 14.9% of the pool is assigned to strains that are not in it (median
##    across samples; total pool frequency retained 0.8505). It scales with the
##    number of absent columns, not with anything subtler: 3.6% at R=150 (48
##    absent), 7.5% at 250, 11.8% at 400, 14.9% at 540 -- about 0.34 per mille
##    per absent candidate throughout.
##
## 2. COVERAGE BARELY MATTERS. Over a 64-fold reduction in depth, leakage at
##    R=540 rises only from 14.9% to 20.1%, and most of that appears in the last
##    two halvings. Low coverage is not what makes an unknown membership list
##    expensive.
##
## 3. IT IS NOT SPREAD EVENLY. Among strain-samples above 0.5% of the pool the
##    median retained fraction is 0.898, but the 5-95% range is 0.396-1.04, and
##    individual strains are displaced outright: ECA36 goes from 4.17 per mille
##    at R=102 to 0.00 at R=540 while JU3226 -- absent from the pool, IBS 0.9898
##    to ECA36 -- absorbs 4.39 per mille. That is a near-exact swap of one
##    strain for its look-alike.
##
## 4. AND YET IT DOES NOT REACH THE PHENOTYPE. Spearman agreement between NNLS
##    growth slopes and the MIP-seq slopes is 0.975 at R=540 and 0.974 at R=102
##    at full depth; RMSE is 0.0001 for both. Restricting the reference to the
##    true pool and renormalising -- what an experiment with a known membership
##    list gets -- changes the third decimal at most. The cost of not knowing
##    the membership, measured as loss in rho against the R=102 reference, is
##    between -0.001 and +0.003 at every depth above 1/64.
##
##    So the two things are both true and neither cancels the other: individual
##    strain identity degrades badly at the extreme of relatedness, and the
##    slope phenotype that Figure 1 and the pooled GWAS are built on does not
##    notice. Displacement is consistent across the samples of a replicate, so
##    it moves a strain's whole trajectory rather than its trend.
##
## 5. RELATEDNESS GRADES THE MAGNITUDE, NOT THE MEMBERSHIP. Mean frequency
##    absorbed rises 9.6-fold across IBS deciles, from 0.081 per mille in the
##    least related tenth of absent candidates to 0.785 in the most related
##    (R=540, full depth). But a whole-set Spearman reports only +0.051, because
##    NNLS is non-negative and its solutions are sparse: 56.8% of absent
##    candidates absorb EXACTLY zero, and the ties flatten the statistic. Split
##    in two, it resolves cleanly -- IBS does not predict which candidates enter
##    the support (logistic slope -1.13, p = 0.71) and does predict how much
##    those that do absorb take (rho +0.35 at 540, +0.46 at 150). The whole-set
##    rho declining from +0.304 at R=150 to +0.051 at R=540 is entirely the zero
##    fraction climbing from 25% to 57%; the among-absorbers correlation is flat
##    at 0.34-0.46 throughout.
##
##    So: relatedness does direct leakage, and the near-duplicate pairs above
##    IBS ~0.985 are where it becomes wholesale displacement. What is NOT
##    established is which columns the solver admits in the first place -- that
##    is a property of the collinearity of the whole design rather than of any
##    one pairwise distance, and nothing here explains it.
##
## MARKERS: the archive's marker set, so this sits on the same scale as
## supplemental_data/deconvolution/baugh_strain_similarity.tsv; sampled to
## N_MARKERS under a fixed seed because 540 columns of the full set does not fit
## comfortably in memory. The IBS block is cached to a TSV, so the 12 GB
## genotype load happens once.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(ggtext)
})

OUT   <- "plots/diagnostics"
PRED  <- file.path(OUT, "reference_size_predictions.rds")
NNLS  <- "supplemental_data/deconvolution/baugh_nnls_with_mipseq.RData"
SIM   <- "supplemental_data/deconvolution/baugh_strain_similarity.tsv"
BOOT  <- "data/baugh/2024bootstrapINPUT.Rdata"
GENO  <- "data/genotypes/processed_genotype_matrix.Rda"
IBSC  <- file.path(OUT, "CACHE_absent_strain_ibs.tsv")

N_MARKERS <- 2e5
SEED      <- 1
REFRESH   <- "--refresh-ibs" %in% commandArgs(TRUE)

say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep=""); flush.console() }

for (f in c(PRED, NNLS, SIM))
  if (!file.exists(f)) stop("need ", f, if (f == PRED)
    " -- run scripts/diagnostic_reference_size.R first" else "", call. = FALSE)

## ---------------------------------------------------------------------------
## 1. predictions, pool membership, sample metadata
## ---------------------------------------------------------------------------
pred <- readRDS(PRED)
stopifnot(all(c("strain","sample","frq","ref_size","fraction") %in% names(pred)))

tot <- pred %>% group_by(ref_size, fraction, sample) %>%
  summarise(t = sum(frq), .groups = "drop")
if (max(abs(tot$t - 1)) > 1e-8)
  stop("predicted frequencies do not sum to 1 per sample (max deviation ",
       signif(max(abs(tot$t - 1)), 3), "); leakage would not be a fraction ",
       "of the pool and the panels below would not be comparable across ",
       "reference sizes.", call. = FALSE)

POOL      <- sort(unique(pred$strain[pred$ref_size == min(pred$ref_size)]))
REF_SIZES <- sort(unique(pred$ref_size))
FRACTIONS <- sort(unique(pred$fraction))
say(length(POOL), " pool strains | reference sizes ",
    paste(REF_SIZES, collapse = ", "), " | ", length(FRACTIONS), " depths")
stopifnot(length(POOL) == min(REF_SIZES))

e <- new.env(); load(NNLS, e)
arch <- e$wgs_mip_results; rm(e)
meta <- arch %>% distinct(sample, replicate, day, baseline)
stopifnot(!anyNA(meta$day), !anyNA(meta$replicate), nrow(meta) == n_distinct(pred$sample))

## The reference is nested -- each size adds strains to the one below -- so a
## strain's tier is the smallest reference it appears in. Taken from the
## predictions rather than re-derived from the seeded shuffle, so the two cannot
## disagree about who was a candidate at which size.
tier <- pred %>% distinct(strain, ref_size) %>%
  group_by(strain) %>% summarise(enters_at = min(ref_size), .groups = "drop")

## ---------------------------------------------------------------------------
## 2. leakage: how much of the pool goes to strains that are not in it
## ---------------------------------------------------------------------------
leak_sample <- pred %>%
  mutate(absent = !strain %in% POOL) %>%
  group_by(ref_size, fraction, sample) %>%
  summarise(leak = sum(frq[absent]),
            n_absent_used = sum(absent & frq > 1e-6),
            n_absent      = sum(absent),
            top_absent    = if (any(absent)) max(frq[absent]) else 0,
            .groups = "drop")

leak <- leak_sample %>%
  group_by(ref_size, fraction) %>%
  summarise(n_samples      = n(),
            leak_median    = median(leak),
            leak_mean      = mean(leak),
            leak_min       = min(leak),
            leak_max       = max(leak),
            n_absent       = first(n_absent),
            used_median    = median(n_absent_used),
            top_absent_max = max(top_absent),
            .groups = "drop")
write_tsv(leak, file.path(OUT, "TABLE_downsample_leakage.tsv"))

cat("\n== off-pool leakage, median across samples (% of the pool) ==\n")
print(as.data.frame(leak %>%
  transmute(fraction, ref_size, pct = round(100 * leak_median, 2)) %>%
  pivot_wider(names_from = ref_size, values_from = pct, names_prefix = "R=")),
  row.names = FALSE)

cat("\n== absent strains given non-trivial frequency (>1e-6), median sample ==\n")
print(as.data.frame(leak %>%
  transmute(fraction, ref_size, used = used_median) %>%
  pivot_wider(names_from = ref_size, values_from = used, names_prefix = "R=")),
  row.names = FALSE)

## ---------------------------------------------------------------------------
## 3. IBS from every absent strain to its closest pool member
## ---------------------------------------------------------------------------
if (!REFRESH && file.exists(IBSC)) {
  ibs <- read_tsv(IBSC, show_col_types = FALSE)
  say("IBS read from cache (", nrow(ibs), " strains); --refresh-ibs to rebuild")
} else {
  for (f in c(BOOT, GENO))
    if (!file.exists(f)) stop("need ", f, " to build the IBS cache -- see ",
                              "DATA_AVAILABILITY.md", call. = FALSE)
  say("loading the archive for its marker set")
  e <- new.env(); load(BOOT, e)
  markers <- rownames(e$flipped_bootstrap_input[[1]])
  rm(e); invisible(gc())
  say("  ", length(markers), " archive markers")

  say("loading the CeNDR genotype matrix (12 GB, ~25 s)")
  g2 <- new.env(); load(GENO, g2)
  cn <- sub("_.*$", "", colnames(g2$g))
  want <- sort(unique(pred$strain))
  stopifnot(all(want %in% cn))
  mi <- match(intersect(markers, rownames(g2$g)), rownames(g2$g))
  say("  archive markers present in CeNDR: ", length(mi))

  set.seed(SEED)
  if (length(mi) > N_MARKERS) {
    mi <- sort(sample(mi, N_MARKERS))
    say("  sampled ", length(mi), " markers under seed ", SEED)
  }
  GT <- g2$g[mi, match(want, cn), drop = FALSE]
  rm(g2); invisible(gc())
  colnames(GT) <- want

  keep <- which(rowSums(is.na(GT)) == 0)
  GT <- GT[keep, , drop = FALSE]
  say("  complete markers: ", nrow(GT))
  vals <- unique(as.vector(GT[seq_len(min(5000, nrow(GT))), , drop = FALSE]))
  if (!all(vals %in% c(0, 1)))
    stop("genotypes are not 0/1 coded -- found ", paste(sort(vals), collapse = ", "),
         ".\n  The crossprod IBS identity below assumes 0/1.", call. = FALSE)
  ## a 0/1 marker varies iff its row sum is neither 0 nor the strain count;
  ## invariant sites inflate every IBS by the same amount and carry no
  ## information about who confounds whom
  rs <- rowSums(GT)
  GT <- GT[rs > 0 & rs < ncol(GT), , drop = FALSE]
  n <- nrow(GT)
  say("  variable markers used: ", n)
  if (n < 1e4) stop("only ", n, " variable markers -- too few for a stable IBS.",
                    call. = FALSE)

  ## mean(x == y) = (both-one + both-zero) / n, the same identity and therefore
  ## the same scale as baugh_strain_similarity.R
  cp <- crossprod(GT)
  cs <- colSums(GT)
  same <- (cp + (n - outer(cs, cs, "+") + cp)) / n
  diag(same) <- NA
  rm(GT, cp); invisible(gc())

  ## nearest POOL member, for every strain. For a pool strain this excludes
  ## itself (the NA diagonal), so its value is comparable with the published
  ## nn_ibs among the 102.
  M <- same[, POOL, drop = FALSE]
  j <- apply(M, 1, function(r) if (all(is.na(r))) NA_integer_ else which.max(r))
  ibs <- tibble(strain      = rownames(M),
                nn_pool_ibs = vapply(seq_len(nrow(M)), function(i)
                                if (is.na(j[i])) NA_real_ else M[i, j[i]], 0),
                nn_pool     = ifelse(is.na(j), NA_character_, colnames(M)[j]),
                mean_pool_ibs = rowMeans(M, na.rm = TRUE))
  rm(same, M); invisible(gc())
  write_tsv(ibs, IBSC)
  say("IBS cached to ", IBSC)
}

stopifnot(all(unique(pred$strain) %in% ibs$strain),
          all(ibs$nn_pool_ibs > 0.4, na.rm = TRUE),
          all(ibs$nn_pool_ibs <= 1, na.rm = TRUE))

## Does the flip-invariance argument hold? Pool strains' nearest-pool-member IBS
## should track the published table computed on the oriented matrix.
pub <- read_tsv(SIM, show_col_types = FALSE)
chk <- ibs %>% filter(strain %in% POOL) %>% inner_join(pub, by = "strain")
r_chk <- cor(chk$nn_pool_ibs, chk$nn_ibs)
say("check vs published nn_ibs on the ", nrow(chk), " pool strains: r = ",
    sprintf("%.4f", r_chk), ", median |diff| = ",
    sprintf("%.4f", median(abs(chk$nn_pool_ibs - chk$nn_ibs))))
if (r_chk < 0.95)
  stop("IBS computed here does not track the published table (r = ",
       signif(r_chk, 3), "). The marker set or the coding differs; fix that ",
       "before reading the relatedness panels.", call. = FALSE)

## ---------------------------------------------------------------------------
## 4. is leakage directed by relatedness to the pool?
## ---------------------------------------------------------------------------
absent_leak <- pred %>%
  filter(!strain %in% POOL) %>%
  group_by(ref_size, fraction, strain) %>%
  summarise(leak = mean(frq), .groups = "drop") %>%
  inner_join(ibs %>% select(strain, nn_pool_ibs, nn_pool), by = "strain")

## A RANK CORRELATION ON THE WHOLE SET IS THE WRONG STATISTIC HERE, and reading
## it alone was an error worth leaving a marker against. NNLS is non-negative,
## so its solutions are sparse: at 540 candidates 56.8% of the absent columns
## get EXACTLY zero, and those ties drag any Spearman toward nothing. The
## response has to be decomposed --
##
##   (i)  does IBS predict WHETHER a candidate enters the support at all? No.
##        Logistic P(absorbs > 0) ~ IBS has slope -1.13, p = 0.71, and the rank
##        correlation with the indicator is -0.021. Which columns the solver
##        admits is a property of the collinearity of the whole design, not of
##        one pairwise distance, and nothing here explains it.
##   (ii) among those that do absorb, does IBS predict HOW MUCH? Yes, and
##        consistently: rho +0.46, +0.39, +0.34, +0.35 at 150, 250, 400, 540.
##
## The whole-set rho falling from +0.304 to +0.051 across reference sizes is
## therefore an artifact of the zero fraction climbing from 25% to 57%, not a
## relationship that weakens. Both parts are computed below so the mixture
## cannot be misread again.
leak_ibs <- absent_leak %>%
  group_by(ref_size, fraction) %>%
  summarise(n = n(),
            zero_frac = mean(leak == 0),
            rho = suppressWarnings(cor(nn_pool_ibs, leak, method = "spearman")),
            p   = suppressWarnings(cor.test(nn_pool_ibs, leak,
                                            method = "spearman")$p.value),
            rho_nonzero = if (sum(leak > 0) > 3)
              suppressWarnings(cor(nn_pool_ibs[leak > 0], leak[leak > 0],
                                   method = "spearman")) else NA_real_,
            rho_support = suppressWarnings(cor(nn_pool_ibs,
                                               as.numeric(leak > 0),
                                               method = "spearman")),
            .groups = "drop")
write_tsv(leak_ibs, file.path(OUT, "TABLE_downsample_leakage_ibs.tsv"))

cat("\n== Spearman of leakage against IBS, whole set / among absorbers ==\n")
print(as.data.frame(leak_ibs %>%
  transmute(fraction, ref_size, v = sprintf("%+.3f / %+.3f", rho, rho_nonzero)) %>%
  pivot_wider(names_from = ref_size, values_from = v, names_prefix = "R=")),
  row.names = FALSE)
cat("\n== fraction of absent candidates absorbing EXACTLY zero ==\n")
print(as.data.frame(leak_ibs %>%
  transmute(fraction, ref_size, z = round(100 * zero_frac, 1)) %>%
  pivot_wider(names_from = ref_size, values_from = z, names_prefix = "R=")),
  row.names = FALSE)

worst <- absent_leak %>%
  filter(ref_size == max(REF_SIZES), fraction == max(FRACTIONS)) %>%
  arrange(desc(leak))
write_tsv(worst, file.path(OUT, "TABLE_absent_strain_leakage.tsv"))
cat("\n== the ten absent strains that absorb the most, R=", max(REF_SIZES),
    " at full depth ==\n", sep = "")
print(as.data.frame(worst %>% slice_head(n = 10) %>%
  transmute(strain, leak_permille = round(1000 * leak, 2),
            nn_pool_ibs = round(nn_pool_ibs, 4), nn_pool)), row.names = FALSE)

## ---------------------------------------------------------------------------
## 5. does it still matter once the reference is restricted to the true pool?
##
## The renormalised series keeps only pool columns and rescales them to sum to
## 1, which is what an experiment with a KNOWN membership list gets. If the
## renormalised slopes recover the accuracy of the R=102 reference, then leakage
## costs precision only when membership is unknown, and the answer to "so what
## is the big deal" is: the big deal is bounded, and bounded by this much.
## ---------------------------------------------------------------------------
ols <- function(x, y) {
  o <- is.finite(x) & is.finite(y); x <- x[o]; y <- y[o]
  if (length(x) < 2) return(NA_real_)
  xc <- x - mean(x); d <- sum(xc^2)
  if (d == 0) return(NA_real_)
  sum(xc * y) / d
}
## identical to Figure 1 and to diagnostic_reference_size.R
slopes_for <- function(d) {
  f <- d %>% left_join(meta, by = "sample") %>% filter(!baseline, day != 17)
  base <- f %>% filter(day == 1) %>% select(replicate, strain, b = frq)
  f %>% left_join(base, by = c("replicate", "strain")) %>%
    mutate(delta = frq - b) %>%
    group_by(replicate, strain) %>%
    summarise(slope = ols(day, delta), .groups = "drop") %>%
    group_by(strain) %>% summarise(slope = mean(slope, na.rm = TRUE), .groups = "drop")
}
mip_slopes <- slopes_for(arch %>% transmute(strain, sample, frq = published_frq)) %>%
  rename(mip = slope)

score <- function(d) {
  s <- slopes_for(d) %>% inner_join(mip_slopes, by = "strain") %>%
    filter(strain != "N2", is.finite(slope), is.finite(mip))
  tibble(n = nrow(s), rho = cor(s$slope, s$mip, method = "spearman"),
         rmse = sqrt(mean((s$slope - s$mip)^2)))
}

say("scoring slopes, as fitted and after restricting to the true pool")
slp <- pred %>% group_by(ref_size, fraction) %>% group_split() %>%
  map_dfr(function(d) {
    R <- d$ref_size[1]; fr <- d$fraction[1]
    asis <- score(d %>% select(strain, sample, frq))
    rn <- d %>% filter(strain %in% POOL) %>%
      group_by(sample) %>% mutate(frq = frq / sum(frq)) %>% ungroup()
    ren <- score(rn %>% select(strain, sample, frq))
    tibble(ref_size = R, fraction = fr,
           rho_asis = asis$rho, rho_renorm = ren$rho,
           rmse_asis = asis$rmse, rmse_renorm = ren$rmse)
  })
write_tsv(slp, file.path(OUT, "TABLE_downsample_leakage_slopes.tsv"))

cat("\n== slope agreement with MIP-seq: as fitted / restricted to the pool ==\n")
print(as.data.frame(slp %>%
  transmute(fraction, ref_size,
            v = sprintf("%.3f / %.3f", rho_asis, rho_renorm)) %>%
  pivot_wider(names_from = ref_size, values_from = v, names_prefix = "R=")),
  row.names = FALSE)

base_rho <- slp %>% filter(ref_size == min(REF_SIZES)) %>% select(fraction, ref_rho = rho_asis)
gap <- slp %>% inner_join(base_rho, by = "fraction") %>%
  filter(ref_size == max(REF_SIZES)) %>%
  transmute(fraction,
            cost_unknown = ref_rho - rho_asis,
            cost_known   = ref_rho - rho_renorm)
cat("\n== cost of not knowing the membership, R=", max(REF_SIZES),
    " vs R=", min(REF_SIZES), " (loss in rho) ==\n", sep = "")
print(as.data.frame(gap %>% mutate(across(starts_with("cost"), ~round(.x, 4)))),
      row.names = FALSE)

## ---------------------------------------------------------------------------
## 6. the figure
## ---------------------------------------------------------------------------
pt <- function(l) paste0("<span style='font-size:13pt;color:#111111'>**", l, "**</span>")
th <- theme_classic(base_size = 11) +
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_markdown(size = 11.5),
        plot.title.position = "plot",
        legend.key.size = grid::unit(9, "pt"))
RCOL <- setNames(colorRampPalette(c("#1B6C7A", "#7A9A3B", "#C08A2E", "#B5623C",
                                    "#9E4257"))(length(REF_SIZES)),
                 as.character(REF_SIZES))
lab_x <- scale_x_log10(breaks = FRACTIONS,
                       labels = c("1/64","1/32","1/16","1/8","1/4","1/2","1"))

## A -- leakage against coverage. R=102 is kept in deliberately: it has no
## absent columns, so it is a flat zero and is the floor the others are read
## against, and keeping it means this one legend also names panel D's lines.
pA <- ggplot(leak, aes(fraction, 100 * leak_median, colour = factor(ref_size))) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.7) +
  lab_x + scale_colour_manual(values = RCOL, name = "Reference\nstrains") +
  labs(x = "Fraction of full-depth reads",
       y = "Off-pool leakage (% of the pool)", title = pt("A")) + th

## B -- where it goes. Decile MEANS, not medians: the distribution is
## zero-inflated -- most absent candidates absorb nothing at all -- so every
## decile median is ~0 and a median line draws a flat trace through data that
## plainly has structure. The mean is the quantity that actually adds up to the
## leakage in panel A, which is what this panel is decomposing. A loess is wrong
## here for the usual reason: the response is bounded below by zero and the
## smoother dips through it.
sub <- absent_leak %>% filter(ref_size == max(REF_SIZES), fraction == max(FRACTIONS))
bins <- sub %>%
  mutate(b = cut(nn_pool_ibs, breaks = quantile(nn_pool_ibs, seq(0, 1, 0.1)),
                 include.lowest = TRUE)) %>%
  group_by(b) %>% summarise(x = median(nn_pool_ibs), y = mean(1000 * leak),
                            share = mean(leak > 5e-4), .groups = "drop")
cat("\n== absorbed per mille by IBS decile, R=", max(REF_SIZES),
    " at full depth ==\n", sep = "")
print(as.data.frame(bins %>% transmute(ibs_median = round(x, 4),
        mean_permille = round(y, 3),
        share_over_half_permille = round(share, 3))), row.names = FALSE)

pB <- ggplot(sub, aes(nn_pool_ibs, 1000 * leak)) +
  geom_point(size = 0.9, alpha = 0.35, colour = "grey40") +
  geom_line(data = bins, aes(x, y), linewidth = 0.7,
            colour = RCOL[[as.character(max(REF_SIZES))]]) +
  geom_point(data = bins, aes(x, y), size = 1.9,
             colour = RCOL[[as.character(max(REF_SIZES))]]) +
  labs(x = "IBS to the closest pool member",
       y = "Frequency absorbed (per mille)",
       subtitle = "points, candidates; line, decile mean",
       title = pt("B")) + th +
  theme(plot.subtitle = element_text(size = 8.5, colour = "grey35"))

## C -- and whether coverage changes that direction. BOTH series are drawn on
## purpose: the whole-set correlation alone reads as a relationship dissolving
## as the reference grows, and it is not -- it is the exact-zero fraction
## climbing. Putting the among-absorbers series beside it means the panel cannot
## be misread without also ignoring the dashed lines.
pC <- ggplot(leak_ibs %>% filter(ref_size > min(REF_SIZES)),
             aes(fraction, colour = factor(ref_size))) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_line(aes(y = rho, linetype = "all absent candidates"), linewidth = 0.6) +
  geom_point(aes(y = rho), size = 1.7) +
  geom_line(aes(y = rho_nonzero, linetype = "those absorbing > 0"),
            linewidth = 0.6) +
  lab_x + scale_colour_manual(values = RCOL, guide = "none") +
  scale_linetype_manual(values = c("all absent candidates" = "solid",
                                   "those absorbing > 0" = "22"), name = NULL) +
  labs(x = "Fraction of full-depth reads",
       y = "ρ, IBS vs frequency absorbed", title = pt("C")) + th +
  theme(legend.position = c(0.66, 0.90),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA))

## D -- the part that decides whether any of it matters
pD <- ggplot(slp, aes(fraction, colour = factor(ref_size))) +
  geom_line(aes(y = rho_asis, linetype = "reference as given"), linewidth = 0.6) +
  geom_point(aes(y = rho_asis), size = 1.5) +
  geom_line(aes(y = rho_renorm, linetype = "restricted to the pool"), linewidth = 0.6) +
  lab_x + scale_colour_manual(values = RCOL, guide = "none") +
  scale_linetype_manual(values = c("reference as given" = "solid",
                                   "restricted to the pool" = "22"), name = NULL) +
  labs(x = "Fraction of full-depth reads", y = "Spearman ρ, slopes vs MIP-seq",
       title = pt("D")) + th + theme(legend.position = c(0.68, 0.22))

fig <- (pA | pB) / (pC | pD)
ggsave(file.path(OUT, "downsample_leakage.pdf"), fig,
       width = 10, height = 7, device = cairo_pdf)
ggsave(file.path(OUT, "downsample_leakage.png"), fig,
       width = 10, height = 7, dpi = 300, bg = "white")
say("wrote downsample_leakage.{pdf,png} and three tables")
