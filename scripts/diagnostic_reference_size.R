## Diagnostic -- slope accuracy against coverage AND reference size ----------
##
##   Rscript scripts/diagnostic_reference_size.R
##     -> plots/diagnostics/reference_size_vs_coverage.{pdf,png}
##     -> plots/diagnostics/TABLE_reference_size_vs_coverage.tsv
##
## NEEDS THE ARCHIVE: data/baugh/2024bootstrapINPUT.Rdata and
## data/genotypes/processed_genotype_matrix.Rda, both Dryad-hosted.
##
## THE QUESTION. In a real pooled experiment you do not know which strains are
## in the pool, so the NNLS reference carries candidates that are absent. Every
## extra candidate is another column the solver can misassign frequency to, and
## columns are near-collinear in proportion to relatedness. So: how does slope
## recovery degrade as coverage falls, and how much worse does it get as the
## reference admits more candidate strains?
##
## THE DESIGN. The 102 pooled strains are always in the reference; larger
## references add strains drawn from the other 438 CeNDR isotypes under a fixed
## seed. Coverage is varied by BINOMIAL THINNING of the archived alt counts,
## which is exactly read-level thinning -- if reads are kept independently with
## probability p, the alt reads retained are Binomial(alt, p) -- and needs no
## knowledge of the reference-allele counts, which the archive does not hold.
##
## DEPTH IS THEREFORE RELATIVE, not absolute. The archive stores alt counts
## only, so the absolute coverage of the full data cannot be recovered and the
## 0.25x-10x labels of baugh_downsampled_slopes.rda cannot be reproduced. Depth
## here is the fraction of full-depth reads retained, and the mean alt count per
## marker is reported alongside it so the two can be related later if the
## reference-allele counts are ever archived.
##
## ORIENTATION. The archive's genotypes are FLIPPED: markers whose alternate
## allele is common in the 102-strain panel were recoded so the counted allele
## is the minor one, and the counts were recoded to match. Strains added to the
## reference must be recoded the same way or the design matrix is inconsistent
## with the counts. The mask is recovered per marker by comparing the archive
## against the CeNDR matrix on the 102 shared strains: agreement 1 means not
## flipped, agreement 0 means flipped, and a marker that agrees partly means the
## two matrices disagree on genotype calls and is dropped rather than guessed.
## On a 20,000-marker probe this was 19,883 identical, 117 complemented and
## ZERO ambiguous, so the mask is exact.
##
## FIDELITY CHECK. At reference size 102 and full depth the pipeline must
## reproduce the archived NNLS frequencies. It is checked, not assumed, and the
## script stops if it fails -- without it, any degradation seen here could be a
## bug in this script rather than a property of the deconvolution.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse); library(RcppML); library(patchwork); library(ggtext)
})

OUT   <- "plots/diagnostics"
BOOT  <- "data/baugh/2024bootstrapINPUT.Rdata"
GENO  <- "data/genotypes/processed_genotype_matrix.Rda"
NNLS  <- "supplemental_data/deconvolution/baugh_nnls_with_mipseq.RData"
SIM   <- "supplemental_data/deconvolution/baugh_strain_similarity.tsv"
for (f in c(BOOT, GENO, NNLS, SIM))
  if (!file.exists(f)) stop("need ", f, " -- see DATA_AVAILABILITY.md", call. = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

REF_SIZES <- c(102, 150, 250, 400, 540)
FRACTIONS <- c(1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1)
SEED      <- 1
say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep=""); flush.console() }

## ---------------------------------------------------------------------------
## 1. the archive
## ---------------------------------------------------------------------------
say("loading the Baugh design matrix and counts")
e <- new.env(); load(BOOT, e)
gt0     <- e$flipped_bootstrap_input[[1]]
counts0 <- e$flipped_bootstrap_input[[2]]
rm(e); invisible(gc())
pool <- sub("_.*$", "", colnames(gt0))
markers <- rownames(gt0)
say("  ", nrow(gt0), " markers x ", ncol(gt0), " pooled strains, ",
    ncol(counts0), " samples")
stopifnot(!anyNA(gt0), !anyNA(counts0), all(counts0 >= 0))

## ---------------------------------------------------------------------------
## 2. orientation mask and the enlarged, correctly oriented genotype matrix
## ---------------------------------------------------------------------------
say("loading the CeNDR genotype matrix (12 GB)")
g2 <- new.env(); load(GENO, g2)
cn  <- sub("_.*$", "", colnames(g2$g))
stopifnot(all(pool %in% cn), all(markers %in% rownames(g2$g)))
mi  <- match(markers, rownames(g2$g))
pi_ <- match(pool, cn)

say("recovering the flip mask on all ", length(markers), " markers")
A <- gt0
B <- g2$g[mi, pi_, drop = FALSE]
agree <- rowMeans(A == B)
flip  <- agree == 0
ok    <- agree == 1 | flip
say("  identical ", sum(agree == 1), " | complemented ", sum(flip),
    " | ambiguous ", sum(!ok), " (dropped)")
rm(A, B); invisible(gc())

extra <- setdiff(cn, pool)
set.seed(SEED)
extra_order <- sample(extra)               # one fixed order; nested references
keep_strains <- c(pool, extra_order[seq_len(max(REF_SIZES) - length(pool))])
say("building the oriented reference for ", length(keep_strains), " strains")
GT <- g2$g[mi[ok], match(keep_strains, cn), drop = FALSE]
rm(g2); invisible(gc())
GT[flip[ok], ] <- 1 - GT[flip[ok], ]        # same recoding the archive used
colnames(GT) <- keep_strains
gt0 <- gt0[ok, , drop = FALSE]
counts0 <- counts0[ok, , drop = FALSE]
## Values only: identical() would also compare dimnames, and GT carries bare
## strain names where the archive carries "STRAIN_STRAIN".
n_mismatch <- sum(GT[, seq_along(pool)] != gt0)
say("  oriented reference vs archive on the pooled columns: ", n_mismatch,
    " mismatching cells of ", format(length(gt0), big.mark = ","))
if (n_mismatch != 0)
  stop("the flip mask did not reproduce the archive's genotypes -- ",
       n_mismatch, " cells differ. Do not trust anything downstream.",
       call. = FALSE)
rm(gt0); invisible(gc())
say("  ", nrow(GT), " markers x ", ncol(GT), " candidate strains")

## ---------------------------------------------------------------------------
## 3. thinned count matrices, generated once so every reference size sees the
##    same reads at a given depth
## ---------------------------------------------------------------------------
say("thinning counts")
set.seed(SEED)
thinned <- lapply(FRACTIONS, function(p) {
  if (p == 1) return(counts0)
  m <- matrix(rbinom(length(counts0), as.vector(counts0), p),
              nrow = nrow(counts0), dimnames = dimnames(counts0))
  m
})
names(thinned) <- sprintf("%.5f", FRACTIONS)
alt_per_marker <- sapply(thinned, function(m) mean(colSums(m)) / nrow(m))
say("  mean alt reads per marker: ",
    paste(sprintf("%.3f", alt_per_marker), collapse = ", "))

## ---------------------------------------------------------------------------
## 4. deconvolve, for every reference size x depth
## ---------------------------------------------------------------------------
solve_freq <- function(gt, ct) {
  GGp <- crossprod(gt)
  Gy  <- crossprod(gt, ct)
  p   <- apply(Gy, 2, function(x) as.vector(RcppML::nnls(GGp, matrix(x), fast_nnls = TRUE)))
  p   <- sweep(p, 2, colSums(p), "/")
  dimnames(p) <- list(colnames(gt), colnames(ct))
  p
}

res <- list()
for (R in REF_SIZES) {
  idx <- seq_len(R)
  gtR <- GT[, idx, drop = FALSE]
  say("reference size ", R, ": crossprod")
  for (k in seq_along(FRACTIONS)) {
    fr <- FRACTIONS[k]
    P <- solve_freq(gtR, thinned[[k]])
    ## built explicitly rather than via as_tibble(as.table(P)): that idiom
    ## needs NAMED dimnames and silently is not what a matrix carries
    res[[length(res) + 1L]] <- tibble(
      strain   = rep(rownames(P), times = ncol(P)),
      sample   = rep(colnames(P), each  = nrow(P)),
      frq      = as.vector(P),
      ref_size = R, fraction = fr)
  }
  say("  done ", R)
  rm(gtR); invisible(gc())
}
pred <- bind_rows(res)
saveRDS(pred, file.path(OUT, "reference_size_predictions.rds"))

## ---------------------------------------------------------------------------
## 5. FIDELITY: reference 102 at full depth must reproduce the archive
## ---------------------------------------------------------------------------
a <- new.env(); load(NNLS, a)
arch <- as_tibble(a$wgs_mip_results) %>% select(strain, sample, frq, published_frq)
chk <- pred %>% filter(ref_size == 102, fraction == 1) %>%
  inner_join(arch, by = c("strain", "sample"), suffix = c("_new", "_arch"))
mx <- max(abs(chk$frq_new - chk$frq), na.rm = TRUE)
say("fidelity at 102 strains, full depth: ", nrow(chk),
    " strain-samples, max |difference| ", signif(mx, 3))
if (!(mx < 1e-6))
  stop("pipeline does not reproduce the archived frequencies (max diff ",
       signif(mx, 3), "). Fix before reading anything below.", call. = FALSE)

## ---------------------------------------------------------------------------
## 6. slopes, exactly as Figure 1 fits them
## ---------------------------------------------------------------------------
ols <- function(x, y) {
  o <- is.finite(x) & is.finite(y); x <- x[o]; y <- y[o]
  if (length(x) < 2) return(NA_real_)
  xc <- x - mean(x); d <- sum(xc^2)
  if (d == 0) return(NA_real_)
  sum(xc * y) / d
}
meta <- tibble(sample = colnames(counts0)) %>%
  mutate(baseline  = grepl("_baseline$", sample),
         replicate = sub("^(rep[0-9]+)_.*$", "\\1", sample),
         day       = as.integer(sub("^rep[0-9]+_d([0-9]+).*$", "\\1", sample)))
stopifnot(!anyNA(meta$day), !anyNA(meta$replicate))

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

grid <- pred %>% group_by(ref_size, fraction) %>% group_split()
acc <- map_dfr(grid, function(d) {
  R <- d$ref_size[1]; fr <- d$fraction[1]
  s <- slopes_for(d %>% select(strain, sample, frq)) %>%
    inner_join(mip_slopes, by = "strain") %>%
    filter(strain != "N2", is.finite(slope), is.finite(mip))
  tibble(ref_size = R, fraction = fr, n = nrow(s),
         rho = cor(s$slope, s$mip, method = "spearman"),
         rmse = sqrt(mean((s$slope - s$mip)^2)))
})
acc <- acc %>% mutate(alt_per_marker = alt_per_marker[sprintf("%.5f", fraction)])
write_tsv(acc, file.path(OUT, "TABLE_reference_size_vs_coverage.tsv"))

cat("\n== Spearman of NNLS against MIP-seq slopes ==\n")
print(as.data.frame(acc %>% select(ref_size, fraction, rho) %>%
  mutate(rho = round(rho, 3)) %>%
  pivot_wider(names_from = ref_size, values_from = rho,
              names_prefix = "R=")), row.names = FALSE)

## ---------------------------------------------------------------------------
## 7. does relatedness bite harder at low coverage?
## ---------------------------------------------------------------------------
sim <- read_tsv(SIM, show_col_types = FALSE) %>%
  transmute(strain, nn_ibs = coalesce(nn_ibs_wild, nn_ibs))
disc <- map_dfr(grid, function(d) {
  R <- d$ref_size[1]; fr <- d$fraction[1]
  x <- d %>% select(strain, sample, frq) %>%
    inner_join(arch %>% select(strain, sample, mip = published_frq),
               by = c("strain", "sample")) %>%
    filter(strain != "N2", is.finite(mip)) %>%
    group_by(strain) %>% summarise(d_abs = mean(abs(frq - mip)), .groups = "drop") %>%
    inner_join(sim, by = "strain")
  ct <- suppressWarnings(cor.test(x$nn_ibs, x$d_abs, method = "spearman"))
  tibble(ref_size = R, fraction = fr, n = nrow(x),
         rho_ibs = unname(ct$estimate), p = ct$p.value,
         median_gap_permille = 1000 * median(x$d_abs))
})
write_tsv(disc, file.path(OUT, "TABLE_reference_size_leakage.tsv"))
cat("\n== relatedness against discrepancy, by coverage and reference size ==\n")
print(as.data.frame(disc %>% select(ref_size, fraction, rho_ibs) %>%
  mutate(rho_ibs = round(rho_ibs, 3)) %>%
  pivot_wider(names_from = ref_size, values_from = rho_ibs,
              names_prefix = "R=")), row.names = FALSE)
cat("\n== median discrepancy, per mille of pool frequency ==\n")
print(as.data.frame(disc %>% select(ref_size, fraction, median_gap_permille) %>%
  mutate(median_gap_permille = round(median_gap_permille, 2)) %>%
  pivot_wider(names_from = ref_size, values_from = median_gap_permille,
              names_prefix = "R=")), row.names = FALSE)

## ---------------------------------------------------------------------------
## 8. the figure
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

pA <- ggplot(acc, aes(fraction, rho, colour = factor(ref_size))) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.7) +
  lab_x + scale_colour_manual(values = RCOL, name = "Reference\nstrains") +
  labs(x = "Fraction of full-depth reads", y = "Spearman ρ, slopes vs MIP-seq",
       title = pt("A")) + th
pB <- ggplot(acc, aes(fraction, rmse, colour = factor(ref_size))) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.7) +
  lab_x + scale_colour_manual(values = RCOL, guide = "none") +
  labs(x = "Fraction of full-depth reads", y = "RMSE of the slope",
       title = pt("B")) + th
pC <- ggplot(disc, aes(fraction, median_gap_permille, colour = factor(ref_size))) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.7) +
  lab_x + scale_y_log10() + scale_colour_manual(values = RCOL, guide = "none") +
  labs(x = "Fraction of full-depth reads",
       y = "Median |NNLS − MIP| (per mille)", title = pt("C")) + th
pD <- ggplot(disc, aes(fraction, rho_ibs, colour = factor(ref_size))) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_line(linewidth = 0.6) + geom_point(size = 1.7) +
  lab_x + scale_colour_manual(values = RCOL, guide = "none") +
  labs(x = "Fraction of full-depth reads",
       y = "ρ, relatedness vs discrepancy", title = pt("D")) + th

fig <- (pA | pB) / (pC | pD)
ggsave(file.path(OUT, "reference_size_vs_coverage.pdf"), fig,
       width = 10, height = 7, device = cairo_pdf)
ggsave(file.path(OUT, "reference_size_vs_coverage.png"), fig,
       width = 10, height = 7, dpi = 300, bg = "white")
say("wrote reference_size_vs_coverage.{pdf,png} and two tables to ", OUT)
