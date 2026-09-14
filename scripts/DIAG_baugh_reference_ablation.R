## Does adding the three non-pool strains really improve agreement, or is any --
## extra column a sponge? -----------------------------------------------------
##
## The objection this answers: adding columns to a least-squares problem almost
## always reduces residual, so "the fit got better" is not evidence that the
## extra strains are real. The guard is that agreement is scored against MIP-seq,
## an INDEPENDENT measurement that the fit never sees. But that guard is only
## convincing with a negative control, so this script fits four references:
##
##   FULL         all 103 strains, the reference as shipped
##   RESTRICTED   the 100 strains the MIP panel measured
##   SYNTH        the 100, plus 3 SYNTHETIC near-twins built by perturbing an
##                existing pool strain to the same IBS the real extras have to
##                their twins. A synthetic twin carries no new information --
##                it is a sponge by construction. If it helps as much as the
##                real extras, the improvement is generic.
##   DROP3        97 strains: the 100 minus the three twins CX11264, PS2025 and
##                NIC256, to show what losing a genuinely present strain costs
##                on the same scale.
##
## Each fit is scored against MIP per strain, and the per-strain change is
## plotted against IBS to the three extras.
##
## REQUIRES a file outside this repository:
##   /Users/Stefan/UCLA/Projects/bulkGWAS/baugh_wgs/cluster_data/
##     20220908_Baugh_BulkL1_Bootstrap_Input_flippedCommon_NAfix.RData
##
## WHAT IT FINDS
##
## The sponge hypothesis is refuted, but the effect is far narrower than the
## aggregate suggests, and both halves matter.
##
## Refuted: the IBS-matched synthetic twins absorb 0.094% of the pool mass
## against the real extras' 3.21%, a 34-fold difference, and the SYNTH fit is
## indistinguishable from RESTRICTED (per-strain RMSD correlation 0.99991, max
## absolute difference 4.1e-4). NNLS does not hand mass to any near-twin column
## on offer; matching a strain's IBS is not enough, the column has to match the
## actual allele counts. So the gain from the real extras is specific to them.
##
## Narrow: the penalty for removing them is not diffuse and is not a function
## of relatedness. Spearman rho between a strain's max IBS to an extra and its
## RMSD penalty is 0.124 (p = 0.22), and 0.043 (p = 0.68) once the three twins
## are set aside. PS2025 alone carries 41% of the total penalty and the three
## twins carry 58%. Excluding those three, mean RMSD rises only 2.9% when the
## extras are dropped, against 6.5% including them, and 45 of 99 strains are
## actually BETTER without them.
##
## So "including the extras improves agreement" is true in aggregate but is
## carried by one strain pair (ECA348/PS2025, IBS 0.984) and secondarily by
## NIC260/NIC256. It is not a general effect of adding related genotypes.
##
## Writes plots/diagnostics/DIAG_baugh_reference_ablation.{pdf,png} and
## baugh_reference_ablation.tsv. Nothing else is modified.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(ggrepel)
})

GT   <- "/Users/Stefan/UCLA/Projects/bulkGWAS/baugh_wgs/cluster_data/20220908_Baugh_BulkL1_Bootstrap_Input_flippedCommon_NAfix.RData"
MIPF <- "supplemental_data/deconvolution/mipseq_frequencies.txt.gz"
DIAG <- "plots/diagnostics"
EXTRA <- c("CX11262", "ECA348", "NIC260")
TWIN  <- c(CX11262 = "CX11264", ECA348 = "PS2025", NIC260 = "NIC256")
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)

e <- new.env(); load(GT, e)
gt <- e$flipped_bootstrap_input[[1]]; ct <- e$flipped_bootstrap_input[[2]]
keep <- which(rowSums(is.na(gt)) == 0); gt <- gt[keep, ]; ct <- ct[keep, ]
strains <- sub("_.*$", "", colnames(gt))
colnames(gt) <- strains
message(sprintf("markers %d, strains %d", nrow(gt), ncol(gt)))

## --- MIP reference ---------------------------------------------------------
ln <- readLines(gzfile(MIPF)); hdr <- strsplit(ln[1], "\t")[[1]]
body <- strsplit(ln[-1], "\t")
mipf <- tibble(strain = sapply(body, `[`, 1),
  !!!setNames(lapply(seq_along(hdr),
      function(j) as.numeric(sapply(body, `[`, j + 1))), hdr)) %>%
  pivot_longer(-strain, names_to = "sample", values_to = "mip") %>%
  mutate(sample = gsub("_BL", "_d1_baseline", sample)) %>%
  filter(is.finite(mip))
pool <- unique(mipf$strain)

## --- IBS of every strain to the three extras -------------------------------
ibs <- function(a, b) mean(gt[, a] == gt[, b])
ibs_tab <- expand_grid(strain = setdiff(strains, EXTRA), extra = EXTRA) %>%
  mutate(ibs = map2_dbl(strain, extra, ibs)) %>%
  group_by(strain) %>%
  summarise(max_ibs = max(ibs), nearest_extra = extra[which.max(ibs)],
            .groups = "drop")

## --- synthetic sponge twins, IBS-matched to the real extras ----------------
real_ibs <- map_dbl(EXTRA, ~ ibs(.x, TWIN[[.x]]))
names(real_ibs) <- EXTRA
message("real extra-to-twin IBS: ",
        paste(sprintf("%s/%s %.3f", EXTRA, TWIN[EXTRA], real_ibs), collapse = "; "))
set.seed(11)
synth <- map(EXTRA, function(x) {
  src <- gt[, TWIN[[x]]]
  n_flip <- round((1 - real_ibs[[x]]) * length(src))
  idx <- sample.int(length(src), n_flip)
  src[idx] <- 1L - src[idx]
  src
})
synth <- do.call(cbind, synth)
colnames(synth) <- paste0("SYNTH_", TWIN[EXTRA])
message("synthetic IBS to source: ",
        paste(sprintf("%.3f", map_dbl(seq_len(3),
              ~ mean(synth[, .x] == gt[, TWIN[EXTRA][.x]]))), collapse = " "))

## --- fits ------------------------------------------------------------------
fit <- function(mat, tag) {
  Gy <- crossprod(mat, ct); GG <- crossprod(mat)
  p <- apply(Gy, 2, function(x) as.vector(RcppML::nnls(GG, matrix(x), fast_nnls = TRUE)))
  p <- apply(p, 2, function(x) x / sum(x)); rownames(p) <- colnames(mat)
  message(sprintf("  %-11s fitted %d strains", tag, nrow(p)))
  as_tibble(p, rownames = "strain") %>%
    pivot_longer(-strain, names_to = "sample", values_to = "frq")
}
pool_idx <- which(strains %in% pool)
fits <- list(
  FULL       = fit(gt, "FULL"),
  RESTRICTED = fit(gt[, pool_idx], "RESTRICTED"),
  SYNTH      = fit(cbind(gt[, pool_idx], synth), "SYNTH"),
  DROP3      = fit(gt[, setdiff(pool_idx, which(strains %in% TWIN))], "DROP3"))

score <- imap_dfr(fits, function(d, tag) {
  d %>% inner_join(mipf, by = c("strain", "sample")) %>%
    filter(strain != "N2") %>%
    group_by(strain) %>%
    summarise(rmsd = sqrt(mean((frq - mip)^2)), bias = mean(frq - mip),
              rho = suppressWarnings(cor(frq, mip, method = "spearman")),
              .groups = "drop") %>% mutate(fit = tag)
})

wide <- score %>% select(strain, fit, rmsd) %>%
  pivot_wider(names_from = fit, values_from = rmsd) %>%
  left_join(ibs_tab, by = "strain") %>%
  mutate(delta_restricted = RESTRICTED - FULL,
         delta_synth      = SYNTH - FULL,
         is_twin = strain %in% TWIN)
write_tsv(wide, file.path(DIAG, "baugh_reference_ablation.tsv"))

ov <- score %>% group_by(fit) %>%
  summarise(mean_rmsd = mean(rmsd), median_rmsd = median(rmsd), n = n())
cat("\n=== overall per-strain RMSD against MIP (N2 excluded) ===\n"); print(ov)
cat(sprintf("\nRESTRICTED vs FULL: %d of %d strains worse, %d better\n",
    sum(wide$delta_restricted > 0), nrow(wide), sum(wide$delta_restricted < 0)))
cat(sprintf("SYNTH      vs FULL: %d of %d strains worse, %d better\n",
    sum(wide$delta_synth > 0), nrow(wide), sum(wide$delta_synth < 0)))

## --- plots -----------------------------------------------------------------
theme_set(theme_bw(9) + theme(
  plot.title = element_text(face = "bold", size = 9.5),
  plot.subtitle = element_text(size = 7.5, colour = "grey30"),
  panel.grid.minor = element_blank(), legend.position = "top",
  legend.title = element_blank(), legend.margin = margin(0,0,-4,0)))
COL_W <- "#C4302B"; COL_B <- "#2E7BB6"; COL_N <- "#8A8A8A"
lab_n <- 7
sc <- function(d, ycol, title, sub) {
  d2 <- d %>% mutate(dir = ifelse(.data[[ycol]] > 0, "worse without", "better without"))
  hi <- d2 %>% slice_max(abs(.data[[ycol]]), n = lab_n)
  ggplot(d2, aes(FULL, .data[[ycol]] + FULL)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
    geom_point(aes(colour = dir), alpha = 0.8, size = 1.9) +
    geom_point(data = d2 %>% filter(is_twin), shape = 21, size = 3.4,
               stroke = 0.7, colour = "black", fill = NA) +
    geom_text_repel(data = hi, aes(label = strain), size = 2.6,
                    min.segment.length = 0, max.overlaps = Inf, seed = 2) +
    scale_colour_manual(values = c("better without" = COL_B, "worse without" = COL_W)) +
    scale_x_log10() + scale_y_log10() +
    labs(title = title, subtitle = sub,
         x = "per-strain RMSD, FULL reference (103)",
         y = "per-strain RMSD, ablated reference")
}
pA <- sc(wide, "delta_restricted",
  "A  Drop the 3 non-pool strains",
  "Above the line = that strain got WORSE when the 3 were removed. Circled: the 3 near-twins.")
pB <- sc(wide, "delta_synth",
  "B  Negative control: 3 SYNTHETIC sponges instead",
  "IBS-matched fake twins carrying no new information. If sponging explains A, this should match it.")

pC <- ggplot(wide, aes(max_ibs, delta_restricted)) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = 2) +
  geom_point(aes(colour = delta_restricted > 0), alpha = 0.85, size = 2) +
  geom_point(data = wide %>% filter(is_twin), shape = 21, size = 3.6,
             stroke = 0.7, colour = "black", fill = NA) +
  geom_text_repel(data = wide %>% slice_max(abs(delta_restricted), n = lab_n),
                  aes(label = strain), size = 2.6, min.segment.length = 0,
                  max.overlaps = Inf, seed = 3) +
  scale_colour_manual(values = c(`TRUE` = COL_W, `FALSE` = COL_B),
                      labels = c("better without", "worse without")) +
  labs(title = "C  Change against genetic proximity to the removed strains",
       subtitle = sprintf("Spearman rho = %.3f (p = %.3g) between max IBS to an extra and the RMSD penalty",
         cor(wide$max_ibs, wide$delta_restricted, method = "spearman"),
         cor.test(wide$max_ibs, wide$delta_restricted, method = "spearman")$p.value),
       x = "max IBS to any of CX11262, ECA348, NIC260",
       y = "RMSD change when the 3 are removed")

long <- score %>% filter(fit != "DROP3") %>%
  mutate(fit = factor(fit, c("FULL", "RESTRICTED", "SYNTH")))
pD <- ggplot(long, aes(fit, rmsd)) +
  geom_line(aes(group = strain), colour = "grey85", linewidth = 0.25) +
  geom_point(colour = COL_N, size = 1.2, alpha = 0.6) +
  geom_line(data = long %>% filter(strain %in% TWIN),
            aes(group = strain, colour = strain), linewidth = 0.8) +
  geom_point(data = long %>% filter(strain %in% TWIN), aes(colour = strain), size = 2.4) +
  scale_y_log10() +
  labs(title = "D  Every strain across the three references",
       subtitle = "Coloured: the three near-twins of the removed strains",
       x = NULL, y = "per-strain RMSD vs MIP")

fig <- (pA | pB) / (pC | pD) + plot_layout(heights = c(1, 1))
ggsave(file.path(DIAG, "DIAG_baugh_reference_ablation.pdf"), fig, width = 11.5, height = 9.5)
ggsave(file.path(DIAG, "DIAG_baugh_reference_ablation.png"), fig, width = 11.5, height = 9.5, dpi = 200)
cat(sprintf("\nwrote %s/DIAG_baugh_reference_ablation.{pdf,png}\n", DIAG))
