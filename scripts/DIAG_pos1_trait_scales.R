## What the pos-1 pilot phenotype looks like on every scale it could be put on -
##
##   Rscript scripts/DIAG_pos1_trait_scales.R
##     -> plots/diagnostics/DIAG_pos1_trait_scales.{pdf,png}
##
## The 231-strain pilot has 81 strains absent from every pos-1 pool. Every way
## of scoring the response handles them differently, and the choice moves the
## mapping more than any modelling decision downstream of it.
##
## DELTA IS DOMINATED BY A HANDFUL OF STRAINS. Raw delta has skew 4.16 and
## kurtosis 37.3, and the FIVE most extreme strains hold 70.5% of its total sum
## of squares. |delta| correlates with control-pool frequency at rho = +0.875:
## on this scale "responded a lot" and "was abundant to begin with" are nearly
## the same statement. Variance-stabilising helps but does not fix it -- vst is
## skew 1.62, kurtosis 9.79, top five 39.0%, rho = +0.725.
##
## AND LOG RATIOS BLOW UP THE ZEROS. log2(pos/ctrl) is -Inf when a strain is
## absent; the shipped log2fc column lands those 81 at -58 to -69 against an
## observed range of -9.9 to +6.8. Scanning that column directly is scanning 81
## outliers, not a ratio. Any ratio-scale analysis has to either drop the
## censored strains or bound them (rank, or a censored model).
##
## Reads only supplemental_data. Exploratory; nothing in the manuscript reads
## this file.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
PH <- "supplemental_data/phenotypes"; DIAG <- "plots/diagnostics"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)

tr   <- fread(cmd = paste("gzcat", shQuote(file.path(PH,"pos1_2023_association_traits.csv.gz"))))
keep <- tr[is.finite(`vst_ctrl_pos-1_T2`)]$strain
raw  <- fread(cmd = paste("gzcat", shQuote(file.path(PH,"pos1_2023_sample_frequencies.csv.gz"))))
d5   <- raw[depth_cutoff == 5 & strain %in% keep][
           , .(frq = sum(frq)), by = .(strain, sample_info, rnai)]
L <- min(d5[frq > 0]$frq)
Z <- merge(d5[rnai == "pos-1", .(pos = mean(frq)), by = strain],
           d5[rnai != "pos-1", .(ctrl = mean(frq)), by = strain], by = "strain")
Z[, cens := pos == 0]
Z <- merge(Z, tr[, .(strain, delta = `delta_ctrl_pos-1_T2`,
                     vst = `vst_ctrl_pos-1_T2`, log2fc = `log2fc_ctrl_pos-1_T2`)], by = "strain")
Z[, ratio := ifelse(cens, NA_real_, log2(pos / ctrl))]
Z[, thresh := log2(L / ctrl)]

## tobit null fit, then E[y | y <= t] for the censored
nll <- function(p) { m <- p[1]; s <- exp(p[2])
  -(sum(dnorm(Z[cens == FALSE]$ratio, m, s, log = TRUE)) +
    sum(pnorm(Z[cens == TRUE]$thresh, m, s, log.p = TRUE))) }
fit <- optim(c(mean(Z[cens==FALSE]$ratio), log(sd(Z[cens==FALSE]$ratio))), nll, method = "BFGS")
mu <- fit$par[1]; sg <- exp(fit$par[2]); a <- (Z$thresh - mu) / sg
Z[, tobit := ifelse(cens, mu - sg * dnorm(a) / pnorm(a), ratio)]
r <- rep(NA_real_, nrow(Z))
r[Z$cens]  <- mean(seq_len(sum(Z$cens)))
r[!Z$cens] <- sum(Z$cens) + rank(Z[cens == FALSE]$ratio, ties.method = "average")
Z[, rankit := qnorm((r - 0.5) / .N)]

LONG <- melt(Z[, .(strain, cens, ctrl,
                   `delta (absolute)` = delta, `vst (absolute)` = vst,
                   `log2fc as shipped` = log2fc, `log2(pos/ctrl), observed only` = ratio,
                   `tobit imputed` = tobit, `rankit, censored tied` = rankit)],
             id.vars = c("strain","cens","ctrl"), variable.name = "scale", value.name = "y")
LONG[, scale := factor(scale, levels = c("delta (absolute)","vst (absolute)","log2fc as shipped",
                                         "log2(pos/ctrl), observed only","tobit imputed","rankit, censored tied"))]
LONG[, grp := fifelse(cens, "absent from every pos-1 pool (81)", "present (150)")]
COL <- c("absent from every pos-1 pool (81)" = "#C4302B", "present (150)" = "#1B3A6B")

sk <- function(x) mean((x-mean(x))^3)/sd(x)^3
top5 <- function(x) 100*sum(sort((x-mean(x))^2, decreasing=TRUE)[1:5])/sum((x-mean(x))^2)
lab <- LONG[is.finite(y), .(txt = sprintf("skew %+.2f | top 5 hold %.0f%% of SS", sk(y), top5(y))), by = scale]

pA <- ggplot(LONG[is.finite(y)], aes(y, fill = grp)) +
  geom_histogram(bins = 55, colour = NA) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE,
            hjust = -0.04, vjust = 1.5, size = 2.5, colour = "grey25") +
  facet_wrap(~ scale, scales = "free", ncol = 2) +
  scale_fill_manual(values = COL, name = NULL) +
  labs(x = "phenotype value", y = "strains",
       title = "The same experiment, scored six ways",
       subtitle = paste("Red is the 81 strains absent from every pos-1 pool. On the absolute scales they sit inside the",
                        "\ndistribution; log2fc as shipped throws them to -58..-69, far outside the observed range of -9.9..6.8;",
                        "\nthe censored codings put them at the responsive end without inventing an order among them.")) +
  theme_bw(9) + theme(legend.position = "top", panel.grid.minor = element_blank())

pB <- ggplot(Z, aes(ctrl, abs(delta), colour = cens)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "grey45", linetype = "dashed") +
  geom_point(size = 0.9, alpha = 0.7) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(`FALSE` = "#1B3A6B", `TRUE` = "#C4302B"), guide = "none") +
  labs(x = "control-pool frequency (log)", y = "|delta| (log)",
       title = "On the absolute scale, response magnitude is mostly starting abundance",
       subtitle = sprintf(paste("Spearman(|delta|, control frequency) = %+.3f. The dashed line is |delta| = control frequency,",
                                "\nwhich is exactly where the 81 absent strains must fall: their pos-1 frequency is 0, so their",
                                "\ndelta IS minus their starting abundance. They are not censored on this scale -- they are exactly",
                                "\nobserved -- but what is observed about them is largely how abundant they were."),
                          cor(abs(Z$delta), Z$ctrl, method = "spearman"))) +
  theme_bw(9) + theme(panel.grid.minor = element_blank())

fig <- pA / pB + plot_layout(heights = c(2.1, 1)) + plot_annotation(tag_levels = "A")
ggsave(file.path(DIAG,"DIAG_pos1_trait_scales.pdf"), fig, width = 9, height = 12)
ggsave(file.path(DIAG,"DIAG_pos1_trait_scales.png"), fig, width = 9, height = 12, dpi = 180)
cat(sprintf("[%s] wrote DIAG_pos1_trait_scales.{pdf,png}\n", format(Sys.time(), "%H:%M:%S")))
print(lab)
