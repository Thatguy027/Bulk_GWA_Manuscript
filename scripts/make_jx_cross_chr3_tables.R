## Stage the JU2466 x XZ1516 chromosome III right-arm data -------------------
##
##   Rscript scripts/make_jx_cross_chr3_tables.R
##     -> supplemental_data/mapping/jx_cross_chr3_profile.tsv.gz
##        supplemental_data/mapping/jx_cross_chr3_peaks.tsv
##        supplemental_data/mapping/jx_cross_sid2_window.tsv
##
## WHY THIS CROSS MATTERS. JU2466 and XZ1516 BOTH carry sid-2 96K, so T96K
## cannot segregate between them -- yet this cross has a chromosome III
## right-arm QTL. Whatever drives it is therefore not T96K, which makes it the
## one cross that can separate the sid-2 T96K story from anything else on that
## arm.
##
## THE EXPERIMENT IS INCOMPLETE and the tables say so. There is NO HT115
## control pool, so every contrast is one RNAi condition against another rather
## than against an unselected baseline, and a peak can be driven by either
## side. Only three timepoint-2 pools exist (pos-1, mig-6, par-1); the
## timepoint-1 samples have no contrast files.
##
## SOURCE IS OUTSIDE THIS REPOSITORY: ../bulkGWAS/xqtl_analysis/NJX_rnai.
## Set JX_PLOTS to override.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse)})

JX  <- Sys.getenv("JX_PLOTS",
       "/Users/Stefan/UCLA/Projects/bulkGWAS/xqtl_analysis/NJX_rnai/plots")
OUT <- "supplemental_data/mapping"
CONTRASTS <- c("pos1-par1","pos1-mig6","mig6-par1")
SID2 <- 13680248L; GEN <- 13305359L      # sid-2 T96K, and the general:4 peak
WIN  <- 25000L
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
if (!dir.exists(JX)) stop("JU2466 x XZ1516 plots not found at ", JX, call. = FALSE)

prof <- list(); peaks <- list(); win <- list()
for (c in CONTRASTS) {
  f <- file.path(JX, sprintf("JU2466_XZ1516_F2-2_contrast_%s_10000_plot_DF.tsv", c))
  if (!file.exists(f)) stop("missing ", f, call. = FALSE)
  d <- read_tsv(f, show_col_types = FALSE) %>% filter(chrom == "III")
  at <- function(p) d$LOD[which.min(abs(d$physical.position - p))]

  ## a thinned profile over the right arm, enough to draw and to re-check
  prof[[c]] <- d %>% filter(physical.position >= 12.4e6, physical.position <= 13.9e6) %>%
    mutate(bin = floor(physical.position / 2000)) %>%
    group_by(bin) %>%
    slice_max(LOD, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(contrast = c, position = physical.position, LOD)

  r <- d %>% filter(physical.position > 12e6)
  i <- which.max(r$LOD)
  peaks[[c]] <- tibble(contrast = c,
                       right.peak = r$physical.position[i], right.LOD = r$LOD[i],
                       LOD.at.general = at(GEN), LOD.at.sid2 = at(SID2),
                       left.peak = { l <- d %>% filter(physical.position < 10e6)
                                     l$physical.position[which.max(l$LOD)] },
                       left.LOD = { l <- d %>% filter(physical.position < 10e6)
                                    max(l$LOD) })

  ## parental allele fractions in a window on sid-2, per pool, from raw counts
  w <- d %>% filter(abs(physical.position - SID2) <= WIN)
  for (sm in unique(sub("^p1_", "", grep("^p1_", names(d), value = TRUE)))) {
    p1 <- sum(w[[paste0("p1_", sm)]], na.rm = TRUE)
    p2 <- sum(w[[paste0("p2_", sm)]], na.rm = TRUE)
    win[[sm]] <- tibble(sample = sm,
                        condition = str_extract(sm, "(pos1|par1|mig6|rpn12|vha5)$"),
                        p1.JU2466 = p1, p2.XZ1516 = p2,
                        frac.JU2466 = p1 / (p1 + p2), n = p1 + p2)
  }
}
write_tsv(bind_rows(prof), file.path(OUT, "jx_cross_chr3_profile.tsv.gz"))
write_tsv(bind_rows(peaks), file.path(OUT, "jx_cross_chr3_peaks.tsv"))
write_tsv(bind_rows(win) %>% distinct(sample, .keep_all = TRUE) %>% arrange(desc(frac.JU2466)),
          file.path(OUT, "jx_cross_sid2_window.tsv"))

cat("\n== chromosome III peaks, JU2466 x XZ1516 ==\n")
print(as.data.frame(bind_rows(peaks) %>%
  transmute(contrast, right.peak.Mb = round(right.peak/1e6, 3),
            right.LOD = round(right.LOD, 1),
            LOD.at.general = round(LOD.at.general, 1),
            LOD.at.sid2 = round(LOD.at.sid2, 1))), row.names = FALSE)
cat("\n== JU2466 allele fraction at sid-2 +/- 25 kb ==\n")
print(as.data.frame(bind_rows(win) %>% distinct(sample, .keep_all = TRUE) %>%
  transmute(condition, frac.JU2466 = round(frac.JU2466, 3), n = round(n)) %>%
  arrange(desc(frac.JU2466))), row.names = FALSE)
msg("wrote three tables to ", OUT)
