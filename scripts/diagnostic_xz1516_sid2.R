## Does XZ1516's extra sid-2 variation explain the JU2466 x XZ1516 chrIII QTL?
##
##   Rscript scripts/diagnostic_xz1516_sid2.R
##     -> plots/diagnostics/TABLE_xz1516_sid2_variants.tsv
##        plots/diagnostics/TABLE_jx_chr3_coincidence.tsv
##
## JU2466 and XZ1516 both carry sid-2 96K, so T96K cannot segregate between
## them, yet their cross has a chromosome III right-arm QTL. Two hypotheses,
## tested here against data already in the repository.
##
##   1  AN ALLELIC SERIES AT sid-2. XZ1516 carries five missense variants
##      JU2466 does not, plus a splice-region variant. If those make its SID-2
##      handle dsRNA differently, sid-2 itself is still the gene. Figure 4's
##      mechanism for T96K is charge -- 96K adds a positive charge in the most
##      positive solvent-exposed pocket of an acidic domain -- so the test is
##      whether XZ1516's extras can act on charge at all, and whether they are
##      anywhere near that pocket.
##
##   2  IT IS THE OTHER CHROMOSOME III LOCUS. The mig-6 census found an
##      independent general-response locus at III:13.31 Mb, ~370 kb proximal to
##      sid-2 and classified general:4 in N2 x XZ1516. If the JU2466 x XZ1516
##      peak sits there rather than on sid-2, the two crosses converge on a
##      locus that is not sid-2 at all.
##
## Reads only staged tables, so it runs from a clone.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse)})
ST <- "supplemental_data/structure"; MAP <- "supplemental_data/mapping"
OUT <- "plots/diagnostics"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SID2 <- 13680248; GEN <- 13305359

## ---- 1. can XZ1516's extra variants act through charge? -------------------
## Formal side-chain charge at neutral pH. Only a substitution that changes
## this can act on the electrostatic mechanism at all; a neutral-to-neutral
## change cannot, wherever it sits.
FORMAL <- c(D = -1, E = -1, K = 1, R = 1, H = 0.1,
            A = 0, V = 0, M = 0, Q = 0, P = 0, T = 0, I = 0, L = 0,
            S = 0, N = 0, G = 0, F = 0, W = 0, Y = 0, C = 0)
q <- read_tsv(file.path(ST, "sid2_local_charge.tsv"), show_col_types = FALSE)
t96 <- q %>% filter(resid == 96)

V <- tribble(
  ~resid, ~from, ~to, ~carrier,
  78,  "D", "A", "XZ1516 only",   141, "M", "V", "XZ1516 only",
  144, "Q", "P", "XZ1516 only",   151, "A", "I", "XZ1516 only",
  209, "L", "M", "XZ1516 only",
  96,  "T", "K", "shared 96K",    153, "P", "T", "shared",
  5,   "V", "L", "JU1793 only") %>%
  left_join(q %>% select(resid, q_local_pH44, x, y, z), by = "resid") %>%
  mutate(substitution = paste0(from, resid, to),
         formal.charge.change = FORMAL[to] - FORMAL[from],
         local.charge = q_local_pH44,
         dist.to.T96 = sqrt((x - t96$x)^2 + (y - t96$y)^2 + (z - t96$z)^2),
         modelled = !is.na(q_local_pH44),
         charge.capable = formal.charge.change != 0) %>%
  select(resid, substitution, carrier, formal.charge.change, local.charge,
         dist.to.T96, modelled, charge.capable) %>%
  arrange(carrier, resid)
write_tsv(V, file.path(OUT, "TABLE_xz1516_sid2_variants.tsv"))

cat("== 1. XZ1516's extra sid-2 variants, and whether charge can be their route ==\n")
print(as.data.frame(V %>% transmute(substitution, carrier,
  `formal Dq` = formal.charge.change,
  `local charge` = ifelse(is.na(local.charge), NA, round(local.charge, 2)),
  `A to T96` = ifelse(is.na(dist.to.T96), NA, round(dist.to.T96, 1)),
  modelled)), row.names = FALSE)
cat(sprintf("\n  ectodomain local charge: median %.2f, range %.2f to %.2f\n",
            median(q$q_local_pH44), min(q$q_local_pH44), max(q$q_local_pH44)))
cat(sprintf("  T96 local charge %.2f, %.0f%% percentile of the ectodomain\n",
            t96$q_local_pH44, 100 * mean(q$q_local_pH44 <= t96$q_local_pH44)))
xz <- V %>% filter(carrier == "XZ1516 only")
cat(sprintf("  of XZ1516's %d extra variants, %d can change charge at all: %s\n",
            nrow(xz), sum(xz$charge.capable),
            paste(xz$substitution[xz$charge.capable], collapse = ", ")))

## ---- 2. where does the JU2466 x XZ1516 peak actually sit? -----------------
pk <- read_tsv(file.path(MAP, "jx_cross_chr3_peaks.tsv"), show_col_types = FALSE)
cqf <- read_tsv("plots/TABLE_cross_qtl_full.tsv", show_col_types = FALSE) %>%
  filter(cross == "N2xXZ1516", contrast == "ht115 vs mig6", chrom == "III",
         peak.rank == 1)
lo <- cqf$lcon.Mb * 1e6; hi <- cqf$rcon.Mb * 1e6
co <- pk %>% mutate(
  peak.in.general.interval = right.peak >= lo & right.peak <= hi,
  closer.to = ifelse(abs(right.peak - GEN) < abs(right.peak - SID2),
                     "general locus 13.31", "sid-2 13.68"),
  kb.from.general = (right.peak - GEN)/1000,
  kb.from.sid2 = (right.peak - SID2)/1000)
write_tsv(co, file.path(OUT, "TABLE_jx_chr3_coincidence.tsv"))

cat("\n== 2. the JU2466 x XZ1516 chrIII right-arm peak against the two candidates ==\n")
cat(sprintf("  N2 x XZ1516 general:4 locus: peak %.3f Mb, interval %.3f-%.3f Mb, LOD %.0f\n",
            cqf$peak.Mb, cqf$lcon.Mb, cqf$rcon.Mb, cqf$peak.LOD))
cat(sprintf("  sid-2 T96K sits at %.3f Mb, %.0f kb distal of that peak\n\n",
            SID2/1e6, (SID2 - GEN)/1000))
print(as.data.frame(co %>% transmute(contrast,
  `peak (Mb)` = round(right.peak/1e6, 3), `peak LOD` = round(right.LOD, 1),
  `LOD at 13.31` = round(LOD.at.general, 1),
  `LOD at sid-2` = round(LOD.at.sid2, 1),
  `in general interval` = peak.in.general.interval, `closer to` = closer.to)),
  row.names = FALSE)

w <- read_tsv(file.path(MAP, "jx_cross_sid2_window.tsv"), show_col_types = FALSE)
cat("\n== direction, from raw counts at sid-2 +/- 25 kb ==\n")
print(as.data.frame(w %>% transmute(condition, frac.JU2466 = round(frac.JU2466, 3),
                                    n = round(n))), row.names = FALSE)
cat("  no HT115 pool exists for this cross, so par-1 is the least-selected\n")
cat("  reference rather than an unselected baseline.\n")
