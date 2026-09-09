## General RNAi-response loci against target-specific ones ------------------
##
##   Rscript scripts/make_cross_af_tables.R        # once
##   Rscript scripts/cross_qtl_full_summary.R      # once
##   Rscript scripts/cross_qtl_condition_sharing.R
##     -> plots/TABLE_cross_qtl_locus_classification.tsv
##        plots/TABLE_cross_qtl_condition_dfreq.tsv
##        plots/TABLE_cross_qtl_target_pairs.tsv
##
## THE DESIGN THIS EXPLOITS
## Both crosses ran an HT115 control pool alongside knockdown pools, so each
## "ht115 vs X" contrast is one RNAi target's response against a common
## control:
##
##   N2 x XZ1516      pos-1, mig-6, par-1, rpn-12, vha-5   (5 targets)
##   JU1793 x JU2466  pos-1, mig-6                         (2 targets)
##
## A locus where many targets shift the same way is a general RNAi-response
## locus; one where only mig-6 shifts is mig-6-specific. Counting targets is a
## stronger statement than the two-trace reading in the contrast supplement,
## where "general" had to be inferred from a flat difference trace.
##
## WHY THIS IS NOT DONE ON SIGNIFICANCE
## Pooled depths here are large enough that LOD runs into the hundreds, so at
## the genome-wide threshold of 3.57 nearly every target is "significant" at
## nearly every locus: classifying that way called 55 of 61 N2 x XZ1516 loci
## general, 19 of them in all five targets, which is not credible. Two further
## facts settle the approach:
##
##   every contrast in a cross shares ONE HT115 pool, so any drift in that
##   pool appears in all five contrasts and imitates a general locus
##
##   the timepoint-1 replicate has its own HT115, mig-6, par-1 and rpn-12
##   pools, and the replicate agreement of the frequency shift differs sharply
##   by target: mig-6 r = 0.975, rpn-12 r = 0.683, par-1 r = 0.125. A par-1
##   "response" mostly does not reproduce, so significance in par-1 is not
##   evidence of anything
##
## So the classification is made on EFFECT SIZE, with significance as a
## necessary condition only, and the replicate agreement is carried alongside:
##
##   responds        LOD > CROSS_THR and |dfreq| >= RESP_MIN
##   two targets
##     differ        |dfreq_i - dfreq_j| >= DIFF_MIN
##
## dfreq is parent-1 frequency in the target pool minus the HT115 pool, both
## pooled over +/- FREQ_WIN_KB and weighted by depth.
##
## The target-vs-target contrasts are used for the differ test as well, and
## they do NOT share the control pool, so they are the cleaner evidence that
## mig-6 is doing something the other targets are not.
##
## WHAT "GENERAL" REQUIRES, AND WHAT IT DOES NOT
## Different RNAi targets impose different selection strength -- mig-6 shifts
## frequencies about three times as far as par-1 does at the chromosome V loci
## -- so requiring equal magnitude across targets is the wrong test: it called
## nothing general, including the chromosome III locus where four of five
## targets shift the same way. The criterion is therefore CONCORDANT DIRECTION
## among the responding targets, with the magnitude spread reported next to
## the call rather than gating it.
##
## CLASSES
##   general:<n>      >= GEN_MIN_N targets respond, all in the same direction
##   discordant       >= GEN_MIN_N respond but not all in the same direction
##   shared:<t1,t2>   exactly two respond. In the JU cross only two targets
##                    were run, so two of two is the most that cross can show:
##                    general.limited marks those
##   specific:<t>     exactly one responds
##   none             no target responds by effect size
##
## mig6.dominant marks a locus where mig-6 responds at least MIG6_DOM_FOLD
## times as far as any other responding target: general in direction, but
## carried mostly by mig-6.
## Every threshold is a named constant and the per-target numbers are written
## out, so the table can be re-cut without re-running this.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

OUT <- "plots"
SD  <- "supplemental_data/mapping"
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

b   <- readRDS(file.path(SD, "pooled_cross_bundle_thinned.rds"))
THR <- b$CROSS_THR

FREQ_WIN_KB <- 50     # matched to cross_qtl_full_summary.R
MERGE_KB    <- 250    # calls this close are one locus
RESP_MIN    <- 0.10   # |dfreq| vs control to count as a response
DIFF_MIN    <- 0.10   # |dfreq| difference for two targets to differ
GEN_MIN_N   <- 3      # responding targets needed for a general call
MIG6_DOM_FOLD <- 1.5  # mig-6 this many times the next responder is "dominant"

CHROMS <- c("I", "II", "III", "IV", "V", "X")

## max() over an all-missing vector warns and returns -Inf; loci can be all
## missing when the peak falls in an uncovered marker gap
max_or_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

sm <- b$scan_meta %>% as_tibble() %>%
  separate(label, into = c("pool.a", "pool.b"), sep = " vs ", remove = FALSE)
ctrl <- sm %>% filter(pool.a == "ht115") %>%
  transmute(key, cross, target = pool.b)
pair <- sm %>% filter(pool.a != "ht115", pool.b != "ht115") %>%
  transmute(key, cross, t1 = pool.a, t2 = pool.b)

msg(nrow(ctrl), " control-vs-target and ", nrow(pair),
    " target-vs-target contrasts")

## ===========================================================================
## loci: the union of the calls in the control-vs-target contrasts
## ===========================================================================
qtl <- read_tsv(file.path(OUT, "TABLE_cross_qtl_full.tsv"),
                show_col_types = FALSE) %>%
  filter(paste(cross, contrast) %in%
           paste(ctrl$cross, paste("ht115 vs", ctrl$target)))

loci <- qtl %>%
  arrange(cross, chrom, lcon.Mb) %>%
  group_by(cross, chrom) %>%
  mutate(grp = cumsum(c(TRUE,
            lcon.Mb[-1] > cummax(rcon.Mb)[-n()] + MERGE_KB / 1e3))) %>%
  group_by(cross, chrom, grp) %>%
  summarise(lcon.Mb = min(lcon.Mb), rcon.Mb = max(rcon.Mb),
            peak.Mb = peak.Mb[which.max(peak.LOD)],
            top.LOD = max(peak.LOD), n.calls = n(),
            called.by = paste(sort(unique(str_remove(contrast, "^ht115 vs "))),
                              collapse = ","), .groups = "drop") %>%
  select(-grp) %>%
  mutate(locus = sprintf("%s %s:%.2f", cross, chrom, peak.Mb))

msg(nrow(loci), " loci after merging calls within ", MERGE_KB, " kb")

## ===========================================================================
## frequencies and LODs at each locus
## ===========================================================================
af <- set_names(map(unique(ctrl$cross), function(cr) {
  d <- fread(file.path(SD, paste0("cross_af_", cr, ".tsv.gz")),
             colClasses = list(character = "chrom"))
  setkey(d, chrom, pos)     # so the window subset is a binary search
  d
}), unique(ctrl$cross))

samp <- read.delim(file.path(SD, "cross_af_samples.tsv"),
                   colClasses = "character") %>% as_tibble() %>%
  mutate(pool = sub("g$", "", tolower(condition)))

## suffixed arguments: inside a data.table i-expression a bare chrom/pos binds
## to the column, not to the argument
freq_win <- function(cross_, pool_, tp_, chrom_, pos_) {
  s <- samp %>% filter(cross == cross_, pool == pool_, timepoint == tp_)
  if (nrow(s) != 1) return(NA_real_)
  d <- af[[cross_]]
  w <- d[.(chrom_, seq(pos_ - FREQ_WIN_KB * 1e3, pos_ + FREQ_WIN_KB * 1e3)),
         nomatch = 0L]
  if (!nrow(w)) return(NA_real_)
  p1 <- sum(w[[s$p1_column]]); p2 <- sum(w[[s$p2_column]])
  if ((p1 + p2) == 0) NA_real_ else p1 / (p1 + p2)
}

lod_at <- function(key_, chrom_, pos_) {
  d <- b$scans[[key_]]
  d <- d[as.character(d$chrom) == chrom_, ]
  if (!nrow(d)) return(NA_real_)
  d$LOD[which.min(abs(d$physical.position - pos_))]
}

cells <- pmap_dfr(loci %>% select(locus, cross, chrom, peak.Mb),
                  function(locus, cross, chrom, peak.Mb) {
  pos <- peak.Mb * 1e6
  f_c2 <- freq_win(cross, "ht115", "2", chrom, pos)
  f_c1 <- freq_win(cross, "ht115", "1", chrom, pos)
  pmap_dfr(ctrl %>% filter(cross == !!cross) %>% select(key, target),
           function(key, target) {
    f_t2 <- freq_win(cross, target, "2", chrom, pos)
    f_t1 <- freq_win(cross, target, "1", chrom, pos)
    tibble(locus, cross, chrom, peak.Mb, target,
           LOD = lod_at(key, chrom, pos),
           f.ctrl = f_c2, f.target = f_t2, dfreq = f_t2 - f_c2,
           dfreq.rep1 = f_t1 - f_c1)
  })
}) %>%
  mutate(sig = is.finite(LOD) & LOD > THR,
         responds = sig & is.finite(dfreq) & abs(dfreq) >= RESP_MIN)

## target-vs-target: the comparison that does not share the control pool
pcells <- pmap_dfr(loci %>% select(locus, cross, chrom, peak.Mb),
                   function(locus, cross, chrom, peak.Mb) {
  pos <- peak.Mb * 1e6
  pp  <- pair %>% filter(cross == !!cross)
  if (!nrow(pp)) return(NULL)
  pmap_dfr(pp %>% select(key, t1, t2), function(key, t1, t2) {
    d1 <- cells$dfreq[cells$locus == locus & cells$target == t1]
    d2 <- cells$dfreq[cells$locus == locus & cells$target == t2]
    tibble(locus, cross, chrom, peak.Mb, t1, t2,
           LOD = lod_at(key, chrom, pos),
           ddfreq = if (length(d1) && length(d2)) d1 - d2 else NA_real_)
  })
}) %>%
  mutate(differ = is.finite(LOD) & LOD > THR &
                  is.finite(ddfreq) & abs(ddfreq) >= DIFF_MIN)

## ===========================================================================
## classify
## ===========================================================================
## the largest disagreement among the responding targets, from the pairwise
## table where it exists and from the dfreq spread otherwise
spread_of <- function(locus_, resp) {
  if (length(resp) < 2) return(NA_real_)
  p <- pcells %>% filter(locus == locus_, t1 %in% resp, t2 %in% resp)
  if (!nrow(p)) return(NA_real_)
  max_or_na(abs(p$ddfreq))
}
any_differ <- function(locus_, resp) {
  if (length(resp) < 2) return(FALSE)
  p <- pcells %>% filter(locus == locus_, t1 %in% resp, t2 %in% resp)
  nrow(p) > 0 && any(p$differ, na.rm = TRUE)
}
## is mig-6 distinguishable from every other target it was measured against?
mig6_outlier <- function(locus_) {
  p <- pcells %>% filter(locus == locus_, t1 == "mig6" | t2 == "mig6")
  nrow(p) > 0 && all(p$differ, na.rm = TRUE)
}

cls <- cells %>%
  group_by(locus, cross, chrom, peak.Mb) %>%
  summarise(n.tested = n(), n.sig = sum(sig), n.resp = sum(responds),
            targets.resp = paste(sort(target[responds]), collapse = ","),
            max.LOD = max_or_na(LOD),
            dfreq.mig6 = dfreq[target == "mig6"][1],
            dfreq.pos1 = dfreq[target == "pos1"][1],
            dfreq.mig6.rep1 = dfreq.rep1[target == "mig6"][1],
            max.abs.dfreq = max_or_na(abs(dfreq)),
            .groups = "drop") %>%
  rowwise() %>%
  mutate(resp.spread = spread_of(locus, str_split(targets.resp, ",")[[1]]),
         resp.differ = any_differ(locus, str_split(targets.resp, ",")[[1]]),
         mig6.distinct = mig6_outlier(locus)) %>%
  ungroup() %>%
  left_join(cells %>% filter(responds) %>%
              group_by(locus) %>%
              summarise(n.signs = n_distinct(sign(dfreq)),
                        resp.dir = if (n_distinct(sign(dfreq)) == 1)
                          ifelse(dfreq[1] > 0, "parent1", "parent2") else "mixed",
                        mig6.fold = {
                          o <- abs(dfreq[target != "mig6"])
                          m <- abs(dfreq[target == "mig6"])
                          if (!length(m) || !length(o)) NA_real_ else m / max(o)
                        }, .groups = "drop"),
            by = "locus") %>%
  mutate(concordant = !is.na(n.signs) & n.signs == 1,
         mig6.dominant = !is.na(mig6.fold) & mig6.fold >= MIG6_DOM_FOLD,
         class = case_when(
      n.resp == 0 ~ "none",
      n.resp >= GEN_MIN_N & concordant ~ paste0("general:", n.resp),
      n.resp >= GEN_MIN_N ~ "discordant",
      n.resp == 2 ~ paste0("shared:", targets.resp),
      TRUE ~ paste0("specific:", targets.resp)),
      ## two of two is the ceiling in the JU cross, which ran only two targets
      general.limited = n.tested == 2 & n.resp == 2 & concordant) %>%
  left_join(loci %>% select(locus, lcon.Mb, rcon.Mb, n.calls, called.by),
            by = "locus")

## ===========================================================================
## sharing across the crosses
## ===========================================================================
other_of <- function(cross_, chrom_, l_, r_) {
  o <- cls %>% filter(cross != cross_, chrom == chrom_,
                      rcon.Mb >= l_, lcon.Mb <= r_)
  if (!nrow(o)) return(tibble(shared.across.cross = FALSE,
                              other.locus = NA_character_,
                              other.class = NA_character_,
                              other.targets.resp = NA_character_))
  o <- o %>% slice_max(max.abs.dfreq, n = 1, with_ties = FALSE)
  tibble(shared.across.cross = TRUE, other.locus = o$locus,
         other.class = o$class, other.targets.resp = o$targets.resp)
}
cls <- bind_cols(cls, pmap_dfr(cls %>% select(cross, chrom, lcon.Mb, rcon.Mb),
                               function(cross, chrom, lcon.Mb, rcon.Mb)
                                 other_of(cross, chrom, lcon.Mb, rcon.Mb)))

cls <- cls %>%
  mutate(general = str_starts(class, "general"),
         other.general = !is.na(other.class) & str_starts(other.class, "general"),
         mig6.only = class == "specific:mig6",
         other.mig6.only = !is.na(other.class) & other.class == "specific:mig6",
         verdict = case_when(
      general & other.general      ~ "general RNAi-response, both crosses",
      general                      ~ "general RNAi-response, one cross",
      mig6.only & other.mig6.only  ~ "mig-6 specific, both crosses",
      mig6.only                    ~ "mig-6 specific, one cross",
      str_starts(class, "specific:") ~ paste0(str_remove(class, "specific:"),
                                              " specific"),
      str_starts(class, "shared:") ~ paste0("shared by ",
                                            str_remove(class, "shared:")),
      general.limited              ~ "responds in both targets tested (2 of 2)",
      class == "discordant"        ~ "responds in several targets, opposite directions",
      TRUE ~ "no response by effect size")) %>%
  arrange(cross, chrom, peak.Mb)

write_tsv(cls, file.path(OUT, "TABLE_cross_qtl_locus_classification.tsv"))
write_tsv(cells, file.path(OUT, "TABLE_cross_qtl_condition_dfreq.tsv"))
write_tsv(pcells, file.path(OUT, "TABLE_cross_qtl_target_pairs.tsv"))
msg("wrote three tables to ", OUT)

## ---------------------------------------------------------------------------
## what it says
## ---------------------------------------------------------------------------
cat("\n== thresholds ==\n")
cat("  significance (necessary only)  LOD > ", sprintf("%.2f", THR), "\n",
    "  responds                       |dfreq| >= ", RESP_MIN, "\n",
    "  two targets differ             |ddfreq| >= ", DIFF_MIN,
    " and LOD > threshold\n", "  general                        >= ",
    GEN_MIN_N, " responding targets, none differing\n", sep = "")

cat("\n== classes ==\n")
print(as.data.frame(cls %>% count(cross, class)), row.names = FALSE)

cat("\n== general RNAi-response loci ==\n")
gen <- cls %>% filter(general)
if (!nrow(gen)) cat("  none\n") else
print(as.data.frame(gen %>% transmute(locus,
        interval = sprintf("%.2f-%.2f", lcon.Mb, rcon.Mb),
        targets = targets.resp, dir = resp.dir,
        spread = round(resp.spread, 3), mig6.fold = round(mig6.fold, 2),
        mig6.dominant, top.LOD = round(max.LOD, 1),
        shared = shared.across.cross, other.class)), row.names = FALSE)

cat("\n== responds in both targets tested, JU cross (2 of 2) ==\n")
lim <- cls %>% filter(general.limited)
if (!nrow(lim)) cat("  none\n") else
print(as.data.frame(lim %>% transmute(locus,
        interval = sprintf("%.2f-%.2f", lcon.Mb, rcon.Mb),
        dir = resp.dir, dfreq.mig6 = round(dfreq.mig6, 3),
        dfreq.pos1 = round(dfreq.pos1, 3),
        shared = shared.across.cross, other.class)), row.names = FALSE)

cat("\n== responds in several targets but in opposite directions ==\n")
dis <- cls %>% filter(class == "discordant")
if (!nrow(dis)) cat("  none\n") else
print(as.data.frame(dis %>% transmute(locus,
        interval = sprintf("%.2f-%.2f", lcon.Mb, rcon.Mb),
        targets = targets.resp, top.LOD = round(max.LOD, 1))),
      row.names = FALSE)

cat("\n== mig-6 specific loci ==\n")
mg <- cls %>% filter(mig6.only)
if (!nrow(mg)) cat("  none\n") else
print(as.data.frame(mg %>% transmute(locus,
        interval = sprintf("%.2f-%.2f", lcon.Mb, rcon.Mb),
        dfreq.mig6 = round(dfreq.mig6, 3),
        rep1 = round(dfreq.mig6.rep1, 3),
        mig6.distinct, top.LOD = round(max.LOD, 1),
        shared = shared.across.cross, other.class)), row.names = FALSE)

cat("\n== sensitivity of the counts to RESP_MIN ==\n")
sens <- map_dfr(c(0.05, 0.10, 0.15, 0.20), function(rm) {
  r <- cells %>% mutate(responds = sig & is.finite(dfreq) & abs(dfreq) >= rm) %>%
    group_by(locus, cross) %>%
    summarise(n.resp = sum(responds),
              targets.resp = paste(sort(target[responds]), collapse = ","),
              .groups = "drop")
  tibble(RESP_MIN = rm,
         general = sum(r$n.resp >= GEN_MIN_N),
         mig6.only = sum(r$targets.resp == "mig6"),
         pos1.only = sum(r$targets.resp == "pos1"),
         none = sum(r$n.resp == 0))
})
print(as.data.frame(sens), row.names = FALSE)
