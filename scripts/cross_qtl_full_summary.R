## Every cross QTL, with effect size and parental frequencies ----------------
##
##   Rscript scripts/make_cross_af_tables.R      # once, stages the AF tables
##   Rscript scripts/cross_qtl_full_summary.R
##     -> plots/TABLE_cross_qtl_full.tsv
##
## WHAT THIS IS NOT
## Figure 2 draws QTL from the same scans but shows a deliberate subset, and
## two constants do the cutting:
##
##   LOD_MIN       <- 100    Figure2.R, drawing only: peaks below it are absent
##   MAX_PEAKS     <- 1      one interval per chromosome per contrast, so a
##                           secondary peak is never called at all
##   PEAK_MIN_FRAC <- 0.30   and a second peak would have needed 30% of the
##                           chromosome maximum to be considered
##
## This table drops the LOD 100 cutoff entirely and calls up to MAX_PEAKS
## peaks per chromosome, keeping only the genome-wide significance threshold
## the scans were designed against: CROSS_THR, which is 3.57, not 100.
##
## The interval rule is unchanged -- the contiguous run of markers within a 5%
## LOD drop of the peak, masked against taller neighbours -- so intervals here
## are comparable with TABLE_cross_qtl_tracks_5pct.tsv from Figure 2.
##
## TWO EFFECT SIZES, BECAUSE THEY ANSWER DIFFERENT QUESTIONS
##   contrast.beta   the model effect at the peak marker, signed by the order
##                   of the contrast label: positive means the FIRST pool named
##                   carries more of parent 1
##   dfreq           the same difference on the allele-frequency scale: parent
##                   1's frequency in pool A minus pool B, each pooled over a
##                   +/- FREQ_WIN_KB window around the peak and weighted by
##                   read depth, which is the estimator load_parent_freq() in
##                   Figure3_common.R uses
##
## f1.a and f1.b are those two window frequencies themselves, so a reader can
## see whether a QTL is a shift between two intermediate frequencies or a pool
## running to fixation.
##
## SHARED BETWEEN THE CROSSES
## Three contrasts were run in both crosses -- HT115 vs mig-6, HT115 vs pos-1,
## and mig-6 vs pos-1 -- and only those can be shared. For each QTL:
##   other.LOD.at.peak   the other cross's LOD for the same contrast at this
##                       exact position, which is the honest comparison: it
##                       does not require the other cross to have called a peak
##   shared              a QTL of the same contrast in the other cross whose
##                       interval overlaps this one
##   shared.suggestive   no overlapping call, but the other cross still clears
##                       CROSS_THR at this position
## The twelve N2 x XZ1516 contrasts with no JU counterpart get NA, not FALSE.
##
## RESOLUTION CAVEAT
## The scans come from the thinned bundle, which keeps the highest-LOD marker
## per 5 kb bin. Peak positions and peak LODs are exact; interval bounds can
## contract by up to one bin per edge, as measured in make_thinned_bundle.R.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

OUT <- "plots"
SD  <- "supplemental_data/mapping"
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

b <- readRDS(file.path(SD, "pooled_cross_bundle_thinned.rds"))

## interval and peak-calling constants: the interval rule of Figure 2, but
## without its display cutoffs
DROP_FRAC     <- 0.05   # interval = contiguous run within 5% of the peak LOD
MIN_GAP       <- 1e6    # buffer masked either side of a called interval
MAX_PEAKS     <- 5      # Figure 2 uses 1
PEAK_MIN_FRAC <- 0.10   # Figure 2 uses 0.30
FREQ_WIN_KB   <- 50     # half-width of the parental-frequency window
THR           <- b$CROSS_THR

CHROMS <- c("I", "II", "III", "IV", "V", "X")

## ===========================================================================
## peaks and intervals
## ===========================================================================
## Faithful to peak_intervals_frac() in Figure2.R, including the masking fix:
## the run is measured only over markers still available, so a secondary
## peak's interval cannot bleed through a taller neighbour.
peak_intervals <- function(d, frac = DROP_FRAC, threshold = THR) {
  d <- d[is.finite(d$LOD), ]
  if (!nrow(d)) return(NULL)
  bind_rows(lapply(split(d, as.character(d$chrom)), function(x) {
    x <- x[order(x$pos), ]
    chrom_max <- max(x$LOD)
    avail <- rep(TRUE, nrow(x))
    out <- list()
    for (k in seq_len(MAX_PEAKS)) {
      if (!any(avail)) break
      i <- which(avail)[which.max(x$LOD[avail])]
      if (x$LOD[i] <= threshold || x$LOD[i] < PEAK_MIN_FRAC * chrom_max) break
      above <- avail & (x$LOD > (x$LOD[i] - x$LOD[i] * frac))
      run <- cumsum(c(TRUE, diff(above) != 0))
      sel <- which(run == run[i])
      lo <- min(sel); hi <- max(sel)
      out[[k]] <- tibble(chrom = as.character(x$chrom[1]),
                         peak.position = x$pos[i], peak.LOD = x$LOD[i],
                         contrast.beta = x$contrast.beta[i], z = x$z[i],
                         lcon = x$pos[lo], rcon = x$pos[hi], peak.rank = k)
      avail[x$pos >= x$pos[lo] - MIN_GAP & x$pos <= x$pos[hi] + MIN_GAP] <- FALSE
    }
    bind_rows(out)
  })) %>% mutate(width.kb = (rcon - lcon) / 1e3)
}

scan_of <- function(key) {
  b$scans[[key]] %>%
    filter(chrom %in% CHROMS) %>%
    transmute(chrom = as.character(chrom), pos = physical.position,
              LOD, contrast.beta, z)
}

## the pools of a contrast, in the order the label names them
pools_of <- function(label) strsplit(label, " vs ", fixed = TRUE)[[1]]
canon_pair <- function(label) paste(sort(pools_of(label)), collapse = "|")

meta <- b$scan_meta %>% as_tibble() %>%
  mutate(pair = map_chr(label, canon_pair),
         pool.a = map_chr(label, ~ pools_of(.x)[1]),
         pool.b = map_chr(label, ~ pools_of(.x)[2]))

qtl <- pmap_dfr(meta %>% select(key, cross, label, pair, pool.a, pool.b),
                function(key, cross, label, pair, pool.a, pool.b) {
  iv <- peak_intervals(scan_of(key))
  if (is.null(iv)) return(NULL)
  iv %>% mutate(cross = cross, contrast = label, pair = pair,
                pool.a = pool.a, pool.b = pool.b, .before = 1)
})

msg(nrow(qtl), " QTL over ", nrow(meta), " contrasts at LOD > ",
    round(THR, 2), " (Figure 2's cutoff of 100 would keep ",
    sum(qtl$peak.LOD > 100), ")")

## ===========================================================================
## parental frequencies in a window around each peak
## ===========================================================================
af <- set_names(map(unique(meta$cross), function(cr) {
  f <- file.path(SD, paste0("cross_af_", cr, ".tsv.gz"))
  if (!file.exists(f))
    stop("missing ", f, "\n  run: Rscript scripts/make_cross_af_tables.R",
         call. = FALSE)
  fread(f, colClasses = list(character = "chrom"))
}), unique(meta$cross))

samp <- read.delim(file.path(SD, "cross_af_samples.tsv"),
                   colClasses = "character") %>% as_tibble() %>%
  ## every contrast is an F2-2 comparison, so the pools are the timepoint-2
  ## samples; the timepoint-1 and OP50 samples are in the AF tables but are
  ## not part of any contrast
  filter(timepoint == "2") %>%
  mutate(pool = sub("g$", "", tolower(condition)))

## depth-weighted parent-1 frequency in the window, the load_parent_freq()
## estimator: sum the counts over markers, then divide
## The arguments are deliberately suffixed: inside a data.table i-expression
## a bare `chrom` or `pos` resolves to the COLUMN, not to the argument, so
## naming them the same silently selects the whole genome as the window.
freq_in_window <- function(cross_, pool_, chrom_, pos_) {
  s <- samp %>% filter(cross == cross_, pool == pool_)
  if (nrow(s) != 1) return(c(f = NA_real_, n = NA_real_))
  d  <- af[[cross_]]
  lo <- pos_ - FREQ_WIN_KB * 1e3
  hi <- pos_ + FREQ_WIN_KB * 1e3
  w  <- d[chrom == chrom_ & pos >= lo & pos <= hi]
  if (!nrow(w)) return(c(f = NA_real_, n = 0, depth = 0))
  p1 <- sum(w[[s$p1_column]]); p2 <- sum(w[[s$p2_column]])
  c(f = if ((p1 + p2) > 0) p1 / (p1 + p2) else NA_real_,
    n = nrow(w), depth = p1 + p2)
}

fr <- pmap_dfr(qtl %>% select(cross, pool.a, pool.b, chrom, peak.position),
               function(cross, pool.a, pool.b, chrom, peak.position) {
  a <- freq_in_window(cross, pool.a, chrom, peak.position)
  bb <- freq_in_window(cross, pool.b, chrom, peak.position)
  tibble(f1.a = a[["f"]], f1.b = bb[["f"]], n.markers.window = a[["n"]],
         depth.a = a[["depth"]], depth.b = bb[["depth"]])
})
qtl <- bind_cols(qtl, fr) %>% mutate(dfreq = f1.a - f1.b)

## ===========================================================================
## is a secondary peak its own QTL, or the shoulder of one sweep?
## ===========================================================================
## Selection in these pools is strong enough that a single sweep can span most
## of a chromosome, and MIN_GAP alone will then call its shoulders as separate
## peaks. The discriminating measurement is the LOD trough between a peak and
## the nearest taller peak on the same chromosome: if the trace never drops
## below the genome-wide threshold between them, they are one signal, not two.
trough_to_taller <- function(key, chrom_, pos_, lod_) {
  d <- scan_of(key)
  d <- d[d$chrom == chrom_, ]
  taller <- qtl %>% filter(contrast == meta$label[meta$key == key],
                           cross == meta$cross[meta$key == key],
                           chrom == chrom_, peak.LOD > lod_)
  if (!nrow(taller)) return(NA_real_)
  near <- taller$peak.position[which.min(abs(taller$peak.position - pos_))]
  seg <- d[d$pos >= min(pos_, near) & d$pos <= max(pos_, near), ]
  if (nrow(seg) < 3) return(NA_real_)
  min(seg$LOD)
}

qtl <- qtl %>%
  mutate(key = meta$key[match(paste(cross, contrast), paste(meta$cross, meta$label))]) %>%
  rowwise() %>%
  mutate(trough.LOD = trough_to_taller(key, chrom, peak.position, peak.LOD)) %>%
  ungroup() %>%
  mutate(separated = is.na(trough.LOD) | trough.LOD < THR) %>%
  select(-key)

## ===========================================================================
## shared between the crosses
## ===========================================================================
lod_at <- function(key, chrom, pos) {
  d <- b$scans[[key]]
  d <- d[as.character(d$chrom) == chrom, ]
  if (!nrow(d)) return(NA_real_)
  d$LOD[which.min(abs(d$physical.position - pos))]
}

shared_cols <- pmap_dfr(qtl %>% select(cross, pair, chrom, peak.position,
                                       lcon, rcon),
                        function(cross, pair, chrom, peak.position, lcon, rcon) {
  om <- meta %>% filter(pair == !!pair, cross != !!cross)
  if (!nrow(om))
    return(tibble(other.cross = NA_character_, other.LOD.at.peak = NA_real_,
                  other.peak.Mb = NA_real_, other.peak.LOD = NA_real_,
                  shared = NA, shared.suggestive = NA))
  ol <- lod_at(om$key[1], chrom, peak.position)
  oq <- qtl %>% filter(cross == om$cross[1], pair == !!pair, chrom == !!chrom,
                       rcon >= !!lcon, lcon <= !!rcon) %>%
    slice_max(peak.LOD, n = 1, with_ties = FALSE)
  tibble(other.cross = om$cross[1], other.LOD.at.peak = ol,
         other.peak.Mb = if (nrow(oq)) oq$peak.position / 1e6 else NA_real_,
         other.peak.LOD = if (nrow(oq)) oq$peak.LOD else NA_real_,
         shared = nrow(oq) > 0,
         shared.suggestive = nrow(oq) == 0 & is.finite(ol) & ol > THR)
})

tab <- bind_cols(qtl, shared_cols) %>%
  transmute(cross, contrast, pair, pool.a, pool.b,
            chrom = factor(chrom, levels = CHROMS),
            peak.Mb = peak.position / 1e6, peak.LOD,
            lcon.Mb = lcon / 1e6, rcon.Mb = rcon / 1e6, width.kb,
            peak.rank, contrast.beta, z,
            f1.a, f1.b, dfreq, n.markers.window, depth.a, depth.b,
            trough.LOD, separated,
            other.cross, other.LOD.at.peak, other.peak.Mb, other.peak.LOD,
            shared, shared.suggestive,
            drawn.in.figure2 = peak.LOD > 100 & peak.rank == 1) %>%
  arrange(cross, contrast, chrom, desc(peak.LOD))

write_tsv(tab, file.path(OUT, "TABLE_cross_qtl_full.tsv"))
msg("wrote ", file.path(OUT, "TABLE_cross_qtl_full.tsv"))

## ---------------------------------------------------------------------------
## what the table says
## ---------------------------------------------------------------------------
cat("\n== QTL counts ==\n")
cat("  genome-wide threshold           LOD ", sprintf("%.2f", THR), "\n", sep = "")
cat("  QTL at that threshold           ", nrow(tab), "\n", sep = "")
cat("  of those, drawn in Figure 2     ", sum(tab$drawn.in.figure2),
    "  (LOD > 100 and the chromosome's top peak)\n", sep = "")
cat("  secondary peaks Figure 2 cannot call (rank > 1)  ",
    sum(tab$peak.rank > 1), "\n", sep = "")
cat("\n  of the ", sum(tab$peak.rank > 1), " secondary peaks, ",
    sum(tab$peak.rank > 1 & tab$separated),
    " are separated from the taller peak by a\n  sub-threshold LOD trough; the other ",
    sum(tab$peak.rank > 1 & !tab$separated),
    " are shoulders of one sweep and should not\n  be counted as independent QTL.\n", sep = "")
cat("  QTL with no covered marker in the +/- ", FREQ_WIN_KB,
    " kb frequency window: ", sum(is.na(tab$f1.a) | is.na(tab$f1.b)),
    "\n", sep = "")

cat("\n== the three contrasts run in both crosses ==\n")
print(as.data.frame(tab %>% filter(!is.na(other.cross)) %>%
  transmute(cross, contrast, chrom, peak.Mb = round(peak.Mb, 2),
            LOD = round(peak.LOD, 1),
            interval = sprintf("%.2f-%.2f", lcon.Mb, rcon.Mb),
            beta = round(contrast.beta, 3),
            f1.a = round(f1.a, 3), f1.b = round(f1.b, 3),
            dfreq = round(dfreq, 3),
            other.LOD = round(other.LOD.at.peak, 1),
            shared, shared.suggestive)), row.names = FALSE)

cat("\n== the twelve N2 x XZ1516 contrasts with no JU counterpart ==\n")
print(as.data.frame(tab %>% filter(is.na(other.cross)) %>%
  transmute(contrast, chrom, peak.Mb = round(peak.Mb, 2),
            LOD = round(peak.LOD, 1), rank = peak.rank,
            beta = round(contrast.beta, 3), dfreq = round(dfreq, 3))),
  row.names = FALSE)
