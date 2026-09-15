## Precompute for the pooled-panel optimiser ---------------------------------
##
##   Rscript scripts/DIAG_pool_optimizer_precompute.R
##     -> .pool_opt_cache/precomputed.rds   (gitignored, ~0.3 GB)
##
## Builds everything the optimiser needs to score a candidate panel in
## milliseconds, so that a swap-based search can afford tens of thousands of
## evaluations.
##
## CANDIDATE UNIVERSE. The CaeNDR 20231213 hard-filter panel, 611 isotypes,
## restricted to the 571 that have per-strain divergent-region calls. The 40
## without them are dropped rather than treated as having none: a strain with no
## divergent calls would have all of its privateness scored as the untaxed kind,
## which is exactly the quantity being optimised, so keeping them would hand
## them an unearned advantage.
##
## TWO MARKER SETS, because the two halves of the problem want opposite things.
## The plink2 extraction steps are run by this script; plink2 must be on PATH
## (or named by PLINK2). Intermediates land in .pool_opt_cache, which is
## gitignored -- nothing here is small enough or stable enough to commit.
##
##   rare    carriers 1-25 in the universe, MtDNA dropped. These are what make
##           a strain identifiable to the deconvolution. MtDNA is dropped
##           because it does not recombine, so its markers are one linked block
##           masquerading as many independent ones.
##   pruned  carriers 30-581, LD-pruned at r2 0.5 in 50 kb windows. These are
##           what association mapping actually tests.
##
## WHAT MAKES THE SEARCH FAST
##   - privateness is evaluated from a (marker, strain) carrier list rather than
##     a genotype matrix, so a panel's private-marker counts cost one pass over
##     ~11 M integer pairs instead of a 2.4 M x n matrix operation
##   - the panel kinship for ANY subset is a double-centred submatrix of one
##     precomputed 571 x 571 cross-product: H G_S H, with H = I - J/n. Centring
##     each marker within the panel is exactly post-multiplication by H, so no
##     marker ever has to be re-centred
##
## SOURCES OUTSIDE THE REPOSITORY -- CENDR_PLINK, CENDR_DIVERGENT. Exploratory,
## on the pool_optimization branch.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(BEDMatrix)
})

PLINK <- Sys.getenv("CENDR_PLINK",
                    "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213/filtered_geno10")
DIV   <- Sys.getenv("CENDR_DIVERGENT",
                    "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213/20231213_c_elegans_divergent_regions_strain.bed")
CACHE   <- ".pool_opt_cache"
CHUNK   <- 100000L
CAR_MAX <- 25L        # carriers in the universe for a marker to count as rare
CAR_MIN <- 30L        # carriers for a marker to be worth LD-pruning as common
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
for (f in c(paste0(PLINK, ".bed"), DIV)) if (!file.exists(f)) stop("missing source: ", f)
dir.create(CACHE, showWarnings = FALSE)
plink2 <- Sys.getenv("PLINK2", "plink2")
if (nzchar(Sys.which(plink2)) == FALSE) stop("plink2 not on PATH", call. = FALSE)
run <- function(...) {
  st <- system2(plink2, c("--allow-extra-chr", "--silent", ...))
  if (st != 0) stop("plink2 failed with status ", st, call. = FALSE)
}

## --- marker and strain selection, done with plink2 --------------------------
## Rerunning is cheap relative to the R work below, so the selection is redone
## every time rather than cached; that keeps the cache from going stale against
## an edited threshold.
if (!file.exists(file.path(CACHE, "freq.acount")))
  run("--bfile", PLINK, "--freq", "counts", "--out", file.path(CACHE, "freq"))
f <- fread(file.path(CACHE, "freq.acount"))
setnames(f, c("CHROM", "ID", "REF", "ALT", "ALT_CTS", "OBS_CT"))
f[, car := ALT_CTS / 2]     # isotypes are homozygous, so plink's 0/2 is 0/1
## MtDNA is dropped from BOTH sets: it does not recombine, so its markers are
## one linked block that would otherwise count as many independent ones.
fwrite(f[CHROM != "MtDNA" & car >= 1 & car <= CAR_MAX, .(ID)],
       file.path(CACHE, "rare.ids"), col.names = FALSE)
fwrite(f[CHROM != "MtDNA" & car >= CAR_MIN & car <= max(f$OBS_CT) / 2 - CAR_MIN, .(ID)],
       file.path(CACHE, "common.ids"), col.names = FALSE)

fam <- fread(paste0(PLINK, ".fam"), header = FALSE)
bed0 <- fread(DIV, col.names = c("chrom", "start", "end", "strain"))
uni <- fam[V1 %in% bed0$strain, .(V1, V2)]
msg(sprintf("candidate universe: %d of %d isotypes have divergent-region calls",
            nrow(uni), nrow(fam)))
fwrite(uni, file.path(CACHE, "universe.fam"), sep = "\t", col.names = FALSE)

run("--bfile", PLINK, "--keep", file.path(CACHE, "universe.fam"),
    "--extract", file.path(CACHE, "rare.ids"), "--make-bed",
    "--out", file.path(CACHE, "rare"))
run("--bfile", PLINK, "--keep", file.path(CACHE, "universe.fam"),
    "--extract", file.path(CACHE, "common.ids"), "--indep-pairwise", "50kb", "1", "0.5",
    "--out", file.path(CACHE, "prune"))
run("--bfile", PLINK, "--keep", file.path(CACHE, "universe.fam"),
    "--extract", file.path(CACHE, "prune.prune.in"), "--make-bed",
    "--out", file.path(CACHE, "pruned"))

## --- rare markers as a carrier list ----------------------------------------
X <- BEDMatrix(file.path(CACHE, "rare"), simple_names = TRUE)
strains <- rownames(X)
bim <- fread(file.path(CACHE, "rare.bim"), header = FALSE,
             col.names = c("chrom", "id", "cm", "pos", "a1", "a2"))
msg(sprintf("rare set: %d markers x %d strains", ncol(X), nrow(X)))

parts <- split(seq_len(ncol(X)), ceiling(seq_len(ncol(X)) / CHUNK))
car <- rbindlist(lapply(parts, function(ix) {
  g <- X[, ix, drop = FALSE]
  w <- which(g > 0, arr.ind = TRUE)      # 0/2 coding; NA is not > 0
  data.table(m = ix[w[, 2]], s = w[, 1])
}))
setkey(car, m)
msg(sprintf("carrier pairs: %s", format(nrow(car), big.mark = ",")))

## carrier counts WITHIN the 571-strain universe; a marker can lose carriers to
## the universe restriction, and one with none left can never be private
nc <- car[, .N, by = m]
car <- car[nc[N >= 1L], on = "m"]

## --- classify each carrier pair by that strain's own divergent regions ------
bed <- fread(DIV, col.names = c("chrom", "start", "end", "strain"))
bed <- bed[strain %in% strains][, start := start + 1L]     # BED is 0-based
bed[, s := match(strain, strains)]
setkey(bed, s, chrom, start, end)

car[, `:=`(chrom = bim$chrom[m], pos = bim$pos[m])]
car[, `:=`(start = pos, end = pos)]
hit <- foverlaps(car, bed, by.x = c("s", "chrom", "start", "end"), nomatch = NULL)
car[, in_div := FALSE]
car[unique(hit[, .(s, chrom, pos = i.start)]), in_div := TRUE, on = .(s, chrom, pos)]
msg(sprintf("carrier pairs inside the carrier's own divergent regions: %.1f%%",
            100 * mean(car$in_div)))
car <- car[, .(m, s, in_div)]

## --- pruned common markers, and the one cross-product the search reuses -----
P <- BEDMatrix(file.path(CACHE, "pruned"), simple_names = TRUE)
stopifnot(identical(rownames(P), strains))
gp <- P[, ]
gp[is.na(gp)] <- 0                   # 0.7% missing; mean-imputing at 0 after
gp <- gp / 2                         # the 0/2 homozygous coding -> 0/1
storage.mode(gp) <- "double"
msg(sprintf("pruned set: %d markers x %d strains", ncol(gp), nrow(gp)))

## G = X'X over pruned markers, strains x strains. Any panel's marker-centred
## kinship is H G[S,S] H, so this is computed once and never again.
G <- tcrossprod(gp)                  # strains x strains
msg(sprintf("gram: %d x %d", nrow(G), ncol(G)))

saveRDS(list(strains = strains, car = car, n_rare = ncol(X),
             gp = gp, G = G, pruned_bim = fread(file.path(CACHE, "pruned.bim"),
               header = FALSE, col.names = c("chrom","id","cm","pos","a1","a2"))),
        file.path(CACHE, "precomputed.rds"), compress = FALSE)
msg("wrote ", file.path(CACHE, "precomputed.rds"))
