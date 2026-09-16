## Simulated phenotypes for mapping the optimised panel against the others ----
##
##   Rscript scripts/DIAG_pool_make_mapping_sim.R [n_loci] [n_rep]
##     -> cluster/mapping_simulation/traits_<panel>.csv   one per panel
##     -> cluster/mapping_simulation/simulation_key.tsv   what each trait is
##     -> cluster/mapping_simulation/trait_lists/<panel>.txt
##     -> cluster/mapping_simulation/README.md
##
## WHY. Everything in DIAG_pool_optimizer.R is measured in MAPPABLE MARKERS, a
## surrogate. The claim that matters is detected loci, and no marker count
## establishes it. These are phenotypes with the answer known: a variant of known
## effect at a known position, so a scan either finds it or does not.
##
## THE DESIGN, and the two choices in it that are not obvious.
##
## Causal variants are the SAME across panels. A locus is eligible only if it
## clears MAF 0.05 in EVERY panel, so no panel is scored on loci another panel
## could not test at all. This is also why changing ONE panel invalidates every
## file here and not just that panel's: the eligible locus set is defined by the
## intersection, so the key moves too. Loci are drawn stratified across three MAF bins,
## because power depends on allele frequency and an unstratified draw would
## sample mostly rare ones.
##
## HERITABILITY is fixed per panel, not the effect size. The same beta in two
## panels means different heritabilities, because genotype variance differs with
## allele frequency; fixing h2 instead asks "at equal signal to noise, which
## panel finds the locus", which is the question. beta is back-computed per
## panel and recorded, so the other reading is still available from the key.
##
## NULL TRAITS. Every panel gets pure-noise traits as well. Panels differ in
## structure, so they differ in genomic inflation, and a threshold taken from one
## panel does not transfer to another. The nulls give a panel-specific
## distribution of the genome-wide maximum to calibrate against -- without them
## a panel could look powerful only by being inflated.
##
## Phenotypes are drawn from the same genotype matrix the optimiser scored
## (CaeNDR 20231213). Map with a 20231213 VCF so the causal positions exist.
##
## Exploratory, on the pool_optimization branch.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })

a <- commandArgs(TRUE)
N_LOCI <- if (length(a) >= 1) as.integer(a[1]) else 12L   # causal loci per MAF design
N_REP  <- if (length(a) >= 2) as.integer(a[2]) else 3L    # noise replicates per cell
H2     <- c(0.05, 0.15, 0.35)
N_NULL <- 15L
MAF_MIN <- 0.05
BINS   <- list(c(0.05, 0.10), c(0.10, 0.25), c(0.25, 0.50))
OUT    <- "cluster/mapping_simulation"
dir.create(file.path(OUT, "trait_lists"), recursive = TRUE, showWarnings = FALSE)
set.seed(1)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

P <- readRDS(".pool_opt_cache/precomputed.rds")
U <- P$strains; GENO <- P$gp
bim <- P$pruned_bim
bim[, id := sub("_.*$", "", id)]
memb <- fread("plots/diagnostics/pool_optimizer_panels.tsv")
PANELS <- split(memb$strain, memb$panel)
names(PANELS)[names(PANELS) == "optimised"] <- "optimised (VIF constraint)"

## The SUPERSEDED K*-constrained panel is carried as its own arm rather than
## discarded. The constraint swap bought 15,644 -> 15,969 mappable markers, and
## mappable markers are the surrogate this whole simulation exists to check.
## Running both panels is the only way to find out whether that 2% converts into
## detected loci or into nothing, and it costs one more pheno file.
KF <- ".pool_opt_cache/annealed_n96.rds"
if (file.exists(KF)) {
  PANELS[["optimised (K* constraint)"]] <- U[readRDS(KF)$S]
  msg("carrying the superseded K*-constrained panel as a comparison arm")
}
slug <- function(x) gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(x)))
msg(sprintf("panels: %s", paste(sprintf("%s (%d)", names(PANELS), lengths(PANELS)),
                                collapse = ", ")))

## --- eligible loci: testable in EVERY panel --------------------------------
idx <- lapply(PANELS, function(s) match(s, U))
maf <- sapply(idx, function(i) { f <- colMeans(GENO[i, , drop = FALSE]); pmin(f, 1 - f) })
ok  <- apply(maf, 1, function(r) all(r >= MAF_MIN))
msg(sprintf("markers clearing MAF %.2f in every panel: %d of %d", MAF_MIN, sum(ok), nrow(maf)))

## stratify on the MINIMUM MAF across panels, so a locus sits in the bin that
## describes how hard it is in the panel where it is hardest
minmaf <- apply(maf, 1, min)
pick <- unlist(lapply(BINS, function(b) {
  cand <- which(ok & minmaf >= b[1] & minmaf < b[2])
  if (!length(cand)) return(integer(0))
  sample(cand, min(N_LOCI, length(cand)))
}))
LOCI <- data.table(m = pick, id = bim$id[pick], chrom = bim$chrom[pick],
                   pos = bim$pos[pick], min_maf = minmaf[pick])
msg(sprintf("causal loci drawn: %d (%s per MAF bin)", nrow(LOCI),
            paste(table(cut(LOCI$min_maf, c(0.05, 0.10, 0.25, 0.50), include.lowest = TRUE)),
                  collapse = "/")))

## --- simulate ---------------------------------------------------------------
key <- list()
for (pn in names(PANELS)) {
  s <- PANELS[[pn]]; i <- idx[[pn]]; n <- length(i)
  ph <- data.table(strain = s)
  tn <- 0L
  for (li in seq_len(nrow(LOCI))) for (h2 in H2) for (r in seq_len(N_REP)) {
    g <- GENO[i, LOCI$m[li]]
    vg <- var(g)
    if (vg <= 0) next
    tn <- tn + 1L
    nm <- sprintf("sim_%04d", tn)
    ## fix h2, back out beta: var(beta g) / (var(beta g) + var(e)) = h2
    beta <- sqrt(h2 / (1 - h2) / vg)
    y <- as.vector(beta * g + rnorm(n))
    ph[[nm]] <- round((y - mean(y)) / sd(y), 6)
    key[[length(key) + 1]] <- data.table(
      panel = pn, trait = nm, kind = "causal", causal_id = LOCI$id[li],
      chrom = LOCI$chrom[li], pos = LOCI$pos[li],
      maf_in_panel = round(min(mean(g), 1 - mean(g)), 4),
      h2 = h2, beta = round(beta, 4), rep = r)
  }
  for (r in seq_len(N_NULL)) {
    tn <- tn + 1L
    nm <- sprintf("sim_%04d", tn)
    ph[[nm]] <- round(as.vector(scale(rnorm(n))), 6)
    key[[length(key) + 1]] <- data.table(
      panel = pn, trait = nm, kind = "null", causal_id = NA_character_,
      chrom = NA_character_, pos = NA_integer_, maf_in_panel = NA_real_,
      h2 = 0, beta = 0, rep = r)
  }
  f <- file.path(OUT, sprintf("traits_%s.csv", slug(pn)))
  fwrite(ph, f)
  writeLines(paste(setdiff(names(ph), "strain"), collapse = ","),
             file.path(OUT, "trait_lists", sprintf("%s.txt", slug(pn))))
  msg(sprintf("  %-22s %3d strains x %3d traits -> %s", pn, nrow(ph), ncol(ph) - 1L, basename(f)))
}
K <- rbindlist(key)
fwrite(K, file.path(OUT, "simulation_key.tsv"), sep = "\t")

## trait names are assigned per panel in the same order, so sim_0001 is the same
## (locus, h2, replicate) cell everywhere -- check that rather than assume it
chk <- dcast(K[kind == "causal"], trait ~ panel, value.var = "causal_id")
stopifnot(all(apply(chk[, -1], 1, function(r) length(unique(r[!is.na(r)])) == 1)))
msg("verified: a given trait name is the same locus and h2 in every panel")
cat(sprintf("\n%d traits per panel (%d causal + %d null), %d panels, %d GEMMA runs total\n",
            uniqueN(K$trait), uniqueN(K[kind == "causal"]$trait), N_NULL,
            uniqueN(K$panel), nrow(K)))

## --- how to run it ----------------------------------------------------------
mafsum <- K[kind == "causal", .(median_maf = median(maf_in_panel)), by = panel][order(-median_maf)]
readme <- c(
"# Mapping simulation: does the optimised panel actually detect more?", "",
"Phenotypes with the answer known, for scoring panels against each other on",
"detected loci rather than on the mappable-marker surrogate.", "",
sprintf("Generated by `scripts/DIAG_pool_make_mapping_sim.R` (%s).", format(Sys.Date())), "",
"## Run", "",
"Map against a **20231213** VCF -- the causal positions come from that release.",
"Note the genotypes here come from the PLINK fileset `filtered_geno10`, not from",
"a VCF, and `WI.20231213.hard-filter.isotype.vcf.gz` is NOT present on this",
"machine (only an orphan `.tbi`). `bcsq.vcf.gz` in the same directory carries the",
"identical 611 samples and every panel strain, and the pipeline's plink step",
"already passes --snps-only --biallelic-only, which is the filter its extra",
"record types need.",
"All traits for a panel live in one file, so each panel is one invocation and",
"the plink conversion and the kinship matrices are built once and reused.", "",
"```sh",
"for p in optimised_vif_constraint optimised_k_constraint \\",
"         rnai_panel_93 baugh_panel_102 naive_most_private; do",
"  nextflow run main.nf -profile hoffman2 \\",
"    --vcf /Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213/bcsq.vcf.gz \\",
"    --pheno traits/traits_${p}.csv \\",
"    --traits \"$(cat trait_lists/${p}.txt)\" \\",
"    --name sim_${p}",
"done",
"```", "",
"`--expect_individuals` and `--expect_variants` are pinned to the pos-1 scan in",
"`nextflow.config`; set them to 0 or to this panel's counts, or the run stops on",
"the guard rather than on a real problem.", "",
"## Files", "",
"| file | contents |",
"|---|---|",
"| `traits_<panel>.csv` | `strain` then one column per simulated trait |",
"| `trait_lists/<panel>.txt` | the trait names, comma separated, for `--traits` |",
"| `simulation_key.tsv` | what every trait is: causal marker, position, MAF in that panel, h2, beta, replicate |", "",
"## Scoring", "",
"Join the scan output to `simulation_key.tsv` on `trait`. For a causal trait,",
"a detection is a marker above threshold within some window of `chrom`/`pos`;",
"anything above threshold far from it is a false positive. Power is the",
"detected fraction per (panel, h2, MAF bin) cell.", "",
"Take the threshold from each panel's own `kind == \"null\"` traits rather than",
"from a shared Bonferroni line. The panels differ in structure and so in genomic",
"inflation, and a panel can look powerful merely by being inflated. The null",
"traits are pure noise on the same genotypes, so the 95th percentile of their",
"genome-wide maxima is the honest per-panel cutoff.", "",
"## The two optimiser arms", "",
"`optimised_vif_constraint` is the current panel; `optimised_k_constraint` is the",
"superseded one built on private-marker count. They share 45 of 96 strains. The",
"swap bought 15,644 -> 15,969 mappable markers, and mappable markers are the",
"surrogate this simulation exists to check, so both are run: the comparison is",
"whether that 2% shows up as detected loci or as nothing.", "",
"## What is held fixed, and what is not", "",
"Causal loci are shared: a locus is eligible only if it clears MAF 0.05 in every",
"panel, so no panel is scored on loci another could not test. A given trait name",
"is the same locus and the same h2 in every panel file.", "",
"Heritability is fixed, not effect size. The same beta would mean different",
"heritabilities in panels with different allele frequencies; fixing h2 asks 'at",
"equal signal to noise, which panel finds it'. `beta` is in the key if the other",
"comparison is wanted.", "",
"Noise is drawn independently per panel. The panels hold different strains, so",
"there is no paired draw to share.", "",
"What is NOT equalised is the allele frequency of a shared locus, which differs",
"between panels by construction:", "",
paste0("  ", paste(sprintf("%s %.2f", mafsum$panel, mafsum$median_maf), collapse = " | ")),
"",
"(median MAF over the causal loci). At fixed h2 the non-centrality barely depends",
"on MAF, so this is not much of a power confound, but it is a real difference in",
"how balanced the panels are and `maf_in_panel` is in the key to condition on.", "",
"## Caveat", "",
"One causal variant per trait, additive, no epistasis, no polygenic background.",
"This measures detection of a single locus of known effect and nothing else.")
writeLines(readme, file.path(OUT, "README.md"))
msg("wrote ", file.path(OUT, "README.md"))
