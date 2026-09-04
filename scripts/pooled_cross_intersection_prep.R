## pooled_cross_intersection_prep.R ----------------------------------------
## Builds the joined data bundle behind the pooled-GWAS x cross-mapping
## figures. Run once; pooled_cross_intersection_figures.R reads the cache.
##
##   Rscript scripts/pooled_cross_intersection_prep.R
##
## Inputs
##   data/pooled_RNAi_expt/data/processed_strain_fq.rda        pooled strain frequencies
##   data/pooled_RNAi_expt/reanalysis/mapping/*_loco_results.csv.gz   GEMMA LOCO GWAS
##   data/pooled_RNAi_expt/reanalysis/mapping/plink/traits.{gen,sample}  panel genotypes
##   data/cross_experiments/N2-XZ_export/                      N2 x XZ1516 xQTL
##   data/cross_experiments/JU1793-JU2466_export/              JU1793 x JU2466 xQTL
##   <CENDR>/Caenorhabditis_elegans.WBcel235.112.gff3          gene models
##   <CENDR>/bcsq.vcf.gz                                       parental genotypes
##
## Output
##   data/pooled_cross_intersection/bundle.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

## -- paths ------------------------------------------------------------------
ROOT     <- "."
POOLED   <- file.path(ROOT, "data/pooled_RNAi_expt")
MAPDIR   <- file.path(POOLED, "reanalysis/mapping")
NXZ      <- file.path(ROOT, "data/cross_experiments/N2-XZ_export")
JUX      <- file.path(ROOT, "data/cross_experiments/JU1793-JU2466_export")
CENDR    <- "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213"
OUT      <- file.path(ROOT, "data/pooled_cross_intersection")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

CHROMS <- c("I", "II", "III", "IV", "V", "X")
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## The two crosses, and which pooled trait each RNAi condition corresponds to.
CROSSES <- tibble::tribble(
  ~cross,            ~dir, ~p1,      ~p2,       ~phase_to,
  "N2xXZ1516",       NXZ,  "N2",     "XZ1516",  "N2",
  "JU1793xJU2466",   JUX,  "JU1793", "JU2466",  "JU1793")

## ===========================================================================
## 1. Pooled per-strain phenotypes
## ===========================================================================
## Convention taken from data/pooled_RNAi_expt/scripts/20240328_calculate_frequencies.R
## ("t2 - t2c" block): the ctrl baseline is the MEAN of ctrl's frequency across
## its replicates per strain x sample_type, subtracted from every individual RNAi
## replicate; the sign is RNAi - ctrl. Metadata fixes follow
## scripts/20260813_delta_freq_correlations.R.
msg("pooled phenotypes")
e <- new.env()
load(file.path(POOLED, "data/processed_strain_fq.rda"), envir = e)

fq <- e$processed_strain_frq %>%
  ## S19 (ctrl pellet) and S20 (ric-3 supernatant) were flagged bad and replaced
  ## by the "_b" re-preps S54/S55; excluded in the legacy script too.
  filter(!sample %in% c("S19", "S20")) %>%
  ## S28 is the true rde-3 t2 supernatant sample, orphaned by a SAMPLE NUMBER
  ## collision with S29 in experiment_meta.tsv.
  mutate(Gene        = ifelse(sample == "S28", "rde-3",       Gene),
         sample_type = ifelse(sample == "S28", "supernatant", sample_type),
         Timepoint   = ifelse(sample == "S28", 2,             Timepoint)) %>%
  filter(!(sample == "S29" & Gene == "rde-3")) %>%
  ## S48 (mig-6 supernatant) failed NNLS deconvolution -- all-NA frq.
  filter(sample != "S48")

ctrl_base <- fq %>%
  filter(Gene == "ctrl", Timepoint == 2) %>%
  group_by(strain, sample_type) %>%
  summarise(ctrlfrq = mean(frq, na.rm = TRUE), .groups = "drop")

pheno <- fq %>%
  filter(Timepoint == 2, Gene != "ctrl") %>%
  left_join(ctrl_base, by = c("strain", "sample_type")) %>%
  mutate(delta_fq = frq - ctrlfrq) %>%
  group_by(strain, Gene, sample_type) %>%
  summarise(pheno = mean(delta_fq, na.rm = TRUE),
            n_rep = sum(!is.na(delta_fq)), .groups = "drop")

msg("  ", n_distinct(pheno$strain), " strains x ", n_distinct(pheno$Gene), " RNAi targets")

## ===========================================================================
## 2. Pooled GWAS (GEMMA LOCO)
## ===========================================================================
msg("pooled GWAS")
gwas <- map_dfr(c("mig-6", "pos-1"), function(g) {
  fread(file.path(MAPDIR, sprintf("vst_ctrl_%s_T2_loco_results.csv.gz", g)),
        select = c("chr", "rs", "ps", "af", "beta", "se", "p_wald")) %>%
    as_tibble() %>% mutate(gene = g, .before = 1)
}) %>%
  filter(chr %in% CHROMS) %>%
  mutate(neglog10p = -log10(p_wald))

## Bonferroni over the tested markers. Reported alongside the eigen-decomposition
## style threshold the cross pipeline uses so the two scales are explicit.
N_MARKER <- gwas %>% filter(gene == "mig-6") %>% nrow()
GWAS_BF  <- -log10(0.05 / N_MARKER)
msg("  ", N_MARKER, " markers; Bonferroni -log10p = ", round(GWAS_BF, 3))

gwas_peaks <- gwas %>%
  group_by(gene, chr) %>% slice_max(neglog10p, n = 1, with_ties = FALSE) %>%
  ungroup() %>% mutate(significant = neglog10p > GWAS_BF) %>%
  arrange(gene, desc(neglog10p))

## ===========================================================================
## 3. Cross scans + support intervals (as shipped by each export)
## ===========================================================================
msg("cross scans")
read_scan <- function(path) {
  fread(path, select = c("chrom", "physical.position", "contrast.beta", "z",
                         "LOD", "neglog10p")) %>%
    as_tibble() %>% filter(chrom %in% CHROMS)
}

scan_files <- CROSSES %>%
  mutate(f = map(dir, ~list.files(file.path(.x, "plot_data"),
                                  pattern = "_plot_DF\\.tsv\\.gz$", full.names = TRUE))) %>%
  unnest(f) %>%
  mutate(contrast = f %>% basename() %>%
           sub("^.*_contrast_", "", .) %>% sub("_[0-9]+_plot_DF\\.tsv\\.gz$", "", .),
         ## "HT115g-MIG6g" and "HT115-mig6" name the same comparison
         label = contrast %>% tolower() %>% gsub("g(-|$)", "\\1", .) %>%
           gsub("-", " vs ", ., fixed = TRUE))

scans <- set_names(map(scan_files$f, read_scan),
                   paste0(scan_files$cross, " | ", scan_files$label))
scan_meta <- scan_files %>% select(cross, contrast, label) %>%
  mutate(key = names(scans), .before = 1)

## The exports ship intervals computed with the same drop rule (10% of each
## chromosome's own peak), so they are read rather than recomputed.
ivs <- CROSSES %>%
  mutate(f = map(dir, ~list.files(file.path(.x, "intervals"),
                                  pattern = "_intervals\\.tsv$", full.names = TRUE))) %>%
  unnest(f) %>%
  mutate(d = map(f, ~as_tibble(fread(.x)))) %>%
  unnest(d) %>%
  mutate(contrast = f %>% basename() %>%
           sub("^.*_contrast_", "", .) %>% sub("_[0-9]+_intervals\\.tsv$", "", .),
         label = contrast %>% tolower() %>% gsub("g(-|$)", "\\1", .) %>%
           gsub("-", " vs ", ., fixed = TRUE),
         key = paste0(cross, " | ", label)) %>%
  select(cross, contrast, label, key, chrom, peak.position, peak.LOD,
         lod.drop.used, lcon, rcon, marker, width.kb, significant)

## Threshold used by the cross pipeline (alpha = 0.05 over 2000 effective tests,
## "package" LOD convention -- see N2-XZ_export/scripts/config.yml).
CROSS_THR <- stats::qchisq(2 * 0.05 / 2000, df = 1, lower.tail = FALSE) / (2 * log(10))
msg("  ", length(scans), " contrasts; cross LOD threshold = ", round(CROSS_THR, 3))

## ===========================================================================
## 4. Gene models + RNAi machinery
## ===========================================================================
msg("gene annotation")
gff <- file.path(CENDR, "Caenorhabditis_elegans.WBcel235.112.gff3")
genes <- fread(cmd = sprintf(
  paste0("awk -F'\\t' '$3==\"gene\" && $1!~/^#/ { n=\".\"; ",
         "if (match($9,/Name=[^;]+/)) n=substr($9,RSTART+5,RLENGTH-5); ",
         "b=\".\"; if (match($9,/biotype=[^;]+/)) b=substr($9,RSTART+8,RLENGTH-8); ",
         "print $1\"\\t\"$4\"\\t\"$5\"\\t\"n\"\\t\"b }' %s"), shQuote(gff)),
  col.names = c("chrom", "start", "end", "gene", "biotype")) %>%
  as_tibble() %>% filter(chrom %in% CHROMS)

## Canonical dsRNA uptake / RNAi machinery genes, plus the RNAi targets assayed.
RNAI_GENES <- c(
  "sid-1", "sid-2", "sid-3", "sid-5",
  "rde-1", "rde-2", "rde-3", "rde-4", "rde-8", "rde-10", "rde-11", "rde-12",
  "mut-2", "mut-7", "mut-14", "mut-15", "mut-16",
  "dcr-1", "drh-1", "drh-3", "ergo-1", "alg-1", "alg-2",
  "sago-1", "sago-2", "ppw-1", "ppw-2",
  "nrde-1", "nrde-2", "nrde-3", "nrde-4", "hrde-1",
  "eri-1", "eri-6", "rrf-1", "rrf-2", "rrf-3",
  "rsd-2", "rsd-3", "rsd-6", "tsn-1")
TARGET_GENES <- c("pos-1", "mig-6", "vha-5", "rpn-12", "par-1")

rnai_genes <- genes %>%
  filter(gene %in% c(RNAI_GENES, TARGET_GENES)) %>%
  mutate(class = case_when(gene %in% c("sid-1", "sid-2", "sid-3", "sid-5") ~ "dsRNA uptake",
                           gene %in% TARGET_GENES                          ~ "RNAi target",
                           TRUE                                            ~ "RNAi machinery")) %>%
  arrange(chrom, start)

## ===========================================================================
## 5. Panel genotypes at loci of interest
## ===========================================================================
msg("panel genotypes at GWAS peaks")
panel_strains <- readLines(file.path(MAPDIR, "plink/traits.sample"))[-(1:2)] %>%
  str_split(" ") %>% map_chr(1)

## traits.gen is Oxford format: chrom rsid pos alleleA alleleB then 3 dosages
## per sample. alleleA is GEMMA's allele1 (the one `af` and `beta` refer to).
gen_at <- function(chr, ps) {
  ln <- system(sprintf("grep -m1 '^%s %s:%d ' %s", chr, chr, ps,
                       shQuote(file.path(MAPDIR, "plink/traits.gen"))), intern = TRUE)
  if (!length(ln)) return(NULL)
  f <- str_split(ln, " ")[[1]]
  m <- matrix(as.numeric(f[-(1:5)]), nrow = 3)
  tibble(strain   = panel_strains,
         genotype = apply(m, 2, function(v)
           if (all(v == 0)) NA_character_ else c(f[4], "het", f[5])[which.max(v)]),
         allele1 = f[4], allele0 = f[5])
}

gwas_sig <- gwas_peaks %>% filter(significant)
panel_gt <- pmap_dfr(gwas_sig %>% select(gene, chr, ps),
                     function(gene, chr, ps) {
                       g <- gen_at(chr, ps)
                       if (is.null(g)) return(NULL)
                       g %>% mutate(gwas_gene = gene, chrom = chr, ps = ps, .before = 1)
                     })

## ===========================================================================
## 6. Parental genotypes + cross informativeness
## ===========================================================================
## Whether a cross can test a locus at all depends on its two parents carrying
## different alleles there. Extracted once from the CeNDR isotype VCF; see the
## sibling shell step in the session notes for how parents_gt.tsv.gz is made.
msg("parental genotypes / cross informativeness")
PGT <- file.path(OUT, "parents_gt.tsv.gz")
if (!file.exists(PGT)) {
  vcf <- file.path(CENDR, "bcsq.vcf.gz")
  stopifnot(file.exists(vcf))
  msg("  extracting from ", basename(vcf), " (~4 min)")
  system(sprintf(
    "bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT]\\n' -s N2,XZ1516,JU1793,JU2466 %s | gzip -1 > %s",
    shQuote(vcf), shQuote(PGT)))
}
pgt <- fread(PGT, col.names = c("chrom", "ps", "ref", "alt", "N2", "XZ1516", "JU1793", "JU2466"))
pgt <- pgt[chrom %in% CHROMS]

hard <- function(x) fifelse(x == "0/0", 0L, fifelse(x == "1/1", 1L, NA_integer_))
pgt[, `:=`(g_N2 = hard(N2), g_XZ = hard(XZ1516), g_J93 = hard(JU1793), g_J66 = hard(JU2466))]
pgt[, `:=`(d_nxz = fifelse(!is.na(g_N2) & !is.na(g_XZ),  as.integer(g_N2  != g_XZ),  NA_integer_),
           d_jux = fifelse(!is.na(g_J93) & !is.na(g_J66), as.integer(g_J93 != g_J66), NA_integer_))]

BIN <- 1e5
pgt[, bin := floor(ps / BIN) * BIN]
inform <- pgt[, .(n_nxz = sum(!is.na(d_nxz)), k_nxz = sum(d_nxz, na.rm = TRUE),
                  n_jux = sum(!is.na(d_jux)), k_jux = sum(d_jux, na.rm = TRUE)),
              by = .(chrom, bin)] %>%
  as_tibble() %>%
  mutate(N2xXZ1516     = ifelse(n_nxz > 0, k_nxz / n_nxz, NA_real_),
         JU1793xJU2466 = ifelse(n_jux > 0, k_jux / n_jux, NA_real_)) %>%
  select(chrom, bin, n_nxz, n_jux, N2xXZ1516, JU1793xJU2466) %>%
  pivot_longer(c(N2xXZ1516, JU1793xJU2466), names_to = "cross", values_to = "discordance")

## parental genotypes at the significant pooled QTL, for the cross-design panel
parent_gt_qtl <- gwas_sig %>%
  select(gwas_gene = gene, chrom = chr, ps, af, beta) %>%
  left_join(pgt[, .(chrom, ps, ref, alt, N2, XZ1516, JU1793, JU2466)],
            by = c("chrom", "ps"))

## ===========================================================================
## 7. Two-way intersection
## ===========================================================================
msg("intersection")

## A. cross LOD at each significant pooled GWAS peak
isect_A <- pmap_dfr(gwas_sig %>% select(gene, chr, ps), function(gene, chr, ps) {
  imap_dfr(scans, function(s, key) {
    d  <- s %>% filter(chrom == chr, is.finite(LOD))
    if (!nrow(d)) return(NULL)
    j  <- which.min(abs(d$physical.position - ps))
    iv <- ivs %>% filter(key == !!key, chrom == chr)
    tibble(gwas_gene = gene, chrom = chr, ps = ps, key = key,
           cross_LOD_at_peak   = d$LOD[j],
           marker_dist_kb      = abs(d$physical.position[j] - ps) / 1e3,
           chrom_peak_LOD      = if (nrow(iv)) iv$peak.LOD else NA_real_,
           chrom_peak_position = if (nrow(iv)) iv$peak.position else NA_real_,
           inside_interval     = if (nrow(iv)) ps >= iv$lcon & ps <= iv$rcon else NA)
  })
}) %>% left_join(scan_meta, by = "key")

## B. best pooled GWAS signal inside each significant cross support interval
gwas_dt <- as.data.table(gwas)
isect_B <- ivs %>% filter(significant) %>%
  pmap_dfr(function(...) {
    iv <- tibble(...)
    map_dfr(c("mig-6", "pos-1"), function(g) {
      sub <- gwas_dt[gene == g & chr == iv$chrom & ps >= iv$lcon & ps <= iv$rcon]
      tibble(key = iv$key, cross = iv$cross, label = iv$label, chrom = iv$chrom,
             marker = iv$marker, cross_LOD = iv$peak.LOD, gwas_gene = g,
             n_snp = nrow(sub),
             best_neglog10p = if (nrow(sub)) max(sub$neglog10p) else NA_real_,
             best_ps        = if (nrow(sub)) sub$ps[which.max(sub$neglog10p)] else NA_integer_,
             gwas_sig       = if (nrow(sub)) max(sub$neglog10p) > GWAS_BF else NA)
    })
  })

## C. RNAi genes under each significant cross interval
isect_C <- ivs %>% filter(significant) %>%
  pmap_dfr(function(...) {
    iv <- tibble(...)
    h  <- rnai_genes %>% filter(chrom == iv$chrom, end >= iv$lcon, start <= iv$rcon)
    if (!nrow(h)) return(NULL)
    h %>% transmute(key = iv$key, cross = iv$cross, label = iv$label, chrom,
                    marker = iv$marker, cross_LOD = iv$peak.LOD,
                    peak.position = iv$peak.position,
                    gene, class, gene_start = start,
                    dist_to_peak_kb = pmin(abs(start - iv$peak.position),
                                           abs(end   - iv$peak.position)) / 1e3)
  })

## ===========================================================================
## 8. Save
## ===========================================================================
bundle <- list(
  pheno = pheno, gwas = gwas, gwas_peaks = gwas_peaks, gwas_sig = gwas_sig,
  GWAS_BF = GWAS_BF, N_MARKER = N_MARKER,
  scans = scans, scan_meta = scan_meta, ivs = ivs, CROSS_THR = CROSS_THR,
  CROSSES = CROSSES, genes = genes, rnai_genes = rnai_genes,
  panel_gt = panel_gt, parent_gt_qtl = parent_gt_qtl, inform = inform,
  isect_A = isect_A, isect_B = isect_B, isect_C = isect_C,
  built = Sys.time())

saveRDS(bundle, file.path(OUT, "bundle.rds"))
msg("wrote ", file.path(OUT, "bundle.rds"))

## -- console summary --------------------------------------------------------
cat("\n== significant pooled GWAS QTL ==\n")
print(as.data.frame(gwas_sig %>% select(gene, chr, ps, af, beta, neglog10p)),
      row.names = FALSE, digits = 4)

cat("\n== parental genotypes there (REF = N2 allele) ==\n")
print(as.data.frame(parent_gt_qtl), row.names = FALSE, digits = 4)

cat("\n== cross intervals containing a genome-wide significant pooled hit ==\n")
print(as.data.frame(isect_B %>% filter(gwas_sig) %>%
                      select(cross, label, marker, cross_LOD, gwas_gene,
                             best_neglog10p, best_ps)), row.names = FALSE, digits = 4)

cat("\n== dsRNA-uptake genes under cross intervals ==\n")
print(as.data.frame(isect_C %>% filter(class == "dsRNA uptake") %>%
                      arrange(desc(cross_LOD)) %>%
                      select(cross, label, marker, cross_LOD, gene, dist_to_peak_kb)),
      row.names = FALSE, digits = 4)
