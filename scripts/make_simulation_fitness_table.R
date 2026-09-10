## Stage the seven published traits the simulation used as fitness ------------
##
##   Rscript scripts/make_simulation_fitness_table.R
##     <SIM_TRAITS>/*.tsv, <SIM_TRAITS>/TableS4_GWAS-phenotypes.csv
##     -> supplemental_data/deconvolution/simulation_fitness_traits.tsv
##
## WHY THIS EXISTS
## METHODS.txt carried a [TO FILL] saying the simulation's fitness values were
## never archived, so the reported r-squared could not be recomputed -- they
## were read back out of text embedded in the 2021 per-trait PDFs by
## scripts/extract_sim_reported_r2.py. That is now fixable. The fitness was not
## a random draw at all: scripts/legacy/haploReg_original.R lines 309-334 read a
## directory of published traits with validated QTL and used the trait values
## THEMSELVES as fitness. Those files survive outside this repository.
##
## The seven traits arrive in six files, which is why the directory looks like
## it holds five:
##
##   ARSENIC_PC1.tsv               PC1
##   GWA_ascr5_dauer.tsv           assay_norm
##   albendazole_q75TOF.tsv        Albendazole_q75.TOF
##   mtdna_ratio.tsv               mtDNA_ratio
##   propionate_L1_survival.tsv    value
##   TableS4_GWAS-phenotypes.csv   amsacrine_f.L1 AND etoposide_median.TOF
##
## What is deposited is the raw full join on strain -- the published values, not
## the fitness. The two transforms the original applies (shift each trait by its
## own minimum so fitness is non-negative, then set NA to 0 so a strain with no
## published value for a trait is absent from that pool) are left to the
## consumer, scripts/simulation_recompute_r2.R, so that the deposit stays a
## record of the traits and the derivation stays visible in code.
##
## SOURCE IS OUTSIDE THIS REPOSITORY. Default location below; set SIM_TRAITS to
## override. The output is small enough to commit, so the recompute runs from a
## clone once this has been staged.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(purrr)
})

SRC <- Sys.getenv("SIM_TRAITS",
                  "/Users/Stefan/UCLA/Projects/bulkGWAS/traits_with_validated_qtl")
OUT <- "supplemental_data/deconvolution/simulation_fitness_traits.tsv"

TRAITS <- c("Albendazole_q75.TOF", "PC1", "assay_norm", "mtDNA_ratio",
            "value", "amsacrine_f.L1", "etoposide_median.TOF")

if (!dir.exists(SRC))
  stop("trait directory not found at '", SRC, "'.\n",
       "  These files are not in the repository. Set SIM_TRAITS to your copy.",
       call. = FALSE)

tsvs <- list.files(SRC, "\\.tsv$", full.names = TRUE)
csv  <- file.path(SRC, "TableS4_GWAS-phenotypes.csv")
if (!length(tsvs) || !file.exists(csv))
  stop("expected per-trait .tsv files and TableS4_GWAS-phenotypes.csv in ", SRC,
       call. = FALSE)

parts <- c(lapply(tsvs, read_tsv, show_col_types = FALSE),
           list(read_csv(csv, show_col_types = FALSE)))
phenos <- reduce(parts, full_join, by = "strain")

missing <- setdiff(TRAITS, names(phenos))
if (length(missing))
  stop("trait column(s) absent from the source files: ",
       paste(missing, collapse = ", "), call. = FALSE)

phenos <- phenos %>% select(strain, all_of(TRAITS)) %>% arrange(strain)

write_tsv(phenos, OUT)
message(sprintf("wrote %s: %d strains x %d traits", OUT, nrow(phenos), length(TRAITS)))
for (t in TRAITS)
  message(sprintf("  %-22s %3d with a value, %3d NA",
                  t, sum(!is.na(phenos[[t]])), sum(is.na(phenos[[t]]))))
