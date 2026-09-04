## sid-2 coding variation in the wild population ----------------------------
##
##   Rscript scripts/sid2_variant_table.R
##     -> data/structure/sid2_variants_cendr.tsv
##
## Builds the table behind the variant panel of Figure 4: every annotated
## protein-altering sid-2 variant, with its allele frequency across the CeNDR
## isotype set.
##
## WHAT IS INCLUDED
## Protein-altering changes only -- missense, stop and frameshift. Synonymous
## variants and the one splice-region variant are dropped, and so is
## everything non-coding.
##
## COVERAGE, WHICH IS THE IMPORTANT CAVEAT
## Amino-acid annotation comes from the bcftools csq (BCSQ) run in
## scripts/pooled_cross_candidate_variation.R, which annotated variants
## segregating between the parents of the two crosses -- N2, XZ1516, JU1793 and
## JU2466. It is therefore NOT an exhaustive catalogue of sid-2 coding
## variation in the species: a variant private to some other wild isolate is
## absent from this table because nothing here annotated it. 81 variants of any
## kind segregate in the 3 kb sid-2 span in CeNDR; this table carries the
## annotated protein-altering subset. The script reports both numbers so the
## figure can state the coverage rather than imply completeness.
##
## Allele frequencies are computed across the whole CeNDR 20210121 isotype set
## (540 isotypes), not across the four parents, so they are population
## frequencies for an incompletely ascertained variant set. Allele orientation
## is checked against the BCSQ reference and alternate alleles and the script
## stops if any variant disagrees.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

ANN   <- "plots/pooled_cross_intersection/TABLE_sid2_parental_variants.tsv"
LDOUT <- "data/structure/sid2_variant_ld.tsv"
PARENTS <- c("JU1793", "JU2466", "N2", "XZ1516")
PLINK <- "data/genotypes/CeNDR20210121_Plink/III"
OUT   <- "data/structure/sid2_variants_cendr.tsv"
SPAN  <- c(13679000L, 13682000L)     # the sid-2 span, from the annotation

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

ann <- read_tsv(ANN, show_col_types = FALSE) %>%
  filter(gene == "sid-2") %>%
  distinct(pos, ref, alt, consequence, aa_change) %>%
  arrange(pos)
msg("annotated sid-2 variants: ", nrow(ann), " (",
    paste(sort(unique(ann$consequence)), collapse = ", "), ")")

## --- allele frequencies across the CeNDR isotype set -----------------------
scratch <- file.path(tempdir(), "sid2_freq")
dir.create(scratch, showWarnings = FALSE, recursive = TRUE)
stem <- file.path(scratch, "sid2")
st <- system2("plink2",
              c("--bfile", PLINK, "--chr", "III",
                "--from-bp", SPAN[1], "--to-bp", SPAN[2],
                "--freq", "--out", stem,
                ## roman-numeral chromosomes are "extra" to plink2
                "--allow-extra-chr"),
              stdout = TRUE, stderr = TRUE)
if (!file.exists(paste0(stem, ".afreq"))) {
  cat(st, sep = "\n"); stop("plink2 produced no .afreq")
}
frq <- read_tsv(paste0(stem, ".afreq"), show_col_types = FALSE) %>%
  rename(chrom = `#CHROM`) %>%
  mutate(pos = as.integer(sub(".*:", "", ID)))
msg("variants of any kind in the span, CeNDR: ", nrow(frq))

## --- join, checking the allele orientation --------------------------------
d <- ann %>% inner_join(frq %>% select(pos, REF, ALT, ALT_FREQS, OBS_CT),
                        by = "pos")
stopifnot(nrow(d) == nrow(ann))
bad <- d %>% filter(ref != REF | alt != ALT)
if (nrow(bad)) {
  print(as.data.frame(bad)); stop("allele orientation disagrees with BCSQ")
}
msg("allele orientation matches BCSQ for all ", nrow(d), " variants")

## --- protein-altering only, with residue and a display label --------------
## aa_change looks like "96T>96K", or "151A>151I;151A>151T" when a variant hits
## more than one transcript. Parse the first record for the residue and the
## reference amino acid, then collapse the distinct alternate amino acids.
parse_aa <- function(s) {
  parts <- str_split(s, ";")[[1]]
  m <- str_match(parts, "^(\\d+)([A-Za-z*])>(\\d+)([A-Za-z*])$")
  if (any(is.na(m[, 1]))) return(tibble(residue = NA_integer_, label = NA_character_))
  tibble(residue = as.integer(m[1, 2]),
         label = paste0(m[1, 3], m[1, 2], "/",
                        paste(unique(m[, 5]), collapse = "/")) %>%
           str_replace("/", "") %>% str_replace("^(\\w)(\\d+)", "\\1\\2"))
}

v <- d %>%
  filter(consequence %in% c("missense", "stop_gained", "stop_lost",
                            "frameshift", "start_lost")) %>%
  rowwise() %>%
  mutate(residue = parse_aa(aa_change)$residue,
         ## label as refAA-position-altAA, alternates joined by "/"
         label = {
           m <- str_match(str_split(aa_change, ";")[[1]],
                          "^(\\d+)([A-Za-z*])>(\\d+)([A-Za-z*])$")
           paste0(m[1, 3], m[1, 2], paste(unique(m[, 5]), collapse = "/"))
         }) %>%
  ungroup() %>%
  transmute(chrom = "III", pos, ref, alt, consequence, aa_change,
            residue, label,
            af = ALT_FREQS, n_alleles = OBS_CT,
            n_isotypes = round(OBS_CT / 2)) %>%
  arrange(residue)

stopifnot(!any(is.na(v$residue)))
msg("protein-altering variants kept: ", nrow(v), " of ", nrow(d),
    " annotated | dropped ",
    paste(setdiff(d$consequence, v$consequence), collapse = ", "))

cat("\n== sid-2 protein-altering variants, CeNDR frequency ==\n")
print(as.data.frame(v %>% transmute(pos, label, residue,
        `AF (%)` = round(100 * af, 1), isotypes = n_isotypes)),
      row.names = FALSE)

## --- genotypes at these sites, and the LD among them ----------------------
## plink2 --export A names each column with the allele it counted, which is
## not always the BCSQ reference, so the counted allele is checked per site
## rather than assumed.
writeLines(paste0("III:", v$pos), paste0(stem, ".snps"))
st2 <- system2("plink2",
               c("--bfile", PLINK, "--extract", paste0(stem, ".snps"),
                 "--export", "A", "--out", paste0(stem, "_gt"),
                 "--allow-extra-chr"),
               stdout = TRUE, stderr = TRUE)
if (!file.exists(paste0(stem, "_gt.raw"))) {
  cat(st2, sep = "\n"); stop("plink2 produced no .raw")
}
raw <- read_tsv(paste0(stem, "_gt.raw"), show_col_types = FALSE)
cols <- grep("^III:", names(raw), value = TRUE)
gt <- as.matrix(raw[, cols]); rownames(gt) <- raw$IID
counted <- sub(".*_", "", cols)
gpos <- as.integer(sub("_.*", "", sub("III:", "", cols)))
stopifnot(setequal(gpos, v$pos))
gt <- gt[, order(match(gpos, v$pos))]
counted <- counted[order(match(gpos, v$pos))]
colnames(gt) <- v$label
## Dosage of the ALT allele, whichever allele plink happened to count.
## Direct column indexing, NOT sweep(..., FUN = ifelse): ifelse returns a
## result the length of its CONDITION, so a scalar per-column flag collapsed
## each column to a single value and silently destroyed the matrix.
flip <- counted != v$alt
alt_dose <- gt
if (any(flip)) alt_dose[, flip] <- 2 - gt[, flip]
stopifnot(identical(dim(alt_dose), dim(gt)),
          identical(colnames(alt_dose), v$label))
msg("  counted alleles: ", paste(counted, collapse = " "),
    " | flipped to ALT dosage where needed")

## LD among the protein-altering variants
R <- cor(alt_dose, use = "pairwise.complete.obs")
ld <- as_tibble(as.table(R), .name_repair = "minimal") %>%
  setNames(c("variant_a", "variant_b", "r")) %>%
  filter(variant_a < variant_b) %>%
  mutate(r2 = r^2) %>% arrange(desc(r2))
write_tsv(ld, LDOUT)
msg("wrote ", LDOUT)
cat("\n== LD among protein-altering sid-2 variants (top pairs) ==\n")
print(as.data.frame(ld %>% filter(r2 > 0.2) %>%
        transmute(variant_a, variant_b, r = round(r, 3), r2 = round(r2, 3))),
      row.names = FALSE)

## amino acid carried by each parent: ALT dosage 2 -> alt residue, 0 -> ref
aa_ref <- str_match(str_split(v$aa_change, ";") %>% map_chr(1),
                    "^\\d+([A-Za-z*])>")[, 2]
aa_alt <- v %>% pull(label) %>% str_replace("^[A-Za-z*]\\d+", "")
pget <- function(strain) {
  if (!strain %in% rownames(alt_dose)) return(rep(NA_character_, nrow(v)))
  d <- alt_dose[strain, ]
  ifelse(is.na(d), NA_character_, ifelse(d == 2, aa_alt, aa_ref))
}
for (sn in PARENTS) v[[paste0(tolower(sn), "_aa")]] <- pget(sn)
v$parents_differ <- v$ju1793_aa != v$ju2466_aa

cat("\n== allele carried by each cross parent ==\n")
print(as.data.frame(v %>% transmute(label, residue,
        `AF (%)` = round(100 * af, 1),
        JU1793 = ju1793_aa, JU2466 = ju2466_aa, N2 = n2_aa, XZ1516 = xz1516_aa,
        differ = parents_differ)), row.names = FALSE)
msg("JU1793 and JU2466 differ at ", sum(v$parents_differ), " of ", nrow(v),
    " protein-altering sites: ",
    paste(v$label[v$parents_differ], collapse = ", "))

writeLines(c(
  sprintf("# sid-2 protein-altering variants with CeNDR allele frequencies"),
  sprintf("# generated %s by scripts/sid2_variant_table.R", Sys.Date()),
  sprintf("# annotation: bcftools csq, variants segregating among N2, XZ1516, JU1793, JU2466"),
  sprintf("# NOT an exhaustive catalogue: %d variants of any kind segregate in the", nrow(frq)),
  sprintf("# %.1f kb sid-2 span in CeNDR; %d are annotated here and %d are protein-altering",
          diff(SPAN) / 1000, nrow(d), nrow(v)),
  sprintf("# frequencies are across the CeNDR 20210121 isotype set (%d isotypes)",
          max(v$n_isotypes))), OUT)
write_tsv(v, OUT, append = TRUE, col_names = TRUE)
msg("wrote ", OUT)
