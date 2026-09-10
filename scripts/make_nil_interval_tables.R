## Stage the gene content and parental variation of the 37 kb NIL interval ----
##
##   Rscript scripts/make_nil_interval_tables.R
##     <BCSQ VCF>, <WS283 csq GFF3>
##     -> supplemental_data/mapping/nil_interval_genes.tsv
##        supplemental_data/mapping/nil_interval_exons.tsv
##        supplemental_data/mapping/nil_interval_parent_variants.tsv
##
## WHY THIS EXISTS
## Figure 3 resolves the chromosome III QTL to 13.6577-13.6950 Mb but does not
## say what is inside it. The question a reader asks next is which genes those
## 37 kb hold and which of them differ between the two parents of the cross in
## a way that could change a protein. Answering it needs the BCSQ-annotated
## CeNDR VCF (10 GB) and a WormBase GFF3, neither of which can live in the
## repository -- so the answer is extracted here into three small tables and
## the figure script reads only those.
##
## SOURCES ARE OUTSIDE THIS REPOSITORY. Defaults below; override with
## CENDR_BCSQ and WS_GFF3. Requires bcftools and tabix on PATH.
##
## IMPACT CLASSES. bcftools csq emits consequence terms, not the HIGH/MODERATE
## impact ladder that snpEff and VEP use, so the mapping is made explicit here
## rather than assumed:
##
##   HIGH      stop_gained, stop_lost, start_lost, frameshift,
##             splice_acceptor, splice_donor
##   MODERATE  missense, inframe_insertion, inframe_deletion
##   LOW       synonymous, splice_region, start_retained, stop_retained
##   MODIFIER  everything else -- intron, UTR, non-coding, unannotated
##
## Every differing site is kept, not only the protein-altering ones, because
## the denominator is the point: "two missense changes" means little without
## the 25 other differences it is two of.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(data.table); library(dplyr); library(readr); library(stringr)})

VCF  <- Sys.getenv("CENDR_BCSQ", "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213/bcsq.vcf.gz")
GFF  <- Sys.getenv("WS_GFF3",    "/Users/Stefan/UCLA/Genomics_Data/Annotations/c_elegans.PRJNA13758.WS283.csq.gff3.gz")
DIV  <- Sys.getenv("CENDR_DIVERGENT",
                   "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213/20231213_c_elegans_divergent_regions_strain.bed")
OUT  <- "supplemental_data/mapping"
CHR  <- "III"; FROM <- 13657700L; TO <- 13695000L     # Figure3_common.R RESOLVED
PARENTS <- c("JU1793", "JU2466")
REGION  <- sprintf("%s:%d-%d", CHR, FROM, TO)

for (f in c(VCF, GFF)) if (!file.exists(f)) stop("missing source: ", f, call. = FALSE)
for (b in c("bcftools", "tabix")) if (!nzchar(Sys.which(b))) stop("need ", b, " on PATH", call. = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## ---------------------------------------------------------------------------
## gene models
## ---------------------------------------------------------------------------
gff_attr <- function(x, key) {
  m <- str_match(x, paste0("(?:^|;)", key, "=([^;]*)"))[, 2]
  m
}
raw <- fread(cmd = sprintf("tabix %s %s", shQuote(GFF), REGION), sep = "\t", header = FALSE,
             col.names = c("chrom","source","feature","start","end","score","strand","frame","attr"))
msg("GFF rows in the interval: ", nrow(raw))

genes <- raw %>% filter(feature == "gene") %>%
  transmute(wbgene = gff_attr(attr, "Name"),
            locus  = gff_attr(attr, "locus"),
            sequence.name = gff_attr(attr, "sequence_name"),
            biotype = gff_attr(attr, "biotype"),
            strand, start, end) %>%
  mutate(label = coalesce(locus, sequence.name)) %>%
  arrange(start)

## the locus name often sits on the mRNA rather than the gene, so fill from there
mrna <- raw %>% filter(feature == "mRNA") %>%
  transmute(wbgene = sub("^gene:", "", gff_attr(attr, "Parent")),
            mrna.locus = gff_attr(attr, "locus"),
            transcript = gff_attr(attr, "Name"))
genes <- genes %>% left_join(mrna %>% distinct(wbgene, mrna.locus), by = "wbgene") %>%
  mutate(label = coalesce(locus, mrna.locus, sequence.name)) %>%
  select(-mrna.locus)

## exons, tagged with the gene they belong to via their transcript
tx2gene <- mrna %>% select(transcript, wbgene)
exons <- raw %>% filter(feature == "exon") %>%
  transmute(transcript = sub("^transcript:", "", gff_attr(attr, "Parent")),
            strand, start, end) %>%
  left_join(tx2gene, by = "transcript") %>%
  left_join(genes %>% select(wbgene, label), by = "wbgene") %>%
  filter(!is.na(wbgene)) %>% arrange(start)

write_tsv(genes, file.path(OUT, "nil_interval_genes.tsv"))
write_tsv(exons, file.path(OUT, "nil_interval_exons.tsv"))
msg("genes: ", nrow(genes), " (", sum(genes$biotype == "protein_coding"), " protein-coding)",
    "; exon rows: ", nrow(exons))

## ---------------------------------------------------------------------------
## the two parents' variation
## ---------------------------------------------------------------------------
q <- sprintf("bcftools query -r %s -s %s -f '%%POS\\t%%REF\\t%%ALT\\t[%%GT\\t]%%INFO/BCSQ\\n' %s",
             REGION, paste(PARENTS, collapse = ","), shQuote(VCF))
v <- fread(cmd = q, sep = "\t", header = FALSE, na.strings = "",
           col.names = c("pos","ref","alt","gt.JU1793","gt.JU2466","bcsq"))
msg("sites in the interval: ", nrow(v))

norm_gt <- function(g) gsub("\\|", "/", g)
v <- v %>% mutate(across(starts_with("gt."), norm_gt),
                  called = gt.JU1793 != "./." & gt.JU2466 != "./.",
                  differs = called & gt.JU1793 != gt.JU2466)

HIGH <- "stop_gained|stop_lost|start_lost|frameshift|splice_acceptor|splice_donor"
MOD  <- "missense|inframe"
LOW  <- "synonymous|splice_region|start_retained|stop_retained"

## BCSQ packs consequence|gene|transcript|biotype[|strand|aa|nt] and can carry
## several comma-separated entries; take the most severe, and keep its gene
sev <- function(s) {
  if (is.na(s) || s == ".") return("MODIFIER")
  if (grepl(HIGH, s)) "HIGH" else if (grepl(MOD, s)) "MODERATE"
  else if (grepl(LOW, s)) "LOW" else "MODIFIER"
}
first_field <- function(s, i) {
  if (is.na(s) || s == ".") return(NA_character_)
  parts <- strsplit(strsplit(s, ",", fixed = TRUE)[[1]][1], "|", fixed = TRUE)[[1]]
  if (length(parts) >= i) parts[i] else NA_character_
}
var <- v %>% filter(differs) %>%
  rowwise() %>%
  mutate(impact      = sev(bcsq),
         consequence = first_field(bcsq, 1),
         gene        = first_field(bcsq, 2),
         transcript  = first_field(bcsq, 3),
         aa.change   = first_field(bcsq, 6)) %>%
  ungroup() %>%
  mutate(consequence = ifelse(is.na(consequence), "unannotated", consequence),
         impact = factor(impact, levels = c("HIGH","MODERATE","LOW","MODIFIER"))) %>%
  select(pos, ref, alt, gt.JU1793, gt.JU2466, impact, consequence, gene,
         transcript, aa.change) %>%
  arrange(pos)
write_tsv(var, file.path(OUT, "nil_interval_parent_variants.tsv"))

cat("\n== sites where the parents differ, by impact ==\n")
print(as.data.frame(var %>% count(impact, .drop = FALSE)), row.names = FALSE)
cat("\n== the protein-altering ones ==\n")
print(as.data.frame(var %>% filter(impact %in% c("HIGH","MODERATE")) %>%
                      select(pos, ref, alt, impact, consequence, gene, aa.change)),
      row.names = FALSE)

## ---------------------------------------------------------------------------
## the two caveats that decide how strong the claim can be
## ---------------------------------------------------------------------------
cat(sprintf("\nno-call in either parent: %d of %d sites\n", sum(!v$called), nrow(v)))
if (file.exists(DIV)) {
  d <- fread(DIV, header = FALSE, select = 1:4,
             col.names = c("chrom","start","end","strain"))
  hit <- d %>% filter(chrom == CHR, end > FROM, start < TO, strain %in% PARENTS)
  cat("divergent regions overlapping the interval in either parent: ", nrow(hit), "\n", sep = "")
  if (nrow(hit)) print(as.data.frame(hit), row.names = FALSE)
} else cat("divergent-region BED not found; overlap not checked\n")
msg("wrote three tables to ", OUT)
