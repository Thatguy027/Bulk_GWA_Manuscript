## Stage gene content and parental variation at the top mig-6 cross QTL ------
##
##   Rscript scripts/make_mig6_locus_tables.R
##     -> supplemental_data/mapping/mig6_locus_{genes,exons,variants,divergent}.tsv
##
## The same census scripts/make_nil_interval_tables.R runs on the 37 kb NIL
## interval, applied to every independent HT115-vs-mig-6 cross QTL above
## LOD 500, in a 100 kb window centred on each peak.
##
## WHICH LOCI. plots/TABLE_cross_qtl_full.tsv, contrast "ht115 vs mig6",
## peak.LOD > 500, restricted to peak.rank == 1 & separated. The rank filter
## matters: eight peaks clear LOD 500 but four are secondary peaks that
## cross_qtl_full_summary.R marks unseparated -- shoulders of one sweep rather
## than independent QTL -- and a 100 kb window on a shoulder is a window on
## the same locus twice.
##
## DIVERGENT REGIONS are staged alongside, per parent, because in these windows
## they are not a footnote: a divergent region is where a wild isolate's
## short reads stop aligning to N2 well enough to call, so an apparent absence
## of coding variation inside one is uninformative rather than reassuring.
##
## Sources are outside the repository -- CENDR_BCSQ, WS_GFF3, CENDR_DIVERGENT.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(data.table); library(dplyr); library(readr); library(stringr); library(purrr)})

VCF <- Sys.getenv("CENDR_BCSQ", "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213/bcsq.vcf.gz")
GFF <- Sys.getenv("WS_GFF3",    "/Users/Stefan/UCLA/Genomics_Data/Annotations/c_elegans.PRJNA13758.WS283.csq.gff3.gz")
DIV <- Sys.getenv("CENDR_DIVERGENT",
                  "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213/20231213_c_elegans_divergent_regions_strain.bed")
OUT   <- "supplemental_data/mapping"
WIN   <- 100000L                      # total window width, centred on the peak
LODMIN <- 500
PARENTS <- list(N2xXZ1516 = c("N2","XZ1516"), JU1793xJU2466 = c("JU1793","JU2466"))
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
for (f in c(VCF, GFF, DIV)) if (!file.exists(f)) stop("missing source: ", f, call. = FALSE)

loci <- read_tsv("plots/TABLE_cross_qtl_full.tsv", show_col_types = FALSE) %>%
  filter(contrast == "ht115 vs mig6", peak.LOD > LODMIN, peak.rank == 1, separated) %>%
  transmute(cross, chrom = as.character(chrom), peak.Mb, peak.LOD, dfreq,
            peak.bp = round(peak.Mb * 1e6),
            from = pmax(1L, as.integer(round(peak.Mb * 1e6) - WIN/2)),
            to   = as.integer(round(peak.Mb * 1e6) + WIN/2),
            locus = sprintf("%s %s:%.2f Mb", cross, chrom, peak.Mb)) %>%
  arrange(desc(peak.LOD))
msg("loci: ", nrow(loci))
print(as.data.frame(loci %>% select(locus, peak.LOD, from, to)), row.names = FALSE)

gff_attr <- function(x, key) str_match(x, paste0("(?:^|;)", key, "=([^;]*)"))[, 2]
HIGH <- "stop_gained|stop_lost|start_lost|frameshift|splice_acceptor|splice_donor"
MOD  <- "missense|inframe"
LOW  <- "synonymous|splice_region|start_retained|stop_retained"
## BCSQ can carry several comma-separated consequences for one site. Scanning
## the whole string for severity but reporting fields from the FIRST entry
## mismatches them -- a site whose first entry is 5_prime_utr but which also
## carries an inframe indel came out MODERATE while being labelled
## "5_prime_utr". Pick the most severe ENTRY, then read its fields.
rank1 <- function(e) if (grepl(HIGH, e)) 1L else if (grepl(MOD, e)) 2L else
  if (grepl(LOW, e)) 3L else 4L
IMPACTS <- c("HIGH","MODERATE","LOW","MODIFIER")
worst_entry <- function(s) {
  if (is.na(s) || s == ".") return(NA_character_)
  e <- strsplit(s, ",", fixed = TRUE)[[1]]
  e <- e[nzchar(e)]
  if (!length(e)) return(NA_character_)
  e[which.min(vapply(e, rank1, 1L))]
}
sev <- function(s) {
  e <- worst_entry(s)
  if (is.na(e)) "MODIFIER" else IMPACTS[rank1(e)]
}
fld <- function(s, i) {
  e <- worst_entry(s)
  if (is.na(e)) return(NA_character_)
  p <- strsplit(e, "|", fixed = TRUE)[[1]]
  if (length(p) >= i) p[i] else NA_character_
}

dv <- fread(DIV, header = FALSE, select = 1:4, col.names = c("chrom","start","end","strain"))

one <- function(r) {
  reg <- sprintf("%s:%d-%d", r$chrom, r$from, r$to)
  par <- PARENTS[[r$cross]]

  raw <- fread(cmd = sprintf("tabix %s %s", shQuote(GFF), reg), sep = "\t", header = FALSE,
               col.names = c("chrom","source","feature","start","end","score","strand","frame","attr"))
  mrna <- raw %>% filter(feature == "mRNA") %>%
    transmute(wbgene = sub("^gene:", "", gff_attr(attr, "Parent")),
              mrna.locus = gff_attr(attr, "locus"), transcript = gff_attr(attr, "Name"))
  genes <- raw %>% filter(feature == "gene") %>%
    transmute(wbgene = gff_attr(attr, "Name"), locus = gff_attr(attr, "locus"),
              sequence.name = gff_attr(attr, "sequence_name"),
              biotype = gff_attr(attr, "biotype"), strand, start, end) %>%
    left_join(mrna %>% distinct(wbgene, mrna.locus), by = "wbgene") %>%
    mutate(label = coalesce(locus, mrna.locus, sequence.name)) %>%
    select(-mrna.locus) %>% arrange(start)
  exons <- raw %>% filter(feature == "exon") %>%
    transmute(transcript = sub("^transcript:", "", gff_attr(attr, "Parent")), strand, start, end) %>%
    left_join(mrna %>% select(transcript, wbgene), by = "transcript") %>%
    inner_join(genes %>% select(wbgene, label), by = "wbgene") %>% arrange(start)

  q <- sprintf("bcftools query -r %s -s %s -f '%%POS\\t%%REF\\t%%ALT\\t[%%GT\\t]%%INFO/BCSQ\\n' %s",
               reg, paste(par, collapse = ","), shQuote(VCF))
  v <- fread(cmd = q, sep = "\t", header = FALSE, na.strings = "",
             col.names = c("pos","ref","alt","gt.p1","gt.p2","bcsq")) %>%
    mutate(across(starts_with("gt."), ~ gsub("\\|", "/", .x)),
           called = gt.p1 != "./." & gt.p2 != "./.")
  var <- v %>% filter(called, gt.p1 != gt.p2) %>% rowwise() %>%
    mutate(impact = sev(bcsq), consequence = fld(bcsq, 1), gene = fld(bcsq, 2),
           transcript = fld(bcsq, 3), aa.change = fld(bcsq, 6)) %>% ungroup() %>%
    mutate(consequence = coalesce(consequence, "unannotated"),
           alt.parent = ifelse(gt.p1 == "0/0", par[2], par[1]))

  d <- dv %>% filter(chrom == r$chrom, end > r$from, start < r$to, strain %in% par) %>%
    mutate(start = pmax(start, r$from), end = pmin(end, r$to))

  tag <- function(x) x %>% mutate(locus = r$locus, cross = r$cross, chrom = r$chrom,
                                  peak.bp = r$peak.bp, from = r$from, to = r$to)
  ## precomputed: tibble() evaluates its arguments in order, so a column named
  ## `genes` would shadow the genes table before the later expressions read it
  n_genes  <- nrow(genes); n_coding <- sum(genes$biotype == "protein_coding")
  n_sites  <- nrow(v);     n_nocall <- sum(!v$called); n_diff <- nrow(var)
  n_high   <- sum(var$impact == "HIGH"); n_mod <- sum(var$impact == "MODERATE")
  nc_p1 <- sum(v$gt.p1 == "./."); nc_p2 <- sum(v$gt.p2 == "./.")
  div_bp   <- if (nrow(d)) sum(d$end - d$start) else 0L
  div_str  <- paste(sort(unique(d$strain)), collapse = ",")
  list(genes = tag(genes), exons = tag(exons), var = tag(var), div = tag(d),
       summ = tibble(locus = r$locus, cross = r$cross, chrom = r$chrom,
                     peak.Mb = r$peak.Mb, peak.LOD = r$peak.LOD, dfreq = r$dfreq,
                     peak.bp = r$peak.bp, from = r$from, to = r$to,
                     genes = n_genes, coding = n_coding,
                     sites = n_sites, nocall = n_nocall, differ = n_diff,
                     HIGH = n_high, MODERATE = n_mod,
                     divergent.bp = div_bp, divergent.strains = div_str,
                     parent1 = par[1], parent2 = par[2],
                     nocall.p1 = nc_p1, nocall.p2 = nc_p2))
}

res <- lapply(seq_len(nrow(loci)), function(i) one(loci[i, ]))
write_tsv(bind_rows(map(res, "genes")), file.path(OUT, "mig6_locus_genes.tsv"))
write_tsv(bind_rows(map(res, "exons")), file.path(OUT, "mig6_locus_exons.tsv"))
write_tsv(bind_rows(map(res, "var")),   file.path(OUT, "mig6_locus_variants.tsv"))
write_tsv(bind_rows(map(res, "div")),   file.path(OUT, "mig6_locus_divergent.tsv"))
summ <- bind_rows(map(res, "summ"))
write_tsv(summ, file.path(OUT, "mig6_locus_summary.tsv"))
write_tsv(summ %>% select(locus, parent1, parent2, nocall.p1, nocall.p2, sites),
          file.path(OUT, "mig6_locus_nocall.tsv"))

cat("\n== per-locus census ==\n")
print(as.data.frame(summ %>% select(locus, peak.LOD, coding, sites, differ,
                                    HIGH, MODERATE, nocall.p1, nocall.p2,
                                    divergent.bp)),
      row.names = FALSE)
cat("\n== every HIGH- or MODERATE-impact difference ==\n")
pa <- bind_rows(map(res, "var")) %>% filter(impact %in% c("HIGH","MODERATE"))
print(as.data.frame(pa %>% select(locus, pos, impact, consequence, gene, aa.change, alt.parent)),
      row.names = FALSE)
msg("wrote five tables to ", OUT)
