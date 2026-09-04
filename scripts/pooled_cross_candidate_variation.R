## pooled_cross_candidate_variation.R --------------------------------------
## Does a gene under a QTL peak actually carry variation between the two
## parents of the cross that mapped it? Positional overlap alone is not
## evidence: a 900 kb support interval routinely contains dozens of genes with
## protein-altering differences, so a candidate is only interesting relative to
## that background.
##
## For each significant support interval this script counts, per gene, the
## number of distinct sites carrying a protein-altering consequence
## (missense / stop / frameshift / splice-acceptor-donor / start-lost /
## inframe indel) at which the two parents of that cross differ, and reports
## where a named gene of interest ranks among them.
##
##   Rscript scripts/pooled_cross_candidate_variation.R
##
## Requires bcftools on PATH and the BCSQ-annotated CeNDR isotype VCF.
##
## Outputs
##   plots/pooled_cross_intersection/TABLE_candidate_coding_variation.tsv
##   plots/pooled_cross_intersection/TABLE_interval_candidate_burden.tsv
##   plots/pooled_cross_intersection/FIG12_candidate_variation.{pdf,png}
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(ggtext)
})

b     <- readRDS("data/pooled_cross_intersection/bundle.rds")
VCF   <- "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213/bcsq.vcf.gz"
OUT   <- "plots/pooled_cross_intersection"
CACHE <- "data/pooled_cross_intersection/candidate_variation.rds"
stopifnot(file.exists(VCF), nzchar(Sys.which("bcftools")))

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x

## which two parents belong to each cross, and their column in the query
PARENTS <- c(N2xXZ1516 = "N2|XZ1516", JU1793xJU2466 = "JU1793|JU2466")

## consequences we count as protein-altering
HI <- paste("missense", "stop_gained", "stop_lost", "frameshift",
            "splice_acceptor", "splice_donor", "start_lost", "inframe",
            sep = "|")

## ---------------------------------------------------------------------------
## one interval -> per-gene protein-altering site counts for that cross's parents
## ---------------------------------------------------------------------------
interval_genes <- function(cross, chrom, start, end) {
  pp  <- strsplit(PARENTS[[cross]], "|", fixed = TRUE)[[1]]
  cmd <- sprintf(
    "bcftools query -f '%%POS\\t%%REF\\t%%ALT\\t%%INFO/BCSQ[\\t%%GT]\\n' -s %s -r %s:%d-%d %s 2>/dev/null",
    paste(pp, collapse = ","), chrom, start, end, shQuote(VCF))
  raw <- system(cmd, intern = TRUE)
  if (!length(raw)) return(NULL)

  f <- tstrsplit(raw, "\t", fixed = TRUE)
  d <- data.table(pos = as.integer(f[[1]]), ref = f[[2]], alt = f[[3]],
                  bcsq = f[[4]], a = f[[5]], b = f[[6]])
  ## keep only sites where the two parents are confidently different
  d <- d[a %in% c("0/0", "1/1") & b %in% c("0/0", "1/1") & a != b]
  n_disc <- nrow(d)
  if (!nrow(d)) return(list(genes = NULL, variants = NULL, n_disc = 0L))

  ## which parent carries the alternate allele at each discordant site
  d[, `:=`(parent1 = pp[1], parent2 = pp[2],
           alt_carrier = fifelse(a == "1/1", pp[1], pp[2]),
           chrom = chrom)]

  ## BCSQ packs one or more consequence records per site, comma separated,
  ## each pipe-delimited: consequence|gene|transcript|biotype|strand|aa|dna
  ann <- d[bcsq != "."]
  recs <- if (nrow(ann))
    ann[, .(rec = unlist(strsplit(bcsq, ",", fixed = TRUE))), by = pos] else
    data.table(pos = integer(), rec = character())
  if (nrow(recs)) {
    parts <- tstrsplit(recs$rec, "|", fixed = TRUE)
    recs[, `:=`(csq = parts[[1]], gene = parts[[2]],
                aa  = if (length(parts) >= 6) parts[[6]] else NA_character_)]
    recs <- recs[!is.na(gene) & nzchar(gene)]
  }

  ## variant-level table: every discordant site, with its worst consequence
  worst <- if (nrow(recs))
    unique(recs[, .(pos, csq, gene, aa)])[
      , .(gene = paste(sort(unique(gene)), collapse = ";"),
          consequence = paste(sort(unique(csq)), collapse = ";"),
          aa_change = paste(sort(unique(aa[nzchar(aa) & !is.na(aa)])), collapse = ";"),
          protein_altering = any(grepl(HI, csq))), by = pos] else
    data.table(pos = integer(), gene = character(), consequence = character(),
               aa_change = character(), protein_altering = logical())

  variants <- merge(d[, .(chrom, pos, ref, alt, gt_parent1 = a, gt_parent2 = b,
                          parent1, parent2, alt_carrier)],
                    worst, by = "pos", all.x = TRUE)
  variants[is.na(protein_altering), protein_altering := FALSE]
  variants[is.na(consequence), consequence := "intergenic/noncoding"]
  setorder(variants, pos)

  hi <- recs[grepl(HI, csq)]
  if (!nrow(hi)) return(list(genes = NULL, variants = variants, n_disc = n_disc))

  ## a site can appear once per transcript; count distinct sites per gene
  g <- unique(hi[, .(gene, pos, csq, aa)])[
    , .(n_sites = uniqueN(pos),
        changes = paste(sort(unique(ifelse(nzchar(aa) & !is.na(aa), aa, csq))),
                        collapse = ", ")),
    by = gene][order(-n_sites)]
  list(genes = g, variants = variants, n_disc = n_disc)
}

## ---------------------------------------------------------------------------
## run over every significant interval of the six focal contrasts
## ---------------------------------------------------------------------------
KEY6 <- c("N2xXZ1516 | ht115 vs mig6", "N2xXZ1516 | ht115 vs pos1",
          "N2xXZ1516 | pos1 vs mig6",
          "JU1793xJU2466 | ht115 vs mig6", "JU1793xJU2466 | ht115 vs pos1",
          "JU1793xJU2466 | mig6 vs pos1")

ivs <- b$ivs %>% filter(key %in% KEY6, significant) %>% arrange(cross, label, chrom)

if (file.exists(CACHE)) {
  res <- readRDS(CACHE)
  msg("loaded cache ", CACHE)
} else {
  msg(nrow(ivs), " intervals to profile")
  res <- pmap(ivs %>% select(cross, label, chrom, lcon, rcon, peak.LOD, marker),
              function(cross, label, chrom, lcon, rcon, peak.LOD, marker) {
    r <- interval_genes(cross, chrom, lcon, rcon)
    msg("  ", cross, " ", label, " ", marker,
        "  discordant=", r$n_disc %||% 0,
        " genes=", if (is.null(r$genes)) 0 else nrow(r$genes))
    list(meta = tibble(cross, label, chrom, lcon, rcon, peak.LOD, marker,
                       n_discordant = r$n_disc %||% 0L,
                       n_genes = if (is.null(r$genes)) 0L else nrow(r$genes)),
         genes = r$genes, variants = r$variants)
  })
  saveRDS(res, CACHE)
  msg("cached ", CACHE)
}

burden <- map_dfr(res, "meta") %>%
  mutate(width.kb = (rcon - lcon) / 1e3,
         genes_per_100kb = n_genes / (width.kb / 100))

per_gene <- map_dfr(seq_along(res), function(i) {
  g <- res[[i]]$genes; if (is.null(g)) return(NULL)
  m <- res[[i]]$meta
  as_tibble(g) %>%
    mutate(cross = m$cross, label = m$label, chrom = m$chrom, marker = m$marker,
           peak.LOD = m$peak.LOD, n_genes_in_interval = m$n_genes,
           rank = rank(-n_sites, ties.method = "min"), .before = 1)
})

## ---------------------------------------------------------------------------
## how the RNAi-pathway genes we annotate actually fare
## ---------------------------------------------------------------------------
rnai <- b$rnai_genes$gene

focus <- per_gene %>% filter(gene %in% rnai) %>%
  select(cross, label, marker, peak.LOD, gene, n_sites, rank,
         n_genes_in_interval, changes)

## RNAi genes that lie in an interval but carry NO protein-altering difference
in_interval <- pmap_dfr(ivs %>% select(cross, label, chrom, lcon, rcon, marker, peak.LOD),
  function(cross, label, chrom, lcon, rcon, marker, peak.LOD) {
    b$rnai_genes %>% filter(chrom == !!chrom, end >= lcon, start <= rcon) %>%
      transmute(cross, label, marker, peak.LOD, gene, class)
  })

silent <- in_interval %>%
  anti_join(focus %>% select(cross, label, marker, gene),
            by = c("cross", "label", "marker", "gene")) %>%
  mutate(n_sites = 0L, rank = NA_integer_, changes = NA_character_) %>%
  left_join(burden %>% select(cross, label, marker, n_genes_in_interval = n_genes),
            by = c("cross", "label", "marker"))

focus_all <- bind_rows(focus %>% left_join(in_interval %>%
                         select(cross, label, marker, gene, class),
                       by = c("cross", "label", "marker", "gene")),
                       silent) %>%
  arrange(desc(peak.LOD), desc(n_sites))

write_tsv(focus_all, file.path(OUT, "TABLE_candidate_coding_variation.tsv"))
write_tsv(burden %>% arrange(desc(peak.LOD)),
          file.path(OUT, "TABLE_interval_candidate_burden.tsv"))
write_tsv(per_gene %>% arrange(desc(peak.LOD), rank),
          file.path(OUT, "TABLE_all_genes_coding_variation.tsv"))

## ---------------------------------------------------------------------------
## variant-level dumps: every parental difference under every cross QTL
## ---------------------------------------------------------------------------
VDIR <- file.path(OUT, "variants_by_qtl")
dir.create(VDIR, showWarnings = FALSE, recursive = TRUE)

all_variants <- map_dfr(res, function(x) {
  if (is.null(x$variants) || !nrow(x$variants)) return(NULL)
  m <- x$meta
  as_tibble(x$variants) %>%
    mutate(cross = m$cross, contrast = m$label, qtl = m$marker,
           peak.LOD = m$peak.LOD, peak.position = NA_integer_, .before = 1)
})
## attach the peak position for each interval
all_variants <- all_variants %>%
  select(-peak.position) %>%
  left_join(ivs %>% select(cross, contrast = label, qtl = marker,
                           peak.position, lcon, rcon),
            by = c("cross", "contrast", "qtl")) %>%
  mutate(dist_to_peak_kb = round((pos - peak.position) / 1e3, 2)) %>%
  relocate(dist_to_peak_kb, .after = pos)

write_tsv(all_variants, file.path(OUT, "TABLE_parental_variants_all_qtl.tsv.gz"))
msg("wrote TABLE_parental_variants_all_qtl.tsv.gz  (",
    nrow(all_variants), " rows across ", n_distinct(all_variants$qtl), " intervals)")

## one file per QTL, and a protein-altering-only companion
all_variants %>% group_split(cross, contrast, qtl) %>% walk(function(v) {
  slug <- sprintf("%s__%s__%s",
                  v$cross[1],
                  gsub("[^A-Za-z0-9]+", "-", v$contrast[1]),
                  gsub("[^A-Za-z0-9]+", "_", v$qtl[1]))
  write_tsv(v, file.path(VDIR, paste0(slug, ".tsv.gz")))
})
write_tsv(all_variants %>% filter(protein_altering),
          file.path(OUT, "TABLE_parental_variants_protein_altering.tsv"))
msg("wrote ", length(list.files(VDIR)), " per-QTL variant files to ", VDIR)

## ---------------------------------------------------------------------------
## the sid-2 variants specifically -- the locus the NIL / fine-mapping work
## follows up, where T96K is the variant that turns out to contribute
## ---------------------------------------------------------------------------
SID2 <- b$rnai_genes %>% filter(gene == "sid-2")
sid2_variants <- map_dfr(names(PARENTS), function(cr) {
  pp <- strsplit(PARENTS[[cr]], "|", fixed = TRUE)[[1]]
  r  <- interval_genes(cr, SID2$chrom, SID2$start - 2000, SID2$end + 2000)
  if (is.null(r$variants) || !nrow(r$variants)) return(NULL)
  as_tibble(r$variants) %>% mutate(cross = cr, .before = 1)
}) %>%
  mutate(in_sid2 = pos >= SID2$start & pos <= SID2$end,
         focal   = grepl("96T>96K", aa_change, fixed = TRUE)) %>%
  arrange(cross, pos)

write_tsv(sid2_variants, file.path(OUT, "TABLE_sid2_parental_variants.tsv"))

cat("\n== sid-2 protein-altering variants, by cross ==\n")
print(as.data.frame(sid2_variants %>%
        filter(protein_altering, grepl("sid-2", gene)) %>%
        transmute(cross, pos, ref, alt, aa_change, consequence,
                  alt_carrier, focal)), row.names = FALSE)

## ---------------------------------------------------------------------------
## FIG 12
## ---------------------------------------------------------------------------
PAL <- c(nxz = "#0B7A75", juju = "#D57A00", hi = "#C4302B", grey = "#8C8C8C")
CROSS_LAB <- c(N2xXZ1516 = "N2 × XZ1516", JU1793xJU2466 = "JU1793 × JU2466")
wrap_txt <- function(x, n = 120) paste(strwrap(x, width = n), collapse = "\n")

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = base_size),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_text(face = "bold", size = base_size + 0.5),
          plot.subtitle = element_text(size = base_size - 1, colour = "grey30"),
          plot.title.position = "plot",
          legend.key.size = grid::unit(9, "pt"))
}

## (a) candidate burden: how many genes could explain each peak
pa <- burden %>%
  mutate(lab = paste0(chrom, ":", round(lcon / 1e6, 1), "-", round(rcon / 1e6, 1)),
         row = paste0(CROSS_LAB[cross], " · ", label),
         cross = factor(cross, levels = names(CROSS_LAB))) %>%
  ggplot(aes(width.kb, n_genes, colour = cross)) +
  geom_point(aes(size = peak.LOD), alpha = 0.8) +
  ## labelling all 36 intervals is illegible; name the strong ones only
  ggrepel::geom_text_repel(aes(label = ifelse(peak.LOD >= 150 | n_genes >= 100, lab, "")),
                           size = 2.4, seed = 5, min.segment.length = 0,
                           box.padding = 0.45, segment.size = 0.25,
                           max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = unname(PAL[c("nxz", "juju")]),
                      labels = unname(CROSS_LAB), name = NULL) +
  scale_size_continuous(range = c(1.2, 5), name = "peak LOD") +
  labs(x = "support-interval width (kb)",
       y = "genes with a protein-altering\nparental difference",
       title = "Every interval holds many candidates",
       subtitle = wrap_txt(paste0(
         "Positional overlap with a QTL peak is not evidence on its own: interval width, ",
         "not peak strength, sets how many genes could explain the signal."), 62)) +
  theme_pub(9)

## (b) the RNAi-pathway genes we annotate, scored against their own interval
pb <- focus_all %>%
  filter(!is.na(n_genes_in_interval)) %>%
  mutate(tag = sprintf("%s\n%s · %s", gene, CROSS_LAB[cross], marker),
         tag = fct_reorder(tag, n_sites),
         status = case_when(n_sites == 0 ~ "no coding difference",
                            rank <= 5    ~ "top-5 in its interval",
                            TRUE         ~ "present but unremarkable")) %>%
  ggplot(aes(n_sites, tag, fill = status)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = ifelse(is.na(rank), "0 of  —",
                               sprintf("%d, rank %d/%d", n_sites, rank, n_genes_in_interval))),
            hjust = -0.12, size = 2.5, colour = "grey20") +
  scale_fill_manual(values = c(`top-5 in its interval` = PAL[["hi"]],
                               `present but unremarkable` = PAL[["grey"]],
                               `no coding difference` = "grey85"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.34))) +
  labs(x = "distinct protein-altering sites between the two parents", y = NULL,
       title = "Only *sid-2* stands out against its interval",
       subtitle = wrap_txt(paste0(
         "Each bar is an RNAi-pathway gene lying inside a significant support interval, ",
         "counted for the two parents of the cross that mapped it. sid-5 carries no coding ",
         "difference between JU1793 and JU2466 at all, and rde-12 is an ordinary member of ",
         "a hyper-divergent interval. mut-15 does not appear because it falls in no support ",
         "interval. sid-2's two JU-cross missense variants are V5L and T96K, the latter the ",
         "variant the NIL fine-mapping follows up."), 96)) +
  theme_pub(9) +
  theme(axis.text.y = element_text(size = 6.8, lineheight = 0.95),
        plot.title = element_markdown(face = "bold", size = 9.5),
        legend.position = "top", legend.margin = margin(b = -6))

fig12 <- (pa | pb) + plot_layout(widths = c(1, 1.35)) +
  plot_annotation(
    title = "Which genes under a peak actually differ between the parents that mapped it",
    subtitle = wrap_txt(paste0(
      "Protein-altering = missense, stop gained/lost, frameshift, splice acceptor/donor, ",
      "start lost, or inframe indel, from BCFtools/csq annotation of the CeNDR isotype VCF. ",
      "Counts are distinct genomic sites, collapsed across transcripts."), 145),
    tag_levels = "a",
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8.5, colour = "grey25"),
                  plot.tag = element_text(face = "bold", size = 12)))

ggsave(file.path(OUT, "FIG12_candidate_variation.pdf"), fig12,
       width = 13, height = 6.4, device = cairo_pdf)
ggsave(file.path(OUT, "FIG12_candidate_variation.png"), fig12,
       width = 13, height = 6.4, dpi = 300, bg = "white")
msg("wrote FIG12_candidate_variation.{pdf,png}")

## ---------------------------------------------------------------------------
cat("\n== interval candidate burden ==\n")
print(as.data.frame(burden %>% arrange(desc(peak.LOD)) %>%
        transmute(cross, label, marker, peak.LOD = round(peak.LOD, 1),
                  width.kb = round(width.kb), n_discordant, n_genes)),
      row.names = FALSE)

cat("\n== RNAi-pathway genes inside an interval, scored ==\n")
print(as.data.frame(focus_all %>%
        transmute(cross, label, marker, gene, n_sites, rank,
                  of = n_genes_in_interval, peak.LOD = round(peak.LOD, 1))),
      row.names = FALSE)
