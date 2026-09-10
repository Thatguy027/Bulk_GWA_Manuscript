## Do the censused-interval genes have eQTL, and do the parents differ? -------
##
##   Rscript scripts/candidate_eqtl_expression.R
##     -> plots/diagnostics/TABLE_candidate_eqtl.tsv
##
## Two independent questions about every gene in the five censused windows (the
## 37 kb NIL interval and the four 100 kb mig-6 QTL windows):
##
##   1. does it have a mapped eQTL in the 207-isolate expression study, and is
##      that eQTL local (near the gene, i.e. plausibly cis) or distant?
##   2. do the two parents of the cross that mapped the locus differ in its
##      expression, and by how much relative to the spread across all 207?
##
## These are complementary to the coding-variant census. A gene with no
## protein-altering difference between the parents can still be the cause
## through expression -- set-25 and lam-2 are exactly that shape -- and a gene
## with missense differences and no expression difference is a coding candidate
## only.
##
## SOURCE IS OUTSIDE THIS REPOSITORY: the 207-isolate expression matrix, its
## eQTL table and its feature table. Set CE207_DIR to override.
##
## WHAT THE PARENTAL COMPARISON CAN AND CANNOT SAY. One expression value per
## strain, so there is no within-strain replication and no p-value on a
## parental difference. The difference is therefore reported as a z-score
## against the 207-strain distribution -- how unusual the gap is, not whether
## it is significant. Expression is on the study's own transformed scale.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(data.table); library(dplyr); library(readr); library(tidyr); library(stringr)})

E    <- Sys.getenv("CE207_DIR", "/Users/Stefan/UCLA/Genomics_Data/CeNDR/expression")
MAP  <- "supplemental_data/mapping"
OUT  <- "plots/diagnostics"
EXPR <- file.path(E, "Ce207expression.csv")
QTL  <- file.path(E, "ce207_qtl.tsv")
for (f in c(EXPR, QTL)) if (!file.exists(f)) stop("missing source: ", f, call. = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## the genes to ask about, and which cross each window belongs to
nil <- read_tsv(file.path(MAP, "nil_interval_genes.tsv"), show_col_types = FALSE) %>%
  transmute(locus = "NIL interval III:13.66-13.70", cross = "JU1793xJU2466",
            wbgene, label, biotype)
mg <- read_tsv(file.path(MAP, "mig6_locus_genes.tsv"), show_col_types = FALSE) %>%
  transmute(locus, cross, wbgene, label, biotype)
genes <- bind_rows(nil, mg) %>% filter(biotype == "protein_coding") %>% distinct()
PAR <- list(N2xXZ1516 = c("N2","XZ1516"), JU1793xJU2466 = c("JU1793","JU2466"))
msg("coding genes across the five windows: ", nrow(genes))

## which of them carry a protein-altering parental difference
pa <- bind_rows(
  read_tsv(file.path(MAP, "nil_interval_parent_variants.tsv"), show_col_types = FALSE) %>%
    transmute(gene, impact),
  read_tsv(file.path(MAP, "mig6_locus_variants.tsv"), show_col_types = FALSE) %>%
    transmute(gene, impact)) %>%
  filter(impact %in% c("HIGH","MODERATE")) %>% count(gene, name = "coding.changes")

## ---- eQTL ------------------------------------------------------------------
q <- fread(QTL, sep = "\t") %>% as_tibble()
eq <- q %>% filter(WormBaseGeneID %in% genes$wbgene) %>%
  group_by(wbgene = WormBaseGeneID) %>%
  summarise(n.eQTL = n(),
            classes = paste(sort(unique(eQTL_classification)), collapse = "/"),
            best.logP = max(logP), max.var.exp = max(var_exp),
            hyperdivergent = paste(sort(unique(HyperDivergent)), collapse = "/"),
            .groups = "drop")

## ---- parental expression ---------------------------------------------------
x <- fread(EXPR, sep = ",") %>% as_tibble()
strain_cols <- setdiff(names(x), c("transcript","WormBaseGeneID","GeneName","biotype","h2","H2"))
xg <- x %>% filter(WormBaseGeneID %in% genes$wbgene)
msg("expression rows matched: ", nrow(xg), " (transcript level)")

long <- xg %>% select(all_of(c("WormBaseGeneID","transcript","h2","H2", strain_cols))) %>%
  pivot_longer(all_of(strain_cols), names_to = "strain", values_to = "expr")
stats <- long %>% group_by(WormBaseGeneID, transcript) %>%
  summarise(h2 = first(h2), H2 = first(H2),
            mean207 = mean(expr, na.rm = TRUE), sd207 = sd(expr, na.rm = TRUE),
            .groups = "drop")
pick <- long %>% inner_join(genes %>% select(wbgene, locus, cross),
                            by = c("WormBaseGeneID" = "wbgene"),
                            relationship = "many-to-many") %>%
  rowwise() %>% filter(strain %in% PAR[[cross]]) %>% ungroup() %>%
  mutate(which = ifelse(strain == sapply(cross, function(c) PAR[[c]][1]), "p1", "p2")) %>%
  select(WormBaseGeneID, transcript, locus, cross, which, expr) %>%
  pivot_wider(names_from = which, values_from = expr) %>%
  inner_join(stats, by = c("WormBaseGeneID","transcript")) %>%
  mutate(delta = p1 - p2, z = delta / sd207)

res <- pick %>%
  inner_join(genes %>% select(wbgene, label), by = c("WormBaseGeneID" = "wbgene")) %>%
  left_join(eq, by = c("WormBaseGeneID" = "wbgene")) %>%
  left_join(pa, by = c("label" = "gene")) %>%
  mutate(coding.changes = tidyr::replace_na(coding.changes, 0L),
         n.eQTL = tidyr::replace_na(n.eQTL, 0L)) %>%
  select(locus, gene = label, transcript, cross, h2, H2,
         p1, p2, delta, sd207, z, n.eQTL, classes, best.logP, max.var.exp,
         hyperdivergent, coding.changes) %>%
  arrange(locus, desc(abs(z)))
write_tsv(res, file.path(OUT, "TABLE_candidate_eqtl.tsv"))

cat(sprintf("\ngenes with an expression measurement: %d of %d coding genes\n",
            dplyr::n_distinct(res$gene), dplyr::n_distinct(genes$label)))
cat(sprintf("genes with at least one mapped eQTL   : %d\n",
            dplyr::n_distinct(res$gene[res$n.eQTL > 0])))
cat("\n== the biggest parental expression differences, |z| >= 1.5 ==\n")
print(as.data.frame(res %>% filter(abs(z) >= 1.5) %>%
  transmute(locus = str_replace(locus, " Mb$",""), gene,
            p1 = round(p1,2), p2 = round(p2,2), z = round(z,2),
            eQTL = ifelse(n.eQTL > 0, classes, "-"),
            var.exp = ifelse(is.na(max.var.exp), NA, round(max.var.exp,2)),
            coding = coding.changes)), row.names = FALSE)

cat("\n== every gene with a mapped eQTL ==\n")
print(as.data.frame(res %>% filter(n.eQTL > 0) %>%
  transmute(locus = str_replace(locus, " Mb$",""), gene, n.eQTL, classes,
            logP = round(best.logP,1), var.exp = round(max.var.exp,2),
            z = round(z,2), coding = coding.changes, hyperdiv = hyperdivergent)),
  row.names = FALSE)
msg("wrote TABLE_candidate_eqtl.tsv")
