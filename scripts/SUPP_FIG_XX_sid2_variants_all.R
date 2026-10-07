## Supplement -- every protein-altering sid-2 variant in the latest CeNDR release
##
##   Rscript scripts/SUPP_FIG_XX_sid2_variants_all.R
##     -> plots/SUPP_FIG_XX_sid2_variants_all.{pdf,png}
##
## Figure 4D draws only the sites where the three mapping parents differ. This is
## the whole catalogue it was cut from: every protein-altering sid-2 change in
## CeNDR 20250625 -- missense, the in-frame deletion and the frameshift -- drawn
## by the same sid2_variant_panel() so the two read as one design.
##
## WHAT DIFFERS FROM 4D, AND WHY
##   Source. 4D reads sid2_variants_cendr.tsv, annotated only at sites
##   segregating among the four cross parents. This reads the 20250625 export,
##   sid2_variants_cendr20250625.tsv, which annotates the population.
##   Bars are COUNTS of isotypes carrying the alternate, not frequencies: the
##   export lists carriers only, so there is no denominator to divide by.
##   Residue 151 is split into its two forms. 13680412 alone gives 151T; with the
##   partner SNV at 13680413 it gives 151I. 4D keeps the single "A151I/T" row
##   because only 151T is parental; here each form has its own count.
##
## The parental columns are read from the same carrier lists. The band still
## marks sites where JU1793 and JU2466 differ, as in 4D.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggtext)
})

OUT <- "plots"
VAR <- "supplemental_data/structure/sid2_variants_cendr20250625.tsv"

source("scripts/figure_palette.R")
source("scripts/figure_theme.R")
source("scripts/sid2_variant_panel.R")
TOPO_COL2 <- c(`Signal peptide` = "grey72", `Extracellular` = "#9EC5DE",
               `TM helix` = "#37474F", `Cytoplasmic` = "grey88")
panel_title <- function(letter) NULL      # one panel; the caption carries the title
## the charge ramp and limits of Figure 4C; QLIM must equal the renderer's
QLIM <- 2
ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(64)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

v <- read_tsv(VAR, comment = "#", show_col_types = FALSE) %>%
  mutate(carr = strsplit(carriers, " "))
C <- setNames(v$carr, v$site)

## residue 151, resolved into its two forms
v151 <- v %>% filter(site == "snp412")
r151 <- bind_rows(
  v151 %>% mutate(label = "A151T", alt_aa = "T",
                  carr = list(setdiff(C$snp412, C$snp413))),
  v151 %>% mutate(label = "A151I", alt_aa = "I",
                  carr = list(intersect(C$snp412, C$snp413))))
stopifnot(all(C$snp413 %in% C$snp412))   # the partner never occurs alone

vr <- bind_rows(v %>% filter(!site %in% c("snp412", "snp413")), r151) %>%
  mutate(n = lengths(carr))
parent_aa <- function(s) ifelse(map_lgl(vr$carr, ~ s %in% .x), vr$alt_aa, vr$ref_aa)
vr <- vr %>%
  mutate(ju1793_aa = parent_aa("JU1793"), ju2466_aa = parent_aa("JU2466"),
         xz1516_aa = parent_aa("XZ1516"),
         parents_differ = ju1793_aa != ju2466_aa,
         focal = label == "T96K") %>%
  arrange(residue, label)

cat("\n== protein-altering sid-2 sites, CeNDR 20250625 ==\n")
print(as.data.frame(vr %>% select(label, consequence, n, ju1793_aa, ju2466_aa,
                                  xz1516_aa)), row.names = FALSE)
on96 <- vr %>% filter(label != "T96K") %>%
  mutate(on_96K = map_int(carr, ~ sum(.x %in% C$T96K)))
cat(sprintf("\nisotypes named: %d | T96K %d, all of them 153T: %s | 153T without 96K: %d\n",
            length(unique(unlist(vr$carr))), length(C$T96K),
            all(C$T96K %in% C$P153T), length(setdiff(C$P153T, C$T96K))))
cat("sites never on 96K:", paste(on96$label[on96$on_96K == 0], collapse = ", "),
    "| only on 96K:", paste(on96$label[on96$on_96K == on96$n], collapse = ", "), "\n")

## ---- pins: the caption quotes these ----------------------------------------
PINNED <- c(V5L = 9, A10_I12del = 4, Q37P = 7, S43G = 27, D78A = 6,
            `N94Lfs*6` = 1, T96K = 317, S111G = 1, G117R = 1, M141V = 6,
            Q144P = 6, A151T = 85, A151I = 62, V152A = 52, P153T = 328,
            T158A = 1, L209M = 147, T223A = 9, A291P = 2, G299E = 13)
got <- setNames(vr$n, vr$label)[names(PINNED)]
if (any(is.na(got)) || any(got != PINNED))
  stop("pinned carrier counts moved: ",
       paste(names(PINNED)[is.na(got) | got != PINNED], collapse = ", "))
stopifnot(nrow(vr) == 20, all(C$T96K %in% C$P153T),
          length(setdiff(C$P153T, C$T96K)) == 11,
          identical(vr$label[vr$parents_differ], c("V5L", "T96K")))
msg("pins agree")

p <- sid2_variant_panel(vr, letter = NULL, base_size = BASE_SIZE,
                        ramp = ramp, qlim = QLIM, bar = "n",
                        bar_max = max(vr$n), bar_text = function(x) as.character(x),
                        bar_header = "Isotypes carrying")

ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_variants_all.pdf"), p, width = 6.2,
       height = 9.2, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_variants_all.png"), p, width = 6.2,
       height = 9.2, dpi = 300, bg = "white")
msg("wrote SUPP_FIG_XX_sid2_variants_all.{pdf,png}")
