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
##   Nothing in the data: both read the CeNDR 20250625 export
##   (sid2_cendr20250625_variants.csv) through sid2_cendr_variants(), and both
##   draw frequencies over the release's 684 isotypes. 4D keeps only the rows
##   where the three mapping parents differ; this keeps all 20.
##   Residue 151 appears as both of its forms here, 151T and 151I; 4D shows
##   only 151T, the one a parent (XZ1516) carries.
##
## The parental columns are read from the same carrier lists. The band still
## marks sites where JU1793 and JU2466 differ, as in 4D.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggtext)
})

OUT <- "plots"

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

vr <- sid2_parent_columns(sid2_cendr_variants())
C  <- setNames(vr$carriers, vr$label)

cat("\n== protein-altering sid-2 sites, CeNDR 20250625 ==\n")
print(as.data.frame(vr %>% transmute(label, consequence, n, af = round(af, 4),
                                     ju1793_aa, ju2466_aa, xz1516_aa)),
      row.names = FALSE)
on96 <- vr %>% filter(label != "T96K") %>%
  mutate(on_96K = map_int(carriers, ~ sum(.x %in% C$T96K)))
cat(sprintf("\nisotypes named: %d | T96K %d, all of them 153T: %s | 153T without 96K: %d\n",
            length(unique(unlist(vr$carriers))), length(C$T96K),
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
                        ramp = ramp, qlim = QLIM, bar_text = sid2_pct)

ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_variants_all.pdf"), p, width = 6.2,
       height = 9.2, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_sid2_variants_all.png"), p, width = 6.2,
       height = 9.2, dpi = 300, bg = "white")
msg("wrote SUPP_FIG_XX_sid2_variants_all.{pdf,png}")
