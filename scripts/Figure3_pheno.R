## Figure 3, version 3 -- genome-wide scan with the pooled phenotype inset ---
##
##   Rscript scripts/Figure3_pheno.R  ->  plots/Figure3_pheno.{pdf,png}
##
##   A  the JU1793 x JU2466 HT115-vs-pos-1 cross scan, all six chromosomes,
##      with the pooled pos-1 phenotype distribution inset over chromosomes
##      I and II, where the scan is flat
##   B  the NIL introgressions on the right arm of chromosome III
##   C  embryo hatching under pos-1 RNAi for the same strains
##
## Everything shared with the other two variants is in Figure3_common.R.
##
## What the inset adds
## -------------------
## It closes the loop back to Figure 2. Both cross parents were phenotyped in
## the pooled panel, and they sit at opposite ends of it -- JU1793 near the
## resistant extreme, JU2466 in the sensitive tail -- which is why they were
## picked as parents. Stated in prose that is an assertion; drawn here it is
## the reason the cross in panel A has any power at all, and it puts the
## pooled experiment, the cross and the NILs in one frame.
##
## The trait is the VST pos-1 response, the same trait mapped by GWAS in
## Figure 2, so the inset is not a separate measurement introduced here.
## Positive = the strain gained frequency under pos-1 RNAi = resistant.
##
## The inset occupies dead space rather than displacing anything: chromosomes
## I and II peak at LOD 14 and 7 against a chrIII peak of 140, so the top ~70%
## of those two panels is empty. It is placed with patchwork's inset_element
## aligned to the panel region, so it sits below the facet strips and above
## the chrI / chrII traces, both of which stay visible.
## ---------------------------------------------------------------------------

source("scripts/Figure3_common.R")

msg("panel A: cross scan, genome-wide, with pooled phenotype inset")
scan_t <- thin_scan(load_scan())
pA  <- panel_A_genome(scan_t)
pIn <- pheno_inset()

## left/right in panel-relative units. Chromosome I ends at 15.3% of the
## genome and II at 30.9%, so 0.02-0.32 covers those two panels.
pA_in <- pA + inset_element(pIn, left = 0.020, bottom = 0.28,
                            right = 0.325, top = 1.00, align_to = "panel")

msg("panel B: NIL introgressions")
pB <- panel_B()

msg("panel C: hatching")
pC <- panel_C()

fig <- three_panel(pA_in, pB, pC, heights = c(1, 1))

ggsave(file.path(OUT, "Figure3_pheno.pdf"), fig, width = 9.6, height = 6.6,
       device = cairo_pdf)
ggsave(file.path(OUT, "Figure3_pheno.png"), fig, width = 9.6, height = 6.6,
       dpi = 300, bg = "white")
msg("wrote Figure3_pheno.{pdf,png}")
