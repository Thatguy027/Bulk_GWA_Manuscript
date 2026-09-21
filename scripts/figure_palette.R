## The manuscript's figure palette, in one place ------------------------------
##
##   source("scripts/figure_palette.R")
##
## Every figure script takes its colours from here rather than defining its own,
## so a palette change is one edit instead of eighty-three scattered ones.
##
## OKABE-ITO. The eight-colour set designed for colour-vision deficiency
## (Okabe & Ito 2008, https://jfly.uni-koeln.de/color/), which is also the
## conventional choice in genetics. It replaces a hand-assembled palette whose
## significant-peak marker was brick red (#C4302B) and whose eigen threshold was
## green (#1A7F5A) -- the two drawn on the same Manhattan panels, and the one
## pairing a red-green deficiency collapses outright.
##
## N2 IS ORANGE BY CONVENTION and nothing else in the palette is, which is why
## JU1793 moved off the orange it previously had.
##
## EIGHT HUES OVER THIRTEEN ROLES, so some are reused. Each reuse is between
## roles that never share a panel, and that is the property to preserve when
## editing this file:
##   COL_XZ / COL_EIG / CROSS_COL[N2xXZ1516]  -- blue #0072B2
##   COL_N2 / CROSS_COL[JU1793xJU2466]         -- orange #E69F00
## The N2 x XZ1516 cross takes XZ1516's blue, which is a parent's colour.
## JU1793 is reddish purple rather than blue because blue against JU2466's
## green was the weakest contrast in the palette.
##
## Greyscale is the one thing this palette does not survive for the categorical
## roles: #0072B2 and #009E73 have nearly the same lightness. Figures that must
## read in black and white need a second channel -- linetype, shape or a direct
## label. The Manhattan panels are the exception and are greyscale by design.

## --- strains and genotypes --------------------------------------------------
COL_N2      <- "#E69F00"   # reference strain, orange by convention
COL_JU1793  <- "#CC79A7"
COL_JU2466  <- "#009E73"
COL_XZ      <- "#0072B2"   # the blue JU1793 vacated; no figure shows both
COL_EDIT    <- "#999999"   # edited lines

## --- RNAi targets -----------------------------------------------------------
COL_MIG     <- "#56B4E9"   # mig-6
COL_POS     <- "#D55E00"   # pos-1

## --- crosses ----------------------------------------------------------------
CROSS_COL   <- c(N2xXZ1516 = "#0072B2", JU1793xJU2466 = "#E69F00")

## --- Manhattan points -------------------------------------------------------
## The scans are drawn in alternating shades of grey by chromosome, the
## conventional Manhattan treatment, so that colour in those panels means
## "this marker cleared Bonferroni" and nothing else. Both shades carry a
## slight cool bias rather than being neutral grey.
COL_GW_A    <- "#2F3438"   # odd-numbered chromosomes
COL_GW_B    <- "#9AA1A6"   # even-numbered chromosomes

## Alternating shade for a scan, as a per-row vector. It is returned as a vector
## rather than mapped through aes() on purpose: the panels that need it already
## spend their colour scale on something else, and ggplot allows one scale per
## aesthetic. `d` needs a `chrom` factor in genome order.
gw_shade <- function(d) ifelse(as.integer(d$chrom) %% 2 == 1, COL_GW_A, COL_GW_B)

## --- data marks -------------------------------------------------------------
## COL_PT is the mass of points in a Manhattan or a scatter; COL_PEAK is drawn
## on top of it, so the two are separated by LIGHTNESS rather than hue -- point
## size reinforces it (0.4 against 1.2). The old palette separated them by hue,
## slate against brick red, which is the pairing this palette exists to remove.
COL_PT      <- "#6E7B85"   # data points; a grey carrying a little of the blue
COL_FIT     <- "#D55E00"   # regression or fit line over those points
COL_HIST    <- "#009E73"   # histogram fill

## The mig-6-against-pos-1 contrast. It was #E08214, an orange, which now reads
## as N2; reddish purple is free in every panel it appears in.
COL_MIGR    <- "#CC79A7"

## --- allele aliases ---------------------------------------------------------
## Residue 96 states ARE their parent strains, so they are not separate colours.
COL_96T     <- COL_JU1793   # the JU1793 allele
COL_96K     <- COL_JU2466   # the JU2466 allele

## --- structure figures -----------------------------------------------------
## The focal residue in the SID-2 panels. It used to be JU1793's colour, since
## the model is the 96T allele, but JU1793 is now blue and would collide with
## the lysines K93/K132 drawn beside it. Vermillion is Okabe-Ito's emphasis
## slot; it is also COL_POS, and the two never share a panel.
COL_FOCAL   <- "#D55E00"   # T96, the residue the figure is about

## --- mapping annotation -----------------------------------------------------
COL_PEAK    <- "#CD2626"   # firebrick3: markers clearing Bonferroni
COL_EIG     <- "#0072B2"   # eigen-decomposition threshold (dashed)
COL_THR     <- "#666666"   # Bonferroni threshold (solid)
COL_REGION  <- "#D9D9D9"   # the interval the NIL series resolves
