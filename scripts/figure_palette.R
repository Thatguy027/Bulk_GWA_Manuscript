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
##   COL_JU1793 / COL_EIG / CROSS_COL[N2xXZ1516]  -- blue #0072B2
##   COL_N2     / CROSS_COL[JU1793xJU2466]        -- orange #E69F00
## The crosses deliberately take their colour from a parent.
##
## Greyscale is the one thing this palette does not survive: #0072B2 and
## #009E73 have nearly the same lightness. Figures that must read in black and
## white need a second channel -- linetype, shape or a direct label.

## --- strains and genotypes --------------------------------------------------
COL_N2      <- "#E69F00"   # reference strain, orange by convention
COL_JU1793  <- "#0072B2"
COL_JU2466  <- "#009E73"
COL_XZ      <- "#CC79A7"
COL_EDIT    <- "#999999"   # edited lines

## --- RNAi targets -----------------------------------------------------------
COL_MIG     <- "#56B4E9"   # mig-6
COL_POS     <- "#D55E00"   # pos-1

## --- crosses ----------------------------------------------------------------
CROSS_COL   <- c(N2xXZ1516 = "#0072B2", JU1793xJU2466 = "#E69F00")

## --- structure figures -----------------------------------------------------
## The focal residue in the SID-2 panels. It used to be JU1793's colour, since
## the model is the 96T allele, but JU1793 is now blue and would collide with
## the lysines K93/K132 drawn beside it. Vermillion is Okabe-Ito's emphasis
## slot; it is also COL_POS, and the two never share a panel.
COL_FOCAL   <- "#D55E00"   # T96, the residue the figure is about

## --- mapping annotation -----------------------------------------------------
COL_PEAK    <- "#111111"   # markers clearing a threshold
COL_EIG     <- "#0072B2"   # eigen-decomposition threshold (dashed)
COL_THR     <- "#666666"   # Bonferroni threshold (solid)
COL_REGION  <- "#D9D9D9"   # the interval the NIL series resolves
