## Locate the inputs baugh_frequencies(refresh = TRUE) needs -----------------
##
## Sourced ONLY on the slow path. Kept out of Figure1_common.R because that is a
## figure script and figure scripts must run from a deposit-only clone; this
## file is allowed to look in data/, as the other archive readers are.
##
## The deposit carries the CACHE, which is what the default path reads, so a
## fresh clone regenerates every figure without this file ever being sourced.
## What it finds is the genotype matrix and counts BEHIND that cache, which the
## deposit does not carry:
##
##   2024bootstrapINPUT.Rdata          data/baugh/ only
##   bootstrap array (ab, bootse)      deposited as baugh_bootstrap_array.rda,
##                                     byte-identical to data/baugh/
##                                     2024baugh_bootstrap_prediction.rda
##
## Returns a list of the two paths, erroring with both candidates named if
## neither exists.
## ---------------------------------------------------------------------------

baugh_refresh_paths <- function(deposit_dir) {
  pick <- function(label, ...) {
    cand <- c(...)
    hit <- cand[file.exists(cand)]
    if (!length(hit))
      stop("cannot find ", label, " -- looked for:\n  ",
           paste(cand, collapse = "\n  "),
           "\nThe deposit does not carry it; baugh_frequencies() reads the ",
           "cache instead.", call. = FALSE)
    hit[1]
  }
  list(
    boot  = pick("the genotype matrix and counts",
                 file.path(deposit_dir, "2024bootstrapINPUT.Rdata"),
                 file.path("data", "baugh", "2024bootstrapINPUT.Rdata")),
    bpred = pick("the bootstrap array",
                 file.path(deposit_dir, "baugh_bootstrap_array.rda"),
                 file.path("data", "baugh", "2024baugh_bootstrap_prediction.rda")))
}
