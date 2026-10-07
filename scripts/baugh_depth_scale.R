## The stored downsampling key is HALF the sequencing depth it names ---------
##
## Sourced, not run. Everything that plots or reports the Baugh downsampling
## series goes through here, so the correction is defined once.
##
## THE BUG, IN THE SCRIPT THAT PRODUCED THE DATA. In
## scripts/baugh_L1_DownSample_Counts.R the count matrix handed to
## downsampleCounts() is built by gathering alt_ct and ref_ct into ONE column:
##
##   ds_input <- ... mutate(ref_ct = dp - alt_ct) %>%
##     select(sample, marker, alt_ct, ref_ct) %>%
##     gather(ref_alt, ct, -sample, -marker) %>%
##     unite(uid, marker, ref_alt, sep = "-") %>% spread(sample, ct)
##
## so ds_input has TWO rows per marker. The draw is then
##
##   downsampleCounts(count.matrix = ds_input, downsample.total = ds_dp * nrow(ds_input))
##
## and each draw is one read's allele observation at a marker. The reads drawn
## per sample are therefore ds_dp * 2 * n_markers, and the mean reads per marker
## -- which is what "depth" means -- is 2 * ds_dp, not ds_dp. Had the intent been
## literal the denominator would have been the marker count, not nrow(ds_input).
##
## Verified rather than assumed. Replaying that reshape on an allele-count table
## with the same columns gives exactly two rows per marker, suffixed alt_ct and
## ref_ct. The Baugh samples carry a mean depth of 24x to 92x (median 41x,
## estimated from the archived alternate counts and the per-marker alternate
## frequency), so the doubled draw is feasible without replacement at every step
## -- the thinnest sample has 24.1x against the 20 reads per marker the top step
## needs, which is tight and passes. A literal reading would also have been
## feasible, so that test does not discriminate; the arithmetic above does.
##
## WHAT IS AND IS NOT CORRECTED. The ds_n stored in
## supplemental_data/deconvolution/baugh_downsampled_slopes.rda is left alone:
## it is the script's own parameter and the key every downstream join uses.
## Only the DISPLAYED and REPORTED depth is doubled. No correlation moves --
## the estimates are what they always were, and only the depth each is
## attributed to changes.
## ---------------------------------------------------------------------------

## ds_dp as stored -> mean reads per marker
DS_DEPTH_FACTOR <- 2
ds_depth <- function(ds_n) ds_n * DS_DEPTH_FACTOR

## axis labeller: takes the stored key, prints the true depth
ds_depth_label <- function(ds_n, suffix = "×") paste0(ds_depth(ds_n), suffix)
