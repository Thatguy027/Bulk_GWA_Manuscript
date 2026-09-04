g_filters <- function(genotypes, 
                      vcf_stats, 
                      vcf_gt_rate = 1000 ){
  
  af_matrix <- data.frame(marker = row.names(genotypes), af = rowSums(genotypes, na.rm = T)/ncol(genotypes))
  
  variant_sites_inpop <- af_matrix %>%
    dplyr::filter(af>0)
  
  genotypes_variantSites <- genotypes[row.names(genotypes) %in% variant_sites_inpop$marker,]
  
  # only consider markers with high genotyping rate in the CeNDR collection
  hi_gt_rate_markers <- vcf_stats %>%
    dplyr::filter(an>vcf_gt_rate, # most strains have genotype information
                  nchar(ref) == nchar(alt),
                  nchar(ref) == 1) %>% # not indels 
    tidyr::unite(marker, chrom, pos, sep = ":", remove = F) %>%
    tidyr::unite(marker, marker, alt, sep = "_") %>%
    dplyr::pull(marker)
  
  genotypes_variantSites_high_gt <- genotypes_variantSites[row.names(genotypes_variantSites) %in% hi_gt_rate_markers, ]
  
  return(genotypes_variantSites_high_gt)
}

flip_common_variants <- function(genotypes, allele_counts){
  
  # find high frequency variants
  af_matrix <- data.frame(marker = row.names(genotypes), af = rowSums(genotypes, na.rm =T)/ncol(genotypes))
  
  high_alt <- af_matrix %>%
    dplyr::filter(af>0.5)
  
  # flip coding of high frequency variants
  g_changed <- genotypes[row.names(genotypes) %in% high_alt$marker,]
  
  g_changed <- abs(g_changed-1)
  
  # extract unchanged genotypes
  g_not_changed <- genotypes[!(row.names(genotypes) %in% high_alt$marker),]
  
  # append unchanged and flipped genotypes
  g_modified <- rbind(g_not_changed,g_changed)
  
  # Flip high frequency counts, spread alt calls, and sort by marker
  flip_high_alt_counts <- pr_alt_df_bias %>%
    dplyr::mutate(alt_count = ifelse(marker %in% high_alt$marker, dp-alt_ct, alt_ct)) %>%
    dplyr::select(marker, sample, alt_count) %>%
    dplyr::filter(marker %in% row.names(g_modified)) %>%
    tidyr::spread(sample, alt_count) %>%
    tidyr::separate(marker, into = c("chrom","pos","alt"), convert = T, remove = F) %>%
    dplyr::arrange(chrom, pos)
  
  # probably redundant, but filter modified matrix to match alt counts
  g_final <- g_modified[row.names(g_modified) %in% (flip_high_alt_counts$marker),]
  
  # resort genotypes and put back into a matrix
  g_resort <- g_final %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="marker") %>%
    tidyr::separate(marker, into = c("chrom","pos","alt"), convert = T, remove = F) %>%
    dplyr::arrange(chrom, pos) %>%
    dplyr::select(-(chrom:alt)) %>%
    tibble::column_to_rownames(var = "marker") %>%
    as.matrix()
  
  if(identical(flip_high_alt_counts$marker, row.names(g_resort))){
    
    t1m <- dplyr::select(flip_high_alt_counts, -marker, -chrom, -pos, -alt) %>%
      as.matrix()
    
    g_resort[is.na(g_resort)] <- 0
    
    g_n_counts <- list(g_resort,t1m)
    
    return(g_n_counts)
  } else{
    print("counts and g dont match")
    return(NULL)
  }
  
}

plot.inflation <- function (x, size = 2) {
  
  # Get the number of p-values.
  n <- length(x)
  
  # Compute the negative log10(p-values), and sort them from largest
  # to smallest.
  y <- rev(sort(-log10(x)))
  
  # Create the q-q plot.
  return(ggplot(data.frame(x = -log10((1:n)/n),y = y),aes(x = x,y = y)) +
           geom_abline(intercept = 0,slope = 1,color = "magenta") +
           geom_point(color = "dodgerblue",shape = 20,size = 2) +
           labs(x = "Expected -log10 p-value",
                y = "Observed -log10 p-value") +
           theme(axis.line = element_blank()))
}

calc_inflation_factor <- function(map_df){
  z_scores = map_df$beta/map_df$se
  chisq <- z_scores^2
  print(median(chisq)/qchisq(0.5,1))
}
