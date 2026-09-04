library(BEDMatrix)
library(RcppML)
library(extraDistr)
library(tidyverse)
library(ggpubr)
library(GGally)

# g is a matrix (n x s , n markers by s strains)
# t1 is a vector alternate allele counts  (n x 1 , n markers)
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load in genotypes
mats=list()
for(chrom in c('I','II','III', 'IV', 'V', 'X')){
  base::print(chrom)
  g=as.matrix(BEDMatrix(paste0('../data/CeNDR20210121_Plink/', chrom, '.bed') ))
  g=t(g)
  # g[is.na(g)]=0
  #to simplify for now
  g[g==2]=1
  mats[[chrom]]=g
}

g=rbind(mats[[1]], mats[[2]], mats[[3]], mats[[4]], mats[[5]], mats[[6]])

rm(mats)

# clean up strain names
isotype_lookup <- data.table::fread("meta_files/strain_isotype_lookup.tsv") 

strains <- data.table::fread("meta_files/BulkCe_strainsets.tsv", col.names = c("strain"), header = F) %>%
  dplyr::left_join(., isotype_lookup, by = "strain") %>%
  na.omit() %>%
  dplyr::select(strain = isotype) 

f_strains <- paste(strains$strain,strains$strain,sep = "_") %>% sort()


g_exp <- g[,f_strains]
g_exp[1:10,1:10]

allele_cts = rowSums(g_exp,na.rm = T)
allele_cts[1:10]

variable_sites <- allele_cts[names(which(allele_cts > 0))]

g_exp <- g_exp[names(variable_sites),]
dim(g_exp)

rm(g)

marker_df <- data.frame(marker = rownames(g_exp)) %>%
  tidyr::separate(marker, into = c("chrom", "pos"), sep = ":", convert = T, remove = F) %>%
  tidyr::separate(pos, into = c("pos", "alt"), sep = "_", convert = T)

# load bulk sample ALT calls
alt_dir <- "data/aser/"

sample_df <- data.frame(sample_file = grep("table",list.files(alt_dir), value = T)) %>%
  dplyr::mutate(sample_id = gsub(pattern = ".table", replacement = "", sample_file)) %>%
  tidyr::separate(sample_id, into = c("sample_name", "lane", "junk")) %>%
  dplyr::select(-sample_name, -junk) %>%
  tidyr::pivot_wider(names_from = lane, values_from = sample_file) %>%
  tidyr::unnest(c(A,B))


for(alt_files in 1:nrow(sample_df)){
  
  lane_a <- sample_df[alt_files,1]$A
  lane_b <- sample_df[alt_files,2]$B
  
  sample_name <- gsub("_A_1.table","",lane_a)
  
  base::print(sample_name)
  
  lane_a_afdf <- data.table::fread(glue::glue("{alt_dir}{lane_a}"))%>%
    dplyr::select(chrom = contig, pos = position, ref = refAllele, alt = altAllele, alt_ct_a = altCount, dp_a = totalCount) %>%
    dplyr::select(-ref)
  
  lane_b_afdf <- data.table::fread(glue::glue("{alt_dir}{lane_b}"))%>%
    dplyr::select(chrom = contig, pos = position, ref = refAllele, alt = altAllele, alt_ct_b = altCount, dp_b = totalCount) %>%
    dplyr::select(-ref)
  
  lane_combined <- dplyr::left_join(marker_df, lane_a_afdf, by = c("chrom", "pos","alt")) %>%
    dplyr::left_join(., lane_b_afdf, by = c("chrom", "pos","alt")) %>%
    dplyr::mutate(dp = dp_a + dp_b,
                  alt_ct = alt_ct_a + alt_ct_b,
                  sample = sample_name) %>%
    dplyr::select(marker:alt,sample,alt_ct, dp)
  
  if(!exists("alt_df_bias")){
    alt_df_bias <- lane_combined
  } else {
    alt_df_bias <- dplyr::bind_rows(alt_df_bias, lane_combined)
  }
}

save(alt_df_bias, file="data/20240320_unprocessed_allele_counts.Rdata")
load(file="data/20240320_unprocessed_allele_counts.Rdata")


# remove samples with crap coverage
average_dp_sample <- alt_df_bias %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(med_dp = median(dp,na.rm = T)) 

ggplot(average_dp_sample)+
  geom_point(aes(x=sample, y = med_dp))+
  theme_bw(18)+
  labs(x = "Sample", y = "Median depth")

good_samples <- average_dp_sample %>%
  dplyr::filter(med_dp >= 10) %>%
  dplyr::pull(sample)

# find low coverage variants across all samples
dp_cut_df <- alt_df_bias %>%
  dplyr::filter(sample %in% good_samples) %>%
  dplyr::select(marker, sample, dp, alt) %>%
  tidyr::spread(sample, dp)

dp_per_site <- apply(dp_cut_df[,4:ncol(dp_cut_df)], 1, function(x){sum(x, na.rm = T)})

dp_cut_df$sum_dp = dp_per_site

ggplot(dp_cut_df) +
  aes(x = sum_dp/(ncol(dp_cut_df)-4))+
  geom_histogram(bins = 100)+
  xlim(0,100)

good_dp_marker <- dp_cut_df %>% # 2919290
  # dplyr::filter(sum_dp/(ncol(dp_cut_df)-4) > 10) %>% # 2918858
  # dplyr::filter(sum_dp/(ncol(dp_cut_df)-4) < 250) %>% # 2914674
  na.omit() %>% # 2914517
  dplyr::pull(marker)

pr_alt_df_bias <- alt_df_bias %>%
  dplyr::filter(marker %in% good_dp_marker)

uniq_mrkr <- unique(pr_alt_df_bias$marker)
g_pruned <- g_exp[row.names(g_exp) %in% uniq_mrkr,colnames(g_exp) %in% f_strains]
# g_pruned is the baseline for all other genotype processing
# pr_alt_df_bias is the baseline for all other allele count processing

#################################################### - apply filters to the genotypes

genotype_calls_vcf <- data.table::fread("../20220909_Ascr_BZ/published_data/WI.20210121.hard-filter.isotype_with_cM_allele_counts.tsv",
                                        col.names = c("chrom","pos","ref","alt","ac","an"))

# remove sites that are not variable in the mapping population
# remove sites that are not genotyped in a user-defined number of CeNDR isolates
# genotypes is a sites (row) x strain (column) matrix
# vcf_stats is a data frame with columns c("chrom","pos","ref","alt","ac","an"), generated by bcftools query
# # # bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n' WI.20210121.hard-filter.isotype_with_cM.vcf.gz > WI.20210121.hard-filter.isotype_with_cM_allele_counts.tsv
# vcf_gt_rate is a number corresponding to the minimum number of allowed genotype calls in the CeNDR VCF to include the site as good for NNLS

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

g_variable_hiGT <- g_filters(genotypes = g_pruned,
                             vcf_stats = genotype_calls_vcf,
                             vcf_gt_rate = 1000)

allele_counts_variable_hiGT <- pr_alt_df_bias %>%
  dplyr::filter(marker %in% row.names(g_variable_hiGT))

#################################################### - flip coding of common alleles in genotype matrix and the alt count data

flip_common_variants <- function(genotypes, allele_counts){
  
  # find high frequency variants
  af_matrix <- data.frame(marker = row.names(genotypes), af = rowSums(genotypes, na.rm = T)/ncol(genotypes))
  
  high_alt <- af_matrix %>%
    dplyr::filter(af>0.5)
  
  # flip coding of high frequency variants
  g_changed <- genotypes[row.names(genotypes) %in% high_alt$marker,]
  
  g_changed[g_changed==1] <- -1
  g_changed[g_changed==0] <- 1
  g_changed[g_changed==-1] <- 0
  
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
    
    g_n_counts <- list(g_resort,t1m)
    
    return(g_n_counts)
  } else{
    print("counts and g dont match")
    return(NULL)
  }
  
  
  
  
}

flipped_bootstrap_input <- flip_common_variants(genotypes = g_variable_hiGT, 
                                                allele_counts = allele_counts_variable_hiGT)

save(flipped_bootstrap_input, file="data/cluster_data/20230410_cleanGT_flippedCommon_n_Counts.Rdata")

