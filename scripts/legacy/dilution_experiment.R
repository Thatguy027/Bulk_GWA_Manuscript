library(BEDMatrix)
library(RcppML)
library(extraDistr)
library(tidyverse)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))

# load functions
source("common_functions.R")

# load genotypes - g
# and 
# genotype calls - genotype_calls_vcf
load("../data/genotypes/processed_genotype_matrix.Rda")

# load and clean up strain names
isotype_lookup <- data.table::fread("../data/meta/strain_isotype_lookup.tsv") 

strains <- data.table::fread("../data/experiments/20211111_BulkGWA/BC_strains.txt", col.names = c("strain","set"), header = F) %>%
  dplyr::left_join(., isotype_lookup, by = "strain") %>%
  na.omit() %>%
  dplyr::select(strain = isotype, set) %>%
  dplyr::filter(set %in% c("B","C"))

f_strains <- paste(strains$strain,strains$strain,sep = "_") %>% sort()

# load bulk sample ALT calls
alt_dir <- "../data/experiments/20211111_BulkGWA/aser/"

for(alt_files in grep("BC",grep("table",list.files(alt_dir), value = T), value = T)){
  
  sample = gsub(".table","",alt_files)
  print(sample)
  afdf <- data.table::fread(glue::glue("{alt_dir}{alt_files}"))%>%
    dplyr::select(chrom = contig, pos = position, ref = refAllele, alt = altAllele, alt_ct = altCount, dp = totalCount)%>%
    dplyr::mutate(af = alt_ct/dp,
                  sample = sample,
                  marker = paste0(chrom,":",pos,"_",alt)) 
  
  if(!exists("alt_df_bias")){
    alt_df_bias <- afdf
  } else {
    alt_df_bias <- dplyr::bind_rows(alt_df_bias, afdf)
  }
}

n_samples <- length(unique(alt_df_bias$sample))

# find low coverage variants across all sampels
dp_cut_df <- alt_df_bias %>%
  dplyr::select(marker, sample, dp, ref, alt) %>%
  tidyr::spread(sample, dp)

dp_per_site <- apply(dp_cut_df[,4:ncol(dp_cut_df)], 1, function(x){sum(x, na.rm = T)})

dp_cut_df$sum_dp = dp_per_site

ggplot(dp_cut_df) +
  aes(x = sum_dp/n_samples)+
  geom_histogram()+
  xlim(0,200)

good_dp_marker <- dp_cut_df %>% # 2919290
  dplyr::filter(sum_dp/n_samples > 5) %>% # 2898656
  dplyr::filter(sum_dp/n_samples < 100) %>% # 2897148
  # na.omit() %>% # 2896940
  dplyr::pull(marker)

pr_alt_df_bias <- alt_df_bias %>%
  dplyr::filter(marker %in% good_dp_marker)

uniq_mrkr <- unique(pr_alt_df_bias$marker)
g_pruned <- g[row.names(g) %in% uniq_mrkr,colnames(g) %in% f_strains]

save(g_pruned, pr_alt_df_bias, file="../data/experiments/20211111_BulkGWA/BC_gpruned_n_counts.rda")

#################################################### apply gentype filters to dataset 
g_variable_hiGT <- g_filters(genotypes = g_pruned,
                             vcf_stats = genotype_calls_vcf,
                             vcf_gt_rate = 1000)

allele_counts_variable_hiGT <- pr_alt_df_bias %>%
  dplyr::filter(marker %in% row.names(g_variable_hiGT))

flipped_bootstrap_input <- flip_common_variants(genotypes = g_variable_hiGT, 
                                                allele_counts = allele_counts_variable_hiGT)

####################################################  - remove sites with any NA

gt_with_na <- apply(flipped_bootstrap_input[[1]], 1, function(x){
  sum(is.na(x))
})

a_ct <- flipped_bootstrap_input[[2]][which(gt_with_na==0),]
gt <- flipped_bootstrap_input[[1]][which(gt_with_na==0),]

####################################################  - predict strain frequencies

fullGy=crossprod(gt,a_ct)
fullGGp=crossprod(gt)

predictions=apply(fullGy,2,function(x) {
  as.vector(RcppML::nnls(fullGGp,matrix(x), fast_nnls=T)) #fast_nnls=F, cd_maxit=10000, cd_tol=1e-10))
})

predictions=apply(predictions, 2, function(x) x/sum(x))

predictions_df <- data.frame(predictions) %>%
  dplyr::mutate(strain = colnames(fullGGp)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(strain = strsplit(strain,split = "_")[[1]][1]) %>%
  tidyr::gather(sample, frq, -strain) %>%
  dplyr::rename(new_frq = frq)

save(predictions_df, file="../data/experiments/20211111_BulkGWA/BC_flipCommon_NArm_predictions.rda")
load(file="../data/experiments/20211111_BulkGWA/BC_flipCommon_NArm_predictions.rda")

plot_df <- predictions_df %>%
  dplyr::left_join(., strains, by = "strain") %>%
  dplyr::group_by(set, sample) %>%
  dplyr::summarise(set_frq = round(sum(new_frq), digits = 2))

library(ggrepel)

ggplot(plot_df)+
  aes(x = sample, y = set_frq, fill = set, label = set_frq) +
  geom_point(size = 4, shape =21)+
  theme_bw(18)+
  scale_fill_manual(values = c("hotpink3", "cadetblue3"))+
  scale_color_manual(values = c("hotpink3", "cadetblue3"))+
  geom_text_repel(aes(color = set),
                  force_pull   = 0, # do not pull toward data points
                  nudge_x      = 0.2,
                  direction    = "y",
                  angle        = 0,
                  hjust        = 0,
                  segment.size = 0.2,
                  max.iter = 1e4, max.time = 1)+
  labs(x = "Sample", y = "Set frequency")+
  theme(legend.position = "none")

ggsave(filename = "../figures/BC_DNA_dilution_prediction_plot.pdf", height = 8, width = 10)




