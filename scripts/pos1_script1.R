library(BEDMatrix)
library(RcppML)
library(extraDistr)
library(tidyverse)
library(smplot2)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))

# load functions
source("common_functions.R")

# load genotypes - g
# and 
# genotype calls - genotype_calls_vcf
load("../data/genotypes/processed_genotype_matrix.Rda")

experiment_folder <- "../data/experiments/20230406_pos1-rnai/"

#################################################### load and clean up strain names
isotype_lookup <- data.table::fread("../data/meta/strain_isotype_lookup.tsv") 

strains <- data.table::fread("../data/meta/strain_sets/BulkCe_strainsets.tsv", col.names = c("strain","set"), header = F) %>%
  dplyr::left_join(., isotype_lookup, by = "strain") %>%
  na.omit() %>%
  dplyr::select(strain = isotype, set) %>%
  dplyr::filter(set %in% c("C","B","E","F","G"))

f_strains <- paste(strains$strain,strains$strain,sep = "_") %>% sort()

####################################################  load bulk sample ALT calls
alt_dir <- glue::glue("{experiment_folder}aser/")

for(alt_files in grep("table",list.files(alt_dir), value = T)){
  
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


# remove samples with crap coverage
average_dp_sample <- alt_df_bias %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(med_dp = median(dp,na.rm = T)) 

ggplot(average_dp_sample)+
  geom_point(aes(x=sample, y = med_dp))+
  theme_bw(18)+
  labs(x = "Sample", y = "Median depth")

good_samples <- average_dp_sample %>%
  dplyr::filter(med_dp > 10) %>%
  dplyr::pull(sample)

# find low coverage variants across all samples
dp_cut_df <- alt_df_bias %>%
  # dplyr::filter(sample %in% good_samples) %>%
  dplyr::select(marker, sample, dp, ref, alt) %>%
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
g_pruned <- g[row.names(g) %in% uniq_mrkr,colnames(g) %in% f_strains]
# g_pruned is the baseline for all other genotype processing
# pr_alt_df_bias is the baseline for all other allele count processing

#################################################### - apply filters to the genotypes

#################################################### apply genotype filters to dataset 

# rerun because NAs were switched to 0s above
g_variable_hiGT <- g_filters(genotypes = g_pruned,
                             vcf_stats = genotype_calls_vcf,
                             vcf_gt_rate = 1000)

allele_counts_variable_hiGT <- pr_alt_df_bias %>%
  dplyr::filter(marker %in% row.names(g_variable_hiGT))

flipped_bootstrap_input <- flip_common_variants(genotypes = g_variable_hiGT, 
                                                allele_counts = allele_counts_variable_hiGT)

# this needs to go onto the cluster for running the bootstrap procedure
save(flipped_bootstrap_input, file = glue::glue("{experiment_folder}2024bootstrapINPUT.Rdata"))
save(allele_counts_variable_hiGT, file = glue::glue("{experiment_folder}2024downsample_alleleCts.Rdata"))

load(glue::glue("{experiment_folder}2024bootstrapINPUT.Rdata"))
load(glue::glue("{experiment_folder}2024downsample_alleleCts.Rdata"))

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

# load("../data/experiments/baugh_wgs/processed_data/20220908_flipCommon_processed_bootstraps.RData")

save(predictions_df, file = glue::glue("{experiment_folder}2024pub_prediction_df.Rdata"))

load(glue::glue("{experiment_folder}2024pub_prediction_df.Rdata"))

####################################################  - load and append meta data
metaDF <- data.table::fread(glue::glue("{experiment_folder}/meta_files/pos1_meta.tsv"))

predictions_df <- dplyr::left_join(predictions_df, metaDF, by = "sample")

####################################################  - load and append bootstrap data from cluster

load(glue::glue("{experiment_folder}2024baugh_bootstrap_prediction.rda"))

strains <- unlist(lapply(strsplit(colnames(flipped_bootstrap_input[[1]]), split = "_"), function(x){x[1]}))

bootme=apply(ab, c(1,2), function(x) mean(x))

#################################################### - convert bootstrapping output to a dataframe
boot_df <- data.frame(strain = strains, bootmean = bootme, bootse = bootse) %>%
  tidyr::gather(stat_sample, value, -strain) %>%
  tidyr::separate(stat_sample, into = c("stat","sample"), convert = T, sep = "\\.") %>%
  tidyr::spread(stat, value)

boot_df %>%
  dplyr::left_join(., predictions_df, by = c("sample", "strain")) %>% 
  dplyr::filter(grepl("baseline", condition)) %>% 
  ggplot()+
  aes(x=new_frq, y = bootmean)+
  geom_point(size = 0.01)+
  geom_pointrange(aes(xmin=new_frq-bootse, xmax=new_frq+bootse))

prediction_boot_df <- boot_df %>%
  dplyr::left_join(., predictions_df, by = c("sample", "strain")) 

save(prediction_boot_df, file = glue::glue("{experiment_folder}2024pub_prediction_bootstrap_combined_df.Rdata"))

#################################################### - construct phenotypes and visualize

t2_ctrl <- prediction_boot_df %>% 
  dplyr::filter(time == "t2", condition %in% c("ctrlA","ctrlB"))%>%
  dplyr::mutate(rep = ifelse(grepl("B",condition), "B", "A")) %>%
  dplyr::select(strain, btm_ctrl = bootmean, bts_ctrl = bootse, ctrl_frq=new_frq, control_cond=condition, time, rep)

t2_pos1 <- prediction_boot_df %>% 
  dplyr::filter(time == "t2", condition %in% c("pos1B_pos1B","pos1A_pos1A")) %>%
  dplyr::mutate(rep = ifelse(grepl("B",condition), "B", "A")) %>%
  dplyr::left_join(., t2_ctrl, by = c("strain", "time", "rep"))

t2_pos1 %>%
  ggplot()+
  aes(x = btm_ctrl, y = bootmean)+
  geom_point()+
  facet_grid(.~rep)

t2_pos1 %>%
  dplyr::mutate(d_ctrl = bootmean - btm_ctrl,
                f_frq = new_frq - ctrl_frq) %>% 
  ggplot()+
  aes(x = d_ctrl, y = f_frq)+
  geom_point()+
  facet_grid(.~rep)

plot_df <- t2_pos1 %>%
  dplyr::select(strain, condition, btm_ctrl, bts_ctrl, bootmean, bootse) %>%
  tidyr::gather(variable, value, -strain, -condition) %>%
  # tidyr::spread(condition, value) %>% 
  dplyr::filter(!variable%in%c("bootse","bts_ctrl"))

ggplot(plot_df)+
  aes(x = variable, y = value, group = strain)+
  geom_line() +
  facet_wrap(~condition)

# visualize replicate correlation  
t2_pos1 %>%
  dplyr::mutate(delta_boot = bootmean-btm_ctrl) %>%
  dplyr::select(strain, sample, delta_boot) %>%
  tidyr::spread(sample, delta_boot) %>%
  dplyr::mutate(cors = cor(A10, A13, method = "spearman")) %>% 
  ggplot()+
  aes(x = A10, y = A13)+
  geom_point()+
  sm_statCorr(color = "firebrick3", corr_method = "spearman",
              linetype = "dashed", size = 1)+
  geom_point(color = "#0072B2", size = 4, data = t2_pos1 %>%
               dplyr::mutate(delta_boot = bootmean-btm_ctrl) %>%
               dplyr::select(strain, sample, delta_boot) %>%
               tidyr::spread(sample, delta_boot) %>%
               dplyr::filter(strain == "ECA738"))+
  geom_point(color = "#CC7AA7", size = 4, data = t2_pos1 %>%
               dplyr::mutate(delta_boot = bootmean-btm_ctrl) %>%
               dplyr::select(strain, sample, delta_boot) %>%
               tidyr::spread(sample, delta_boot) %>%
               dplyr::filter(strain == "ECA760"))+
  theme_bw(18)+
  labs(x = "Replicate 1", y= "Replicate 2")

ggsave(filename = "../figures/original_pos1_rep_correlation.pdf", height = 6, width=8)


t2pos1_traits <- t2_pos1 %>%
  dplyr::mutate(delta_boot = bootmean-btm_ctrl) %>%
  dplyr::select(strain, sample, delta_boot) %>%
  tidyr::spread(sample, delta_boot) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(mean_rep = (A10+A13)/2) %>% 
  dplyr::filter(!(A10 == 0 & A13 == 0))

write.table(t2pos1_traits, file = "../data/experiments/20230406_pos1-rnai/2024_t2pos1_traits.tsv",
            quote = F, col.names = T, row.names = F)

