library(BEDMatrix)
library(RcppML)
library(extraDistr)
library(tidyverse)
library(GGally)
library(broom)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# function to downsample reads ------------------------------------------------------------------------------------------
downsampleCounts=function(count.matrix, downsample.total, with.replacement=FALSE) {
  downsample.count.matrix=matrix(0, nrow=nrow(count.matrix), ncol=ncol(count.matrix))
  rownames(downsample.count.matrix)=rownames(count.matrix)
  colnames(downsample.count.matrix)=colnames(count.matrix)
  total=ncol(count.matrix)
  # create progress bar
  pb =txtProgressBar(min = 1, max = ncol(count.matrix), style = 3)
  for(i in 1:ncol(count.matrix)){
    setTxtProgressBar(pb, i)
    
    if(with.replacement){
      dsgenes=sample(rep(rownames(count.matrix), as.integer(count.matrix[,i])), downsample.total,replace=T)
      
    }else{
      dsgenes=sample(rep(rownames(count.matrix), as.integer(count.matrix[,i])), downsample.total)
    }
    dcounts=plyr::count(dsgenes)
    downsample.count.matrix[match(as.character(dcounts$x), rownames(count.matrix)),i]=dcounts$freq
  }
  close(pb)
  return(downsample.count.matrix)
}
#--------------------------------------------------------------------------------------------------------------------------
load(file = "data/baugh/2024bootstrapINPUT.Rdata") # matrix formatted genotypes and allele counts
# save(file = "../data/experiments/baugh_wgs/2024downsample_alleleCts.Rdata") # formatted allele counts
#--------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------- prediction on full data set

####################################################  - remove sites with any NA

gt_with_na <- apply(flipped_bootstrap_input[[1]], 1, function(x){
  sum(is.na(x))
})

a_ct <- flipped_bootstrap_input[[2]][which(gt_with_na==0),]
gt <- flipped_bootstrap_input[[1]][which(gt_with_na==0),]

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
  tidyr::gather(sample, frq, -strain)

#-------------------------------------------------------------------------------------------------------------------------- generate downsample input

# allele_counts_variable_hiGT <- as.data.frame(flipped_bootstrap_input[[2]]) %>%
#   dplyr::mutate(marker= row.names(flipped_bootstrap_input[[1]]))

# dplyr::select(marker, sample, alt_count, dp) 

average_dp_sample <- allele_counts_variable_hiGT %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(med_dp = median(dp))

ds_input <- allele_counts_variable_hiGT %>%
  dplyr::filter(marker %in% row.names(gt)) %>%
  dplyr::mutate(ref_ct = dp - alt_ct) %>%
  dplyr::select(sample, marker, alt_ct, ref_ct) %>%
  tidyr::gather(ref_alt, ct, -sample, -marker) %>%
  tidyr::unite(uid,  marker, ref_alt, sep ="-") %>%
  tidyr::spread(sample, ct) %>%
  tibble::column_to_rownames(var = "uid") %>%
  as.matrix()

dim(ds_input)
dim(gt)
dim(a_ct)

#--------------------------------------------------------------------------------------------------------------------------  downsample counts

for(ds_dp in c(0.25, 0.5, 1, 3, 5, 10)){
  print(ds_dp)
  test_ct <- data.frame(downsampleCounts(count.matrix = ds_input, downsample.total = ds_dp*nrow(ds_input))) %>%
    dplyr::mutate(downsample_x = ds_dp) %>%
    tibble::rownames_to_column(var = "marker") %>%
    tidyr::separate(marker, into = c("marker","alt_ref"), sep = "-") %>%
    dplyr::filter(alt_ref == "alt_ct") %>%
    dplyr::select(-alt_ref)
  
  if(!exists("full_ds_df")){
    full_ds_df <- test_ct
  } else {
    full_ds_df <- dplyr::bind_rows(full_ds_df, test_ct)
  }
  print(nrow(test_ct))
}

#-------------------------------------------------------------------------------------------------------------------------- process downsample counts
# prepare and run strain frequency predictions for downsampled data sets
df2predict <- full_ds_df %>%
  tidyr::gather(sample, ct, -marker, -downsample_x) %>%
  tidyr::unite(sample_ds, sample, downsample_x, sep = "_") %>%
  tidyr::spread(sample_ds, ct) %>%
  tidyr::separate(marker, into = c("chrom", "pos"), remove = F, convert = T) %>%
  dplyr::arrange(chrom, pos) 

identical(df2predict$marker, row.names(gt))

#-------------------------------------------------------------------------------------------------------------------------- downsample counts strain inference
ds_predict_in <- df2predict %>%
  dplyr::select(-(marker:pos)) %>%
  as.matrix()

fullGy=crossprod(gt,ds_predict_in)

ds_predictions=apply(fullGy,2,function(x) {
  as.vector(RcppML::nnls(fullGGp,matrix(x), fast_nnls=T)) #fast_nnls=F, cd_maxit=10000, cd_tol=1e-10))
  
})

ds_predictions=apply(ds_predictions, 2, function(x) x/sum(x))

ds_predictions_df <- data.frame(ds_predictions) %>%
  dplyr::mutate(strain = colnames(fullGGp)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(strain = strsplit(strain,split = "_")[[1]][1]) %>%
  tidyr::gather(sample, ds_frq, -strain) %>%
  # dplyr::rowwise() %>%
  dplyr::mutate(ds_n = str_split(sample, "_")) %>%
  tidyr::unnest(ds_n) %>%
  dplyr::filter(ds_n %in% c("0.25","0.5","1","3","5","10")) %>%
  dplyr::mutate(ds_rep = paste0("_",ds_n)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(n_sample = gsub(ds_rep, "", sample)) %>%
  dplyr::select(strain, sample = n_sample, ds_frq, ds_n) %>%
  dplyr::mutate(ds_n = as.numeric(ds_n)) %>%
  dplyr::left_join(., predictions_df, by = c("sample","strain"))

#-------------------------------------------------------------------------------------------------------------------------- plot downsampled frq estimates per sample

dir.create("../data/experiments/baugh_wgs/2024Downsample_plots")

for(sample_name in unique(ds_predictions_df$sample)){
  
  sample_dp <- average_dp_sample %>%
    dplyr::filter(sample == sample_name) %>%
    dplyr::pull(med_dp)
  
  plot_to_save <- ds_predictions_df %>%
    dplyr::filter(sample == sample_name, strain != "N2") %>%
    ggplot()+
    aes(x = frq , y = ds_frq)+
    geom_point(size = 0.5)+
    theme_bw(12)+
    scale_color_manual(values = c("gray50", "firebrick3"))+
    # geom_errorbar(aes(xmin=frq-bootse, xmax=frq+bootse), height=.02)+
    ggpubr::stat_cor(method = "spearman", label.x = 0.01, label.y = 0.06,inherit.aes = T,size = 2, color = "firebrick3")+
    facet_grid(sample~ds_n) +
    labs(x = "Full coverage frq", y = "Downsampled frq", title = glue::glue("{sample_name}; {sample_dp}X Depth"))
  
  ggsave(plot_to_save, filename = glue::glue("../data/experiments/baugh_wgs/2024Downsample_plots/{sample_name}_downsample_cor_flipped.pdf"),height = 6, width = 10)
}



save(ds_predictions_df, file = "data/baugh/2026downsampled_dataset.rda")

load("../data/experiments/baugh_wgs/downsampled_dataset.rda")
load("../data/experiments/baugh_wgs/2024_processedBOOTs_with_MIP.RData")

library(broom)

baseline_df <- ds_predictions_df %>% 
  dplyr::filter(strain!="N2") %>%
  tidyr::separate(sample, into = c("replicate", "day"), sep = "_", remove = F) %>%
  dplyr::mutate(baseline = ifelse(grepl("baseline", sample), T, F)) %>%
  dplyr::mutate(day = as.numeric(gsub("d", "", day))) %>%
  dplyr::filter(!baseline, day == 1) %>%
  dplyr::select(replicate, strain, ds_n, base_frq = ds_frq)

delta_df <- ds_predictions_df %>% 
  dplyr::filter(strain!="N2") %>%
  tidyr::separate(sample, into = c("replicate", "day"), sep = "_", remove = F) %>%
  dplyr::mutate(baseline = ifelse(grepl("baseline", sample), T, F)) %>%
  dplyr::mutate(day = as.numeric(gsub("d", "", day))) %>%
  dplyr::filter(!baseline, day != 17) %>%
  dplyr::select(strain:ds_n) %>%
  dplyr::left_join(., baseline_df, by = c("replicate", "strain", "ds_n")) %>%
  dplyr::mutate(delta_d1_wgs = ds_frq - base_frq) %>%
  dplyr::select(strain:day, ds_n, delta_d1_wgs)


wgs_slope <- delta_df %>% 
  group_by(replicate, strain,ds_n) %>% 
  nest() %>% 
  mutate(wgs_model = map(data, ~lm(delta_d1_wgs ~ day, data = .x) %>% 
                           tidy)) %>% 
  unnest(wgs_model) %>% 
  filter(term == 'day') %>%
  dplyr::select(replicate, strain, ds_n, ds_wgs_slope = estimate, wgs_p = p.value) 

full_baseline_df <- wgs_mip_results %>%
  dplyr::filter(strain!="N2") %>%
  dplyr::filter(!baseline, day == 1) %>%
  dplyr::select(replicate, strain, base_frq = frq, base_pubfrq = published_frq)

full_delta_df <- wgs_mip_results %>%
  dplyr::filter(strain!="N2") %>%
  dplyr::filter(!baseline, day != 17) %>%
  dplyr::select(sample:strain, frq, published_frq) %>%
  dplyr::left_join(., full_baseline_df, by = c("replicate", "strain")) %>%
  dplyr::mutate(delta_d1_wgs = frq - base_frq,
                delta_d1_mip = published_frq - base_pubfrq) %>%
  dplyr::select(sample:strain, delta_d1_wgs, delta_d1_mip)

mip_slope <- full_delta_df %>% 
  group_by(replicate, strain) %>% 
  na.omit() %>%
  nest() %>% 
  mutate(wgs_model = map(data, ~lm(delta_d1_mip ~ day, data = .x) %>% 
                           tidy)) %>% 
  unnest(wgs_model) %>% 
  filter(term == 'day') %>%
  dplyr::select(replicate, strain, mip_slope = estimate, mip_p = p.value)

mip_slop_df <- mip_slope %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(mip_slope_mean = mean(mip_slope)) 

# pubtraits <- data.table::fread("published_data/eLife_traits.txt")

wgs_slope%>%
  dplyr::group_by(strain, ds_n) %>%
  dplyr::summarise(mean_slope = mean(ds_wgs_slope)) %>%
  dplyr::left_join(., mip_slop_df, by = c("strain")) %>%
  ggplot()+
  aes(x = mip_slope_mean, y = mean_slope)+
  geom_point(size = 2, alpha = 0.7)+
  facet_grid(.~ds_n)+
  ggpubr::stat_cor(method = "spearman", label.x = 0.002, label.y = 0,inherit.aes = T,size = 2, color = "firebrick3")+
  theme_bw(26)+
  labs(x = "MIP-seq slope", y = "NNLS slope")

ggsave(filename = "../figures/downsample_SLOPEcor.pdf",height = 8, width = 16)

wgs_slope%>%
  dplyr::group_by(strain, ds_n) %>%
  dplyr::summarise(mean_slope = mean(ds_wgs_slope)) %>%
  dplyr::left_join(., mip_slop_df, by = c("strain")) %>%
  dplyr::filter(ds_n == 1) %>%
  ggplot()+
  aes(x = mip_slope_mean, y = mean_slope)+
  geom_point(size = 2, alpha = 0.7)+
  facet_grid(.~ds_n)+
  ggpubr::stat_cor(method = "spearman", label.x = 0.002, label.y = 0,inherit.aes = T,size = 2, color = "firebrick3")+
  theme_bw(26)+
  labs(x = "MIP-seq slope", y = "NNLS slope")+
  geom_abline(slope=1, intercept = 0, linetype=2, color = "firebrick3")

ggsave(filename = "../figures/downsample_SLOPEcor_1x.pdf",height = 8, width = 10)
