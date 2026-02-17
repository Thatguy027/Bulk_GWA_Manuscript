library(tidyverse)
library(RcppML)
library(broom)
library(cowplot)
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("data/baugh/2024bootstrapINPUT.Rdata")
load("data/baugh/2024baugh_bootstrap_prediction.rda")

#################################################### 

gt_with_na <- apply(flipped_bootstrap_input[[1]], 1, function(x){
  sum(is.na(x))
})

a_ct <- flipped_bootstrap_input[[2]][which(gt_with_na==0),]
gt <- flipped_bootstrap_input[[1]][which(gt_with_na==0),]

#################################################### - run baseline predictions for combining with bootstrap outputs
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

strains <- unlist(lapply(strsplit(colnames(flipped_bootstrap_input[[1]]), split = "_"), function(x){x[1]}))

bootme=apply(ab, c(1,2), function(x) mean(x))

#################################################### - convert bootstrapping output to a dataframe
boot_df <- data.frame(strain = strains, bootmean = bootme, bootse = bootse) %>%
  tidyr::gather(stat_sample, value, -strain) %>%
  tidyr::separate(stat_sample, into = c("stat","sample"), convert = T, sep = "\\.") %>%
  tidyr::spread(stat, value)
#################################################### - combine baseline predictions with bootstrapping output
pr_boot_df <- boot_df %>%
  dplyr::left_join(., predictions_df, by = c("sample", "strain")) %>%
  tidyr::separate(sample, into = c("replicate", "day"), sep = "_", remove = F) %>%
  dplyr::mutate(baseline = ifelse(grepl("baseline", sample), T, F)) %>%
  dplyr::select(sample, replicate, day, baseline, strain, bootmean, bootse, frq) %>%
  dplyr::mutate(day = as.numeric(gsub("d", "", day)))
#################################################### - load published mip frequencies and combine with wgs data set

elife_fq <- data.table::fread(file = "data/baugh/MIPseq_frequencies.txt") %>%
  dplyr::rename(strain = V1) %>%
  tidyr::gather(sample, published_frq, -strain) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(sample = ifelse(grepl("_BL", sample), gsub("_BL", "_d1_baseline", sample), sample))

wgs_mip_results <- pr_boot_df %>%
  dplyr::left_join(., elife_fq, by  = c("sample","strain"))

wgs_mip_results %>%
  dplyr::filter(grepl("baseline", sample), strain!="N2") %>% 
  ggplot()+
  aes(x = published_frq, y = frq)+
  geom_point()+
  facet_wrap(~sample) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.01, label.y = 0.07,inherit.aes = T)+
  theme_bw(18)

save(wgs_mip_results, file="data/baugh/2024_processedBOOTs_with_MIP.RData")

load(file="data/baugh/2024_processedBOOTs_with_MIP.RData")
#################################################### - calculate slopes for WGS and MIP frequencies for comparison

baseline_df <- wgs_mip_results %>%
  dplyr::filter(!baseline, day == 1) %>%
  dplyr::select(replicate, strain, base_frq = frq, base_pubfrq = published_frq)

delta_df <- wgs_mip_results %>%
  dplyr::filter(!baseline, day != 17) %>%
  dplyr::select(sample:strain, frq, published_frq) %>%
  dplyr::left_join(., baseline_df, by = c("replicate", "strain")) %>%
  dplyr::mutate(delta_d1_wgs = frq - base_frq,
                delta_d1_mip = published_frq - base_pubfrq) %>%
  dplyr::select(sample:strain, delta_d1_wgs, delta_d1_mip)


wgs_slope <- delta_df %>% 
  group_by(replicate, strain) %>% 
  nest() %>% 
  mutate(wgs_model = map(data, ~lm(delta_d1_wgs ~ day, data = .x) %>% 
                           tidy)) %>% 
  unnest(wgs_model) %>% 
  filter(term == 'day') %>%
  dplyr::select(replicate, strain, wgs_slope = estimate, wgs_p = p.value)

mip_slope <- delta_df %>% 
  group_by(replicate, strain) %>% 
  na.omit() %>%
  nest() %>% 
  mutate(wgs_model = map(data, ~lm(delta_d1_mip ~ day, data = .x) %>% 
                           tidy)) %>% 
  unnest(wgs_model) %>% 
  filter(term == 'day') %>%
  dplyr::select(replicate, strain, mip_slope = estimate, mip_p = p.value)

pubtraits <- data.table::fread("data/baugh/eLife_traits.txt")

slope_df <- dplyr::left_join(wgs_slope, mip_slope, by = c("replicate", "strain"))

slope_plot <- slope_df %>%
  dplyr::filter(strain!="N2") %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(mip_slope_mean = mean(mip_slope, na.rm = T),
                   wgs_slope_mean = mean(wgs_slope, na.rm = T)) %>%
  dplyr::left_join(., pubtraits, by = "strain") %>%
  ggplot()+
  aes(x = mip_slope_mean, y = wgs_slope_mean)+
  geom_point(size =3, alpha = 0.7)+
  ggpubr::stat_cor(method = "spearman", label.x = 0.001, label.y = 0,inherit.aes = T, color = "firebrick3",p.accuracy = 0.0001, r.accuracy = 0.01)+
  theme_bw(26)+
  labs(x = "MIP-seq slope", y = "NNLS slope")+
  geom_abline(slope=1, intercept = 0, linetype=2, color = "firebrick3")

ggsave(slope_plot, filename = "figures/mip-nnls_comparison.pdf", height = 8, width = 10)

#################################################### downsampling 

load("data/baugh/2026downsampled_dataset.rda")
load("data/baugh/2024_processedBOOTs_with_MIP.RData")

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

wgs_slope%>%
  dplyr::group_by(strain, ds_n) %>%
  dplyr::summarise(mean_slope = mean(ds_wgs_slope)) %>%
  dplyr::left_join(., mip_slop_df, by = c("strain")) %>%
  ggplot()+
  aes(x = mip_slope_mean, y = mean_slope)+
  geom_point(size = 2, alpha = 0.7)+
  facet_grid(.~ds_n)+
  ggpubr::stat_cor(method = "spearman", label.x = 0.002, label.y = 0,inherit.aes = T,size = 2, color = "firebrick3",p.accuracy = 0.0001, r.accuracy = 0.01)+
  theme_bw(26)+
  labs(x = "MIP-seq slope", y = "NNLS slope")

# ggsave(filename = "../figures/downsample_SLOPEcor.pdf",height = 8, width = 16)

ds_cor <- wgs_slope%>%
  dplyr::group_by(strain, ds_n) %>%
  dplyr::summarise(mean_slope = mean(ds_wgs_slope)) %>%
  dplyr::left_join(., mip_slop_df, by = c("strain")) %>%
  dplyr::group_by(ds_n) %>% 
  na.omit() %>%
  dplyr::summarise(dscor = cor(mean_slope,mip_slope_mean, method = "spearman")) %>%
  dplyr::filter(ds_n!=0.5)%>%
  ggplot()+
  aes(x = factor(ds_n), y = dscor)+
  geom_point(size = 5)+
  theme_bw(26)+
  theme(panel.grid.minor.x = element_blank())+
  scale_x_discrete(labels = c("0.25", "1", "3", "5", "10") # What they say
  )+
  labs(x = "Dowsampling depth", y = expression(paste("Spearman's ", rho)))+
  ylim(0.8,0.95)

cowplot::plot_grid(slope_plot, ds_cor, nrow = 2, labels = c("A", "B"), label_size = 26)

ggsave(filename = "plots/figure1.pdf", height = 16, width = 10)  
ggsave(filename = "plots/figure1.png", height = 16, width = 10)  


#################################################### get correlations for each sample and plot distribution of correlations

sample_correlations <- wgs_mip_results %>%
  group_by(sample) %>%  # Or group_by(day, replicate) if sample isn't unique
  summarize(
    correlation = cor(frq, published_frq, method = "spearman", use = "complete.obs")
  )

# 2. Plot the histogram of these correlations
corplt <- ggplot(sample_correlations, aes(x = correlation)) +
  geom_histogram(binwidth = 0.02, fill = "#69b3a2", color = NA) +
  geom_vline(aes(xintercept = median(correlation)), color = "firebrick3", linetype = "dashed") +
  labs(
    # title = "Distribution of Correlations across Samples",
    # subtitle = paste("Median Correlation:", round(median(sample_correlations$correlation), 3)),
    x = expression(paste("Spearman's ", rho)),
    y = "Count"
  ) +
  theme_bw(26)

# The syntax is just "Top / (Bottom_Left + Bottom_Right)"
slope_plot / (corplt + ds_cor) + plot_layout(heights = c(3, 2)) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 26))

ggsave(filename = "plots/figure1_wcor.pdf", height = 12, width = 12)  
ggsave(filename = "plots/figure1_wcor.png", height = 12, width = 12)  
