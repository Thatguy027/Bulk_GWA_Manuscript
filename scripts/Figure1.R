library(tidyverse)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("data/baugh/2024bootstrapINPUT.Rdata")
load("data/baugh/2024baugh_bootstrap_prediction.rda")

#################################################### 
#################################################### 
####################################################  - Switch NA to 0

# flipped_bootstrap_input[[1]][is.na(flipped_bootstrap_input[[1]])]<-0
# 
# fullGy=crossprod(flipped_bootstrap_input[[1]],flipped_bootstrap_input[[2]])
# fullGGp=crossprod(flipped_bootstrap_input[[1]])

#################################################### 
#################################################### 
####################################################  - remove sites with any NA

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


boot_df %>%
  dplyr::filter(grepl("baseline", sample)) %>% 
  dplyr::left_join(., predictions_df, by = c("sample", "strain")) %>% 
  ggplot()+
  aes(x=frq, y = bootmean)+
  geom_point()+
  facet_wrap(~sample, scales = "free")+
  geom_pointrange(aes(xmin=frq-bootse, xmax=frq+bootse))

#################################################### - combine baseline predictions with bootstrapping output

pr_boot_df <- boot_df %>%
  dplyr::left_join(., predictions_df, by = c("sample", "strain")) %>%
  tidyr::separate(sample, into = c("replicate", "day"), sep = "_", remove = F) %>%
  dplyr::mutate(baseline = ifelse(grepl("baseline", sample), T, F)) %>%
  dplyr::select(sample, replicate, day, baseline, strain, bootmean, bootse, frq) %>%
  dplyr::mutate(day = as.numeric(gsub("d", "", day)))

library(ggrepel)

plot_means <- pr_boot_df %>%
  dplyr::filter(!baseline) %>%
  dplyr::group_by(strain, day) %>%
  dplyr::mutate(mean_reps = mean(frq)) %>%
  dplyr::distinct(strain, day, mean_reps)

ggplot(plot_means)+
  aes(x = day, y = mean_reps, group = strain, label = strain)+
  geom_line()+
  geom_line(aes(color = strain), data = plot_means %>% dplyr::filter((mean_reps > 0.015)))+
  geom_text_repel(aes(color = strain),
                  data = plot_means %>% dplyr::filter(day == 17, (mean_reps > 0.015)),
                  force_pull   = 0, # do not pull toward data points
                  nudge_x      = 0.1,
                  direction    = "y",
                  angle        = 0,
                  hjust        = 0,
                  segment.size = 0.2,
                  max.iter = 1e4, max.time = 1)+
  theme(legend.position = "none")

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

wgs_mip_results %>%
  dplyr::filter(grepl("baseline", sample), strain!="N2") %>% 
  dplyr::select(sample, strain, frq, published_frq) %>%
  tidyr::gather(frq_type, value, -sample, -strain) %>%
  ggplot()+
  aes(x = value, fill = frq_type)+
  geom_density(alpha = 0.5)+
  facet_wrap(~sample) +
  theme_bw(18)

save(wgs_mip_results, file="data/baugh/2024_processedBOOTs_with_MIP.RData")

load(file="data/baugh/2024_processedBOOTs_with_MIP.RData")

wgs_mip_results %>%
  ggplot()+
  aes(x = published_frq, y = frq)+
  geom_point()+
  facet_wrap(~sample) +
  xlim(0, 0.075) +
  ylim(0,0.075) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.01, label.y = 0.07,inherit.aes = T)+
  theme_bw(18)

#################################################### - run PCA for WGS and MIP frequencies for comparison

pca_input_wgs <- wgs_mip_results %>%
  dplyr::filter(strain!= "N2") %>%
  na.omit() %>%
  dplyr::select(sample, strain, frq) %>%
  tidyr::spread(strain, frq)

wgs_pca <- prcomp(pca_input_wgs[,2:ncol(pca_input_wgs)], center = T)$rotation

pca_input_mip <- wgs_mip_results %>%
  dplyr::filter(strain!= "N2") %>%
  na.omit() %>%
  dplyr::select(sample, strain, published_frq) %>%
  tidyr::spread(strain, published_frq)

mip_pca <- prcomp(pca_input_mip[,2:ncol(pca_input_mip)], center = T)$rotation

cor(mip_pca[,1], wgs_pca[,1], method = "spearman")
plot(mip_pca[,1], wgs_pca[,1])

# plot the correlation between PC1s 
pcplot <- data.frame(strain = names(mip_pca[,1]),
           pc1mip = mip_pca[,1], 
           pc1wgs = wgs_pca[,1]) %>%
  ggplot()+
  aes(x = pc1mip, y = pc1wgs)+
  geom_point(size=2, alpha = 0.7)+
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0,inherit.aes = T, color = "firebrick3")+
  theme_bw(26)+
  labs(x = "MIP-seq PC1", y = "NNLS PC1")

ggsave(pcplot, filename = "figures/2024Baugh_mip-nnls_comparison_PC1.pdf", height = 8, width = 10)

pca_cor_df <- data.frame(pc = names(diag(cor(mip_pca, wgs_pca, method = "spearman"))), pca_cor = diag(cor(mip_pca, wgs_pca, method = "spearman"))) %>%
  dplyr::mutate(pc = as.numeric(gsub("PC", "", pc)))

pca_cor_df %>%
  ggplot()+
  aes(x = pc, y = abs(pca_cor), color = factor(sign(pca_cor)), label = round(pca_cor,3))+
  scale_color_manual(values=c("red","black"))+
  geom_point()+
  theme_bw(18)+
  theme(legend.position = "none")+
  ylim(0.825, 0.99)+
  xlim(0,9)+
  geom_text_repel(data = pca_cor_df %>% dplyr::filter(pc <= 9),
                  force_pull   = 0, # do not pull toward data points
                  nudge_x      = 0.1,
                  direction    = "y",
                  angle        = 0,
                  hjust        = 0,
                  segment.size = 0.2,
                  max.iter = 1e4, max.time = 1)


n2_shared_sites <- data.frame(strain = names(fullGGp[,"N2_N2"]),
                              shared = (fullGGp[,"N2_N2"])) %>%
  tidyr::separate(strain, c("strain","tmp"), sep = "_") %>%
  dplyr::select(-tmp)

wgs_mip_results %>%
  dplyr::mutate(delta_pub = abs(frq - published_frq),
                sign_delta = sign(frq - published_frq)) %>%
  na.omit() %>%
  dplyr::arrange(delta_pub) %>%
  dplyr::group_by(sample) %>%
  dplyr::arrange(delta_pub) %>%
  dplyr::mutate(delta_rank = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(sum_delta_rank = sum(delta_rank)) %>% 
  dplyr::left_join(., n2_shared_sites, by = "strain")%>%
  ggplot()+
  aes(x = sum_delta_rank, y = shared)+
  geom_point()

#################################################### - run PCA for WGS and MIP frequencies for comparison - as defined by publication

baseline_df <- wgs_mip_results %>%
  dplyr::filter(!baseline, day == 1) %>%
  dplyr::select(replicate, strain, base_frq = frq, base_pubfrq = published_frq)

####### - doesnt work bc lots of inf
delta_df <- wgs_mip_results %>%
  dplyr::filter(!baseline) %>%
  dplyr::select(sample:strain, frq, published_frq) %>%
  dplyr::left_join(., baseline_df, by = c("replicate", "strain")) %>%
  dplyr::mutate(delta_d1_wgs = frq/base_frq,
                delta_d1_mip = published_frq/base_pubfrq) %>%
  dplyr::select(sample:strain, delta_d1_wgs, delta_d1_mip) %>%
  dplyr::mutate(log2_wgs = log2(delta_d1_wgs),
                log2_mip = log2(delta_d1_mip)) 

pca_input_wgs <- delta_df %>%
  dplyr::filter(strain!= "N2") %>%
  na.omit() %>%
  dplyr::select(sample, strain, log2_wgs) %>%
  tidyr::spread(strain, log2_wgs)

wgs_pca <- prcomp(pca_input_wgs[,2:ncol(pca_input_wgs)], center = T)$rotation

pca_input_mip <- wgs_mip_results %>%
  dplyr::filter(strain!= "N2") %>%
  na.omit() %>%
  dplyr::select(sample, strain, published_frq) %>%
  tidyr::spread(strain, published_frq)

mip_pca <- prcomp(pca_input_mip[,2:ncol(pca_input_mip)], center = T)$rotation

cor(mip_pca[,1], wgs_pca[,1], method = "spearman")
plot(mip_pca[,1], wgs_pca[,1])

# plot the correlation between PC1s 
pcplot <- data.frame(strain = names(mip_pca[,1]),
                     pc1mip = mip_pca[,1], 
                     pc1wgs = wgs_pca[,1]) %>%
  ggplot()+
  aes(x = pc1mip, y = pc1wgs)+
  geom_point(size=2, alpha = 0.7)+
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0,inherit.aes = T, color = "firebrick3")+
  theme_bw(26)+
  labs(x = "MIP-seq PC1", y = "NNLS PC1")

ggsave(pcplot, filename = "figures/2024Baugh_mip-nnls_comparison_PC1.pdf", height = 8, width = 10)

pca_cor_df <- data.frame(pc = names(diag(cor(mip_pca, wgs_pca, method = "spearman"))), pca_cor = diag(cor(mip_pca, wgs_pca, method = "spearman"))) %>%
  dplyr::mutate(pc = as.numeric(gsub("PC", "", pc)))

pca_cor_df %>%
  ggplot()+
  aes(x = pc, y = abs(pca_cor), color = factor(sign(pca_cor)), label = round(pca_cor,3))+
  scale_color_manual(values=c("red","black"))+
  geom_point()+
  theme_bw(18)+
  theme(legend.position = "none")+
  ylim(0.825, 0.99)+
  xlim(0,9)+
  geom_text_repel(data = pca_cor_df %>% dplyr::filter(pc <= 9),
                  force_pull   = 0, # do not pull toward data points
                  nudge_x      = 0.1,
                  direction    = "y",
                  angle        = 0,
                  hjust        = 0,
                  segment.size = 0.2,
                  max.iter = 1e4, max.time = 1)


n2_shared_sites <- data.frame(strain = names(fullGGp[,"N2_N2"]),
                              shared = (fullGGp[,"N2_N2"])) %>%
  tidyr::separate(strain, c("strain","tmp"), sep = "_") %>%
  dplyr::select(-tmp)

wgs_mip_results %>%
  dplyr::mutate(delta_pub = abs(frq - published_frq),
                sign_delta = sign(frq - published_frq)) %>%
  na.omit() %>%
  dplyr::arrange(delta_pub) %>%
  dplyr::group_by(sample) %>%
  dplyr::arrange(delta_pub) %>%
  dplyr::mutate(delta_rank = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(sum_delta_rank = sum(delta_rank)) %>% 
  dplyr::left_join(., n2_shared_sites, by = "strain")%>%
  ggplot()+
  aes(x = sum_delta_rank, y = shared)+
  geom_point()

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

library(broom)
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

slope_df %>%
  dplyr::filter(strain!="N2") %>%
  ggplot()+
  aes(x = mip_slope, y = wgs_slope)+
  geom_point()+
  facet_wrap(~replicate, scales = "free")+
  ggpubr::stat_cor(method = "spearman", label.x = 0.001, label.y = 0,inherit.aes = T, color = "firebrick3")+
  theme_bw(18)

slope_plot <- slope_df %>%
  dplyr::filter(strain!="N2") %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(mip_slope_mean = mean(mip_slope, na.rm = T),
                   wgs_slope_mean = mean(wgs_slope, na.rm = T)) %>%
  dplyr::left_join(., pubtraits, by = "strain") %>%
  ggplot()+
  aes(x = mip_slope_mean, y = wgs_slope_mean)+
  geom_point(size =3, alpha = 0.7)+
  ggpubr::stat_cor(method = "spearman", label.x = 0.001, label.y = 0,inherit.aes = T, color = "firebrick3")+
  theme_bw(26)+
  labs(x = "MIP-seq slope", y = "NNLS slope")+
  geom_abline(slope=1, intercept = 0, linetype=2, color = "firebrick3")

ggsave(slope_plot, filename = "figures/mip-nnls_comparison.pdf", height = 8, width = 10)

#################################################### - Figure one plot

fig1 <- ggpubr::ggarrange(slope_plot, pcplot, nrow = 1, widths = c(1,1), align = "h")

ggsave(fig1, filename = "figures/FIG1_baugh_comparison.pdf", height = 6, width = 12)
ggsave(fig1, filename = "figures/FIG1_baugh_comparison.png", height = 6, width = 12, dpi = 300)

#################################################### downsampling 

load("data/baugh/2026downsampled_dataset.rda")
load("data/baugh/2024_processedBOOTs_with_MIP.RData")

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
  dplyr::group_by(ds_n) %>% 
  na.omit() %>%
  dplyr::summarise(dscor = cor(mean_slope,mip_slope_mean, method = "spearman")) %>%
  dplyr::filter(ds_n!=0.5)%>%
  ggplot()+
  aes(x = factor(ds_n), y = dscor)+
  geom_point(size = 3)+
  theme_bw(18)+
  theme(panel.grid.minor.x = element_blank())+
  scale_x_discrete(labels = c("0.25", "1", "3", "5", "10") # What they say
  )+
  labs(x = "Dowsampling depth", y = expression(paste("Spearman's ", rho)))+
  ylim(0.8,0.95)


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
