library(tidyverse)
library("GGally")
library(ggrepel)
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("data/cluster_data/20230410_cleanGT_flippedCommon_n_Counts.Rdata")
# load("cluster_data/pos1RNAi_bootstrap_prediction.rda")

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

meta_df <- data.table::fread(file = "meta_files/experiment_meta.tsv") %>%
  dplyr::mutate(sample = paste0("S",`SAMPLE NUMBER`))


processed_strain_frq <- predictions_df %>%
  dplyr::left_join(., meta_df, by = "sample")

save(processed_strain_frq, file = "data/processed_strain_fq.rda")
load(file = "data/processed_strain_fq.rda")

t0 <- processed_strain_frq %>%
  dplyr::filter(Timepoint == 0) 

t0 %>%
  dplyr::select(strain, sample, frq) %>%
  tidyr::pivot_wider(names_from = sample, values_from = frq) %>%
  ggplot()+
  aes(x = S1, y = S2)+
  geom_point() +
  theme_bw(18)

mt0 <- t0 %>%
  dplyr::select(strain, sample, frq) %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(t0frq = mean(frq))

# t2 - t0   ---------     # t2 - t0   ---------     # t2 - t0   ---------     # t2 - t0   ---------     # t2 - t0   ---------     # t2 - t0

t2 <- processed_strain_frq %>%
  dplyr::filter(Timepoint == 2) %>%
  dplyr::left_join(., mt0, by = "strain") %>%
  dplyr::mutate(delta_fq = frq - t0frq)

delta_t0_traits <- t2 %>%
  tidyr::unite(trait, sample_type, Gene, sample, sep = "_") %>%
  dplyr::select(strain, trait, delta_fq) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(nzero = length(which(delta_fq==0))) %>% 
  dplyr::filter(nzero < 30) %>% # remove strains that are 0 in more than 30 of 39 traits tested
  tidyr::pivot_wider(names_from = trait, values_from = delta_fq)  %>%
  dplyr::select(-nzero,-`supernatant_mig-6_S48`)

write.table(delta_t0_traits, file = "cluster_data/delta_t0_traits.tsv", quote = F, col.names = T, row.names = F, sep = "\t")

dir.create("plots")

unique(t2$Gene)
unique(t2$sample_type)
outlier_strains <- c("XZ1516","JU2466","EG4946","PX179","JU2578","QX1791","DL226","NIC262","CX11276","MY2741","ECA396","XZ1513")
ctrl_outliers_high <- c("XZ1516","JU2578","ED3048","JU1088","DL226","ED3049")
ctrl_outliers_low <- c("ED3011","JU2001","JU3127","JU3125","XZ1513","NIC260","NIC251","NIC262")

# t2 - t0

for(clones in unique(t2$Gene)){
  for(sample_prep in unique(t2$sample_type)){
    plt_df <- t2 %>%
      dplyr::filter(Gene == clones, sample_type == sample_prep) %>%
      dplyr::select(strain, sample, delta_fq) %>%
      dplyr::filter(!sample%in%c("S19","S20"))%>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(plot_axes = as.character(cur_group_id())) %>%
      dplyr::ungroup() %>%
      dplyr::select(-sample) %>%
      na.omit()
    
    if(nrow(plt_df) > 100) {
      
      replicate_cor <- plt_df %>%
        tidyr::pivot_wider(names_from = plot_axes, values_from = delta_fq) %>%
        dplyr::mutate(cor_rep = cor(`1`,`2`,method = "spearman")) %>%
        dplyr::distinct(cor_rep) %>%
        dplyr::pull(cor_rep) %>%
        round(digits=2)
      
      outliers <- plt_df %>% dplyr::group_by(plot_axes) %>% 
        dplyr::filter(quantile(abs(delta_fq), 0.8) < abs(delta_fq)) %>% 
        tidyr::pivot_wider(names_from = plot_axes, values_from = delta_fq) 
      
      plot_df_1 <- plt_df %>%
        tidyr::pivot_wider(names_from = plot_axes, values_from = delta_fq) 
      
      plot_save <- plot_df_1 %>%
        ggplot()+
        geom_point(aes(x = `1`, y = `2`))+
        geom_point(aes(x = `1`, y = `2`), color = "darkorchid", data = plot_df_1 %>% dplyr::filter(strain%in%ctrl_outliers_high), size = 3)+
        geom_point(aes(x = `1`, y = `2`), color = "cadetblue3", data = plot_df_1 %>% dplyr::filter(strain%in%ctrl_outliers_low), size = 3)+
        # geom_text_repel(aes(x = `1`, y = `2`, label = strain), data = plot_df_1 %>% dplyr::filter(strain%in%outlier_strains) , color = "blue",
        #                 min.segment.length = 0, seed = 42, box.padding = 0.5) +
        geom_text_repel(aes(x = `1`, y = `2`, label = strain), data = outliers,min.segment.length = 0, seed = 42, box.padding = 0.5, max.overlaps = 20) +
        geom_point(aes(x = `1`, y = `2`), color = "red", data = outliers, size = 2)+
        theme_bw(18)+
        labs(title = glue::glue("{clones} {sample_prep} Delta Frq: spearman's rho = {replicate_cor}"))
      
      ggsave(plot_save, filename = glue::glue("plots/{clones}_{sample_prep}_delta_repCor_t2-t0.pdf", height = 6, width=6))
    }

  }
}

# t2 - t2c   ---------     # t2 - t2c   ---------     # t2 - t2c   ---------     # t2 - t2c   ---------     # t2 - t2c   ---------     # t2 - t2c

timepoint = 2

t2c <- processed_strain_frq %>%
  dplyr::filter(Gene == "ctrl", Timepoint == timepoint) 

# 
# t2c %>%
#   dplyr::select(strain, sample, sample_type, frq) %>%
#   tidyr::pivot_wider(names_from = sample, values_from = frq)
#   ggplot()+
#   aes(x = sample)

mt2c <- t2c %>%
  dplyr::select(strain, sample, sample_type, frq) %>%
  dplyr::group_by(strain, sample_type) %>%
  dplyr::summarise(ctrlfrq = mean(frq))

outliers_h <- mt2c %>% dplyr::group_by(sample_type) %>% 
  dplyr::filter(quantile(abs(ctrlfrq), 0.85) < abs(ctrlfrq)) %>% 
  tidyr::pivot_wider(names_from = sample_type, values_from = ctrlfrq) 

outliers_l <- mt2c %>% dplyr::group_by(sample_type) %>% 
  dplyr::filter(quantile(abs(ctrlfrq), 0.15) > abs(ctrlfrq)) %>% 
  tidyr::pivot_wider(names_from = sample_type, values_from = ctrlfrq) 

mt2c %>%
  tidyr::spread(sample_type, ctrlfrq) %>%
  ggplot()+
  aes(x = pellet, y= supernatant)+
  geom_text_repel(aes(x = pellet, y = supernatant, label = strain), data = outliers_h,min.segment.length = 0, seed = 42, box.padding = 0.5, max.overlaps = 20) +
  geom_point()

t2 <- processed_strain_frq %>%
  dplyr::filter(Timepoint == timepoint, Gene!="ctrl") %>%
  dplyr::left_join(., mt2c, by = c("strain", "sample_type")) %>%
  dplyr::mutate(delta_fq = frq - ctrlfrq)

pres_plot <- t2 %>%
  dplyr::filter(Gene %in% c("pos-1", "mig-6")) %>%
  dplyr::group_by(strain, Gene) %>%
  dplyr::mutate(mean_delta = mean(delta_fq, na.rm =T)) %>%
  dplyr::distinct(strain, Gene, mean_delta)

pres_plot%>%
  dplyr::filter(Gene == "mig-6") %>%
  ggplot()+
  aes(x = mean_delta )+
  geom_histogram() +
  geom_vline(aes(xintercept = mean_delta), data = pres_plot %>% dplyr::filter(Gene == "mig-6", strain == "JU1793"), color = "#9794E7", size = 2) +
  geom_vline(aes(xintercept = mean_delta), data = pres_plot %>% dplyr::filter(Gene == "mig-6", strain == "JU2466"), color = "#91D29C", size = 2) +
  theme_bw(18)+
  labs(x = "Change in strain frequency", y = "Count")

ggsave(filename = "mig6_JuHighlight.pdf", height = 5, width = 6)

pres_plot%>%
  dplyr::filter(Gene == "pos-1") %>%
  ggplot()+
  aes(x = mean_delta )+
  geom_histogram() +
  geom_vline(aes(xintercept = mean_delta), data = pres_plot %>% dplyr::filter(Gene == "pos-1", strain == "JU1793"), color = "#9794E7", size = 2) +
  geom_vline(aes(xintercept = mean_delta), data = pres_plot %>% dplyr::filter(Gene == "pos-1", strain == "JU2466"), color = "#91D29C", size = 2) +
  theme_bw(18)+
  labs(x = "Change in strain frequency", y = "Count")

ggsave(filename = "pos1_JuHighlight.pdf", height = 5, width = 6)


delta_t2_traits <- t2 %>%
  tidyr::unite(trait, sample_type, Gene, sample, sep = "_") %>%
  dplyr::select(strain, trait, delta_fq) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(nzero = length(which(delta_fq==0))) %>% 
  dplyr::filter(nzero < 30) %>% # remove strains that are 0 in more than 30 of 39 traits tested
  tidyr::pivot_wider(names_from = trait, values_from = delta_fq) %>%
  dplyr::select(-nzero,-`supernatant_mig-6_S48`)


pos1focus <- t2 %>%
  dplyr::filter(`SAMPLE NUMBER`%in%c(17,38,49)) %>%
  tidyr::unite(trait, sample_type, Gene, sample, sep = "_") %>%
  dplyr::select(strain, trait, delta_fq) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(nzero = length(which(delta_fq==0))) %>% 
  dplyr::filter(nzero < 30) %>% # remove strains that are 0 in more than 30 of 39 traits tested
  tidyr::pivot_wider(names_from = trait, values_from = delta_fq) %>%
  dplyr::select(-nzero)

t2_pos <-data.frame(strain = pos1focus$strain, 
                   prcomp(as.matrix(pos1focus[,2:ncol(pos1focus)]))$x)

pos_traits <- pos1focus %>%
  dplyr::left_join(., t2_pos, by = "strain")

poscor <- cor(pos_traits[,2:ncol(pos_traits)])
library(corrplot)
corrplot(poscor, method="circle")


write.table(pos_traits, file = "cluster_data/pos1_focus_traits.tsv", quote = F, col.names = T, row.names = F, sep = "\t")

# in case there is a general RNAi trait - guessing it might be on chromosome II
t2_pc <-data.frame(strain = delta_t2_traits$strain, 
                   prcomp(as.matrix(delta_t2_traits[,2:ncol(delta_t2_traits)]))$x)

t2_pc %>%
  ggplot()+
  aes(x = PC1, y = PC)+
  geom_point()


write.table(t2_pc[,1:4], file = "cluster_data/delta_t2_PCtraits_fix.tsv", quote = F, col.names = T, row.names = F, sep = "\t")

write.table(delta_t2_traits, file = "cluster_data/delta_t2_traits_fix.tsv", quote = F, col.names = T, row.names = F, sep = "\t")


# t2 - t2pos   ---------     # t2 - t2pos   ---------     # t2 - t2pos   ---------     # t2 - t2pos   ---------     # t2 - t2pos   ---------     # t2 - t2pos

timepoint = 2

t2c <- processed_strain_frq %>%
  dplyr::filter(Gene == "pos-1", Timepoint == timepoint) 

# 
# t2c %>%
#   dplyr::select(strain, sample, sample_type, frq) %>%
#   tidyr::pivot_wider(names_from = sample, values_from = frq)
#   ggplot()+
#   aes(x = sample)

mt2c <- t2c %>%
  dplyr::select(strain, sample, sample_type, frq) %>%
  dplyr::group_by(strain, sample_type) %>%
  dplyr::summarise(ctrlfrq = mean(frq))

outliers_h <- mt2c %>% dplyr::group_by(sample_type) %>% 
  dplyr::filter(quantile(abs(ctrlfrq), 0.85) < abs(ctrlfrq)) %>% 
  tidyr::pivot_wider(names_from = sample_type, values_from = ctrlfrq) 

outliers_l <- mt2c %>% dplyr::group_by(sample_type) %>% 
  dplyr::filter(quantile(abs(ctrlfrq), 0.15) > abs(ctrlfrq)) %>% 
  tidyr::pivot_wider(names_from = sample_type, values_from = ctrlfrq) 

mt2c %>%
  tidyr::spread(sample_type, ctrlfrq) %>%
  ggplot()+
  aes(x = pellet, y= supernatant)+
  geom_text_repel(aes(x = pellet, y = supernatant, label = strain), data = outliers_h,min.segment.length = 0, seed = 42, box.padding = 0.5, max.overlaps = 20) +
  geom_point()

t2 <- processed_strain_frq %>%
  dplyr::filter(Timepoint == timepoint, ! Gene%in%c("ctrl", "pos-1")) %>%
  dplyr::left_join(., mt2c, by = c("strain", "sample_type")) %>%
  dplyr::mutate(delta_fq = frq - ctrlfrq)

delta_t2_traits <- t2 %>%
  tidyr::unite(trait, sample_type, Gene, sample, sep = "_") %>%
  dplyr::select(strain, trait, delta_fq) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(nzero = length(which(delta_fq==0))) %>% 
  dplyr::filter(nzero < 30) %>% # remove strains that are 0 in more than 30 of 39 traits tested
  tidyr::pivot_wider(names_from = trait, values_from = delta_fq) %>%
  dplyr::select(-nzero,-`supernatant_mig-6_S48`)


# in case there is a general RNAi trait - guessing it might be on chromosome II
t2_pc <-data.frame(strain = delta_t2_traits$strain, 
                   prcomp(as.matrix(delta_t2_traits[,2:ncol(delta_t2_traits)]))$x)

t2_pc %>%
  ggplot()+
  aes(x = PC1, y = PC2)+
  geom_point()


t2pos <- dplyr::left_join(delta_t2_traits, t2_pc[,1:4], by = "strain")

write.table(t2pos, file = "cluster_data/delta_t2pos_traits_fix.tsv", quote = F, col.names = T, row.names = F, sep = "\t")








ctrl_outliers_high <- outliers_h %>% dplyr::pull(strain)
ctrl_outliers_low <- outliers_l %>% dplyr::pull(strain)

for(clones in unique(t2$Gene)){
  for(sample_prep in unique(t2$sample_type)){
    plt_df <- t2 %>%
      dplyr::filter(Gene == clones, sample_type == sample_prep) %>%
      dplyr::select(strain, sample, delta_fq) %>%
      dplyr::filter(!sample%in%c("S19","S20"))%>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(plot_axes = as.character(cur_group_id())) %>%
      dplyr::ungroup() %>%
      dplyr::select(-sample) %>%
      na.omit()
    
    if(nrow(plt_df) > 100) {
      
      replicate_cor <- plt_df %>%
        tidyr::pivot_wider(names_from = plot_axes, values_from = delta_fq) %>%
        dplyr::mutate(cor_rep = cor(`1`,`2`,method = "spearman")) %>%
        dplyr::distinct(cor_rep) %>%
        dplyr::pull(cor_rep) %>%
        round(digits=2)
      
      outliers <- plt_df %>% dplyr::group_by(plot_axes) %>% 
        dplyr::filter(quantile(abs(delta_fq), 0.9) < abs(delta_fq)) %>% 
        tidyr::pivot_wider(names_from = plot_axes, values_from = delta_fq) 
      
      plot_df_1 <- plt_df %>%
        tidyr::pivot_wider(names_from = plot_axes, values_from = delta_fq) 
      
      outlier_df <- plot_df_1 %>%
         dplyr::filter(strain %in% c(outliers$strain, "JU1793", "AB1","JU2466"))
      
      plot_save <- plot_df_1 %>%
        ggplot()+
        geom_point(aes(x = `1`, y = `2`))+
        geom_point(aes(x = `1`, y = `2`), color = "darkorchid", data = plot_df_1 %>% dplyr::filter(strain%in%ctrl_outliers_high), size = 3)+
        geom_point(aes(x = `1`, y = `2`), color = "cadetblue3", data = plot_df_1 %>% dplyr::filter(strain%in%ctrl_outliers_low), size = 3)+
        # geom_text_repel(aes(x = `1`, y = `2`, label = strain), data = plot_df_1 %>% dplyr::filter(strain%in%outlier_strains) , color = "blue",
        #                 min.segment.length = 0, seed = 42, box.padding = 0.5) +
        geom_text_repel(aes(x = `1`, y = `2`, label = strain), data = outlier_df,min.segment.length = 0, seed = 42, box.padding = 0.5, max.overlaps = 20) +
        geom_point(aes(x = `1`, y = `2`), color = "red", data = outlier_df, size = 2)+
        theme_bw(18)+
        labs(title = glue::glue("{clones} {sample_prep} Delta Frq: spearman's rho = {replicate_cor}"))
      
      ggsave(plot_save, filename = glue::glue("plots/{clones}_{sample_prep}_delta_control_repCor_delta_T2.pdf", height = 6, width=6))
    }
    
  }
}

rnaiclone <- "vha-5"

pos1_values <- t2 %>% 
  dplyr::filter(sample_type == "pellet") %>%
  dplyr::filter(Gene == "pos-1") %>%
  dplyr::select(strain, sample_type, Timepoint, pos1d = delta_fq, pos1s = `SAMPLE NUMBER`)

# vha-5
vha5 <- t2 %>% 
  dplyr::filter(sample_type == "pellet") %>%
  dplyr::filter(Gene == rnaiclone) %>%
  dplyr::left_join(., pos1_values, by = c("strain", "sample_type", "Timepoint"))


vha5_bigd <- vha5 %>%
  dplyr::filter(sample_type == "pellet") %>%
  dplyr::filter(!(frq == 0 & ctrlfrq == 0)) %>%
  dplyr::filter(quantile((delta_fq), 0.87) < (delta_fq) ) 

vha5 %>%
  dplyr::filter(sample_type == "pellet") %>%
  ggplot()+
  aes(x = ctrlfrq, y = frq)+
  geom_point()+
  facet_grid(pos1s~`SAMPLE NUMBER`)+
  geom_text_repel(aes(label = strain), data = vha5_bigd,min.segment.length = 0, seed = 42, box.padding = 0.5, max.overlaps = 20) 
  
vha5 %>%
  dplyr::filter(sample_type == "pellet") %>%
  ggplot()+
  aes(x = pos1d, y = delta_fq)+
  geom_point()+
  facet_grid(pos1s~`SAMPLE NUMBER`)+
  geom_text_repel(aes(label = strain), data = vha5_bigd,min.segment.length = 0, seed = 42, box.padding = 0.5, max.overlaps = 20) 


processed_strain_frq %>% 
  dplyr::filter(Gene == "vha-5", `SAMPLE NUMBER` %in% c(14, 25, 35, 46)) %>%
  dplyr::select(strain, frq, sample_type) %>%
  tidyr::pivot_wider(names_from = sample_type, values_from = frq) %>%
  ggplot()+
  aes(x = pellet, y = supernatant)+
  geom_point()


focal_strain <- t2 %>%
  dplyr::filter(sample_type == "pellet") %>%
  dplyr::mutate(faceting = paste(Gene, sample_type, sample, sep = "_")) %>%
  dplyr::select(strain, sample, delta_fq,faceting) %>%
  dplyr::filter(strain == "JU1793")

t2 %>%
  dplyr::filter(sample_type == "pellet") %>%
  na.omit() %>%
  dplyr::mutate(faceting = paste(Gene, sample_type, sample, sep = "_")) %>%
  dplyr::select(strain, sample, delta_fq,faceting) %>%
  ggplot()+
  aes(x = delta_fq) +
  geom_histogram()+
  geom_vline(aes(xintercept = delta_fq), data = focal_strain, color = "red")+
  facet_wrap(vars(faceting), scales = "free")


# mig-6
processed_strain_frq %>% 
  dplyr::filter(Gene == "mig-6", `SAMPLE NUMBER` %in% c(37, 48)) %>%
  dplyr::select(strain, frq, sample_type) %>%
  tidyr::pivot_wider(names_from = sample_type, values_from = frq) %>%
  ggplot()+
  aes(x = pellet, y = supernatant)+
  geom_point()


judf <- t2 %>%
  dplyr::mutate(faceting = paste(Gene, sample_type, sample, sep = "_")) %>%
  dplyr::select(strain, sample, delta_fq,faceting) %>%
  dplyr::filter(strain == "JU1793")

t2 %>%
  na.omit() %>%
  dplyr::mutate(faceting = paste(Gene, sample_type, sample, sep = "_")) %>%
  dplyr::select(strain, sample, delta_fq,faceting) %>%
  ggplot()+
  aes(x = delta_fq) +
  geom_histogram()+
  geom_vline(aes(xintercept = delta_fq), data = judf, color = "red")+
  facet_wrap(vars(faceting), scales = "free")




t2 %>%
  dplyr::group_by(sample) %>%
  dplyr::arrange(desc(delta_fq)) %>% View()
  dplyr::



contrast_samples <- c("S2","S20")

processed_strain_frq %>%
  dplyr::filter(sample %in% contrast_samples) %>% 
  ggplot()+
  aes(x = sample, y = frq, group = strain)+
  geom_line()

delta_df <- processed_strain_frq %>%
  dplyr::filter(sample %in% contrast_samples) %>% 
  dplyr::select(strain:frq) %>%
  tidyr::spread(sample, frq)

colnames(delta_df) <- c("strain", "t0", "t1")

delta_df %>%
  dplyr::mutate(delta_fq = t1-t0) %>%
  dplyr::arrange(desc(delta_fq))

contrast_samples <- c("S1","S51") # many outliers
contrast_samples <- c("S1","S50") # two outliers


# strains missing in all timepoints and experiments
# missin_strains <- processed_strain_frq %>%
#   dplyr::filter(sample!="A1") %>%
#   dplyr::group_by(strain) %>%
#   dplyr::mutate(sum_fq = sum(frq)) %>% 
#   dplyr::filter(sum_fq == 0) %>%
#   dplyr::pull(strain) %>%
#   unique()
# 
# processed_strain_frq <- processed_strain_frq %>%
#   dplyr::filter(!(strain %in% missin_strains)) # remove strains that are not detected in any experiment/timepoint

bootme=apply(ab, c(1,2), function(x) mean(x))

boot_df <- data.frame(strain = colnames(fullGGp), bootmean = bootme, bootse = bootse) %>%
  tidyr::gather(stat_sample, value, -strain) %>%
  tidyr::separate(stat_sample, into = c("stat","sample"), convert = T, sep = "\\.") %>%
  tidyr::spread(stat, value) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(strain = strsplit(strain,split = "_")[[1]][1]) 

processed_strain_frq <- processed_strain_frq %>%
  dplyr::left_join(., boot_df, by = c("sample","strain"))

write.table(processed_strain_frq, file = "processed_data/pr_strain_frq.tsv", sep = "\t", col.names = F, row.names = F)

processed_strain_frq <- data.table::fread("processed_data/pr_strain_frq.tsv",col.names = c("strain","sample","frq","time","condition", "bootmean","bootse"))



pltdf <- processed_strain_frq %>% dplyr::filter(sample == "A1" | sample == "A13") %>% dplyr::select(strain, sample, frq, bootmean, bootse)

#highest increasing strains
incstr <- pltdf %>% dplyr::select(strain, sample, frq ) %>% 
  tidyr::spread(sample, frq) %>% 
  dplyr::mutate(delta = A13 - A1) %>% 
  dplyr::arrange(desc(delta)) %>% 
  dplyr::filter(delta>quantile(delta,0.9)) %>%
  dplyr::pull(strain)

# pos-1
pltdf%>% 
  ggplot() +
  aes(x = sample, y = frq, group = strain, label = strain)+
  geom_line(alpha = 0.5)+
  geom_line(color = "blue", data = pltdf %>% dplyr::filter(strain %in% incstr))+
  geom_point(alpha = 0.5)+
  geom_errorbar(aes(ymin=frq-bootse, ymax=frq+bootse), width=0, alpha = 0.5)+
  geom_text_repel(
    nudge_x = 0.1, direction = "y", hjust = "left", data = pltdf %>% dplyr::filter(sample == "A13", strain %in% incstr)
  )+
  theme_bw()

# ctrl

pltdf <- processed_strain_frq %>% dplyr::filter(sample == "A1" | sample == "A7") %>% dplyr::select(strain, sample, frq, bootmean, bootse)

pltdf%>% 
  ggplot() +
  aes(x = sample, y = frq, group = strain, label = strain)+
  geom_line(alpha = 0.5)+
  geom_line(color = "blue", data = pltdf %>% dplyr::filter(strain %in% incstr))+
  geom_point(alpha = 0.5)+
  geom_errorbar(aes(ymin=frq-bootse, ymax=frq+bootse), width=0, alpha = 0.5)+
  # geom_text_repel(
  #   nudge_x = 0.1, direction = "y", hjust = "left", data = pltdf %>% dplyr::filter(sample == "A7", strain %in% incstr)
  # )+
  theme_bw()



samp=c("A10")
pltdf <- processed_strain_frq %>% 
  dplyr::filter(sample == "A1" | sample %in% samp) %>% 
  dplyr::select(strain, sample, frq, bootmean, bootse) 
pltdf%>% 
  dplyr::filter(strain %in% c("ECA760", "ECA738", "ECA1255", "JU2522")) %>%
  ggplot() +
  aes(x = sample, y = frq, group = strain, label = strain)+
  geom_line(alpha = 0.5)+
  # geom_line(color = "blue", data = pltdf %>% dplyr::filter(strain %in% incstr))+
  geom_point(alpha = 0.5)+
  geom_errorbar(aes(ymin=frq-bootse, ymax=frq+bootse), width=0, alpha = 0.5)+
  geom_text_repel(
    nudge_x = 0.1, direction = "y", hjust = "left", data = pltdf %>% dplyr::filter(sample == samp, strain %in% c("ECA760", "ECA738", "ECA1255", "JU2522"))
  )+
  theme_bw()

samp=c("A4", "A10")
pltdf <- processed_strain_frq %>% 
  dplyr::filter(sample == "A1" | sample %in% samp) %>% 
  dplyr::select(strain, sample, frq, bootmean, bootse) %>%
  dplyr::mutate(sample_n = as.numeric(gsub("A","",sample)))

pltdf%>% 
  # dplyr::filter(strain %in% c("ECA760", "ECA738", "ECA1255", "JU2522")) %>%
  ggplot() +
  aes(x = sample_n, y = frq, group = strain, label = strain)+
  geom_line(alpha = 0.5)+
  # geom_line(color = "blue", data = pltdf %>% dplyr::filter(strain %in% incstr))+
  geom_point(alpha = 0.5)+
  theme_bw()+
  ylim(0,0.05)



repa <- c("A10", "A11","A12","A13", "A6", "A7")

processed_strain_frq %>%
  dplyr::filter(sample %in% repa) %>%
  dplyr::select(strain:frq) %>%
  tidyr::spread(sample, frq) %>%
  ggpairs(2:7)

# sample correlation high for pos1 reps
repa <- c("A10", "A11","A12","A13")
processed_strain_frq %>%
  dplyr::filter(sample %in% repa) %>%
  dplyr::select(strain:frq) %>%
  tidyr::spread(sample, frq) %>%
  ggpairs(2:5)
  
  aes(x = sample, y = frq, color = sample)+
  geom_point()+
  geom_errorbar(aes(ymin=frq-bootse, ymax=frq+bootse), width=0,
                position=position_dodge(0.05))+
  facet_grid(sample~.)


repa <- c("A10", "A11","A12","A13", "A1")

mean_rnai <- processed_strain_frq %>%
  dplyr::filter(sample %in% repa) %>%
  dplyr::select(strain:frq) %>%
  tidyr::spread(sample, frq) %>%
  dplyr::rowwise() %>%
  dplyr::filter(sum(A1,A10,A11,A12,A13)!=0) %>%
  dplyr::mutate(mean_rnai = sum(A10,A11,A12,A13)/4) %>%
  dplyr::mutate(delta_baseline = mean_rnai-A1) %>%
  dplyr::select(strain,mean_rnai, delta_baseline)
  
write.table(mean_rnai, file = "cluster_data/mean_rnai.tsv", quote = F, col.names = T, row.names = F, sep = "\t")

mean_rnai %>%
  ggplot()+
  aes(x = mean_rnai, y = delta_baseline)+
  geom_point()


processed_strain_frq %>% dplyr::filter(sample == "A2" | sample == "A3") %>% dplyr::select(strain, sample, frq, bootmean, bootse) %>% 
  ggplot() +
  aes(x = bootmean, y = frq, color = sample)+
  geom_point()+
  geom_errorbar(aes(ymin=frq-bootse, ymax=frq+bootse), width=0,
                position=position_dodge(0.05))+
  facet_grid(sample~.)


# maintain replicates of t0 for baseline frequency
t0_replicate_df <- processed_strain_frq %>%
  dplyr::filter(sample == "A1") %>%
  dplyr::select(strain, t0_frq = frq) %>%
  dplyr::distinct()

normalized_frq_df <- processed_strain_frq %>% 
  dplyr::left_join(., t0_replicate_df, by = c("strain")) %>% 
  dplyr::mutate(delta_t0 = frq - t0_frq)

# t1 comparison
ctrl_frq <- normalized_frq_df %>%
  dplyr::filter(time == "t1", grepl("ctrl", condition)) %>%
  dplyr::select(strain, frq, condition) %>%
  tidyr::spread(condition, frq)

pos1_frq <- normalized_frq_df %>%
  dplyr::filter(time == "t1", grepl("pos", condition)) %>%
  dplyr::select(strain, frq, condition) %>%
  tidyr::spread(condition, frq) %>%
  dplyr::left_join(., ctrl_frq, by = c("strain")) 

pos1_frq %>%
  dplyr::mutate(pos_mean = (pos1A_outgrowth+ pos1B_outgrowth)/2,
                ctrl_mean = (ctrlA+ ctrlB)/2) %>%
  dplyr::mutate(posA_Aratio = pos1A_outgrowth/ctrlA,
                posA_Bratio = pos1A_outgrowth/ctrlB,
                posB_Aratio = pos1B_outgrowth/ctrlA,
                posB_Bratio = pos1B_outgrowth/ctrlB,
                mean_ratio = pos_mean/ctrl_mean) %>%
  ggplot()+
  aes(y = posA_Bratio, x= mean_ratio)+
  geom_point()




# t2 comparison
ctrla_frq <- normalized_frq_df %>%
  dplyr::filter(time == "t2", grepl("ctrlA", condition)) %>%
  dplyr::select(strain, ctrla_frq = frq, ctrla_mean = bootmean, ctrla_se = bootse)

ctrlb_frq <- normalized_frq_df %>%
  dplyr::filter(time == "t2", "ctrlB" == condition) %>%
  dplyr::select(strain, ctrlb_frq = frq, ctrlb_mean = bootmean, ctrlb_se = bootse)

pos1_frq <- normalized_frq_df %>%
  dplyr::filter(time == "t2", grepl("pos", condition), !grepl("ctr1", condition),!grepl("ctrl", condition)) %>%
  dplyr::select(strain, frq, condition, t0_frq, bootmean, bootse) %>%
  # tidyr::spread(condition, frq) %>%
  dplyr::left_join(., ctrla_frq, by = c("strain")) %>%
  dplyr::left_join(., ctrlb_frq, by = c("strain"))

# t0 v t2
t0_frq <- normalized_frq_df %>%
  dplyr::filter(time == "t0") %>%
  dplyr::select(strain, t0_frq = frq, t0_mean = bootmean, t0_se = bootse)

pos1_frq <- normalized_frq_df %>%
  dplyr::filter(time == "t2", grepl("pos", condition), !grepl("ctr1", condition),!grepl("ctrl", condition)) %>%
  dplyr::select(strain, frq, condition, bootmean, bootse)%>%
  dplyr::left_join(., t0_frq, by = c("strain"))

pos1_frq %>%
  dplyr::mutate(delta_t0= t0_frq-frq) %>% 
  ggplot()+
  aes(x = t0_frq, y = frq, color = condition)+
  geom_point()+
  geom_errorbar(aes(ymin=frq-bootse, ymax=frq+bootse), width=0)+
  geom_errorbar(aes(xmin=t0_frq-t0_se, xmax=t0_frq+t0_se), width=0)+
  geom_abline(intercept = 0)
  
pos1_frq %>%
  dplyr::mutate(delta_t0= t0_frq-frq) %>% 
  ggplot()+
  aes(y = delta_t0, x = t0_frq)+
  geom_point()

rnai_resistant <- pos1_frq %>%
  dplyr::mutate(delta_t0= t0_frq-frq)%>%
  dplyr::filter(delta_t0 < -0.005) %>%
  dplyr::distinct(strain) %>%
  dplyr::pull(strain)


pos1_frq_long <- pos1_frq %>%
  tidyr::gather(t0_data, value, -(strain:bootse))

mapping_df <- pos1_frq %>%
  dplyr::mutate(delta_frq=t0_frq-frq) %>%
  dplyr::select(strain, condition, delta_frq) %>%
  dplyr::group_by(strain) %>%
  dplyr::filter(!sum(delta_frq)==0) %>% 
  tidyr::spread(condition, delta_frq)

write_tsv(mapping_df, file = "cluster_data/pos1rnai.tsv", col_names = T)


for(res_strain in unique(processed_strain_frq$strain)){
  pos1_frq_long %>%
    dplyr::filter(strain == res_strain) %>%
    tidyr::spread(t0_data, value)%>%
    ggplot()+
    geom_point(size =2, aes(x = condition, y = frq))+
    geom_point(size =2, aes(x = condition, y = t0_frq), color = "red")+
    geom_errorbar(aes(x = condition, y = frq, ymin=frq-bootse, ymax=frq+bootse), width=0)+
    geom_errorbar(aes(x = condition, y = t0_frq, ymin=t0_frq-t0_se, ymax=t0_frq+t0_se), width=0, color = "red")+
    theme_bw(18)
  
  ggsave(glue::glue("plots/{res_strain}_frq_t2.png", height= 6, width = 10))
}



resistant_strains <- c()
for(res_strain in list.files("plots/maybe/")){
  resistant_strains <- c(resistant_strains, strsplit(res_strain,split = "_")[[1]][1])
  
  # temp_df <- bi_df %>%
  #   dplyr::mutate(resistant = ifelse(strain%in%resistant_strains,T,F))
  
}

bi_df <- pos1_frq_long %>%
  dplyr::distinct(strain) %>%
  dplyr::mutate(bi_pheno = ifelse(strain %in% resistant_strains, 1, 0))

write_tsv(bi_df, file = "cluster_data/binary_pos1rnai.tsv", col_names = T)

library(milorGWAS)

geno.x <- BEDMatrix::BEDMatrix("/Users/Stefan/UCLA/Genomics_Data/CeNDR/20220216/binary_rnai/binaryRNAi")

fam_df <- read.table(file = "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20220216/binary_rnai/binaryRNAi.fam")
bim_df <- read.table(file = "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20220216/binary_rnai/binaryRNAi.bim")

K <- as.matrix(read.table(file = "/Users/Stefan/UCLA/Genomics_Data/CeNDR/20220216/binary_rnai/gemmeGRM.X.sXX.txt"))

bi_pheno <- read.table("~/UCLA/Projects/bulkGWAS/20230406_pos1-rnai/cluster_data/binary_pos1rnai.tsv",header = T)
y <- bi_pheno$bi_pheno
  
gwas.binary = milorGWAS::association.test.logistic(x=as.bed.matrix(as.matrix(geno.x)),
                                                       Y=y,
                                                       X=matrix(1,length(y),1), K=K, algorithm='amle')



pos1_frq %>%
  ggplot()+
  aes(x = strain, y = frq-t0_frq, color = condition)+
  geom_point()+
  geom_errorbar(aes(ymin=frq-bootse, ymax=frq+bootse), width=0,
                position=position_dodge(0.05))


pos1_frq %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(strain_mean = mean(frq),
                ctrl_mean = (ctrlA+ctrlB)/2) %>%
  dplyr::mutate(delta_a_mean = ctrl_mean - strain_mean,
                delta_b_mean = ctrl_mean - strain_mean,
                delta_a = ctrlA - frq,
                delta_b = ctrlB - frq) %>%
  dplyr::mutate(ratio_a_mean = strain_mean/ctrl_mean,
                ratio_b_mean = strain_mean/ctrl_mean,
                ratio_a = frq/ctrlA,
                ratio_b = frq/ctrlB) %>%
  ggplot()+
  aes(y = ratio_a_mean, x = ratio_b_mean, color = pos_cond)+
  geom_point()

  
pos1_frq %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(strain_mean = mean(frq),
                ctrl_mean = (ctrlA+ctrlB)/2) %>%
  dplyr::mutate(delta_mean = ctrl_mean - strain_mean,
                delta_a = ctrlA - frq,
                delta_b = ctrlB - frq) %>%
  dplyr::mutate(ratio_mean = strain_mean/ctrl_mean,
                ratio_a = frq/ctrlA,
                ratio_b = frq/ctrlB) %>%
  dplyr::filter(ratio_mean > 4) %>%
  dplyr::distinct(strain) %>%
  dplyr::pull(strain)

rnai_res <- pos1_frq %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(strain_mean = mean(frq),
                ctrl_mean = (ctrlA+ctrlB)/2) %>%
  dplyr::mutate(delta_mean = ctrl_mean - strain_mean,
                delta_a = ctrlA - frq,
                delta_b = ctrlB - frq) %>%
  dplyr::mutate(ratio_mean = strain_mean/ctrl_mean,
                ratio_a = frq/ctrlA,
                ratio_b = frq/ctrlB) %>%
  dplyr::filter(ratio_a > 4, ratio_b > 4) %>%
  dplyr::distinct(strain) %>%
  dplyr::pull(strain)


for(res_strain in rnai_res){
  normalized_frq_df %>%
    dplyr::filter(time %in% c("t2"), grepl("pos", condition), !grepl("ctr1", condition),!grepl("ctrl", condition)) %>%
    dplyr::select(strain, frq, condition, t0_frq) %>%
    tidyr::spread(condition, frq) %>%
    dplyr::left_join(., ctrl_frq, by = c("strain")) %>%
    tidyr::gather(pos_cond, frq, -strain) %>%
    dplyr::filter(strain == res_strain) %>%
    ggplot()+
    aes(x = pos_cond, y = frq)+
    geom_point(size =2)+
    theme_bw(18)+
    ylim(c(0,NA))
  
  ggsave(glue::glue("plots/{res_strain}_frq_t2.png", height= 6, width = 10))
  
  normalized_frq_df %>%
    dplyr::filter(time %in% c("t1"), grepl("pos", condition), !grepl("ctr1", condition),!grepl("ctrl", condition)) %>%
    dplyr::select(strain, frq, condition, t0_frq) %>%
    tidyr::spread(condition, frq) %>%
    tidyr::gather(pos_cond, frq, -strain) %>%
    dplyr::filter(strain == res_strain) %>%
    ggplot()+
    aes(x = pos_cond, y = frq)+
    geom_point(size =2)+
    theme_bw(18)+
    ylim(c(0,NA))
  
  ggsave(glue::glue("plots/{res_strain}_frq_t1.png", height= 6, width = 10))
}




t0_replicate_df <- processed_strain_frq %>%
  dplyr::filter(timepoint == 0) %>%
  dplyr::select(strain, replicate, t0_frq = frq) %>%
  dplyr::distinct() %>%
  dplyr::filter(replicate != "L4")

normalized_frq_df <- normalized_frq_df %>% 
  dplyr::left_join(., t0_replicate_df, by = c("strain", "replicate")) %>% 
  dplyr::mutate(delta_t0 = frq - t0_frq)

normalized_frq_df %>%
  dplyr::filter(grepl("L", sample)) %>%
  ggplot()+
  aes(x = timepoint, y = delta_t0, group = strain) %>%
  geom_line()+
  facet_wrap(~replicate)



top_strains <- normalized_frq_df %>%
  dplyr::filter(grepl("L", sample), timepoint == 13) %>%
  dplyr::group_by(strain) %>%
  dplyr::filter(timepoint == max(timepoint, na.rm = T)) %>%
  dplyr::arrange(desc(delta_t0)) %>%
  dplyr::distinct(strain, .keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::top_n(50) %>%
  dplyr::arrange(desc(delta_t0)) 

bot_strains <- normalized_frq_df %>%
  dplyr::filter(grepl("L", sample), timepoint == 13) %>%
  dplyr::group_by(strain) %>%
  dplyr::filter(timepoint == max(timepoint, na.rm = T)) %>%
  dplyr::arrange(desc(delta_t0)) %>%
  dplyr::distinct(strain, .keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::top_n(-1) %>%
  dplyr::arrange(desc(delta_t0)) 

outliers <- dplyr::bind_rows(top_strains, bot_strains) %>%
  dplyr::mutate(rank_order = factor(1:n())) %>%
  dplyr::select(strain, rank_order)

normalized_frq_df %>%
  dplyr::filter(grepl("L", sample)) %>%
  dplyr::left_join(., outliers, by = "strain") %>%
  dplyr::filter(!is.na(rank_order)) %>%
  dplyr::arrange(rank_order) %>%
  dplyr::mutate(strain = factor(strain, levels = unique(strain))) %>% 
  ggplot()+
  aes(x = timepoint, y = delta_t0, group = replicate) %>%
  geom_line()+
  facet_wrap(.~strain)

normalized_frq_df %>%
  dplyr::filter(grepl("_", sample) | sample == "W7") %>%
  dplyr::left_join(., outliers, by = "strain") %>%
  dplyr::filter(!is.na(rank_order)) %>%
  dplyr::arrange(rank_order) %>%
  dplyr::mutate(strain = factor(strain, levels = unique(strain))) %>% 
  ggplot()+
  aes(x = timepoint, y = delta_t1, group = replicate) %>%
  geom_line()+
  facet_wrap(.~strain)

save(normalized_frq_df, file="processed_data/20220914_wgsL1_normalized frequencies.RData")

load("processed_data/20220914_wgsL1_normalized frequencies.RData")

# run PCA
pca_input_wgs <- normalized_frq_df %>%
  dplyr::filter(!grepl("W", sample)) %>%
  dplyr::select(strain, sample, delta_t1) %>%
  dplyr::group_by(strain) %>%
  dplyr::filter(sum(delta_t1) > 0 ) %>%
  tidyr::spread(strain, delta_t1)

wgs_pca <- prcomp(pca_input_wgs[,2:ncol(pca_input_wgs)], center = T)$rotation

wgs_pca %>%
  as.data.frame() %>%
  tibble::rowid_to_column(var = "strain") %>%
  ggplot()+
  aes(x = PC1, y = PC4)+
  geom_point()

# extract strain slopes
library(broom)

wgs_slope <- normalized_frq_df %>%
  dplyr::filter(!grepl("W", sample), timepoint!=0) %>%
  dplyr::select(strain:timepoint, delta_t1) %>%
  dplyr::select(-frq) %>%
  dplyr::group_by(strain) %>%
  dplyr::filter(sum(delta_t1) > 0 ) %>%
  group_by(replicate, strain) %>% 
  nest() %>% 
  mutate(wgs_model = map(data, ~lm(delta_t1 ~ timepoint, data = .x) %>% 
                           tidy)) %>% 
  unnest(wgs_model) %>% 
  filter(term == 'timepoint') %>%
  dplyr::select(replicate, strain, wgs_slope = estimate, wgs_p = p.value)


pr_slopes <- wgs_slope %>%
  dplyr::select(-wgs_p) %>%
  tidyr::spread(replicate, wgs_slope) %>%
  dplyr::mutate(slope_mean = (L1+L2+L3)/3) 

pr_slopes%>%
  ggplot()+
  aes(x = L3, y = slope_mean)+
  geom_point()
  
delta_df <- normalized_frq_df %>%
  dplyr::filter(!grepl("W", sample)) %>%
  dplyr::select(strain, sample, delta_t1) %>%
  dplyr::group_by(strain) %>%
  dplyr::filter(sum(delta_t1) > 0 ) %>%
  tidyr::spread(sample, delta_t1) %>%
  dplyr::rename(L1_d = L1, L2_d = L2, L3_d = L3)

l1_wgs_traits <- pr_slopes %>%
  dplyr::left_join(., wgs_pca %>%
                     as.data.frame() %>%
                     tibble::rownames_to_column(var = "strain") %>%
                     dplyr::select(strain:PC5), by = "strain") %>%
  dplyr::left_join(., delta_df, by = "strain")


write.table(l1_wgs_traits, file = "cluster_data/20220921_WGS_L1_traits.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

data.table::fread("cluster_data/20220921_WGS_L1_traits.tsv") %>%
  tidyr::gather(trait, value, -strain) %>%
  dplyr::group_by(strain) %>%
  dplyr::filter(sum(value) > 0 ) %>%
  tidyr::spread(trait, value)

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

elife_fq <- data.table::fread(file = "published_data/MIPseq_frequencies.txt") %>%
  dplyr::rename(strain = V1) %>%
  tidyr::gather(sample, published_frq, -strain) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(sample = ifelse(grepl("_BL", sample), gsub("_BL", "_d1_baseline", sample), sample))

wgs_mip_results <- pr_boot_df %>%
  dplyr::left_join(., elife_fq, by  = c("sample","strain"))

# save(pr_boot_df, file = "Processed_bootstrap_strain_frequencies.rda")

wgs_mip_results %>%
  dplyr::filter(grepl("baseline", sample), strain!="N2") %>% 
  ggplot()+
  aes(x = published_frq, y = frq)+
  geom_point()+
  facet_wrap(~sample) +
  # xlim(0, 0.075) +
  # ylim(0,0.075) +
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
  # xlim(0, 0.075) +
  # ylim(0,0.075) +
  # ggpubr::stat_cor(method = "spearman", label.x = 0.01, label.y = 0.07,inherit.aes = T)+
  theme_bw(18)

save(wgs_mip_results, file="processed_data/20220908_flipCommon_processed_bootstraps.RData")

load(file="processed_data/20220908_flipCommon_processed_bootstraps.RData")

wgs_mip_results %>%
  ggplot()+
  aes(x = published_frq, y = frq)+
  geom_point()+
  facet_wrap(~sample) +
  xlim(0, 0.075) +
  ylim(0,0.075) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.01, label.y = 0.07,inherit.aes = T)+
  theme_bw(18)

pca_input_wgs <- wgs_mip_results %>%
  dplyr::filter(strain!= "N2") %>%
  na.omit() %>%
  dplyr::select(sample, strain, frq) %>%
  tidyr::spread(strain, frq)



pca_input_mip <- wgs_mip_results %>%
  dplyr::filter(strain!= "N2") %>%
  na.omit() %>%
  dplyr::select(sample, strain, published_frq) %>%
  tidyr::spread(strain, published_frq)

mip_pca <- prcomp(pca_input_mip[,2:ncol(pca_input_mip)], center = T)$rotation

cor(mip_pca[,1], wgs_pca[,1], method = "spearman")
plot(mip_pca[,1], wgs_pca[,1])

plot(abs(diag(cor(mip_pca, wgs_pca))))

heatmap(cor(mip_pca, wgs_pca), Rowv=NA, Colv=NA)

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

wgs_mip_results %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(sumfq = sum(frq, na.rm = T),
                   sumpub = sum(published_frq, na.rm = T)) %>% View()



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

pubtraits <- data.table::fread("published_data/eLife_traits.txt")

slope_df <- dplyr::left_join(wgs_slope, mip_slope, by = c("replicate", "strain"))

slope_df %>%
  dplyr::filter(strain!="N2") %>%
  ggplot()+
  aes(x = mip_slope, y = wgs_slope)+
  geom_point()+
  facet_wrap(~replicate, scales = "free")+
  ggpubr::stat_cor(method = "spearman", label.x = 0.001, label.y = 0,inherit.aes = T, color = "firebrick3")+
  theme_bw(18)

slope_df %>%
  dplyr::filter(strain!="N2") %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(mip_slope_mean = mean(mip_slope, na.rm = T),
                   wgs_slope_mean = mean(wgs_slope, na.rm = T)) %>%
  dplyr::left_join(., pubtraits, by = "strain") %>%
  ggplot()+
  aes(x = Slope, y = wgs_slope_mean)+
  geom_point()+
  ggpubr::stat_cor(method = "spearman", label.x = 0.001, label.y = 0,inherit.aes = T, color = "firebrick3")+
  theme_bw(18)

delta_df %>%
  ggplot()+
  aes(x = day, y = delta_d1_wgs, group = replicate)+
  geom_line()+
  facet_wrap(~strain)

delta_df %>%
  ggplot()+
  aes(x = day, y = delta_d1_mip, group = replicate)+
  geom_line()+
  facet_wrap(~strain)



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