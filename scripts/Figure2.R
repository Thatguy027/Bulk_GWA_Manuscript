library(tidyverse)
library(smplot2)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("data/pos1_rnai_sensitive/processed_strain_fq.rda")

# get delta traits 
timepoint = 2

# time 2 control
t2c <- processed_strain_frq %>%
  dplyr::filter(Gene == "ctrl", Timepoint == timepoint) 

# mean of time 2 control replicates
mt2c <- t2c %>%
  dplyr::select(strain, sample, sample_type, frq) %>%
  dplyr::group_by(strain, sample_type) %>%
  dplyr::summarise(ctrlfrq = mean(frq))

# get t2 traits for each rnai condition
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
  dplyr::filter(Gene %in% c("pos-1", "mig-6")) %>%
  tidyr::unite(trait, sample_type, Gene, sample, sep = "_") %>%
  dplyr::select(strain, trait, delta_fq) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(nzero = length(which(delta_fq==0))) %>% 
  dplyr::filter(nzero < 5) %>% # remove strains that are 0 in more than 30 of 39 traits tested
  tidyr::pivot_wider(names_from = trait, values_from = delta_fq) %>%
  dplyr::select(-nzero,-`supernatant_mig-6_S48`) # bad sample

# in case there is a general RNAi trait - guessing it might be on chromosome II
t2_pc <-data.frame(strain = delta_t2_traits$strain, 
                   prcomp(as.matrix(delta_t2_traits[,2:ncol(delta_t2_traits)]))$x)

delta_t2_traits_addmeans <- delta_t2_traits %>%
  dplyr::mutate(pel_pos = (`pellet_pos-1_S17`+`pellet_pos-1_S38`)/2,
                sup_pos = (`supernatant_pos-1_S27`+`supernatant_pos-1_S49`)/2,
                pel_mig = (`pellet_mig-6_S16`+`pellet_mig-6_S37`)/2)

write.table(delta_t2_traits_addmeans, file = "../data/pos1_rnai_sensitive/2026_delta_t2_traits.tsv", quote = F, col.names = T, row.names = F, sep = "\t")

# run gemma 
base_dir <- "data/pos1_rnai_sensitive/"
trait_file <- "2026_delta_t2_traits.tsv"

trait_values <- data.table::fread(glue::glue("{base_dir}{trait_file}")) %>%
  tidyr::gather(trait, value, -strain)

# the chrx output uses all snps for K
# loco_files <- grep("I\\.assoc|II\\.assoc|III\\.assoc|IV\\.assoc|V\\.assoc|X\\.assoc",list.files(glue::glue("{base_dir}/gemma/")), value=T)
loco_files <- grep("X\\.assoc",list.files(glue::glue("{base_dir}/gemma/")), value=T)

trait_id <- colnames(data.table::fread(glue::glue("{base_dir}{trait_file}")))
trait_id  <- data.frame(trait = trait_id[2:length(trait_id)], id = 1:(length(trait_id)-1))


for(locochr in loco_files){
  
  # print(locochr)
  
  left_out_chr = strsplit(locochr, split = "\\.")[[1]][3]
  trait_n = trait_id %>%
    dplyr::filter(id == as.numeric(strsplit(locochr, split = "\\.")[[1]][2])) %>%
    dplyr::pull(trait)
  
  if(left_out_chr == "X"){
    # chr_res <- data.table::fread(glue::glue("{base_dir}/gemma/{locochr}"))%>%
    #   dplyr::filter(chr=="23") %>%
    #   dplyr::mutate(trait = trait_n,
    #                 chr="X")
    chr_res <- data.table::fread(glue::glue("{base_dir}/gemma/{locochr}"))%>%
      dplyr::mutate(trait = trait_n)
      } else {
    chr_res <- data.table::fread(glue::glue("{base_dir}/gemma/{locochr}"))%>%
      dplyr::filter(chr==left_out_chr) %>%
      dplyr::mutate(trait = trait_n)
  }
  
  
  if(!exists("loco_res")){
    loco_res <- chr_res
  } else {
    loco_res <- dplyr::bind_rows(loco_res, chr_res)
  }
}

facet_s <- loco_res %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(min_chr = min(ps),
                   max_chr = max(ps)) %>%
  tidyr::gather(bla, ps, -chr)


don <- loco_res %>%  
  # filter(-log10(p_wald)>1) %>%
  dplyr::filter(trait == "pel_pos") %>%
  dplyr::mutate(chr = ifelse(chr == "23", "X", chr)) %>%
  dplyr::filter(chr != "MtDNA") %>%
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(ps)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(loco_res %>% dplyr::mutate(chr = ifelse(chr == "23", "X", chr)) %>% dplyr::filter(chr != "MtDNA") , ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, ps) %>%
  mutate( BPcum=ps+tot) %>%
  dplyr::mutate(factrait = factor(trait, levels = c("pel_pos","pel_mig"), labels = c("pos-1","mig-6")))

axisdf <- don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

ggplot(don, aes(x=BPcum, y=-log10(p_wald))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "grey10"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0.1, 0.1) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  facet_grid(factrait~., scales = "free_y") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )+
  labs(x = "Genomic Position (Mb)",y = expression(-log[10](italic(p))))

# loco_res %>%
#   # dplyr::filter(trait == mkpl) %>%
#   dplyr::mutate(factrait = factor(trait, levels = c("pel_pos","pel_mig"), labels = c("pos-1","mig-6"))) %>%
#   ggplot()+
#   aes(x = ps/1e6, y = -log10(p_wald))+
#   geom_point(size = 0.5, alpha = 0.5)+
#   facet_grid(factrait~chr, scales = "free", space = "free") +
#   scale_y_continuous(limits = c(0, 10), expand = c(0.01, 0.01))+
#   theme_bw(18)+
#   labs(x = "Genomic Position (Mb)",y = expression(-log[10](italic(p))))

ggsave(filename = "plots/Figure2.png", height = 8, width = 12)
ggsave(filename = "plots/Figure2.pdf", height = 8, width = 12)
