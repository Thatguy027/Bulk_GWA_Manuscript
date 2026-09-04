library(tidyverse)
library(smplot2)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("data/pos1_rnai_sensitive/processed_strain_fq.rda")

# get delta traits 
timepoint = 2
base_size=18

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

mig6_ph_plot <- pres_plot%>%
  dplyr::filter(Gene == "mig-6") %>%
  ggplot()+
  aes(x = mean_delta )+
  geom_histogram() +
  # scale_fill_manual(values = c(JU1793 = "#F34C00", JU2466 = "#40B4AB")) +
  geom_vline(aes(xintercept = mean_delta), data = pres_plot %>% dplyr::filter(Gene == "mig-6", strain == "JU1793"), color = "#F34C00", size = 1,linetype = "dashed") +
  geom_vline(aes(xintercept = mean_delta), data = pres_plot %>% dplyr::filter(Gene == "mig-6", strain == "JU2466"), color = "#40B4AB", size = 1,linetype = "dashed") +
  theme_bw(base_size)+
  labs(x = "Change in strain frequency", y = "Count")

ggsave(filename = "mig6_JuHighlight.pdf", height = 5, width = 6)

pres_plot%>%
  dplyr::filter(Gene == "pos-1") %>%
  ggplot()+
  aes(x = mean_delta )+
  geom_histogram() +
  geom_vline(aes(xintercept = mean_delta), data = pres_plot %>% dplyr::filter(Gene == "pos-1", strain == "JU1793"), color = "#9794E7", size = 1,linetype = "dashed") +
  geom_vline(aes(xintercept = mean_delta), data = pres_plot %>% dplyr::filter(Gene == "pos-1", strain == "JU2466"), color = "#91D29C", size = 1,linetype = "dashed") +
  theme_bw(base_size)+
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

write.table(delta_t2_traits_addmeans, file = "data/pos1_rnai_sensitive/2026_delta_t2_traits.tsv", quote = F, col.names = T, row.names = F, sep = "\t")

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
  dplyr::filter(trait == "pel_pos") %>%
  arrange(chr, ps) %>%
  mutate( BPcum=ps+tot) %>%
  dplyr::mutate(factrait = factor(trait, levels = c("pel_pos","pel_mig"), labels = c("pos-1","mig-6")))

axisdf <- don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

pos_man_plot <- ggplot(don, aes(x=BPcum, y=-log10(p_wald))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=0.7) +
  scale_color_manual(values = rep(c("grey40", "grey10"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0.1, 0.1) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw(base_size) +
  # facet_grid(factrait~., scales = "free_y") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )+
  labs(y = expression(-log[10](italic(p))))

# mig6-specific qtl
result_file <- "data/cross_experiments/Nov2024_JU_cross_pos_mig/plots/JU1793_JU2466_F2-2_contrast_MIG6g-POS1g_10000_plot_DF.tsv"

result_df <- data.table::fread(result_file) %>%
  dplyr::select(chrom, physical.position, p)

effective.n.tests <- 2000
sig_line <- -log10(0.05 / effective.n.tests)

don_cross <- result_df %>%  
  dplyr::filter(chrom != "MtDNA") %>%
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(physical.position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(result_df, ., by=c("chrom"="chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrom, physical.position) %>%
  mutate( BPcum=physical.position+tot) 

axisdf <- don_cross %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

p_manhattan <- ggplot(don_cross, aes(x=BPcum, y=-log10(p))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chrom)), alpha=0.8, size=0.7) +
  geom_hline(yintercept = sig_line, linewidth = 0.6, linetype = "dashed", color = "firebrick3") +
  scale_color_manual(values = rep(c("grey40", "grey10"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0.1, 0.1) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw(base_size) +
  # facet_grid(factrait~., scales = "free_y") +
  theme(
    legend.position="none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )+
  labs(x = "", y = expression(-log[10](italic(p))))


library(patchwork)

mig6_ph_plot <- mig6_ph_plot + theme(legend.position = "none")
pos_man_plot <- pos_man_plot + theme(legend.position = "none")
p_manhattan <- p_manhattan + theme(legend.position = "none")

# assemble with A (top), B (bottom-left), C (bottom-right)
(pos_man_plot / ((mig6_ph_plot | p_manhattan) + plot_layout(widths = c(1, 2)))) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = base_size))

ggsave(filename = "plots/Figure2.png", height = 8, width = 12)
ggsave(filename = "plots/Figure2.pdf", height = 8, width = 12)
