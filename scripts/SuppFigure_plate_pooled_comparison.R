library(tidyverse)
library(smplot2)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("../data/pos1_original/2024pub_prediction_bootstrap_combined_df.Rdata")


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

# load in plate based pos1 phenotyping and combine with change in frequency plot
plate_phenotypes <- data.table::fread("../data/pos1_original/cendr_pos1_plate_pheno.tsv") %>%
  dplyr::select(strain, plate_pos1=og_trait)

plate_pooled_pos1 <- dplyr::full_join(t2_pos1, plate_phenotypes, by= "strain") %>%
  dplyr::mutate(delta_boot = bootmean-btm_ctrl) %>%
  dplyr::group_by(condition, rep) %>%
  dplyr::arrange(desc(delta_boot)) %>%
  dplyr::mutate(strain_rank = 1:n()) %>%
  na.omit()

# significance testing
# cor.test(plate_pooled_pos1$plate_pos1, plate_pooled_pos1$strain_rank, method = "kendall")
# cor.test(plate_pooled_pos1$delta_boot, plate_pooled_pos1$plate_pos1, method = "spearman")
# 
# clinfun::jonckheere.test(plate_pooled_pos1$delta_boot, as.numeric(plate_pooled_pos1$plate_pos1), nperm = 1000)

ggplot(plate_pooled_pos1, aes(x = factor(plate_pos1), y = delta_boot)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +
  geom_jitter(width = 0.1, alpha = 0.6, color = "darkblue", size = 2) +
  # facet_grid(rep~.)+
  labs(
    x = "Individual Strain Phenotype",
    y = "Pooled ∆Frequency"
  ) +
  theme_bw(26) +
  # Adds the Spearman correlation and p-value to the plot automatically
  ggpubr::stat_cor(aes(group = 1), method = "spearman", label.x.npc = "left", label.y.npc = "top")


ggplot(plate_pooled_pos1, aes(x = factor(plate_pos1), 
                    y = delta_boot, 
                    fill = rep,    # Fills the boxes by replicate
                    color = rep)) + # Colors the points and text by replicate
  
  # 1. The Boxplots: 'dodge' places them side-by-side
  geom_boxplot(
    position = position_dodge(width = 0.8), 
    outlier.shape = NA, # Hides default outliers so we don't plot them twice with jitter
    alpha = 0.4,
    color = "black"     # Keeps the box borders black for a cleaner look
  ) +
  
  # 2. The Points: 'jitterdodge' aligns points to their specific box
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
    alpha = 0.6, size =2
  ) +
  
  # 3. The Statistics: 'group' calculates a separate correlation for each rep
  ggpubr::stat_cor(
    aes(group = rep), 
    method = "spearman", 
    label.x.npc = "left", 
    label.y.npc = "top" 
    # ggpubr will automatically stack the two equations and color code them!
  ) +
  
  # 4. Clean up labels and theme
  labs(
    x = "Individual Strain Phenotype",
    y = "Pooled ∆Frequency",
    fill = "Replicate",
    color = "Replicate"
  ) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02")) +  # Custom colors (optional)
  scale_color_manual(values = c("#1b9e77", "#d95f02")) +
  theme_bw(26) +
  theme(legend.position = "top")
