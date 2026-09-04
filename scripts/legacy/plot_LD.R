
library(tidyverse)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))

III_LD <- data.table::fread("focal_variant_ld_chrIII_filtered.ld") %>%
  dplyr::mutate(focal_dist = abs(BP_A - BP_B))

genome_LD <- data.table::fread("focal_variant_ld_genome_wide.ld") %>%
  dplyr::mutate(focal_dist = abs(BP_A - BP_B))

gplt <- ggplot(genome_LD)+
  aes(x = BP_B/1e6, y = R2)+
  geom_point()+
  theme_bw(18)+
  facet_grid(.~CHR_B, scales = "free", space = "free")+
  xlab("Genomic Position")

ggsave(gplt, filename = "sid2LD_genome.png", height = 6, width = 12)

sid2v = III_LD %>%
  dplyr::filter(BP_A==BP_B)

ggplot(III_LD)+
  aes(x = BP_B/1e6, y = R2)+
  geom_point()+
  theme_bw(18)+
  geom_point(color = "red", data=sid2v)+
  facet_grid(.~CHR_B, scales = "free", space = "free")+
  xlim(12,13.780142)+
  xlab("Genomic Position")

ggsave(filename = "sid2LD_zoom.pdf", height = 6, width = 10)
ggsave(filename = "sid2LD_zoom.png", height = 6, width = 10)
