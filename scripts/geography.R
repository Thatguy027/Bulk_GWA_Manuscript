
library(tidyverse)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))

sid2variant <- data.table::fread("sid2T96K.tsv") 

strain_df <- data.table::fread("20250625_c_elegans_strain_data.csv") %>%
  dplyr::filter(isotype == strain, 
                isotype %in% c(sid2variant$V7)) 

sid2variant <- sid2variant %>%
  dplyr::select(isotype = V7, gt = V8) %>%
  dplyr::mutate(GT = ifelse(grepl("1", gt), "ALT",
                            ifelse(grepl("0", gt), "REF", NA))) %>%
  dplyr::filter(isotype %in% c(strain_df$isotype)) 

isolation_info <- sid2variant %>%
  dplyr::left_join(., strain_df %>% dplyr::select(isotype, lat = latitude, long = longitude, landscape, substrate), by = "isotype")

isolation_info$lat <- as.numeric(isolation_info$lat)
isolation_info$long <- as.numeric(isolation_info$long)

world <- map_data(map = "world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

# world pop
ggplot()+ geom_map(data=world, map=world,
                   aes(x=long, y=lat, map_id=region),
                   color="white", fill="#7f7f7f", size=0.05, alpha=1/4)+
  geom_point(data = isolation_info, aes(x=long, y=lat, fill=GT), shape =21, alpha = 0.7) + 
  scale_fill_manual(values = c("blue", "orange"),name = "Thr96Lys")+
  theme_map()
