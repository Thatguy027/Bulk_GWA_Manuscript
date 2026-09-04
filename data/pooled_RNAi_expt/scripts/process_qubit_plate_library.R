library(tidyverse)
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

ratios <- data.table::fread("data/Library_RNAi_Qubit_SZLC.csv", header = T) %>%
  tidyr::gather(col, ratio, -Well) %>%
  dplyr::rename(row = Well) %>%
  dplyr::mutate(well = paste0(row, col))

standards_wells <- c(paste0("F",7:12),
               paste0("G",7:12),
               paste0("H", 7:12))

standard_conc <- rep(c(0,1,2,4,8,10),3) 

standard_df <- data.frame(conc=standard_conc, well=standards_wells) %>%
  dplyr::left_join(., ratios, by = "well")

standard <- lm(conc ~ 0+ratio, data = standard_df)

samples <- data.table::fread("data/Library_RNAi_Qubit_SZLC.csv", header = T) %>%
  dplyr::rename(row = Well) %>%
  tidyr::gather(col, ratio, -row) 

conc = predict( standard, newdata = samples)*10

samples_conc <- data.frame(samples, conc) %>% 
  dplyr::select(-ratio) %>%
  tidyr::spread(col, conc) %>%
  dplyr::select(row, "1", "2","3","4","5","6","7","8","9","10","11","12")


write.table(samples_conc, "data/Library_RNAi_Qubit_SZLC_conc.csv", col.names = T, row.names = F, quote = F, sep = ",")
