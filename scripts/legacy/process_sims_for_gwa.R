library(tidyverse)


setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))

for(trait_file in list.files("../traits_with_validated_qtl/")){
  temp_ori_trait <- data.table::fread(glue::glue("../traits_with_validated_qtl/","{trait_file}"))
  
  if(!exists("original_traits")){
    original_traits <- temp_ori_trait
  } else {
    original_traits <- dplyr::full_join(original_traits,temp_ori_trait, by = "strain")
  }
  
}

sim_freq <- readRDS("plots/nnls_estFreqs.RDS")


for(trait in names(sim_freq)){
  
  temp_df <- t(sim_freq[[trait]])
  df_cols <- paste(trait,colnames(temp_df),sep = "_")
  strains <- sapply(str_split(row.names(temp_df), pattern = "_"), function(x){x[1]})
  temp_df <- data.frame(strain = strains, temp_df)
  colnames(temp_df) <- c("strain",df_cols)
  
  phenotyped_strains <- dplyr::select(original_traits, strain, `trait`) %>%
    na.omit() %>%
    dplyr::pull(strain)
  
  temp_df <- dplyr::filter(temp_df, strain %in% phenotyped_strains)
  
  if(!exists("trait_df")){
    trait_df <- temp_df
  } else {
    trait_df <- dplyr::full_join(trait_df, temp_df, by = "strain")
  }
  
  
}

write.table(trait_df, file = "processed_simFreq_traits.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

#NXF_VER=20.04.0 nextflow run andersenlab/nemascan -profile mappings --vcf 20210121 --traitfile processed_simFreq_traits.tsv


