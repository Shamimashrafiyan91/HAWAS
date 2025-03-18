#converting output of sig regions to a format that is acceptable for PASTAA running
library(dplyr)
library(tidyr)
library(stringr)
sig_regions <- read.table("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/sig_unique_regions_after_FDR_smaller_0.05_Luekemia_Corrected_scale_CNN_gens_byDSEQ_COV_SEX_no_filtering_before_FDR.txt", header = TRUE)
sig_regions <- sig_regions %>%mutate (rank = row_number())
sig_regions_up <- sig_regions %>%select(names, rank)

sig_regions_up_tm <- sig_regions_up %>%
separate(names, into = c("chr","start","end"), sep= "\\.") %>%
mutate(new_format = paste0(chr, ":", start, "-", end)
,rank = row_number())%>%
select(new_format, rank)

write.table(sig_regions_up_tm, "/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/sig_region_CNN_after_scale_correction_input_PASTAA.txt",row.names = FALSE, col.names = FALSE, quote = FALSE)
