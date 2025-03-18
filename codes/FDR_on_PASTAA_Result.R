# finding adjusted p-value for TFs that I giot from PASTAA:

CNN_PASTAA <- read.table("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/Pastaa_output_STARE_sorted_CNN_corrected_scale.txt")
colnames(CNN_PASTAA) = c("TFs", "p-value", "V3", "V4","V5","V6","V7")
CNN_PASTAA$adj_p_value <- p.adjust(CNN_PASTAA$`p-value`, method = "BH")

CNN_PASTAA_up = CNN_PASTAA[CNN_PASTAA$adj_p_value <= 0.01,]
write.table(CNN_PASTAA_up, "/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/PASTTA_CNN_scaled_with_adj_p-value.txt",col.names = TRUE)
