# HAWAS
A novel concept for conducting histone-acetylome-wide association studies
First: finding significant genes for each method and also real RNA data: run script "Leukemia_diff_expression_IHEC_data.R", you can change the paths on the script and run it, I provided all input data in "data_folder". Mentioned script also plots two scatterplots for fold change (RF vs real) and (CNN vs real) for ~15,000 significant genes from real data.

Second: To compare Leukemia related genes from Disgenet with those that we found in the previous step use script "Disgenet_vs_Leukemia_genes_from_RF_and_CNN.R". Disgenet Leukemia genes are downloaded from Disgenet old portal and its in "data_folder".

____________________________________________________________________________________________________________________________________________

Finding significant regions (related to Leukemia)
Third: we need to find ISP score for each region. To do that use script "ISM_Leukemia_RF_all0_application.R", afterward we have a file for each gene which has the score for all regions (features) in that gene.

Fourth: Doing T-test on all regions from all significant genes. To do that use script ""
