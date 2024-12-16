# HAWAS: A Novel Framework for Histone-Acetylome-Wide Association Studies  

This repository provides scripts and data for conducting a Histone-Acetylome-Wide Association Study (HAWAS) to identify significant genes and regions related to leukemia. Below is an overview of the analysis pipeline and instructions for using the provided scripts.

---

## Step 1: Identifying Significant Genes  
To identify significant genes for leukemia using multiple methods, including Random Forest (RF), Convolutional Neural Networks (CNN), and real RNA-seq data:  

1. Run the script **`Leukemia_diff_expression_IHEC_data.R`**.  
   - Modify the paths in the script to match your local setup.  
   - Input data is available in the provided **`data_folder`**.  

2. This script uses ~15,000 significant genes from real RNA-seq data and generate two scatter plots:  
   - Fold change comparison between **RF** and **real RNA data**.  
   - Fold change comparison between **CNN** and **real RNA data**.  

---

## Step 2: Comparing Leukemia Genes with DisGeNET  
To compare leukemia-related genes identified in Step 1 with known leukemia genes from **DisGeNET**:  

1. Run the script **`Disgenet_vs_Leukemia_genes_from_RF_and_CNN.R`**.  
   - Input data for DisGeNET leukemia genes is located in **`data_folder`** (downloaded from the old DisGeNET portal).  

2. This script will produce a dot plot to evaluate overlaps and unique findings.  

---

## Step 3: Finding Significant Regions Related to Leukemia  
To identify significant regions (features) within each significant gene:  

1. **Calculate ISP Scores**  
   - Use the script **`ISM_Leukemia_RF_all0_application.R`**.  
   - This script calculates ISP scores for each region in all significant genes.  
   - The output includes a file for each gene, containing the ISP scores for all its regions.  

2. **Perform T-tests on All Regions from All Genes**  
   - Conduct T-tests on regions from all significant genes using the script **`T_Test_on_ISP_to_find_sig_regions.R`**.  
   - This step performs a T-test on all regions and saves all result regardles of their p-value.  
   - make sure to provide a path to save the result for "T_test_results" parameter.

3. **Perform FDR on result of previous step**

   - Perform FDR  using **`FDR_on_all_regions.R`** on result of previous step and filter the region based on adjusted p-value = 0.05
   - If there is redundant regions, pick the one with smallest p-value

## Step 4: Identify transcription factors (TFs) associated with the significant regions  
   - 

---

## Input Data  
All necessary input files are provided in the **`data_folder`** directory. Ensure that paths in the scripts are updated to correctly point to this folder before running the analysis.  

---

## Output Overview  
- **Step 1**: List of ~15,000 significant genes and scatter plots comparing fold changes.  
- **Step 2**: Analysis of overlaps and unique genes between the identified significant genes and DisGeNET leukemia genes.  
- **Step 3**: Files with ISP scores for regions and results of T-tests identifying significant regions.  

---

This pipeline serves as a comprehensive approach for studying leukemia-related genes and regions using HAWAS. Feel free to explore the scripts and adapt them as needed for your research.
