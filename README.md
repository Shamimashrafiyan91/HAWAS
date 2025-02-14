# HAWAS: A Novel Framework for Histone-Acetylome-Wide Association Studies  

The HAWAS-gene approach leverages ML-based predictions of gene expression to identify genes with significant alterations between healthy and disease states. By comparing the predicted expression profiles of thousands of genes, this approach enables the discovery of potential disease-associated genes driven by epigenetic regulation.

## HAWAS-gene Test

In the first step, pre-trained ML models, such as CRE-RF and Binned-CNN, are used to predict gene expression counts based on H3K27ac signal data from both control and disease samples. For \( z \) genes, there are \( z \) corresponding models, where each model \( M_i \) (\( i \in \{1, 2, \dots, z\} \)) is associated with an input matrix \( G_{n,m} \). Here, \( n \) represents the total number of samples, including both control and disease groups, and \( m \) represents the number of genomic features (regions) of the gene model \( M_i \).

For each gene \( i \), the model produces the predicted expression count and stores it in \( E_{s,i} \). This process is repeated for all \( z \) genes, and the final results are stored in the output matrix \( E_{n,z} \), where \( n \) (rows) corresponds to the input samples, and \( z \) (columns) represents the predicted expression values for all genes.

In the second step, differentially expressed genes between control and disease samples are identified using DESeq2. Specifically, a design matrix is constructed to model the two conditions: control and disease. DESeq2 is then applied to the predicted expression count matrix \( E_{n,z} \). The test returns a list of disease-associated genes based on statistically significant expression changes.

### Algorithm: HAWAS-gene Test


Input: H3K27ac signal data for control and disease samples, z gene models
Output: A list of disease-associated genes

Step 1: Predict Gene Expression Using Pre-trained Models
1. Initialize matrix E^{n × z} for storing predicted expression values
2. For each gene i in {1,2,...,z}:
   a. Load pre-trained model M_i
   b. Load gene matrix G^{n × m} for gene i
   c. For each sample s in {1,2,...,n} in G_{s,m}:
      i. E_{s,i} = Predict expression of G_{s,*} with model M_i

Step 2: Identify Differentially Expressed Genes Using DESeq2
3. Define design matrix to distinguish control and disease samples
4. Apply DESeq2 on E_{n,z} matrix to identify disease-associated genes
5. Output: A list of disease-associated genes

---------

## HAWAS-gene Test

In the first step, pre-trained ML models, such as CRE-RF and Binned-CNN, are used to predict gene expression counts based on H3K27ac signal data from both control and disease samples. For \( z \) genes, there are \( z \) corresponding models, where each model \( M_i \) (\( i \in \{1, 2, \dots, z\} \)) is associated with an input matrix \( G_{n,m} \). Here, \( n \) represents the total number of samples, including both control and disease groups, and \( m \) represents the number of genomic features (regions) of the gene model \( M_i \).

For each gene \( i \), the model produces the predicted expression count and stores it in \( E_{s,i} \). This process is repeated for all \( z \) genes, and the final results are stored in the output matrix \( E_{n,z} \), where \( n \) (rows) corresponds to the input samples, and \( z \) (columns) represents the predicted expression values for all genes.

In the second step, differentially expressed genes between control and disease samples are identified using DESeq2. Specifically, a design matrix is constructed to model the two conditions: control and disease. DESeq2 is then applied to the predicted expression count matrix \( E_{n,z} \). The test returns a list of disease-associated genes based on statistically significant expression changes.

### Algorithm: HAWAS-gene Test


Input: H3K27ac signal data for control and disease samples, z gene models
Output: A list of disease-associated genes

Step 1: Predict Gene Expression Using Pre-trained Models
1. Initialize matrix E^{n × z} for storing predicted expression values
2. For each gene i in {1,2,...,z}:
   a. Load pre-trained model M_i
   b. Load gene matrix G^{n × m} for gene i
   c. For each sample s in {1,2,...,n} in G_{s,m}:
      i. E_{s,i} = Predict expression of G_{s,*} with model M_i

Step 2: Identify Differentially Expressed Genes Using DESeq2
3. Define design matrix to distinguish control and disease samples
4. Apply DESeq2 on E_{n,z} matrix to identify disease-associated genes
5. Output: A list of disease-associated genes


---
# Running the code
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

---
## Step 4: Identify transcription factors (TFs) associated with the significant regions  
Before starting I suggest you to visit the PASTAA website and learn running it in detail, I also tried to elaborat everything here
   - We used PASTAA tool, it has three scripts:  PSCM_to_PSEM.cpp, TRAP.cpp, and  PASTAA.cpp.
   - We used STARE for first two steps to get the affinity file.
   - For the second step we need to provid DNA sequence for significant regions.
      - use **`Finding_DNA_sed_for_sig_regions.R`** to find the DNA sequence for significant regions, dont forget to provide the path to save the DNA seq fasta file.
   - last part needs significant regions and affinity file from previous step.
     - affinity matrix is named **`Leukemia_sig_regions_TRAP_Affinity.txt`** and you can find it in "data_folder".
     - The significant regions file should be formatted as shown below. Each line represents a genomic region and its corresponding rank or significance value:

    
   - **Format Explanation**:
      - Each line consists of:
        - The **chromosomal region** in the format `chr<chromosome>:<start>-<end>`.
        - A **rank or significance value** (integer) separated by a tab.

   - **Example**:
      - `chr1:203817821-203818036    1`
      - `chr3:10273780-10274017    2`
      - `chr15:29310344-29310612    2`
       

Ensure that the file adheres to this format for downstream analyses.



---

## Input Data  
All necessary input files are provided in the **`data_folder`** directory. Ensure that paths in the scripts are updated to correctly point to this folder before running the analysis.  



