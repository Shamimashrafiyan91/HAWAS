# HAWAS: A Novel Framework for Histone-Acetylome-Wide Association Studies  

The HAWAS-gene approach leverages machine learning (ML) models to predict gene expression from epigenetic signals, enabling the identification of genes and genomic regions associated with disease. It systematically compares predicted gene expression profiles between control and disease samples to pinpoint significant epigenetic regulatory changes.

HAWAS consists of two main components:

- **HAWAS-gene test** – Identifies disease-associated genes based on predicted gene expression.
- **HAWAS-region test** – Identifies significant regulatory regions within disease-associated genes.

## HAWAS-gene Test

Histone modifications, particularly H3K27ac, play a crucial role in gene regulation and are closely associated with active transcription. The HAWAS-gene approach leverages machine learning (ML)-based predictions of gene expression to identify genes with significant alterations between healthy and disease states. By comparing the predicted expression profiles of thousands of genes, this approach enables the discovery of potential disease-associated genes driven by epigenetic regulation.

### Input:
- **Gene Models**: `z` gene models.
- **Gene Matrix**: For each gene `i`, the input is a matrix `G_{n,m}` representing H3K27ac signal data:
   - `n` rows: Samples.
   - `m` columns: Regulatory regions.


### Algorithm:

1. **Predict Gene Expression Using Pre-trained Models**:
    Initialize matrix `E^{n × z}` for storing predicted expression values.
    For each gene `i` in `{1,2,...,z}`:
      - Load pre-trained model `M_i`.
      - Load gene matrix `G^{n × m}` for gene `i`.
      - For each sample `s` in `{1,2,...,n}` in `G_{s,m}`:
        - Predict expression of `G_{s,*}` using model `M_i` and store in `E_{s,i}`.

2. **Identify Differentially Expressed Genes Using DESeq2**:
   Define a design matrix to distinguish control and disease samples.
   Apply DESeq2 on `E_{n,z}` matrix to identify disease-associated genes.

### Output:
- A list of disease-associated genes


#### References
- Love, M.I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8
---


## HAWAS-region Test

While identifying HAWAS genes provides critical insights, understanding which regulatory elements contribute to these expression changes is essential for deciphering disease mechanisms. The HAWAS-region test extends the analysis to the level of gene regulatory regions by quantifying the contribution of individual regions (features) to gene expression using **ISP (Importance Score for Prediction)**. This method identifies regulatory elements whose influence on gene expression differs between healthy and disease conditions.

### Input
- **Gene Models**: `z` gene models.
- **Gene Matrix**: For each gene `i`, the input is a matrix `G_{n,m}` representing H3K27ac signal data:
  - `n` rows: Samples.
  - `m` columns: Regulatory regions.

### Algorithm

1. **Compute ISP Score**:
   - For each gene `i`:
     - For each sample `s` in `{1, 2, ..., n}`:
       - For each region `r` in `{1, 2, ..., m}`:
         - Compute the ISP score for sample `s` and region `r` using the ISP formula (refer to Eq.~\ref{eq:ISP}).
         - Store the result in matrix `I_{s,r}`.
   - **Note**: For simplicity, it is assumed that each gene `i` has the same number of regions `m`.

2. **Perform t-test**:
   - For each gene `i`:
     - For each region `r` in `{1, 2, ..., m}`:
       - Perform a two-sided t-test between the control and disease groups using the ISP scores from `I_{*,r}`.
       - Store the resulting p-value in `L_{i,m}`.
   - **Note**: A t-test is chosen because it is well-suited for comparing continuous values (e.g., region importance scores).

3. **Apply FDR Correction**:
   - After processing all `z` genes, there are `z x m` p-values stored in list `L`.
   - Apply False Discovery Rate (FDR) correction to adjust for multiple testing ([Benjamini & Hochberg, 1995](#)).
   - The output is a list `W` containing significant regions with their adjusted p-values.

### Output
- **List `W`**: Contains significant regulatory regions with their adjusted p-values.

#### References
- Benjamini, Y., & Hochberg, Y. (1995). Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. *Journal of the Royal Statistical Society: Series B (Methodological)*, 57(1), 289–300.

---
# Case study: Chronic Lymphocytic Leukemia (CLL) data

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
      - `chr15:29310344-29310612    3`
       

Ensure that the file adheres to this format for downstream analyses.



---

## Input Data  
All necessary input files are provided in the **`data_folder`** directory. Ensure that paths in the scripts are updated to correctly point to this folder before running the analysis.  



