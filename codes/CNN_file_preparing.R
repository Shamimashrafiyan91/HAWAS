# #finding ISP score for CNN genes
# 
# input_path <- "/projects/apog/work/IHEC/ValidateInteractions/CNN_bestwarm_os_ism_all_leukemia_expand_len10"
# output_path <- "/projects/apog/work/models/1MB/Luekemia/ISP_CNN_SEX_Cov"
# desired_samples <-readLines("/projects/apog/work/models/1MB/results/samples_leukemia_adult_removed_outliers.txt") 
# sig_genes_CNN <- readLines("/projects/apog/work/models/1MB/results/sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_CNN.txt"
# )

# Load necessary libraries
library(foreach)
library(iterators)
library(doParallel)
library(data.table)

# Set input and output paths
input_path <- "/projects/apog/work/IHEC/ValidateInteractions/CNN_bestwarm_os_ism_all_leukemia_expand_len10"
output_path <- "/projects/apog/work/models/1MB/Luekemia/ISP_CNN_SEX_Cov/"

# Load the sample and gene information
desired_samples <- readLines("/projects/apog/work/models/1MB/results/samples_leukemia_adult_removed_outliers.txt")
print(head(desired_samples))
sig_genes_CNN <- readLines("/projects/apog/work/models/1MB/results/sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_CNN.txt")
sig_genes_CNN <- sig_genes_CNN[1:3]
# Function to process each gene file for a given sample
# process_gene_file <- function(gene_file) {
#   # Read the gene file, skip the first two comment lines, and load the data using data.table for speed
#   gene_data <- fread(gzfile(gene_file), skip = 2, col.names = c("chr", "start", "end", "predicted", "ISM"))
#   
#   # Create region names from the first three columns
#   gene_data[, region := paste0("chr", chr, ".", start, ".", end)]
#   
#   # Compute the formula for each region
#   gene_data[, value := -1 * (log2((ISM + 1) / (predicted + 1)))]
#   
#   # Return only region names and the computed values
#   return(gene_data[, .(region, value)])
# }


process_gene_file <- function(gene_file) {
  print(paste("Processing file:", gene_file))  # Add this to debug
  # Check if file exists
  if (!file.exists(gene_file)) {
    stop(paste("File does not exist:", gene_file))
  }
  
  # Try reading the file
  gene_data <- fread(gzfile(gene_file), skip = 2, col.names = c("chr", "start", "end", "predicted", "ISM"))
  
  # Create region names from the first three columns
  gene_data[, region := paste0("chr", chr, ".", start, ".", end)]
  
  # Compute the formula for each region
  gene_data[, value := -1 * (log2((ISM + 1) / (predicted + 1)))]
  
  # Return only region names and the computed values
  return(gene_data[, .(region, value)])
}


# Set up parallel backend to use the number of cores available
num_cores <- 58 # leave 1 core for the system
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Loop through significant genes in parallel
foreach(gene = sig_genes_CNN, .packages = c("data.table"), .export = c("process_gene_file")) %dopar% {
  # Initialize an empty list to store data for each sample
  sample_values_list <- list()
  
  # Loop through the samples
  for (sample_name in desired_samples) {
     print("________________________________")
     print(sample_name)
     # sample_folder <- file.path(input_path, sample_name)
     sample_folder <- paste0(input_path, "/", sample_name)
    print("________________________________")
    print(sample_folder)
    # gene_file <- file.path(sample_folder, paste0(gene, ".txt.gz"))
    gene_file <- paste0(sample_folder, "/", paste0(gene))
    
    # Check if the gene file exists for the sample
    if (file.exists(gene_file)) {
      # Process the gene file and store the resulting data
      gene_data <- process_gene_file(gene_file)
      
      # Store the gene data in the list, keyed by sample name
      sample_values_list[[sample_name]] <- gene_data
    }
  }
  
  # Combine data from all samples for this gene
  if (length(sample_values_list) > 0) {
    # Merge all sample data based on regions
    combined_df <- Reduce(function(x, y) merge(x, y, by = "region", all = TRUE), sample_values_list)
    
    # Set the sample names as rownames
    rownames(combined_df) <- names(sample_values_list)
    
    # Remove the 'region' column and transpose the data frame (regions become columns)
    combined_df <- as.data.table(t(combined_df[, -1, with = FALSE]))
    colnames(combined_df) <- combined_df$region  # Set regions as column names
    
    # Save the dataframe as a text file and RDS file
    write.table(combined_df, file = file.path(output_path, paste0(gene, ".txt")), sep = "\t", quote = FALSE, col.names = NA)
    saveRDS(combined_df, file = file.path(output_path, paste0(gene, ".rds")))
  }
}

# Stop the cluster
stopCluster(cl)
print("Done!")
