# Load necessary libraries
library(foreach)
library(iterators)
library(doParallel)
library(data.table)

# Set input and output paths
input_path <- "/projects/apog/work/IHEC/ValidateInteractions/CNN_bestwarm_os_ism_all_leukemia_expand_len10"
output_path <- "/projects/apog/work/models/1MB/Luekemia/ISP_CNN_SEX_Cov_corrected_scale/"

# Load the sample and gene information
desired_samples <- readLines("/projects/apog/work/models/1MB/results/samples_leukemia_adult_removed_outliers.txt")
sig_genes_CNN <- readLines("/projects/apog/work/models/1MB/results/sig_genes_CNN_corrected_scaled.txt")

# sig_genes_CNN <- sig_genes_CNN[1:10]  # Limit to 10 genes for testing

# Function to process each gene file for a given sample
process_gene_file <- function(gene_file, sample_name) {
  if (!file.exists(gene_file)) {
    warning(paste("File does not exist:", gene_file))
    return(NULL)
  }
  
  # Read the file and check for errors
  tryCatch({
    gene_data <- fread(cmd = paste("zcat", gene_file), skip = 2, col.names = c("chr", "start", "end", "predicted", "ISM"))
    
    if (nrow(gene_data) == 0) {
      warning(paste("No valid data in file:", gene_file))
      return(NULL)
    }
    
    # Create unique region names
    gene_data[, region := paste0("chr", chr, ".", start, ".", end)]
    
    # Compute the formula
    gene_data[, (sample_name) := -1 * (log2((ISM + 1) / (predicted + 1)))]
    
    # Return only the region names and computed values
    return(gene_data[, .(region, get(sample_name))])
    
  }, error = function(e) {
    warning(paste("Error reading file:", gene_file, "Error:", conditionMessage(e)))
    return(NULL)
  })
}

# Set up parallel backend
num_cores <- 58
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Ensure function is properly exported
clusterExport(cl, varlist = c("process_gene_file", "input_path", "desired_samples"), envir = environment())

# Process genes in parallel
results <- foreach(gene = sig_genes_CNN, .packages = c("data.table")) %dopar% {
  sample_values_list <- list()
  
  for (sample_name in desired_samples) {
    sample_folder <- file.path(input_path, sample_name)
    gene_file <- file.path(sample_folder, paste0(gene, ".txt.gz"))
    
    if (file.exists(gene_file)) {
      gene_data <- process_gene_file(gene_file, sample_name)
      
      if (!is.null(gene_data) && nrow(gene_data) > 0) {
        # Rename the second column to include the sample name
        # setnames(gene_data, old = names(gene_data)[2], new = paste0("sample_", sample_name))
        setnames(gene_data, old = names(gene_data)[2], new = paste0(sample_name))
        sample_values_list[[sample_name]] <- gene_data
      }
    }
  }
  
  # Combine data from all samples for this gene
  if (length(sample_values_list) > 0) {
    # Merge while ensuring unique column names
    combined_df <- Reduce(function(x, y) merge(x, y, by = "region", all = TRUE), sample_values_list)
    
    row_names = combined_df$region
    combined_df = combined_df[,-1]
    rownames(combined_df)= row_names
    combined_df = t(combined_df)
    colnames(combined_df) = row_names
    combined_df <- as.data.frame(combined_df)
    
    combined_df$region <- rownames(combined_df)
    rownames(combined_df) <- NULL
    combined_df <- combined_df[, c("region", setdiff(names(combined_df), "region"))]
    # Save the dataframe
    write.table(combined_df, file = file.path(output_path, paste0(gene, ".txt")), sep = "\t", quote = FALSE, col.names = NA)
    saveRDS(combined_df, file = file.path(output_path, paste0(gene, ".rds")))
  }
  
}

# Stop the cluster
stopCluster(cl)
print("Done!")