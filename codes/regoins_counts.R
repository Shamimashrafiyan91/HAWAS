library(parallel)

# Define file path and get list of files
file_path <- "/projects/apog/work/models/1MB/Luekemia/ISP_CNN_SEX_Cov_corrected_scale/" # Path to RDS files
files <- list.files(file_path, full.names = TRUE, pattern = "\\.rds$")  # Only pick .rds files

# Set number of cores (leave one free for stability)
num_cores <- min(42, detectCores() - 1)

# Function to process each file
process_files <- function(file) {
  gene_name <- sub("\\.rds$", "", basename(file))  # Extract gene name from filename

  # Try to read file, return NA if fails
  gene <- tryCatch({
    readRDS(file)
  }, error = function(e) return(NA))

  # Check if data is valid before processing
  if (!is.data.frame(gene)) return(setNames(NA, gene_name))

  # Count columns
  num_columns <- ncol(gene) - 1

  # Return named vector (gene_name -> num_columns)
  return(setNames(num_columns, gene_name))
}

# Apply parallel processing
results <- mclapply(files, process_files, mc.cores = num_cores)

# Convert list to a named numeric vector, removing NAs
col_vec <- unlist(results)
col_vec <- col_vec[!is.na(col_vec)]  # Remove NA values

# Save results to file
output_file <- "/projects/apog/work/models/1MB/results/CNN_corrected_scale_regions_counts.txt"
write.table(data.frame(Gene=names(col_vec), Num_Regions=col_vec),
            file = output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Generate histogram only if col_vec is non-empty
if (length(col_vec) > 0) {
  pdf(file = "/projects/apog/work/models/1MB/results/hist_Sig_genes_Binned-CNN_corrected_scale.pdf", width = 5, height = 5)
  hist(col_vec, breaks = 30, main = "Number of regions per Sig Binned-CNN Gene",
       xlab = "Number of regions", col = "black", border = "black")
  dev.off()
} else {
  message("⚠️ Warning: No valid data for histogram!")
}

# Print a few results for verification
print(head(col_vec))









################################################ RF ################################################ 

# library(parallel)
# 
# # Define file path and get list of files
# file_path <- "/projects/apog/work/models/1MB/Luekemia/ISP_RF_SEX_covariate/" # Path to RDS files
# files <- list.files(file_path, full.names = TRUE, pattern = "\\.RDS$")  # Only pick .rds files
# print(length(files))
# # Set number of cores (leave one free for stability)
# num_cores <- min(42, detectCores() - 1)
# 
# # Function to process each file
# process_files <- function(file) {
#   gene_name <- sub("\\.RDS$", "", basename(file))  # Extract gene name from filename
#   # print(gene_name)
#   # Try to read file, return NA if fails
#   gene <- tryCatch({
#     readRDS(file)
#   }, error = function(e) return(NA))
#   
#   if (is.matrix(gene)) {
#     gene <- as.data.frame(gene)
#   }
#   
#   # Check if data is valid before processing
#   if (!is.data.frame(gene)) return(setNames(NA, gene_name))
#   
#   # Count columns
#   print(dim(gene))
#   num_columns <- ncol(gene)
#   
#   # Return named vector (gene_name -> num_columns)
#   return(setNames(num_columns, gene_name))
# }
# 
# # Apply parallel processing
# results <- mclapply(files, process_files, mc.cores = num_cores)
# print(length(results))
# 
# # Convert list to a named numeric vector, removing NAs
# col_vec <- unlist(results)
# print(head(col_vec))
# col_vec <- col_vec[!is.na(col_vec)]  # Remove NA values
# 
# # Save results to file
# output_file <- "/projects/apog/work/models/1MB/results/CRE_RF_regions_counts.txt"
# write.table(data.frame(Gene=names(col_vec), Num_Regions=col_vec),
#             file = output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
# 
# # Generate histogram only if col_vec is non-empty
# if (length(col_vec) > 0) {
#   pdf(file = "/projects/apog/work/models/1MB/results/hist_Sig_genes_CRE_RF_scale.pdf", width = 5, height = 5)
#   hist(col_vec, breaks = 30, main = "Number of regions per Sig CRE_RF Gene",
#        xlab = "Number of regions", col = "black", border = "black")
#   dev.off()
# } else {
#   message("⚠️ Warning: No valid data for histogram!")
# }
# 
# # Print a few results for verification
# print(head(col_vec))
# 
