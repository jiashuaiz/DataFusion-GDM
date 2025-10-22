# Dataset Preparation Script - 01-BESMI_Prepare_datasets.r
# This script prepares data for Bootstrapping Evaluation for Structural Missingness Imputation (BESMI)
# V1.0.0
# Author: Zhu, Jiashuai (The University of Melbourne; Agriculture Victoria Research)

# Load necessary libraries
library(tidyverse)

# 1. Input Data Processing ################

# Function to load and prepare the genetic distance matrix
prepare_full_dataset <- function(input_path) {
  # Check if input file exists
  if (!file.exists(input_path)) {
    stop(paste("Input file not found:", input_path))
  }
  
  # Check file extension and load accordingly
  if (grepl("\\.csv$", input_path)) {
    # Load CSV file
    gdm_data <- read.csv(input_path, row.names = 1, check.names = FALSE)
    gdm <- as.matrix(gdm_data)
  } else if (grepl("\\.RData$", input_path)) {
    # Load RData file
    load(input_path)
    
    # Assuming the matrix is loaded as 'gdm'
    if (!exists("gdm")) {
      stop("GDM not found in the loaded data. Please ensure the RData file contains a 'gdm' object.")
    }
  } else {
    stop("Unsupported file format. Please provide a CSV or RData file.")
  }
  
  # Validate that gdm is a matrix
  if (!is.matrix(gdm)) {
    gdm <- as.matrix(gdm)
  }
  
  # Verify the matrix is symmetric
  if (!isSymmetric(gdm)) {
    stop("The genetic distance matrix is not symmetric.")
  }
  
  # Check for missing values
  if (any(is.na(gdm))) {
    stop("The genetic distance matrix already contains missing values.")
  }
  
  # Return the prepared matrix
  return(gdm)
}

# 2. Parameter Setup ################

# Function to determine the number of bootstrap sampling sizes based on k
determine_sampling_sizes <- function(k) {
  if (k <= 3) {
    return(15)  # More samplings for small k (10-20 range)
  } else if (k <= 6) {
    return(8)   # Moderate samplings for medium k (5-10 range)
  } else {
    return(4)   # Fewer samplings for large k (3-5 range)
  }
}

# 3. BESMI Core Process ################

# Function to create masked matrices
create_masked_matrices <- function(full_matrix, k, bs_i) {
  # Get all population names
  pop_names <- rownames(full_matrix)
  n_pops <- length(pop_names)
  
  # Randomly select k populations to mask (group U)
  group_u_indices <- sample(1:n_pops, k)
  group_u <- pop_names[group_u_indices]
  group_s <- pop_names[-group_u_indices]
  
  # Create a copy of the full matrix
  masked_matrix <- full_matrix
  
  # Create mask_position matrix (TRUE = masked, FALSE = not masked)
  mask_position <- matrix(FALSE, nrow = nrow(masked_matrix), ncol = ncol(masked_matrix))
  rownames(mask_position) <- pop_names
  colnames(mask_position) <- pop_names
  
  # Mask only cells where both row and column correspond to populations in group U
  for (i in group_u) {
    for (j in group_u) {
      # Set cell to NA in masked matrix
      masked_matrix[i, j] <- NA
      
      # Set corresponding cell to TRUE in mask_position
      mask_position[i, j] <- TRUE
    }
  }
  
  # Calculate percentage of matrix masked
  total_cells <- nrow(masked_matrix) * ncol(masked_matrix)
  masked_cells <- sum(mask_position)
  masked_percentage <- round(100 * masked_cells / total_cells, 1)
  
  return(list(
    masked_matrix = masked_matrix,
    mask_position = mask_position,
    group_u = group_u,
    group_s = group_s,
    masked_percentage = masked_percentage
  ))
}

# 4. Output Management ################

# Function to ensure directories exist
create_directories <- function() {
  dirs <- c(
    "data",
    "data/training_set",
    "data/position_set", 
    "data/imputation_set",
    "data/summary_tables"
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat(sprintf("Created directory: %s\n", dir))
    }
  }
}

# Function to check and create input directory if necessary
ensure_input_path <- function(input_path) {
  # Extract directory path from the input file path
  dir_path <- dirname(input_path)
  
  # If directory doesn't exist, create it
  if (!dir.exists(dir_path) && dir_path != ".") {
    dir.create(dir_path, recursive = TRUE)
    cat(sprintf("Created directory: %s\n", dir_path))
  }
  
  return(input_path)
}

# Function to process data and generate matrices for each k and bootstrap sampling
process_and_save_data <- function(full_matrix, k_values) {
  create_directories()
  saveRDS(full_matrix, "data/full_dataset.rds")
  
  # Initialize summary table
  summary_table <- data.frame(
    k_value = integer(),
    bootstrap_sample_i = integer(),
    group_u = character(),
    group_s = character(),
    masked_percentage = numeric(),
    training_path = character(),
    mask_position_path = character(),
    stringsAsFactors = FALSE
  )
  
  # Process each k value
  for (k in k_values) {
    # Determine number of bootstrap samples for this k
    n_samples <- determine_sampling_sizes(k)
    
    # Process each bootstrap sample
    for (bs_i in 1:n_samples) {
      # Create masked matrices
      result <- create_masked_matrices(full_matrix, k, bs_i)
      
      # Define file paths
      training_path <- sprintf("data/training_set/masked_k%d_bs%d.rds", k, bs_i)
      mask_position_path <- sprintf("data/position_set/mask_k%d_bs%d.rds", k, bs_i)
      
      # Save matrices
      saveRDS(result$masked_matrix, training_path)
      saveRDS(result$mask_position, mask_position_path)
      
      # Update summary table
      summary_table <- rbind(summary_table, data.frame(
        k_value = k,
        bootstrap_sample_i = bs_i,
        group_u = paste(result$group_u, collapse = ";"),
        group_s = paste(result$group_s, collapse = ";"),
        masked_percentage = result$masked_percentage,
        training_path = training_path,
        mask_position_path = mask_position_path,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Save summary table
  write.csv(summary_table, "data/summary_tables/summary_table_datasets.csv", row.names = FALSE)  
  return(summary_table)
}

# Main execution function
run_besmi_preparation <- function(
    input_path,
    max_k = NULL,
    custom_k_values = NULL
) {
  
  # Check if input file exists
  if (!file.exists(input_path)) {
    cat(sprintf("Input file not found: %s\n", input_path))
    cat("Calling external script to generate the missing matrix...\n")
    
    # Call external script
    system("Rscript ../Simulate_structured_distM/01-Simulate_structured_distM.r")
    
    # Check if the file is now created
    if (!file.exists(input_path)) {
      stop("Error: The input file was not created by the external script.")
    }
  }
  
  # Load and prepare the full dataset
  full_matrix <- prepare_full_dataset(input_path)
  
  # Get total number of populations
  n_pops <- nrow(full_matrix)
  
  # Determine k values to use
  if (!is.null(custom_k_values)) {
    k_values <- custom_k_values
  } else {
    max_k <- ifelse(is.null(max_k), n_pops - 1, min(max_k, n_pops - 1))
    k_values <- 1:max_k
  }
  
  # Process data and get summary table
  summary_table <- process_and_save_data(full_matrix, k_values)
  
  # Return the prepared data and summary information
  return(list(
    full_matrix = full_matrix,
    summary_table = summary_table,
    k_values = k_values
  ))
}


# RUN
if (T) {
  set.seed(233) 
  
  result <- run_besmi_preparation(
    input_path = "./data/GDM_simulated.csv",
    max_k = 5  # maximum of populations to mask
  )
  
  # Display summary information
  cat("\nBESMI Dataset Preparation Complete\n")
  cat("------------------------------------\n")
  cat("Total populations:", nrow(result$full_matrix), "\n")
  cat("k values processed:", paste(result$k_values, collapse = ", "), "\n")
  cat("Total scenarios created:", nrow(result$summary_table), "\n")
  cat("Output directory: ./data/\n")
}
