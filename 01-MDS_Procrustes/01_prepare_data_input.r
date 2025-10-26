# 01_prepare_data_input.r
# Purpose: Prepare distance matrices from input CSV file with controlled shared information
# V1.0.0
# Author: Zhu, Jiashuai (The University of Melbourne; Agriculture Victoria Research)

# Function to prepare matrices from an input CSV file
prepare_matrices_from_csv <- function(input_file, k, bias_mu = 0.1, noise_sd = 0.05, randomize_shared = TRUE, seed = 42) {
  full_matrix <- read.csv(input_file, header = TRUE, row.names = 1)
  full_matrix <- as.matrix(full_matrix)
  
  if (nrow(full_matrix) != ncol(full_matrix)) {
    stop("Input matrix must be square")
  }
  
  # Get all population names
  all_pop_names <- rownames(full_matrix)
  np <- length(all_pop_names)
  
  # Check if k is valid
  if (k < 1 || k > np) {
    stop("k must be between 1 and np (number of populations in input matrix)")
  }
  
  set.seed(seed)
  
  # Select k populations to be shared
  if (randomize_shared) {
    # Randomly select shared populations
    shared_pop_names <- sample(all_pop_names, k)
  } else {
    # Use the first k populations as shared
    shared_pop_names <- all_pop_names[1:k]
  }
  
  # Remaining populations are unique
  unique_pop_names <- setdiff(all_pop_names, shared_pop_names)
  
  matrix_A <- full_matrix
  matrix_B <- full_matrix
  
  # Add noise
  noise_A <- matrix(rnorm(nrow(matrix_A) * ncol(matrix_A), mean = 0, sd = noise_sd), nrow = nrow(matrix_A))
  matrix_A <- matrix_A + noise_A
  matrix_A[matrix_A < 0] <- abs(matrix_A[matrix_A < 0])
  
  noise_B <- matrix(rnorm(nrow(matrix_B) * ncol(matrix_B), mean = bias_mu, sd = noise_sd), nrow = nrow(matrix_B))
  
  matrix_B <- matrix_B + noise_B
  matrix_B[matrix_B < 0] <- abs(matrix_B[matrix_B < 0])
  
  # Ensure zeros on diagonal
  diag(matrix_A) <- 0
  diag(matrix_B) <- 0
  
  # Create new labels for each matrix
  label_A <- rep("", np)
  label_B <- rep("", np)
  names(label_A) <- all_pop_names
  names(label_B) <- all_pop_names
  
  # Set labels for shared populations
  for (pop in shared_pop_names) {
    label_A[pop] <- paste0("S_", pop)
    label_B[pop] <- paste0("S_", pop)
  }
  
  # Set labels for unique populations
  for (pop in unique_pop_names) {
    label_A[pop] <- paste0("A_", pop)
    label_B[pop] <- paste0("B_", pop)
  }
  
  # Update row and column names
  rownames(matrix_A) <- label_A[rownames(matrix_A)]
  colnames(matrix_A) <- label_A[colnames(matrix_A)]
  
  rownames(matrix_B) <- label_B[rownames(matrix_B)]
  colnames(matrix_B) <- label_B[colnames(matrix_B)]
  
  # debugging information
  cat("Total populations:", np, "\n")
  cat("Shared populations:", k, "\n")
  cat("Unique populations:", length(unique_pop_names), "\n")
  
  # Return the matrices and parameters
  return(list(A = matrix_A, B = matrix_B, np = np, k = k))
}

# Debug
if (FALSE) {
  
  # User-configurable parameters
  input_file <- "GDM_simulated.csv"
  k <- 25           # Number of shared populations
  noise_sd <- 0.05  # Standard deviation for noise
  output_dir <- "data"  # Output directory
  
  cat(sprintf("Preparing matrices from input file with k=%d\n", k))
  
  matrices <- prepare_matrices_from_csv(input_file, k=k, noise_sd=noise_sd)
  heatmap(matrices[["A"]])
  heatmap(matrices[["B"]])
  print(matrices[["np"]])
  print(matrices[["k"]])
}
