# 01_prepare_data_random.r
# Purpose: Simulate random distance matrices
# V1.0.0
# Author: Zhu, Jiashuai (The University of Melbourne; Agriculture Victoria Research)

# Function to create a random symmetric distance matrix
create_random_symmetric_matrix <- function(n, min_val=0, max_val=1, bias_mu=0.05, noise_sd=0.001) {
  matrix_base <- matrix(runif(n*n, min=min_val, max=max_val), nrow=n, ncol=n)
  random_bias <- matrix(rnorm(n*n, mean=bias_mu, sd=noise_sd), nrow=n, ncol=n)
  matrix_out = matrix_base + random_bias
  matrix_out <- (matrix_out + t(matrix_out)) / 2  # symmetric
  diag(matrix_out) <- 0
  return(matrix_out)
}

# Function to simulate genetic distance matrices with shared base values
simulate_distance_matrices <- function(np, k, bias_mu, noise_sd) {
  if (k < 1 || k > np) {
    stop("k must be between 1 and np")
  }
  
  # Generate population names with leading zeros
  shared_pop_names <- sprintf(paste0("S_%0", 3, "d"), 1:k)
  unique_A_pop_names <- sprintf(paste0("A_%0", 3, "d"), (k+1):np)
  unique_B_pop_names <- sprintf(paste0("B_%0", 3, "d"), (k+1):np)
  
  A_pop_names <- c(shared_pop_names, unique_A_pop_names)
  B_pop_names <- c(shared_pop_names, unique_B_pop_names)
  
  # Create a single large base matrix that we'll use for both A and B
  # This ensures shared populations have the same base values
  full_base_A <- matrix(runif(np*np, min=0, max=1), nrow=np, ncol=np)
  full_base_A <- (full_base_A + t(full_base_A)) / 2  # make symmetric
  diag(full_base_A) <- 0
  
  # Create a separate base matrix for B's unique populations
  full_base_B <- matrix(0, nrow=np, ncol=np)
  
  # Copy the shared part from A to B (upper left k×k quadrant)
  full_base_B[1:k, 1:k] <- full_base_A[1:k, 1:k]
  
  # For the unique part of B, generate new random values
  unique_B_matrix_base <- matrix(runif((np-k)*(np-k), min=0, max=1), nrow=np-k, ncol=np-k)
  unique_B_matrix_base <- (unique_B_matrix_base + t(unique_B_matrix_base)) / 2  # make symmetric
  diag(unique_B_matrix_base) <- 0
  full_base_B[(k+1):np, (k+1):np] <- unique_B_matrix_base
  
  # Generate new random values for the cross-section between shared and unique-B
  B_shared_unique_base <- matrix(runif(k*(np-k), min=0, max=1), nrow=k, ncol=np-k)
  full_base_B[1:k, (k+1):np] <- B_shared_unique_base
  full_base_B[(k+1):np, 1:k] <- t(B_shared_unique_base)
  
  # Add noise to both matrices
  A_noise <- matrix(rnorm(np*np, mean=bias_mu, sd=noise_sd), nrow=np, ncol=np)
  A_noise <- (A_noise + t(A_noise)) / 2  # make symmetric
  diag(A_noise) <- 0
  A <- full_base_A + A_noise
  
  B_noise <- matrix(rnorm(np*np, mean=bias_mu, sd=noise_sd), nrow=np, ncol=np)
  B_noise <- (B_noise + t(B_noise)) / 2  # make symmetric
  diag(B_noise) <- 0
  B <- full_base_B + B_noise
  
  # Set row and column names
  rownames(A) <- colnames(A) <- A_pop_names
  rownames(B) <- colnames(B) <- B_pop_names
  
  # Return matrices along with the np and k values
  return(list(A=A, B=B, np=np, k=k))
}

# debug
if (F) {
  
  library(stats)
  
  set.seed(66) # seed for reproducibility
  
  # User-configurable parameters
  np <- 50          # Number of populations
  k <- 25           # Number of shared populations (1 ≤ k ≤ np)
  bias_mu <- 0.3    # bias
  noise_sd <- 0.05  # Standard deviation for noise
  output_dir <- "data"  # Output directory
  
  cat(sprintf("Generating matrices with np=%d and k=%d\n", np, k))
  
  matrices <- simulate_distance_matrices(np=np, k=k, bias_mu=bias_mu, noise_sd=noise_sd)
  
  # Check if shared populations have the same base values (before noise)
  shared_indices_A <- 1:k
  shared_indices_B <- 1:k
  
  # Extract shared submatrices
  shared_A <- matrices$A[shared_indices_A, shared_indices_A]
  shared_B <- matrices$B[shared_indices_B, shared_indices_B]
  
  # Verify similar base values (they won't be exactly equal due to noise)
  cat("Mean absolute difference between shared parts:", mean(abs(shared_A - shared_B)), "\n")
  cat("Correlation between shared parts:", cor(as.vector(shared_A), as.vector(shared_B)), "\n")
  
  # Visualize
  par(mfrow=c(1,2))
  heatmap(matrices[["A"]], main="Matrix A")
  heatmap(matrices[["B"]], main="Matrix B")
  print(matrices[["np"]])
  print(matrices[["k"]])
}