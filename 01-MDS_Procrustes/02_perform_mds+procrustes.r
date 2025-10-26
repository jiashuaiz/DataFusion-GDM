# 02_perform_mds.r
# Purpose: Perform MDS on distance matrices to obtain coordinate spaces
# V1.0.0
# Author: Zhu, Jiashuai (The University of Melbourne; Agriculture Victoria Research)

# MDS----
# Double-center the matrices
double_center <- function(D) {
  n <- nrow(D)
  J <- diag(n) - (1/n) * matrix(1, n, n)
  B <- -0.5 * J %*% (D^2) %*% J
  return(B)
}

# Perform MDS on a pair of distance matrices
perform_mds <- function(A, B) {
  # Double-center the matrices A B
  A_centered <- double_center(A)
  B_centered <- double_center(B)
  
  # Perform eigen-decomposition
  eigen_A <- eigen(A_centered, symmetric = TRUE)
  eigen_B <- eigen(B_centered, symmetric = TRUE)
  
  # Determine optimal dimensionality
  pos_eigen_A <- sum(eigen_A$values > 0)
  pos_eigen_B <- sum(eigen_B$values > 0)
  d_max_pos <- min(pos_eigen_A, pos_eigen_B)
  
  # Calculate proportion of variance explained
  var_explained_A <- cumsum(eigen_A$values[1:d_max_pos]) / sum(eigen_A$values[eigen_A$values > 0])
  var_explained_B <- cumsum(eigen_B$values[1:d_max_pos]) / sum(eigen_B$values[eigen_B$values > 0])
  
  # Find dimensions needed for 100% variance in each matrix
  d_100_A <- min(which(var_explained_A >= 1))
  d_100_B <- min(which(var_explained_B >= 1))

  # Choose optimal dimensionality
  d_opt <- min(d_100_A, d_100_B)
  
  # Calculate MDS coordinates
  X <- eigen_A$vectors[, 1:d_opt] %*% diag(sqrt(eigen_A$values[1:d_opt]))
  Y <- eigen_B$vectors[, 1:d_opt] %*% diag(sqrt(eigen_B$values[1:d_opt]))

  # Name the rows of the coordinate matrices
  rownames(X) <- rownames(A)
  rownames(Y) <- rownames(B)
  
  # Create variance information structure
  var_info <- list(
    eigenvalues_A = eigen_A$values,
    eigenvalues_B = eigen_B$values,
    var_explained_A = var_explained_A,
    var_explained_B = var_explained_B
  )
  
  # Create a single list with all MDS results
  mds_results <- list(
    X = X,
    Y = Y,
    d_opt = d_opt,
    variance_info = var_info
  )
  
  return(mds_results)
}


# procrustes analysis -----
# function to perform Procrustes transformation
apply_procrustes <- function(X_base, Y_base, Y) {
  # Perform Procrustes transformation
  proc <- vegan::procrustes(X_base, Y_base, scale = TRUE, symmetric = FALSE)
  
  # Extract transformation parameters
  translation <- proc$translation
  dilation <- proc$scale
  rotation <- proc$rotation
  
  # Apply transformation to the full Y matrix (Rotate → Scale → Translate)
  Y_transformed <- (Y %*% rotation) * dilation + 
    matrix(translation, nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
  
  return(Y_transformed)
}

# Function to convert MDS coordinates back to distance matrix
coords_to_distances <- function(coords) {
  # Use base R dist() function to compute Euclidean distances
  dist_M <- as.matrix(dist(coords))
  rownames(dist_M) = colnames(dist_M) <- rownames(coords)
  return(dist_M)
}


# debug
if (F) {
  
  # parameters
  np <- 50          # Number of populations
  k <- 25           # Number of shared populations
  output_dir <- "data"  # Output directory
  
  # Load the distance matrices
  M <- readRDS(file.path(output_dir, sprintf("M_np%d_k%d.rds", np, k)))
  A <- M$A
  B <- M$B
  
  # Perform MDS
  cat("Performing MDS analysis...\n")
  mds_results <- perform_mds(A, B)
  
  # Save MDS results
  saveRDS(mds_results, file = file.path(output_dir, sprintf("MDS_np%d_k%d.rds", np, k)))
  cat("MDS analysis completed. Results saved.\n")
  
  # Extract coordinates
  X <- mds_results$X
  Y <- mds_results$Y
  d_opt <- mds_results$d_opt
  
  # Perform Procrustes analysis
  cat("Performing Procrustes analysis...\n")
  
  # Extract common populations
  pop_common <- intersect(rownames(A), rownames(B))
  n_common <- length(pop_common)
  cat("Using", n_common, "common populations for Procrustes analysis\n")
  
  # Extract sub-matrices containing only common populations
  X_sub <- X[pop_common, 1:d_opt]
  Y_sub <- Y[pop_common, 1:d_opt]
  
  Y_transformed <- apply_procrustes(X_sub, Y_sub, Y)
  
  # Map Y_transformed coordinates back to genetic distance space
  B_calibrated <- coords_to_distances(Y_transformed)
  
  # # Extract only the common populations for comparison
  # A_common <- A[pop_common, pop_common]
  # B_common <- B[pop_common, pop_common]
  # B_cal_common <- B_calibrated[pop_common, pop_common]
  
  the_prior <- mean((A - B)^2) #!
  the_posterior <- mean((A - B_calibrated)^2) #!
}