# 03_procrustes_sensitivity.r
# Purpose: Test sensitivity of Procrustes alignment for different numbers of shared populations
# V1.0.0
# Author: Zhu, Jiashuai (The University of Melbourne; Agriculture Victoria Research)

rm(list = ls())

# Check for required packages
required_packages <- c("vegan", "ggplot2")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing package: %s\n", pkg))
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Choose which data preparation method to use
data_source <- "input"  # Options: "random" or "input"

# Configure parameters based on data source
if (data_source == "random") {
  source("01_prepare_data_random.r")
  np <- 50          # Number of populations for random data
  bias_mu <- 0.1    # Bias for the 2 datasets
  noise_sd <- 0.05  # noise level
  prepare_data_func <- function(k, seed) {
    simulate_distance_matrices(np=np, k=k, 
                               bias_mu=bias_mu, noise_sd=noise_sd)
  }
  file_prefix <- "random"
} else if (data_source == "input") {
  source("01_prepare_data_input.r")
  input_file <- "GDM_simulated.csv"
  bias_mu <- 0.1    # Bias for the 2 datasets
  noise_sd <- 0.05  # Standard deviation for noise
  
  # Read the input file to determine np
  temp_matrix <- read.csv(input_file, header = TRUE, row.names = 1)
  np <- nrow(temp_matrix)
  
  prepare_data_func <- function(k, seed) {
    prepare_matrices_from_csv(input_file, k=k, 
                              bias_mu=bias_mu, noise_sd=noise_sd, 
                              randomize_shared=TRUE, seed=seed)
  }
  file_prefix <- "input"
} else {
  stop("Invalid data_source. Choose 'random' or 'input'")
}

source("02_perform_mds+procrustes.r")

# User-configurable parameters
k_values <- seq(2, min(np, 30), by = 2)  # Ensure k doesn't exceed np
num_tests <- 10  # Number of tests per k
output_dir <- "data"  # Directory to save results
seed_base <- 42   # Base seed for reproducibility

dir.create(output_dir, showWarnings = FALSE)

results <- data.frame()

for (k in k_values) {
  cat(sprintf("Testing sensitivity for k = %d of %d populations\n", k, np))
  
  for (test_id in 1:num_tests) {
    test_seed <- seed_base + test_id * 100 + k
    
    matrices <- prepare_data_func(k, test_seed)
    
    A <- matrices$A
    B <- matrices$B
    
    # Perform MDS
    mds_results <- perform_mds(A, B)
    X <- mds_results$X
    Y <- mds_results$Y
    d_opt <- mds_results$d_opt
    
    # Perform Procrustes transformation
    pop_common <- intersect(rownames(A), rownames(B))
    X_sub <- X[pop_common, 1:d_opt]
    Y_sub <- Y[pop_common, 1:d_opt]
    Y_transformed <- apply_procrustes(X_sub, Y_sub, Y)
    
    # Map Y_transformed coordinates back to genetic distance space
    B_calibrated <- coords_to_distances(Y_transformed)
    
    # Compute the values before and after calibration
    the_prior <- mean((A - B)^2)
    the_posterior <- mean((A - B_calibrated)^2)
    
    improvement <- (the_prior - the_posterior) / the_prior * 100
    
    results <- rbind(results, data.frame(
      k = k,
      test_id = test_id,
      the_prior = the_prior,
      the_posterior = the_posterior,
      improvement = improvement,
      data_source = data_source
    ))
    
    # Save the first test of each k for detailed inspection
    if (F) { # test_id == 1
      detailed_results <- list(
        matrices = matrices,
        mds_results = mds_results,
        the_prior = the_prior,
        the_posterior = the_posterior,
        improvement = improvement
      )
      saveRDS(detailed_results, 
              file.path(output_dir, sprintf("detailed/detailed_%s_k%d.rds", file_prefix, k)))
    }
  }
}

# Save results
results_filename <- sprintf("sensitivity_results_%s.rds", file_prefix)
saveRDS(results, file = file.path(output_dir, results_filename))

# Calculate summary statistics for each k
summary_results <- aggregate(cbind(the_prior, the_posterior, improvement) ~ k, data = results, 
                             FUN = function(x) c(mean = mean(x), sd = sd(x), min = min(x), max = max(x)))
summary_results <- do.call(data.frame, summary_results)

# Save summary
summary_filename <- sprintf("sensitivity_summary_%s.rds", file_prefix)
saveRDS(summary_results, file = file.path(output_dir, summary_filename))

cat("Sensitivity analysis completed. Results saved.\n")