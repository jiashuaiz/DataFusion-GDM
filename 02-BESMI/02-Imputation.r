# BESMI Framework - 02-Imputation.r
# - For a single dataset in BESMI framework
# - Tails-Chain Approach
# V1.0.0
# Author: Zhu, Jiashuai (The University of Melbourne; Agriculture Victoria Research)

# Load necessary libraries
library(tidyverse)
library(reshape2)

# Matrix Symmetry Function
make_symmetric <- function(M) {
  M_sym <- (M + t(M)) / 2
  diag(M_sym) <- 0
  return(M_sym)
}

# Evaluation Metrics
calculate_distance <- function(a, b, method = "mae") {
  if (length(a) != length(b)) {
    stop("Vectors must have the same length")
  }
  
  # Remove NA values
  valid <- !is.na(a) & !is.na(b)
  if (sum(valid) < 3) {  # Need at least 3 points for meaningful calculation
    return(NA)
  }
  
  if (method == "ssd") {
    # Sum of Squared Differences
    return(sum((a - b)^2, na.rm = TRUE))
  } else if (method == "mae") {
    # Mean Absolute Error
    return(mean(abs(a - b), na.rm = TRUE))
  } else if (method == "correlation") {
    # Pearson correlation - note: higher is better, unlike other metrics
    return(cor(a, b, method = "pearson", use = "pairwise.complete.obs"))
  } else if (method == "rmse") {
    # Root Mean Squared Error
    return(sqrt(mean((a - b)^2, na.rm = TRUE)))
  } else {
    stop("Unknown distance method. Use 'ssd', 'mae', 'correlation', or 'rmse'")
  }
}

# Initialized Imputation (iteration = 0)
initialize_M <- function(M) {
  df <- as.data.frame(M)
  
  # mean imputation (Baseline)
  for (col in colnames(df)) {
    df[[col]][is.na(df[[col]])] <- mean(df[[col]], na.rm = TRUE)
  }
  
  M_init <- as.matrix(df)
  
  M_init <- make_symmetric(M_init)
  
  return(M_init)
}

# Tails-Chain Iterative Imputation Process ---------------------------------------
library(mice)
iterative_imputation <- function(M_input,
                                 M_mask,
                                 M_real = NULL, 
                                 method = "lasso.norm", 
                                 max_iterations = 5, 
                                 imputation_convergence_threshold = 1e-3,
                                 propagation_convergence_threshold = 1e-3,
                                 distance_metric = "mae",
                                 k = NA, bs_i = NA) {
  
  # Initialize
  iteration <- 1
  converged <- FALSE
  
  all_metrics <- data.frame(
    k = integer(),
    bs = integer(),
    iteration = integer(),
    imputation_dis = numeric(),
    propagation_dis = numeric(),
    runtime = numeric(),
    improvement_pct = numeric(),
    converged = logical(),
    averaged = logical(),
    stringsAsFactors = FALSE
  )
  
  tails_chain <- list()
  
  M_initialized <- initialize_M(M_input)
  tails_chain[[1]] <- M_initialized
  
  cat(sprintf("Starting imputation (k=%s, bs=%s, metric=%s)\n", k, bs_i, distance_metric))
  
  # iteration
  while (!converged && iteration <= max_iterations) {
    # Record iteration start time
    start_time_iter <- Sys.time()
    
    cat(".")
    
    # Convert matrix to data frame for MICE
    input_df <- as.data.frame(M_input)
    
    # Apply MICE imputation
    # Note: Each run is independent with progressively longer internal iterations
    mice_result <- tryCatch({
      mice(input_df, m = 1,
           maxit = iteration,  # Progressive internal iterations 
           method = method, 
           where = M_mask,
           printFlag = FALSE)
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(mice_result)) {
      cat("\nImputation failed at iteration", iteration, "\n")
      break
    }
    
    # Extract imputed matrix and make symmetric
    M_imputed <- as.matrix(complete(mice_result))
    M_imputed <- make_symmetric(M_imputed)
    
    # Calculate imputation distance
    if (!is.null(M_real)) {
      true_values <- M_real[M_mask]
      imputed_values <- M_imputed[M_mask]
      imputation_dis <- calculate_distance(true_values, imputed_values, method = distance_metric)
    } else {
      # No true values available
      imputation_dis <- NA
    }
    
    # Calculate propagation distance (between consecutive elements in tails-chain)
    # This is calculated BEFORE any averaging is applied
    if (iteration > 1) {
      prev_matrix <- tails_chain[[iteration]]
      prev_values <- prev_matrix[M_mask]
      curr_values <- M_imputed[M_mask]
      propagation_dis <- calculate_distance(prev_values, curr_values, method = distance_metric)
    } else {
      # First iteration - compare with initialized matrix
      init_values <- tails_chain[[1]][M_mask]
      curr_values <- M_imputed[M_mask]
      propagation_dis <- calculate_distance(init_values, curr_values, method = distance_metric)
    }
    
    # Conditional matrix handling based on propagation distance
    averaged <- FALSE
    if (iteration > 1) {
      # Check if propagation distance exceeds threshold
      if (distance_metric == "correlation") {
        # For correlation, higher is better - we average if NOT converged
        if (propagation_dis < (1 - propagation_convergence_threshold)) {
          # Not converged, so average with previous matrix
          M_final <- (M_imputed + tails_chain[[iteration]]) / 2
          M_final <- make_symmetric(M_final)  # Ensure symmetry
          averaged <- TRUE
        } else {
          # Converged, so use current matrix directly
          M_final <- M_imputed
        }
      } else {
        # For distance metrics, lower is better - we average if NOT converged
        if (propagation_dis >= propagation_convergence_threshold) {
          # Not converged, so average with previous matrix
          M_final <- (M_imputed + tails_chain[[iteration]]) / 2
          M_final <- make_symmetric(M_final)  # Ensure symmetry
          averaged <- TRUE
        } else {
          # Converged, so use current matrix directly
          M_final <- M_imputed
        }
      }
    } else {
      # First iteration always uses the direct result
      M_final <- M_imputed
    }
    
    # Add to tails-chain
    tails_chain[[iteration + 1]] <- M_final
    
    # Calculate improvement percentage in imputation distance
    if (iteration > 1) {
      prev_imputation_dis <- all_metrics$imputation_dis[nrow(all_metrics)]
      
      if (distance_metric == "correlation") {
        # For correlation, higher is better
        improvement_pct <- ifelse(!is.na(imputation_dis) && !is.na(prev_imputation_dis),
                                  100 * (imputation_dis - prev_imputation_dis) / abs(prev_imputation_dis), 0)
      } else {
        # For distance metrics, lower is better
        improvement_pct <- ifelse(!is.na(imputation_dis) && !is.na(prev_imputation_dis) && prev_imputation_dis > 0,
                                  100 * (prev_imputation_dis - imputation_dis) / prev_imputation_dis, 0)
      }
    } else {
      improvement_pct <- NA
    }
    
    # Calculate runtime for this iteration
    end_time_iter <- Sys.time()
    runtime <- as.numeric(difftime(end_time_iter, start_time_iter, units = "secs"))
    
    # Check convergence criteria based on BOTH distances
    is_converged <- FALSE
    
    # First check imputation distance if we have real values
    if (!is.na(imputation_dis) && !is.null(M_real)) {
      if (distance_metric == "correlation") {
        # For correlation, check if we're close to 1.0 (perfect correlation)
        if (imputation_dis > (1 - imputation_convergence_threshold)) {
          is_converged <- TRUE
        }
      } else {
        # For other metrics, check if below threshold
        if (imputation_dis < imputation_convergence_threshold) {
          is_converged <- TRUE
        }
      }
    }
    
    # Then check propagation distance
    if (!is.na(propagation_dis)) {
      if (distance_metric == "correlation") {
        # For correlation, check if change is minimal (close to 1.0 is good)
        if (propagation_dis > (1 - propagation_convergence_threshold)) {
          is_converged <- TRUE; 
        }
      } else {
        # For distance metrics, check if below threshold
        if (propagation_dis < propagation_convergence_threshold) {
          is_converged <- TRUE;
        }
      }
    }
    
    # Store metrics for this iteration
    new_row <- data.frame(
      k = k,
      bs = bs_i,
      iteration = iteration,
      imputation_dis = imputation_dis,
      propagation_dis = propagation_dis,
      runtime = runtime,
      improvement_pct = improvement_pct,
      converged = is_converged,
      averaged = averaged,
      stringsAsFactors = FALSE
    )
    
    all_metrics <- rbind(all_metrics, new_row)
    
    # Increment iteration counter
    iteration <- iteration + 1
  }
  
  # Minimal completion status
  cat("\nComplete\n")

  # Return results - metrics dataframe, final matrix, and all matrices in the tails-chain
  return(list(
    final_matrix = tails_chain[[length(tails_chain)]],
    metrics = all_metrics,
    tails_chain = tails_chain
  ))
}

# KNN Imputation Process --------------------------------------
library(DMwR2)
knn_impute <- function(M_input,
                       M_mask,
                       M_real = NULL,
                       distance_metric = "mae",
                       k = NA, bs_i = NA) {
  num_cols <- ncol(M_input)
  best_score <- ifelse(distance_metric == "correlation", -Inf, Inf)  # Initial value based on metric
  best_matrix <- NULL
  best_k <- NA
  
  # Convert to data frame for knnImputation()
  input_df <- as.data.frame(M_input)
  
  # Initialize storage for metrics with ALL columns matching iterative_imputation
  all_metrics <- data.frame(
    k = integer(),
    bs = integer(),
    iteration = integer(),      # This will be the k value in KNN
    imputation_dis = numeric(),
    propagation_dis = numeric(),
    runtime = numeric(),
    improvement_pct = numeric(),
    converged = logical(),
    averaged = logical(),
    stringsAsFactors = FALSE
  )
  
  # Store matrices for all k
  tails_chain <- list()
  
  # Initial matrix (similar to initialize_M in iterative_imputation)
  initialized_matrix <- initialize_M(M_input)
  tails_chain[[1]] <- initialized_matrix
  
  # Maximum k should not exceed number of complete rows
  max_k <- min(num_cols, sum(rowSums(is.na(input_df)) == 0))
  if (max_k < 1) max_k <- 1  # Ensure at least k=1 is tried
  
  # Previous score for improvement calculation
  prev_score <- NA
  
  for (curr_k in 1:max_k) {
    # Start time
    start_time <- Sys.time()
    
    # Perform KNN imputation with error handling
    imputed_df <- tryCatch({
      knnImputation(input_df, k = curr_k)
    }, error = function(e) {
      cat("Error with k =", curr_k, ":", e$message, "\n")
      return(NULL)
    })
    
    # Skip to next iteration if imputation failed
    if (is.null(imputed_df)) {
      next
    }
    
    # Convert to matrix and make symmetric
    imputed_matrix <- as.matrix(imputed_df)
    imputed_matrix <- make_symmetric(imputed_matrix)  # Added symmetry enforcement
    
    # Store the matrix in tails_chain
    tails_chain[[curr_k + 1]] <- imputed_matrix
    
    # Evaluate quality if real values are available
    if (!is.null(M_real)) {
      true_values <- M_real[M_mask]
      imputed_values <- imputed_matrix[M_mask]
      
      # Compute distance metric
      score <- calculate_distance(true_values, imputed_values, method = distance_metric)
      
      # Update best_k based on the metric
      if (!is.na(score)) {
        if ((distance_metric == "correlation" && score > best_score) ||  # Higher is better
            (distance_metric != "correlation" && score < best_score)) {  # Lower is better
          best_k <- curr_k
          best_score <- score
          best_matrix <- imputed_matrix
        }
      }
    } else {
      score <- NA
    }
    
    # Update previous score for next iteration
    prev_score <- score
    
    # Compute runtime
    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Store metrics (with all columns matching iterative_imputation)
    new_row <- data.frame(
      k = k,                      # Dataset k parameter
      bs = bs_i,                  # Dataset bootstrap index
      iteration = curr_k,         # KNN k parameter (stored as iteration)
      imputation_dis = score,
      propagation_dis = NA,
      runtime = runtime,
      improvement_pct = NA,
      converged = NA,
      averaged = FALSE,           # KNN never uses averaging
      stringsAsFactors = FALSE
    )
    
    all_metrics <- rbind(all_metrics, new_row)
  }
  
  cat("Best k selected:", best_k, "with", distance_metric, "score:", best_score, "\n")
  
  # Return matching structure to iterative_imputation
  return(list(
    final_matrix = best_matrix,
    metrics = all_metrics,
    tails_chain = tails_chain
  ))
}

# Single Dataset Processing --------------------------------------
impute_single_dataset <- function(input_path, method = "lasso.norm", 
                                  max_iterations = 5, 
                                  imputation_convergence_threshold = 1e-3,
                                  propagation_convergence_threshold = 1e-3,
                                  distance_metric = "mae",
                                  output_dir = "data/imputation_set") {
  
  # Extract dataset identifiers from the input path
  k <- as.numeric(gsub(".*masked_k([0-9]+)_bs.*", "\\1", input_path))
  bs_i <- as.numeric(gsub(".*masked_k[0-9]+_bs([0-9]+)\\.rds", "\\1", input_path))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create output paths
  method_safe <- gsub("\\.", "_", method)
  filled_path <- sprintf("%s/filled_k%d_bs%d_%s.rds", output_dir, k, bs_i, method_safe)
  intermediates_path <- sprintf("%s/intermediates_k%d_bs%d_%s.rds", output_dir, k, bs_i, method_safe)
  
  # Initialize
  empty_metrics <- data.frame(
    k = k,
    bs = bs_i,
    iteration = NA_integer_,
    imputation_dis = NA_real_,
    propagation_dis = NA_real_,
    runtime = NA_real_,
    improvement_pct = NA_real_,
    converged = FALSE,
    averaged = FALSE,
    masked_cells = NA_integer_,
    masked_percent = NA_real_,
    distance_metric = distance_metric,
    impt_method = method,
    imputation_convergence_threshold = imputation_convergence_threshold,
    propagation_convergence_threshold = propagation_convergence_threshold,
    input_path = input_path,
    output_path = filled_path,
    intermediates_path = intermediates_path,
    stringsAsFactors = FALSE
  )
  
  # Load datasets
  M_input <- tryCatch({
    readRDS(input_path)
  }, error = function(e) {
    cat(paste("Failed to load input matrix:", input_path, "-", e$message), "\n")
    return(empty_metrics)
  })
  
  mask_path <- gsub("training_set/masked_", "position_set/mask_", input_path)
  M_mask <- tryCatch({
    readRDS(mask_path)
  }, error = function(e) {
    cat(paste("Failed to load mask matrix:", mask_path, "-", e$message), "\n")
    return(empty_metrics)
  })
  
  # Calculate masked cells information
  masked_cells <- sum(M_mask)
  masked_percent <- 100 * masked_cells / (nrow(M_input) * ncol(M_input))
  
  empty_metrics$masked_cells <- masked_cells
  empty_metrics$masked_percent <- masked_percent
  
  # Load the true values
  full_dataset_path <- "data/full_dataset.rds"
  M_real <- tryCatch({
    if (file.exists(full_dataset_path)) {
      readRDS(full_dataset_path)
    } else {
      NULL
    }
  }, error = function(e) {
    NULL
  })
  
  # Run the imputation
  imputation_result <- tryCatch({
    if (toupper(method) == "KNN") {
      # Use KNN imputation
      cat(sprintf("Starting KNN imputation (k=%s, bs=%s, metric=%s)\n", k, bs_i, distance_metric))
      knn_impute(
        M_input = M_input,
        M_mask = M_mask,
        M_real = M_real,
        distance_metric = distance_metric,
        k = k,
        bs_i = bs_i
      )
    } else {
      # Use iterative imputation with MICE
      iterative_imputation(
        M_input = M_input,
        M_mask = M_mask,
        M_real = M_real,
        method = method,
        max_iterations = max_iterations,
        imputation_convergence_threshold = imputation_convergence_threshold,
        propagation_convergence_threshold = propagation_convergence_threshold,
        distance_metric = distance_metric,
        k = k,
        bs_i = bs_i
      )
    }
  }, error = function(e) {
    cat("Imputation failed:", e$message, "\n")
    return(list(
      final_matrix = NULL,
      metrics = empty_metrics,
      tails_chain = NULL
    ))
  })
  
  # Process metrics from imputation
  if (!is.null(imputation_result$metrics) && nrow(imputation_result$metrics) > 0) {
    # Ensure all required columns exist in the output
    metrics_df <- imputation_result$metrics
    result_metrics <- empty_metrics[0, ]
    
    # Add all rows from imputation metrics
    for (i in 1:nrow(metrics_df)) {
      row <- metrics_df[i, ]
      
      new_row <- empty_metrics[1, ]
      
      for (col in names(metrics_df)) {
        if (col %in% names(new_row)) {
          new_row[[col]] <- row[[col]]
        }
      }
      
      result_metrics <- rbind(result_metrics, new_row)
    }
    
    # Update metadata
    result_metrics$masked_cells <- masked_cells
    result_metrics$masked_percent <- masked_percent
    result_metrics$distance_metric <- distance_metric
    result_metrics$impt_method <- method
    result_metrics$imputation_convergence_threshold <- imputation_convergence_threshold
    result_metrics$propagation_convergence_threshold <- propagation_convergence_threshold
    result_metrics$input_path <- input_path
    result_metrics$output_path <- filled_path
    result_metrics$intermediates_path <- intermediates_path
  } else {
    # Use the empty metrics template if no metrics were produced
    result_metrics <- empty_metrics
  }
  
  # Save results with error handling
  tryCatch({
    if (!is.null(imputation_result$final_matrix)) {
      saveRDS(imputation_result$final_matrix, filled_path)

      if (!exists("result_metrics")) {
        result_metrics <- empty_metrics
      }
      
      # Save all matrices in the tails-chain and metrics
      intermediates <- list(
        tails_chain = imputation_result$tails_chain,
        metrics = result_metrics,
        runtime = sum(result_metrics$runtime, na.rm = TRUE),
        converged = any(result_metrics$converged, na.rm = TRUE),
        iterations = max(result_metrics$iteration, na.rm = TRUE)
      )
      saveRDS(intermediates, intermediates_path)
    }
  }, error = function(e) {
    cat("Error saving results:", e$message, "\n")
  })
  
  return(result_metrics)
}

# debug
if (F) {
  result_df <- impute_single_dataset(
    input_path = "data/training_set/masked_k5_bs1.rds",
    method = "lasso.norm",
    max_iterations = 5,
    imputation_convergence_threshold = 1e-3,
    propagation_convergence_threshold = 1e-3,
    distance_metric = "mae",
    output_dir = "data/imputation_set"
  )
}