# BESMI Framework - 03-Batch_Processing.r
# Batch Processing for BESMI Imputation Framework
# V1.0.0
# Author: Zhu, Jiashuai (The University of Melbourne; Agriculture Victoria Research)

# Load required libraries
library(tidyverse)
library(mice)
library(reshape2)
library(parallel)
library(memoise)
library(pryr)

source("02-Imputation.r")

# Set up output directories
if (!dir.exists("data/imputation_set")) {
  dir.create("data/imputation_set", recursive = TRUE)
}
if (!dir.exists("data/summary_tables")) {
  dir.create("data/summary_tables", recursive = TRUE)
}

# Function to log memory usage
log_memory <- function(stage, dataset_name, log_file = "data/batch.log") {
  mem_usage <- as.numeric(pryr::mem_used()) / 1024^2  # Convert to MB
  message <- sprintf("%s - Memory usage for %s: %.2f MB", 
                     format(Sys.time(), "%H:%M:%S"), 
                     dataset_name, mem_usage)
  write(message, log_file, append = TRUE)
  return(mem_usage)
}

# Main function for batch imputation processing
run_batch_imputation <- function(method = "lasso.norm", 
                                 max_iter = 5,
                                 imputation_convergence_threshold = 1e-6,
                                 propagation_convergence_threshold = 1e-6,
                                 distance_metric = "mae",
                                 output_dir = "data/imputation_set",
                                 k_filter = NULL) {
  
  # method-safe string for filenames
  method_safe <- gsub("\\.", "_", method)
  
  # Dataset -------------------------------------------------
  
  # Start logging
  log_file <- "data/batch.log"
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  write(paste("=== BESMI Batch Processing Started at", timestamp, "==="), log_file)
  write(paste("Method:", method), log_file, append = TRUE)
  write(paste("Max iterations:", max_iter), log_file, append = TRUE)
  write(paste("Distance metric:", distance_metric), log_file, append = TRUE)
  write(paste("Imputation convergence threshold:", imputation_convergence_threshold), log_file, append = TRUE)
  write(paste("Propagation convergence threshold:", propagation_convergence_threshold), log_file, append = TRUE)
  
  # Load dataset summary table
  summary_path <- "data/summary_tables/summary_table_datasets.csv"
  tryCatch({
    summary_table <- read.csv(summary_path, stringsAsFactors = FALSE)
  }, error = function(e) {
    message <- paste("Failed to load dataset summary table:", e$message)
    write(message, log_file, append = TRUE)
    stop(message)
  })
  
  # Extract input paths
  dataset_paths <- summary_table$training_path
  
  # Validate that all files exist
  missing_files <- dataset_paths[!file.exists(dataset_paths)]
  if (length(missing_files) > 0) {
    message <- paste("Some dataset files are missing:", paste(missing_files, collapse = ", "))
    write(message, log_file, append = TRUE)
    stop(message)
  }
  
  # Filter by k value if specified
  if (!is.null(k_filter)) {
    k_values <- as.numeric(gsub(".*masked_k([0-9]+)_bs.*", "\\1", dataset_paths))
    dataset_paths <- dataset_paths[k_values == k_filter]
    if (length(dataset_paths) == 0) {
      message <- paste("No datasets found with k =", k_filter)
      write(message, log_file, append = TRUE)
      stop(message)
    }
    write(paste("Filtered to", length(dataset_paths), "datasets with k =", k_filter), 
          log_file, append = TRUE)
  }
  
  # Initialize -----------------------------------------
  Metrics_df <- data.frame(
    k = integer(),
    bs = integer(),
    iteration = integer(),
    imputation_dis = numeric(),
    propagation_dis = numeric(),
    runtime = numeric(),
    improvement_pct = numeric(),
    converged = logical(),
    averaged = logical(),
    masked_cells = integer(),
    masked_percent = numeric(),
    distance_metric = character(),
    impt_method = character(),
    imputation_convergence_threshold = numeric(),
    propagation_convergence_threshold = numeric(),
    input_path = character(),
    output_path = character(),
    intermediates_path = character(),
    stringsAsFactors = FALSE
  )
  
  # Checkpointing Mechanism --------------------------------------------
  checkpoint_file <- paste0("data/batch_checkpoint_", method_safe, ".rds")
  processed_datasets <- c()
  
  # previous metrics
  previous_metrics <- data.frame()
  metrics_filename <- paste0("data/metrics_df_", method_safe, ".csv")
  
  if (file.exists(metrics_filename)) {
    tryCatch({
      previous_metrics <- read.csv(metrics_filename, stringsAsFactors = FALSE)
      write(paste("Loaded previous metrics with", nrow(previous_metrics), "rows"), 
            log_file, append = TRUE)
    }, error = function(e) {
      write(paste("Error loading previous metrics:", e$message), 
            log_file, append = TRUE)
    })
  }
  
  # Check if checkpoint file exists
  if (file.exists(checkpoint_file)) {
    tryCatch({
      processed_datasets <- readRDS(checkpoint_file)
      if (length(processed_datasets) > 0) {
        write(paste("Resuming from checkpoint for method", method, 
                    "- Already processed", length(processed_datasets), "datasets"), 
              log_file, append = TRUE)
      }
    }, error = function(e) {
      write(paste("Error loading checkpoint file:", e$message), log_file, append = TRUE)
      processed_datasets <- c()
    })
  } else {
    write(paste("No existing checkpoint found for method", method, 
                "- Starting fresh processing"), log_file, append = TRUE)
  }
  
  # Filter out already processed datasets
  datasets_to_process <- dataset_paths[!dataset_paths %in% processed_datasets]
  
  # Set up progress tracking
  total_datasets <- length(datasets_to_process)
  write(paste("Processing", total_datasets, "datasets"), log_file, append = TRUE)
  
  # *Batch Processing ---------------------------------------------------
  
  # Call batch_impute for dataset processing
  result <- batch_impute(
    dataset_paths = datasets_to_process,
    the_method = method,
    max_iter = max_iter,
    imputation_convergence_threshold = imputation_convergence_threshold,
    propagation_convergence_threshold = propagation_convergence_threshold,
    distance_metric = distance_metric,
    output_dir = output_dir,
    log_file = log_file,
    processed_datasets = processed_datasets,
    checkpoint_file = checkpoint_file
  )
  
  # Combine with any previous metrics
  if (nrow(previous_metrics) > 0) {
    Metrics_df <- rbind(previous_metrics, result)
    write(paste("Combined with previous metrics. Total rows:", nrow(Metrics_df)), 
          log_file, append = TRUE)
  } else {
    Metrics_df <- result
  }
  
  # Generate Summary Statistics ----------------------------------------
  
  # Group the combined metrics by dataset
  dataset_groups <- Metrics_df %>%
    group_by(k, bs) %>%
    summarize(
      chain_length = max(iteration, na.rm = TRUE),
      chain_converged = any(converged, na.rm = TRUE),
      total_runtime = sum(runtime, na.rm = TRUE),
      final_impt_dis = imputation_dis[which.max(iteration)],
      avg_impt_dis = mean(imputation_dis, na.rm = TRUE),
      final_ppgt_dis = propagation_dis[which.max(iteration)],
      avg_ppgt_dis = mean(propagation_dis, na.rm = TRUE),
      best_imputation_dis = if(first(distance_metric) == "correlation") {
        max(imputation_dis, na.rm = TRUE)
      } else {
        min(imputation_dis, na.rm = TRUE)
      },
      # Transfer metadata that is dataset-specific
      masked_cells = first(masked_cells),
      masked_percent = first(masked_percent),
      distance_metric = first(distance_metric),
      impt_method = first(impt_method),
      imputation_convergence_threshold = first(imputation_convergence_threshold),
      propagation_convergence_threshold = first(propagation_convergence_threshold),
      input_path = first(input_path),
      output_path = first(output_path),
      intermediates_path = first(intermediates_path)
    )
  
  summary_table_impt <- as.data.frame(dataset_groups)
  
  # Results Consolidation ----------------------------------------------
  
  metrics_filename <- paste0("data/metrics_df_", method_safe, ".csv")
  summary_filename <- paste0("data/summary_tables/summary_table_imputation_", method_safe, ".csv")
  
  # Save summary table and metrics
  saveRDS(Metrics_df, paste0("data/metrics_df_", method_safe, ".rds"))
  write.csv(Metrics_df, metrics_filename, row.names = FALSE)
  write.csv(summary_table_impt, summary_filename, row.names = FALSE)
  
  # Log filenames
  write(paste("Saved metrics to:", metrics_filename), log_file, append = TRUE)
  write(paste("Saved summary to:", summary_filename), log_file, append = TRUE)
  
  # Log completion summary
  timestamp_end <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  write(paste("=== BESMI Batch Processing Completed at", timestamp_end, "==="), 
        log_file, append = TRUE)
  write(paste("Total datasets processed in this run:", nrow(result)), 
        log_file, append = TRUE)
  write(paste("Total datasets processed overall:", nrow(summary_table_impt)), 
        log_file, append = TRUE)
  write(paste("Total converged datasets:", sum(summary_table_impt$chain_converged)), 
        log_file, append = TRUE)
  write(paste("Average runtime per dataset:", 
              mean(summary_table_impt$total_runtime, na.rm = TRUE), "seconds"), 
        log_file, append = TRUE)
}

# Function to execute batch imputation on datasets
batch_impute <- function(dataset_paths, 
                         the_method = "lasso.norm",
                         max_iter = 5,
                         imputation_convergence_threshold = 1e-6,
                         propagation_convergence_threshold = 1e-6,
                         distance_metric = "mae",
                         output_dir = "data/imputation_set",
                         log_file = "data/batch.log",
                         processed_datasets = c(),
                         checkpoint_file = "data/batch_checkpoint.rds",
                         k_filter = NULL) {
  
  # Initialize result dataframe
  combined_metrics <- data.frame()
  
  # Current progress count
  current <- 0
  total <- length(dataset_paths)
  
  # Progress bar
  cat("\nBatch Processing Progress:\n")
  progress_bar <- function(current, total, width = 50) {
    percent <- round(current / total * 100)
    filled <- round(width * current / total)
    bar <- paste(rep("=", filled), collapse = "")
    if (filled < width) {
      bar <- paste0(bar, ">", paste(rep(" ", width - filled - 1), collapse = ""))
    }
    cat(sprintf("\r[%s] %d%% (%d/%d)", bar, percent, current, total))
    if (current == total) cat("\n")
  }
  
  # If there are no datasets to process, return empty dataframe
  if (length(dataset_paths) == 0) {
    write("No datasets to process in this run.", log_file, append = TRUE)
    return(data.frame())
  }
  
  # Process each dataset sequentially
  for (dataset_path in dataset_paths) {
    # Extract k and bs values from filename
    k <- as.numeric(gsub(".*masked_k([0-9]+)_bs.*", "\\1", dataset_path))
    bs_i <- as.numeric(gsub(".*masked_k[0-9]+_bs([0-9]+)\\.rds", "\\1", dataset_path))
    
    # Only process if k_filter matches or is NULL
    if (!is.null(k_filter) && k != k_filter) {
      next
    }
    
    # Log start of processing for this dataset
    dataset_name <- basename(dataset_path)
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    write(paste(timestamp, "- Starting dataset:", dataset_name), log_file, append = TRUE)
    
    # Log initial memory usage
    initial_mem <- log_memory("Before", dataset_name, log_file)
    
    # Process the dataset with error handling
    result <- tryCatch({
      # Call the imputation function from 02-Imputation.r
      impute_single_dataset(
        input_path = dataset_path,
        method = the_method,
        max_iterations = max_iter,
        imputation_convergence_threshold = imputation_convergence_threshold,
        propagation_convergence_threshold = propagation_convergence_threshold,
        distance_metric = distance_metric,
        output_dir = output_dir
      )
    }, error = function(e) {
      # Log error
      error_msg <- paste("Error processing", dataset_name, ":", e$message)
      write(error_msg, log_file, append = TRUE)
      
      # Return empty dataframe with error information
      data.frame(
        k = k,
        bs = bs_i,
        iteration = NA,
        imputation_dis = NA,
        propagation_dis = NA,
        runtime = NA,
        improvement_pct = NA,
        converged = FALSE,
        averaged = NA,
        masked_cells = NA,
        masked_percent = NA,
        distance_metric = distance_metric,
        impt_method = the_method,
        imputation_convergence_threshold = imputation_convergence_threshold,
        propagation_convergence_threshold = propagation_convergence_threshold,
        input_path = dataset_path,
        output_path = NA,
        intermediates_path = NA,
        error = e$message,
        stringsAsFactors = FALSE
      )
    })
    
    # Log memory usage after processing
    final_mem <- log_memory("After", dataset_name, log_file)
    mem_diff <- final_mem - initial_mem
    write(sprintf("%s - Memory change for %s: %.2f MB", 
                  format(Sys.time(), "%H:%M:%S"), 
                  dataset_name, mem_diff), log_file, append = TRUE)
    
    # Combine results
    combined_metrics <- rbind(combined_metrics, result)
    
    # Update processed datasets list and save checkpoint
    processed_datasets <- c(processed_datasets, dataset_path)
    saveRDS(processed_datasets, checkpoint_file)
    
    # Update progress
    current <- current + 1
    progress_bar(current, total)
    
    # Log completion of this dataset
    timestamp_end <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    iterations <- max(result$iteration, na.rm = TRUE)
    total_runtime <- sum(result$runtime, na.rm = TRUE)
    converged <- any(result$converged, na.rm = TRUE)
    write(paste(timestamp_end, "- Completed dataset:", dataset_name, 
                "- Iterations:", iterations, 
                "- Runtime:", total_runtime, "sec", 
                "- Converged:", converged), 
          log_file, append = TRUE)
    
    # garbage collection to free memory
    gc()
  }
  
  return(combined_metrics)
}

# RUN
if (T) {
  run_batch_imputation(
    method = "midastouch",
    max_iter = 5,
    # k_filter = 5,  # Only process datasets with k=5
    imputation_convergence_threshold = 1e-3,
    propagation_convergence_threshold = 1e-3,
    distance_metric = "mae",
  )
}