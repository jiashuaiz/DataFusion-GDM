# BESMI Batch Processing helpers

#' Run BESMI imputation for a list of dataset paths
#' @param dataset_paths Character vector of RDS paths to masked matrices
#' @param the_method Imputation method (e.g., 'lasso.norm' or 'KNN')
#' @param max_iter Maximum iterations for iterative methods
#' @param imputation_convergence_threshold Convergence threshold for imputation metric
#' @param propagation_convergence_threshold Convergence threshold for propagation metric
#' @param distance_metric Distance metric for evaluation ('mae','ssd','rmse','correlation')
#' @param output_dir Output directory for imputed matrices (defaults to a temporary location)
#' @param k_filter Optional numeric filter for k value
#' @param full_dataset_path Optional path to a full matrix RDS used as ground truth
#' @return Data frame of metrics for all datasets
#' @export
besmi_batch_impute <- function(dataset_paths,
                               the_method = "lasso.norm",
                               max_iter = 5,
                               imputation_convergence_threshold = 1e-6,
                               propagation_convergence_threshold = 1e-6,
                               distance_metric = "mae",
                               output_dir = file.path(tempdir(), "DataFusionGDM_imputation"),
                               k_filter = NULL,
                               full_dataset_path = NULL) {
  combined_metrics <- data.frame()
  if (length(dataset_paths) == 0) return(combined_metrics)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  for (dataset_path in dataset_paths) {
    k <- as.numeric(gsub(".*masked_k([0-9]+)_bs.*", "\\1", dataset_path))
    bs_i <- as.numeric(gsub(".*masked_k[0-9]+_bs([0-9]+)\\.rds", "\\1", dataset_path))
    if (!is.null(k_filter) && k != k_filter) next
    res <- tryCatch({
      besmi_impute_single_dataset(dataset_path, method = the_method, max_iterations = max_iter,
                                  imputation_convergence_threshold = imputation_convergence_threshold,
                                  propagation_convergence_threshold = propagation_convergence_threshold,
                                  distance_metric = distance_metric, output_dir = output_dir,
                                  full_dataset_path = full_dataset_path)
    }, error = function(e) {
      data.frame(k = k, bs = bs_i, iteration = NA, imputation_dis = NA, propagation_dis = NA,
                 runtime = NA, improvement_pct = NA, converged = FALSE, averaged = NA,
                 masked_cells = NA, masked_percent = NA, distance_metric = distance_metric,
                 impt_method = the_method, imputation_convergence_threshold = imputation_convergence_threshold,
                 propagation_convergence_threshold = propagation_convergence_threshold,
                 input_path = dataset_path, output_path = NA, intermediates_path = NA,
                 error = e$message, stringsAsFactors = FALSE)
    })
    combined_metrics <- rbind(combined_metrics, res)
    gc()
  }
  combined_metrics
}

#' Impute a single dataset from masked matrix path
#' @param input_path Path to masked matrix RDS
#' @param method Imputation method ('lasso.norm' or 'KNN')
#' @param max_iterations Maximum iterations for iterative methods
#' @param imputation_convergence_threshold Convergence threshold for imputation metric
#' @param propagation_convergence_threshold Convergence threshold for propagation metric
#' @param distance_metric Distance metric name
#' @param output_dir Output directory for results (defaults to a temporary location)
#' @param full_dataset_path Optional path to a full matrix RDS used as ground truth
#' @return Data frame of per-iteration metrics
#' @export
besmi_impute_single_dataset <- function(input_path, method = "lasso.norm", 
                                  max_iterations = 5, 
                                  imputation_convergence_threshold = 1e-3,
                                  propagation_convergence_threshold = 1e-3,
                                  distance_metric = "mae",
                                  output_dir = file.path(tempdir(), "DataFusionGDM_imputation"),
                                  full_dataset_path = NULL) {
  k <- as.numeric(gsub(".*masked_k([0-9]+)_bs.*", "\\1", input_path))
  bs_i <- as.numeric(gsub(".*masked_k[0-9]+_bs([0-9]+)\\.rds", "\\1", input_path))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  method_safe <- gsub("\\.", "_", method)
  filled_path <- sprintf("%s/filled_k%d_bs%d_%s.rds", output_dir, k, bs_i, method_safe)
  intermediates_path <- sprintf("%s/intermediates_k%d_bs%d_%s.rds", output_dir, k, bs_i, method_safe)
  empty_metrics <- data.frame(
    k = k, bs = bs_i, iteration = NA_integer_, imputation_dis = NA_real_, propagation_dis = NA_real_,
    runtime = NA_real_, improvement_pct = NA_real_, converged = FALSE, averaged = FALSE,
    masked_cells = NA_integer_, masked_percent = NA_real_, distance_metric = distance_metric, impt_method = method,
    imputation_convergence_threshold = imputation_convergence_threshold,
    propagation_convergence_threshold = propagation_convergence_threshold,
    input_path = input_path, output_path = filled_path, intermediates_path = intermediates_path,
    stringsAsFactors = FALSE
  )
  M_input <- tryCatch({ readRDS(input_path) }, error = function(e) return(empty_metrics))
  mask_path <- gsub("training_set/masked_", "position_set/mask_", input_path)
  M_mask <- tryCatch({ readRDS(mask_path) }, error = function(e) return(empty_metrics))
  masked_cells <- sum(M_mask)
  masked_percent <- 100 * masked_cells / (nrow(M_input) * ncol(M_input))
  empty_metrics$masked_cells <- masked_cells
  empty_metrics$masked_percent <- masked_percent
  M_real <- NULL
  if (!is.null(full_dataset_path)) {
    M_real <- tryCatch({ if (file.exists(full_dataset_path)) readRDS(full_dataset_path) else NULL }, error = function(e) NULL)
  }
  imputation_result <- tryCatch({
    if (toupper(method) == "KNN") {
      besmi_knn_impute(M_input = M_input, M_mask = M_mask, M_real = M_real, distance_metric = distance_metric, k = k, bs_i = bs_i)
    } else {
      besmi_iterative_imputation(M_input = M_input, M_mask = M_mask, M_real = M_real, method = method,
        max_iterations = max_iterations, imputation_convergence_threshold = imputation_convergence_threshold,
        propagation_convergence_threshold = propagation_convergence_threshold, distance_metric = distance_metric,
        k = k, bs_i = bs_i)
    }
  }, error = function(e) list(final_matrix = NULL, metrics = empty_metrics, tails_chain = NULL))
  if (!is.null(imputation_result$final_matrix)) {
    saveRDS(imputation_result$final_matrix, filled_path)
    intermediates <- list(tails_chain = imputation_result$tails_chain, metrics = imputation_result$metrics)
    saveRDS(intermediates, intermediates_path)
  }
  if (!is.null(imputation_result$metrics) && nrow(imputation_result$metrics) > 0) {
    metrics_df <- imputation_result$metrics
    result_metrics <- empty_metrics[0, ]
    for (i in seq_len(nrow(metrics_df))) {
      row <- metrics_df[i, ]; new_row <- empty_metrics[1, ]
      for (col in names(metrics_df)) if (col %in% names(new_row)) new_row[[col]] <- row[[col]]
      result_metrics <- rbind(result_metrics, new_row)
    }
    result_metrics$masked_cells <- masked_cells
    result_metrics$masked_percent <- masked_percent
    result_metrics$distance_metric <- distance_metric
    result_metrics$impt_method <- method
    result_metrics$imputation_convergence_threshold <- imputation_convergence_threshold
    result_metrics$propagation_convergence_threshold <- propagation_convergence_threshold
    result_metrics$input_path <- input_path
    result_metrics$output_path <- filled_path
    result_metrics$intermediates_path <- intermediates_path
  } else result_metrics <- empty_metrics
  result_metrics
}
