# BESMI: Bootstrapping Evaluation for Structural Missingness Imputation

#' Prepare full GDM dataset from CSV or RData
#' @param input_path Path to CSV or RData file containing the full distance matrix
#' @return Symmetric numeric matrix
#' @export
besmi_prepare_full_dataset <- function(input_path) {
  if (!file.exists(input_path)) stop(paste("Input file not found:", input_path))
  if (grepl("\\.csv$", input_path)) {
    gdm_data <- utils::read.csv(input_path, row.names = 1, check.names = FALSE)
    gdm <- as.matrix(gdm_data)
  } else if (grepl("\\.RData$", input_path)) {
    e <- new.env(parent = emptyenv())
    load(input_path, envir = e)
    if (!exists("gdm", envir = e)) stop("GDM not found in the loaded data. Expect object 'gdm'.")
    gdm <- get("gdm", envir = e)
  } else stop("Unsupported file format. Provide CSV or RData.")
  if (!is.matrix(gdm)) gdm <- as.matrix(gdm)
  if (!isSymmetric(gdm)) stop("The genetic distance matrix is not symmetric.")
  if (any(is.na(gdm))) stop("The genetic distance matrix already contains missing values.")
  gdm
}

#' Determine bootstrap sample count for a given k
#' @keywords internal
.besmi_determine_sampling_sizes <- function(k) {
  if (k <= 3) 15 else if (k <= 6) 8 else 4
}

#' Create masked matrices for BESMI
#' @param full_matrix Full symmetric matrix
#' @param k Number of populations to mask (as U)
#' @param seed Optional seed for reproducibility
#' @return List with masked_matrix, mask_position, group_u, group_s, masked_percentage
#' @export
besmi_create_masked_matrices <- function(full_matrix, k, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pop_names <- rownames(full_matrix)
  n_pops <- length(pop_names)
  group_u_indices <- sample(seq_len(n_pops), k)
  group_u <- pop_names[group_u_indices]
  group_s <- pop_names[-group_u_indices]
  masked_matrix <- full_matrix
  mask_position <- matrix(FALSE, nrow = nrow(masked_matrix), ncol = ncol(masked_matrix))
  rownames(mask_position) <- pop_names
  colnames(mask_position) <- pop_names
  for (i in group_u) for (j in group_u) { masked_matrix[i, j] <- NA; mask_position[i, j] <- TRUE }
  total_cells <- nrow(masked_matrix) * ncol(masked_matrix)
  masked_cells <- sum(mask_position)
  masked_percentage <- round(100 * masked_cells / total_cells, 1)
  list(masked_matrix = masked_matrix, mask_position = mask_position, group_u = group_u, group_s = group_s, masked_percentage = masked_percentage)
}

#' Initialize matrix by column means
#' @keywords internal
.besmi_initialize_M <- function(M) {
  df <- as.data.frame(M)
  for (col in colnames(df)) df[[col]][is.na(df[[col]])] <- mean(df[[col]], na.rm = TRUE)
  M_init <- as.matrix(df)
  (M_init + t(M_init)) / 2
}

#' Distance metrics
#' @keywords internal
.besmi_calculate_distance <- function(a, b, method = "mae") {
  if (length(a) != length(b)) stop("Vectors must have the same length")
  valid <- !is.na(a) & !is.na(b)
  if (sum(valid) < 3) return(NA)
  if (method == "ssd") sum((a - b)^2, na.rm = TRUE)
  else if (method == "mae") mean(abs(a - b), na.rm = TRUE)
  else if (method == "correlation") stats::cor(a, b, method = "pearson", use = "pairwise.complete.obs")
  else if (method == "rmse") sqrt(mean((a - b)^2, na.rm = TRUE))
  else stop("Unknown distance method")
}

#' Iterative imputation with MICE (tails-chain)
#' @param M_input Matrix with NAs to impute
#' @param M_mask Logical mask matrix (TRUE indicates masked positions)
#' @param M_real Optional ground truth matrix
#' @param method MICE method (e.g., 'lasso.norm')
#' @param max_iterations Max outer iterations
#' @param imputation_convergence_threshold Threshold for imputation distance
#' @param propagation_convergence_threshold Threshold for propagation distance
#' @param distance_metric Distance metric name
#' @param k Dataset parameter k (for logging)
#' @param bs_i Bootstrap index (for logging)
#' @return List with final_matrix, metrics, tails_chain
#' @export
besmi_iterative_imputation <- function(M_input,
                                 M_mask,
                                 M_real = NULL, 
                                 method = "lasso.norm", 
                                 max_iterations = 5, 
                                 imputation_convergence_threshold = 1e-3,
                                 propagation_convergence_threshold = 1e-3,
                                 distance_metric = "mae",
                                 k = NA, bs_i = NA) {
  iteration <- 1; all_metrics <- data.frame(); tails_chain <- list()
  M_initialized <- .besmi_initialize_M(M_input)
  tails_chain[[1]] <- M_initialized
  while (iteration <= max_iterations) {
    input_df <- as.data.frame(M_input)
    mice_result <- tryCatch({
      mice::mice(input_df, m = 1, maxit = iteration, method = method, where = M_mask, printFlag = FALSE)
    }, error = function(e) NULL)
    if (is.null(mice_result)) break
    M_imputed <- as.matrix(mice::complete(mice_result))
    M_imputed <- (M_imputed + t(M_imputed)) / 2
    if (!is.null(M_real)) {
      true_values <- M_real[M_mask]
      imputed_values <- M_imputed[M_mask]
      imputation_dis <- .besmi_calculate_distance(true_values, imputed_values, method = distance_metric)
    } else imputation_dis <- NA
    if (iteration > 1) {
      prev_matrix <- tails_chain[[iteration]]
      prev_values <- prev_matrix[M_mask]
      curr_values <- M_imputed[M_mask]
      propagation_dis <- .besmi_calculate_distance(prev_values, curr_values, method = distance_metric)
    } else {
      init_values <- tails_chain[[1]][M_mask]
      curr_values <- M_imputed[M_mask]
      propagation_dis <- .besmi_calculate_distance(init_values, curr_values, method = distance_metric)
    }
    averaged <- FALSE
    if (iteration > 1) {
      if (distance_metric == "correlation") {
        if (propagation_dis < (1 - propagation_convergence_threshold)) {
          M_final <- (M_imputed + tails_chain[[iteration]]) / 2; M_final <- (M_final + t(M_final)) / 2; averaged <- TRUE
        } else M_final <- M_imputed
      } else {
        if (!is.na(propagation_dis) && propagation_dis >= propagation_convergence_threshold) {
          M_final <- (M_imputed + tails_chain[[iteration]]) / 2; M_final <- (M_final + t(M_final)) / 2; averaged <- TRUE
        } else M_final <- M_imputed
      }
    } else M_final <- M_imputed
    tails_chain[[iteration + 1]] <- M_final
    new_row <- data.frame(k = k, bs = bs_i, iteration = iteration,
                          imputation_dis = imputation_dis, propagation_dis = propagation_dis,
                          runtime = NA_real_, improvement_pct = NA_real_,
                          converged = FALSE, averaged = averaged, stringsAsFactors = FALSE)
    all_metrics <- rbind(all_metrics, new_row)
    iteration <- iteration + 1
  }
  list(final_matrix = tails_chain[[length(tails_chain)]], metrics = all_metrics, tails_chain = tails_chain)
}

#' KNN imputation sweep (uses VIM::kNN)
#' @param M_input Matrix with NAs
#' @param M_mask Logical mask matrix
#' @param M_real Optional ground truth
#' @param distance_metric Distance metric name
#' @param k Dataset parameter k
#' @param bs_i Bootstrap index
#' @return List with final_matrix, metrics, tails_chain
#' @export
besmi_knn_impute <- function(M_input, M_mask, M_real = NULL, distance_metric = "mae", k = NA, bs_i = NA) {
  if (!requireNamespace("VIM", quietly = TRUE)) stop("VIM package is required for KNN imputation. Please install VIM from CRAN.")
  num_cols <- ncol(M_input)
  best_score <- ifelse(distance_metric == "correlation", -Inf, Inf)
  best_matrix <- NULL; best_k <- NA
  input_df <- as.data.frame(M_input)
  all_metrics <- data.frame(); tails_chain <- list()
  initialized_matrix <- .besmi_initialize_M(M_input)
  tails_chain[[1]] <- initialized_matrix
  max_k <- min(num_cols, sum(rowSums(is.na(input_df)) == 0)); if (max_k < 1) max_k <- 1
  for (curr_k in seq_len(max_k)) {
    imputed_df <- tryCatch({ VIM::kNN(input_df, k = curr_k, imp_var = FALSE) }, error = function(e) NULL)
    if (is.null(imputed_df)) next
    imputed_matrix <- as.matrix(imputed_df)
    imputed_matrix <- (imputed_matrix + t(imputed_matrix)) / 2
    tails_chain[[curr_k + 1]] <- imputed_matrix
    if (!is.null(M_real)) {
      true_values <- M_real[M_mask]; imputed_values <- imputed_matrix[M_mask]
      score <- .besmi_calculate_distance(true_values, imputed_values, method = distance_metric)
      if (!is.na(score)) {
        if ((distance_metric == "correlation" && score > best_score) || (distance_metric != "correlation" && score < best_score)) {
          best_k <- curr_k; best_score <- score; best_matrix <- imputed_matrix
        }
      }
    } else score <- NA
    new_row <- data.frame(k = k, bs = bs_i, iteration = curr_k, imputation_dis = score,
                          propagation_dis = NA, runtime = NA, improvement_pct = NA,
                          converged = NA, averaged = FALSE, stringsAsFactors = FALSE)
    all_metrics <- rbind(all_metrics, new_row)
  }
  if (is.null(M_real)) {
    best_score <- NA_real_
  }
  list(
    final_matrix = best_matrix,
    metrics = all_metrics,
    tails_chain = tails_chain,
    best_k = best_k,
    best_score = best_score
  )
}
