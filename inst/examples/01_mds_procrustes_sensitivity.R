# Example: MDS + Procrustes sensitivity using package APIs

library(DataFusionGDM)
library(ggplot2)
library(vegan)

set.seed(42)

# Prefer packaged example data
pkg_file <- system.file("extdata", "GDM_simulated.csv", package = "DataFusionGDM")
if (nzchar(pkg_file)) {
  input_file <- pkg_file
} else {
  if (!file.exists("data/GDM_simulated.csv")) {
    dir.create("data", showWarnings = FALSE)
    export_simulated_gdm("data/GDM_simulated.csv", scenario = "default", n_pops = 40)
  }
  input_file <- "data/GDM_simulated.csv"
}

# Parameters
k_values <- seq(2, 30, by = 2)
num_tests <- 10
seed_base <- 42

# Data prep function equivalent
prepare_matrices_from_csv <- function(input_file, k, bias_mu = 0.1, noise_sd = 0.05, randomize_shared = TRUE, seed = 42) {
  full_matrix <- as.matrix(read.csv(input_file, header = TRUE, row.names = 1))
  all_pop_names <- rownames(full_matrix)
  np <- length(all_pop_names)
  if (k < 1 || k > np) stop("k must be between 1 and np")
  set.seed(seed)
  shared_pop_names <- if (randomize_shared) sample(all_pop_names, k) else all_pop_names[1:k]
  unique_pop_names <- setdiff(all_pop_names, shared_pop_names)
  matrix_A <- full_matrix
  matrix_B <- full_matrix
  noise_A <- matrix(rnorm(nrow(matrix_A) * ncol(matrix_A), mean = 0, sd = noise_sd), nrow = nrow(matrix_A))
  matrix_A <- abs(matrix_A + noise_A)
  noise_B <- matrix(rnorm(nrow(matrix_B) * ncol(matrix_B), mean = bias_mu, sd = noise_sd), nrow = nrow(matrix_B))
  matrix_B <- abs(matrix_B + noise_B)
  diag(matrix_A) <- 0; diag(matrix_B) <- 0
  label_A <- label_B <- rep("", np); names(label_A) <- names(label_B) <- all_pop_names
  for (pop in shared_pop_names) { label_A[pop] <- paste0("S_", pop); label_B[pop] <- paste0("S_", pop) }
  for (pop in unique_pop_names) { label_A[pop] <- paste0("A_", pop); label_B[pop] <- paste0("B_", pop) }
  rownames(matrix_A) <- label_A[rownames(matrix_A)]; colnames(matrix_A) <- label_A[colnames(matrix_A)]
  rownames(matrix_B) <- label_B[rownames(matrix_B)]; colnames(matrix_B) <- label_B[colnames(matrix_B)]
  list(A = matrix_A, B = matrix_B, np = np, k = k)
}

results <- data.frame()
for (k in k_values) {
  for (test_id in 1:num_tests) {
    test_seed <- seed_base + test_id * 100 + k
    matrices <- prepare_matrices_from_csv(input_file, k = k, seed = test_seed)
    A <- matrices$A; B <- matrices$B
    mds_results <- perform_mds(A, B)
    X <- mds_results$X; Y <- mds_results$Y; d_opt <- mds_results$d_opt
    pop_common <- intersect(rownames(A), rownames(B))
    X_sub <- X[pop_common, 1:d_opt]
    Y_sub <- Y[pop_common, 1:d_opt]
    Y_transformed <- apply_procrustes(X_sub, Y_sub, Y)
    B_calibrated <- coords_to_distances(Y_transformed)
    the_prior <- mean((A - B)^2)
    the_posterior <- mean((A - B_calibrated)^2)
    improvement <- (the_prior - the_posterior) / the_prior * 100
    results <- rbind(results, data.frame(k = k, test_id = test_id, the_prior = the_prior,
                                         the_posterior = the_posterior, improvement = improvement,
                                         data_source = "input"))
  }
}

summary_results <- aggregate(cbind(the_prior, the_posterior, improvement) ~ k, data = results, 
                             FUN = function(x) c(mean = mean(x), sd = sd(x), min = min(x), max = max(x)))
summary_results <- do.call(data.frame, summary_results)

# Plot
results_df <- summary_results
results_df$prop <- results_df$k / max(results_df$k)

p1 <- ggplot(results_df, aes(x = prop)) +
  geom_line(aes(y = the_prior.mean, color = "Prior"), linewidth = 1) +
  geom_point(aes(y = the_prior.mean, color = "Prior"), size = 2) +
  geom_ribbon(aes(ymin = the_prior.mean - the_prior.sd, ymax = the_prior.mean + the_prior.sd, fill = "Prior"), alpha = 0.2) +
  geom_line(aes(y = the_posterior.mean, color = "Posterior"), linewidth = 1) +
  geom_point(aes(y = the_posterior.mean, color = "Posterior"), size = 2) +
  geom_ribbon(aes(ymin = the_posterior.mean - the_posterior.sd, ymax = the_posterior.mean + the_posterior.sd, fill = "Posterior"), alpha = 0.2) +
  scale_x_continuous(labels = scales::percent) +
  labs(title = "the Prior vs Posterior", x = "Proportion of Shared Populations", y = "Differences") +
  scale_color_manual(values = c("Prior" = "red", "Posterior" = "blue")) +
  scale_fill_manual(values = c("Prior" = "red", "Posterior" = "blue")) +
  theme_minimal() + theme(legend.title = element_blank())

print(p1)
