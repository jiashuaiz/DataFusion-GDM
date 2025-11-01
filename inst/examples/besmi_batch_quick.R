library(DataFusionGDM)

full <- simulate_genetic_distances(n_pops = 40, verbose = FALSE, seed = 233)$distance_matrix

summary_table <- data.frame()
for (k in 1:2) {
  for (bs_i in 1:2) {
    res <- besmi_create_masked_matrices(full, k, seed = 100 + 10 * k + bs_i)
    base_dir <- tempdir()
    dir.create(file.path(base_dir, "training_set"), recursive = TRUE, showWarnings = FALSE)
    training_path <- file.path(base_dir, sprintf("training_set/masked_k%d_bs%d.rds", k, bs_i))
    saveRDS(res$masked_matrix, training_path)
    summary_table <- rbind(summary_table, data.frame(k = k, bs = bs_i, training_path = training_path))
  }
}

paths <- summary_table$training_path
output_dir <- file.path(tempdir(), "imputation_set")
metrics <- besmi_batch_impute(paths, the_method = "midastouch", max_iter = 2,
  imputation_convergence_threshold = 1e-3, propagation_convergence_threshold = 1e-3,
  distance_metric = "mae", output_dir = output_dir)
print(head(metrics))

