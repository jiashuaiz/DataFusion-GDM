## Example: BESMI prepare datasets and batch imputation (with quick visualization)

library(DataFusionGDM)
library(tidyverse)

# Ensure input GDM exists
if (!file.exists("data/GDM_simulated.csv")) {
  dir.create("data", showWarnings = FALSE)
  export_simulated_gdm("data/GDM_simulated.csv", scenario = "default", n_pops = 40)
}

full_matrix <- besmi_prepare_full_dataset("data/GDM_simulated.csv")

# Prepare masked datasets for k = 1..5
k_values <- 1:5
summary_table <- data.frame()
for (k in k_values) {
  n_samples <- if (k <= 3) 15 else if (k <= 6) 8 else 4
  for (bs_i in 1:n_samples) {
    res <- besmi_create_masked_matrices(full_matrix, k, seed = 233 + bs_i + 100 * k)
    dir.create("data/training_set", recursive = TRUE, showWarnings = FALSE)
    dir.create("data/position_set", recursive = TRUE, showWarnings = FALSE)
    training_path <- sprintf("data/training_set/masked_k%d_bs%d.rds", k, bs_i)
    mask_position_path <- sprintf("data/position_set/mask_k%d_bs%d.rds", k, bs_i)
    saveRDS(res$masked_matrix, training_path)
    saveRDS(res$mask_position, mask_position_path)
    summary_table <- rbind(summary_table, data.frame(
      k_value = k, bootstrap_sample_i = bs_i,
      group_u = paste(res$group_u, collapse = ";"),
      group_s = paste(res$group_s, collapse = ";"),
      masked_percentage = res$masked_percentage,
      training_path = training_path,
      mask_position_path = mask_position_path,
      stringsAsFactors = FALSE
    ))
  }
}
dir.create("data/summary_tables", recursive = TRUE, showWarnings = FALSE)
write.csv(summary_table, "data/summary_tables/summary_table_datasets.csv", row.names = FALSE)

# Batch process imputation (example method)
dataset_paths <- summary_table$training_path
metrics <- besmi_batch_impute(dataset_paths = dataset_paths, the_method = "midastouch", max_iter = 5,
                              imputation_convergence_threshold = 1e-3,
                              propagation_convergence_threshold = 1e-3,
                              distance_metric = "mae",
                              output_dir = "data/imputation_set")

write.csv(metrics, "data/metrics_df_midastouch.csv", row.names = FALSE)

# Quick visualization: imputation distance by iteration
if (nrow(metrics) > 0) {
  library(ggplot2)
  gg <- ggplot(metrics, aes(x = iteration, y = imputation_dis, group = interaction(k, bs))) +
    geom_line(alpha = 0.3) +
    geom_smooth(se = FALSE, color = "blue") +
    labs(title = "Imputation distance by iteration", x = "Iteration", y = "Distance") +
    theme_minimal()
  print(gg)
}
