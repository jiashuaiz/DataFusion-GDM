# Example: BESMI end-to-end (concise)

cat("[DataFusionGDM] Preparing input...\n")
library(DataFusionGDM)

if (!file.exists("data/GDM_simulated.csv")) {
  dir.create("data", showWarnings = FALSE)
  export_simulated_gdm("data/GDM_simulated.csv", scenario = "default", n_pops = 40)
}

full <- besmi_prepare_full_dataset("data/GDM_simulated.csv")

cat("[DataFusionGDM] Creating masked datasets...\n")
summary_table <- data.frame()
for (k in 1:2) {
  for (bs_i in 1:2) {
    res <- besmi_create_masked_matrices(full, k, seed = 100 + 10 * k + bs_i)
    dir.create("data/training_set", recursive = TRUE, showWarnings = FALSE)
    dir.create("data/position_set", recursive = TRUE, showWarnings = FALSE)
    tr <- sprintf("data/training_set/masked_k%d_bs%d.rds", k, bs_i)
    mp <- sprintf("data/position_set/mask_k%d_bs%d.rds", k, bs_i)
    saveRDS(res$masked_matrix, tr)
    saveRDS(res$mask_position, mp)
    summary_table <- rbind(summary_table, data.frame(k = k, bs = bs_i, training_path = tr, mask_path = mp))
  }
}

cat("[DataFusionGDM] Running batch imputation...\n")
paths <- summary_table$training_path
metrics <- besmi_batch_impute(paths, the_method = "midastouch", max_iter = 3,
  imputation_convergence_threshold = 1e-3, propagation_convergence_threshold = 1e-3,
  distance_metric = "mae", output_dir = "data/imputation_set")

cat("[DataFusionGDM] Preview metrics:\n")
print(head(metrics))

cat("[DataFusionGDM] Done.\n")


