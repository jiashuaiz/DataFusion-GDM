# End-to-end functional test for DataFusionGDM
suppressPackageStartupMessages({
  library(DataFusionGDM)
})

message("Running scenario (default)...")
ans <- run_genetic_scenario("default", n_pops = 20, verbose = FALSE)
stopifnot(is.list(ans), !is.null(ans$results), !is.null(ans$plots))

# Validate distance matrix
DM <- ans$results$distance_matrix
if (!(is.matrix(DM) || (is.numeric(DM) && length(dim(DM)) == 2))) stop("Distance matrix not a matrix")

# Plot to PDF device to avoid interactive dependence
pdf(file = tempfile(fileext = ".pdf"))
ans$plots$mds()
ans$plots$heatmap()
dev.off()

message("Testing export_simulated_gdm and re-loading...")
ret <- export_simulated_gdm(file.path(tempdir(), "GDM_sim_default.csv"), scenario = "simple", n_pops = 10, verbose = FALSE)
full <- besmi_prepare_full_dataset(ret)
stopifnot(is.matrix(full), nrow(full) == 10, ncol(full) == 10)

message("Testing BESMI masking + iterative imputation...")
mk <- besmi_create_masked_matrices(full, k = 3, seed = 1)
imp_iter <- besmi_iterative_imputation(mk$masked_matrix, mk$mask_position, M_real = full, max_iterations = 2)
stopifnot(is.matrix(imp_iter$final_matrix))

message("Testing BESMI KNN imputation (VIM)...")
imp_knn <- besmi_knn_impute(mk$masked_matrix, mk$mask_position, M_real = full)
stopifnot(is.list(imp_knn))

message("Testing BESMI batch (single dataset)...")
base_dir <- file.path(tempdir(), "besmi_batch")
tr_dir <- file.path(base_dir, "data", "training_set")
pos_dir <- file.path(base_dir, "data", "position_set")
imp_dir <- file.path(base_dir, "data", "imputation_set")
dir.create(tr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pos_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(mk$masked_matrix, file.path(tr_dir, "masked_k3_bs1.rds"))
saveRDS(mk$mask_position, file.path(pos_dir, "mask_k3_bs1.rds"))

paths <- list.files(tr_dir, pattern = "\\.rds$", full.names = TRUE)
metrics <- besmi_batch_impute(paths, the_method = "KNN", max_iter = 1, output_dir = imp_dir)
stopifnot(is.data.frame(metrics), nrow(metrics) >= 1)

message("All end-to-end tests passed.")
