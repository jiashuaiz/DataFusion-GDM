# Example: MDS + Procrustes calibration (concise)

cat("[DataFusionGDM] Loading example distance matrix...\n")
library(DataFusionGDM)

pkg_file <- system.file("extdata", "GDM_simulated.csv", package = "DataFusionGDM")
if (!nzchar(pkg_file)) {
  export_simulated_gdm("GDM_simulated.csv", scenario = "default", n_pops = 40)
  pkg_file <- "GDM_simulated.csv"
}

mat <- as.matrix(read.csv(pkg_file, row.names = 1, check.names = FALSE))
A <- mat; B <- mat * 1.05

cat("[DataFusionGDM] Performing MDS and Procrustes...\n")
mds <- perform_mds(A, B)
Yt <- apply_procrustes(mds$X, mds$Y, mds$Y)
B2 <- coords_to_distances(Yt)

rmse_before <- sqrt(mean((A - B)^2))
rmse_after  <- sqrt(mean((A - B2)^2))
cat(sprintf("[DataFusionGDM] RMSE before: %.4f | after: %.4f\n", rmse_before, rmse_after))

cat("[DataFusionGDM] Plotting absolute differences (before vs after)...\n")
op <- par(mfrow = c(1, 2))
heatmap(abs(A - B), Rowv = NA, Colv = NA, symm = TRUE, main = "Abs diff (before)")
heatmap(abs(A - B2), Rowv = NA, Colv = NA, symm = TRUE, main = "Abs diff (after)")
par(op)

cat("[DataFusionGDM] Done.\n")


