# Local smoke test for DataFusionGDM
suppressPackageStartupMessages({
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos = "https://cloud.r-project.org")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "https://cloud.r-project.org")
  if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan", repos = "https://cloud.r-project.org")
  if (!requireNamespace("mice", quietly = TRUE)) install.packages("mice", repos = "https://cloud.r-project.org")
})

message("Loading package...")
devtools::load_all(".")

message("Running simulation...")
res <- simulate_genetic_distances(n_pops = 12, n_major_groups = 3, n_subgroups = 6,
                                  geo_dims = 2, genetic_dims = 2, use_noise = FALSE, verbose = FALSE)
stopifnot(is.matrix(res$distance_matrix), nrow(res$distance_matrix) == 12)

message("Testing plot creation...")
p <- create_mds_plot(res$distance_matrix, res$population_info)
invisible(ggplot2::ggplot_build(p))

message("Testing MDS + Procrustes...")
A <- res$distance_matrix
set.seed(123)
B <- A + matrix(rnorm(length(A), 0, 0.01), nrow = nrow(A))
diag(B) <- 0
B <- (B + t(B)) / 2
mds <- perform_mds(A, B)
rn <- rownames(A)
d <- max(1, min(ncol(mds$X), ncol(mds$Y)))
Xb <- mds$X[rn, 1:d, drop = FALSE]
Yb <- mds$Y[rn, 1:d, drop = FALSE]
Yt <- apply_procrustes(Xb, Yb, mds$Y[, 1:d, drop = FALSE])
Bcal <- coords_to_distances(Yt)
prior <- mean((A - B)^2)
post <- mean((A - Bcal)^2)
message(sprintf("MDS/Procrustes improvement: %.4f -> %.4f", prior, post))

message("Testing BESMI iterative imputation...")
set.seed(1)
mask <- matrix(FALSE, nrow = nrow(A), ncol = ncol(A))
idx <- sample(seq_len(nrow(A)), 3)
mask[idx, idx] <- TRUE
Min <- A
Min[mask] <- NA
imp <- besmi_iterative_imputation(Min, M_mask = mask, M_real = A, max_iterations = 2)
stopifnot(is.matrix(imp$final_matrix))

message("Smoke tests passed.")
