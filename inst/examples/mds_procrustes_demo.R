library(DataFusionGDM)

set.seed(42)
res <- simulate_genetic_distances(n_pops = 30, verbose = FALSE)
G <- res$distance_matrix
A <- G + matrix(rnorm(length(G), 0, 0.02), nrow = nrow(G)); diag(A) <- 0
B <- G + matrix(rnorm(length(G), 0.03, 0.02), nrow = nrow(G)); diag(B) <- 0

mds <- perform_mds(A, B)
Yt <- apply_procrustes(mds$X, mds$Y, mds$Y)
B_cal <- coords_to_distances(Yt)

cat(sprintf("MSE prior: %.4f\n", mean((A - B)^2)))
cat(sprintf("MSE post:  %.4f\n", mean((A - B_cal)^2)))

