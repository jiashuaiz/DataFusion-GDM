library(DataFusionGDM)

res <- simulate_genetic_distances(n_pops = 30, verbose = FALSE, seed = 42)
G <- res$distance_matrix
idx <- seq_len(length(G))
A_noise <- matrix(sin(idx) * 0.02, nrow = nrow(G))
B_noise <- matrix(cos(idx) * 0.02 + 0.03, nrow = nrow(G))
A <- G + A_noise; diag(A) <- 0
B <- G + B_noise; diag(B) <- 0

mds <- perform_mds(A, B)
Yt <- apply_procrustes(mds$X, mds$Y, mds$Y)
B_cal <- coords_to_distances(Yt)

cat(sprintf("MSE prior: %.4f\n", mean((A - B)^2)))
cat(sprintf("MSE post:  %.4f\n", mean((A - B_cal)^2)))

