test_that("package loads", {
  expect_silent(library(DataFusionGDM))
})

test_that("simulate_genetic_distances returns consistent structure", {
  res <- simulate_genetic_distances(n_pops = 12, n_major_groups = 3, n_subgroups = 6,
                                    geo_dims = 2, genetic_dims = 2, use_noise = FALSE, verbose = FALSE)
  expect_true(is.matrix(res$distance_matrix))
  expect_equal(nrow(res$distance_matrix), 12)
  expect_equal(ncol(res$distance_matrix), 12)
})

test_that("MDS + Procrustes runs on small matrices", {
  set.seed(1)
  res <- simulate_genetic_distances(n_pops = 10, n_major_groups = 3, n_subgroups = 6,
                                    geo_dims = 2, genetic_dims = 2, use_noise = FALSE, verbose = FALSE)
  A <- res$distance_matrix
  B <- (A + matrix(rnorm(length(A), 0, 0.01), nrow(A)))
  diag(B) <- 0
  B <- (B + t(B)) / 2
  mds <- perform_mds(A, B)
  rn <- rownames(A)
  d <- max(1, min(ncol(mds$X), ncol(mds$Y)))
  Xb <- mds$X[rn, 1:d, drop = FALSE]
  Yb <- mds$Y[rn, 1:d, drop = FALSE]
  Yt <- apply_procrustes(Xb, Yb, mds$Y[, 1:d, drop = FALSE])
  Bcal <- coords_to_distances(Yt)
  expect_true(is.matrix(Bcal))
  expect_equal(dim(Bcal), dim(A))
})

test_that("BESMI iterative imputation runs with tiny settings", {
  set.seed(2)
  res <- simulate_genetic_distances(n_pops = 10, n_major_groups = 3, n_subgroups = 6,
                                    geo_dims = 2, genetic_dims = 2, use_noise = FALSE, verbose = FALSE)
  A <- res$distance_matrix
  mask <- matrix(FALSE, nrow = nrow(A), ncol = ncol(A))
  idx <- sample(seq_len(nrow(A)), 2)
  mask[idx, idx] <- TRUE
  Min <- A
  Min[mask] <- NA
  imp <- besmi_iterative_imputation(Min, M_mask = mask, M_real = A, max_iterations = 1)
  expect_true(is.matrix(imp$final_matrix))
})


