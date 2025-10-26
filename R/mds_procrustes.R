# MDS and Procrustes alignment utilities

#' Double-center a distance matrix
#' @keywords internal
.double_center <- function(D) {
  n <- nrow(D)
  J <- diag(n) - (1/n) * matrix(1, n, n)
  B <- -0.5 * J %*% (D^2) %*% J
  return(B)
}

#' Perform MDS on a pair of distance matrices
#' @param A First distance matrix
#' @param B Second distance matrix
#' @return A list with coordinates `X`, `Y`, optimal dimension `d_opt`, and variance info
#' @export
perform_mds <- function(A, B) {
  A_centered <- .double_center(A)
  B_centered <- .double_center(B)
  eigen_A <- eigen(A_centered, symmetric = TRUE)
  eigen_B <- eigen(B_centered, symmetric = TRUE)
  # Count numerically positive eigenvalues
  pos_eigen_A <- sum(eigen_A$values > 1e-12)
  pos_eigen_B <- sum(eigen_B$values > 1e-12)
  d_max_pos <- max(1, min(pos_eigen_A, pos_eigen_B))
  # Proportion of variance for positive eigen subspace
  pos_vals_A <- eigen_A$values[eigen_A$values > 1e-12]
  pos_vals_B <- eigen_B$values[eigen_B$values > 1e-12]
  denom_A <- sum(pos_vals_A)
  denom_B <- sum(pos_vals_B)
  var_explained_A <- if (d_max_pos > 0 && denom_A > 0) cumsum(eigen_A$values[1:d_max_pos]) / denom_A else rep(0, d_max_pos)
  var_explained_B <- if (d_max_pos > 0 && denom_B > 0) cumsum(eigen_B$values[1:d_max_pos]) / denom_B else rep(0, d_max_pos)
  # Find dimensions reaching ~100% with tolerance
  idx_A <- which(var_explained_A >= 0.999)
  idx_B <- which(var_explained_B >= 0.999)
  d_100_A <- if (length(idx_A)) min(idx_A) else d_max_pos
  d_100_B <- if (length(idx_B)) min(idx_B) else d_max_pos
  # Choose optimal dimensionality
  d_opt <- max(1, min(d_100_A, d_100_B, d_max_pos))
  # Calculate MDS coordinates
  X <- eigen_A$vectors[, 1:d_opt, drop = FALSE] %*% diag(sqrt(pmax(eigen_A$values[1:d_opt], 0)))
  Y <- eigen_B$vectors[, 1:d_opt, drop = FALSE] %*% diag(sqrt(pmax(eigen_B$values[1:d_opt], 0)))
  rownames(X) <- rownames(A)
  rownames(Y) <- rownames(B)
  list(
    X = X,
    Y = Y,
    d_opt = d_opt,
    variance_info = list(
      eigenvalues_A = eigen_A$values,
      eigenvalues_B = eigen_B$values,
      var_explained_A = var_explained_A,
      var_explained_B = var_explained_B
    )
  )
}

#' Procrustes alignment and mapping back to distances
#' @param X_base Base coordinates for target alignment
#' @param Y_base Base coordinates for source alignment
#' @param Y Full source coordinates to transform
#' @return Transformed coordinates matrix
#' @export
apply_procrustes <- function(X_base, Y_base, Y) {
  proc <- vegan::procrustes(X_base, Y_base, scale = TRUE, symmetric = FALSE)
  translation <- proc$translation
  dilation <- proc$scale
  rotation <- proc$rotation
  Y_transformed <- (Y %*% rotation) * dilation + 
    matrix(translation, nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
  return(Y_transformed)
}

#' Convert coordinate matrix to distance matrix
#' @param coords Numeric coordinate matrix
#' @return Symmetric distance matrix
#' @export
coords_to_distances <- function(coords) {
  dist_M <- as.matrix(stats::dist(coords))
  rownames(dist_M) = colnames(dist_M) <- rownames(coords)
  return(dist_M)
}
