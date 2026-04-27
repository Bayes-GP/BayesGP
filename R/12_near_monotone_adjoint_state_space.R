adj_SS_cov_mat_c2 <- function(s, t, c = 1){
  M <- matrix(nrow = 2, ncol = 2)
  M[1, 1] <- 1 / 12 * (
    4 * c * s^3 - 3 * s^4 - 4 * c * t^3 + t^4 +
      6 * (2 * c * s - s^2) * t^2 -
      4 * (3 * c * s^2 - 2 * s^3) * t
  ) / (c - t)
  M[1, 2] <- 1 / 2 * (
    c * s^2 - s^3 + (c - s) * t^2 -
      2 * (c * s - s^2) * t
  ) / (c - t)
  M[2, 2] <- (c * s - s^2 - (c - s) * t) / (c - t)
  Matrix::forceSymmetric(M)
}

adj_SS_cov_mat_c5 <- function(s, t, c = 1){
  M <- matrix(nrow = 2, ncol = 2)
  M[1, 1] <- 2 / 15 * (
    24 * c^4 - 12 * c^3 * s - 3 * c^2 * s^2 + c * s^3 -
      10 * c * t^3 + 15 * (3 * c^2 - c * s) * t^2 -
      24 * (c^3 - 2 * c^2 * t + c * t^2) * sqrt(c - s) * sqrt(c - t) -
      30 * (2 * c^3 - c^2 * s) * t
  ) / c
  M[1, 2] <- 1 / 5 * (
    (c^2 - 2 * c * s + s^2) * sqrt(c - s) * sqrt(c) -
      5 * (c^2 - 2 * c * t + t^2) * sqrt(c - s) * sqrt(c) +
      4 * (c^2 - 2 * c * t + t^2) * sqrt(c - t) * sqrt(c)
  ) / (sqrt(c - s) * sqrt(c))
  M[2, 2] <- 1 / 2 * (2 * c * s - s^2 - 2 * c * t + t^2) / (c - s)
  Matrix::forceSymmetric(M)
}

adj_R_trans_matrix_c1 <- function(s, x, c = 1, alpha = 2){
  R <- matrix(nrow = 2, ncol = 2)
  R[1, 1] <- 1
  R[2, 1] <- 0
  R[1, 2] <- (
    alpha * (c - s)^(1 / alpha) * c -
      alpha * (c - s)^(1 / alpha) * x -
      (alpha * c - alpha * s) * (c - x)^(1 / alpha)
  ) / ((alpha - 1) * (c - s)^(1 / alpha))
  R[2, 2] <- (c - x)^(1 / alpha) / (c - s)^(1 / alpha)
  R
}

adj_R_trans_matrix_c2 <- function(s, x, c = 1){
  R <- matrix(nrow = 2, ncol = 2)
  R[1, 1] <- 1
  R[2, 1] <- 0
  R[1, 2] <- 1 / 2 * (2 * c * s - s^2 - 2 * c * x + x^2) / (c - x)
  R[2, 2] <- (c - s) / (c - x)
  R
}

adj_SS_cov_mat <- function(s, t, c = 1, alpha = 2){
  curvature <- resolve_curvature_parameter(a = alpha)

  if(isTRUE(all.equal(curvature, 1))){
    return(adj_SS_cov_mat_c2(s = s, t = t, c = c))
  }
  if(isTRUE(all.equal(curvature, -1))){
    return(adj_SS_cov_mat_c2(s = s, t = t, c = c))
  }
  if(isTRUE(all.equal(curvature, 2))){
    return(adj_SS_cov_mat_c5(s = s, t = t, c = c))
  }

  stop(
    "Reverse exact mGP state-space currently supports non-left references ",
    "only for `a` equal to 1, 2, or -1. Use `computation = \"fem\"` for ",
    "other curvatures."
  )
}

adj_R_trans_matrix <- function(s, x, c = 1, alpha = 2){
  curvature <- resolve_curvature_parameter(a = alpha)

  if(isTRUE(all.equal(curvature, 1))){
    return(adj_R_trans_matrix_c2(s = s, x = x, c = c))
  }
  if(isTRUE(all.equal(curvature, -1))){
    return(adj_R_trans_matrix_c2(s = s, x = x, c = c))
  }
  if(isTRUE(all.equal(curvature, 2))){
    return(adj_R_trans_matrix_c1(s = s, x = x, c = c, alpha = 2))
  }

  stop(
    "Reverse exact mGP state-space currently supports non-left references ",
    "only for `a` equal to 1, 2, or -1. Use `computation = \"fem\"` for ",
    "other curvatures."
  )
}

adj_mGP_joint_prec <- function(t_vec, a = NULL, alpha = NULL, c){
  curvature <- resolve_curvature_parameter(a = a, alpha = alpha)
  t_vec <- sort(unique(as.numeric(t_vec)))

  if(length(t_vec) == 0){
    return(Matrix::Matrix(0, nrow = 0, ncol = 0, sparse = TRUE))
  }
  if(any(!is.finite(t_vec)) || any(t_vec <= 0)){
    stop("Reverse exact mGP support must contain finite positive distances.")
  }
  if(any(t_vec >= c)){
    stop("Reverse exact mGP support requires all distances to be smaller than `c`.")
  }

  compute_ci <- function(grid, i){
    as.matrix(Matrix::forceSymmetric(solve(adj_SS_cov_mat(
      s = grid[i + 1],
      t = grid[i],
      alpha = curvature,
      c = c
    ))))
  }
  compute_ai <- function(grid, i, ci){
    transition <- adj_R_trans_matrix(
      x = grid[i],
      s = grid[i + 1],
      alpha = curvature,
      c = c
    )
    t(transition) %*% ci %*% transition
  }
  compute_bi <- function(grid, i, ci){
    transition <- adj_R_trans_matrix(
      x = grid[i],
      s = grid[i + 1],
      alpha = curvature,
      c = c
    )
    -t(transition) %*% ci
  }

  n <- length(t_vec)
  if(n == 1){
    return(as(as.matrix(Matrix::forceSymmetric(compute_ci(c(0, t_vec[1]), 1))), "dgCMatrix"))
  }

  block_b <- list()
  block_ac <- list()
  block_c <- list()

  for(i in seq_len(n - 1)){
    block_c[[i]] <- compute_ci(t_vec, i)
  }
  for(i in 2:(n - 1)){
    block_ac[[i]] <- compute_ai(t_vec, i, block_c[[i]]) + block_c[[i - 1]]
  }
  block_ac[[1]] <- compute_ai(t_vec, 1, block_c[[1]]) + compute_ci(c(0, t_vec[1]), 1)
  block_ac[[n]] <- compute_ci(t_vec, n - 1)

  for(i in seq_len(n - 1)){
    block_b[[i]] <- compute_bi(t_vec, i, block_c[[i]])
  }

  precision <- matrix(0, nrow = 0, ncol = 2 * n)
  for(i in seq_len(n - 1)){
    precision <- rbind(
      precision,
      cbind(
        matrix(0, nrow = 2, ncol = 2 * (i - 1)),
        block_ac[[i]],
        block_b[[i]],
        matrix(0, nrow = 2, ncol = 2 * (n - i - 1))
      )
    )
  }
  precision <- rbind(
    precision,
    cbind(matrix(0, nrow = 2, ncol = 2 * (n - 1)), block_ac[[n]])
  )

  as(as.matrix(Matrix::forceSymmetric(precision)), "dgCMatrix")
}
