IWP_joint_prec_local <- getFromNamespace("IWP_joint_prec", "BayesGP")
local_poly_local <- getFromNamespace("local_poly", "BayesGP")
compute_weights_precision_local <- getFromNamespace("compute_weights_precision", "BayesGP")

build_exact_observation_matrix_test <- function(observed_x, support_x){
  n_obs <- length(observed_x)
  n_support <- length(support_x)
  observation_matrix <- matrix(0, nrow = n_obs, ncol = 2 * n_support)

  if(n_support == 0){
    return(observation_matrix)
  }

  matched_index <- match(observed_x, support_x)
  positive_rows <- which(!is.na(matched_index))

  if(length(positive_rows) > 0){
    observation_matrix[cbind(positive_rows, 2 * matched_index[positive_rows] - 1)] <- 1
  }

  observation_matrix
}

build_exact_iwp_prior_covariance <- function(x){
  shifted_x <- sort(x - min(x))
  support_x <- shifted_x[shifted_x > 0]

  if(length(support_x) == 0){
    return(matrix(0, nrow = length(shifted_x), ncol = length(shifted_x)))
  }

  exact_B <- build_exact_observation_matrix_test(
    observed_x = shifted_x,
    support_x = support_x
  )
  exact_P <- as.matrix(IWP_joint_prec_local(t_vec = support_x, p = 2))

  as.matrix(exact_B %*% solve(exact_P, t(exact_B)))
}

build_ospline_iwp_prior_covariance <- function(x, k){
  shifted_x <- sort(x - min(x))
  knot_sequence <- seq(0, max(shifted_x), length.out = k)
  temp_data <- data.frame(z = shifted_x)

  temp_instance <- methods::new(
    "iwp",
    response_var = as.name("y"),
    smoothing_var = as.name("z"),
    order = 2,
    knots = knot_sequence,
    k = k,
    observed_x = shifted_x,
    sd.prior = list(),
    psd.prior = list(),
    boundary.prior = list(),
    data = temp_data,
    X = matrix(0, nrow = 0, ncol = 0),
    B = matrix(0, nrow = 0, ncol = 0),
    P = matrix(0, nrow = 0, ncol = 0),
    initial_location = 0,
    region = range(c(0, shifted_x))
  )

  ospline_B <- as.matrix(local_poly_local(temp_instance))
  ospline_P <- as.matrix(compute_weights_precision_local(temp_instance))

  as.matrix(ospline_B %*% solve(ospline_P, t(ospline_B)))
}

frobenius_relative_error <- function(approximate_cov, exact_cov){
  sqrt(sum((approximate_cov - exact_cov)^2) / sum(exact_cov^2))
}

max_relative_error <- function(approximate_cov, exact_cov){
  max(abs(approximate_cov - exact_cov)) / max(abs(exact_cov))
}

test_that("native O-spline IWP covariance approaches the exact state-space covariance", {
  evaluation_grids <- list(
    balanced = c(0.10, 0.35, 0.75, 1.30, 2.05, 2.95, 4.00, 5.20, 6.60, 8.15),
    clustered = c(0.20, 0.24, 0.31, 0.52, 1.00, 1.80, 3.00, 4.70, 7.00, 8.50)
  )

  for(x in evaluation_grids){
    exact_cov <- build_exact_iwp_prior_covariance(x)
    coarse_cov <- build_ospline_iwp_prior_covariance(x, k = 20)
    refined_cov <- build_ospline_iwp_prior_covariance(x, k = 160)

    coarse_frob <- frobenius_relative_error(coarse_cov, exact_cov)
    refined_frob <- frobenius_relative_error(refined_cov, exact_cov)
    refined_max <- max_relative_error(refined_cov, exact_cov)

    expect_lt(refined_frob, coarse_frob)
    expect_lt(refined_frob, 2e-5)
    expect_lt(refined_max, 1.2e-5)
  }
})
