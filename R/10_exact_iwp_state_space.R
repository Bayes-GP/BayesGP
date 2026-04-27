build_iwp_state_space_instance <- function(rand_effect, response_var, data, envir = parent.frame()){
  smoothing_var <- extract_term_argument(rand_effect, "smoothing_var", envir = envir)
  if(is.null(smoothing_var)){
    smoothing_var <- extract_term_argument(rand_effect, "x", envir = envir)
  }
  if(is.null(smoothing_var)){
    smoothing_var <- toString(rand_effect[[2]])
  }

  if(is.null(data[[smoothing_var]])){
    stop("Cannot find smoothing variable `", smoothing_var, "` in `data`.")
  }

  order <- extract_term_argument(rand_effect, "order", envir = envir, default = 2)
  if(length(order) != 1 || order != 2){
    stop("Exact state-space `iwp` terms currently require `order = 2`.")
  }

  observed_x <- data[[smoothing_var]]
  grid <- extract_term_argument(rand_effect, "grid", envir = envir, default = NULL)
  domain_x <- if(is.null(grid)) observed_x else c(observed_x, grid)

  initial_location <- resolve_left_boundary(
    x = domain_x,
    initial_location = extract_term_argument(rand_effect, "initial_location", envir = envir, default = "left")
  )

  shifted_observed_x <- observed_x - initial_location
  shifted_domain_x <- domain_x - initial_location

  sd.prior <- normalize_sd_prior(extract_term_argument(rand_effect, "sd.prior", envir = envir))
  boundary.prior <- normalize_boundary_prior(
    extract_term_argument(rand_effect, "boundary.prior", envir = envir),
    dimension = 1
  )

  support_grid <- sort(unique(shifted_domain_x[shifted_domain_x > 0]))
  B <- build_exact_observation_matrix(
    observed_x = shifted_observed_x,
    support_x = support_grid
  )
  P <- build_exact_precision_with_diagnostics(
    support_grid = support_grid,
    model_name = "iwp",
    precision_builder = function(grid) IWP_joint_prec(t_vec = grid, p = order)
  )

  instance <- methods::new(
    "iwp",
    response_var = as.name(response_var),
    smoothing_var = as.name(smoothing_var),
    order = order,
    knots = numeric(),
    k = length(support_grid),
    observed_x = sort(shifted_observed_x),
    sd.prior = sd.prior,
    psd.prior = list(),
    boundary.prior = boundary.prior,
    data = data,
    X = matrix(shifted_observed_x, ncol = 1),
    B = as.matrix(B),
    P = P,
    initial_location = initial_location,
    region = range(domain_x)
  )

  list(
    instance = instance,
    meta = list(
      component = smoothing_var,
      model_name = "iwp",
      computation_method = "state-space",
      support_grid = support_grid,
      order = order
    )
  )
}

exact_iwp_samples <- function(object, component, refined_x = NULL,
                              include_intercept = FALSE, deriv = 0){
  instance <- NULL
  metadata <- object$exact_iwp_meta[[component]]

  if(is.null(metadata)){
    stop("Unknown component `", component, "`.")
  }

  for(candidate in object$instances){
    if(as.character(candidate@smoothing_var) == component){
      instance <- candidate
      break
    }
  }

  if(is.null(instance)){
    stop("Unable to locate the stored exact `iwp` instance for `", component, "`.")
  }

  if(!deriv %in% c(0, 1)){
    stop("Exact state-space `iwp` prediction currently supports only `deriv = 0` or `deriv = 1`.")
  }

  coef_samps <- extract_random_samples(object, component = component)
  boundary_samps <- extract_boundary_samples(object, component = component)
  intercept_samps <- if(include_intercept && deriv == 0) {
    extract_intercept_samples(object)
  } else {
    matrix(0, nrow = 1, ncol = ncol(coef_samps))
  }

  available_x <- c(instance@initial_location, instance@initial_location + metadata$support_grid)
  eval_x <- if(is.null(refined_x)) available_x else refined_x
  available_lookup <- match(eval_x, available_x)

  if(anyNA(available_lookup)){
    stop("Exact state-space `iwp` evaluation currently requires `refined_x` to be included in the term support. Supply `grid = ...` in `f()` when fitting if extra points are needed.")
  }

  fitted_draws <- matrix(0, nrow = length(available_x), ncol = ncol(coef_samps))
  if(length(metadata$support_grid) > 0){
    fitted_draws[-1, ] <- coef_samps[seq(deriv + 1, 2 * length(metadata$support_grid), by = 2), , drop = FALSE]
  }

  if(nrow(boundary_samps) > 0){
    boundary_eval <- if(deriv == 0) {
      matrix(eval_x - instance@initial_location, ncol = 1)
    } else {
      matrix(1, nrow = length(eval_x), ncol = 1)
    }
    fitted_draws <- fitted_draws[available_lookup, , drop = FALSE] + boundary_eval %*% boundary_samps
  } else {
    fitted_draws <- fitted_draws[available_lookup, , drop = FALSE]
  }

  if(deriv == 0){
    fitted_draws <- fitted_draws + matrix(rep(intercept_samps, each = nrow(fitted_draws)), nrow = nrow(fitted_draws))
  }

  data.frame(x = eval_x, as.data.frame(as.matrix(fitted_draws)))
}
