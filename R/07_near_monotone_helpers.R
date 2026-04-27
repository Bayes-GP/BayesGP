normalize_sd_prior <- function(prior){
  if(is.null(prior)){
    return(list(prior = "exp", param = list(u = 1, alpha = 0.5)))
  }

  if(is.numeric(prior) && length(prior) == 1){
    return(list(prior = "exp", param = list(u = as.numeric(prior), alpha = 0.5)))
  }

  if(!is.list(prior)){
    stop("`sd.prior` must be NULL, a scalar, or a list.")
  }

  if(is.null(prior$prior)){
    prior$prior <- "exp"
  }

  if(prior$prior != "exp"){
    stop("BayesGP currently supports only exponential PC priors for `sd.prior`.")
  }

  if(is.null(prior$param)){
    stop("`sd.prior` lists must include a `param` entry.")
  }

  if(is.numeric(prior$param) && length(prior$param) == 1){
    prior$param <- list(u = as.numeric(prior$param), alpha = 0.5)
  }

  if(is.null(prior$param$u)){
    stop("`sd.prior$param$u` must be supplied.")
  }

  if(is.null(prior$param$alpha)){
    prior$param$alpha <- 0.5
  }

  prior
}

resolve_step_prior_value <- function(prior, label = "`sd.prior`"){
  has_h <- !is.null(prior$h)
  has_step <- !is.null(prior$step)

  if(!has_h && !has_step){
    return(NULL)
  }

  if(has_h && has_step && !isTRUE(all.equal(as.numeric(prior$h), as.numeric(prior$step)))){
    stop(label, " should supply at most one of `h` and `step`, or keep them identical.")
  }

  step_value <- if(has_h) prior$h else prior$step

  if(!is.numeric(step_value) || length(step_value) != 1 || !is.finite(step_value) || step_value <= 0){
    stop(label, "`h` / `step` must be a single positive number.")
  }

  as.numeric(step_value)
}

convert_nearmono_sd_prior <- function(sd.prior, model_name, curvature,
                                      c_original, c_shift,
                                      initial_location){
  step_value <- resolve_step_prior_value(sd.prior)

  if(is.null(step_value)){
    return(list(sd.prior = sd.prior, psd.prior = list()))
  }

  x_value <- sd.prior$x
  if(is.null(x_value)){
    stop(
      "Near-monotone PSD priors require `sd.prior$x` when `sd.prior$h` or ",
      "`sd.prior$step` is supplied."
    )
  }

  if(!is.numeric(x_value) || length(x_value) != 1 || !is.finite(x_value)){
    stop("`sd.prior$x` must be a single finite number for near-monotone PSD priors.")
  }
  if(x_value + c_original <= 0 || x_value + step_value + c_original <= 0){
    stop("Near-monotone PSD priors require `sd.prior$x + c` and `sd.prior$x + h + c` to be positive on the original scale.")
  }
  centered_x_value <- as.numeric(x_value) - initial_location

  psd.prior <- sd.prior
  converted_prior <- sd.prior
  converted_prior$h <- NULL
  converted_prior$step <- step_value
  converted_prior$x <- as.numeric(x_value)
  converted_prior$centered_x <- centered_x_value

  converted_prior$param <- switch(
    model_name,
    mgp = prior_conversion_mgp(
      prior = sd.prior$param,
      h = step_value,
      c = c_shift,
      a = curvature,
      x = centered_x_value
    ),
    tiwp2 = prior_conversion_tiwp2(
      prior = sd.prior$param,
      h = step_value,
      c = c_shift,
      a = curvature,
      x = centered_x_value,
      normalized = TRUE
    ),
    stop("Near-monotone PSD prior conversion is available only for `mgp` and `tiwp2`.")
  )

  list(sd.prior = converted_prior, psd.prior = psd.prior)
}

normalize_boundary_prior <- function(prior, dimension = 1){
  if(is.null(prior)){
    prior <- list(mean = 0, prec = 0.001)
  }

  if(is.null(prior$mean)){
    prior$mean <- 0
  }

  if(is.null(prior$prec)){
    prior$prec <- 0.001
  }

  prior$mean <- rep(prior$mean[1], length.out = dimension)
  prior$prec <- rep(prior$prec[1], length.out = dimension)
  prior
}

normalize_family_control <- function(control.family, family){
  family <- tolower(family)

  if(missing(control.family) || is.null(control.family)){
    control.family <- list()
  }

  if(is.numeric(control.family) && length(control.family) == 1){
    control.family <- list(sd.prior = control.family)
  }

  if(family == "gaussian"){
    sd_known <- NULL

    if(!is.null(control.family$sd)){
      if(is.numeric(control.family$sd) && length(control.family$sd) == 1){
        sd_known <- as.numeric(control.family$sd)
      } else {
        stop("`control.family$sd` must be a scalar when the Gaussian SD is fixed.")
      }
    }

    if(is.null(sd_known)){
      if(is.null(control.family$sd.prior) && !is.null(control.family$prior)){
        control.family$sd.prior <- control.family$prior
      }
      control.family$sd.prior <- normalize_sd_prior(control.family$sd.prior)
      control.family$sd_known <- FALSE
      control.family$sd <- NULL
    } else {
      if(sd_known <= 0){
        stop("`control.family$sd` must be strictly positive.")
      }
      control.family$sd.prior <- NULL
      control.family$sd_known <- TRUE
      control.family$sd <- sd_known
    }
  }

  control.family
}

normalize_fixed_control <- function(control.fixed, fixed_effects_names){
  if(missing(control.fixed) || is.null(control.fixed)){
    control.fixed <- list(intercept = list(prec = 0.001, mean = 0))
  }

  if(is.null(control.fixed$intercept)){
    control.fixed$intercept <- list(prec = 0.001, mean = 0)
  }

  if(is.null(control.fixed$intercept$prec)){
    control.fixed$intercept$prec <- 0.001
  }

  if(is.null(control.fixed$intercept$mean)){
    control.fixed$intercept$mean <- 0
  }

  for(name in fixed_effects_names){
    if(is.null(control.fixed[[name]])){
      control.fixed[[name]] <- list(prec = 0.001, mean = 0)
    }
    if(is.null(control.fixed[[name]]$prec)){
      control.fixed[[name]]$prec <- 0.001
    }
    if(is.null(control.fixed[[name]]$mean)){
      control.fixed[[name]]$mean <- 0
    }
  }

  control.fixed
}

resolve_computation_method <- function(computation = NULL, method = NULL,
                                       default = "state-space",
                                       allowed = c("state-space", "fem"),
                                       label = "Smooth terms"){
  if(!is.null(method)){
    stop(
      "Use `computation` inside `f()` for ", tolower(label), ", not `method`. ",
      "`model_fit(method = ...)` is reserved for the inference algorithm."
    )
  }

  if(is.null(computation)){
    return(default)
  }

  computation <- tolower(computation)

  if(!computation %in% allowed){
    stop(label, " require `computation` to be one of: ", paste(sprintf("'%s'", allowed), collapse = ", "), ".")
  }

  computation
}

get_support_grid_diagnostics <- function(support_grid){
  support_grid <- as.numeric(support_grid)
  support_grid <- sort(unique(support_grid))

  list(
    n = length(support_grid),
    span = if(length(support_grid) <= 1) 0 else diff(range(support_grid)),
    min_gap = if(length(support_grid) <= 1) NA_real_ else min(diff(support_grid))
  )
}

validate_exact_state_space_support <- function(support_grid, model_name,
                                               scale_label = "support grid"){
  diagnostics <- get_support_grid_diagnostics(support_grid)
  gap_scale <- max(1, diagnostics$span)
  near_duplicate_tol <- 1e-10 * gap_scale

  if(any(!is.finite(support_grid))){
    stop(
      "Exact state-space `", model_name, "` received non-finite values on the ",
      scale_label, "."
    )
  }

  if(!is.na(diagnostics$min_gap) && diagnostics$min_gap <= near_duplicate_tol){
    stop(
      "Exact state-space `", model_name, "` received nearly overlapping locations on the ",
      scale_label, " (min gap = ", signif(diagnostics$min_gap, 4), "). ",
      "This usually happens when the support grid is too dense or contains near-duplicates. ",
      "Try a sparser `grid`, remove nearly duplicated locations, or switch to ",
      "`computation = \"fem\"`."
    )
  }

  diagnostics
}

build_exact_precision_with_diagnostics <- function(support_grid, model_name,
                                                   precision_builder,
                                                   scale_label = "support grid"){
  diagnostics <- validate_exact_state_space_support(
    support_grid = support_grid,
    model_name = model_name,
    scale_label = scale_label
  )

  if(diagnostics$n == 0){
    return(matrix(0, nrow = 0, ncol = 0))
  }

  precision_matrix <- tryCatch(
    as.matrix(precision_builder(support_grid)),
    error = function(e){
      stop(
        "Unable to construct the exact state-space precision matrix for `",
        model_name, "` on the ", scale_label, ". ",
        "This usually happens when the support grid is too dense or contains near-duplicates. ",
        "Try a sparser `grid`, remove nearly duplicated locations, or switch to ",
        "`computation = \"fem\"`. Original error: ", conditionMessage(e)
      )
    }
  )

  if(any(!is.finite(precision_matrix))){
    stop(
      "Exact state-space `", model_name, "` produced non-finite entries in the precision matrix. ",
      "Try a sparser `grid`, remove nearly duplicated locations, or switch to ",
      "`computation = \"fem\"`."
    )
  }

  reciprocal_condition <- tryCatch(
    as.numeric(rcond(precision_matrix)),
    error = function(e) NA_real_
  )

  chol_error <- tryCatch(
    {
      chol(precision_matrix)
      NULL
    },
    error = function(e) conditionMessage(e)
  )

  if(!is.null(chol_error)){
    stop(
      "Exact state-space `", model_name, "` precision matrix is numerically too close to singular ",
      "(rcond = ", format(reciprocal_condition, digits = 4), ", min gap = ",
      if(is.na(diagnostics$min_gap)) "NA" else format(diagnostics$min_gap, digits = 4),
      ", support points = ", diagnostics$n, " on the ", scale_label, "). ",
      "This usually happens when the support grid is too dense or contains near-duplicates. ",
      "Try a sparser `grid`, remove nearly duplicated locations, or switch to ",
      "`computation = \"fem\"`. Cholesky failed with: ", chol_error
    )
  }

  precision_matrix
}

resolve_left_boundary <- function(x, initial_location){
  left_boundary <- min(x)

  if(is.null(initial_location)){
    return(left_boundary)
  }

  if(is.character(initial_location)){
    if(length(initial_location) != 1 || tolower(initial_location) != "left"){
      stop("Near-monotone smooth terms currently require `initial_location = \"left\"`.")
    }
    return(left_boundary)
  }

  if(!isTRUE(all.equal(as.numeric(initial_location), left_boundary))){
    stop("Near-monotone smooth terms currently require `initial_location` to equal the left boundary of the domain.")
  }

  left_boundary
}

resolve_reference_location <- function(x, initial_location){
  domain_range <- range(x)

  if(is.null(initial_location)){
    return(domain_range[1])
  }

  if(is.character(initial_location)){
    if(length(initial_location) != 1){
      stop("`initial_location` must be a scalar.")
    }
    initial_location <- tolower(initial_location)
    if(initial_location == "left"){
      return(domain_range[1])
    }
    if(initial_location == "middle"){
      return(mean(domain_range))
    }
    if(initial_location == "right"){
      return(domain_range[2])
    }
    stop("`initial_location` must be numeric or one of 'left', 'middle', or 'right'.")
  }

  if(!is.numeric(initial_location) || length(initial_location) != 1 || !is.finite(initial_location)){
    stop("`initial_location` must be numeric or one of 'left', 'middle', or 'right'.")
  }

  tolerance <- sqrt(.Machine$double.eps)
  if(initial_location < domain_range[1] - tolerance || initial_location > domain_range[2] + tolerance){
    stop("`initial_location` must lie inside the near-monotone domain.")
  }

  as.numeric(initial_location)
}

make_boxcox_transform <- function(a = NULL, alpha = NULL, c = 1,
                                  normalized = TRUE, ref_location = 0){
  curvature <- resolve_curvature_parameter(a = a, alpha = alpha)
  lambda <- (curvature - 1) / curvature
  reference <- c + ref_location

  if(normalized){
    if(reference <= 0){
      stop("Normalized Box-Cox transforms require `c + ref_location > 0`.")
    }
    if(isTRUE(all.equal(lambda, 0))){
      return(function(x){
        reference * (log(x + c) - log(reference))
      })
    }

    return(function(x){
      ((x + c)^lambda - reference^lambda) / (lambda * reference^(lambda - 1))
    })
  }

  if(isTRUE(all.equal(lambda, 0))){
    return(function(x){
      log(x + c)
    })
  }

  function(x){
    (x + c)^lambda
  }
}

make_boxcox_transform_derivative <- function(x, a = NULL, alpha = NULL, c = 1,
                                             normalized = TRUE, ref_location = 0){
  curvature <- resolve_curvature_parameter(a = a, alpha = alpha)
  lambda <- (curvature - 1) / curvature
  reference <- c + ref_location

  if(normalized){
    if(reference <= 0){
      stop("Normalized Box-Cox transforms require `c + ref_location > 0`.")
    }
    return(((x + c) / reference)^(lambda - 1))
  }

  if(isTRUE(all.equal(lambda, 0))){
    return(1 / (x + c))
  }

  lambda * (x + c)^(lambda - 1)
}

normalize_psd_model <- function(model) {
  model <- tolower(model)

  if (model %in% c("mgp", "tiwp2", "iwp", "sgp")) {
    return(model)
  }

  if (model %in% c("tiwp", "tiwp-2", "t-iwp2")) {
    return("tiwp2")
  }

  stop(
    "`model` must be one of 'mgp', 'tiwp2', 'iwp', or 'sgp'."
  )
}

normalize_smooth_term_model <- function(model) {
  model <- tolower(as.character(model))
  supported <- c("iid", "iwp", "sgp", "mgp", "tiwp2")

  if(length(model) != 1 || is.na(model) || !nzchar(model)){
    stop("`model` must be a single smooth-term model name.")
  }

  if(model %in% supported){
    return(model)
  }

  if(model == "iwp2"){
    stop(
      "Unknown smooth-term model `iwp2`. Use `model = \"tiwp2\"` for ",
      "transformed IWP2, or use `model = \"iwp\", order = 2` for native IWP2."
    )
  }

  stop(
    "Unknown smooth-term model `", model, "`. Use supported_models() to list ",
    "available model names."
  )
}

PSD_compute_tiwp2 <- function(x, h, a = NULL, alpha = NULL,
                              c = 1, sd = 1, normalized = TRUE){
  transform <- make_boxcox_transform(
    a = a,
    alpha = alpha,
    c = c,
    normalized = normalized
  )
  transformed_step <- transform(x + h) - transform(x)
  abs(sd) * sqrt((transformed_step^3) / 3)
}

PSD_compute <- function(model, h, x = NULL, sd = 1, ...) {
  model <- normalize_psd_model(model)

  if (model %in% c("mgp", "tiwp2") && is.null(x)) {
    stop("`x` must be supplied when `model` is 'mgp' or 'tiwp2'.")
  }

  switch(
    model,
    mgp = PSD_compute_mgp(x = x, h = h, sd = sd, ...),
    tiwp2 = PSD_compute_tiwp2(x = x, h = h, sd = sd, ...),
    iwp = PSD_compute_iwp(h = h, sd = sd, ...),
    sgp = PSD_compute_sgp(h = h, sd = sd, ...),
    stop("Unsupported PSD model `", model, "`.")
  )
}

simulate_tiwp_boxcox <- function(x, a = NULL, alpha = NULL, c = 1, sd = 1,
                                 initial_level = 0, initial_slope = 1,
                                 normalized_transform = TRUE){
  transform <- make_boxcox_transform(
    a = a,
    alpha = alpha,
    c = c,
    normalized = normalized_transform,
    ref_location = x[1]
  )
  transformed_x <- transform(x)
  transformed_initial_slope <- initial_slope / make_boxcox_transform_derivative(
    x = x[1],
    a = a,
    alpha = alpha,
    c = c,
    normalized = normalized_transform,
    ref_location = x[1]
  )

  sim_IWp_Var(
    t = transformed_x,
    p = 2,
    sd = sd,
    initial_vec = c(initial_level, transformed_initial_slope)
  )[, 2]
}

boundary_basis_matrix <- function(x, a = NULL, alpha = NULL, c = 1,
                                  initial_location = min(x),
                                  normalized = TRUE){
  curvature <- resolve_curvature_parameter(a = a, alpha = alpha)
  shifted_x <- x - initial_location
  transformed_arg <- shifted_x + c

  if(any(transformed_arg <= 0)){
    stop("All evaluation points must satisfy `x + c > 0` on the original scale.")
  }

  transform <- make_boxcox_transform(
    a = curvature,
    c = c,
    normalized = normalized,
    ref_location = 0
  )

  matrix(transform(shifted_x), ncol = 1)
}

build_exact_observation_matrix <- function(observed_x, support_x){
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

extract_term_argument <- function(rand_effect, name, envir = parent.frame(), default = NULL){
  if(name %in% names(rand_effect)){
    return(eval(rand_effect[[name]], envir = envir))
  }

  default
}

build_fixed_design <- function(parse_result, data, family){
  fixed_effects <- parse_result$fixed_effects
  offset_effects <- parse_result$offset_effects
  family <- tolower(family)
  family_is_coxph <- family %in% c("cox", "coxph")
  family_is_cc <- family %in% c("casecrossover", "cc")

  design_mat_fixed <- list()
  fixed_effects_names <- character()

  if(!family_is_coxph && !family_is_cc){
    design_mat_fixed[[length(design_mat_fixed) + 1]] <- matrix(1, nrow = nrow(data), ncol = 1)
  }

  for(fixed_effect in fixed_effects){
    fixed_name <- toString(fixed_effect)

    if(is.null(data[[fixed_name]])){
      stop("Fixed effects must be plain column names available in `data`.")
    }

    if(is.numeric(data[[fixed_name]])){
      design_mat_fixed[[length(design_mat_fixed) + 1]] <- matrix(data[[fixed_name]], ncol = 1)
      fixed_effects_names <- c(fixed_effects_names, fixed_name)
    } else {
      level <- data[[fixed_name]]
      Xf <- stats::model.matrix(~level)[, -1, drop = FALSE]
      for(col_index in seq_len(ncol(Xf))){
        design_mat_fixed[[length(design_mat_fixed) + 1]] <- Xf[, col_index, drop = FALSE]
      }
      fixed_effects_names <- c(fixed_effects_names, colnames(Xf))
    }
  }

  offset_sum <- 0
  for(offset_effect in offset_effects){
    offset_name <- toString(offset_effect[[2]])
    offset_sum <- offset_sum + data[[offset_name]]
  }

  list(
    design_mat_fixed = design_mat_fixed,
    fixed_effects_names = fixed_effects_names,
    offset_sum = offset_sum
  )
}

collect_extended_terms <- function(parse_result){
  rand_effects <- parse_result$rand_effects
  model_names <- vapply(rand_effects, function(term){
    tolower(as.character(term$model))
  }, character(1))

  list(
    extended = rand_effects[model_names %in% c("mgp", "tiwp2")],
    native = rand_effects[!model_names %in% c("mgp", "tiwp2")]
  )
}

extract_random_samples <- function(object, component){
  component_index <- object$random_samp_indexes[[component]]
  if(is.null(component_index)){
    stop("Unknown component `", component, "`.")
  }

  as.matrix(object$samps$samps[component_index, , drop = FALSE])
}

extract_boundary_samples <- function(object, component){
  component_index <- object$boundary_samp_indexes[[component]]

  if(is.null(component_index) || length(component_index) == 0){
    return(matrix(0, nrow = 0, ncol = ncol(object$samps$samps)))
  }

  as.matrix(object$samps$samps[component_index, , drop = FALSE])
}

extract_intercept_samples <- function(object){
  intercept_index <- object$fixed_samp_indexes[["intercept"]]

  if(is.null(intercept_index)){
    return(matrix(0, nrow = 1, ncol = ncol(object$samps$samps)))
  }

  as.matrix(object$samps$samps[intercept_index, , drop = FALSE])
}

summarize_function_samples <- function(samples, probs = c(0.025, 0.975)){
  data.frame(
    mean = rowMeans(samples),
    lower = apply(samples, 1, quantile, probs = probs[1]),
    upper = apply(samples, 1, quantile, probs = probs[2])
  )
}

evaluate_train_prediction_draws <- function(samples, truth, n_train){
  safe_mean <- function(values){
    if(length(values) == 0){
      return(NA_real_)
    }
    mean(values)
  }

  compute_interval_metrics <- function(probability){
    interval <- apply(
      samples,
      1,
      quantile,
      probs = c((1 - probability) / 2, 1 - (1 - probability) / 2)
    )

    interpolation_index <- seq_len(n_train)
    prediction_index <- setdiff(seq_along(truth), interpolation_index)

    list(
      inter = safe_mean(
        truth[interpolation_index] >= interval[1, interpolation_index] &
          truth[interpolation_index] <= interval[2, interpolation_index]
      ),
      pred = safe_mean(
        truth[prediction_index] >= interval[1, prediction_index] &
          truth[prediction_index] <= interval[2, prediction_index]
      )
    )
  }

  mean_draw <- rowMeans(samples)
  interpolation_index <- seq_len(n_train)
  prediction_index <- setdiff(seq_along(truth), interpolation_index)

  coverage_50 <- compute_interval_metrics(0.5)
  coverage_80 <- compute_interval_metrics(0.8)
  coverage_95 <- compute_interval_metrics(0.95)

  list(
    data_inter = data.frame(
      coverage_50 = coverage_50$inter,
      coverage_80 = coverage_80$inter,
      coverage_95 = coverage_95$inter,
      rmse = sqrt(safe_mean((mean_draw[interpolation_index] - truth[interpolation_index])^2)),
      rmle = safe_mean(abs((mean_draw[interpolation_index] - truth[interpolation_index]) / truth[interpolation_index]))
    ),
    data_pred = data.frame(
      coverage_50 = coverage_50$pred,
      coverage_80 = coverage_80$pred,
      coverage_95 = coverage_95$pred,
      rmse = sqrt(safe_mean((mean_draw[prediction_index] - truth[prediction_index])^2)),
      rmle = safe_mean(abs((mean_draw[prediction_index] - truth[prediction_index]) / truth[prediction_index]))
    )
  )
}

evaluate_global_interval_draws <- function(samples, truth, alpha = 0.2){
  lower <- apply(samples, 1, quantile, probs = alpha / 2)
  upper <- apply(samples, 1, quantile, probs = 1 - alpha / 2)

  data.frame(
    coverage = mean(truth >= lower & truth <= upper),
    interval_width = mean(upper - lower)
  )
}

evaluate_regionwise_interval_draws <- function(samples, truth, x,
                                               alpha = 0.2,
                                               breaks = quantile(x, probs = c(0, 1 / 3, 2 / 3, 1)),
                                               labels = c("left", "middle", "right")){
  lower <- apply(samples, 1, quantile, probs = alpha / 2)
  upper <- apply(samples, 1, quantile, probs = 1 - alpha / 2)
  region_index <- findInterval(x, vec = breaks, all.inside = TRUE, rightmost.closed = TRUE)

  result <- data.frame(
    coverage = vapply(seq_along(labels), function(i){
      selected <- region_index == i
      mean(truth[selected] >= lower[selected] & truth[selected] <= upper[selected])
    }, numeric(1)),
    interval_width = vapply(seq_along(labels), function(i){
      selected <- region_index == i
      mean(upper[selected] - lower[selected])
    }, numeric(1))
  )

  rownames(result) <- labels
  result
}
