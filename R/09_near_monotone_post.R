#' Extract posterior smooth draws
#'
#' @param object Fitted object returned by \code{\link{model_fit}}.
#' @param component Smooth-term name.
#' @param refined_x Optional evaluation grid. Exact state-space terms can only
#'   be evaluated on points included in the term's support grid.
#' @param include_intercept Whether to add the model intercept draw.
#'
#' @return A data frame with `x` in the first column and posterior draws in the
#'   remaining columns.
#' @export
smooth_samples <- function(object, component, refined_x = NULL, include_intercept = FALSE){
  instance <- NULL
  metadata <- object$near_mono_meta[[component]]

  if(is.null(metadata)){
    stop("Unknown component `", component, "`.")
  }

  for(candidate in object$instances){
    if(as.character(candidate@smoothing_var) == component){
      instance <- candidate
      break
    }
  }

  coef_samps <- extract_random_samples(object, component = component)
  boundary_samps <- extract_boundary_samples(object, component = component)
  intercept_samps <- if(include_intercept) extract_intercept_samples(object) else matrix(0, nrow = 1, ncol = ncol(coef_samps))

  if(metadata$computation_method == "state-space" && metadata$model_name == "mgp" && !isTRUE(metadata$reference_is_left)){
    eval_x <- if(is.null(refined_x)) {
      c(
        instance@initial_location - rev(metadata$raw_negative_support),
        instance@initial_location,
        instance@initial_location + metadata$raw_positive_support
      )
    } else {
      refined_x
    }
    shifted_eval_x <- eval_x - instance@initial_location
    validate_two_sided_exact_evaluation(
      x = shifted_eval_x,
      positive_support = metadata$positive_support,
      negative_support = metadata$negative_support
    )
    design_eval <- build_two_sided_exact_observation_matrix(
      observed_x = shifted_eval_x,
      positive_support = metadata$positive_support,
      negative_support = metadata$negative_support,
      state_size = 2
    )
    boundary_eval <- boundary_basis_matrix(
      x = eval_x,
      a = metadata$curvature,
      c = metadata$shift,
      initial_location = instance@initial_location,
      normalized = metadata$normalized_boundary
    )
    fitted_draws <- design_eval %*% coef_samps
    if(nrow(boundary_samps) > 0){
      fitted_draws <- fitted_draws + boundary_eval %*% boundary_samps
    }
    fitted_draws <- fitted_draws + matrix(rep(intercept_samps, each = nrow(fitted_draws)), nrow = nrow(fitted_draws))
  } else if(metadata$computation_method == "state-space" && metadata$model_name == "tiwp2"){
    eval_x <- if(is.null(refined_x)) {
      c(
        instance@initial_location - rev(metadata$raw_negative_support),
        instance@initial_location,
        instance@initial_location + metadata$raw_positive_support
      )
    } else {
      refined_x
    }
    shifted_eval_x <- eval_x - instance@initial_location
    transform <- make_boxcox_transform(
      a = metadata$curvature,
      c = metadata$shift,
      normalized = TRUE,
      ref_location = 0
    )
    transformed_eval_x <- transform(shifted_eval_x)
    validate_two_sided_exact_evaluation(
      x = transformed_eval_x,
      positive_support = metadata$positive_support,
      negative_support = metadata$negative_support
    )
    design_eval <- build_two_sided_exact_observation_matrix(
      observed_x = transformed_eval_x,
      positive_support = metadata$positive_support,
      negative_support = metadata$negative_support,
      state_size = 2
    )
    boundary_eval <- boundary_basis_matrix(
      x = eval_x,
      a = metadata$curvature,
      c = metadata$shift,
      initial_location = instance@initial_location,
      normalized = metadata$normalized_boundary
    )
    fitted_draws <- design_eval %*% coef_samps
    if(nrow(boundary_samps) > 0){
      fitted_draws <- fitted_draws + boundary_eval %*% boundary_samps
    }
    fitted_draws <- fitted_draws + matrix(rep(intercept_samps, each = nrow(fitted_draws)), nrow = nrow(fitted_draws))
  } else if(metadata$computation_method == "state-space"){
    available_x <- c(instance@initial_location, instance@initial_location + metadata$support_grid)
    eval_x <- if(is.null(refined_x)) available_x else refined_x

    available_lookup <- match(eval_x, available_x)
    if(anyNA(available_lookup)){
      stop("Exact state-space evaluation currently requires `refined_x` to be included in the term support. Supply `grid = ...` in `f()` when fitting if extra points are needed.")
    }

    process_draws <- matrix(0, nrow = length(available_x), ncol = ncol(coef_samps))
    if(length(metadata$support_grid) > 0){
      process_draws[-1, ] <- coef_samps[seq(1, 2 * length(metadata$support_grid), by = 2), , drop = FALSE]
    }

    boundary_eval <- boundary_basis_matrix(
      x = eval_x,
      a = metadata$curvature,
      c = metadata$shift,
      initial_location = instance@initial_location,
      normalized = metadata$normalized_boundary
    )
    fitted_draws <- process_draws[available_lookup, , drop = FALSE]
    if(nrow(boundary_samps) > 0){
      fitted_draws <- fitted_draws + boundary_eval %*% boundary_samps
    }
    fitted_draws <- fitted_draws + matrix(rep(intercept_samps, each = nrow(fitted_draws)), nrow = nrow(fitted_draws))
  } else if(metadata$computation_method == "fem" && metadata$model_name == "mgp"){
    eval_x <- if(is.null(refined_x)) sort(unique(instance@data[[as.character(instance@smoothing_var)]])) else refined_x
    shifted_eval_x <- eval_x - instance@initial_location
    if(!isTRUE(metadata$reference_is_left)){
      shifted_domain <- instance@region - instance@initial_location
      design_eval <- cbind(
        build_side_design(
          x = shifted_eval_x,
          side = "positive",
          k = instance@k,
          region = range(c(0, shifted_domain[shifted_domain > 0])),
          boundary = metadata$fem_boundary,
          basis_builder = function(x, k, region, boundary){
            as.matrix(Compute_Design(x = x, k = k, region = region, boundary = boundary))
          }
        ),
        build_side_design(
          x = shifted_eval_x,
          side = "negative",
          k = instance@k,
          region = range(c(0, -shifted_domain[shifted_domain < 0])),
          boundary = metadata$fem_boundary,
          basis_builder = function(x, k, region, boundary){
            as.matrix(Compute_Design(x = x, k = k, region = region, boundary = boundary))
          }
        )
      )
    } else {
      design_eval <- as.matrix(
        Compute_Design(
          shifted_eval_x,
          k = instance@k,
          region = range(instance@region - instance@initial_location),
          boundary = metadata$fem_boundary
        )
      )
    }
    boundary_eval <- boundary_basis_matrix(
      x = eval_x,
      a = metadata$curvature,
      c = metadata$shift,
      initial_location = instance@initial_location,
      normalized = metadata$normalized_boundary
    )
    fitted_draws <- design_eval %*% coef_samps
    if(nrow(boundary_samps) > 0){
      fitted_draws <- fitted_draws + boundary_eval %*% boundary_samps
    }
    fitted_draws <- fitted_draws + matrix(rep(intercept_samps, each = nrow(fitted_draws)), nrow = nrow(fitted_draws))
  } else if(metadata$computation_method == "fem" && metadata$model_name == "tiwp2"){
    eval_x <- if(is.null(refined_x)) sort(unique(instance@data[[as.character(instance@smoothing_var)]])) else refined_x
    shifted_eval_x <- eval_x - instance@initial_location
    transform <- make_boxcox_transform(
      a = metadata$curvature,
      c = metadata$shift,
      normalized = TRUE,
      ref_location = 0
    )
    transformed_eval_x <- transform(shifted_eval_x)
    if(!isTRUE(metadata$reference_is_left)){
      design_eval <- cbind(
        build_tiwp_side_design(
          x = transformed_eval_x,
          side = "positive",
          knots = metadata$transformed_knots_positive,
          response_var = "y",
          smoothing_var = as.character(instance@smoothing_var),
          sd.prior = instance@sd.prior,
          psd.prior = instance@psd.prior,
          boundary.prior = instance@boundary.prior
        ),
        build_tiwp_side_design(
          x = transformed_eval_x,
          side = "negative",
          knots = metadata$transformed_knots_negative,
          response_var = "y",
          smoothing_var = as.character(instance@smoothing_var),
          sd.prior = instance@sd.prior,
          psd.prior = instance@psd.prior,
          boundary.prior = instance@boundary.prior
        )
      )
      boundary_eval <- boundary_basis_matrix(
        x = eval_x,
        a = metadata$curvature,
        c = metadata$shift,
        initial_location = instance@initial_location,
        normalized = metadata$normalized_boundary
      )
      fitted_draws <- design_eval %*% coef_samps
      if(nrow(boundary_samps) > 0){
        fitted_draws <- fitted_draws + boundary_eval %*% boundary_samps
      }
      fitted_draws <- fitted_draws + matrix(rep(intercept_samps, each = nrow(fitted_draws)), nrow = nrow(fitted_draws))
    } else {
      fitted_draws <- compute_post_fun_iwp(
        samps = coef_samps,
        global_samps = if(nrow(boundary_samps) == 0) NULL else boundary_samps,
        knots = metadata$transformed_knots,
        refined_x = transformed_eval_x,
        p = 2,
        intercept_samps = intercept_samps
      )
      fitted_draws <- as.matrix(fitted_draws[, -1, drop = FALSE])
    }
  } else {
    stop("Unsupported component type.")
  }

  data.frame(x = eval_x, as.data.frame(as.matrix(fitted_draws)))
}

#' Summarize posterior smooth draws
#'
#' @inheritParams smooth_samples
#' @param level Pointwise interval level.
#'
#' @return Posterior mean and interval summary.
#' @export
smooth_summary <- function(object, component, refined_x = NULL,
                           include_intercept = FALSE, level = 0.95){
  samps <- smooth_samples(
    object = object,
    component = component,
    refined_x = refined_x,
    include_intercept = include_intercept
  )
  extract_mean_interval_given_samps(samps = samps, level = level)
}

validate_two_sided_exact_evaluation <- function(x, positive_support, negative_support){
  positive_x <- x[x > 0]
  negative_x <- abs(x[x < 0])

  if(length(positive_x) > 0 && anyNA(match(positive_x, positive_support))){
    stop("Exact state-space evaluation currently requires `refined_x` to be included in the term support. Supply `grid = ...` in `f()` when fitting if extra points are needed.")
  }
  if(length(negative_x) > 0 && anyNA(match(negative_x, negative_support))){
    stop("Exact state-space evaluation currently requires `refined_x` to be included in the term support. Supply `grid = ...` in `f()` when fitting if extra points are needed.")
  }

  invisible(TRUE)
}
