#' @export
#' @method summary FitResult
summary.FitResult <- function(object, ...) {
  output <- list(
    fit = fit_result_info(object),
    smooth_terms = fit_result_smooth_terms(object),
    parameters = post_table(object)
  )
  class(output) <- c("summary.FitResult", "list")
  return(output)
}

#' @export
#' @method print FitResult
print.FitResult <- function(x, ...) {
  info <- fit_result_info(x)
  smooth_terms <- fit_result_smooth_terms(x)

  cat("BayesGP FitResult\n")
  cat("Family:", info$family, "\n")
  if(isTRUE(info$gaussian_sd_known)){
    cat("Gaussian SD:", format(info$gaussian_sd_value, digits = 6), "(fixed)\n")
  }
  if(fit_result_is_case_crossover(info)){
    cat("Strata:", info$strata, "(n =", info$n_strata, ")\n")
  }
  if(is.finite(info$log_marginal_likelihood)){
    cat("Log marginal likelihood:", format(info$log_marginal_likelihood, digits = 8), "\n")
  } else {
    cat("Log marginal likelihood: unavailable\n")
  }
  cat("Smooth terms:", info$n_smooth_terms, "\n")
  cat("Posterior samples:", info$n_posterior_samples, "\n")

  if(nrow(smooth_terms) > 0){
    cat("\nSmooth term specification:\n")
    print(smooth_terms, row.names = FALSE)
  }

  cat("\nUse summary(object) for posterior parameter summaries.\n")
  invisible(x)
}

#' @export
#' @method print summary.FitResult
print.summary.FitResult <- function(x, ...) {
  cat("BayesGP FitResult Summary\n\n")
  if(is.finite(x$fit$log_marginal_likelihood)){
    cat("Log marginal likelihood:", format(x$fit$log_marginal_likelihood, digits = 8), "\n\n")
  }
  print(x$fit, row.names = FALSE)

  if(nrow(x$smooth_terms) > 0){
    cat("\nSmooth term specification:\n")
    print(x$smooth_terms, row.names = FALSE)
  }

  cat("\nPosterior/prior summaries for parameters:\n")
  output <- capture.output(print(as.data.frame(x$parameters), row.names = FALSE))
  cat(output, sep = "\n")
  cat("For Normal prior, P1 is its mean and P2 is its variance. \n")
  cat("For Exponential prior, prior is specified as P(theta > P1) = P2. \n")
  invisible(x)
}

fit_result_info <- function(object){
  samps <- tryCatch(object$samps$samps, error = function(e) NULL)
  info <- data.frame(
    family = if(is.null(object$family)) NA_character_ else object$family,
    log_marginal_likelihood = fit_result_log_marginal(object),
    n_smooth_terms = length(object$instances),
    n_fixed_effects = length(object$fixed_samp_indexes),
    n_posterior_samples = if(is.null(samps)) NA_integer_ else ncol(samps),
    stringsAsFactors = FALSE
  )
  family <- tolower(info$family[[1]])

  if(identical(family, "gaussian")){
    info$gaussian_sd_known <- isTRUE(object$family_sd_known)
    info$gaussian_sd_value <- if(is.null(object$family_sd_value)) NA_real_ else object$family_sd_value
    return(info)
  }

  family_info <- fit_result_family_info(object)
  if(family %in% c("casecrossover", "cc")){
    info$strata <- family_info_string(family_info$strata)
    info$n_strata <- family_info_integer(family_info$n_strata)
    info$max_stratum_size <- family_info_integer(family_info$max_stratum_size)
    info$count_source <- family_info_string(family_info$count_source)
    info$total_count <- family_info_numeric(family_info$total_count)
  } else if(identical(family, "binomial")){
    info$size <- family_info_string(family_info$size)
    info$total_trials <- family_info_numeric(family_info$total_trials)
  } else if(family %in% c("cox", "coxph")){
    info$censoring <- family_info_string(family_info$censoring)
  }

  info
}

fit_result_log_marginal <- function(object){
  value <- tryCatch(object$mod$normalized_posterior$lognormconst, error = function(e) NULL)
  if(is.null(value) || length(value) == 0){
    return(NA_real_)
  }
  as.numeric(value[[1]])
}

fit_result_is_case_crossover <- function(info){
  family <- tolower(info$family[[1]])
  family %in% c("casecrossover", "cc") && "strata" %in% names(info)
}

fit_result_family_info <- function(object){
  if(is.null(object$family_info)){
    return(list())
  }
  object$family_info
}

build_family_fit_info <- function(family, data, response_var, size = NULL, cens = NULL,
                                  weight = NULL, strata = NULL){
  family <- tolower(family)

  if(family %in% c("casecrossover", "cc")){
    strata_name <- family_info_name(strata)
    strata_values <- family_info_column(data, strata_name)
    count_source <- family_info_name(weight, default = response_var)
    count_values <- family_info_column(data, count_source)

    return(list(
      strata = strata_name,
      n_strata = unique_value_count(strata_values),
      max_stratum_size = max_group_size(strata_values),
      count_source = count_source,
      total_count = numeric_sum_or_na(count_values)
    ))
  }

  if(identical(family, "binomial")){
    size_name <- family_info_name(size, default = "1")
    size_values <- family_info_column(data, size_name)
    if(is.null(size_values)){
      size_values <- rep(1, nrow(data))
    }

    return(list(
      size = size_name,
      total_trials = numeric_sum_or_na(size_values)
    ))
  }

  if(family %in% c("cox", "coxph")){
    return(list(
      censoring = family_info_name(cens, default = "all observed")
    ))
  }

  list()
}

family_info_name <- function(x, default = NA_character_){
  default <- family_info_scalar_character(default)
  if(is.null(x) || length(x) == 0){
    return(default)
  }
  value <- family_info_scalar_character(x)
  if(length(value) == 0 || is.na(value) || identical(value, "NULL")){
    return(default)
  }
  value
}

family_info_column <- function(data, column){
  column <- family_info_name(column)
  if(is.null(column) || length(column) == 0 || is.na(column) || !column %in% names(data)){
    return(NULL)
  }
  data[[column]]
}

family_info_scalar_character <- function(x){
  if(is.null(x) || length(x) == 0){
    return(NA_character_)
  }
  as.character(x)[[1]]
}

unique_value_count <- function(x){
  if(is.null(x)){
    return(NA_integer_)
  }
  length(unique(x))
}

max_group_size <- function(x){
  if(is.null(x) || length(x) == 0){
    return(NA_integer_)
  }
  as.integer(max(table(x)))
}

numeric_sum_or_na <- function(x){
  if(is.null(x)){
    return(NA_real_)
  }
  sum(as.numeric(x))
}

family_info_string <- function(x){
  if(is.null(x) || length(x) == 0){
    return(NA_character_)
  }
  as.character(x[[1]])
}

family_info_integer <- function(x){
  if(is.null(x) || length(x) == 0){
    return(NA_integer_)
  }
  as.integer(x[[1]])
}

family_info_numeric <- function(x){
  if(is.null(x) || length(x) == 0){
    return(NA_real_)
  }
  as.numeric(x[[1]])
}

fit_result_smooth_terms <- function(object){
  if(length(object$instances) == 0){
    return(data.frame(
      component = character(),
      model = character(),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(object$instances, function(instance){
    component <- as.character(slot_or_default(instance, "smoothing_var", NA_character_))
    instance_class <- class(instance)[1]
    metadata <- object$near_mono_meta[[component]]
    exact_iwp_metadata <- object$exact_iwp_meta[[component]]

    if(!is.null(metadata)){
      return(near_mono_smooth_term_row(component, instance, metadata))
    }

    if(!is.null(exact_iwp_metadata)){
      return(exact_iwp_smooth_term_row(component, instance, exact_iwp_metadata))
    }

    if(instance_class == "iid"){
      return(iid_smooth_term_row(component, instance))
    }

    if(instance_class == "iwp"){
      return(iwp_smooth_term_row(component, instance))
    }

    if(instance_class == "sgp"){
      return(sgp_smooth_term_row(component, instance))
    }

    data.frame(
      component = component,
      model = instance_class,
      stringsAsFactors = FALSE
    )
  })

  bind_smooth_term_rows(rows)
}

iid_smooth_term_row <- function(component, instance){
  data.frame(
    component = component,
    model = "iid",
    n_levels = instance_level_count(instance, component),
    stringsAsFactors = FALSE
  )
}

iwp_smooth_term_row <- function(component, instance){
  data.frame(
    component = component,
    model = "iwp",
    computation = "fem",
    order = scalar_or_default(slot_or_default(instance, "order", NA_real_), NA_real_),
    basis_size = scalar_or_default(slot_or_default(instance, "k", NA_real_), NA_real_),
    domain = format_numeric_range(slot_or_default(instance, "region", numeric())),
    reference = scalar_or_default(slot_or_default(instance, "initial_location", NA_real_), NA_real_),
    stringsAsFactors = FALSE
  )
}

exact_iwp_smooth_term_row <- function(component, instance, metadata){
  data.frame(
    component = component,
    model = "iwp",
    computation = metadata$computation_method,
    order = scalar_or_default(metadata$order, NA_real_),
    basis_size = ncol_or_default(slot_or_default(instance, "B", NULL)),
    support_size = length(metadata$support_grid),
    domain = format_numeric_range(slot_or_default(instance, "region", numeric())),
    reference = scalar_or_default(slot_or_default(instance, "initial_location", NA_real_), NA_real_),
    stringsAsFactors = FALSE
  )
}

sgp_smooth_term_row <- function(component, instance){
  angular_frequency <- scalar_or_default(slot_or_default(instance, "a", NA_real_), NA_real_)
  data.frame(
    component = component,
    model = "sgp",
    computation = "fem",
    basis_size = scalar_or_default(slot_or_default(instance, "k", NA_real_), NA_real_),
    domain = format_numeric_range(slot_or_default(instance, "region", numeric())),
    reference = scalar_or_default(slot_or_default(instance, "initial_location", NA_real_), NA_real_),
    angular_frequency = angular_frequency,
    frequency = if(is.finite(angular_frequency)) angular_frequency / (2 * pi) else NA_real_,
    period = if(is.finite(angular_frequency) && angular_frequency != 0) (2 * pi) / angular_frequency else NA_real_,
    harmonics = scalar_or_default(slot_or_default(instance, "m", NA_real_), NA_real_),
    stringsAsFactors = FALSE
  )
}

near_mono_smooth_term_row <- function(component, instance, metadata){
  is_state_space <- identical(metadata$computation_method, "state-space")
  data.frame(
    component = component,
    model = metadata$model_name,
    computation = metadata$computation_method,
    basis_size = scalar_or_default(slot_or_default(instance, "k", NA_real_), NA_real_),
    support_size = if(is_state_space) length(metadata$support_grid) else NA_integer_,
    domain = format_numeric_range(slot_or_default(instance, "region", numeric())),
    reference = scalar_or_default(metadata$reference_location, NA_real_),
    a = scalar_or_default(metadata$curvature, NA_real_),
    c = scalar_or_default(if(!is.null(metadata$c)) metadata$c else metadata$shift, NA_real_),
    stringsAsFactors = FALSE
  )
}

bind_smooth_term_rows <- function(rows){
  all_columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  preferred_order <- c(
    "component", "model", "computation", "order", "basis_size",
    "support_size", "n_levels", "domain", "reference", "a", "c",
    "angular_frequency", "frequency", "period", "harmonics"
  )
  all_columns <- c(preferred_order[preferred_order %in% all_columns], setdiff(all_columns, preferred_order))

  rows <- lapply(rows, function(row){
    for(column in setdiff(all_columns, names(row))){
      row[[column]] <- NA
    }
    row[all_columns]
  })

  output <- do.call(rbind, rows)
  row.names(output) <- NULL
  drop_empty_smooth_term_columns(output)
}

drop_empty_smooth_term_columns <- function(x){
  keep <- c("component", "model")
  informative <- vapply(names(x), function(column){
    if(column %in% keep){
      return(TRUE)
    }
    values <- x[[column]]
    !all(is.na(values))
  }, logical(1))

  x[, informative, drop = FALSE]
}

instance_level_count <- function(instance, component){
  data <- slot_or_default(instance, "data", NULL)
  if(is.data.frame(data) && component %in% names(data)){
    return(length(unique(data[[component]])))
  }

  ncol_or_default(slot_or_default(instance, "B", NULL))
}

slot_or_default <- function(object, slot, default = NULL){
  if(isS4(object) && slot %in% methods::slotNames(object)){
    return(methods::slot(object, slot))
  }
  default
}

scalar_or_default <- function(x, default = NA_real_){
  if(is.null(x) || length(x) == 0){
    return(default)
  }
  x[[1]]
}

ncol_or_default <- function(x, default = NA_integer_){
  if(is.null(x) || is.null(dim(x))){
    return(default)
  }
  ncol(x)
}

format_numeric_range <- function(x){
  if(!is.numeric(x) || length(x) == 0 || any(!is.finite(x))){
    return(NA_character_)
  }
  x <- range(x)
  paste0("[", format(x[1], digits = 6), ", ", format(x[2], digits = 6), "]")
}


#' To predict the GP component in the fitted model, at the locations specified in `newdata`. 
#' @param object The fitted object from the function `model_fit`.
#' @param newdata The dataset that contains the locations to be predicted for the specified GP. Its column names must include `variable`.
#' @param variable The name of the variable to be predicted, should be in the `newdata`.
#' @param deriv The degree of derivative that the user specifies for inference. Only applicable for a GP in the `iwp` type.
#' @param include.intercept A logical variable specifying whether the intercept should be accounted when doing the prediction. The default is TRUE. For Coxph model, this 
#' variable will be forced to FALSE.
#' @param only.samples A logical variable indicating whether only the posterior samples are required. The default is FALSE, and the summary of posterior samples will be reported.
#' @param quantiles A numeric vector of quantiles that predict.FitResult will produce, the default is c(0.025, 0.5, 0.975).
#' @param boundary.condition A string specifies whether the boundary.condition should be considered in the prediction, should be one of c("yes", "no", "only"). The default option is "Yes".
#' @param ... Additional arguments. Currently unused.
#'
#' @details
#' For near-monotone terms (`model = "mgp"` or `model = "tiwp2"`) and for
#' exact `iwp` state-space terms (`model = "iwp", computation =
#' "state-space"`), `predict()` also supports `only.samples = TRUE` to return
#' posterior draws directly. Exact
#' state-space fits can only be evaluated on the support grid used at fit
#' time; provide those locations through `grid = ...` inside `f()` when
#' fitting if predictions are needed at additional points.
#' @export
predict.FitResult <- function(object, newdata = NULL, variable, deriv = 0, include.intercept = TRUE, only.samples = FALSE, quantiles = c(0.025, 0.5, 0.975), boundary.condition = "Yes", ...) {
  if(object$family == "Coxph" || object$family == "coxph"| object$family == "cc" | object$family == "casecrossover" | object$family == "CaseCrossover"){
    include.intercept = FALSE ## No intercept for coxph model
  }
  if(!is.null(object$near_mono_meta[[variable]])){
    eval_x <- if(is.null(newdata)) NULL else sort(newdata[[variable]])
    f <- smooth_samples(
      object = object,
      component = variable,
      refined_x = eval_x,
      include_intercept = include.intercept
    )
    if(only.samples){
      names(f)[-1] <- paste0("samp", seq_len(ncol(f) - 1))
      return(f)
    }
    fpos <- extract_mean_interval_given_samps(f, quantiles = quantiles)
    names(fpos)[names(fpos) == "x"] <- variable
    return(fpos)
  }
  if(!is.null(object$exact_iwp_meta[[variable]])){
    eval_x <- if(is.null(newdata)) NULL else sort(newdata[[variable]])
    f <- exact_iwp_samples(
      object = object,
      component = variable,
      refined_x = eval_x,
      include_intercept = include.intercept,
      deriv = deriv
    )
    if(only.samples){
      names(f)[-1] <- paste0("samp", seq_len(ncol(f) - 1))
      return(f)
    }
    fpos <- extract_mean_interval_given_samps(f, quantiles = quantiles)
    names(fpos)[names(fpos) == "x"] <- variable
    return(fpos)
  }
  samps <- object$samps
  for (instance in object$instances) {
    if(sum(names(object$random_samp_indexes) == variable) >= 2){
      stop("There are more than one variables with the name `variable`: please refit the model with different names for them.")
    }else if(sum(names(object$random_samp_indexes) == variable) == 0){
      stop("The specified variable cannot be found in the fitted model, please check the name.")
    }
    global_samps <- samps$samps[object$boundary_samp_indexes[[variable]], , drop = F]
    if(boundary.condition == "no"){
      global_samps <- NULL
    }
    coefsamps <- samps$samps[object$random_samp_indexes[[variable]], ]
    if(boundary.condition == "only"){
      coefsamps <- 0 * coefsamps
    }
    if (instance@smoothing_var == variable && class(instance) == "iwp") {
      IWP <- instance
      ## Step 2: Initialization
      if (is.null(newdata)) {
        refined_x_final <- seq(from = min(IWP@observed_x), to = max(IWP@observed_x), length.out = 3000)
      } else {
        refined_x_final <- sort(newdata[[variable]] - IWP@initial_location) # initialize according to `initial_location`
      }
      if(include.intercept){
        intercept_samps <- samps$samps[object$fixed_samp_indexes[["intercept"]], , drop = F]
      } else{
        intercept_samps <- NULL
      }
      ## Step 3: apply `compute_post_fun_iwp` to samps
      f <- compute_post_fun_iwp(
        samps = coefsamps, global_samps = global_samps,
        knots = IWP@knots,
        refined_x = refined_x_final,
        p = IWP@order,
        degree = deriv,
        intercept_samps = intercept_samps
      )
      f[,1] <- f[,1] + IWP@initial_location
    }
    else if (instance@smoothing_var == variable && class(instance) == "sgp") {
      sGP <- instance
      ## Step 2: Initialization
      if (is.null(newdata)) {
        refined_x_final <- seq(from = min(sGP@observed_x), to = max(sGP@observed_x), length.out = 3000)
      } else {
        refined_x_final <- sort(newdata[[variable]] - sGP@initial_location) # initialize according to `initial_location`
      }
      if(include.intercept){
        intercept_samps <- samps$samps[object$fixed_samp_indexes[["intercept"]], , drop = F]
      } else{
        intercept_samps <- NULL
      }
      ## Step 3: apply `compute_post_fun_sgp` to samps
      f <- compute_post_fun_sgp(
        samps = coefsamps, global_samps = global_samps,
        k = sGP@k,
        refined_x = refined_x_final,
        a = sGP@a, 
        m = sGP@m,
        region = (sGP@region - sGP@initial_location),
        intercept_samps = intercept_samps,
        boundary = sGP@boundary,
        initial_location = 0
      )
      f[,1] <- f[,1] + sGP@initial_location
    }
  }
  if(only.samples){
    names(f)[-1] <- paste0("samp", 1:(ncol(f)-1))
    return(f)
  }
  ## Step 4: summarize the prediction
  fpos <- extract_mean_interval_given_samps(f, quantiles = quantiles)
  names(fpos)[names(fpos) == "x"] <- variable
  return(fpos)
}

#' @export
plot.FitResult <- function(x, ...) {
  object <- x
  ### Step 1: predict with newdata = NULL
  for (instance in object$instances) { ## for each variable in model_fit
    if (class(instance) == "iwp") {
      predict_result <- predict(object, variable = as.character(instance@smoothing_var))
      matplot(
        x = predict_result[,1], y = predict_result[, c("q0.5", "q0.025", "q0.975")], lty = c(1, 2, 2), lwd = c(2, 1, 1),
        col = "black", type = "l",
        ylab = "effect", xlab = as.character(instance@smoothing_var)
      )
    }
    if (class(instance) == "sgp") {
      predict_result <- predict(object, variable = as.character(instance@smoothing_var))
      matplot(
        x = predict_result[,1], y = predict_result[, c("q0.5", "q0.025", "q0.975")], lty = c(1, 2, 2), lwd = c(2, 1, 1),
        col = "black", type = "l",
        ylab = "effect", xlab = as.character(instance@smoothing_var)
      )
    }
  }

  ### Step 2: then plot the original aghq object
  # plot(object$mod)
}


#' Extract the posterior samples from the fitted model for the target fixed variables.
#' 
#' @param model_fit The result from model_fit().
#' @param variables A vector of names of the target fixed variables to sample.
#' @export
sample_fixed_effect <- function(model_fit, variables){
  samps <- model_fit$samps$samps
  index <- model_fit$fixed_samp_indexes[variables]
  selected_samps <- t(samps[unlist(index), ,drop = FALSE])
  colnames(selected_samps) <- variables
  return(selected_samps)
} 


#' Computing the posterior samples of the function or its derivative using the posterior samples
#' of the basis coefficients for iwp
#'
#' @param samps A matrix that consists of posterior samples for the O-spline basis coefficients. Each column
#' represents a particular sample of coefficients, and each row is associated with one basis function. This can
#' be extracted using `sample_marginal` function from `aghq` package.
#' @param global_samps A matrix that consists of posterior samples for the global basis coefficients. If NULL,
#' assume there will be no global polynomials and the boundary conditions are exactly zero.
#' @param intercept_samps A matrix that consists of posterior samples for the intercept parameter. If NULL, assume
#' the function evaluated at zero is zero.
#' @param knots A vector of knots used to construct the O-spline basis, first knot should be viewed as "0",
#' the reference starting location. These k knots will define (k-1) basis function in total.
#' @param refined_x A vector of locations to evaluate the O-spline basis
#' @param p An integer value indicates the order of smoothness
#' @param degree The order of the derivative to take, if zero, implies to consider the function itself.
#' @return A data.frame that contains different samples of the function or its derivative, with the first column
#' being the locations of evaluations x = refined_x.
#' @examples
#' knots <- c(0, 0.2, 0.4, 0.6)
#' samps <- matrix(rnorm(n = (3 * 10)), ncol = 10)
#' result <- compute_post_fun_iwp(samps = samps, knots = knots, refined_x = seq(0, 1, by = 0.1), p = 2)
#' plot(result[, 2] ~ result$x, type = "l", ylim = c(-0.3, 0.3))
#' for (i in 1:9) {
#'   lines(result[, (i + 1)] ~ result$x, lty = "dashed", ylim = c(-0.1, 0.1))
#' }
#' global_samps <- matrix(rnorm(n = 10, sd = 0.1), ncol = 10)
#' result <- compute_post_fun_iwp(global_samps = global_samps, samps = samps, knots = knots, refined_x = seq(0, 1, by = 0.1), p = 2)
#' plot(result[, 2] ~ result$x, type = "l", ylim = c(-0.3, 0.3))
#' for (i in 1:9) {
#'   lines(result[, (i + 1)] ~ result$x, lty = "dashed", ylim = c(-0.1, 0.1))
#' }
#' @export
compute_post_fun_iwp <- function(samps, global_samps = NULL, knots, refined_x, p, degree = 0, intercept_samps = NULL) {
  if (p <= degree) {
    return(message("Error: The degree of derivative to compute is not defined. Should consider higher order smoothing model or lower order of the derivative degree."))
  }
  if (is.null(global_samps)) {
    global_samps <- matrix(0, nrow = (p - 1), ncol = ncol(samps))
  }
  if (nrow(global_samps) != (p - 1)) {
    return(message("Error: Incorrect dimension of global_samps. Check whether the choice of p is consistent with the fitted model."))
  }
  if (ncol(samps) != ncol(global_samps)) {
    return(message("Error: The numbers of posterior samples do not match between the O-splines and global polynomials."))
  }
  if (is.null(intercept_samps)) {
    intercept_samps <- matrix(0, nrow = 1, ncol = ncol(samps))
  }
  if (nrow(intercept_samps) != (1)) {
    return(message("Error: Incorrect dimension of intercept_samps."))
  }
  if (ncol(samps) != ncol(intercept_samps)) {
    return(message("Error: The numbers of posterior samples do not match between the O-splines and the intercept."))
  }

  ### Augment the global_samps to also consider the intercept
  global_samps <- rbind(intercept_samps, global_samps)

  ## Design matrix for the spline basis weights
  B <- dgTMatrix_wrapper(local_poly_helper(knots, refined_x = refined_x, p = (p - degree), neg_sign_order = degree))

  if ((p - degree) >= 1) {
    X <- global_poly_helper(refined_x, p = p)
    X <- as.matrix(X[, 1:(p - degree), drop = FALSE])
    for (i in 1:ncol(X)) {
      X[, i] <- (factorial(i + degree - 1) / factorial(i - 1)) * X[, i]
    }
    fitted_samps_deriv <- X %*% global_samps[(1 + degree):(p), , drop = FALSE] + B %*% samps
  } else {
    fitted_samps_deriv <- B %*% samps
  }
  result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_deriv)))
  result
}



#' Computing the posterior samples of the function using the posterior samples
#' of the basis coefficients for sGP
#'
#' @param samps A matrix that consists of posterior samples for the O-spline basis coefficients. Each column
#' represents a particular sample of coefficients, and each row is associated with one basis function. This can
#' be extracted using `sample_marginal` function from `aghq` package.
#' @param global_samps A matrix that consists of posterior samples for the global basis coefficients. If NULL,
#' assume there will be no global polynomials and the boundary conditions are exactly zero.
#' @param k The number of the sB basis.
#' @param region The region to define the sB basis
#' @param refined_x A vector of locations to evaluate the sB basis
#' @param a The frequency of sGP.
#' @param boundary Logical; whether the seasonal boundary basis is included.
#' @param m The number of harmonics to consider
#' @param intercept_samps Optional posterior samples for an intercept to add to each function draw.
#' @param initial_location The initial location of the sGP.
#' @return A data.frame that contains different samples of the function, with the first column
#' being the locations of evaluations x = refined_x.
#' @export
compute_post_fun_sgp <- function(samps, global_samps = NULL, k, refined_x, a, region, boundary = TRUE, m, intercept_samps = NULL, initial_location = NULL) {
  ## Design matrix for the spline basis weights
  B <- dgTMatrix_wrapper(Compute_B_sB_helper(refined_x = refined_x, k = k, a = a, region = region, boundary = boundary, initial_location = initial_location, m = m))
  X <- cbind(1,global_poly_helper_sgp(refined_x = refined_x, a = a, m = m, initial_location = initial_location))
  if (is.null(intercept_samps)) {
    intercept_samps <- matrix(0, nrow = 1, ncol = ncol(samps))
  }
  if(is.null(global_samps)){
    global_samps <- matrix(0, nrow = (2*m), ncol = ncol(samps))
  }
  global_samps <- rbind(intercept_samps, global_samps)
  f_samps <- X %*% global_samps + B %*% samps

  result <- cbind(x = refined_x, data.frame(as.matrix(f_samps)))
  result
}


#' Construct posterior inference given samples
#'
#' @param samps Posterior samples of f or its derivative, with the first column being evaluation
#' points x. This can be yielded by `compute_post_fun_iwp` function.
#' @param level The level to compute the pointwise interval. Ignored when quantiles are provided.
#' @param quantiles A numeric vector of quantiles to be computed.
#' @return A dataframe with a column for evaluation locations x, and posterior mean and pointwise
#' intervals at that set of locations.
#' @export
extract_mean_interval_given_samps <- function(samps, level = 0.95, quantiles = NULL) {
  x <- samps[, 1]
  samples <- samps[, -1]
  result <- data.frame(x = x)
  if(is.null(quantiles)){
    alpha <- 1 - level
    quantiles <- c((alpha / 2), 0.5, (level + (alpha / 2)))
  }
  for (q in quantiles) {
    result[[paste0("q",q)]] <- as.numeric(apply(samples, MARGIN = 1, quantile, p = q))
  }
  result$mean <- as.numeric(apply(samples, MARGIN = 1, mean))
  result
}




aghq_density_interpolation <- function(theta_marg){
  if(!is.null(theta_marg) && nrow(theta_marg) < 4){
    return("polynomial")
  }
  "spline"
}

#' Obtain the posterior density of a variance parameter in the fitted model
#' 
#' @param object The fitted object from the function `model_fit`.
#' @param component The component of the variance parameter that you want to show. By default this is `NULL`, indicating the family.sd is of interest.
#' @param h For PSD, the unit of predictive step to consider, by default is set to `NULL`, indicating the result is using the same `h` as in the model fitting.
#' @param theta_logprior The log prior function used on the selected variance parameter. By default is `NULL`, and the current Exponential prior will be used.
#' @param MCMC_samps_only For model fitted with MCMC, whether only the posterior samples are needed.
#' @export
var_density <- function(object, component = NULL, h = NULL, theta_logprior = NULL, MCMC_samps_only = FALSE){
  postsigma <- NULL
  if(is.null(component) && isTRUE(object$control.family$sd_known)){
    stop("The Gaussian observation SD is fixed in this fitted model, so there is no family SD posterior density.")
  }
  
  if(is.null(theta_logprior)){
    theta_logprior <- function(theta, prior_alpha, prior_u) {
      lambda <- -log(prior_alpha)/prior_u
      log(lambda/2) - lambda * exp(-theta/2) - theta/2
    }
  }
  priorfunc <- function(x, prior_alpha, prior_u) exp(theta_logprior(x, prior_alpha, prior_u))
  priorfuncsigma <- function(x, prior_alpha, prior_u) (2/x) * exp(theta_logprior(-2*log(x), prior_alpha, prior_u))
  
  if(any(class(object$mod) == "aghq")){
    if(is.null(component)){
      if(object$family != "gaussian"){
        stop("There is no family SD in the fitted model. Please indicate which component of the var-parameter that you want to show in `component`.")
      }
      theta_marg <- object$mod$marginals[[length(object$instances) + 1]]
      if(nrow(theta_marg) <= 2){
        stop("The number of quadrature points is too small, please use aghq_k >= 3.")
      }
      logpostsigma <- aghq::compute_pdf_and_cdf(theta_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = aghq_density_interpolation(theta_marg))
      postsigma <- data.frame(SD = logpostsigma$transparam, 
                                  post = logpostsigma$pdf_transparam,
                                  prior = priorfuncsigma(logpostsigma$transparam, prior_alpha = object$control.family$sd.prior$param$alpha, prior_u = object$control.family$sd.prior$param$u))
    }
    else{
      for (i in 1:length(object$instances)) {
        instance <- object$instances[[i]]
        if(instance@smoothing_var == component){
          theta_marg <- object$mod$marginals[[i]]
          if(nrow(theta_marg) <= 2){
            stop("The number of quadrature points is too small, please use aghq_k >= 3.")
          }
          logpostsigma <- aghq::compute_pdf_and_cdf(theta_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = aghq_density_interpolation(theta_marg))
          postsigma <- data.frame(SD = logpostsigma$transparam, 
                                  post = logpostsigma$pdf_transparam,
                                  prior = priorfuncsigma(logpostsigma$transparam, prior_alpha = object$instances[[i]]@sd.prior$param$alpha, prior_u = object$instances[[i]]@sd.prior$param$u))
          
          if(is.null(h)){
            if(!is.null(instance@sd.prior$h)){
              h <- instance@sd.prior$h
            }
          }
          if(!is.null(h)){
            if(class(instance) == "iwp"){
              p <- instance@order
              correction <- sqrt((h^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2)))
            }
            else if(class(instance) == "sgp"){
              correction <- 0
              for (j in 1:instance@m) {
                correction <- correction + compute_d_step_sgpsd(d = h, a = (j*instance@a))
              }
            }
            else{
              stop("PSD is currently on defined on iwp and sGP, please specify h = NULL for other type of random effect")
            }
            postsigmaPSD <- data.frame(PSD = postsigma$SD * correction, 
                                       post.PSD = postsigma$post / correction,
                                       prior.PSD = postsigma$prior / correction)
            
            postsigma <- cbind(postsigma, postsigmaPSD)
          }
          
        }
      }
      
    }
  }
  
  else if(any(class(object$mod) == "stanfit")){
    sigmaPSD_marg_samps <- NULL
    if(is.null(component)){
      if(object$family != "gaussian"){
        stop("There is no family SD in the fitted model. Please indicate which component of the var-parameter that you want to show in `component`.")
      }
      theta_marg_samps <- object$samps$thet[[length(object$instances) + 1]]
      sigma_marg_samps <- exp(-0.5*theta_marg_samps)
      sigma_marg_density <- density(sigma_marg_samps)
      postsigma <- data.frame(SD = sigma_marg_density$x, 
                              post = sigma_marg_density$y,
                              prior = priorfuncsigma(sigma_marg_density$x, prior_alpha = object$control.family$sd.prior$param$alpha, prior_u = object$control.family$sd.prior$param$u))
    }
    else{
      for (i in 1:length(object$instances)) {
        instance <- object$instances[[i]]
        if(instance@smoothing_var == component){
          theta_marg_samps <- object$samps$thet[[i]]
          sigma_marg_samps <- exp(-0.5*theta_marg_samps)
          sigma_marg_density <- density(sigma_marg_samps)
          postsigma <- data.frame(SD = sigma_marg_density$x, 
                                  post = sigma_marg_density$y,
                                  prior = priorfuncsigma(sigma_marg_density$x, prior_alpha = object$instances[[i]]@sd.prior$param$alpha, prior_u = object$instances[[i]]@sd.prior$param$u))
          if(is.null(h)){
            if(!is.null(instance@sd.prior$h)){
              h <- instance@sd.prior$h
            }
          }
          if(!is.null(h)){
            if(class(instance) == "iwp"){
              p <- instance@order
              correction <- sqrt((h^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2)))
            }
            else if(class(instance) == "sgp"){
              correction <- 0
              for (j in 1:instance@m) {
                correction <- correction + compute_d_step_sgpsd(d = h, a = (j*instance@a))
              }
            }
            else{
              stop("PSD is currently on defined on iwp and sGP, please specify h = NULL for other type of random effect")
            }
            sigmaPSD_marg_samps <- sigma_marg_samps * correction
            postsigmaPSD <- data.frame(PSD = postsigma$SD * correction, 
                                       post.PSD = postsigma$post / correction,
                                       prior.PSD = postsigma$prior / correction)
            
            postsigma <- cbind(postsigma, postsigmaPSD)
          }
          
        }
      }
    }
    if(MCMC_samps_only == TRUE){
      return(list(sigmaPSD_marg_samps = sigmaPSD_marg_samps, sigma_marg_samps = sigma_marg_samps))
    }
  }
  
  else{
    stop("The function `var_density` currently only supports model fittd with `method = aghq`.")
  }
  postsigma <- postsigma[order(postsigma$SD),]
  return(postsigma)
}

#' Plot the posterior density of a variance parameter in the fitted model
#' 
#' @param object The fitted object from the function `model_fit`.
#' @param component The component of the variance parameter that you want to show. By default this is `NULL`, indicating the family.sd is of interest.
#' @param h For PSD, the unit of predictive step to consider, by default is set to `NULL`, indicating the result is using the same `h` as in the model fitting.
#' @param theta_logprior The log prior function used on the selected variance parameter. By default is `NULL`, and the current Exponential prior will be used.
#' @export
var_plot <- function(object, component = NULL, h = NULL, theta_logprior = NULL){
  hyper_result <- var_density(object = object, component = component, h = h, theta_logprior = theta_logprior, MCMC_samps_only = FALSE)
  if("PSD" %in% names(hyper_result)){
    matplot(hyper_result[,'PSD'], hyper_result[,c('post.PSD','prior.PSD')], lty=c(1,2), type = "l", xlab = "PSD", ylab = "Density")
    legend('topright', lty=c(1,2), col=c('black','red'), legend=c('post','prior'), bty = "n")
  }
  else{
    matplot(hyper_result[,'SD'], hyper_result[,c('post','prior')], lty=c(1,2), type = "l", xlab = "SD", ylab = "Density")
    legend('topright', lty=c(1,2), col=c('black','red'), legend=c('post','prior'), bty = "n")
  }
}


#' Obtain the posterior and prior density of all the parameters in the fitted model
#' 
#' @param object The fitted object from the function `model_fit`.
#' @export
para_density <- function(object){
  result_list <- list()
  for (fixed_name in names(object$fixed_samp_indexes)) {
    samps <- sample_fixed_effect(object, variables = fixed_name)
    fixed_density <- density(samps)
    result_list[[fixed_name]] <- data.frame(effect = fixed_density$x, post = fixed_density$y)
  }
  
  for (random_name in names(object$random_samp_indexes)) {
    result_list[[random_name]] <- var_density(object = object, component = random_name)
  }
  
  if(object$family == "gaussian" && !isTRUE(object$control.family$sd_known)){
    result_list[["family_sd"]] <- var_density(object = object)
  }
  
  return(result_list)
}


#' Obtain the posterior summary table for all the parameters in the fitted model
#' @param object The fitted object from the function `model_fit`.
#' @param quantiles The specified quantile to display the posterior summary, default is c(0.025, 0.975).
#' @param digits The significant digits to be kept in the result, default is 3.
#' @export
post_table <- function(object, quantiles = c(0.025, 0.975), digits = 3){
  all_density <- para_density(object)
  all_cdf <- list()
  compute_cdf <- function(x, y){
    new_data_frame <- list(x = x)
    new_data_frame$cdf <- cumsum(y * c(diff(x), 0))
    new_data_frame
  }
  result_table <- c("name", "median", paste0("q", unlist(strsplit(toString(quantiles), split = ", "))), "prior", "prior:P1", "prior:P2")
  for (name in names(object$fixed_samp_indexes)) {
    all_cdf[[name]] <- compute_cdf(x = all_density[[name]]$effect, y = all_density[[name]]$post)
    to_add <- c(name, all_cdf[[name]]$x[max(which(all_cdf[[name]]$cdf <= 0.5))])
    for (q in quantiles) {
      to_add <- c(to_add, all_cdf[[name]]$x[max(which(all_cdf[[name]]$cdf <= q))])
    }
    to_add <- c(to_add, "Normal", object$control.fixed[[name]]$mean, 1/object$control.fixed[[name]]$prec)
    result_table <- rbind(result_table, to_add)
  } 
  
  for (name in names(object$random_samp_indexes)) {
    if("PSD" %in% names(all_density[[name]])){
      all_cdf[[name]] <- compute_cdf(x = all_density[[name]]$PSD, y = all_density[[name]]$post.PSD)
      print_name <- paste0(name, " (PSD)")
    }
    else{
      all_cdf[[name]] <- compute_cdf(x = all_density[[name]]$SD, y = all_density[[name]]$post)
      print_name <- paste0(name, " (SD)")
    }
    to_add <- c(print_name, all_cdf[[name]]$x[max(which(all_cdf[[name]]$cdf <= 0.5))])
    for (q in quantiles) {
      to_add <- c(to_add, all_cdf[[name]]$x[max(which(all_cdf[[name]]$cdf <= q))])
    }
    for (i in 1:length(object$instances)) {
      if(name == object$instances[[i]]@smoothing_var){
        if(!is.null(object$instances[[i]]@sd.prior$h) ||!is.null(object$instances[[i]]@sd.prior$step)){
          param <- object$instances[[i]]@psd.prior$param
        }
        else{
          param <- object$instances[[i]]@sd.prior$param
        }
      }
    }
    object$instances
    to_add <- c(to_add, "Exponential", param$u, param$alpha)
    result_table <- rbind(result_table, to_add)
  }
  
  if("family_sd" %in% names(all_density)){
    all_cdf[["family_sd"]] <- compute_cdf(x = all_density[["family_sd"]]$SD, y = all_density[["family_sd"]]$post)
    to_add <- c("family_sd", all_cdf[["family_sd"]]$x[max(which(all_cdf[["family_sd"]]$cdf <= 0.5))])
    for (q in quantiles) {
      to_add <- c(to_add, all_cdf[["family_sd"]]$x[max(which(all_cdf[["family_sd"]]$cdf <= q))])
    }
    to_add <- c(to_add, "Exponential", object$control.family$sd.prior$param$u, object$control.family$sd.prior$param$alpha)
    result_table <- rbind(result_table, to_add)
  }
  result_table <- data.frame(unname(result_table))
  result_table_names <- unlist(unname(result_table)[1,])
  result_table <- (type.convert(result_table[-1,], as.is = TRUE))
  result_table <- data.frame(lapply(result_table, function(y) if(is.numeric(y)) round(y, digits) else y)) 
  colnames(result_table) <- result_table_names
  result_table
}
