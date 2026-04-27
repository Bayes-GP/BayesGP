get_result_by_method <- function(response_var, data, instances, design_mat_fixed, family, control.family, control.fixed, fixed_effects, fixed_effects_names, offset_sum, aghq_k, size, cens, weight, strata, method, M, customized_template, option_list, envir = parent.frame(), extra_theta_num) {
  family <- tolower(family)
  if(is.null(customized_template)){
    cpp = "BayesGP"
  }
  else{
    cpp = customized_template
  }
  # Family types: gaussian - 0, poisson - 1, binomial - 2, coxph - 3, casecrossover - 4, prior.only - -2.
  if (family == "gaussian") {
    family_type <- 0
  } else if (family == "poisson") {
    family_type <- 1
  } else if (family == "binomial") {
    family_type <- 2
  } else if (family == "coxph") {
    family_type <- 3
  } else if(family == "casecrossover" || family == "cc"){
    family_type <- 4
  }else if(family == "customized"){
    if(is.null(customized_template)){
      stop("In order to use the customized family: please input the name of the complied cpp template as `customized_template`.")
    }
    family_type <- -1
  }
  else if(family == "none"){
    message("The family option is set to `none`, and only prior samples will be produced.")
    family_type <- -2
  }
  
  # Containers for random effects
  X <- list()
  B <- list()
  P <- list()
  logPdet <- list()
  u <- list()
  alpha <- list()
  betaprec <- list()
  betamean <- list()
  
  # Containers for fixed effects
  beta_fixed_prec <- list()
  beta_fixed_mean <- list()
  Xf <- list()
  
  w_count <- 0
  family_sd_known <- identical(family_type, 0) && isTRUE(control.family$sd_known)
  family_sd_value <- if(family_sd_known) as.numeric(control.family$sd) else NA_real_
  # Need a theta for the Gaussian variance, so
  # theta_count starts at 1 if Gaussian
  theta_count <- 0 + (family_type == 0 && !family_sd_known)
  
  for (instance in instances) {
    # For each random effects
    if (class(instance) == "iwp") {
      X[[length(X) + 1]] <- dgTMatrix_wrapper(instance@X)
      if(instance@order != 1){
      for (jj in 1:length(instance@boundary.prior$prec)) {
        betaprec[[length(betaprec) + 1]] <- instance@boundary.prior$prec[jj]
        betamean[[length(betamean) + 1]] <- instance@boundary.prior$mean[jj]
        if(instance@boundary.prior$prec[jj] == Inf){
          betaprec[[length(betaprec)]] = 100000
          X[[length(X)]][,jj] <- 0
        }
      }
      }
      w_count <- w_count + ncol(instance@X)
    }
    else if(class(instance) == "sgp"){
      X[[length(X) + 1]] <- dgTMatrix_wrapper(instance@X)
      for (jj in 1:length(instance@boundary.prior$prec)) {
        betaprec[[length(betaprec) + 1]] <- instance@boundary.prior$prec[jj]
        betamean[[length(betamean) + 1]] <- instance@boundary.prior$mean[jj]
        if(instance@boundary.prior$prec[jj] == Inf){
          betaprec[[length(betaprec)]] = 100000
          X[[length(X)]][,jj] <- 0
        }
      }
      w_count <- w_count + ncol(instance@X)
    }
    B[[length(B) + 1]] <- dgTMatrix_wrapper(instance@B)
    P[[length(P) + 1]] <- dgTMatrix_wrapper(instance@P)
    logPdet[[length(logPdet) + 1]] <- as.numeric(determinant(instance@P, logarithm = TRUE)$modulus)
    u[[length(u) + 1]] <- instance@sd.prior$param$u
    alpha[[length(alpha) + 1]] <- instance@sd.prior$param$alpha
    w_count <- w_count + ncol(instance@B)
    theta_count <- theta_count + 1
  }
  
  # For the variance of the Gaussian family
  # From control.family, if applicable
  if (family_type == 0 && !family_sd_known) {
    u[[length(u) + 1]] <- control.family$sd.prior$param$u
    alpha[[length(alpha) + 1]] <- control.family$sd.prior$param$alpha
  }
  if(family_type == 3 || family_type == 4){
    if(length(design_mat_fixed) >= 1){
      for (i in 1:length(design_mat_fixed)) {
        # For each fixed effects
        beta_fixed_prec[[i]] <- control.fixed[[fixed_effects_names[[i]]]]$prec
        beta_fixed_mean[[i]] <- control.fixed[[fixed_effects_names[[i]]]]$mean
        Xf[[length(Xf) + 1]] <- dgTMatrix_wrapper(design_mat_fixed[[i]])
        w_count <- w_count + ncol(design_mat_fixed[[i]])
      }
    }
    
  }
  
  else{
    for (i in 1:length(design_mat_fixed)) {
      # For each fixed effects
      if (i == 1) {
        beta_fixed_prec[[i]] <- control.fixed$intercept$prec
        beta_fixed_mean[[i]] <- control.fixed$intercept$mean
        
      } else {
        beta_fixed_prec[[i]] <- control.fixed[[fixed_effects_names[[i - 1]]]]$prec
        beta_fixed_mean[[i]] <- control.fixed[[fixed_effects_names[[i - 1]]]]$mean
      }
      Xf[[length(Xf) + 1]] <- dgTMatrix_wrapper(design_mat_fixed[[i]])
      w_count <- w_count + ncol(design_mat_fixed[[i]])
    }
  }
  
  tmbdat <- list(
    # For Random effects
    X = X,
    B = B,
    P = P,
    logPdet = logPdet,
    u = u,
    alpha = alpha,
    betaprec = betaprec,
    betamean = betamean,
    
    # For Fixed Effects:
    beta_fixed_prec = beta_fixed_prec,
    beta_fixed_mean = beta_fixed_mean,
    Xf = Xf,
    
    # Response
    y = data[[response_var]],
    offset_sum = as.numeric(offset_sum),
    
    # Family type
    family_type = family_type,
    gaussian_sd_known = as.numeric(family_sd_known),
    gaussian_sd_value = if(is.na(family_sd_value)) 0 else family_sd_value
  )
  
  if(is.null(size)){
    tmbdat$size <- numeric(length = length(tmbdat$y)) + 1 # A vector of 1s being default
  } else if(is.null(data[[size]])){
    tmbdat$size <- numeric(length = length(tmbdat$y)) + 1 # A vector of 1s being default
  } else{
    tmbdat$size <- data[[size]]
  }
  
  # If Family == "coxph", check whether cens is defined in user's input
  if (family_type == 3) {
    tmbdat$ranks = rank(tmbdat$y, ties.method = "min")
    n <- length(tmbdat$y)
    tmbdat$D <- cbind(Matrix::Matrix(1,n-1,1),Matrix::Diagonal(n-1,-1))
    if(is.null(data[[cens]])){
      tmbdat$cens <- numeric(length = length(tmbdat$y)) + 1 # A vector of 1s being default
    }
    else{
      tmbdat$cens <- data[[cens]]
    }
  }
  
  # If Family == "cc", build grouped multinomial likelihood data.
  if (family_type == 4) {
    if(is.null(weight)){
      count <- data[[response_var]]
    }
    else{
      count <- data[[weight]]
    }
    if(is.null(data[[strata]])){
      stop("The specified names for strata are not correct.")
    }
    else{
      strata <- data[[strata]]
    }

    if(!is.numeric(count) || length(count) != nrow(data)){
      stop("Case-crossover counts must be a numeric vector with one value per row.")
    }
    if(any(!is.finite(count))){
      stop("Case-crossover counts must be finite.")
    }
    if(any(count < 0)){
      stop("Case-crossover counts must be nonnegative.")
    }
    if(any(abs(count - round(count)) > sqrt(.Machine$double.eps))){
      stop("Case-crossover counts must be whole numbers.")
    }

    count <- as.numeric(round(count))
    tmbdat$y <- count

    unique_strata <- unique(strata)
    max_N <- max(sapply(unique_strata, function(s) sum(strata == s)))
    stratum_members <- matrix(0L, nrow = length(unique_strata), ncol = max_N)
    stratum_count <- numeric(length(unique_strata))

    for (i in seq_along(unique_strata)) {
      rows <- which(strata == unique_strata[[i]])
      stratum_members[i, seq_along(rows)] <- rows
      stratum_count[i] <- sum(count[rows])
    }

    tmbdat$stratum_members <- stratum_members
    tmbdat$stratum_count <- stratum_count
  }
  
  tmbparams <- list(
    W = c(rep(0, w_count)), # recall W is everything in the model (RE or FE)
    theta = c(rep(0, theta_count))
  )
  
  if(!is.null(customized_template)){
    if(is.null(extra_theta_num)){
      num_theta_to_add = 0
    }
    else if(is.numeric(extra_theta_num)){
      num_theta_to_add = extra_theta_num
    }
    tmbparams$theta <- c(tmbparams$theta,rep(0,num_theta_to_add))
    theta_count <- length(tmbparams$theta)
  }
  
  if(theta_count == 0 & method != "nlminb"){
    stop("For model with no hyper-parameter, the method cannot be aghq or MCMC.")
  }
  
  if(method == "nlminb"){
    if(theta_count != 0){
      stop("For model with hyper-parameter, the method should be aghq or MCMC.")
    }
    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      DLL = cpp,
      silent = TRUE
    )
    ff$he <- function(w) numDeriv::jacobian(ff$gr, w)
    
    opt <- nlminb(start = ff$par, objective = ff$fn, gradient = ff$gr, hessian = ff$he, 
                  control = list(eval.max = 20000, iter.max = 20000))
    prec_matrix <- Matrix::forceSymmetric(ff$he(opt$par))
    mod = list(mean = opt$par, prec = as.matrix(prec_matrix), opt = opt)
    class(mod) <- method
  }
  else if(method == "aghq"){
    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      random = "W",
      DLL = cpp,
      silent = TRUE
    )
    ff$he <- function(w) numDeriv::jacobian(ff$gr, w)
    mod <- aghq::marginal_laplace_tmb(ff, aghq_k, c(rep(0, theta_count))) # The default value of aghq_k is 4
  }
  else if(method == "MCMC"){
    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      random = "W",
      DLL = cpp,
      silent = TRUE
    )
    ff$he <- function(w) numDeriv::jacobian(ff$gr, w)
    default_option_list <- get_default_option_list_MCMC(option_list = option_list)
    if(default_option_list$silent == TRUE){
      sink(tempfile())
      mod <- tmbstan::tmbstan(
        ff,
        chains = default_option_list$chains,
        cores = default_option_list$cores,
        iter = default_option_list$warmup + M,
        warmup = default_option_list$warmup,
        seed = default_option_list$seed,
        silent = default_option_list$silent,
        laplace = default_option_list$laplace
      )
      sink()
      
    }
    else{
      mod <- tmbstan::tmbstan(
        ff,
        chains = default_option_list$chains,
        cores = default_option_list$cores,
        iter = default_option_list$warmup + M,
        warmup = default_option_list$warmup,
        seed = default_option_list$seed,
        silent = default_option_list$silent,
        laplace = default_option_list$laplace
      )
    }

  }
  return(list(
    mod = mod,
    w_count = w_count,
    family_sd_known = family_sd_known,
    family_sd_value = family_sd_value
  ))
}


#' @title
#' Model fitting with random effects/fixed effects
#'
#' @description
#' Fitting a hierarchical model based on the provided formula, data and parameters such as type of method and family of response.
#' Returning the S4 objects for the random effects, concatenated design matrix for the intercepts and fixed effects, fitted model,
#' indexes to partition the posterior samples.
#'
#' @param formula A formula that contains one response variable, and covariates with either random or fixed effect.
#' @param data A dataframe that contains the response variable and other covariates mentioned in the formula.
#' @param method The inference method used in the model. By default, the method is set to be "aghq".
#' @param family The family of response used in the model. By default, the family is set to be "gaussian".
#' @param control.family Parameters used to specify the family-level behavior.
#'   For Gaussian models, `control.family = list(sd.prior = ...)` keeps the
#'   prior-based observation SD pathway, while `control.family = list(sd = value)`
#'   fixes the Gaussian observation SD at a known value.
#' @param control.fixed Parameters used to specify the priors for the fixed effects. For example control.fixed = list(intercept = list(prec = 0.001, mean = 0)) will setup the prior N(0,1/0.001) for the intercept parameter. When left unspecified, all fixed effect parameters will be assigned independent N(0,1/0.001) priors. 
#' @param aghq_k Number of quadrature points for adaptive Gauss-Hermite quadrature when `method = "aghq"`.
#' @param size The name of the size variable, should be one of the variables in `data`. The default value is "NULL", corresponding to a vector of 1s. This is only used for the Binomial family, and denotes the number of binomial trails.
#' @param cens The name of the right-censoring indicator, should be one of the variables in `data`. The default value is "NULL". This is only used for the CoxPH family.
#' @param weight Optional name of a count/weight variable in `data`. For case-crossover models, this supplies multinomial counts within strata; when omitted, the response variable is used as the count.
#' @param strata Name of the strata variable in `data`. Required for case-crossover models.
#' @param M The number of posterior samples to be taken, by default is 3000.
#' @param customized_template The name of the customized cpp template that the user wants to use instead. By default this is NULL, and the cpp template `BayesGP` will be used.
#' @param Customized_RE The list that contains the compute_B and compute_P functions for the customized random effect. By default, this is NULL and there is not customized random effect in the model.
#' @param option_list A list that controls the details of the inference algorithm, by default is an empty list.
#' @param envir The environment in which the formula and other expressions are to be evaluated. 
#'   Defaults to `parent.frame()`, which refers to the environment from which the function was called.
#'   This allows the function to access variables that are defined in the calling function's scope.
#' @param extra_theta_num An integer number to indicate the extra number of parameters required in the 'theta' vector. Should only be specified when using customized template.
#' @details
#'   BayesGP supports smooth terms through `f()`. Near-monotone `mgp` and
#'   `tiwp2` terms default to
#'   `computation = "fem"` and can use `computation = "state-space"` when an
#'   exact support-grid representation is desired. An exact state-space
#'   representation is also available for native `iwp` terms when `order = 2`.
#'   Exact state-space fits can only be predicted on the support grid supplied
#'   during fitting; use `grid = ...` inside `f()` when predictions are needed
#'   at additional locations. For near-monotone FEM terms, use
#'   `region = c(left, right)` inside `f()` to set the continuous basis domain
#'   for prediction outside the observed range. For near-monotone `mgp` and
#'   `tiwp2` terms, `c` is specified on the original input scale in the base
#'   model `m(x) = (x + c)^((a - 1) / a)` or `m(x) = log(x + c)` when `a = 1`;
#'   any centering around `initial_location` is handled internally. If `c` is
#'   omitted, BayesGP chooses a domain-aware default that makes `x + c`
#'   strictly positive over the term domain.
#'
#' @return A list that contains following items: the S4 objects for the random
#'   effects (instances), concatenated design matrix for the fixed effects
#'   (design_mat_fixed), fitted aghq (mod) and indexes to partition the
#'   posterior samples (boundary_samp_indexes, random_samp_indexes and
#'   fixed_samp_indexes).
#'
#' @export
model_fit <- function(formula, data, method = "aghq", family = "gaussian", control.family, control.fixed, aghq_k = 4, size = NULL, cens = NULL, weight = NULL, strata = NULL, M = 3000, customized_template = NULL, Customized_RE = NULL, option_list = list(), envir = parent.frame(), extra_theta_num = NULL) {
  # parse the input formula
  parse_result <- parse_formula(formula)
  response_var <- parse_result$response
  rand_effects <- parse_result$rand_effects
  fixed_effects <- parse_result$fixed_effects
  offset_effects <- parse_result$offset_effects
  
  instances <- list()
  near_mono_meta <- list()
  exact_iwp_meta <- list()
  design_mat_fixed <- list()
  family = tolower(family)
  family_is_coxph <- FALSE
  if(family == "cox" || family == "coxph"){
    family_is_coxph <- TRUE
    data <- data[order(data[[response_var]]), ]
  }
  
  family_is_cc <- FALSE
  if(family == "casecrossover" || family == "cc"){
    family_is_cc <- TRUE
  }
  control.family <- normalize_family_control(control.family, family = family)

  # For random effects
  for (rand_effect in rand_effects) {
    smoothing_var <- rand_effect$smoothing_var
    if(is.null(smoothing_var)){
      smoothing_var <- rand_effect$x
      if(is.null(smoothing_var)){
        if(toString(rand_effect[[2]]) %in% names(data)){
          smoothing_var <- rand_effect[[2]]
        }
        else{
          stop("Error: Cannot find the specified the covaraite name. The covariate name should be specified as smoothing_var or x, or the first argument in f().")
        }
      }
    }
    model_class <- normalize_smooth_term_model(
      extract_term_argument(rand_effect, "model", envir = envir, default = "iid")
    )

    if(model_class %in% c("mgp", "tiwp2")){
      built_term <- build_nearmono_instance(
        rand_effect = rand_effect,
        response_var = response_var,
        data = data,
        envir = envir
      )
      instances[[length(instances) + 1]] <- built_term$instance
      near_mono_meta[[built_term$meta$component]] <- built_term$meta
      next
    }

    if(model_class == "iwp"){
      requested_iwp_method <- resolve_computation_method(
        computation = extract_term_argument(rand_effect, "computation", envir = envir, default = NULL),
        method = extract_term_argument(rand_effect, "method", envir = envir, default = NULL),
        default = NULL,
        allowed = c("state-space", "fem"),
        label = "IWP smooth terms"
      )

      if(identical(requested_iwp_method, "state-space")){
        built_term <- build_iwp_state_space_instance(
          rand_effect = rand_effect,
          response_var = response_var,
          data = data,
          envir = envir
        )
        instances[[length(instances) + 1]] <- built_term$instance
        exact_iwp_meta[[built_term$meta$component]] <- built_term$meta
        next
      }
    }
    
    sd.prior <- eval(rand_effect$sd.prior, envir = envir)
    if (is.null(sd.prior)) {
      sd.prior <- eval(rand_effect$prior, envir = envir)
      if(is.null(sd.prior)){
        sd.prior <- list(prior = "exp", param = list(u = 1, alpha = 0.5))
      }
    }
    if (length(sd.prior) == 1){
      if(is.numeric(sd.prior)){
        sd.prior <- list(prior = "exp", param = list(u = as.numeric(sd.prior), alpha = 0.5))
      }
    }
    if(!"prior" %in% names(sd.prior)){
      sd.prior$prior <- "exp"
    }
    if(!"param" %in% names(sd.prior)){
      stop("If sd.prior is provided as a list, it must contains a list called param.")
    }else{
      if(length(sd.prior$param) == 1){
        sd.prior$param <- list(u = sd.prior$param[[1]], alpha = 0.5)
      }
      else{
        sd.prior$param <- list(u = sd.prior$param$u, alpha = sd.prior$param$alpha)
        if(is.null(sd.prior$param$alpha)){
          warning("The value of alpha is not provided in sd.prior$param: automatically filled with 0.5.")
          sd.prior$param$alpha <- 0.5
        }
        if(is.null(sd.prior$param$u)){
          stop("Error: The value of u is not provided in sd.prior$param.")
        }
      }
    }
    
    if (sd.prior$prior != "exp" & sd.prior$prior != "Exp" & sd.prior$prior != "exponential" & sd.prior$prior != "Exponential" & sd.prior$prior != "customized") {
      stop("Error: For each random effect, sd.prior currently only supports 'exp' (exponential), or 'customized' as prior.")
    }
    if(sd.prior$param$alpha > 1 | sd.prior$param$alpha < 0){
      if(sd.prior$prior != "customized"){
        stop("Error: The value of sd.prior$param$alpha is not specified as a probability.")
      }
    }
    
    if (model_class == "iwp") {
      order <- eval(rand_effect$order, envir = envir)
      knots <- eval(rand_effect$knots, envir = envir)
      if("k" %in% names(rand_effect)){
        k <- eval(rand_effect$k, envir = envir)
      }
      else if("knots" %in% names(rand_effect)){
        k <- length(knots)
      }
      else{
        k <- NULL
      }
      initial_location <- if("initial_location" %in% names(rand_effect)){
        eval(rand_effect$initial_location, envir = envir)[1]
      } else {
        NULL
      }
      if (!(is.null(k)) && k < 3) {
        stop("Error: parameter <k> in the random effect part should be >= 3.")
      }
      if (is.null(order) | order < 1) {
        stop("Error: Parameter <order> in the random effect part should be >= 1.")
      }
      boundary.prior <- eval(rand_effect$boundary.prior, envir = envir)
      
      if("region" %in% names(rand_effect)){
        region <- eval(rand_effect$region, envir = envir)
      }
      else{
        region <- range(data[[smoothing_var]])
      }
      
      # If the user does not specify initial_location, compute initial_location with
      # the mean of region
      if (is.null(initial_location)) {
        initial_location <- mean(region)
      }
      if (!is.numeric(initial_location)){
        if (initial_location == "middle") {
          initial_location <- mean(region)
        }
        else if (initial_location == "left") {
          initial_location <- min(region)
        }
        else if (initial_location == "right") {
          initial_location <- max(region)
        }
        else{
          stop("initial_location should be either numeric or one of 'left', 'right', 'middle'.")
        }
      }
      # If the user does not specify knots, compute knots with
      # the parameter k
      initialized_smoothing_var <- data[[smoothing_var]] - initial_location
      region_initialized <- region - initial_location
      if (is.null(knots)) {
        default_k <- 30
        if (is.null(k)) {
          knots <- unique(sort(seq(from = min(region_initialized), to = max(region_initialized), length.out = default_k))) # should be length.out
        } else {
          knots <- unique(sort(seq(from = min(region_initialized), to = max(region_initialized), length.out = k)))
        }
      }
      else{
        knots <- knots - initial_location
      }
      observed_x <- sort(initialized_smoothing_var) # initialized_smoothing_var: initialized observed covariate values
      if (is.null(boundary.prior)) {
        if(order != 1){
          boundary.prior <- list(prec = 0.001, mean = 0)
        }
        else{
          boundary.prior <- list()
        }
      }
      if(is.null(boundary.prior$prec)){
        boundary.prior$prec <- 0.001
      }
      if(is.null(boundary.prior$mean)){
        boundary.prior$mean <- 0
      }
      if(length(boundary.prior$mean) != (order-1)){
        boundary.prior$mean <- rep(boundary.prior$mean[1], length.out = (order-1))
      }
      if(length(boundary.prior$prec) != (order-1)){
        boundary.prior$prec <- rep(boundary.prior$prec[1], length.out = (order-1))
      }
      instance <- new(model_class,
        response_var = response_var,
        smoothing_var = smoothing_var, order = order,
        knots = knots, k = length(knots),
        observed_x = observed_x, sd.prior = sd.prior, boundary.prior = boundary.prior, data = data
      )
      # Case for iwp
      instance@initial_location <- initial_location
      instance@region <- region
      instance@X <- global_poly(instance)[, -1, drop = FALSE]
      instance@B <- local_poly(instance)
      instance@P <- compute_weights_precision(instance)
      if(is.numeric(sd.prior$h)){
        instance@psd.prior <- instance@sd.prior
        instance@sd.prior$param <- prior_conversion_iwp(d = sd.prior$h, prior = sd.prior$param, p = eval(rand_effect$order, envir = envir))
      } else if(is.numeric(sd.prior$step)){
        instance@sd.prior$h <- sd.prior$step
        instance@psd.prior <- instance@sd.prior
        instance@sd.prior$param <- prior_conversion_iwp(d = sd.prior$step, prior = sd.prior$param, p = eval(rand_effect$order, envir = envir))
      }
      instances[[length(instances) + 1]] <- instance
    } 
    else if (model_class == "iid") {
      instance <- new(model_class,
        response_var = response_var,
        smoothing_var = smoothing_var, sd.prior = sd.prior, data = data
      )
      # Case for iid
      instance@B <- compute_B(instance)
      instance@P <- compute_P(instance)
      instances[[length(instances) + 1]] <- instance
    }
    else if (model_class == "customized") {
      instance <- new(model_class,
                      response_var = response_var,
                      smoothing_var = smoothing_var, sd.prior = sd.prior, data = data,
                      compute_B = Customized_RE$compute_B, compute_P = Customized_RE$compute_P
      )
      # Case for customized
      instance@B <- compute_B(instance)
      instance@P <- compute_P(instance)
      instances[[length(instances) + 1]] <- instance
    }
    else if (model_class == "sgp"){
      if(!"a" %in% names(rand_effect)){
        if("freq" %in% names(rand_effect)){
          rand_effect$a <- (2*pi)*eval(rand_effect$freq, envir = envir)
        }
        else if("period" %in% names(rand_effect)){
          rand_effect$a <- (2*pi)/eval(rand_effect$period, envir = envir)
        }
      }
      a <- eval(rand_effect$a, envir = envir)
      k <- eval(rand_effect$k, envir = envir)
      # initial_location for sGP should be "left"
      initial_location <- "left"
      if("m" %in% names(rand_effect)){
        m <- eval(rand_effect$m, envir = envir)
      }
      else{
        m <- 1
      }
      if(is.null(k)){
        k <- 30
      }else if (k < 3) {
        stop("Error: parameter <k> in the random effect part should be >= 3.")
      }
      if (a < 0) {
        stop("Error: Parameter <a> in the random effect part should be positive.")
      }
      boundary.prior <- eval(rand_effect$boundary.prior, envir = envir)
      if("region" %in% names(rand_effect)){
        region <- eval(rand_effect$region, envir = envir)
      }
      else{
        region <- range(data[[smoothing_var]])
      }
      # If the user does not specify initial_location, compute initial_location with
      # the min of region
      if (is.null(initial_location)) {
        initial_location <- min(region)
      }
      if (!is.numeric(initial_location)){
        if (initial_location == "left") {
          initial_location <- min(region)
        }
        else{
          stop("initial_location should be either numeric or one of 'left', 'right', 'middle'.")
        }
      }
      initialized_smoothing_var <- data[[smoothing_var]] - initial_location
      observed_x <- sort(initialized_smoothing_var) # initialized_smoothing_var: initialized observed covariate values
      if("accuracy" %in% names(rand_effect)){
        accuracy <- eval(rand_effect$accuracy, envir = envir)
      }
      else{
        accuracy <- 5000
      }
      if (is.null(boundary.prior)) {
        boundary.prior <- list(prec = 0.001, mean = 0)
      }
      if(is.null(boundary.prior$prec)){
        boundary.prior$prec <- 0.001
      }
      if(is.null(boundary.prior$mean)){
        boundary.prior$mean <- 0
      }
      if(length(boundary.prior$mean) != 2){
        boundary.prior$mean <- rep(boundary.prior$mean[1], length.out = (2*m))
      }
      if(length(boundary.prior$prec) != 2){
        boundary.prior$prec <- rep(boundary.prior$prec[1], length.out = (2*m))
      }
      boundary <- TRUE
      if ("boundary" %in% names(rand_effect)){
        boundary <- eval(rand_effect$boundary, envir = envir)
      }
      instance <- new(model_class,
                      response_var = response_var,
                      smoothing_var = smoothing_var, a = a, m = m,
                      k = k, observed_x = observed_x, sd.prior = sd.prior, boundary.prior = boundary.prior, data = data,
                      region = region, accuracy = accuracy, boundary = boundary)
      # Case for sgp
      instance@initial_location <- initial_location
      instance@X <- global_poly(instance)
      instance@B <- compute_B(instance)
      instance@P <- compute_P(instance)
      if(is.numeric(sd.prior$h)){
        instance@psd.prior <- instance@sd.prior
        instance@sd.prior$param <- prior_conversion_sgp(d = sd.prior$h, prior = sd.prior$param, a = eval(rand_effect$a, envir = envir), m = m)
      } else if(is.numeric(sd.prior$step)){
        instance@sd.prior$h <- sd.prior$step
        instance@psd.prior <- instance@sd.prior
        instance@sd.prior$param <- prior_conversion_sgp(d = sd.prior$step, prior = sd.prior$param, a = eval(rand_effect$a, envir = envir), m = m)
      }
      instances[[length(instances) + 1]] <- instance
    }
  }

  fixed_effects_names <- c()
  # For the intercept
  if (!family_is_coxph & !family_is_cc) {
    Xf0 <- matrix(1, nrow = nrow(data), ncol = 1)
    design_mat_fixed[[length(design_mat_fixed) + 1]] <- Xf0
    fixed_effects_names <- c(fixed_effects_names, "intercept")
  }
  # For fixed effects
  for (fixed_effect in fixed_effects) {
    if(is.numeric(data[[fixed_effect]])){
      Xf <- matrix(data[[fixed_effect]], nrow = nrow(data), ncol = 1)
      design_mat_fixed[[length(design_mat_fixed) + 1]] <- Xf
      fixed_effects_names <- c(fixed_effects_names, fixed_effect)
    }
    else{
      Level <- data[[fixed_effect]]
      Xf <- model.matrix(~ Level)[,-1]
      for (col_index in 1:ncol(Xf)) {
        design_mat_fixed[[length(design_mat_fixed) + 1]] <- Xf[,col_index, drop = FALSE]
      }
      fixed_effects_names <- c(fixed_effects_names, colnames(Xf))
    }
  }

  if (missing(control.fixed)) {
    control.fixed <- list(intercept = list(prec = 0.001, mean = 0))
    for (fixed_effect in fixed_effects_names) {
      control.fixed[[fixed_effect]] <- list(prec = 0.001, mean = 0)
    }
  }
  if(!"intercept" %in% names(control.fixed)){
    control.fixed$intercept <- list(prec = 0.001, mean = 0)
  }
  if(is.null(control.fixed$intercept$prec)){
    control.fixed$intercept$prec <- 0.001
  }
  if(is.null(control.fixed$intercept$mean)){
    control.fixed$intercept$mean <- 0
  }
  
  for (fixed_effect in fixed_effects_names) {
    if(!as.character(fixed_effect) %in% names(control.fixed)){
      control.fixed[[fixed_effect]] <- list(prec = 0.001, mean = 0)
    }
    if(is.null(control.fixed[[fixed_effect]]$prec)){
      control.fixed[[fixed_effect]]$prec <- 0.001
    }
    if(is.null(control.fixed[[fixed_effect]]$mean)){
      control.fixed[[fixed_effect]]$mean <- 0
    }
  }
  
  # For offset effects
  offset_sum <- 0
  for (offset_effect in offset_effects) {
    offset_vec <- data[[offset_effect[[2]]]]
    offset_sum <- offset_sum + offset_vec
  }

  result_by_method <- get_result_by_method(response_var = response_var, data = data, instances = instances, design_mat_fixed = design_mat_fixed, family = family, 
                                           control.family = control.family, control.fixed = control.fixed, 
                                           fixed_effects = fixed_effects, fixed_effects_names = fixed_effects_names, aghq_k = aghq_k, size = size, cens = cens,
                                           weight = weight, strata = strata, offset_sum = offset_sum,
                                           method = method, M = M, option_list = option_list, customized_template = customized_template,
                                           envir = envir, extra_theta_num = extra_theta_num)
  mod <- result_by_method$mod
  w_count <- result_by_method$w_count

  global_samp_indexes <- list()
  coef_samp_indexes <- list()
  fixed_samp_indexes <- list()
  rand_effects_names <- c()
  global_effects_names <- c()
  sum_col_ins <- 0
  for (instance in instances) {
    sum_col_ins <- sum_col_ins + ncol(instance@B)
    rand_effects_names <- c(rand_effects_names, instance@smoothing_var)
    if (class(instance) == "iwp") {
      global_effects_names <- c(global_effects_names, instance@smoothing_var)
    }
    else if (class(instance) == "sgp") {
      global_effects_names <- c(global_effects_names, instance@smoothing_var)
    }
  }

  cur_start <- sum_col_ins + 1
  cur_end <- sum_col_ins
  cur_coef_start <- 1
  cur_coef_end <- 0
  for (instance in instances) {
    if (class(instance) == "iwp") {
      cur_end <- cur_end + ncol(instance@X)
      if (instance@order == 1) {
        global_samp_indexes[[length(global_samp_indexes) + 1]] <- numeric()
      } else if (instance@order > 1) {
        global_samp_indexes[[length(global_samp_indexes) + 1]] <- (cur_start:cur_end)
      }
    }
    else if (class(instance) == "sgp") {
      cur_end <- cur_end + ncol(instance@X)
      global_samp_indexes[[length(global_samp_indexes) + 1]] <- (cur_start:cur_end)
    }
    
    cur_coef_end <- cur_coef_end + ncol(instance@B)
    coef_samp_indexes[[length(coef_samp_indexes) + 1]] <- (cur_coef_start:cur_coef_end)
    cur_start <- cur_end + 1
    cur_coef_start <- cur_coef_end + 1
  }
  names(global_samp_indexes) <- global_effects_names
  names(coef_samp_indexes) <- rand_effects_names

  if((cur_end + 1) <= w_count){
    for (fixed_samp_index in ((cur_end + 1):w_count)) {
      fixed_samp_indexes[[length(fixed_samp_indexes) + 1]] <- fixed_samp_index
    }
  }
  names(fixed_samp_indexes) <- fixed_effects_names
  
  fit_result <- list(
    instances = instances, design_mat_fixed = design_mat_fixed, mod = mod,
    boundary_samp_indexes = global_samp_indexes,
    random_samp_indexes = coef_samp_indexes,
    fixed_samp_indexes = fixed_samp_indexes,
    family = family,
    family_info = build_family_fit_info(
      family = family,
      data = data,
      response_var = response_var,
      size = size,
      cens = cens,
      weight = weight,
      strata = strata
    ),
    control.family = control.family,
    control.fixed = control.fixed,
    near_mono_meta = near_mono_meta,
    exact_iwp_meta = exact_iwp_meta,
    family_sd_known = result_by_method$family_sd_known,
    family_sd_value = result_by_method$family_sd_value
  )
  
  if(any(class(fit_result$mod) == "aghq")){
    fit_result$samps <- aghq::sample_marginal(fit_result$mod, M = M)
  }
  else if(any(class(fit_result$mod) == "nlminb")){
    fit_result$samps <- LaplacesDemon::rmvnp(n = M, mu = fit_result$mod$mean, Omega = as.matrix(fit_result$mod$prec))
  }
  else if(any(class(fit_result$mod) == "stanfit")){
    mc_samps <- rstan::extract(fit_result$mod)
    fit_result$samps <- list(samps = t(mc_samps$W), theta = as.list(data.frame(mc_samps$theta)))
  }

  class(fit_result) <- "FitResult"

  return(fit_result)
}




#' @title Repeated fitting Bayesian Hierarchical Models for a sequence of values of the looping variable.
#'
#' @description
#' Performs repeated model fitting over a sequence of values for a specified variable within a hierarchical model.
#' This function repeatedly fits a model for each value of the looping variable, compiles the log marginal likelihoods,
#' and calculates the posterior probabilities for the variable's values.
#'
#' @param loop_holder A string specifying the name of the variable to loop over. The default value is `LOOP`.
#' @param loop_values A numeric vector containing the values to loop over for the specified variable.
#' @param prior_func A function that takes the specified loop_values and returns the values of the prior for the loop variable.
#' By default, it is a uniform prior which returns a constant value, indicating equal probability for all values.
#' @param parallel Logical, indicating whether or not to run the model fitting in parallel (default is FALSE).
#' @param cores The number of cores to use for parallel execution (default is detected cores - 1).
#' @param ... Additional arguments passed to the model fitting function `model_fit`.
#'
#' @return A data frame containing the values of the looping variable, their corresponding log marginal likelihoods,
#' and posterior probabilities.
#' 
#' @export
model_fit_loop <- function(loop_holder = "LOOP", loop_values, prior_func = function(x){1}, parallel = FALSE, cores = (parallel::detectCores() - 1), ...) {
  update_arg <- function(arg_list, arg_name, arg_value) {
    if (arg_name %in% names(arg_list)) {
      arg_list[[arg_name]] <- arg_value
    } else {
      for (name in names(arg_list)) {
        if (is.function(arg_list[[name]])) {
          formals(arg_list[[name]])[[arg_name]] <- arg_value
        }
        if (is.list(arg_list[[name]])) {
          arg_list[[name]] <- update_arg(arg_list[[name]], arg_name, arg_value)
        }
      }
    }
    return(arg_list)
  }
  
  args_list <- list(...)
  log_ml <- c()
  
  run_model <- function(loop_var) {
    safe_env <- new.env(parent = parent.frame())
    assign(loop_holder, loop_var, envir = safe_env)
    args_list$envir <- safe_env
    args_list <- update_arg(args_list, loop_holder, loop_var)
    model_fit_result <- do.call("model_fit", args_list)
    model_fit_result$mod$normalized_posterior$lognormconst
  }
  
  if (!parallel) {
    for (loop_var in loop_values) {
      log_ml <- c(log_ml, run_model(loop_var))
    }
  } else {
    cl <- parallel::makeCluster(cores)
    parallel::clusterEvalQ(cl, {
      required_packages <- c("Matrix", "aghq", "TMB")  # List the required packages
      lapply(required_packages, function(pkg) {
        library(pkg, character.only = TRUE)
      })
    })
    parallel::clusterExport(cl, list("update_arg", "model_fit", "args_list", "loop_holder", "prior_func"), envir = environment())
    log_ml <- parallel::parLapply(cl, loop_values, run_model)
    parallel::stopCluster(cl)
  }
  
  log_joint <- unlist(log_ml) + log(prior_func(loop_values))
  log_joint <- log_joint - max(log_joint)
  post <- exp(log_joint)
  integral_result <- sfsmisc::integrate.xy(loop_values, post)
  post <- post / integral_result
  data.frame(var = loop_values, post = post, log_ml = unlist(log_ml))
  
}
