get_result_by_method <- function(response_var, data, instances, design_mat_fixed, family, control.family, control.fixed, fixed_effects, aghq_k, size, cens, weight, strata, method, M, customized_template, option_list, envir = parent.frame()) {
  if(is.null(customized_template)){
    cpp = "BayesGP"
  }
  else{
    cpp = customized_template
  }
  # Family types: Gaussian - 0, Poisson - 1, Binomial - 2, Coxph - 3, CaseCrossover - 4
  if (family == "Gaussian") {
    family_type <- 0
  } else if (family == "Poisson") {
    family_type <- 1
  } else if (family == "Binomial") {
    family_type <- 2
  } else if (family == "Coxph" || family == "coxph") {
    family_type <- 3
  } else if(family == "casecrossover" || family == "cc" || family == "CaseCrossover"){
    family_type <- 4
  }else if(family == "Customized"){
    if(is.null(customized_template)){
      stop("In order to use the customized family: please input the name of the complied cpp template as `customized_template`.")
    }
    family_type <- -1
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
  # Need a theta for the Gaussian variance, so
  # theta_count starts at 1 if Gaussian
  theta_count <- 0 + (family_type == 0)
  
  for (instance in instances) {
    # For each random effects
    if (class(instance) == "IWP") {
      X[[length(X) + 1]] <- dgTMatrix_wrapper(instance@X)
      betaprec[[length(betaprec) + 1]] <- instance@boundary.prior$prec
      betamean[[length(betamean) + 1]] <- instance@boundary.prior$mean
      w_count <- w_count + ncol(instance@X)
    }
    else if(class(instance) == "sGP"){
      X[[length(X) + 1]] <- dgTMatrix_wrapper(instance@X)
      betaprec[[length(betaprec) + 1]] <- instance@boundary.prior$prec
      betamean[[length(betamean) + 1]] <- instance@boundary.prior$mean
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
  if (family_type == 0) {
    
    if (is.null(control.family$sd_prior)) {
      control.family$sd_prior <- eval(control.family$prior, envir = envir)
      if(is.null(control.family$sd_prior)){
        control.family$sd_prior <- list(prior = "exp", param = list(u = 1, alpha = 0.5))
      }
    }
    if (length(control.family$sd_prior) == 1){
      if(is.numeric(control.family$sd_prior)){
        control.family$sd_prior <- list(prior = "exp", param = list(u = as.numeric(control.family$sd_prior), alpha = 0.5))
      }
    } 
    
    if(!"prior" %in% names(control.family$sd_prior)){
      control.family$sd_prior$prior <- "exp"
    }
    if(!"param" %in% names(control.family$sd_prior)){
      stop("If sd.prior is provided as a list, it must contains a list called param.")
    }else{
      if(length(control.family$sd_prior$param) == 1){
        control.family$sd_prior$param <- list(u = control.family$sd_prior$param[[1]], alpha = 0.5)
      }
      else{
        control.family$sd_prior$param <- list(u = control.family$sd_prior$param$u, alpha = control.family$sd_prior$param$alpha)
        if(is.null(control.family$sd_prior$param$alpha)){
          warnings("The value of alpha is not provided in control.family$sd_prior$param: automatically filled with 0.5.")
          control.family$sd_prior$param$alpha <- 0.5
        }
        if(is.null(control.family$sd_prior$param$u)){
          stop("Error: The value of u is not provided in control.family$sd_prior$param.")
        }
      }
    }
    
    if (control.family$sd_prior$prior != "exp" & control.family$sd_prior$prior != "Exp" & control.family$sd_prior$prior != "exponential" & control.family$sd_prior$prior != "Exponential") {
      stop("Error: For each random effect, control.family$sd_prior currently only supports 'exp' (exponential) as prior.")
    }
    if(control.family$sd_prior$param$alpha > 1 | control.family$sd_prior$param$alpha < 0){
      stop("Error: The value of control.family$sd_prior$param$alpha is not specified as a probability.")
    }
    
    
    u[[length(u) + 1]] <- control.family$sd_prior$param$u
    alpha[[length(alpha) + 1]] <- control.family$sd_prior$param$alpha
  }
  if(family_type == 3 || family_type == 4){
    if(length(design_mat_fixed) >= 1){
      for (i in 1:length(design_mat_fixed)) {
        # For each fixed effects
        beta_fixed_prec[[i]] <- control.fixed[[fixed_effects[[i]]]]$prec
        beta_fixed_mean[[i]] <- control.fixed[[fixed_effects[[i]]]]$mean
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
        beta_fixed_prec[[i]] <- control.fixed[[fixed_effects[[i - 1]]]]$prec
        beta_fixed_mean[[i]] <- control.fixed[[fixed_effects[[i - 1]]]]$mean
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
    
    # Family type
    family_type = family_type
  )
  
  # If Family == "Binomial", check whether size is defined in user's input
  if (family_type == 2) {
    if(is.null(data[[size]])){
      tmbdat$size <- numeric(length = length(tmbdat$y)) + 1 # A vector of 1s being default
    }
    else{
      tmbdat$size <- data[[size]]
    }
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
  
  # If Family == "cc", check whether strata and weight is defined in user's input
  if (family_type == 4) {
    case <- data[[response_var]]
    if(is.null(weight)){
      weight <- data[[response_var]]
    }
    else{
      weight <- data[[weight]]
    }
    if(is.null(data[[strata]])){
      stop("The specified names for strata are not correct.")
    }
    else{
      strata <- data[[strata]]
    }
    case_day <- which(case > 0)
    count <- weight[case_day] 
    
    # Initialize a list to store the control_day for each strata
    unique_strata <- unique(strata)
    control_day_list <- list()
    
    max_N <- max(sapply(unique_strata, function(s) sum(strata == s & case == 0)))
    
    # Loop through each unique strata
    for (s in unique_strata) {
      # Get the indices for case days and control days in the current strata
      case_day_strata <- which(strata == s & case > 0)
      control_day_strata <- which(strata == s & case == 0)
      
      # Initialize the control_day matrix for this strata
      M <- length(case_day_strata)
      N <- max_N
      control_day_matrix <- matrix(0, nrow = M, ncol = N + 1)
      
      # Populate the control_day matrix
      for (i in seq_along(case_day_strata)) {
        control_day_matrix[i, 1] <- case_day_strata[i]
        control_day_matrix[i, 2:(length(control_day_strata) + 1)] <- control_day_strata
      }
      
      # Add this strata's control_days matrix to the list
      control_day_list[[as.character(s)]] <- control_day_matrix
    }
    control_days <- do.call(rbind, control_day_list)
    
    tmbdat$control_days <- control_days
    tmbdat$case_day <- case_day
    tmbdat$count <- count
  }
  
  tmbparams <- list(
    W = c(rep(0, w_count)), # recall W is everything in the model (RE or FE)
    theta = c(rep(0, theta_count))
  )
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
    mod <- tmbstan::tmbstan(
      ff,
      chains = default_option_list$chains,
      cores = default_option_list$cores,
      iter = default_option_list$warmup + M,
      warmup = default_option_list$warmup,
      seed = default_option_list$seed
    )
  }
  return(list(mod = mod, w_count = w_count))
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
#' @param family The family of response used in the model. By default, the family is set to be "Gaussian".
#' @param control.family Parameters used to specify the priors for the family parameters, such as the standard deviation parameter of Gaussian family.
#' @param control.fixed Parameters used to specify the priors for the fixed effects.
#' @param cens The name of the right-censoring indicator, should be one of the variables in `data`. The default value is "NULL".
#' @param M The number of posterior samples to be taken, by default is 3000.
#' @param customized_template The name of the customized cpp template that the user wants to use instead. By default this is NULL, and the cpp template `BayesGP` will be used.
#' @param Customized_RE The list that contains the compute_B and compute_P functions for the customized random effect. By default, this is NULL and there is not customized random effect in the model.
#' @param option_list A list that controls the details of the inference algorithm, by default is an empty list.
#' @param envir The environment in which the formula and other expressions are to be evaluated. 
#'   Defaults to `parent.frame()`, which refers to the environment from which the function was called.
#'   This allows the function to access variables that are defined in the calling function's scope.
#' @return A list that contains following items: the S4 objects for the random effects (instances), concatenated design matrix for
#' the fixed effects (design_mat_fixed), fitted aghq (mod) and indexes to partition the posterior samples
#' (boundary_samp_indexes, random_samp_indexes and fixed_samp_indexes).
#'
#' @export
model_fit <- function(formula, data, method = "aghq", family = "Gaussian", control.family, control.fixed, aghq_k = 4, size = NULL, cens = NULL, weight = NULL, strata = NULL, M = 3000, customized_template = NULL, Customized_RE = NULL, option_list = list(), envir = parent.frame()) {
  # parse the input formula
  parse_result <- parse_formula(formula)
  response_var <- parse_result$response
  rand_effects <- parse_result$rand_effects
  fixed_effects <- parse_result$fixed_effects

  instances <- list()
  design_mat_fixed <- list()
  
  family_is_coxph <- FALSE
  if(family == "Coxph" || family == "coxph"){
    family_is_coxph <- TRUE
    data <- data[order(data[[response_var]]), ]
  }
  
  family_is_cc <- FALSE
  if(family == "casecrossover" || family == "cc" || family == "CaseCrossover"){
    family_is_cc <- TRUE
  }

  # For random effects
  for (rand_effect in rand_effects) {
    smoothing_var <- rand_effect$smoothing_var
    if(is.null(smoothing_var)){
      smoothing_var <- rand_effect$x
      if(is.null(smoothing_var)){
        stop("Error: The covariate name must be specified as smoothing_var or x.")
      }
    }
    model_class <- rand_effect$model
    
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
    
    if (sd.prior$prior != "exp" & sd.prior$prior != "Exp" & sd.prior$prior != "exponential" & sd.prior$prior != "Exponential") {
      stop("Error: For each random effect, sd.prior currently only supports 'exp' (exponential) as prior.")
    }
    if(sd.prior$param$alpha > 1 | sd.prior$param$alpha < 0){
      stop("Error: The value of sd.prior$param$alpha is not specified as a probability.")
    }
    
    if (model_class == "IWP") {
      order <- eval(rand_effect$order, envir = envir)
      knots <- eval(rand_effect$knots, envir = envir)
      k <- eval(rand_effect$k, envir = envir)
      initial_location <- eval(rand_effect$initial_location, envir = envir)
      if (!(is.null(k)) && k < 3) {
        stop("Error: parameter <k> in the random effect part should be >= 3.")
      }
      if (is.null(order) | order < 1) {
        stop("Error: Parameter <order> in the random effect part should be >= 1.")
      }
      boundary.prior <- eval(rand_effect$boundary.prior, envir = envir)
      # If the user does not specify initial_location, compute initial_location with
      # the min of data[[smoothing_var]]
      if (is.null(initial_location)) {
        initial_location <- min(data[[smoothing_var]])
      }
      # If the user does not specify knots, compute knots with
      # the parameter k
      initialized_smoothing_var <- data[[smoothing_var]] - initial_location
      if (is.null(knots)) {
        default_k <- 5
        if (is.null(k)) {
          knots <- unique(sort(seq(from = min(initialized_smoothing_var), to = max(initialized_smoothing_var), length.out = default_k))) # should be length.out
        } else {
          knots <- unique(sort(seq(from = min(initialized_smoothing_var), to = max(initialized_smoothing_var), length.out = k)))
        }
      }
      observed_x <- sort(initialized_smoothing_var) # initialized_smoothing_var: initialized observed covariate values
      if (is.null(boundary.prior)) {
        boundary.prior <- list(prec = 0.01, mean = 0)
      }
      if(is.null(boundary.prior$prec)){
        boundary.prior$prec <- 0.01
      }
      if(is.null(boundary.prior$mean)){
        boundary.prior$mean <- 0
      }
      instance <- new(model_class,
        response_var = response_var,
        smoothing_var = smoothing_var, order = order,
        knots = knots, observed_x = observed_x, sd.prior = sd.prior, boundary.prior = boundary.prior, data = data
      )
      # Case for IWP
      instance@initial_location <- initial_location
      instance@X <- global_poly(instance)[, -1, drop = FALSE]
      instance@B <- local_poly(instance)
      instance@P <- compute_weights_precision(instance)
      instances[[length(instances) + 1]] <- instance
      
      if(is.numeric(sd.prior$h)){
        sd.prior$param <- prior_conversion_IWP(d = sd.prior$h, prior = sd.prior$param, p = eval(rand_effect$order, envir = envir))
      } else if(is.numeric(sd.prior$step)){
        sd.prior$param <- prior_conversion_IWP(d = sd.prior$step, prior = sd.prior$param, p = eval(rand_effect$order, envir = envir))
      }
      
    } 
    else if (model_class == "IID") {
      instance <- new(model_class,
        response_var = response_var,
        smoothing_var = smoothing_var, sd.prior = sd.prior, data = data
      )
      # Case for IID
      instance@B <- compute_B(instance)
      instance@P <- compute_P(instance)
      instances[[length(instances) + 1]] <- instance
    }
    else if (model_class == "Customized") {
      instance <- new(model_class,
                      response_var = response_var,
                      smoothing_var = smoothing_var, sd.prior = sd.prior, data = data,
                      compute_B = Customized_RE$compute_B, compute_P = Customized_RE$compute_P
      )
      # Case for Customized
      instance@B <- compute_B(instance)
      instance@P <- compute_P(instance)
      instances[[length(instances) + 1]] <- instance
    }
    else if (model_class == "sGP"){
      a <- eval(rand_effect$a, envir = envir)
      k <- eval(rand_effect$k, envir = envir)
      initial_location <- eval(rand_effect$initial_location, envir = envir)
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
      # If the user does not specify initial_location, compute initial_location with
      # the min of data[[smoothing_var]]
      if (is.null(initial_location)) {
        initial_location <- min(data[[smoothing_var]])
      }
      initialized_smoothing_var <- data[[smoothing_var]] - initial_location
      observed_x <- sort(initialized_smoothing_var) # initialized_smoothing_var: initialized observed covariate values
      if("region" %in% names(rand_effect)){
        region <- eval(rand_effect$region, envir = envir)
      }
      else{
        region <- range(observed_x)
      }
      if("accuracy" %in% names(rand_effect)){
        accuracy <- eval(rand_effect$accuracy, envir = envir)
      }
      else{
        accuracy <- 0.01
      }
      if (is.null(boundary.prior)) {
        boundary.prior <- list(prec = 0.01, mean = 0)
      }
      if(is.null(boundary.prior$prec)){
        boundary.prior$prec <- 0.01
      }
      if(is.null(boundary.prior$mean)){
        boundary.prior$mean <- 0
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
      # Case for sGP
      instance@initial_location <- initial_location
      instance@X <- global_poly(instance)
      instance@B <- compute_B(instance)
      instance@P <- compute_P(instance)
      instances[[length(instances) + 1]] <- instance
      
      if(is.numeric(sd.prior$h)){
        sd.prior$param <- prior_conversion_sGP(d = sd.prior$h, prior = sd.prior$param, a = eval(rand_effect$a, envir = envir), m = m)
      } else if(is.numeric(sd.prior$step)){
        sd.prior$param <- prior_conversion_sGP(d = sd.prior$step, prior = sd.prior$param, a = eval(rand_effect$a, envir = envir), m = m)
      }
    }
  }

  fixed_effects_names <- c()
  # For the intercept
  if (!family_is_coxph & !family_is_cc) {
    Xf0 <- matrix(1, nrow = nrow(data), ncol = 1)
    design_mat_fixed[[length(design_mat_fixed) + 1]] <- Xf0
    fixed_effects_names <- c(fixed_effects_names, "Intercept")
  }
  # For fixed effects
  for (fixed_effect in fixed_effects) {
    Xf <- matrix(data[[fixed_effect]], nrow = nrow(data), ncol = 1)
    design_mat_fixed[[length(design_mat_fixed) + 1]] <- Xf
    fixed_effects_names <- c(fixed_effects_names, fixed_effect)
  }

  if (missing(control.family)) {
    control.family <- list(sd_prior = list(prior = "exp", param = list(u = 1, alpha = 0.5)))
  }

  if (missing(control.fixed)) {
    control.fixed <- list(intercept = list(prec = 0.01, mean = 0))
    for (fixed_effect in fixed_effects) {
      control.fixed[[fixed_effect]] <- list(prec = 0.01, mean = 0)
    }
  }
  if(!"intercept" %in% names(control.fixed)){
    control.fixed$intercept <- list(prec = 0.01, mean = 0)
  }
  if(is.null(control.fixed$intercept$prec)){
    control.fixed$intercept$prec <- 0.01
  }
  if(is.null(control.fixed$intercept$mean)){
    control.fixed$intercept$mean <- 0
  }
  
  for (fixed_effect in fixed_effects) {
    if(!as.character(fixed_effect) %in% names(control.fixed)){
      control.fixed[[fixed_effect]] <- list(prec = 0.01, mean = 0)
    }
    if(is.null(control.fixed[[fixed_effect]]$prec)){
      control.fixed[[fixed_effect]]$prec <- 0.01
    }
    if(is.null(control.fixed[[fixed_effect]]$mean)){
      control.fixed[[fixed_effect]]$mean <- 0
    }
  }

  result_by_method <- get_result_by_method(response_var = response_var, data = data, instances = instances, design_mat_fixed = design_mat_fixed, family = family, 
                                           control.family = control.family, control.fixed = control.fixed, 
                                           fixed_effects = fixed_effects, aghq_k = aghq_k, size = size, cens = cens,
                                           weight = weight, strata = strata,
                                           method = method, M = M, option_list = option_list, customized_template = customized_template,
                                           envir = envir)
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
    if (class(instance) == "IWP") {
      global_effects_names <- c(global_effects_names, instance@smoothing_var)
    }
    else if (class(instance) == "sGP") {
      global_effects_names <- c(global_effects_names, instance@smoothing_var)
    }
  }

  cur_start <- sum_col_ins + 1
  cur_end <- sum_col_ins
  cur_coef_start <- 1
  cur_coef_end <- 0
  for (instance in instances) {
    if (class(instance) == "IWP") {
      cur_end <- cur_end + ncol(instance@X)
      if (instance@order == 1) {
        global_samp_indexes[[length(global_samp_indexes) + 1]] <- numeric()
      } else if (instance@order > 1) {
        global_samp_indexes[[length(global_samp_indexes) + 1]] <- (cur_start:cur_end)
      }
    }
    else if (class(instance) == "sGP") {
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
    family = family
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




#' @title Bayesian Model Averaging for Hierarchical Models
#'
#' @description
#' Performs Bayesian Model Averaging over a sequence of values for a specified variable within a hierarchical model.
#' This function repeatedly fits a model for each value of the looping variable, compiles the log marginal likelihoods,
#' and calculates the posterior probabilities for the variable's values.
#'
#' @param loop_var_name A string specifying the name of the variable to loop over.
#' @param loop_var_values A numeric vector containing the values to loop over for the specified variable.
#' @param prior_func A function that takes the specified loop_var_values and returns the values of the prior for the loop variable.
#' By default, it is a uniform prior which returns a constant value, indicating equal probability for all values.
#' @param parallel Logical, indicating whether or not to run the model fitting in parallel (default is FALSE).
#' @param cores The number of cores to use for parallel execution (default is detected cores - 1).
#' @param ... Additional arguments passed to the model fitting function `model_fit`.
#'
#' @return A data frame containing the values of the looping variable, their corresponding log marginal likelihoods,
#' and posterior probabilities.
#' 
#' @export
model_fit_BMA <- function(loop_var_name, loop_var_values, prior_func = function(x){1}, parallel = FALSE, cores = (parallel::detectCores() - 1), ...) {
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
    assign(loop_var_name, loop_var, envir = safe_env)
    args_list$envir <- safe_env
    args_list <- update_arg(args_list, loop_var_name, loop_var)
    model_fit_result <- do.call("model_fit", args_list)
    model_fit_result$mod$normalized_posterior$lognormconst
  }
  
  if (!parallel) {
    for (loop_var in loop_var_values) {
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
    parallel::clusterExport(cl, list("update_arg", "model_fit", "args_list", "loop_var_name", "prior_func"), envir = environment())
    log_ml <- parallel::parLapply(cl, loop_var_values, run_model)
    parallel::stopCluster(cl)
  }
  
  log_joint <- unlist(log_ml) + log(prior_func(loop_var_values))
  log_joint <- log_joint - max(log_joint)
  post <- exp(log_joint)
  integral_result <- sfsmisc::integrate.xy(loop_var_values, post)
  post <- post / integral_result
  data.frame(var = loop_var_values, post = post, log_ml = unlist(log_ml))
  
}


