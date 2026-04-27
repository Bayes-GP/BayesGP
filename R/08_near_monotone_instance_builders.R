build_nearmono_instance <- function(rand_effect, response_var, data, envir = parent.frame()){
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

  model_name <- tolower(extract_term_argument(rand_effect, "model", envir = envir))
  computation_method <- resolve_computation_method(
    computation = extract_term_argument(rand_effect, "computation", envir = envir),
    method = extract_term_argument(rand_effect, "method", envir = envir),
    default = "fem",
    label = "Near-monotone smooth terms"
  )
  curvature <- resolve_curvature_parameter(
    a = extract_term_argument(rand_effect, "a", envir = envir),
    alpha = extract_term_argument(rand_effect, "alpha", envir = envir)
  )
  c_argument <- extract_term_argument(
    rand_effect,
    "c",
    envir = envir,
    default = NULL
  )
  if(!is.null(c_argument)){
    if(length(c_argument) != 1){
      stop("The shift parameter `c` must be scalar.")
    }
    if(!is.numeric(c_argument) || !is.finite(c_argument)){
      stop("The shift parameter `c` must be a finite scalar.")
    }
    c_argument <- as.numeric(c_argument)
  }

  observed_x <- data[[smoothing_var]]
  grid <- extract_term_argument(rand_effect, "grid", envir = envir, default = NULL)
  region <- extract_term_argument(rand_effect, "region", envir = envir, default = NULL)
  range_arg <- extract_term_argument(rand_effect, "range", envir = envir, default = NULL)

  if(!is.null(region) && !is.null(range_arg) && !isTRUE(all.equal(as.numeric(region), as.numeric(range_arg)))){
    stop("Supply only one of `region` and `range`, or keep them identical.")
  }
  if(is.null(region)){
    region <- range_arg
  }
  if(!is.null(region)){
    if(computation_method != "fem"){
      stop("Near-monotone `region` / `range` is supported only for `computation = \"fem\"`; use `grid` for state-space support.")
    }
    if(!is.numeric(region) || length(region) != 2 || any(!is.finite(region)) || diff(range(region)) <= 0){
      stop("Near-monotone `region` / `range` must be a finite numeric vector of length 2 with distinct endpoints.")
    }
    region <- range(as.numeric(region))
  }

  domain_x <- if(is.null(grid)) observed_x else c(observed_x, grid)
  if(!is.null(region)){
    domain_span <- range(domain_x)
    tolerance <- sqrt(.Machine$double.eps)
    if(domain_span[1] < region[1] - tolerance || domain_span[2] > region[2] + tolerance){
      stop("Near-monotone `region` / `range` must span the observed values and any supplied `grid` values.")
    }
    domain_x <- c(domain_x, region)
  }
  c_original <- resolve_nearmono_shift(
    c = c_argument,
    domain_x = domain_x
  )
  initial_location <- resolve_reference_location(
    x = domain_x,
    initial_location = extract_term_argument(rand_effect, "initial_location", envir = envir, default = "left")
  )
  c_shift <- initial_location + c_original
  shifted_observed_x <- observed_x - initial_location
  shifted_domain_x <- domain_x - initial_location
  reference_is_left <- isTRUE(all.equal(initial_location, min(domain_x)))
  reference_is_right <- isTRUE(all.equal(initial_location, max(domain_x)))
  two_sided_reference <- !reference_is_left && !reference_is_right

  normalized_boundary <- extract_term_argument(
    rand_effect,
    "normalized_boundary",
    envir = envir,
    default = TRUE
  )
  if(!is.null(extract_term_argument(rand_effect, "allow_zero_c", envir = envir, default = NULL))){
    stop(
      "`allow_zero_c` is no longer supported. Omit `c` to use the ",
      "automatic domain-aware shift, or provide a finite `c` such that ",
      "`x + c > 0` over the near-monotone domain."
    )
  }

  validate_nearmono_shift_domain(
    shifted_domain_x = shifted_domain_x,
    c_shift = c_shift
  )

  sd.prior <- normalize_sd_prior(extract_term_argument(rand_effect, "sd.prior", envir = envir))
  converted_sd_prior <- convert_nearmono_sd_prior(
    sd.prior = sd.prior,
    model_name = model_name,
    curvature = curvature,
    c_original = c_original,
    c_shift = c_shift,
    initial_location = initial_location
  )
  sd.prior <- converted_sd_prior$sd.prior
  psd.prior <- converted_sd_prior$psd.prior
  boundary.prior <- normalize_boundary_prior(
    extract_term_argument(rand_effect, "boundary.prior", envir = envir),
    dimension = 1
  )

  X <- boundary_basis_matrix(
    x = observed_x,
    a = curvature,
    c = c_shift,
    initial_location = initial_location,
    normalized = normalized_boundary
  )

  transformed_support_grid <- numeric()
  transformed_knots <- numeric()
  support_grid <- numeric()
  stored_k <- 0
  fem_accuracy <- extract_term_argument(rand_effect, "accuracy", envir = envir, default = 0.01)
  fem_boundary <- isTRUE(extract_term_argument(rand_effect, "boundary", envir = envir, default = TRUE))

  if(computation_method == "state-space"){
    if(model_name == "mgp"){
      if(curvature %in% c(-2, -0.5)){
        stop("The exact state-space implementation is not available for `a = -2` or `a = -1/2`.")
      }
      positive_support <- sort(unique(shifted_domain_x[shifted_domain_x > 0]))
      negative_support <- sort(unique(abs(shifted_domain_x[shifted_domain_x < 0])))
      support_grid <- sort(unique(abs(shifted_domain_x[shifted_domain_x != 0])))
      P <- build_mgp_state_space_precision(
        positive_support = positive_support,
        negative_support = negative_support,
        curvature = curvature,
        c_shift = c_shift,
        include_negative = !reference_is_left
      )
      B <- build_two_sided_exact_observation_matrix(
        observed_x = shifted_observed_x,
        positive_support = positive_support,
        negative_support = if(reference_is_left) numeric() else negative_support,
        state_size = 2
      )
      attr(B, "raw_positive_support") <- positive_support
      attr(B, "raw_negative_support") <- if(reference_is_left) numeric() else negative_support
      transformed_support_grid <- support_grid
    } else if(model_name == "tiwp2"){
      transform <- make_boxcox_transform(
        a = curvature,
        c = c_shift,
        normalized = TRUE,
        ref_location = 0
      )
      transformed_observed_x <- transform(shifted_observed_x)
      transformed_domain_x <- transform(shifted_domain_x)
      raw_positive_support <- sort(unique(shifted_domain_x[shifted_domain_x > 0]))
      raw_negative_support <- sort(unique(abs(shifted_domain_x[shifted_domain_x < 0])))
      transformed_positive_support <- sort(unique(transformed_domain_x[transformed_domain_x > 0]))
      transformed_negative_support <- sort(unique(abs(transformed_domain_x[transformed_domain_x < 0])))
      transformed_support_grid <- sort(unique(abs(transformed_domain_x[transformed_domain_x != 0])))
      support_grid <- sort(unique(abs(shifted_domain_x[shifted_domain_x != 0])))
      P <- build_block_precision(
        positive_grid = transformed_positive_support,
        negative_grid = transformed_negative_support,
        precision_builder = function(grid) build_exact_precision_with_diagnostics(
          support_grid = grid,
          model_name = "tiwp2",
          precision_builder = function(support) IWP_joint_prec(t_vec = support, p = 2),
          scale_label = "transformed support grid"
        )
      )
      B <- build_two_sided_exact_observation_matrix(
        observed_x = transformed_observed_x,
        positive_support = transformed_positive_support,
        negative_support = transformed_negative_support,
        state_size = 2
      )
      attr(B, "transformed_positive_support") <- attr(B, "positive_support")
      attr(B, "transformed_negative_support") <- attr(B, "negative_support")
      attr(B, "raw_positive_support") <- raw_positive_support
      attr(B, "raw_negative_support") <- raw_negative_support
    } else {
      stop("Unknown near-monotone model `", model_name, "`.")
    }
  } else {
    if(model_name == "mgp"){
      k_input <- extract_term_argument(rand_effect, "k", envir = envir, default = NULL)
      knots_input <- extract_term_argument(rand_effect, "knots", envir = envir, default = NULL)
      if(!is.null(knots_input) && length(knots_input) > 1){
        stop("Explicit knot sequences are not currently supported for `mgp` FEM. Please supply `k` instead.")
      }
      k <- if(!is.null(k_input)) k_input else if(!is.null(knots_input)) knots_input else 30
      stored_k <- k
      support_grid <- sort(unique(shifted_domain_x))
      fem_mgp <- build_two_sided_mgp_fem(
        observed_x = shifted_observed_x,
        domain_x = shifted_domain_x,
        k = k,
        a = curvature,
        c = c_shift,
        accuracy = fem_accuracy,
        boundary = fem_boundary
      )
      B <- fem_mgp$B
      P <- fem_mgp$P
      transformed_support_grid <- support_grid
    } else if(model_name == "tiwp2"){
      k_input <- extract_term_argument(rand_effect, "k", envir = envir, default = NULL)
      knots_input <- extract_term_argument(rand_effect, "knots", envir = envir, default = NULL)
      transform <- make_boxcox_transform(
        a = curvature,
        c = c_shift,
        normalized = TRUE,
        ref_location = 0
      )
      transformed_observed_x <- transform(shifted_observed_x)
      transformed_domain_x <- transform(shifted_domain_x)
      transformed_support_grid <- sort(unique(transformed_domain_x))
      k <- if(!is.null(k_input)) k_input else if(!is.null(knots_input) && length(knots_input) == 1) knots_input else 30
      stored_k <- k
      if(!is.null(knots_input) && length(knots_input) > 1){
        stop("For `tiwp2` FEM, explicit knot vectors are not currently supported. Supply `k` so knots can be generated equally on the transformed scale.")
      }
      fem_tiwp <- build_two_sided_tiwp2_fem(
        observed_x = transformed_observed_x,
        domain_x = transformed_domain_x,
        k = k,
        response_var = response_var,
        smoothing_var = smoothing_var,
        sd.prior = sd.prior,
        psd.prior = psd.prior,
        boundary.prior = boundary.prior
      )
      transformed_knots <- fem_tiwp$knots
      B <- fem_tiwp$B
      P <- fem_tiwp$P
      support_grid <- sort(unique(shifted_domain_x))
    } else {
      stop("Unknown near-monotone model `", model_name, "`.")
    }
  }

  if(computation_method == "state-space"){
    stored_k <- ncol(B)
  }

  instance <- methods::new(
    "iwp",
    response_var = as.name(response_var),
    smoothing_var = as.name(smoothing_var),
    order = 2,
    knots = numeric(),
    k = stored_k,
    observed_x = sort(shifted_observed_x),
    sd.prior = sd.prior,
    psd.prior = psd.prior,
    boundary.prior = boundary.prior,
    data = data,
    X = as.matrix(X),
    B = as.matrix(B),
    P = as.matrix(P),
    initial_location = initial_location,
    region = range(domain_x)
  )

  list(
    instance = instance,
    meta = list(
      component = smoothing_var,
      model_name = model_name,
      computation_method = computation_method,
      curvature = curvature,
      c = c_original,
      c_is_default = is.null(c_argument),
      shift = c_shift,
      internal_shift = c_shift,
      support_grid = support_grid,
      transformed_support_grid = transformed_support_grid,
      transformed_knots = transformed_knots,
      transformed_knots_positive = attr(B, "transformed_knots_positive"),
      transformed_knots_negative = attr(B, "transformed_knots_negative"),
      normalized_boundary = normalized_boundary,
      fem_accuracy = fem_accuracy,
      fem_boundary = fem_boundary,
      reference_location = initial_location,
      reference_is_left = reference_is_left,
      reference_is_right = reference_is_right,
      two_sided_reference = two_sided_reference,
      positive_support = attr(B, "positive_support"),
      negative_support = attr(B, "negative_support"),
      transformed_positive_support = attr(B, "transformed_positive_support"),
      transformed_negative_support = attr(B, "transformed_negative_support"),
      raw_positive_support = attr(B, "raw_positive_support"),
      raw_negative_support = attr(B, "raw_negative_support")
    )
  )
}

resolve_nearmono_shift <- function(c, domain_x){
  if(!is.null(c)){
    return(as.numeric(c))
  }

  if(!is.numeric(domain_x) || length(domain_x) == 0 || any(!is.finite(domain_x))){
    stop("Near-monotone smooth terms require a finite numeric domain to choose the default `c`.")
  }

  lower_endpoint <- min(domain_x)
  if(lower_endpoint > 0){
    return(0)
  }

  -lower_endpoint + nearmono_shift_margin(domain_x)
}

nearmono_shift_margin <- function(domain_x){
  domain_range <- range(domain_x)
  scale <- max(1, diff(domain_range), max(abs(domain_x)))
  1e-6 * scale
}

validate_nearmono_shift_domain <- function(shifted_domain_x, c_shift){
  transformed_arg <- shifted_domain_x + c_shift

  if(any(transformed_arg <= 0)){
    stop(
      "The chosen `c` and input domain must satisfy `x + c > 0` ",
      "on the original scale. Omit `c` to use the automatic ",
      "domain-aware shift."
    )
  }

  invisible(TRUE)
}

build_block_precision <- function(positive_grid, negative_grid, precision_builder){
  blocks <- list()

  if(length(positive_grid) > 0){
    blocks[[length(blocks) + 1]] <- as.matrix(precision_builder(positive_grid))
  }
  if(length(negative_grid) > 0){
    blocks[[length(blocks) + 1]] <- as.matrix(precision_builder(negative_grid))
  }

  if(length(blocks) == 0){
    return(matrix(0, nrow = 0, ncol = 0))
  }
  if(length(blocks) == 1){
    return(blocks[[1]])
  }

  as.matrix(Matrix::bdiag(blocks))
}

build_mgp_state_space_precision <- function(positive_support, negative_support,
                                            curvature, c_shift,
                                            include_negative = TRUE){
  blocks <- list()

  if(length(positive_support) > 0){
    blocks[[length(blocks) + 1]] <- as.matrix(build_exact_precision_with_diagnostics(
      support_grid = positive_support,
      model_name = "mgp",
      precision_builder = function(grid) mGP_joint_prec(t_vec = grid, a = curvature, c = c_shift)
    ))
  }

  if(include_negative && length(negative_support) > 0){
    blocks[[length(blocks) + 1]] <- as.matrix(build_exact_precision_with_diagnostics(
      support_grid = negative_support,
      model_name = "mgp",
      precision_builder = function(grid) adj_mGP_joint_prec(t_vec = grid, a = curvature, c = c_shift),
      scale_label = "reverse support grid"
    ))
  }

  if(length(blocks) == 0){
    return(matrix(0, nrow = 0, ncol = 0))
  }
  if(length(blocks) == 1){
    return(blocks[[1]])
  }

  as.matrix(Matrix::bdiag(blocks))
}

build_two_sided_exact_observation_matrix <- function(observed_x, positive_support,
                                                     negative_support,
                                                     state_size = 2){
  positive_support <- sort(unique(positive_support))
  negative_support <- sort(unique(negative_support))
  n_positive <- length(positive_support)
  n_negative <- length(negative_support)

  observation_matrix <- matrix(
    0,
    nrow = length(observed_x),
    ncol = state_size * (n_positive + n_negative)
  )

  positive_rows <- which(observed_x > 0)
  if(length(positive_rows) > 0 && n_positive > 0){
    positive_match <- match(observed_x[positive_rows], positive_support)
    matched_rows <- positive_rows[!is.na(positive_match)]
    matched_index <- positive_match[!is.na(positive_match)]
    observation_matrix[cbind(matched_rows, state_size * matched_index - (state_size - 1))] <- 1
  }

  negative_rows <- which(observed_x < 0)
  if(length(negative_rows) > 0 && n_negative > 0){
    negative_match <- match(abs(observed_x[negative_rows]), negative_support)
    matched_rows <- negative_rows[!is.na(negative_match)]
    matched_index <- negative_match[!is.na(negative_match)]
    offset <- state_size * n_positive
    observation_matrix[cbind(matched_rows, offset + state_size * matched_index - (state_size - 1))] <- 1
  }

  attr(observation_matrix, "positive_support") <- positive_support
  attr(observation_matrix, "negative_support") <- negative_support
  observation_matrix
}

build_side_design <- function(x, side = c("positive", "negative"), k, region,
                              boundary = TRUE, basis_builder){
  side <- match.arg(side)
  side_x <- if(side == "positive") pmax(x, 0) else pmax(-x, 0)
  active_rows <- if(side == "positive") which(x > 0) else which(x < 0)

  if(length(region) == 0 || max(region) <= 0){
    return(matrix(0, nrow = length(x), ncol = 0))
  }

  if(length(active_rows) > 0){
    design_active <- basis_builder(side_x[active_rows], k = k, region = c(0, max(region)), boundary = boundary)
  } else {
    design_active <- basis_builder(0, k = k, region = c(0, max(region)), boundary = boundary)
  }
  design <- matrix(0, nrow = length(x), ncol = ncol(design_active))
  if(length(active_rows) > 0){
    design[active_rows, ] <- design_active
  }
  design
}

build_two_sided_mgp_fem <- function(observed_x, domain_x, k, a, c,
                                    accuracy = 0.01, boundary = TRUE){
  positive_region <- range(c(0, domain_x[domain_x > 0]))
  negative_region <- range(c(0, -domain_x[domain_x < 0]))

  B_pos <- build_side_design(
    x = observed_x,
    side = "positive",
    k = k,
    region = positive_region,
    boundary = boundary,
    basis_builder = function(x, k, region, boundary){
      as.matrix(Compute_Design(x = x, k = k, region = region, boundary = boundary))
    }
  )
  B_neg <- build_side_design(
    x = observed_x,
    side = "negative",
    k = k,
    region = negative_region,
    boundary = boundary,
    basis_builder = function(x, k, region, boundary){
      as.matrix(Compute_Design(x = x, k = k, region = region, boundary = boundary))
    }
  )

  precision_blocks <- list()
  if(max(positive_region) > 0){
    precision_blocks[[length(precision_blocks) + 1]] <- as.matrix(
      Compute_Prec(a = a, c = c, k = k, region = positive_region,
                   accuracy = accuracy, boundary = boundary)
    )
  }
  if(max(negative_region) > 0){
    precision_blocks[[length(precision_blocks) + 1]] <- as.matrix(
      Compute_Prec_rev(a = a, c = c, k = k, region = negative_region,
                       accuracy = accuracy, boundary = boundary)
    )
  }

  B <- cbind(B_pos, B_neg)
  P <- if(length(precision_blocks) == 0){
    matrix(0, nrow = 0, ncol = 0)
  } else if(length(precision_blocks) == 1) precision_blocks[[1]] else as.matrix(Matrix::bdiag(precision_blocks))

  attr(B, "positive_support") <- sort(unique(domain_x[domain_x > 0]))
  attr(B, "negative_support") <- sort(unique(abs(domain_x[domain_x < 0])))
  list(B = B, P = P)
}

make_iwp_side_instance <- function(x, knots, response_var, smoothing_var,
                                   sd.prior, psd.prior, boundary.prior){
  data <- data.frame(.x = x, .y = 0)
  names(data) <- c(as.character(smoothing_var), as.character(response_var))

  methods::new(
    "iwp",
    response_var = as.name(response_var),
    smoothing_var = as.name(smoothing_var),
    order = 2,
    knots = knots,
    k = length(knots),
    observed_x = sort(x),
    sd.prior = sd.prior,
    psd.prior = psd.prior,
    boundary.prior = boundary.prior,
    data = data,
    X = matrix(0, nrow = 0, ncol = 0),
    B = matrix(0, nrow = 0, ncol = 0),
    P = matrix(0, nrow = 0, ncol = 0),
    initial_location = 0,
    region = range(knots)
  )
}

build_tiwp_side_design <- function(x, side = c("positive", "negative"), knots,
                                   response_var, smoothing_var, sd.prior,
                                   psd.prior, boundary.prior){
  side <- match.arg(side)
  side_x <- if(side == "positive") pmax(x, 0) else pmax(-x, 0)
  active_rows <- if(side == "positive") which(x > 0) else which(x < 0)

  if(length(knots) == 0 || max(knots) <= 0){
    return(matrix(0, nrow = length(x), ncol = 0))
  }

  active_instance <- make_iwp_side_instance(
    x = if(length(active_rows) > 0) side_x[active_rows] else 0,
    knots = knots,
    response_var = response_var,
    smoothing_var = smoothing_var,
    sd.prior = sd.prior,
    psd.prior = psd.prior,
    boundary.prior = boundary.prior
  )
  design_active <- as.matrix(local_poly(active_instance))
  design <- matrix(0, nrow = length(x), ncol = ncol(design_active))
  if(length(active_rows) > 0){
    design[active_rows, ] <- design_active
  }
  design
}

build_two_sided_tiwp2_fem <- function(observed_x, domain_x, k, response_var,
                                      smoothing_var, sd.prior, psd.prior,
                                      boundary.prior){
  positive_domain <- sort(unique(domain_x[domain_x > 0]))
  negative_domain <- sort(unique(abs(domain_x[domain_x < 0])))
  positive_knots <- if(length(positive_domain) > 0) unique(sort(seq(0, max(positive_domain), length.out = k))) else numeric()
  negative_knots <- if(length(negative_domain) > 0) unique(sort(seq(0, max(negative_domain), length.out = k))) else numeric()

  B_pos <- build_tiwp_side_design(
    x = observed_x,
    side = "positive",
    knots = positive_knots,
    response_var = response_var,
    smoothing_var = smoothing_var,
    sd.prior = sd.prior,
    psd.prior = psd.prior,
    boundary.prior = boundary.prior
  )
  B_neg <- build_tiwp_side_design(
    x = observed_x,
    side = "negative",
    knots = negative_knots,
    response_var = response_var,
    smoothing_var = smoothing_var,
    sd.prior = sd.prior,
    psd.prior = psd.prior,
    boundary.prior = boundary.prior
  )

  precision_blocks <- list()
  if(length(positive_knots) > 0 && max(positive_knots) > 0){
    precision_blocks[[length(precision_blocks) + 1]] <- as.matrix(
      compute_weights_precision(make_iwp_side_instance(
        x = positive_domain,
        knots = positive_knots,
        response_var = response_var,
        smoothing_var = smoothing_var,
        sd.prior = sd.prior,
        psd.prior = psd.prior,
        boundary.prior = boundary.prior
      ))
    )
  }
  if(length(negative_knots) > 0 && max(negative_knots) > 0){
    precision_blocks[[length(precision_blocks) + 1]] <- as.matrix(
      compute_weights_precision(make_iwp_side_instance(
        x = negative_domain,
        knots = negative_knots,
        response_var = response_var,
        smoothing_var = smoothing_var,
        sd.prior = sd.prior,
        psd.prior = psd.prior,
        boundary.prior = boundary.prior
      ))
    )
  }

  B <- cbind(B_pos, B_neg)
  P <- if(length(precision_blocks) == 0){
    matrix(0, nrow = 0, ncol = 0)
  } else if(length(precision_blocks) == 1) precision_blocks[[1]] else as.matrix(Matrix::bdiag(precision_blocks))
  attr(B, "transformed_knots_positive") <- positive_knots
  attr(B, "transformed_knots_negative") <- negative_knots
  attr(B, "positive_support") <- positive_domain
  attr(B, "negative_support") <- negative_domain
  list(B = B, P = P, knots = sort(unique(c(positive_knots, negative_knots))))
}
