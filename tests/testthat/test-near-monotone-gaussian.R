make_near_mono_data <- function(){
  x <- seq(0.1, 2.0, length.out = 12)
  truth <- sqrt(x + 1)
  y <- truth + c(
    -0.03, 0.04, 0.01, -0.02, 0.03, -0.01,
    0.02, -0.04, 0.01, 0.03, -0.02, 0.02
  )

  data_full <- data.frame(x = x, y = y, truth = truth)
  data_train <- data_full[1:8, c("x", "y"), drop = FALSE]

  list(data_full = data_full, data_train = data_train)
}

make_zero_left_near_mono_data <- function(){
  x <- seq(0, 2.0, length.out = 12)
  truth <- sqrt(x + 1e-6)
  y <- truth + c(
    -0.03, 0.04, 0.01, -0.02, 0.03, -0.01,
    0.02, -0.04, 0.01, 0.03, -0.02, 0.02
  )

  data_full <- data.frame(x = x, y = y, truth = truth)
  data_train <- data_full[1:8, c("x", "y"), drop = FALSE]

  list(data_full = data_full, data_train = data_train)
}

fit_near_mono_model <- function(model_name, computation_method,
                                normalized_boundary = TRUE,
                                sd_prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))){
  sim_data <- make_near_mono_data()

  model_fit(
    formula = y ~ f(
      x,
      model = model_name,
      computation = computation_method,
      a = 2,
      c = 1,
      k = 12,
      grid = sim_data$data_full$x,
      normalized_boundary = normalized_boundary,
      sd.prior = sd_prior
    ),
    data = sim_data$data_train,
    family = "gaussian",
    control.family = list(sd = 0.1),
    M = 20
  )
}

test_that("known-sd gaussian path works for exact and fem near-monotone terms", {
  sim_data <- make_near_mono_data()
  fit_specs <- list(
    list(name = "mgp_exact", model = "mgp", computation = "state-space", normalized_boundary = FALSE),
    list(name = "tiwp_exact", model = "tiwp2", computation = "state-space", normalized_boundary = TRUE),
    list(name = "mgp_fem", model = "mgp", computation = "fem", normalized_boundary = FALSE),
    list(name = "tiwp_fem", model = "tiwp2", computation = "fem", normalized_boundary = TRUE)
  )

  for(spec in fit_specs){
    fit <- fit_near_mono_model(
      model_name = spec$model,
      computation_method = spec$computation,
      normalized_boundary = spec$normalized_boundary
    )

    pred <- predict(fit, newdata = sim_data$data_full, variable = "x", only.samples = TRUE)
    summary_df <- predict(fit, newdata = sim_data$data_full, variable = "x", only.samples = FALSE)
    smooth_df <- smooth_summary(fit, component = "x", refined_x = sim_data$data_full$x)

    expect_true(isTRUE(fit$family_sd_known), info = spec$name)
    expect_true(isTRUE(fit$control.family$sd_known), info = spec$name)
    expect_equal(fit$family_sd_value, 0.1, tolerance = 1e-12, info = spec$name)
    expect_equal(dim(pred), c(nrow(sim_data$data_full), 21), info = spec$name)
    expect_equal(summary_df$x, sim_data$data_full$x, info = spec$name)
    expect_equal(smooth_df$x, sim_data$data_full$x, info = spec$name)
    expect_false(anyNA(pred[, -1, drop = FALSE]), info = spec$name)
    expect_false(anyNA(summary_df$mean), info = spec$name)
    expect_false(anyNA(smooth_df$mean), info = spec$name)
    expect_error(var_density(fit), "fixed", info = spec$name)
  }
})

test_that("state-space predictions require the support grid when extrapolating exact fits", {
  sim_data <- make_near_mono_data()
  fit <- model_fit(
    formula = y ~ f(
      x,
      model = "mgp",
      computation = "state-space",
      a = 2,
      c = 1,
      sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
    ),
    data = sim_data$data_train,
    family = "gaussian",
    control.family = list(sd = 0.1),
    M = 20
  )

  expect_error(
    predict(fit, newdata = sim_data$data_full, variable = "x", only.samples = TRUE),
    "support"
  )
})

test_that("near-monotone FEM region supports prediction outside the observed range", {
  sim_data <- make_near_mono_data()
  fit_specs <- list(
    list(name = "mgp_region", model = "mgp", normalized_boundary = FALSE, argument_name = "region"),
    list(name = "tiwp2_range", model = "tiwp2", normalized_boundary = TRUE, argument_name = "range")
  )

  for(spec in fit_specs){
    fit <- if(spec$argument_name == "region"){
      model_fit(
        formula = y ~ f(
          x,
          model = "mgp",
          computation = "fem",
          a = 2,
          c = 1,
          k = 12,
          normalized_boundary = FALSE,
          sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5)),
          region = range(sim_data$data_full$x)
        ),
        data = sim_data$data_train,
        family = "gaussian",
        control.family = list(sd = 0.1),
        M = 20
      )
    } else {
      model_fit(
        formula = y ~ f(
          x,
          model = "tiwp2",
          computation = "fem",
          a = 2,
          c = 1,
          k = 12,
          normalized_boundary = TRUE,
          sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5)),
          range = range(sim_data$data_full$x)
        ),
        data = sim_data$data_train,
        family = "gaussian",
        control.family = list(sd = 0.1),
        M = 20
      )
    }

    pred <- predict(fit, newdata = sim_data$data_full, variable = "x", only.samples = TRUE)

    expect_equal(pred$x, sim_data$data_full$x, info = spec$name)
    expect_false(anyNA(pred[, -1, drop = FALSE]), info = spec$name)
  }
})

test_that("FitResult print and summary expose fit diagnostics", {
  fit <- fit_near_mono_model(
    model_name = "mgp",
    computation_method = "fem",
    normalized_boundary = FALSE
  )

  printed <- paste(capture.output(print(fit)), collapse = "\n")
  fit_summary <- summary(fit)
  summary_printed <- paste(capture.output(print(fit_summary)), collapse = "\n")

  expect_match(printed, "Log marginal likelihood")
  expect_match(summary_printed, "Log marginal likelihood")
  expect_true(is.finite(fit_summary$fit$log_marginal_likelihood))
  expect_true("parameters" %in% names(fit_summary))
  expect_true(nrow(fit_summary$smooth_terms) >= 1)
})

test_that("near-monotone FEM and state-space log marginal likelihoods are comparable", {
  sim_data <- make_near_mono_data()
  fit_specs <- expand.grid(
    model = c("mgp", "tiwp2"),
    reference = c("left", "middle"),
    stringsAsFactors = FALSE
  )

  fit_one <- function(model_name, computation_method, reference_kind){
    c_original <- 1.25

    if(computation_method == "state-space"){
      return(model_fit(
        formula = y ~ f(
          x,
          model = model_name,
          computation = "state-space",
          a = 2,
          c = c_original,
          initial_location = reference_kind,
          grid = sim_data$data_full$x,
          sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
        ),
        data = sim_data$data_train,
        family = "gaussian",
        control.family = list(sd = 0.1),
        aghq_k = 3,
        M = 10
      ))
    }

    model_fit(
      formula = y ~ f(
        x,
        model = model_name,
        computation = "fem",
        a = 2,
        c = c_original,
        initial_location = reference_kind,
        k = 60,
        accuracy = 0.4 / 60,
        region = range(sim_data$data_full$x),
        sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
      ),
      data = sim_data$data_train,
      family = "gaussian",
      control.family = list(sd = 0.1),
      aghq_k = 3,
      M = 10
    )
  }

  for(row_index in seq_len(nrow(fit_specs))){
    spec <- fit_specs[row_index, ]
    exact_fit <- fit_one(spec$model, "state-space", spec$reference)
    fem_fit <- fit_one(spec$model, "fem", spec$reference)

    lml_diff <- abs(
      exact_fit$mod$normalized_posterior$lognormconst -
        fem_fit$mod$normalized_posterior$lognormconst
    )

    expect_true(lml_diff < 0.02, info = paste(spec$model, spec$reference))
  }
})

test_that("near-monotone terms reject zero original-scale transformed arguments by default", {
  sim_data <- make_zero_left_near_mono_data()
  fit_specs <- expand.grid(
    model = c("mgp", "tiwp2"),
    computation = c("state-space", "fem"),
    stringsAsFactors = FALSE
  )

  for(row_index in seq_len(nrow(fit_specs))){
    spec <- fit_specs[row_index, ]
    fit_call <- if(spec$computation == "state-space"){
      quote(model_fit(
        formula = y ~ f(
          x,
          model = spec$model,
          computation = spec$computation,
          a = 2,
          c = 0,
          initial_location = "left",
          grid = sim_data$data_full$x,
          sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
        ),
        data = sim_data$data_train,
        family = "gaussian",
        control.family = list(sd = 0.1),
        M = 10
      ))
    } else {
      quote(model_fit(
        formula = y ~ f(
          x,
          model = spec$model,
          computation = spec$computation,
          a = 2,
          c = 0,
          initial_location = "left",
          region = range(sim_data$data_full$x),
          k = 12,
          sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
        ),
        data = sim_data$data_train,
        family = "gaussian",
        control.family = list(sd = 0.1),
        M = 10
      ))
    }

    expect_error(eval(fit_call), "must satisfy `x \\+ c > 0`",
                 info = paste(spec$model, spec$computation))
  }
})

test_that("c equal to zero is allowed when the original-scale domain stays positive", {
  sim_data <- make_near_mono_data()
  fit <- model_fit(
    formula = y ~ f(
      x,
      model = "mgp",
      computation = "state-space",
      a = 2,
      c = 0,
      initial_location = "middle",
      grid = sim_data$data_full$x,
      sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
    ),
    data = sim_data$data_train,
    family = "gaussian",
    control.family = list(sd = 0.1),
    M = 10
  )

  expect_equal(fit$near_mono_meta$x$c, 0)
  expect_gt(fit$near_mono_meta$x$internal_shift, 0)
})

test_that("explicit c is an original-scale additive base-model shift", {
  sim_data <- make_near_mono_data()
  explicit_c <- 10
  reference <- mean(range(sim_data$data_full$x))
  fit <- model_fit(
    formula = y ~ f(
      x,
      model = "mgp",
      computation = "state-space",
      a = 2,
      c = explicit_c,
      initial_location = reference,
      grid = sim_data$data_full$x,
      normalized_boundary = TRUE,
      sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
    ),
    data = sim_data$data_train,
    family = "gaussian",
    control.family = list(sd = 0.1),
    M = 10
  )

  metadata <- fit$near_mono_meta$x
  expected_basis <- boundary_basis_matrix(
    x = sim_data$data_train$x,
    a = 2,
    c = reference + explicit_c,
    initial_location = reference,
    normalized = TRUE
  )

  expect_false(isTRUE(metadata$c_is_default))
  expect_equal(metadata$c, explicit_c)
  expect_equal(metadata$internal_shift, reference + explicit_c)
  expect_equal(fit$instances[[1]]@X, expected_basis, tolerance = 1e-12)
})

test_that("near-monotone terms choose a domain-aware c when c is omitted", {
  zero_data <- make_zero_left_near_mono_data()
  negative_data_full <- data.frame(
    x = seq(-1, 1, length.out = 12),
    y = seq(-0.3, 0.4, length.out = 12)
  )
  negative_data_train <- negative_data_full[1:8, , drop = FALSE]

  exact_fit <- model_fit(
    formula = y ~ f(
      x,
      model = "mgp",
      computation = "state-space",
      a = 2,
      normalized_boundary = FALSE,
      initial_location = "left",
      grid = zero_data$data_full$x,
      sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
    ),
    data = zero_data$data_train,
    family = "gaussian",
    control.family = list(sd = 0.1),
    M = 10
  )
  fem_fit <- model_fit(
    formula = y ~ f(
      x,
      model = "tiwp2",
      computation = "fem",
      a = 2,
      normalized_boundary = TRUE,
      initial_location = "left",
      region = range(negative_data_full$x),
      k = 12,
      accuracy = 0.01,
      sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
    ),
    data = negative_data_train,
    family = "gaussian",
    control.family = list(sd = 0.1),
    M = 10
  )

  exact_meta <- exact_fit$near_mono_meta$x
  fem_meta <- fem_fit$near_mono_meta$x
  exact_pred <- predict(exact_fit, newdata = zero_data$data_full, variable = "x")
  fem_pred <- predict(fem_fit, newdata = negative_data_full, variable = "x")

  expect_true(isTRUE(exact_meta$c_is_default))
  expect_true(isTRUE(fem_meta$c_is_default))
  expect_gt(exact_meta$c, 0)
  expect_gt(fem_meta$c, 1)
  expect_gt(min(zero_data$data_full$x) + exact_meta$c, 0)
  expect_gt(min(negative_data_full$x) + fem_meta$c, 0)
  expect_equal(exact_meta$internal_shift, exact_meta$reference_location + exact_meta$c)
  expect_equal(fem_meta$internal_shift, fem_meta$reference_location + fem_meta$c)
  expect_equal(ncol(exact_fit$instances[[1]]@X), 1)
  expect_equal(ncol(fem_fit$instances[[1]]@X), 1)
  expect_false(anyNA(exact_pred$mean))
  expect_false(anyNA(fem_pred$mean))
})

test_that("allow_zero_c is no longer a supported near-monotone option", {
  sim_data <- make_zero_left_near_mono_data()

  expect_error(
    model_fit(
      formula = y ~ f(
        x,
        model = "tiwp2",
        computation = "state-space",
        a = 2,
        allow_zero_c = TRUE,
        initial_location = "left",
        grid = sim_data$data_full$x,
        sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
      ),
      data = sim_data$data_train,
      family = "gaussian",
      control.family = list(sd = 0.1),
      M = 10
    ),
    "`allow_zero_c` is no longer supported"
  )
  expect_error(
    model_fit(
      formula = y ~ f(
        x,
        model = "mgp",
        computation = "state-space",
        a = 2,
        allow_zero_c = FALSE,
        initial_location = "left",
        grid = sim_data$data_full$x,
        sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
      ),
      data = sim_data$data_train,
      family = "gaussian",
      control.family = list(sd = 0.1),
      M = 10
    ),
    "`allow_zero_c` is no longer supported"
  )
})

test_that("supported model helper includes near-monotone models", {
  models <- supported_models()

  expect_true("mgp" %in% models$model)
  expect_true("tiwp2" %in% models$model)
  expect_match(
    models$description[models$model == "mgp"],
    "FEM is the default",
    fixed = TRUE
  )
  expect_match(
    models$computation[models$model == "tiwp2"],
    "fem, state-space",
    fixed = TRUE
  )
  expect_match(
    models$domain_argument[models$model == "tiwp2"],
    "grid for state-space; region for fem",
    fixed = TRUE
  )
})

test_that("near-monotone default c is zero when the original-scale domain is positive", {
  sim_data <- make_near_mono_data()

  mgp_fit <- model_fit(
    formula = y ~ f(
      x,
      model = "mgp",
      computation = "state-space",
      a = 2,
      grid = sim_data$data_full$x,
      sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
    ),
    data = sim_data$data_train,
    family = "gaussian",
    control.family = list(sd = 0.1),
    M = 10
  )
  tiwp_fit <- model_fit(
    formula = y ~ f(
      x,
      model = "tiwp2",
      computation = "state-space",
      a = 2,
      grid = sim_data$data_full$x,
      sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
    ),
    data = sim_data$data_train,
    family = "gaussian",
    control.family = list(sd = 0.1),
    M = 10
  )

  expect_equal(mgp_fit$near_mono_meta$x$c, 0)
  expect_equal(tiwp_fit$near_mono_meta$x$c, 0)
  expect_true(isTRUE(mgp_fit$near_mono_meta$x$c_is_default))
  expect_true(isTRUE(tiwp_fit$near_mono_meta$x$c_is_default))
  expect_equal(
    mgp_fit$near_mono_meta$x$internal_shift,
    mgp_fit$near_mono_meta$x$reference_location
  )
  expect_equal(
    tiwp_fit$near_mono_meta$x$internal_shift,
    tiwp_fit$near_mono_meta$x$reference_location
  )
})

test_that("normalized near-monotone boundary coefficients are level and slope at the reference", {
  sim_data <- make_near_mono_data()
  c_original <- 1.25
  curvature <- 2
  finite_difference_step <- 1e-5
  reference_values <- c(
    min(sim_data$data_full$x),
    mean(range(sim_data$data_full$x)),
    max(sim_data$data_full$x)
  )
  fit_specs <- expand.grid(
    model = c("mgp", "tiwp2"),
    computation = c("state-space", "fem"),
    reference_index = seq_along(reference_values),
    stringsAsFactors = FALSE
  )

  normalized_original_basis <- function(x, reference){
    lambda <- (curvature - 1) / curvature
    reference_arg <- reference + c_original
    ((x + c_original)^lambda - reference_arg^lambda) /
      (lambda * reference_arg^(lambda - 1))
  }

  for(row_index in seq_len(nrow(fit_specs))){
    spec <- fit_specs[row_index, ]
    reference <- reference_values[spec$reference_index]
    fit_call <- if(spec$computation == "state-space"){
      quote(model_fit(
        formula = y ~ f(
          x,
          model = spec$model,
          computation = spec$computation,
          a = curvature,
          c = c_original,
          initial_location = reference,
          normalized_boundary = TRUE,
          grid = sim_data$data_full$x,
          boundary.prior = list(mean = 0.7, prec = 10),
          sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
        ),
        data = sim_data$data_train,
        family = "gaussian",
        control.family = list(sd = 0.1),
        M = 10
      ))
    } else {
      quote(model_fit(
        formula = y ~ f(
          x,
          model = spec$model,
          computation = spec$computation,
          a = curvature,
          c = c_original,
          initial_location = reference,
          normalized_boundary = TRUE,
          region = range(sim_data$data_full$x),
          k = 20,
          boundary.prior = list(mean = 0.7, prec = 10),
          sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
        ),
        data = sim_data$data_train,
        family = "gaussian",
        control.family = list(sd = 0.1),
        M = 10
      ))
    }
    fit <- eval(fit_call)
    meta <- fit$near_mono_meta$x
    boundary_at_reference <- boundary_basis_matrix(
      x = reference,
      a = curvature,
      c = meta$internal_shift,
      initial_location = reference,
      normalized = TRUE
    )
    boundary_derivative <- (
      boundary_basis_matrix(
        x = reference + finite_difference_step,
        a = curvature,
        c = meta$internal_shift,
        initial_location = reference,
        normalized = TRUE
      ) -
        boundary_basis_matrix(
          x = reference - finite_difference_step,
          a = curvature,
          c = meta$internal_shift,
          initial_location = reference,
          normalized = TRUE
        )
    ) / (2 * finite_difference_step)

    expect_equal(meta$c, c_original, info = paste(spec$model, spec$computation, reference))
    expect_equal(meta$internal_shift, reference + c_original, info = paste(spec$model, spec$computation, reference))
    expect_equal(as.numeric(boundary_at_reference), 0, tolerance = 1e-10,
                 info = paste(spec$model, spec$computation, reference))
    expect_equal(as.numeric(boundary_derivative), 1, tolerance = 1e-8,
                 info = paste(spec$model, spec$computation, reference))
    expect_equal(
      as.numeric(fit$instances[[1]]@X),
      normalized_original_basis(sim_data$data_train$x, reference),
      tolerance = 1e-10,
      info = paste(spec$model, spec$computation, reference)
    )
  }
})

test_that("near-monotone FEM supports non-left reference locations", {
  sim_data <- make_near_mono_data()
  fit_specs <- list(
    list(name = "mgp", model = "mgp", computation = "fem", normalized_boundary = FALSE),
    list(name = "tiwp2", model = "tiwp2", computation = "fem", normalized_boundary = TRUE)
  )

  for(spec in fit_specs){
    fit <- model_fit(
      formula = y ~ f(
        x,
        model = spec$model,
        computation = spec$computation,
        a = 2,
        c = 1,
        k = 12,
        region = range(sim_data$data_full$x),
        initial_location = "middle",
        normalized_boundary = spec$normalized_boundary,
        sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
      ),
      data = sim_data$data_train,
      family = "gaussian",
      control.family = list(sd = 0.1),
      M = 20
    )

    pred <- predict(fit, newdata = sim_data$data_full, variable = "x", only.samples = TRUE)

    expect_true(isTRUE(fit$near_mono_meta$x$two_sided_reference), info = spec$name)
    expect_equal(pred$x, sim_data$data_full$x, info = spec$name)
    expect_false(anyNA(pred[, -1, drop = FALSE]), info = spec$name)
  }
})

test_that("exact tiwp2 supports non-left reference locations through two-sided support", {
  sim_data <- make_near_mono_data()

  fit <- model_fit(
    formula = y ~ f(
      x,
      model = "tiwp2",
      computation = "state-space",
      a = 2,
      c = 1,
      grid = sim_data$data_full$x,
      initial_location = "middle",
      normalized_boundary = TRUE,
      sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
    ),
    data = sim_data$data_train,
    family = "gaussian",
    control.family = list(sd = 0.1),
    M = 20
  )

  pred <- predict(fit, newdata = sim_data$data_full, variable = "x", only.samples = TRUE)

  expect_true(isTRUE(fit$near_mono_meta$x$two_sided_reference))
  expect_equal(pred$x, sim_data$data_full$x)
  expect_false(anyNA(pred[, -1, drop = FALSE]))
})

test_that("exact mgp supports non-left reference locations when reverse precision is available", {
  sim_data <- make_near_mono_data()

  fit <- model_fit(
    formula = y ~ f(
      x,
      model = "mgp",
      computation = "state-space",
      a = 2,
      c = 1,
      grid = sim_data$data_full$x,
      initial_location = "middle",
      normalized_boundary = FALSE,
      sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
    ),
    data = sim_data$data_train,
    family = "gaussian",
    control.family = list(sd = 0.1),
    M = 20
  )

  pred <- predict(fit, newdata = sim_data$data_full, variable = "x", only.samples = TRUE)

  expect_true(isTRUE(fit$near_mono_meta$x$two_sided_reference))
  expect_equal(pred$x, sim_data$data_full$x)
  expect_false(anyNA(pred[, -1, drop = FALSE]))
})

test_that("exact mgp non-left references fail clearly for unsupported reverse curvatures", {
  sim_data <- make_near_mono_data()

  expect_error(
    model_fit(
      formula = y ~ f(
        x,
        model = "mgp",
        computation = "state-space",
        a = 1.5,
        c = 1,
        grid = sim_data$data_full$x,
        initial_location = "middle",
        normalized_boundary = FALSE,
        sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
      ),
      data = sim_data$data_train,
      family = "gaussian",
      control.family = list(sd = 0.1),
      M = 20
    ),
    "Reverse exact mGP state-space currently supports"
  )
})

test_that("near-monotone terms convert PSD-scale sd.prior specifications internally", {
  prior_u <- 0.7
  prior_alpha <- 0.2
  prior_h <- 1.3
  prior_x <- 0.4

  fit_specs <- list(
    list(name = "mgp", model = "mgp", computation = "state-space", normalized_boundary = FALSE),
    list(name = "tiwp2", model = "tiwp2", computation = "state-space", normalized_boundary = TRUE)
  )

  for(spec in fit_specs){
    fit <- fit_near_mono_model(
      model_name = spec$model,
      computation_method = spec$computation,
      normalized_boundary = spec$normalized_boundary,
      sd_prior = list(
        prior = "exp",
        param = list(u = prior_u, alpha = prior_alpha),
        h = prior_h,
        x = prior_x
      )
    )

    instance <- fit$instances[[1]]
    metadata <- fit$near_mono_meta$x
    expected_u <- prior_u / PSD_compute(
      model = spec$model,
      x = prior_x - metadata$reference_location,
      h = prior_h,
      a = 2,
      c = metadata$internal_shift
    )

    expect_equal(instance@psd.prior$param$u, prior_u, info = spec$name)
    expect_equal(instance@psd.prior$param$alpha, prior_alpha, info = spec$name)
    expect_equal(instance@sd.prior$param$u, expected_u, tolerance = 1e-12, info = spec$name)
    expect_equal(instance@sd.prior$param$alpha, prior_alpha, info = spec$name)
    expect_equal(instance@sd.prior$step, prior_h, tolerance = 1e-12, info = spec$name)
    expect_equal(instance@sd.prior$x, prior_x, tolerance = 1e-12, info = spec$name)
    expect_equal(
      instance@sd.prior$centered_x,
      prior_x - metadata$reference_location,
      tolerance = 1e-12,
      info = spec$name
    )
    expect_null(instance@sd.prior$h, info = spec$name)
  }
})

test_that("near-monotone PSD-scale priors require an x location", {
  expect_error(
    fit_near_mono_model(
      model_name = "mgp",
      computation_method = "state-space",
      normalized_boundary = FALSE,
      sd_prior = list(
        prior = "exp",
        param = list(u = 0.7, alpha = 0.2),
        h = 1.3
      )
    ),
    "sd\\.prior\\$x"
  )
})

test_that("exact near-monotone state-space errors early on nearly singular grids", {
  dense_data <- data.frame(
    x = c(0.1, 0.2, 0.200000000001, 0.6),
    y = c(1.0, 1.1, 1.1, 1.3)
  )

  expect_error(
    model_fit(
      formula = y ~ f(
        x,
        model = "mgp",
        computation = "state-space",
        a = 2,
        c = 1,
        sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
      ),
      data = dense_data,
      family = "gaussian",
      control.family = list(sd = 0.1),
      M = 20
    ),
    "sparser `grid`|near-duplicates|computation = \"fem\""
  )
})

test_that("exact near-monotone state-space accepts well-spaced regular grids", {
  support_grid <- seq(0.1, 39.5, by = 0.1)

  precision_matrix <- expect_silent(
    build_exact_precision_with_diagnostics(
      support_grid = support_grid,
      model_name = "mgp",
      precision_builder = function(grid) mGP_joint_prec(t_vec = grid, a = 2, c = 0.3831318282944436)
    )
  )

  expect_true(is.matrix(precision_matrix))
  expect_false(anyNA(precision_matrix))
  expect_false(inherits(try(chol(precision_matrix), silent = TRUE), "try-error"))
})
