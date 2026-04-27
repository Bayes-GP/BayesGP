make_usability_data <- function(){
  x <- seq(0.1, 1.2, length.out = 8)
  data.frame(
    x = x,
    z = seq(0.2, 1.6, length.out = 8),
    group = factor(rep(letters[1:4], each = 2)),
    y = sqrt(x + 1) + c(-0.03, 0.04, 0.01, -0.02, 0.03, -0.01, 0.02, -0.04)
  )
}

known_sd_control <- function(){
  list(sd = 0.1)
}

exp_sd_prior <- function(){
  list(prior = "exp", param = list(u = 1, alpha = 0.5))
}

expect_no_stale_nearmono_columns <- function(fit){
  terms <- summary(fit)$smooth_terms
  printed <- paste(capture.output(print(fit)), collapse = "\n")
  summary_printed <- paste(capture.output(print(summary(fit))), collapse = "\n")

  expect_false("curvature" %in% names(terms))
  expect_false(grepl("curvature", printed, fixed = TRUE))
  expect_false(grepl("curvature", summary_printed, fixed = TRUE))
}

test_that("iid FitResult output is tailored to level random effects", {
  data <- make_usability_data()
  fit <- model_fit(
    y ~ f(group, model = "iid", sd.prior = exp_sd_prior()),
    data = data,
    family = "gaussian",
    control.family = known_sd_control(),
    M = 5
  )

  terms <- summary(fit)$smooth_terms

  expect_equal(names(terms), c("component", "model", "n_levels"))
  expect_equal(terms$model, "iid")
  expect_equal(terms$n_levels, 4)
  expect_no_stale_nearmono_columns(fit)
})

test_that("FitResult summary keeps Gaussian SD fields only for Gaussian fits", {
  data <- make_usability_data()
  fit <- model_fit(
    y ~ f(group, model = "iid", sd.prior = exp_sd_prior()),
    data = data,
    family = "gaussian",
    control.family = known_sd_control(),
    M = 5
  )

  fit_info <- summary(fit)$fit

  expect_true(all(c("gaussian_sd_known", "gaussian_sd_value") %in% names(fit_info)))
  expect_true(isTRUE(fit_info$gaussian_sd_known))
  expect_equal(fit_info$gaussian_sd_value, 0.1)
})

test_that("case-crossover FitResult summary reports strata instead of Gaussian SD fields", {
  data_cc <- data.frame(
    y = c(1, 0, 0, 1, 0, 0),
    x = c(0.2, 0.5, 0.9, 1.2, 1.5, 1.8),
    strata = c("A", "A", "A", "B", "B", "B")
  )

  fit <- model_fit(
    formula = y ~ f(
      x,
      model = "iwp",
      order = 2,
      sd.prior = list(prior = "exp", param = list(u = 0.2, alpha = 0.5))
    ),
    data = data_cc,
    family = "cc",
    strata = "strata",
    aghq_k = 3,
    M = 5,
    method = "aghq"
  )

  fit_summary <- expect_warning(summary(fit), NA)
  fit_info <- fit_summary$fit
  summary_printed <- paste(capture.output(print(fit_summary)), collapse = "\n")

  expect_false(any(grepl("gaussian_sd", names(fit_info), fixed = TRUE)))
  expect_false(grepl("gaussian_sd", summary_printed, fixed = TRUE))
  expect_true(all(c("strata", "n_strata", "max_stratum_size", "count_source", "total_count") %in% names(fit_info)))
  expect_equal(fit_info$strata, "strata")
  expect_equal(fit_info$n_strata, 2)
  expect_equal(fit_info$max_stratum_size, 3)
  expect_equal(fit_info$count_source, "y")
  expect_equal(fit_info$total_count, 2)
})

test_that("native IWP FEM FitResult output omits near-monotone-only fields", {
  data <- make_usability_data()
  fit <- model_fit(
    y ~ f(x, model = "iwp", order = 2, k = 8, sd.prior = exp_sd_prior()),
    data = data,
    family = "gaussian",
    control.family = known_sd_control(),
    M = 5
  )

  terms <- summary(fit)$smooth_terms

  expect_equal(terms$model, "iwp")
  expect_equal(terms$computation, "fem")
  expect_equal(terms$order, 2)
  expect_equal(terms$basis_size, 8)
  expect_true("domain" %in% names(terms))
  expect_false("a" %in% names(terms))
  expect_false("c" %in% names(terms))
  expect_no_stale_nearmono_columns(fit)
})

test_that("exact IWP state-space FitResult output reports support size", {
  data <- make_usability_data()
  fit <- model_fit(
    y ~ f(
      x,
      model = "iwp",
      computation = "state-space",
      order = 2,
      grid = data$x,
      sd.prior = exp_sd_prior()
    ),
    data = data,
    family = "gaussian",
    control.family = known_sd_control(),
    M = 5
  )

  terms <- summary(fit)$smooth_terms

  expect_equal(terms$model, "iwp")
  expect_equal(terms$computation, "state-space")
  expect_equal(terms$order, 2)
  expect_equal(terms$support_size, length(unique(data$x)) - 1)
  expect_false("a" %in% names(terms))
  expect_false("c" %in% names(terms))
  expect_no_stale_nearmono_columns(fit)
})

test_that("SGP FitResult output reports seasonal settings", {
  data <- make_usability_data()
  fit <- model_fit(
    y ~ f(x, model = "sgp", period = 1, m = 2, k = 8, sd.prior = exp_sd_prior()),
    data = data,
    family = "gaussian",
    control.family = known_sd_control(),
    M = 5
  )

  terms <- summary(fit)$smooth_terms

  expect_equal(terms$model, "sgp")
  expect_equal(terms$computation, "fem")
  expect_equal(terms$harmonics, 2)
  expect_equal(terms$period, 1)
  expect_true(all(c("angular_frequency", "frequency", "period") %in% names(terms)))
  expect_false("c" %in% names(terms))
  expect_no_stale_nearmono_columns(fit)
})

test_that("near-monotone FitResult output uses a and c, not curvature", {
  data <- make_usability_data()
  fit_specs <- list(
    list(model = "mgp", computation = "state-space", normalized_boundary = FALSE),
    list(model = "tiwp2", computation = "fem", normalized_boundary = TRUE)
  )

  for(spec in fit_specs){
    fit <- model_fit(
      y ~ f(
        x,
        model = spec$model,
        computation = spec$computation,
        a = 2,
        c = 1,
        k = 8,
        grid = data$x,
        normalized_boundary = spec$normalized_boundary,
        sd.prior = exp_sd_prior()
      ),
      data = data,
      family = "gaussian",
      control.family = known_sd_control(),
      M = 5
    )

    terms <- summary(fit)$smooth_terms

    expect_equal(terms$model, spec$model, info = spec$model)
    expect_equal(terms$computation, spec$computation, info = spec$model)
    expect_equal(terms$a, 2, info = spec$model)
    expect_equal(terms$c, 1, info = spec$model)
    expect_false("normalized_boundary" %in% names(terms), info = spec$model)
    expect_false("boundary" %in% names(terms), info = spec$model)
    expect_false("allow_zero_c" %in% names(terms), info = spec$model)
    expect_no_stale_nearmono_columns(fit)
  }
})

test_that("near-monotone terms default to FEM computation", {
  data <- make_usability_data()
  fit_specs <- list(
    list(model = "mgp", normalized_boundary = FALSE),
    list(model = "tiwp2", normalized_boundary = TRUE)
  )

  for(spec in fit_specs){
    fit <- model_fit(
      y ~ f(
        x,
        model = spec$model,
        a = 2,
        c = 1,
        k = 8,
        normalized_boundary = spec$normalized_boundary,
        sd.prior = exp_sd_prior()
      ),
      data = data,
      family = "gaussian",
      control.family = known_sd_control(),
      M = 5
    )

    terms <- summary(fit)$smooth_terms

    expect_equal(fit$near_mono_meta$x$computation_method, "fem", info = spec$model)
    expect_equal(terms$computation, "fem", info = spec$model)
  }
})

test_that("mixed smooth-term summaries keep only columns used by the fit", {
  data <- make_usability_data()
  fit <- model_fit(
    y ~ f(
      x,
      model = "mgp",
      computation = "state-space",
      a = 2,
      c = 1,
      grid = data$x,
      normalized_boundary = FALSE,
      sd.prior = exp_sd_prior()
    ) + f(group, model = "iid", sd.prior = exp_sd_prior()),
    data = data,
    family = "gaussian",
    control.family = known_sd_control(),
    M = 5
  )

  terms <- summary(fit)$smooth_terms

  expect_equal(nrow(terms), 2)
  expect_true(all(c("component", "model", "a", "c", "n_levels") %in% names(terms)))
  expect_false("curvature" %in% names(terms))
  expect_no_stale_nearmono_columns(fit)
})

test_that("model_arguments is no longer part of the public or internal API", {
  expect_false("model_arguments" %in% getNamespaceExports("BayesGP"))
  expect_false(exists("model_arguments", envir = asNamespace("BayesGP"), inherits = FALSE))
})
