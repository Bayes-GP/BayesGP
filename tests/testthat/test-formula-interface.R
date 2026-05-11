test_that("parse_formula keeps near-monotone computation arguments", {
  parse_result <- parse_formula(
    y ~ f(
      x,
      model = "mgp",
      computation = "state-space",
      a = curvature_value,
      c = shift_value
    ) + z
  )

  expect_equal(parse_result$response, as.symbol("y"))
  expect_equal(parse_result$rand_effects[[1]][[2]], as.symbol("x"))
  expect_equal(parse_result$rand_effects[[1]]$model, "mgp")
  expect_equal(parse_result$rand_effects[[1]]$computation, "state-space")
  expect_equal(parse_result$rand_effects[[1]]$a, as.symbol("curvature_value"))
  expect_equal(parse_result$rand_effects[[1]]$c, as.symbol("shift_value"))
  expect_equal(parse_result$fixed_effects[[1]], as.symbol("z"))
})

test_that("f stores computation for term-level representation choice", {
  term_computation <- f(x, model = "mgp", computation = "state-space", a = 2, c = 1)
  expect_equal(term_computation$model, "mgp")
  expect_equal(term_computation$computation, "state-space")
})

test_that("f does not override model-specific initial_location defaults", {
  term_default <- f(x, model = "mgp", computation = "state-space", a = 2, c = 1)
  term_explicit <- f(
    x,
    model = "mgp",
    computation = "state-space",
    a = 2,
    c = 1,
    initial_location = "middle"
  )

  expect_false("initial_location" %in% names(term_default))
  expect_equal(term_explicit$initial_location, "middle")
})

test_that("model_fit fails early on missing response values", {
  data_sim <- data.frame(x = c(0.1, 0.3, 0.6), y = c(1, NA, 1.2))

  expect_error(
    model_fit(
      y ~ f(x, model = "iid", sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))),
      data = data_sim,
      family = "gaussian",
      control.family = list(sd = 0.1),
      M = 5
    ),
    "column 'y'"
  )
})

test_that("model_fit fails early on non-finite smooth covariates", {
  data_sim <- data.frame(x = c(0.1, NaN, 0.6), y = c(1, 1.1, 1.2))

  expect_error(
    model_fit(
      y ~ f(x, model = "iid", sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))),
      data = data_sim,
      family = "gaussian",
      control.family = list(sd = 0.1),
      M = 5
    ),
    "column 'x'"
  )
})

test_that("model_fit fails early on missing fixed effects and offsets", {
  fixed_data <- data.frame(
    x = c(0.1, 0.3, 0.6),
    z = c(1, NA, 3),
    offset_value = c(0, 0.1, 0.2),
    y = c(1, 1.1, 1.2)
  )
  offset_data <- fixed_data
  offset_data$z <- c(1, 2, 3)
  offset_data$offset_value[2] <- NA

  expect_error(
    model_fit(
      y ~ f(x, model = "iid", sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))) +
        z + offset(offset_value),
      data = fixed_data,
      family = "gaussian",
      control.family = list(sd = 0.1),
      M = 5
    ),
    "column 'z'"
  )

  expect_error(
    model_fit(
      y ~ f(x, model = "iid", sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))) +
        z + offset(offset_value),
      data = offset_data,
      family = "gaussian",
      control.family = list(sd = 0.1),
      M = 5
    ),
    "column 'offset_value'"
  )
})

test_that("model_fit fails early on missing case-crossover family columns", {
  base_data <- data.frame(
    y = c(1, 0, 0, 1),
    x = c(0.1, 0.2, 0.3, 0.4),
    strata = c("A", "A", "B", "B"),
    weight = c(1, 0, 1, 0)
  )
  strata_data <- base_data
  strata_data$strata[2] <- NA
  weight_data <- base_data
  weight_data$weight[2] <- NaN

  expect_error(
    model_fit(
      y ~ f(x, model = "iid", sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))),
      data = strata_data,
      family = "cc",
      strata = "strata",
      weight = "weight",
      M = 5
    ),
    "column 'strata'"
  )

  expect_error(
    model_fit(
      y ~ f(x, model = "iid", sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))),
      data = weight_data,
      family = "cc",
      strata = "strata",
      weight = "weight",
      M = 5
    ),
    "column 'weight'"
  )
})

test_that("term-level method now errors and users must use computation", {
  data_sim <- data.frame(x = c(0.1, 0.3, 0.6), y = c(1, 1.1, 1.2))

  expect_error(
    model_fit(
      y ~ f(
        x,
        model = "mgp",
        method = "fem",
        a = 2,
        c = 1,
        sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
      ),
      data = data_sim,
      family = "gaussian",
      control.family = list(sd = 0.1),
      M = 20
    ),
    "Use `computation` inside `f\\(\\)`"
  )
})

test_that("unknown smooth-term models fail before fitting", {
  data_sim <- data.frame(x = c(0.1, 0.3, 0.6), y = c(1, 1.1, 1.2))

  expect_error(
    model_fit(
      y ~ f(
        x,
        model = "iwp2",
        computation = "fem",
        a = 2,
        c = 0
      ),
      data = data_sim,
      family = "gaussian",
      control.family = list(sd = 0.1),
      M = 20
    ),
    "Use `model = \"tiwp2\"`"
  )
})
