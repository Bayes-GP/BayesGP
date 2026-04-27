make_exact_iwp_data <- function(){
  x <- seq(0.1, 2.0, length.out = 12)
  truth <- sqrt(x + 1)
  y <- truth + c(
    -0.03, 0.04, 0.01, -0.02, 0.03, -0.01,
    0.02, -0.04, 0.01, 0.03, -0.02, 0.02
  )

  data_full <- data.frame(x = x, y = y)
  data_train <- data_full[1:8, , drop = FALSE]

  list(data_full = data_full, data_train = data_train)
}

test_that("exact state-space iwp supports posterior samples and derivatives", {
  sim_data <- make_exact_iwp_data()

  fit <- model_fit(
    formula = y ~ f(
      x,
      model = "iwp",
      computation = "state-space",
      order = 2,
      initial_location = "left",
      grid = sim_data$data_full$x,
      sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
    ),
    data = sim_data$data_train,
    family = "gaussian",
    control.family = list(sd = 0.1),
    M = 20
  )

  pred <- predict(fit, newdata = sim_data$data_full, variable = "x", only.samples = TRUE)
  deriv_pred <- predict(fit, newdata = sim_data$data_full, variable = "x", deriv = 1, only.samples = TRUE)

  expect_equal(fit$exact_iwp_meta$x$model_name, "iwp")
  expect_equal(fit$exact_iwp_meta$x$computation_method, "state-space")
  expect_equal(dim(pred), c(nrow(sim_data$data_full), 21))
  expect_equal(dim(deriv_pred), c(nrow(sim_data$data_full), 21))
  expect_equal(pred$x, sim_data$data_full$x)
  expect_equal(deriv_pred$x, sim_data$data_full$x)
  expect_false(anyNA(pred[, -1, drop = FALSE]))
  expect_false(anyNA(deriv_pred[, -1, drop = FALSE]))
})

test_that("exact state-space iwp requires the support grid for out-of-sample prediction", {
  sim_data <- make_exact_iwp_data()

  fit <- model_fit(
    formula = y ~ f(
      x,
      model = "iwp",
      computation = "state-space",
      order = 2,
      initial_location = "left",
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

test_that("exact state-space iwp errors early on nearly singular grids", {
  dense_data <- data.frame(
    x = c(0.1, 0.2, 0.200000000001, 0.6),
    y = c(1.0, 1.1, 1.1, 1.3)
  )

  expect_error(
    model_fit(
      formula = y ~ f(
        x,
        model = "iwp",
        computation = "state-space",
        order = 2,
        initial_location = "left",
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
