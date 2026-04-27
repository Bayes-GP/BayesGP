test_that("native mgp FEM terms work in case-crossover models", {
  data_cc <- data.frame(
    y = c(1, 0, 0, 1, 0, 0),
    x = c(0.2, 0.5, 0.9, 1.2, 1.5, 1.8),
    strata = c("A", "A", "A", "B", "B", "B")
  )

  fit <- model_fit(
    formula = y ~ f(
      x,
      model = "mgp",
      computation = "fem",
      a = 2,
      c = 1,
      k = 8,
      normalized_boundary = FALSE,
      sd.prior = list(prior = "exp", param = list(u = 0.2, alpha = 0.5))
    ),
    data = data_cc,
    family = "cc",
    strata = "strata",
    M = 20,
    method = "aghq"
  )

  pred <- predict(fit, newdata = data_cc, variable = "x", only.samples = TRUE)
  summary_df <- predict(fit, newdata = data_cc, variable = "x", only.samples = FALSE)

  expect_equal(dim(pred), c(nrow(data_cc), 21))
  expect_equal(summary_df$x, data_cc$x)
  expect_false(anyNA(pred[, -1, drop = FALSE]))
  expect_false(anyNA(summary_df$mean))
})

test_that("native tiwp2 FEM terms work in case-crossover models", {
  data_cc <- data.frame(
    y = c(1, 0, 0, 1, 0, 0),
    x = c(0.2, 0.5, 0.9, 1.2, 1.5, 1.8),
    strata = c("A", "A", "A", "B", "B", "B")
  )

  fit <- model_fit(
    formula = y ~ f(
      x,
      model = "tiwp2",
      computation = "fem",
      initial_location = 1,
      a = 2,
      c = 0,
      k = 8,
      normalized_boundary = TRUE,
      sd.prior = list(prior = "exp", param = list(u = 0.2, alpha = 0.5))
    ),
    data = data_cc,
    family = "cc",
    strata = "strata",
    M = 20,
    method = "aghq"
  )

  pred <- predict(fit, newdata = data_cc, variable = "x", only.samples = TRUE)
  summary_df <- predict(fit, newdata = data_cc, variable = "x", only.samples = FALSE)

  expect_equal(dim(pred), c(nrow(data_cc), 21))
  expect_equal(summary_df$x, data_cc$x)
  expect_false(anyNA(pred[, -1, drop = FALSE]))
  expect_false(anyNA(summary_df$mean))
})

test_that("case-crossover models support multinomial counts within strata", {
  data_cc <- data.frame(
    y = c(2, 1, 3, 1, 2, 2),
    x = c(0.2, 0.5, 0.9, 1.2, 1.5, 1.8),
    strata = c("A", "A", "A", "B", "B", "B")
  )

  fit <- model_fit(
    formula = y ~ f(
      x,
      model = "mgp",
      computation = "fem",
      a = 2,
      c = 1,
      k = 8,
      normalized_boundary = FALSE,
      sd.prior = list(prior = "exp", param = list(u = 0.2, alpha = 0.5))
    ),
    data = data_cc,
    family = "cc",
    strata = "strata",
    M = 20,
    method = "aghq"
  )

  pred <- predict(fit, newdata = data_cc, variable = "x", only.samples = TRUE)

  expect_equal(dim(pred), c(nrow(data_cc), 21))
  expect_false(anyNA(pred[, -1, drop = FALSE]))
})

test_that("native IWP FEM case-crossover fits print with default knots", {
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
    M = 20,
    method = "aghq"
  )

  printed <- paste(capture.output(print(fit)), collapse = "\n")

  expect_match(printed, "BayesGP FitResult")
  expect_match(printed, "Smooth term specification")
  expect_equal(fit$instances[[1]]@k, length(fit$instances[[1]]@knots))
})

test_that("case-crossover fixed effects use all positive counts in a stratum", {
  data_cc <- data.frame(
    y = c(1, 2, 8),
    x = c(0, 1, 2),
    strata = c("A", "A", "A")
  )

  fit <- model_fit(
    formula = y ~ x,
    data = data_cc,
    family = "cc",
    strata = "strata",
    method = "nlminb",
    M = 20,
    control.fixed = list(x = list(mean = 0, prec = 0.001))
  )

  expect_gt(as.numeric(fit$mod$mean[1]), 0.5)
})
