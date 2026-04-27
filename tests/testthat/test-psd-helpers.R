test_that("PSD_compute dispatches to the model-specific helpers", {
  expect_equal(
    PSD_compute(model = "mgp", x = 0.8, h = 1.7, c = 1.2, a = 2, sd = 1.4),
    PSD_compute_mgp(x = 0.8, h = 1.7, c = 1.2, a = 2, sd = 1.4)
  )

  expect_equal(
    PSD_compute(model = "tiwp2", x = 0.8, h = 1.7, c = 1.2, a = 2, sd = 1.4),
    PSD_compute_tiwp2(x = 0.8, h = 1.7, c = 1.2, a = 2, sd = 1.4)
  )

  expect_equal(
    PSD_compute(model = "iwp", h = 1.7, p = 2, sd = 1.4),
    PSD_compute_iwp(h = 1.7, p = 2, sd = 1.4)
  )

  expect_equal(
    PSD_compute(model = "sgp", h = 1.7, a = 0.3, m = 2, sd = 1.4),
    PSD_compute_sgp(h = 1.7, a = 0.3, m = 2, sd = 1.4)
  )
})

test_that("PSD_compute validates model-specific inputs", {
  expect_error(
    PSD_compute(model = "mgp", h = 2, a = 2, c = 1),
    "`x` must be supplied"
  )

  expect_error(
    PSD_compute(model = "unknown", h = 2),
    "`model` must be one of"
  )
})

test_that("PSD helpers stay aligned with prior conversion utilities", {
  mgp_prior <- list(alpha = 0.1, u = 0.7)
  mgp_factor <- PSD_compute_mgp(x = 0.4, h = 1.3, c = 1.1, a = 2, sd = 1)
  expect_equal(
    prior_conversion_mgp(prior = mgp_prior, x = 0.4, h = 1.3, c = 1.1, a = 2)$u,
    mgp_prior$u / mgp_factor
  )
  expect_equal(
    prior_conversion_mgp(prior = mgp_prior, x = 0.4, h = 1.3, c = 1.1, a = 2)$alpha,
    mgp_prior$alpha
  )

  tiwp_prior <- list(alpha = 0.1, u = 0.7)
  tiwp_factor <- PSD_compute_tiwp2(x = 0.4, h = 1.3, c = 1.1, a = 2, sd = 1)
  expect_equal(
    prior_conversion_tiwp2(prior = tiwp_prior, x = 0.4, h = 1.3, c = 1.1, a = 2)$u,
    tiwp_prior$u / tiwp_factor
  )
  expect_equal(
    prior_conversion_tiwp2(prior = tiwp_prior, x = 0.4, h = 1.3, c = 1.1, a = 2)$alpha,
    tiwp_prior$alpha
  )

  iwp_prior <- list(alpha = 0.1, u = 0.7)
  iwp_factor <- PSD_compute_iwp(h = 1.3, p = 2, sd = 1)
  expect_equal(
    prior_conversion_iwp(d = 1.3, prior = iwp_prior, p = 2)$u,
    iwp_prior$u / iwp_factor
  )

  sgp_prior <- list(alpha = 0.1, u = 0.7)
  sgp_factor <- PSD_compute_sgp(h = 1.3, a = 0.3, m = 2, sd = 1)
  expect_equal(
    prior_conversion_sgp(d = 1.3, prior = sgp_prior, a = 0.3, m = 2)$u,
    sgp_prior$u / sgp_factor
  )
})
