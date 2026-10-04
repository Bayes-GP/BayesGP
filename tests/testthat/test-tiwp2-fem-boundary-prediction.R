test_that("tiwp2 FEM predictions reconstruct the fitted boundary basis", {
  data_sim <- data.frame(x = seq(0.1, 2, length.out = 12))
  data_sim$y <- sqrt(data_sim$x + 1)
  fit_specs <- expand.grid(
    a = c(2, 1),
    normalized_boundary = c(TRUE, FALSE),
    initial_location = c("left", "middle", "right"),
    stringsAsFactors = FALSE
  )

  for(row_index in seq_len(nrow(fit_specs))){
    spec <- fit_specs[row_index, ]
    label <- paste("a =", spec$a, spec$initial_location, spec$normalized_boundary)
    fit <- model_fit(
      y ~ f(
        x,
        model = "tiwp2",
        computation = "fem",
        a = spec$a,
        c = 1,
        k = 12,
        initial_location = spec$initial_location,
        normalized_boundary = spec$normalized_boundary,
        boundary.prior = list(mean = 0.7, prec = 10)
      ),
      data = data_sim,
      family = "gaussian",
      control.family = list(sd = 0.1),
      M = 10
    )

    instance <- fit$instances[[1]]
    expected <- as.matrix(
      instance@B %*% extract_random_samples(fit, "x") +
        instance@X %*% extract_boundary_samples(fit, "x")
    )
    smooth_draws <- smooth_samples(fit, "x", refined_x = data_sim$x)
    expect_equal(
      unname(as.matrix(smooth_draws[, -1, drop = FALSE])),
      unname(expected),
      tolerance = 1e-10,
      info = label
    )

    for(include_intercept in c(FALSE, TRUE)){
      expected_prediction <- expected
      if(include_intercept){
        expected_prediction <- sweep(
          expected_prediction,
          2,
          as.numeric(extract_intercept_samples(fit)),
          "+"
        )
      }
      prediction <- predict(
        fit,
        newdata = data_sim,
        variable = "x",
        include.intercept = include_intercept,
        only.samples = TRUE
      )
      expect_equal(
        unname(as.matrix(prediction[, -1, drop = FALSE])),
        unname(expected_prediction),
        tolerance = 1e-10,
        info = paste(label, "include.intercept =", include_intercept)
      )
    }

    # Isolate the boundary contribution and check a new evaluation grid
    # against the analytic original-scale base function.
    boundary_only_fit <- fit
    boundary_only_fit$samps$samps[fit$random_samp_indexes$x, ] <- 0
    boundary_only_fit$samps$samps[fit$boundary_samp_indexes$x, ] <- 1
    evaluation_x <- seq(min(data_sim$x), max(data_sim$x), length.out = 17)
    reference_arg <- instance@initial_location + 1
    if(spec$a == 1){
      expected_boundary <- log(evaluation_x + 1)
      if(spec$normalized_boundary){
        expected_boundary <- reference_arg *
          (expected_boundary - log(reference_arg))
      }
    } else {
      lambda <- (spec$a - 1) / spec$a
      expected_boundary <- (evaluation_x + 1)^lambda
      if(spec$normalized_boundary){
        expected_boundary <- (expected_boundary - reference_arg^lambda) /
          (lambda * reference_arg^(lambda - 1))
      }
    }
    boundary_draws <- predict(
      boundary_only_fit,
      newdata = data.frame(x = evaluation_x),
      variable = "x",
      include.intercept = FALSE,
      only.samples = TRUE
    )
    expect_equal(
      unname(as.matrix(boundary_draws[, -1, drop = FALSE])),
      matrix(expected_boundary, nrow = length(evaluation_x), ncol = ncol(expected)),
      tolerance = 1e-10,
      info = paste(label, "new evaluation grid")
    )
  }
})
