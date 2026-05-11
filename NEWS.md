# BayesGP (development version)

* Bug fix:
- `model_fit()` now validates model input columns before fitting and reports
  missing or non-finite values with the offending column name, avoiding opaque
  downstream optimizer failures.

# BayesGP 0.1.2 (2024-06-16)

* Bug fix: 
- Fixed a bug in the `global_poly_helper_sgp` function and subsequently in the `predict` method that caused the prediction to fail when the smoothing variable is not pre-initialized manually. This should not affect models and predictions with manually initialized smoothing variables.

* Minor changes:
- The definitions of `initial_location` being `left`, `right` and `middle` in IWP and sGP are not defined based on the region rather than the provided smoothing variable.
