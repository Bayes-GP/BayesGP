# BayesGP 0.1.3

* Bug fix:
- `model_fit()` now validates model input columns before fitting and reports
  missing or non-finite values with the offending column name, avoiding opaque
  downstream optimizer failures.
- Near-monotone `mgp` and `tiwp2` smooth terms now use model-specific S4
  classes instead of being reported as native `iwp` objects.
- `var_density()` and `var_plot()` now compute model-specific PSD densities for
  `mgp` and `tiwp2` terms using the original-scale starting location `x`; when
  omitted, `x` defaults to the PSD location stored in the fitted `sd.prior`.

# BayesGP 0.1.2 (2024-06-16)

* Bug fix: 
- Fixed a bug in the `global_poly_helper_sgp` function and subsequently in the `predict` method that caused the prediction to fail when the smoothing variable is not pre-initialized manually. This should not affect models and predictions with manually initialized smoothing variables.

* Minor changes:
- The definitions of `initial_location` being `left`, `right` and `middle` in IWP and sGP are not defined based on the region rather than the provided smoothing variable.
