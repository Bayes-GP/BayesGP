#' List supported smooth-term models
#'
#' Returns the smooth-term model names accepted by \code{\link{f}}, the computation
#' methods available for each model, and the argument used to define the
#' prediction or basis domain.
#'
#' @return A data frame describing the model names accepted by \code{\link{f}} and their
#'   available computation methods.
#' @export
supported_models <- function(){
  data.frame(
    model = c("iid", "iwp", "sgp", "mgp", "tiwp2"),
    computation = c(NA, "fem, state-space", "fem", "fem, state-space", "fem, state-space"),
    description = c(
      "Independent random effect over factor levels.",
      "Integrated Wiener process smooth; FEM/O-spline is the default and exact state-space is available for order 2.",
      "Seasonal Gaussian process smooth.",
      "Near-monotone mGP smooth; FEM is the default and state-space is available explicitly.",
      "Transformed IWP2 near-monotone smooth; FEM is the default and state-space is available explicitly."
    ),
    domain_argument = c(NA, "knots / region", "region", "grid for state-space; region for fem", "grid for state-space; region for fem"),
    stringsAsFactors = FALSE
  )
}
