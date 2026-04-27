#' @keywords internal
#' @import Matrix
#' @importFrom graphics legend matplot
#' @importFrom methods as new setGeneric setMethod slot slotNames
#' @importFrom stats density model.matrix nlminb predict quantile terms
#' @importFrom utils capture.output type.convert
"_PACKAGE"

#' Declare a BayesGP smooth or random-effect term
#'
#' Helper used inside \code{\link{model_fit}} formulas to declare smooth or
#' random-effect terms. Use this page for the shared interface and common prior
#' syntax; use the model-specific help pages for model arguments and examples.
#'
#' @param smoothing_var Column in `data` defining the random-effect or smooth
#'   input.
#' @param model Smooth-term model name. Supported values are `"iid"`, `"iwp"`,
#'   `"sgp"`, `"mgp"`, and `"tiwp2"`.
#' @param computation Optional representation choice. FEM is the default for
#'   models with FEM support; exact state-space representations are available
#'   for selected `iwp`, `mgp`, and `tiwp2` terms.
#' @param sd.prior Prior specification for the smooth standard deviation. The
#'   common exponential PC prior form is `list(prior = "exp", param = list(u =
#'   1, alpha = 0.5))`, meaning `P(SD > u) = alpha`.
#' @param boundary.prior Prior specification for boundary/global coefficients
#'   when the model has such coefficients. A common form is `list(mean = 0,
#'   prec = 0.001)`.
#' @param initial_location Optional reference location. Supported values depend
#'   on the model and include numeric values and named locations such as
#'   `"left"`, `"middle"`, and `"right"`.
#' @param ... Additional model-specific arguments.
#'
#' @section Supported models:
#' \describe{
#'   \item{`model = "iid"`}{Independent random effects over factor levels. See
#'   \code{\link{f_iid}}.}
#'   \item{`model = "iwp"`}{Integrated Wiener process smooth using the BayesGP
#'   FEM/O-spline implementation by default, with an added exact state-space
#'   representation for `order = 2`. See \code{\link{f_iwp}}.}
#'   \item{`model = "sgp"`}{Seasonal Gaussian process smooth using the BayesGP
#'   FEM implementation. See \code{\link{f_sgp}}.}
#'   \item{`model = "mgp"`}{Near-monotone mGP smooth. FEM is the default;
#'   state-space computation is available explicitly. See \code{\link{f_mgp}}.}
#'   \item{`model = "tiwp2"`}{Transformed IWP2 near-monotone smooth. FEM is the
#'   default; state-space computation is available explicitly. See
#'   \code{\link{f_tiwp2}}.}
#' }
#'
#' @section Common argument examples:
#' \preformatted{
#' # Standard exponential PC prior for a smooth SD:
#' sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))
#'
#' # Calibrate the prior through h-step predictive SD where supported:
#' sd.prior = list(
#'   prior = "exp",
#'   param = list(u = 0.5, alpha = 0.1),
#'   h = 1
#' )
#'
#' # Boundary/global coefficient prior:
#' boundary.prior = list(mean = 0, prec = 0.001)
#'
#' # FEM basis domain and exact state-space support:
#' region = c(0, 10)
#' grid = seq(0, 10, length.out = 50)
#' }
#'
#' @examples
#' f(group, model = "iid")
#' f(x, model = "iwp", order = 2, k = 30)
#' f(time, model = "sgp", period = 12, m = 2)
#' f(x, model = "mgp", a = 2, initial_location = "left")
#' f(x, model = "tiwp2", a = 2, computation = "state-space", grid = seq(0, 1, length.out = 20))
#'
#' @seealso \code{\link{f_iid}}, \code{\link{f_iwp}},
#'   \code{\link{f_sgp}}, \code{\link{f_mgp}}, \code{\link{f_tiwp2}},
#'   \code{\link{model_fit}}, \code{\link{supported_models}}
#' @export
f <- function(smoothing_var, model = "iid", computation = NULL,
              sd.prior = NULL, boundary.prior = NULL,
              initial_location = NULL, ...) {
  # Capture the full call
  mc <- match.call(expand.dots = TRUE)
  
  # Replace the first argument with the function name
  mc[[1]] <- as.name("f")
  
  # Ensure all default arguments are included
  mc$model <- model
  mc$computation <- computation
  mc$sd.prior <- sd.prior
  mc$boundary.prior <- boundary.prior
  if(!is.null(initial_location)){
    mc$initial_location <- initial_location[1]
  }
  
  # Replace smoothing_var with its unevaluated form
  mc$smoothing_var <- substitute(smoothing_var)
  
  # Return the modified call
  return(mc)
}

#' Model-specific `f()` documentation for independent random effects
#'
#' Use `f(group, model = "iid")` for exchangeable random effects over observed
#' levels of a grouping variable.
#'
#' @section Required arguments:
#' \describe{
#'   \item{`smoothing_var`}{Grouping variable in `data`. Character, factor, and
#'   numeric inputs are converted to observed levels.}
#'   \item{`model`}{Use `model = "iid"`.}
#' }
#'
#' @section Optional arguments:
#' \describe{
#'   \item{`sd.prior`}{Prior for the random-effect standard deviation. The
#'   default is `list(prior = "exp", param = list(u = 1, alpha = 0.5))`.}
#' }
#'
#' @examples
#' f(group, model = "iid")
#' f(group, model = "iid",
#'   sd.prior = list(prior = "exp", param = list(u = 0.5, alpha = 0.1))
#' )
#'
#' @seealso \code{\link{f}}, \code{\link{supported_models}}
#' @name f_iid
#' @aliases f.iid
NULL

#' Model-specific `f()` documentation for integrated Wiener process smooths
#'
#' Use `f(x, model = "iwp")` for BayesGP integrated Wiener process smooths. The
#' default is the BayesGP FEM/O-spline representation. An exact state-space
#' representation is also available for `order = 2`.
#'
#' @section Required arguments:
#' \describe{
#'   \item{`smoothing_var`}{Numeric smoothing variable in `data`.}
#'   \item{`model`}{Use `model = "iwp"`.}
#'   \item{`order`}{IWP order. The FEM representation supports positive
#'   orders; exact state-space currently requires `order = 2`.}
#' }
#'
#' @section Computation and domain arguments:
#' \describe{
#'   \item{`computation`}{Use `"fem"` for the default FEM/O-spline
#'   representation, or `"state-space"` for the exact order-2 representation.}
#'   \item{`k`}{Number of FEM knots when explicit `knots` are not supplied.
#'   Defaults to 30.}
#'   \item{`knots`}{Explicit FEM knot sequence on the original input scale.}
#'   \item{`region`}{Continuous FEM domain. Defaults to the observed range.}
#'   \item{`grid`}{Discrete support for exact state-space fitting and later
#'   prediction. Defaults to observed `smoothing_var` values.}
#'   \item{`initial_location`}{Reference location. FEM defaults to `"middle"`;
#'   exact state-space defaults to `"left"`.}
#' }
#'
#' @section Prior arguments:
#' \describe{
#'   \item{`sd.prior`}{Prior for the smooth standard deviation. Add `h` or
#'   `step` for predictive-SD calibration.}
#'   \item{`boundary.prior`}{Gaussian prior for boundary/global coefficients,
#'   such as `list(mean = 0, prec = 0.001)`.}
#' }
#'
#' @examples
#' f(x, model = "iwp", order = 2, k = 30)
#' f(x, model = "iwp", order = 2, computation = "state-space", grid = seq(0, 1, length.out = 20))
#' f(x, model = "iwp", order = 3, region = c(0, 10),
#'   sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5)),
#'   boundary.prior = list(mean = 0, prec = 0.001)
#' )
#'
#' @seealso \code{\link{f}}, \code{\link{f_sgp}},
#'   \code{\link{f_mgp}}, \code{\link{f_tiwp2}}
#' @name f_iwp
#' @aliases f.iwp
NULL

#' Model-specific `f()` documentation for seasonal Gaussian process smooths
#'
#' Use `f(x, model = "sgp")` for BayesGP seasonal Gaussian process smooths.
#' Seasonal GP terms use the BayesGP FEM representation.
#'
#' @section Required arguments:
#' \describe{
#'   \item{`smoothing_var`}{Numeric smoothing variable in `data`.}
#'   \item{`model`}{Use `model = "sgp"`.}
#'   \item{`period`, `freq`, or `a`}{Seasonal frequency. Supply exactly one of
#'   period length, cycles-per-unit frequency, or angular frequency `a`.}
#' }
#'
#' @section Computation and domain arguments:
#' \describe{
#'   \item{`m`}{Number of harmonics. Defaults to 1.}
#'   \item{`k`}{Number of B-spline basis functions. Defaults to 30.}
#'   \item{`region`}{Continuous FEM domain. Defaults to the observed range.}
#'   \item{`accuracy`}{Numerical integration resolution. Defaults to 5000.}
#'   \item{`boundary`}{Whether to drop boundary B-spline basis functions.
#'   Defaults to `TRUE`.}
#' }
#'
#' @section Prior arguments:
#' \describe{
#'   \item{`sd.prior`}{Prior for the smooth standard deviation. Add `h` or
#'   `step` for predictive-SD calibration.}
#'   \item{`boundary.prior`}{Gaussian prior for seasonal boundary/global
#'   coefficients, such as `list(mean = 0, prec = 0.001)`.}
#' }
#'
#' @examples
#' f(month, model = "sgp", period = 12, m = 2)
#' f(time, model = "sgp", freq = 1 / 365, k = 40, region = c(0, 365))
#' f(time, model = "sgp", a = 2 * pi / 12,
#'   sd.prior = list(prior = "exp", param = list(u = 0.5, alpha = 0.1)),
#'   boundary.prior = list(mean = 0, prec = 0.001)
#' )
#'
#' @seealso \code{\link{f}}, \code{\link{f_iwp}}
#' @name f_sgp
#' @aliases f.sgp
NULL

#' Model-specific `f()` documentation for near-monotone mGP smooths
#'
#' Use `f(x, model = "mgp")` for near-monotone mGP smooths. FEM is the default
#' computation method; exact state-space computation is available explicitly.
#'
#' @section Required arguments:
#' \describe{
#'   \item{`smoothing_var`}{Numeric smoothing variable in `data`.}
#'   \item{`model`}{Use `model = "mgp"`.}
#'   \item{`a` or `alpha`}{Curvature-control parameter for the base model.
#'   Supply one name.}
#' }
#'
#' @section Computation and domain arguments:
#' \describe{
#'   \item{`computation`}{Use `"fem"` by default, or `"state-space"` for exact
#'   support-grid computation.}
#'   \item{`c`}{Original-scale additive shift in
#'   `m(x) = (x + c)^((a - 1) / a)`, or `m(x) = log(x + c)` when `a = 1`. If
#'   omitted, BayesGP chooses a domain-aware shift so `x + c > 0`.}
#'   \item{`initial_location`}{Reference location: `"left"`, `"middle"`,
#'   `"right"`, or a numeric value inside the domain. Defaults to `"left"`.}
#'   \item{`region` or `range`}{Continuous FEM domain. Defaults to the observed
#'   range.}
#'   \item{`grid`}{Discrete support for exact state-space fitting and later
#'   prediction.}
#'   \item{`k`}{Number of FEM basis functions. Defaults to 30.}
#'   \item{`accuracy`}{FEM integration step. Defaults to 0.01.}
#'   \item{`normalized_boundary`}{Whether the shared boundary basis is scaled
#'   so its derivative is one at the reference. Defaults to `TRUE`.}
#' }
#'
#' @section Prior arguments:
#' \describe{
#'   \item{`sd.prior`}{Prior for the smooth standard deviation. Add `h` or
#'   `step` plus original-scale `x` for predictive-SD calibration.}
#'   \item{`boundary.prior`}{Gaussian prior for the shared boundary
#'   coefficient. The model intercept supplies the shared level when the
#'   family includes an intercept.}
#' }
#'
#' @examples
#' f(pm25, model = "mgp", a = 2)
#' f(pm25, model = "mgp", a = 2, c = 10, initial_location = "middle")
#' f(pm25, model = "mgp", a = 2, region = c(0, 40), k = 40)
#' f(pm25, model = "mgp", a = 2, computation = "state-space", grid = seq(0, 40, length.out = 50))
#' f(pm25, model = "mgp", a = 2,
#'   sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5), h = 1, x = 10),
#'   boundary.prior = list(mean = 0, prec = 0.001)
#' )
#'
#' @seealso \code{\link{f}}, \code{\link{f_tiwp2}}
#' @name f_mgp
#' @aliases f.mgp
NULL

#' Model-specific `f()` documentation for transformed IWP2 near-monotone smooths
#'
#' Use `f(x, model = "tiwp2")` for transformed IWP2 near-monotone smooths.
#' FEM is the default computation method; exact state-space computation is
#' available explicitly.
#'
#' @section Required arguments:
#' \describe{
#'   \item{`smoothing_var`}{Numeric smoothing variable in `data`.}
#'   \item{`model`}{Use `model = "tiwp2"`.}
#'   \item{`a` or `alpha`}{Curvature-control parameter for the Box-Cox base
#'   model. Supply one name.}
#' }
#'
#' @section Computation and domain arguments:
#' \describe{
#'   \item{`computation`}{Use `"fem"` by default, or `"state-space"` for exact
#'   support-grid computation.}
#'   \item{`c`}{Original-scale additive shift in the Box-Cox base model. If
#'   omitted, BayesGP chooses a domain-aware shift so `x + c > 0`.}
#'   \item{`initial_location`}{Reference location: `"left"`, `"middle"`,
#'   `"right"`, or a numeric value inside the domain. Defaults to `"left"`.}
#'   \item{`region` or `range`}{Continuous FEM domain. Defaults to the observed
#'   range.}
#'   \item{`grid`}{Discrete support for exact state-space fitting and later
#'   prediction.}
#'   \item{`k`}{Number of FEM basis functions. Defaults to 30.}
#'   \item{`accuracy`}{FEM integration step. Defaults to 0.01.}
#'   \item{`normalized_boundary`}{Whether the shared boundary basis is scaled
#'   so its derivative is one at the reference. Defaults to `TRUE`.}
#' }
#'
#' @section Prior arguments:
#' \describe{
#'   \item{`sd.prior`}{Prior for the smooth standard deviation. Add `h` or
#'   `step` plus original-scale `x` for predictive-SD calibration.}
#'   \item{`boundary.prior`}{Gaussian prior for the shared boundary
#'   coefficient. The model intercept supplies the shared level when the
#'   family includes an intercept.}
#' }
#'
#' @examples
#' f(pm25, model = "tiwp2", a = 2)
#' f(pm25, model = "tiwp2", a = 2, c = 10, initial_location = "middle")
#' f(pm25, model = "tiwp2", a = 2, region = c(0, 40), k = 40)
#' f(pm25, model = "tiwp2", a = 2, computation = "state-space", grid = seq(0, 40, length.out = 50))
#' f(pm25, model = "tiwp2", a = 2,
#'   sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5), h = 1, x = 10),
#'   boundary.prior = list(mean = 0, prec = 0.001)
#' )
#'
#' @seealso \code{\link{f}}, \code{\link{f_mgp}}
#' @name f_tiwp2
#' @aliases f.tiwp2
NULL

parse_formula <- function(formula) {
  components <- as.list(attributes(terms(formula))$ variables)
  fixed_effects <- list()
  rand_effects <- list()
  offset_effects <- list()
  # Index starts as 3 since index 1 represents "list" and
  # index 2 represents the response variable
  for (i in 3:length(components)) {
    if (startsWith(toString(components[[i]]), "f,")) {
      rand_effects[[length(rand_effects) + 1]] <- components[[i]]
    } 
    else if(startsWith(toString(components[[i]]), "offset,")){
      offset_effects[[length(offset_effects) + 1]] <- components[[i]]
    }
    else {
      fixed_effects[[length(fixed_effects) + 1]] <- components[[i]]
    }
  }
  return(list(response = components[[2]], fixed_effects = fixed_effects, rand_effects = rand_effects, offset_effects = offset_effects))
}

# Create a class for iwp using S4
setClass("iwp", slots = list(
  response_var = "name", smoothing_var = "name", order = "numeric",
  knots = "numeric", k = "numeric",
  observed_x = "numeric", sd.prior = "list",
  psd.prior = "list",
  boundary.prior = "list", data = "data.frame", X = "matrix",
  B = "matrix", P = "matrix", initial_location = "ANY",
  region = "numeric"
))

# Create a class for sgp using S4
setClass("sgp", slots = list(
  response_var = "name", smoothing_var = "name", 
  a = "numeric", freq = "numeric", period = "numeric",
  m = "numeric", k = "numeric",
  knots = "numeric", observed_x = "numeric", sd.prior = "list",
  psd.prior = "list",
  boundary.prior = "list", data = "data.frame", X = "matrix",
  B = "matrix", P = "matrix", initial_location = "ANY", region = "numeric", 
  accuracy = "numeric", boundary = "logical"
))

# Create a class for iid using S4
setClass("iid", slots = list(
  response_var = "name", smoothing_var = "name", sd.prior = "list",
  data = "data.frame", B = "matrix", P = "matrix"
))

# Create a class for customized using S4
setClass("customized", slots = list(
  response_var = "name", smoothing_var = "name", sd.prior = "list",
  data = "data.frame", B = "matrix", P = "matrix",
  compute_B = "function", compute_P = "function"
))

#' @import Matrix
#' @importClassesFrom Matrix dgCMatrix
Compute_Q_sB <- function(a,k,region, accuracy = 5000, boundary = TRUE){
  ss <- function(M) {Matrix::forceSymmetric(M + Matrix::t(M))}
  if((accuracy %% 1) == 0){
    x <- seq(min(region),max(region), length.out = accuracy)
  }
  else{
    x <- seq(min(region),max(region), by = accuracy)
  }
  if(boundary == TRUE){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4,
                                                          dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4))
  }
  Bmatrix <- fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE)
  B1matrix <-  fda::eval.basis(x, B_basis, Lfdobj=1, returnMatrix=TRUE)
  B2matrix <-  fda::eval.basis(x, B_basis, Lfdobj=2, returnMatrix=TRUE)
  cos_matrix <- cos(a*x)
  sin_matrix <- sin(a*x)
  Bcos <- as(apply(Bmatrix, 2, function(x) x*cos_matrix), "dgCMatrix")
  B1cos <- as(apply(B1matrix, 2, function(x) x*cos_matrix), "dgCMatrix")
  B2cos <- as(apply(B2matrix, 2, function(x) x*cos_matrix), "dgCMatrix")
  Bsin <- as(apply(Bmatrix, 2, function(x) x*sin_matrix), "dgCMatrix")
  B1sin <- as(apply(B1matrix, 2, function(x) x*sin_matrix), "dgCMatrix")
  B2sin <- as(apply(B2matrix, 2, function(x) x*sin_matrix), "dgCMatrix")
  
  ### Compute I, L, T:
  Numerical_I <- as(diag(c(diff(c(0,x)))), "dgCMatrix")
  
  ### T
  T00 <- Matrix::t(Bcos) %*% Numerical_I %*% Bcos
  T10 <- Matrix::t(B1cos) %*% Numerical_I %*% Bcos
  T11 <- Matrix::t(B1cos) %*% Numerical_I %*% B1cos
  T20 <- Matrix::t(B2cos) %*% Numerical_I %*% Bcos
  T21 <- Matrix::t(B2cos) %*% Numerical_I %*% B1cos
  T22 <- Matrix::t(B2cos) %*% Numerical_I %*% B2cos
  
  ### L
  L00 <- Matrix::t(Bsin) %*% Numerical_I %*% Bsin
  L10 <- Matrix::t(B1sin) %*% Numerical_I %*% Bsin
  L11 <- Matrix::t(B1sin) %*% Numerical_I %*% B1sin
  L20 <- Matrix::t(B2sin) %*% Numerical_I %*% Bsin
  L21 <- Matrix::t(B2sin) %*% Numerical_I %*% B1sin
  L22 <- Matrix::t(B2sin) %*% Numerical_I %*% B2sin
  
  ### I
  I00 <- Matrix::t(Bsin) %*% Numerical_I %*% Bcos
  I10 <- Matrix::t(B1sin) %*% Numerical_I %*% Bcos
  I11 <- Matrix::t(B1sin) %*% Numerical_I %*% B1cos
  I20 <- Matrix::t(B2sin) %*% Numerical_I %*% Bcos
  I21 <- Matrix::t(B2sin) %*% Numerical_I %*% B1cos
  I22 <- Matrix::t(B2sin) %*% Numerical_I %*% B2cos
  
  
  ### Inner product involving B spline:
  Bmatrix <- as(Bmatrix, "dgCMatrix")
  B1matrix <- as(B1matrix, "dgCMatrix")
  B2matrix <- as(B2matrix, "dgCMatrix")
  
  BB <- Matrix::t(Bmatrix) %*% Numerical_I %*% Bmatrix
  B2B2 <- Matrix::t(B2matrix) %*% Numerical_I %*% B2matrix
  BB2 <- Matrix::t(Bmatrix) %*% Numerical_I %*% B2matrix
  BS <- Matrix::t(Bmatrix) %*% Numerical_I %*% Bsin
  BC <- Matrix::t(Bmatrix) %*% Numerical_I %*% Bcos
  BS1 <- Matrix::t(Bmatrix) %*% Numerical_I %*% B1sin
  BC1 <- Matrix::t(Bmatrix) %*% Numerical_I %*% B1cos
  BS2 <- Matrix::t(Bmatrix) %*% Numerical_I %*% B2sin
  BC2 <- Matrix::t(Bmatrix) %*% Numerical_I %*% B2cos
  B2S <- Matrix::t(B2matrix) %*% Numerical_I %*% Bsin
  B2C <- Matrix::t(B2matrix) %*% Numerical_I %*% Bcos
  B2S1 <- Matrix::t(B2matrix) %*% Numerical_I %*% B1sin
  B2C1 <- Matrix::t(B2matrix) %*% Numerical_I %*% B1cos
  B2S2 <- Matrix::t(B2matrix) %*% Numerical_I %*% B2sin
  B2C2 <- Matrix::t(B2matrix) %*% Numerical_I %*% B2cos
  
  ## G = <phi,phj>
  G <- rbind(cbind(T00, Matrix::t(I00), Matrix::t(BC)), cbind(I00,L00, Matrix::t(BS)), cbind(BC,BS,BB))
  
  ## C = <D^2phi,D^2phj>
  C11 <- T22 - 2*a*ss(I21) - (a^2)*ss(T20) + 2*(a^3)*ss(I10) + 4 * (a^2) * L11 + (a^4)*T00
  C22 <- L22 + 2*a*ss(I21) - (a^2)*ss(L20) - 2*(a^3)*ss(I10) + 4 * (a^2) * T11 + (a^4)*L00
  C12 <- I22 + 2*a*T21 - (a^2)* ss(I20) - 2*a*Matrix::t(L21) - 4*(a^2)*I11 + 2*(a^3)*L10 - 2*(a^3)*Matrix::t(T10) + (a^4)*I00
  
  C13 <- Matrix::t(B2C2) - 2*a*Matrix::t(B2S1) - (a^2)*Matrix::t(B2C)
  C23 <- Matrix::t(B2S2) + 2*a*Matrix::t(B2C1) - (a^2)*Matrix::t(B2S)
  C33 <- B2B2
  C <- rbind(cbind(C11,C12,C13), cbind(Matrix::t(C12), C22, C23), cbind(Matrix::t(C13), Matrix::t(C23), C33))
  
  
  ## M = <phi,D^2phj>
  M11 <- Matrix::t(T20) - (2*a)*Matrix::t(I10) - (a^2)*T00
  M12 <- Matrix::t(I20) + (2*a)*Matrix::t(T10) - (a^2)*I00
  M21 <- Matrix::t(I20) - (2*a)*Matrix::t(L10) - (a^2)*I00
  M22 <- Matrix::t(L20) + (2*a)*Matrix::t(I10) - (a^2)*L00
  
  M13 <- Matrix::t(B2C)
  M23 <- Matrix::t(B2S)
  M31 <- BC2 - (2*a)*BS1 - (a^2)*BC
  M32 <- BS2 + (2*a)*BC1 - (a^2)*BS
  M33 <- BB2
  
  M <- rbind(cbind(M11,M12,M13), cbind(M21,M22,M23), cbind(M31,M32,M33))
  
  
  ### Compute the final precision matrix: Q
  Q <- (a^4)*G + C + (a^2)*ss(M)
  Matrix::forceSymmetric(Q)
}


Compute_B_sB <- function(x, a, k, region, boundary = TRUE){
  if(boundary){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4,
                                                          dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4))
  }
  Bmatrix <- fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE)
  cos_matrix <- cos(a*x)
  sin_matrix <- sin(a*x)
  Bcos <- apply(Bmatrix, 2, function(x) x*cos_matrix)
  Bsin <- apply(Bmatrix, 2, function(x) x*sin_matrix)
  cbind(Bcos, Bsin,Bmatrix)
}


Compute_B_sB_helper <- function(refined_x, a, k, m, region, boundary = TRUE, initial_location = NULL){
  if(is.null(initial_location)){
    initial_location <- min(refined_x)
  }
  refined_x <- refined_x - initial_location
  region <- region - initial_location
  B <- NULL
  for (i in 1:m) {
    B <- cbind(B, Compute_B_sB(x = refined_x, a = (i*a), region = region, k = k, boundary = boundary))
  }
  B
}


setGeneric("compute_B", function(object) {
  standardGeneric("compute_B")
})
setMethod("compute_B", signature = "iid", function(object) {
  smoothing_var <- object@smoothing_var
  x <- as.factor((object@data)[[smoothing_var]])
  B <- model.matrix(~ -1 + x)
  B
})
setMethod("compute_B", signature = "customized", function(object) {
  smoothing_var <- object@smoothing_var
  object@compute_B((object@data)[[smoothing_var]])
})
setMethod("compute_B", signature = "sgp", function(object) {
  smoothing_var <- object@smoothing_var
  a <- object@a
  k <- object@k
  m <- object@m
  initial_location <- object@initial_location
  region <- object@region - initial_location
  accuracy <- object@accuracy
  boundary <- object@boundary
  refined_x <- (object@data)[[smoothing_var]] - initial_location
  B <- NULL
  for (i in 1:m) {
    B <- cbind(B, Compute_B_sB(x = refined_x, a = (i*a), region = region, k = k))
  }
  B
})


setGeneric("compute_P", function(object) {
  standardGeneric("compute_P")
})
setMethod("compute_P", signature = "iid", function(object) {
  smoothing_var <- object@smoothing_var
  x <- (object@data)[[smoothing_var]]
  num_factor <- length(unique(x))
  diag(nrow = num_factor, ncol = num_factor)
})
setMethod("compute_P", signature = "customized", function(object) {
  x <- (object@data)[[object@smoothing_var]]
  object@compute_P(x)
})
setMethod("compute_P", signature = "sgp", function(object) {
  smoothing_var <- object@smoothing_var
  initial_location <- object@initial_location
  a <- object@a
  k <- object@k
  m <- object@m
  region <- object@region - initial_location
  accuracy <- object@accuracy
  boundary <- object@boundary
  refined_x <- (object@data)[[smoothing_var]] - initial_location
  Q <- Compute_Q_sB(a = a, k = k, region = region, accuracy = accuracy)
  if(m >= 2){
    for (i in 2:m) {
      Q <- bdiag(Q, Compute_Q_sB(a = (i*a), k = k, region = region, accuracy = accuracy))
    }
  }
  as(Q, "matrix")
})


setGeneric("local_poly", function(object) {
  standardGeneric("local_poly")
})
setMethod("local_poly", signature = "iwp", function(object) {
  knots <- object@knots
  initial_location <- object@initial_location
  smoothing_var <- object@smoothing_var
  refined_x <- (object@data)[[smoothing_var]] - initial_location
  p <- object@order
  D <- local_poly_helper(knots = knots, refined_x = refined_x, p = p, neg_sign_order = 0)
  D # Local poly design matrix
})

setGeneric("global_poly", function(object) {
  standardGeneric("global_poly")
})
setMethod("global_poly", signature = "iwp", function(object) {
  smoothing_var <- object@smoothing_var
  x <- (object@data)[[smoothing_var]] - object@initial_location
  p <- object@order
  result <- NULL
  for (i in 1:p) {
    result <- cbind(result, x^(i - 1))
  }
  result
})
setMethod("global_poly", signature = "sgp", function(object) {
  smoothing_var <- object@smoothing_var
  a <- object@a
  m <- object@m
  initial_location <- object@initial_location
  refined_x <- (object@data)[[smoothing_var]] - initial_location
  X <- NULL
  for (i in 1:m) {
    X <- cbind(X, cos(i*a*refined_x),sin(i*a*refined_x))
  }
  X 
})


#' Constructing the precision matrix given the knot sequence
#'
#' @param object An `iwp` smooth-term object.
#' @return A precision matrix of the corresponding basis function, should be diagonal matrix with
#' size (k-1) by (k-1).
#' @export
setGeneric("compute_weights_precision", function(object) {
  standardGeneric("compute_weights_precision")
})
#' @rdname compute_weights_precision
#' @aliases compute_weights_precision,iwp-method
setMethod("compute_weights_precision", signature = "iwp", function(object) {
  knots <- object@knots
  if (min(knots) >= 0) {
    as(diag(diff(knots)), "matrix")
  } else if (max(knots) < 0) {
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    as(diag(diff(knots_neg)), "matrix")
  } else {
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    knots_pos <- knots
    knots_pos <- unique(sort(ifelse(knots > 0, knots, 0)))
    d1 <- diff(knots_neg)
    d2 <- diff(knots_pos)
    Precweights1 <- diag(d1)
    Precweights2 <- diag(d2)
    as(Matrix::bdiag(Precweights1, Precweights2), "matrix")
  }
})

#' Constructing the precision matrix given the knot sequence (helper)
#'
#' @param x A vector of knots used to construct the O-spline basis, first knot should be viewed as "0",
#' the reference starting location. These k knots will define (k-1) basis function in total.
#' @return A precision matrix of the corresponding basis function, should be diagonal matrix with
#' size (k-1) by (k-1).
#' @examples
#' compute_weights_precision_helper(x = c(0,0.2,0.4,0.6,0.8))
#' @export
compute_weights_precision_helper <- function(x){
  d <- diff(x)
  Precweights <- diag(d)
  Precweights
}



get_local_poly <- function(knots, refined_x, p) {
  dif <- diff(knots)
  nn <- length(refined_x)
  n <- length(knots)
  D <- matrix(0, nrow = nn, ncol = n - 1)
  for (j in 1:nn) {
    for (i in 1:(n - 1)) {
      if (refined_x[j] <= knots[i]) {
        D[j, i] <- 0
      } else if (refined_x[j] <= knots[i + 1] & refined_x[j] >= knots[i]) {
        D[j, i] <- (1 / factorial(p)) * (refined_x[j] - knots[i])^p
      } else {
        k <- 1:p
        D[j, i] <- sum((dif[i]^k) * ((refined_x[j] - knots[i + 1])^(p - k)) / (factorial(k) * factorial(p - k)))
      }
    }
  }
  D # Local poly design matrix
}

#' Constructing and evaluating the local O-spline basis (design matrix)
#'
#' @param knots A vector of knots used to construct the O-spline basis, first knot should be viewed as "0",
#' the reference starting location. These k knots will define (k-1) basis function in total.
#' @param refined_x A vector of locations to evaluate the O-spline basis
#' @param p An integer value indicates the order of smoothness
#' @param neg_sign_order An integer value N such that D = ((-1)^N)*D for the splines at negative knots. Default is 0.
#' @return A matrix with i,j component being the value of jth basis function
#' value at ith element of refined_x, the ncol should equal to number of knots minus 1, and nrow
#' should equal to the number of elements in refined_x.
#' @examples
#' local_poly_helper(knots = c(0, 0.2, 0.4, 0.6, 0.8), refined_x = seq(0, 0.8, by = 0.1), p = 2)
#' @export
local_poly_helper <- function(knots, refined_x, p = 2, neg_sign_order = 0) {
  if (min(knots) >= 0) {
    # The case of all-positive knots
    D <- get_local_poly(knots, refined_x, p)
  } else if (max(knots) <= 0) {
    # Handle the negative part only
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    D <- get_local_poly(knots_neg, refined_x_neg, p)
    D <- D * ((-1)^neg_sign_order)
  } else {
    # Handle the negative part
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    D1 <- get_local_poly(knots_neg, refined_x_neg, p)
    D1 <- D1 * ((-1)^neg_sign_order)
    
    # Handle the positive part
    refined_x_pos <- ifelse(refined_x > 0, refined_x, 0)
    knots_pos <- unique(sort(ifelse(knots > 0, knots, 0)))
    D2 <- get_local_poly(knots_pos, refined_x_pos, p)
    
    D <- cbind(D1, D2)
  }
  D # Local poly design matrix
}

#' Constructing and evaluating the global polynomials, to account for boundary conditions (design matrix)
#'
#' @param x A vector of locations to evaluate the global polynomials
#' @param p An integer value indicates the order of smoothness
#' @return A matrix with i,j componet being the value of jth basis function
#' value at ith element of x, the ncol should equal to p, and nrow
#' should equal to the number of elements in x
#' @examples
#' global_poly_helper(x = c(0, 0.2, 0.4, 0.6, 0.8), p = 2)
#' @export
global_poly_helper <- function(x, p = 2) {
  result <- NULL
  for (i in 1:p) {
    result <- cbind(result, x^(i - 1))
  }
  result
}

#' Constructing and evaluating the global polynomials, to account for boundary conditions (design matrix) of sgp
#'
#' @param refined_x A vector of locations to evaluate the sB basis
#' @param a The frequency of sgp.
#' @param m The number of harmonics to consider
#' @param initial_location Optional reference location for centering the seasonal basis.
#' @return A matrix with i,j componet being the value of jth basis function
#' value at ith element of x, the ncol should equal to (2*m), and nrow
#' should equal to the number of elements in x
#' @export
global_poly_helper_sgp <- function(refined_x, a, m, initial_location = NULL) {
  if(is.null(initial_location)){
    initial_location <- min(refined_x)
  }
  refined_x <- refined_x - initial_location
  X <- NULL
  for (i in 1:m) {
    X <- cbind(X, cos(i*a*refined_x),sin(i*a*refined_x))
  }
  X 
}

#' Construct prior based on d-step prediction SD (for iwp)
#'
#' @param prior A list that contains alpha and u. This specifies the target prior on the d-step SD \eqn{\sigma(d)}, such that \eqn{P(\sigma(d) > u) = alpha}.
#' @param d A numeric value for the prediction step.
#' @param p An integer for the order of iwp.
#' @return A list that contains alpha and u. The prior for the smoothness parameter \eqn{\sigma} such that \eqn{P(\sigma > u) = alpha}, that yields the ideal prior on the d-step SD.
#' @export
prior_conversion_iwp <- function(d, prior, p) {
  correction_factor <- PSD_compute_iwp(h = d, p = p, sd = 1)
  prior_q <- list(alpha = prior$alpha, u = (prior$u / correction_factor))
  prior_q
}


#' Compute the SD correction factor for sgp
#' @param d A numeric value for the prediction step.
#' @param a The frequency parameter of the sgp.
#' @return The correction factor c that should be used to compute the d-step PSD as c*SD.
compute_d_step_sgpsd <- function(d,a){
  sqrt((1/(a^2))*((d/2) - (sin(2*a*d)/(4*a))))
}

PSD_compute_iwp <- function(h, p, sd = 1) {
  Cp <- (h^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2))
  abs(sd) * sqrt(Cp)
}

PSD_compute_sgp <- function(h, a, m = 1, sd = 1) {
  correction_factor <- 0
  for (i in 1:m) {
    correction_factor <- correction_factor + compute_d_step_sgpsd(d = h, a = (i * a))
  }
  abs(sd) * correction_factor
}


#' Construct prior based on d-step prediction SD (for sgp)
#'
#' @param prior A list that contains alpha and u. This specifies the target prior on the d-step SD \eqn{\sigma(d)}, such that \eqn{P(\sigma(d) > u) = alpha}.
#' @param d A numeric value for the prediction step.
#' @param a The frequency parameter of the sgp.
#' @param m The number of harmonics that should be considered, by default m = 1 represents only the sgp.
#' @return A list that contains alpha and u. The prior for the smoothness parameter \eqn{\sigma} such that \eqn{P(\sigma > u) = alpha}, that yields the ideal prior on the d-step SD.
#' @export
prior_conversion_sgp <- function(d, prior, a, m = 1) {
  correction_factor <- PSD_compute_sgp(h = d, a = a, m = m, sd = 1)
  prior_SD <- list(u = prior$u/correction_factor, alpha = prior$alpha)
  prior_SD
}

#' @import Matrix
#' @importClassesFrom Matrix dgCMatrix
dgTMatrix_wrapper <- function(matrix) {
  # result <- as(as(as(matrix, "dMatrix"), "generalMatrix"), "TsparseMatrix")
  result <- as(matrix, "dgCMatrix")
  result
}

#' Construct default MCMC options
#'
#' @param option_list Optional named list overriding default MCMC options.
#' @return A named list of MCMC options.
#' @export
get_default_option_list_MCMC <- function(option_list = list()){
  default_options <- list(chains = 1, cores = 1, init = "random", seed = 123, warmup = 10000, silent = TRUE, laplace = FALSE)
  required_names <- names(default_options)
  for (name in required_names) {
    if(! name %in% names(option_list)){
      option_list[[name]] <- default_options[[name]]
    }
  }
  option_list
}


#' Custom Template Function
#'
#' This function allows for the dynamic modification of a C++ template
#' within the BayesGP package. Users can specify custom content for the 
#' log-likelihood as well as the log-prior of the variance parameter in the template.
#' The copied template matches the package's compiled template, including fixed
#' Gaussian observation-SD handling and all Gaussian prior normalizing constants
#' used in log marginal likelihood calculations.
#' 
#' 
#' @param SETUP A character string or vector containing the 
#'   lines of C++ code to be inserted before the computation of log-likelihood.
#' @param LOG_LIKELIHOOD A character string or vector containing the 
#'   lines of C++ code to be inserted in the log-likelihood section of 
#'   the template. Should be NULL if no changes are to be made to this section.
#' @param LOG_PRIOR A character string or vector containing the 
#'   lines of C++ code to be inserted in the log-prior (of the variance parameter) section of 
#'   the template. Should be NULL if no changes are to be made to this section.
#' @param compile_template A indicator of whether the new template should be compiled. default is FALSE.
#' @return A string representing the path to the temporary .so (or .dll) file
#'   containing the compiled custom C++ code.
#'
#' @export
custom_template <- function(SETUP = NULL, LOG_LIKELIHOOD = NULL, LOG_PRIOR = NULL, compile_template = FALSE){
  template_path <- system.file("extsrc", "BayesGP.cpp", package = "BayesGP")
  file_content <- readLines(template_path)
  if(!is.null(SETUP)){
    start_line = which(file_content == "  // START OF YOUR SETUP")
    end_line = which(file_content == "  // END OF YOUR SETUP")
    if (length(start_line) == 1 && length(end_line) == 1 && start_line < end_line) {
      # Replace the content
      file_content <- c(file_content[1:(start_line)], SETUP, file_content[(end_line+1):length(file_content)])
    } else {
      warning("Markers not found or in incorrect order")
    }
  }
  if(!is.null(LOG_LIKELIHOOD)){
    start_line = which(file_content == "    // START OF YOUR LOG_LIKELIHOOD: ll")
    end_line = which(file_content == "    // END OF YOUR LOG_LIKELIHOOD")
    if (length(start_line) == 1 && length(end_line) == 1 && start_line < end_line) {
      # Replace the content
      file_content <- c(file_content[1:(start_line)], LOG_LIKELIHOOD, file_content[(end_line+1):length(file_content)])
    } else {
      warning("Markers not found or in incorrect order")
    }
  }
  if(!is.null(LOG_PRIOR)){
    start_line = which(file_content == "  // START OF YOUR LOG_PRIOR: lpT")
    end_line = which(file_content == "  // END OF YOUR LOG_PRIOR")
    if (length(start_line) == 1 && length(end_line) == 1 && start_line < end_line) {
      # Replace the content
      file_content <- c(file_content[1:(start_line)], LOG_PRIOR, file_content[(end_line+1):length(file_content)])
    } else {
      warning("Markers not found or in incorrect order")
    }
  }
  
  # Create a temporary file
  temp_file <- tempfile(fileext = ".cpp") # .cpp extension is added for clarity
  
  # Write the modified content to the temporary file
  writeLines(file_content, temp_file)
  temp_so_file <- sub("\\.cpp$", "", temp_file)
  
  if(compile_template == TRUE){
    suppressWarnings(suppressMessages(invisible(TMB::compile(file = temp_file)) ))
    dyn.load(TMB::dynlib(temp_so_file))
  }
  
  return(temp_so_file)
}


#' Roxygen commands
#'
#' @useDynLib BayesGP
#'
dummy <- function() {
  return(NULL)
}
