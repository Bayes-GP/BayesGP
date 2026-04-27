## let a(.) be a given function
Compute_Prec <- function(a = NULL, alpha = NULL, c, k, region, accuracy = 0.01, boundary = TRUE){
  a <- resolve_curvature_parameter(a = a, alpha = alpha)
  ss <- function(M) {Matrix::forceSymmetric(M + Matrix::t(M))}
  x <- fem_integration_grid(region = region, accuracy = accuracy, c = c)
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
  Bmatrix <- as.matrix(fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE))
  B1matrix <- as.matrix(fda::eval.basis(x, B_basis, Lfdobj=1, returnMatrix=TRUE))
  B2matrix <- as.matrix(fda::eval.basis(x, B_basis, Lfdobj=2, returnMatrix=TRUE))
  a_func <- function(x) {-1/(a*(x+c))}
  a_matrix <- a_func(x)
  B1a <- Matrix::Matrix(apply(B1matrix, 2, function(x) x * a_matrix), sparse = TRUE)
  B2matrix <- Matrix::Matrix(B2matrix, sparse = TRUE)

  ### Compute I, L, T:
  Numerical_I <- Matrix::Diagonal(x = c(diff(c(0, x))))

  G <- Matrix::t(B1a) %*% Numerical_I %*% B1a
  C <- Matrix::t(B2matrix) %*% Numerical_I %*% B2matrix
  M <- ss(Matrix::t(B1a) %*% Numerical_I %*% B2matrix)
  Q <- G - M + C
  Matrix::forceSymmetric(Q)
}

fem_integration_grid <- function(region, accuracy, c){
  region <- range(region)
  lower <- min(region)
  upper <- max(region)

  if(isTRUE(all.equal(c, 0)) && isTRUE(all.equal(lower, 0))){
    x <- seq(lower + accuracy, upper, by = accuracy)
  } else {
    x <- seq(lower, upper, by = accuracy)
  }
  if(length(x) == 0 || tail(x, 1) < upper){
    x <- c(x, upper)
  }
  x
}

Compute_Design <- function(x, k, region, boundary = TRUE){
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
  as.matrix(fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE))
}

Compute_Prec_rev <- function(a = NULL, alpha = NULL, c, k, region, accuracy = 0.01, boundary = TRUE){
  a <- resolve_curvature_parameter(a = a, alpha = alpha)
  ss <- function(M) {Matrix::forceSymmetric(M + Matrix::t(M))}
  x <- seq(min(region),max(region),by = accuracy)
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
  Bmatrix <- as.matrix(fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE))
  B1matrix <- as.matrix(fda::eval.basis(x, B_basis, Lfdobj=1, returnMatrix=TRUE))
  B2matrix <- as.matrix(fda::eval.basis(x, B_basis, Lfdobj=2, returnMatrix=TRUE))
  a_func <- function(x) {1/(a*(c-x))}
  a_matrix <- a_func(x)
  B1a <- Matrix::Matrix(apply(B1matrix, 2, function(x) x * a_matrix), sparse = TRUE)
  B2matrix <- Matrix::Matrix(B2matrix, sparse = TRUE)

  ### Compute I, L, T:
  Numerical_I <- Matrix::Diagonal(x = c(diff(c(0, x))))

  G <- Matrix::t(B1a) %*% Numerical_I %*% B1a
  C <- Matrix::t(B2matrix) %*% Numerical_I %*% B2matrix
  M <- ss(Matrix::t(B1a) %*% Numerical_I %*% B2matrix)
  Q <- G - M + C
  Matrix::forceSymmetric(Q)
}
