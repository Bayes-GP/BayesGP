sampling_from_FEM <- function(x, a, c, k, region, boundary = TRUE, n = 1, accuracy = 0.01){
  Prec <- Compute_Prec(a = a, c = c, k = k, region = region, boundary = boundary, accuracy = accuracy)
  B <- Compute_Design(x, k, region, boundary = boundary)
  if(boundary){
    coefs_samps <- LaplacesDemon::rmvnp(n = n, mu = rep(0,((k-2))), Omega = as.matrix(Prec))
  }
  else{
    coefs_samps <- LaplacesDemon::rmvnp(n = n, mu = rep(0,(k)), Omega = as.matrix(Prec))
  }
  splfd <- B %*% t(coefs_samps)
  (splfd)
}

mGP_sim <- function(t = NULL, mesh_size = 0.01, max_t = 10,
                    a = NULL, alpha = 2, sd = 1, c = 1, initial_vec = NULL){
  curvature <- resolve_curvature_parameter(a = a, alpha = alpha)
  if(is.null(initial_vec)){
    initial_vec <- c(0,0)
  }
  result <- matrix(data = initial_vec, nrow = 1, ncol = 2)
  if (is.null(t)) {
    t <- seq(0, max_t, by = mesh_size)
  }
  for (i in 1:(length(t)-1)) {
    Trans <- R_trans_matrix(s = t[i+1], x = t[i], c = c, a = curvature)
    Sig <- (sd^2) * as.matrix(SS_cov_mat(s = t[i+1], t = t[i], c = c, a = curvature))
    if(sd > 0){
      result_new <- Trans %*% result[i, ] + as.numeric(LaplacesDemon::rmvn(mu = rep(0,2), Sigma = Sig))
    }
    else{
      result_new <- Trans %*% result[i, ]
    }
    result <- rbind(result, matrix(data = result_new, nrow = 1))
  }
  result <- cbind(t, result)
  colnames(result) <- c("t", "func", "func_1st")
  result
}
sim_IWp_Var <- function(t = NULL, mesh_size = 0.01, max_t = 10, p, sd = 1, initial_vec = NULL){
  ### For Precision matrix of the exact method: augmented space
  Compute_Ti <- function(svec,p = 2,i){
    Ti <- matrix(0,nrow = p, ncol = p)
    # delta <- diff(c(0,svec))
    delta <- diff(c(svec))
    denom <- factorial(c(0:(p-1)))
    # numr <- delta[i+1]^(0:(p-1))
    numr <- delta[i]^(0:(p-1))
    Ti[1,] <- numr/denom
    for (i in 2:p) {
      Ti[i,] <- c(rep(0,(i-1)),Ti[(i-1),((i-1):(p-1))])
    }
    Ti
  }
  Compute_Ci <- function(svec, p = 2, i, is.cov = FALSE){
    # delta <- diff(c(0,svec))
    delta <- diff(c(svec))
    Result <- matrix(0,nrow = p, ncol = p)
    # index <- i+1
    index <- i
    for (i in 1:p) {
      for (j in i:p) {
        Result[i,j] <- (delta[index]^(2*p + 1 - i - j))/((2*p + 1 - i - j)*factorial(p-i)*factorial(p-j))
      }
    }
    Result <- Matrix::forceSymmetric(Result)
    if(is.cov == T){
      return(Result)
    }
    else{
      solve(Result)
    }
  }
  Compute_Ai <- function(svec, p = 2, i){
    Ci <- Compute_Ci(svec,p,i)
    Ti <- Compute_Ti(svec,p,i)
    Ai <- t(Ti) %*% Ci
    Ai <- Ai %*% Ti
    Ai
  }
  Compute_Bi <- function(svec, p = 2, i){
    Ci <- Compute_Ci(svec,p,i)
    Ti <- Compute_Ti(svec,p,i)
    Bi <- -t(Ti) %*% Ci
    Bi
  }

  if(is.null(t)){
    t <- seq(0, max_t, by = mesh_size)
  }

  # matT <- Compute_Ti(svec = t[-1], p = p, i = 1)
  # Sig <- (sd^2) * Compute_Ci(svec = t[-1], p = p, i = 1, is.cov = T)

  if(is.null(initial_vec)){
    initial_vec <- c(0,0)
  }

  result <- matrix(data = initial_vec, nrow = 1, ncol = 2)

  for (i in 1:(length(t)-1)) {

    # matT <- Compute_Ti(svec = t[-1], p = p, i = i)
    matT <- Compute_Ti(svec = t, p = p, i = i)
    # Sig <- (sd^2) * Compute_Ci(svec = t[-1], p = p, i = i, is.cov = T)
    Sig <- (sd^2) * as.matrix(Compute_Ci(svec = t, p = p, i = i, is.cov = T))

    if(sd > 0){
      result_new <- matT %*% result[i, ] + as.numeric(LaplacesDemon::rmvn(mu = rep(0,2), Sigma = Sig))
    }
    else{
      result_new <- matT %*% result[i, ]
    }
    result <- rbind(result, matrix(data = result_new, nrow = 1))

  }

  result <- cbind(t, result)
  colnames(result) <- c("t", "func", "func_1st")
  result
}
