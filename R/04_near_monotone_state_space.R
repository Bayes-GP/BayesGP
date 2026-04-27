# case1 : alpha being positive, != 1 or 0. s > t = x:
mspline_cov_c1 <- function(s,t,c = 1, alpha = 2){
  1/3*(3*(alpha^4 + 2*alpha^3)*(c + s)^(1/alpha)*c^((3*alpha + 1)/alpha) - 3*(2*alpha^4 + alpha^3)*c^(2*(alpha + 1)/alpha)*s + 3*((alpha^4 - alpha^3)*c^3 + (alpha^4 - alpha^3)*c^2*s + ((alpha^4 - alpha^3)*c + (alpha^4 - alpha^3)*s)*t^2 + 2*((alpha^4 - alpha^3)*c^2 + (alpha^4 - alpha^3)*c*s)*t)*(c + t)^(2/alpha) - (3*(alpha^4 + 2*alpha^3)*(c + s)^(1/alpha)*c^3 + 3*(alpha^4 + alpha^3 - 2*alpha^2)*(c + s)^(1/alpha)*c^2*t + 3*(alpha^4 + alpha^3 - 2*alpha^2)*(c + s)^(1/alpha)*c*t^2 + (alpha^4 + alpha^3 - 2*alpha^2)*(c + s)^(1/alpha)*t^3 - 3*(alpha^4 + 2*alpha^3)*c^((2*alpha + 1)/alpha)*s - 3*(alpha^4 + 2*alpha^3)*c^((3*alpha + 1)/alpha))*(c + t)^(1/alpha) - 3*(2*alpha^4 + alpha^3)*c^((3*alpha + 2)/alpha) + 3*((alpha^4 + 2*alpha^3)*(c + s)^(1/alpha)*c^((2*alpha + 1)/alpha) - (2*alpha^4 + alpha^3)*c^((alpha + 2)/alpha)*s - (2*alpha^4 + alpha^3)*c^(2*(alpha + 1)/alpha))*t)/((2*alpha^4 + alpha^3 - 6*alpha^2 + alpha + 2)*(c + s)^(1/alpha)*(c + t)^(1/alpha))
}
mspline_deriv_cov_c1 <- function(s,t,c = 1, alpha = 2){
  ((alpha*c + alpha*t)*(c + t)^(2/alpha) - alpha*c^((alpha + 2)/alpha))/((alpha + 2)*(c + s)^(1/alpha)*(c + t)^(1/alpha))
}
mspline_cross_cov_c1 <- function(s,t,c = 1, alpha = 2){
  ((alpha^3 + 2*alpha^2)*(c + s)^(1/alpha)*c^((2*alpha + 1)/alpha) - (2*alpha^3 + alpha^2)*c^((alpha + 2)/alpha)*s + ((2*alpha^3 + alpha^2)*c^2 + (2*alpha^3 + alpha^2)*c*s + ((2*alpha^3 + alpha^2)*c + (2*alpha^3 + alpha^2)*s)*t)*(c + t)^(2/alpha) - ((alpha^3 + 2*alpha^2)*(c + s)^(1/alpha)*c^2 + 2*(alpha^3 + 2*alpha^2)*(c + s)^(1/alpha)*c*t + (alpha^3 + 2*alpha^2)*(c + s)^(1/alpha)*t^2)*(c + t)^(1/alpha) - (2*alpha^3 + alpha^2)*c^(2*(alpha + 1)/alpha))/((2*alpha^3 + 3*alpha^2 - 3*alpha - 2)*(c + s)^(1/alpha)*(c + t)^(1/alpha))
}
SS_cov_c1 <- function(s,t,c = 1, alpha = 2){
  -1/3*(3*(2*alpha^4 + 5*alpha^3 + 2*alpha^2)*(c + s)^(2/alpha)*c^2*t + 3*(2*alpha^4 + 5*alpha^3 + 2*alpha^2)*(c + s)^(2/alpha)*c*t^2 + (2*alpha^4 + 5*alpha^3 + 2*alpha^2)*(c + s)^(2/alpha)*t^3 + (9*alpha^3*c^3 - 6*(alpha^4 - 2*alpha^3 + alpha^2)*c^2*s - 6*(alpha^4 - 2*alpha^3 + alpha^2)*c*s^2 - 2*(alpha^4 - 2*alpha^3 + alpha^2)*s^3)*(c + s)^(2/alpha) + 3*((2*alpha^4 + alpha^3)*c^3 + 2*(2*alpha^4 + alpha^3)*c^2*s + (2*alpha^4 + alpha^3)*c*s^2 + ((2*alpha^4 + alpha^3)*c^2 + 2*(2*alpha^4 + alpha^3)*c*s + (2*alpha^4 + alpha^3)*s^2)*t)*(c + t)^(2/alpha) - 6*(((alpha^4 + 2*alpha^3)*c + (alpha^4 + 2*alpha^3)*s)*(c + s)^(1/alpha)*t^2 + 2*((alpha^4 + 2*alpha^3)*c^2 + (alpha^4 + 2*alpha^3)*c*s)*(c + s)^(1/alpha)*t + ((alpha^4 + 2*alpha^3)*c^3 + (alpha^4 + 2*alpha^3)*c^2*s)*(c + s)^(1/alpha))*(c + t)^(1/alpha))/((2*alpha^4 + alpha^3 - 6*alpha^2 + alpha + 2)*(c + s)^(2/alpha))
}
SS_deriv_cov_c1 <- function(s,t,c = 1, alpha = 2){
  ((alpha*c + alpha*s)*(c + s)^(2/alpha) - (alpha*c + alpha*t)*(c + t)^(2/alpha))/((alpha + 2)*(c + s)^(2/alpha))
}
SS_cross_cov_c1 <- function(s,t,c = 1, alpha = 2){
  (((alpha^3 - alpha^2)*c^2 + 2*(alpha^3 - alpha^2)*c*s + (alpha^3 - alpha^2)*s^2)*(c + s)^(2/alpha) - ((2*alpha^3 + alpha^2)*c^2 + (2*alpha^3 + alpha^2)*c*s + ((2*alpha^3 + alpha^2)*c + (2*alpha^3 + alpha^2)*s)*t)*(c + t)^(2/alpha) + ((alpha^3 + 2*alpha^2)*(c + s)^(1/alpha)*c^2 + 2*(alpha^3 + 2*alpha^2)*(c + s)^(1/alpha)*c*t + (alpha^3 + 2*alpha^2)*(c + s)^(1/alpha)*t^2)*(c + t)^(1/alpha))/((2*alpha^3 + 3*alpha^2 - 3*alpha - 2)*(c + s)^(2/alpha))
}
SS_cov_mat_c1 <- function(s,t,c = 1, alpha = 2){
  M <- matrix(nrow = 2, ncol = 2)
  M[1,1] <- SS_cov_c1(s = s, t = t, c = c, alpha = alpha)
  M[1,2] <- SS_cross_cov_c1(s = s, t = t, c = c, alpha = alpha)
  M[2,2] <- SS_deriv_cov_c1(s = s, t = t, c = c, alpha = alpha)
  Matrix::forceSymmetric(M)
}
R_trans_matrix_c1 <- function(s,x,c = 1, alpha = 2){
  R <- matrix(nrow = 2, ncol = 2)
  R[1,1] <- 1
  R[2,1] <- 0
  R[1,2] <--(alpha*(c + s)^(1/alpha)*c + alpha*(c + s)^(1/alpha)*x - (alpha*c + alpha*s)*(c + x)^(1/alpha))/((alpha - 1)*(c + s)^(1/alpha))
  R[2,2] <- (c + x)^(1/alpha)/(c + s)^(1/alpha)
  R
}



# case 2: alpha = 1. s > t = x:
mspline_cov_c2 <- function(s,t,c = 1){
  1/3*c^3*log(c + s)*log(c) - 1/3*c^3*log(c)^2 + 1/27*t^3*(3*log(c + s) + 2) + 2/9*c^3*log(c) + 1/9*(3*c*log(c + s) + 2*c)*t^2 + 1/9*(3*c^2*log(c + s) + 2*c^2)*t - 1/9*(3*c^3*log(c + s) - 3*c^3*log(c) + 2*c^3 + 3*c^2*t + 3*c*t^2 + t^3)*log(c + t)
}
mspline_deriv_cov_c2 <- function(s,t,c = 1){
  1/3*(3*c^2*t + 3*c*t^2 + t^3)/(c^2 + c*s + (c + s)*t)
}
mspline_cross_cov_c2 <- function(s,t,c = 1){
  1/9*(t^3*(3*log(c + s) + 1) + 3*c^3*log(c) + 3*(3*c*log(c + s) + c)*t^2 + 3*(3*c^2*log(c + s) + c^2)*t - 3*(c^3 + 3*c^2*t + 3*c*t^2 + t^3)*log(c + t))/(c + t)
}
SS_cov_c2 <- function(s,t,c = 1){
  -1/3*c^3*log(c + s)^2 - 1/27*(9*log(c + s)^2 + 6*log(c + s) + 2)*t^3 - 2/9*c^3*log(c + s) + 2/9*c^2*s + 2/9*c*s^2 + 2/27*s^3 - 1/9*(9*c*log(c + s)^2 + 6*c*log(c + s) + 2*c)*t^2 - 1/3*(c^3 + 3*c^2*t + 3*c*t^2 + t^3)*log(c + t)^2 - 1/9*(9*c^2*log(c + s)^2 + 6*c^2*log(c + s) + 2*c^2)*t + 2/9*(t^3*(3*log(c + s) + 1) + 3*c^3*log(c + s) + c^3 + 3*(3*c*log(c + s) + c)*t^2 + 3*(3*c^2*log(c + s) + c^2)*t)*log(c + t)

}
SS_deriv_cov_c2 <- function(s,t,c = 1){
  1/3*(3*c^2*s + 3*c*s^2 + s^3 - 3*c^2*t - 3*c*t^2 - t^3)/(c^2 + 2*c*s + s^2)

}
SS_cross_cov_c2 <- function(s,t,c = 1){
  -1/9*(t^3*(3*log(c + s) + 1) + 3*c^3*log(c + s) - 3*c^2*s - 3*c*s^2 - s^3 + 3*(3*c*log(c + s) + c)*t^2 + 3*(3*c^2*log(c + s) + c^2)*t - 3*(c^3 + 3*c^2*t + 3*c*t^2 + t^3)*log(c + t))/(c + s)

}
SS_cov_mat_c2 <- function(s,t,c = 1){
  M <- matrix(nrow = 2, ncol = 2)
  M[1,1] <- SS_cov_c2(s = s, t = t, c = c)
  M[1,2] <- SS_cross_cov_c2(s = s, t = t, c = c)
  M[2,2] <- SS_deriv_cov_c2(s = s, t = t, c = c)
  Matrix::forceSymmetric(M)
}
R_trans_matrix_c2 <- function(s,x,c = 1){
  R <- matrix(nrow = 2, ncol = 2)
  R[1,1] <- 1
  R[2,1] <- 0
  R[1,2] <- c*log(c + s) + x*log(c + s) - (c + x)*log(c + x)
  R[2,2] <- (c + x)/(c + s)
  R
}

# case 3: alpha < 0 but alpha != -1 or -2. s > t = x:
mspline_cov_c3 <- function(s,t,c = 1, alpha = -3){
  alpha <- abs(alpha)
  1/3*(3*((alpha^4 + alpha^3)*c^(2/alpha)*s + (alpha^4 + alpha^3)*c^((alpha + 2)/alpha))*(c + s)^(1/alpha)*t^2 + 6*((alpha^4 + alpha^3)*c^((alpha + 2)/alpha)*s + (alpha^4 + alpha^3)*c^(2*(alpha + 1)/alpha))*(c + s)^(1/alpha)*t + 3*((alpha^4 + alpha^3)*c^(2*(alpha + 1)/alpha)*s + (alpha^4 + alpha^3)*c^((3*alpha + 2)/alpha))*(c + s)^(1/alpha) - 3*(((2*alpha^4 - alpha^3)*c^3 + (2*alpha^4 - alpha^3)*c^2*s)*(c + s)^(1/alpha) - (alpha^4 - 2*alpha^3)*c^((3*alpha + 1)/alpha) + (((2*alpha^4 - alpha^3)*c^2 + (2*alpha^4 - alpha^3)*c*s)*(c + s)^(1/alpha) - (alpha^4 - 2*alpha^3)*c^((2*alpha + 1)/alpha))*t)*(c + t)^(2/alpha) - ((alpha^4 - alpha^3 - 2*alpha^2)*c^(2/alpha)*t^3 + 3*(alpha^4 - alpha^3 - 2*alpha^2)*c^((alpha + 2)/alpha)*t^2 + 3*(alpha^4 - alpha^3 - 2*alpha^2)*c^(2*(alpha + 1)/alpha)*t - 3*((alpha^4 - 2*alpha^3)*c^((2*alpha + 1)/alpha)*s + (alpha^4 - 2*alpha^3)*c^((3*alpha + 1)/alpha))*(c + s)^(1/alpha) + 3*(alpha^4 - 2*alpha^3)*c^((3*alpha + 2)/alpha))*(c + t)^(1/alpha))/((2*alpha^4 - alpha^3 - 6*alpha^2 - alpha + 2)*(c + t)^(1/alpha)*c^(2/alpha))
}

mspline_deriv_cov_c3 <- function(s,t,c = 1, alpha = -3){
  alpha <- abs(alpha)
  -(alpha*(c + s)^(1/alpha)*(c + t)^(2/alpha)*c - alpha*(c + s)^(1/alpha)*c^(2/alpha)*t - alpha*(c + s)^(1/alpha)*c^((alpha + 2)/alpha))/((alpha - 2)*(c + t)^(1/alpha)*c^(2/alpha))
}
mspline_cross_cov_c3 <- function(s,t,c = 1, alpha = -3){
  alpha <- abs(alpha)
  (((2*alpha^3 - alpha^2)*c^(2/alpha)*s + (2*alpha^3 - alpha^2)*c^((alpha + 2)/alpha))*(c + s)^(1/alpha)*t + ((2*alpha^3 - alpha^2)*c^((alpha + 2)/alpha)*s + (2*alpha^3 - alpha^2)*c^(2*(alpha + 1)/alpha))*(c + s)^(1/alpha) - (((2*alpha^3 - alpha^2)*c^2 + (2*alpha^3 - alpha^2)*c*s)*(c + s)^(1/alpha) - (alpha^3 - 2*alpha^2)*c^((2*alpha + 1)/alpha))*(c + t)^(2/alpha) - ((alpha^3 - 2*alpha^2)*c^(2/alpha)*t^2 + 2*(alpha^3 - 2*alpha^2)*c^((alpha + 2)/alpha)*t + (alpha^3 - 2*alpha^2)*c^(2*(alpha + 1)/alpha))*(c + t)^(1/alpha))/((2*alpha^3 - 3*alpha^2 - 3*alpha + 2)*(c + t)^(1/alpha)*c^(2/alpha))
}
SS_cov_c3 <- function(s,t,c = 1, alpha = -3){
  alpha <- abs(alpha)
  -1/3*(3*((2*alpha^4 - alpha^3)*c^2 + 2*(2*alpha^4 - alpha^3)*c*s + (2*alpha^4 - alpha^3)*s^2)*(c + s)^(2/alpha)*t + 3*((2*alpha^4 - alpha^3)*c^3 + 2*(2*alpha^4 - alpha^3)*c^2*s + (2*alpha^4 - alpha^3)*c*s^2)*(c + s)^(2/alpha) - (9*alpha^3*c^3 + 6*(alpha^4 + 2*alpha^3 + alpha^2)*c^2*s + 6*(alpha^4 + 2*alpha^3 + alpha^2)*c*s^2 + 2*(alpha^4 + 2*alpha^3 + alpha^2)*s^3 - 3*(2*alpha^4 - 5*alpha^3 + 2*alpha^2)*c^2*t - 3*(2*alpha^4 - 5*alpha^3 + 2*alpha^2)*c*t^2 - (2*alpha^4 - 5*alpha^3 + 2*alpha^2)*t^3)*(c + t)^(2/alpha) - 6*(((alpha^4 - 2*alpha^3)*c + (alpha^4 - 2*alpha^3)*s)*(c + s)^(1/alpha)*t^2 + 2*((alpha^4 - 2*alpha^3)*c^2 + (alpha^4 - 2*alpha^3)*c*s)*(c + s)^(1/alpha)*t + ((alpha^4 - 2*alpha^3)*c^3 + (alpha^4 - 2*alpha^3)*c^2*s)*(c + s)^(1/alpha))*(c + t)^(1/alpha))/((2*alpha^4 - alpha^3 - 6*alpha^2 - alpha + 2)*(c + t)^(2/alpha))
}
SS_deriv_cov_c3 <- function(s,t,c = 1, alpha = -3){
  alpha <- abs(alpha)
  -(alpha*(c + s)^(2/alpha)*c + alpha*(c + s)^(2/alpha)*t - (alpha*c + alpha*s)*(c + t)^(2/alpha))/((alpha - 2)*(c + t)^(2/alpha))
}
SS_cross_cov_c3 <- function(s,t,c = 1, alpha = -3){
  alpha <- abs(alpha)
  -(((2*alpha^3 - alpha^2)*c + (2*alpha^3 - alpha^2)*s)*(c + s)^(2/alpha)*t + ((2*alpha^3 - alpha^2)*c^2 + (2*alpha^3 - alpha^2)*c*s)*(c + s)^(2/alpha) - ((alpha^3 + alpha^2)*c^2 + 2*(alpha^3 + alpha^2)*c*s + (alpha^3 + alpha^2)*s^2)*(c + t)^(2/alpha) - ((alpha^3 - 2*alpha^2)*(c + s)^(1/alpha)*c^2 + 2*(alpha^3 - 2*alpha^2)*(c + s)^(1/alpha)*c*t + (alpha^3 - 2*alpha^2)*(c + s)^(1/alpha)*t^2)*(c + t)^(1/alpha))/((2*alpha^3 - 3*alpha^2 - 3*alpha + 2)*(c + t)^(2/alpha))
}
SS_cov_mat_c3 <- function(s,t,c = 1, alpha = -3){
  alpha <- abs(alpha)
  M <- matrix(nrow = 2, ncol = 2)
  M[1,1] <- SS_cov_c3(s = s, t = t, c = c, alpha = alpha)
  M[1,2] <- SS_cross_cov_c3(s = s, t = t, c = c, alpha = alpha)
  M[2,2] <- SS_deriv_cov_c3(s = s, t = t, c = c, alpha = alpha)
  Matrix::forceSymmetric(M)
}
R_trans_matrix_c3 <- function(s,x,c = 1, alpha = -3){
  alpha <- abs(alpha)
  R <- matrix(nrow = 2, ncol = 2)
  R[1,1] <- 1
  R[2,1] <- 0
  R[1,2] <- ((alpha*c + alpha*s)*(c + s)^(1/alpha) - (alpha*c + alpha*x)*(c + x)^(1/alpha))/((alpha + 1)*(c + x)^(1/alpha))
  R[2,2] <- (c + s)^(1/alpha)/(c + x)^(1/alpha)
  R
}

# case 4: alpha = -1
mspline_cov_c4 <- function(s,t,c = 1){
  -1/12*(2*c*t^3 - 3*(2*c*s + s^2)*t^2)/c
}
mspline_deriv_cov_c4 <- function(s,t,c = 1){
  (c + s)*t/c
}
mspline_cross_cov_c4 <- function(s,t,c = 1){
  -1/2*(c*t^2 - (2*c*s + s^2)*t)/c
}
SS_cov_c4 <- function(s,t,c = 1){
  1/12*(4*c*s^3 + 3*s^4 - 4*c*t^3 - t^4 + 6*(2*c*s + s^2)*t^2 - 4*(3*c*s^2 + 2*s^3)*t)/(c + t)
}
SS_deriv_cov_c4 <- function(s,t,c = 1){
  (c*s + s^2 - (c + s)*t)/(c + t)
}
SS_cross_cov_c4 <- function(s,t,c = 1){
  1/2*(c*s^2 + s^3 + (c + s)*t^2 - 2*(c*s + s^2)*t)/(c + t)
}
SS_cov_mat_c4 <- function(s,t,c = 1){
  M <- matrix(nrow = 2, ncol = 2)
  M[1,1] <- SS_cov_c4(s = s, t = t, c = c)
  M[1,2] <- SS_cross_cov_c4(s = s, t = t, c = c)
  M[2,2] <- SS_deriv_cov_c4(s = s, t = t, c = c)
  Matrix::forceSymmetric(M)
}
R_trans_matrix_c4 <- function(s,x,c = 1){
  R <- matrix(nrow = 2, ncol = 2)
  R[1,1] <- 1
  R[2,1] <- 0
  R[1,2] <- 1/2*(2*c*s + s^2 - 2*c*x - x^2)/(c + x)
  R[2,2] <- (c + s)/(c + x)
  R
}

# case 5: alpha = -2 (root # 1)
mspline_cov_c5 <- function(s,t,c = 1){
  -8/27*c^3 - 4/9*c^2*t - 4/9*c*t^2 - 4/27*t^3 + 8/27*(c^2 + c*s)*sqrt(c + s)*sqrt(c) - 4/27*((3*c^2*log(c) + 2*c^2 + (3*c*log(c) + 2*c)*s + (s*(3*log(c) + 2) + 3*c*log(c) + 2*c)*t - 3*(c^2 + c*s + (c + s)*t)*log(c + t))*sqrt(c + s) - 2*(c^2 + c*t)*sqrt(c))*sqrt(c + t)
}
mspline_deriv_cov_c5 <- function(s,t,c = 1){
  sqrt(c + s)*sqrt(c + t)*(log(c + t) - log(c))
}
mspline_cross_cov_c5 <- function(s,t,c = 1){
  2/9*((3*((c + s)*log(c + t) - c*log(c) - s*log(c))*sqrt(c + s)*sqrt(c) + 2*c^2)*sqrt(c + t) - 2*(c^2 + 2*c*t + t^2)*sqrt(c))/sqrt(c)
}
SS_cov_c5 <- function(s,t,c = 1){
  -16/27*c^3 - 4/3*c^2*s - 4/3*c*s^2 - 4/9*s^3 - 4/9*c^2*t - 4/9*c*t^2 - 4/27*t^3 + 16/27*(c^2 + c*s + (c + s)*t)*sqrt(c + s)*sqrt(c + t) + 4/9*(c^3 + 3*c^2*s + 3*c*s^2 + s^3)*log(c + s) - 4/9*(c^3 + 3*c^2*s + 3*c*s^2 + s^3)*log(c + t)
}
SS_deriv_cov_c5 <- function(s,t,c = 1){
  (c + s)*log(c + s) - (c + s)*log(c + t)
}
SS_cross_cov_c5 <- function(s,t,c = 1){
  4/9*sqrt(c + s)*(c + t)^(3/2) - 4/9*c^2 - 8/9*c*s - 4/9*s^2 + 2/3*(c^2 + 2*c*s + s^2)*log(c + s) - 2/3*(c^2 + 2*c*s + s^2)*log(c + t)
}
SS_cov_mat_c5 <- function(s,t,c = 1){
  M <- matrix(nrow = 2, ncol = 2)
  M[1,1] <- SS_cov_c5(s = s, t = t, c = c)
  M[1,2] <- SS_cross_cov_c5(s = s, t = t, c = c)
  M[2,2] <- SS_deriv_cov_c5(s = s, t = t, c = c)
  Matrix::forceSymmetric(M)
}


# case 6: alpha = -1/2 (root # 2)
mspline_cov_c6 <- function(s,t,c = 1){
  1/27*(6*c^6*log(c) + 9*c^5*s*log(c) + 9*c^4*s^2*log(c) + 3*c^3*s^3*log(c) + (3*c^3*log(c) + 2*c^3 + 3*c^2*s + 3*c*s^2 + s^3)*t^3 + 3*(3*c^4*log(c) + 2*c^4 + 3*c^3*s + 3*c^2*s^2 + c*s^3)*t^2 + 3*(3*c^5*log(c) + 2*c^5 + 3*c^4*s + 3*c^3*s^2 + c^2*s^3)*t - 3*(2*c^6 + 3*c^5*s + 3*c^4*s^2 + c^3*s^3 + 3*c^5*t + 3*c^4*t^2 + c^3*t^3)*log(c + t))/c^3
}
mspline_deriv_cov_c6 <- function(s,t,c = 1){
  1/3*((c^2 + 2*c*s + s^2)*t^3 + 3*(c^3 + 2*c^2*s + c*s^2)*t^2 + 3*(c^4 + 2*c^3*s + c^2*s^2)*t)/(c^4 + c^3*t)
}
mspline_cross_cov_c6 <- function(s,t,c = 1){
  1/9*(3*c^6*log(c) + (3*c^3*log(c) + c^3 + 3*c^2*s + 3*c*s^2 + s^3)*t^3 + 3*(3*c^4*log(c) + c^4 + 3*c^3*s + 3*c^2*s^2 + c*s^3)*t^2 + 3*(3*c^5*log(c) + c^5 + 3*c^4*s + 3*c^3*s^2 + c^2*s^3)*t - 3*(c^6 + 3*c^5*t + 3*c^4*t^2 + c^3*t^3)*log(c + t))/(c^4 + c^3*t)
}
SS_cov_c6 <- function(s,t,c = 1){
  1/27*(6*c^5*s + 15*c^4*s^2 + 20*c^3*s^3 + 15*c^2*s^4 + 6*c*s^5 + s^6 - 15*c^2*t^4 - 6*c*t^5 - t^6 - 2*(10*c^3 + 3*(c^3 + 3*c^2*s + 3*c*s^2 + s^3)*log(c + s))*t^3 - 3*(5*c^4 + 6*(c^4 + 3*c^3*s + 3*c^2*s^2 + c*s^3)*log(c + s))*t^2 - 6*(c^5 + 3*(c^5 + 3*c^4*s + 3*c^3*s^2 + c^2*s^3)*log(c + s))*t - 6*(c^6 + 3*c^5*s + 3*c^4*s^2 + c^3*s^3)*log(c + s) + 6*(c^6 + 3*c^5*s + 3*c^4*s^2 + c^3*s^3 + (c^3 + 3*c^2*s + 3*c*s^2 + s^3)*t^3 + 3*(c^4 + 3*c^3*s + 3*c^2*s^2 + c*s^3)*t^2 + 3*(c^5 + 3*c^4*s + 3*c^3*s^2 + c^2*s^3)*t)*log(c + t))/(c^3 + 3*c^2*t + 3*c*t^2 + t^3)
}
SS_deriv_cov_c6 <- function(s,t,c = 1){
  1/3*(3*c^3*s + 6*c^2*s^2 + 4*c*s^3 + s^4 - (c + s)*t^3 - 3*(c^2 + c*s)*t^2 - 3*(c^3 + c^2*s)*t)/(c^3 + 3*c^2*t + 3*c*t^2 + t^3)
}
SS_cross_cov_c6 <- function(s,t,c = 1){
  1/9*(3*c^4*s + 9*c^3*s^2 + 10*c^2*s^3 + 5*c*s^4 + s^5 - (c^2 + 2*c*s + s^2 + 3*(c^2 + 2*c*s + s^2)*log(c + s))*t^3 - 3*(c^3 + 2*c^2*s + c*s^2 + 3*(c^3 + 2*c^2*s + c*s^2)*log(c + s))*t^2 - 3*(c^4 + 2*c^3*s + c^2*s^2 + 3*(c^4 + 2*c^3*s + c^2*s^2)*log(c + s))*t - 3*(c^5 + 2*c^4*s + c^3*s^2)*log(c + s) + 3*(c^5 + 2*c^4*s + c^3*s^2 + (c^2 + 2*c*s + s^2)*t^3 + 3*(c^3 + 2*c^2*s + c*s^2)*t^2 + 3*(c^4 + 2*c^3*s + c^2*s^2)*t)*log(c + t))/(c^3 + 3*c^2*t + 3*c*t^2 + t^3)
}
SS_cov_mat_c6 <- function(s,t,c = 1){
  M <- matrix(nrow = 2, ncol = 2)
  M[1,1] <- SS_cov_c6(s = s, t = t, c = c)
  M[1,2] <- SS_cross_cov_c6(s = s, t = t, c = c)
  M[2,2] <- SS_deriv_cov_c6(s = s, t = t, c = c)
  Matrix::forceSymmetric(M)
}

# Resolve the curvature parameter without relying on partial argument matching.
resolve_curvature_parameter <- function(a = NULL, alpha = NULL){
  if(!is.null(a) && !is.null(alpha) && !isTRUE(all.equal(a, alpha))){
    stop("Please supply only one curvature parameter or keep `a` and `alpha` identical.")
  }

  curvature <- if(!is.null(a)) a else alpha

  if(is.null(curvature)){
    stop("A curvature parameter must be supplied via `a` or `alpha`.")
  }

  curvature
}

# general case: automatically check alpha
mspline_cov <- function(s,t,c = 1, a = NULL, alpha = NULL){
  alpha <- resolve_curvature_parameter(a = a, alpha = alpha)
  if(alpha != 1 & alpha > 0){
    mspline_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    mspline_cov_c2(s = s, t = t, c = c)
  }
  else if(alpha < 0 & alpha != -1 & alpha != -2 & alpha != -(1/2)){
    mspline_cov_c3(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == -1){
    mspline_cov_c4(s = s, t = t, c = c)
  }
  else if(alpha == -2){
    mspline_cov_c5(s = s, t = t, c = c)
  }
  else if(alpha == -1/2){
    mspline_cov_c6(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
mspline_deriv_cov <- function(s,t,c = 1, a = NULL, alpha = NULL){
  alpha <- resolve_curvature_parameter(a = a, alpha = alpha)
  if(alpha != 1 & alpha > 0){
    mspline_deriv_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    mspline_deriv_cov_c2(s = s, t = t, c = c)
  }
  else if(alpha < 0 & alpha != -1 & alpha != -2 & alpha != -(1/2)){
    mspline_deriv_cov_c3(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == -1){
    mspline_deriv_cov_c4(s = s, t = t, c = c)
  }
  else if(alpha == -2){
    mspline_deriv_cov_c5(s = s, t = t, c = c)
  }
  else if(alpha == -1/2){
    mspline_deriv_cov_c6(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
mspline_cross_cov <- function(s,t,c = 1, a = NULL, alpha = NULL){
  alpha <- resolve_curvature_parameter(a = a, alpha = alpha)
  if(alpha != 1 & alpha > 0){
    mspline_cross_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    mspline_cross_cov_c2(s = s, t = t, c = c)
  }
  else if(alpha < 0 & alpha != -1 & alpha != -2 & alpha != -(1/2)){
    mspline_cross_cov_c3(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == -1){
    mspline_cross_cov_c4(s = s, t = t, c = c)
  }
  else if(alpha == -2){
    mspline_cross_cov_c5(s = s, t = t, c = c)
  }
  else if(alpha == -1/2){
    mspline_cross_cov_c6(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
SS_cov <- function(s,t,c = 1, a = NULL, alpha = NULL){
  alpha <- resolve_curvature_parameter(a = a, alpha = alpha)
  if(alpha != 1 & alpha > 0){
    SS_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    SS_cov_c2(s = s, t = t, c = c)
  }
  else if(alpha < 0 & alpha != -1 & alpha != -2 & alpha != -(1/2)){
    SS_cov_c3(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == -1){
    SS_cov_c4(s = s, t = t, c = c)
  }
  else if(alpha == -2){
    SS_cov_c5(s = s, t = t, c = c)
  }
  else if(alpha == -1/2){
    SS_cov_c6(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
SS_deriv_cov <- function(s,t,c = 1, a = NULL, alpha = NULL){
  alpha <- resolve_curvature_parameter(a = a, alpha = alpha)
  if(alpha != 1 & alpha > 0){
    SS_deriv_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    SS_deriv_cov_c2(s = s, t = t, c = c)
  }
  else if(alpha < 0 & alpha != -1 & alpha != -2 & alpha != -(1/2)){
    SS_deriv_cov_c3(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == -1){
    SS_deriv_cov_c4(s = s, t = t, c = c)
  }
  else if(alpha == -2){
    SS_deriv_cov_c5(s = s, t = t, c = c)
  }
  else if(alpha == -1/2){
    SS_deriv_cov_c6(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
SS_cross_cov <- function(s,t,c = 1, a = NULL, alpha = NULL){
  alpha <- resolve_curvature_parameter(a = a, alpha = alpha)
  if(alpha != 1 & alpha > 0){
    SS_cross_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    SS_cross_cov_c2(s = s, t = t, c = c)
  }
  else if(alpha < 0 & alpha != -1 & alpha != -2 & alpha != -(1/2)){
    SS_cross_cov_c3(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == -1){
    SS_cross_cov_c4(s = s, t = t, c = c)
  }
  else if(alpha == -2){
    SS_cross_cov_c5(s = s, t = t, c = c)
  }
  else if(alpha == -1/2){
    SS_cross_cov_c6(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
SS_cov_mat <- function(s,t,c = 1, a = NULL, alpha = NULL){
  alpha <- resolve_curvature_parameter(a = a, alpha = alpha)
  if(alpha != 1 & alpha > 0){
    SS_cov_mat_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    SS_cov_mat_c2(s = s, t = t, c = c)
  }
  else if(alpha < 0 & alpha != -1 & alpha != -2 & alpha != -(1/2)){
    SS_cov_mat_c3(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == -1){
    SS_cov_mat_c4(s = s, t = t, c = c)
  }
  else if(alpha == -2){
    SS_cov_mat_c5(s = s, t = t, c = c)
  }
  else if(alpha == -1/2){
    SS_cov_mat_c6(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
R_trans_matrix <- function(s,x,c = 1, a = NULL, alpha = NULL){
  alpha <- resolve_curvature_parameter(a = a, alpha = alpha)
  if(alpha != 1 & alpha > 0){
    R_trans_matrix_c1(s = s, x = x, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    R_trans_matrix_c2(s = s, x = x, c = c)
  }
  else if(alpha < 0 & alpha != -1){
    R_trans_matrix_c3(s = s, x = x, c = c, alpha = alpha)
  }
  else if(alpha == -1){
    R_trans_matrix_c4(s = s, x = x, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}


# Computing predictive SD for exact mGP:
PSD_compute_mgp <- function(x, h, c = 1, sd = 1, a = NULL, alpha = NULL){
  alpha <- resolve_curvature_parameter(a = a, alpha = alpha)
  abs(sd) * sqrt(SS_cov(s = (x + h), t = (x), c = c, alpha = alpha))
}

# Convert Prior on PSD to SD:
prior_conversion_mgp <- function(prior, d, h = NULL, c = 0.1, a = NULL, alpha = NULL, x = 0){
  alpha <- resolve_curvature_parameter(a = a, alpha = alpha)
  if(!is.null(h)){
    d <- h
  }
  correction_factor <- PSD_compute_mgp(x = x, h = d, c = c, a = alpha, sd = 1)
  prior_q <- list(alpha = prior$alpha, u = (prior$u/correction_factor))
  prior_q
}

prior_conversion_tiwp2 <- function(prior, d, h = NULL, c = 1, a = NULL,
                                   alpha = NULL, x = 0, normalized = TRUE){
  curvature <- resolve_curvature_parameter(a = a, alpha = alpha)
  if(!is.null(h)){
    d <- h
  }
  correction_factor <- PSD_compute_tiwp2(
    x = x,
    h = d,
    c = c,
    a = curvature,
    sd = 1,
    normalized = normalized
  )
  prior_q <- list(alpha = prior$alpha, u = (prior$u / correction_factor))
  prior_q
}


# Computing boundary
boundary_compute <- function(x, c, a = NULL, alpha = NULL){
  alpha <- resolve_curvature_parameter(a = a, alpha = alpha)
  if(alpha == 1){
    cbind(1, log(x+c))
  }
  else{
    cbind(1, (x+c)^((alpha-1)/(alpha)))
  }
}

# Computing joint precision:
mGP_joint_prec <- function(t_vec, a = NULL, alpha = NULL, c){
  alpha <- resolve_curvature_parameter(a = a, alpha = alpha)
  n <- length(t_vec)
  Blist <- list()
  AClist <- list()
  Clist <- list()

  Compute_Ci <- function(t_vec, i){
    Ci <- Matrix::forceSymmetric(solve(SS_cov_mat(t = t_vec[i], s = t_vec[i+1], alpha = alpha, c = c)))
    Ci
  }

  Compute_Ai <- function(t_vec, i, Ci) {
    Ti <- R_trans_matrix(x = t_vec[i], s = t_vec[i+1], alpha = alpha, c = c)
    t(Ti) %*% Ci %*% Ti
  }

  Compute_Bi <- function(t_vec, i, Ci) {
    Ti <- R_trans_matrix(x = t_vec[i], s = t_vec[i+1], alpha = alpha, c = c)
    -t(Ti) %*% Ci
  }

  for (i in 1:(n - 1)) {
    Clist[[i]] <- Compute_Ci(t_vec = t_vec, i = i)
  }

  for (i in 2:(n - 1)) {
    AClist[[i]] <- Compute_Ai(t_vec = t_vec, i = i, Ci = Clist[[i]]) + Clist[[i-1]]
  }

  AClist[[1]] <- Compute_Ai(t_vec = t_vec, i = 1, Ci = Clist[[1]]) + Compute_Ci(t_vec = c(0,t_vec[1]), i = 1)
  AClist[[n]] <- Compute_Ci(t_vec = t_vec, i = (n - 1))

  for (i in 1:(n - 1)) {
    Blist[[i]] <- Compute_Bi(t_vec = t_vec, i = i, Ci = Clist[[i]])
  }

  Qlist <- list()
  Q <- matrix(0, nrow = 0, ncol = n*2)
  for (i in 1:(n-1)) {
    Qlist[[i]] <- cbind(matrix(0,nrow = 2, ncol = 2 * (i-1)), AClist[[i]], Blist[[i]], matrix(0,nrow = 2, ncol = (2 * (n-i-1))) )
    Q <- rbind(Q,Qlist[[i]])
  }
  Q <- rbind(Q,cbind(matrix(0,nrow = 2, ncol = 2 * (n-1)), AClist[[n]]))
  Q <- Matrix::forceSymmetric(Q)
  as(as.matrix(Q), "dgCMatrix")
}

IWP_joint_prec <- function(t_vec, p = 2) {
  Compute_Ti <- function(svec,p = 2,i){
    Ti <- matrix(0,nrow = p, ncol = p)
    delta <- diff(c(0,svec))
    denom <- factorial(c(0:(p-1)))
    numr <- delta[i+1]^(0:(p-1))
    Ti[1,] <- numr/denom
    for (i in 2:p) {
      Ti[i,] <- c(rep(0,(i-1)),Ti[(i-1),((i-1):(p-1))])
    }
    Ti
  }
  Compute_Ci <- function(svec, p = 2, i, is.cov = FALSE){
    delta <- diff(c(0,svec))
    Result <- matrix(0,nrow = p, ncol = p)
    index <- i+1
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
  svec <- t_vec
  n <- length(svec)
  Blist <- list()
  AClist <- list()
  for (i in 1:(n - 1)) {
    AClist[[i]] <- Compute_Ai(svec = svec, i = i , p = p) + Compute_Ci(svec = svec, i = (i-1), p = p)
  }
  AClist[[n]] <- Compute_Ci(svec = svec, i = n - 1, p = p)
  for (i in 1:(n - 1)) {
    Blist[[i]] <- Compute_Bi(svec = svec, i = i, p = p)
  }
  Mlist <- list()
  M <- matrix(0, nrow = 0, ncol = n*p)
  for (i in 1:(n-1)) {
    Mlist[[i]] <- cbind(matrix(0,nrow = p, ncol = p * (i-1)), AClist[[i]], Blist[[i]], matrix(0,nrow = p, ncol = (p * (n-i-1))) )
    M <- rbind(M,Mlist[[i]])
  }
  M <- rbind(M,cbind(matrix(0,nrow = p, ncol = p * (n-1)), AClist[[n]]))
  M <- Matrix::forceSymmetric(M)
  as(as.matrix(M), "dgCMatrix")
}
