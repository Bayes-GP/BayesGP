#define TMB_LIB_INIT R_init_BayesGP
#include <TMB.hpp>
#include <string>
using namespace tmbutils;

template<class Type>
struct list_SparseMatrix_from_R : vector<SparseMatrix<Type> > {
  list_SparseMatrix_from_R(SEXP x){
    int len_x = LENGTH(x);
    (*this).resize(len_x);
    for (int index = 0; index < len_x; index++){
      SEXP cur = VECTOR_ELT(x, index);
      (*this)(index) = asSparseMatrix<Type>(cur);
    }
  }
};

template<class Type>
struct list_Scalar_from_R : vector<Type> {
  list_Scalar_from_R(SEXP x){
    int len_x = LENGTH(x);
    (*this).resize(len_x);
    for (int index = 0; index < len_x; index++){
      SEXP cur = VECTOR_ELT(x, index);
      (*this)(index) = *REAL(cur);
    }
  }
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  // For random effects
  DATA_STRUCT(X, list_SparseMatrix_from_R); // Design matrix
  DATA_STRUCT(B, list_SparseMatrix_from_R); // Design matrix
  DATA_STRUCT(P, list_SparseMatrix_from_R); // Penalty matrix
  DATA_STRUCT(logPdet, list_Scalar_from_R); // Determinant of (fixed) penalty matrix

  // The last element of u and the last element of alpha are for the family
  // (For example: the Gaussian noise variance when family = "Gaussian")
  DATA_STRUCT(u, list_Scalar_from_R);       // Pc prior, u param
  DATA_STRUCT(alpha, list_Scalar_from_R);   // Pc prior, alpha param
  DATA_STRUCT(betaprec, list_Scalar_from_R);// For boundary, beta ~iid N(betamean,1/betaprec)
  DATA_STRUCT(betamean, list_Scalar_from_R);// For boundary, beta ~iid N(betamean,1/betaprec)

  // For fixed effects
  DATA_STRUCT(beta_fixed_prec, list_Scalar_from_R);
  DATA_STRUCT(beta_fixed_mean, list_Scalar_from_R);
  DATA_STRUCT(Xf, list_SparseMatrix_from_R);
  DATA_VECTOR(offset_sum);
  DATA_VECTOR(y); //response variable
  DATA_SCALAR(family_type); // Family types: Gaussian - 0, Poisson - 1, Binomial - 2, Coxph - 3, Case-Crossover - 4
  DATA_SCALAR(gaussian_sd_known);
  DATA_SCALAR(gaussian_sd_value);

  vector<int> betadim(X.size());
  int sum_betadim = 0;
  for (int i = 0; i < X.size(); i++){
    betadim(i) = X(i).cols(); // Number of boundary conditions in RE #i
    sum_betadim += betadim(i);
  }

  vector<int> d(P.size());
  int sum_d = 0;
  for (int i = 0; i < P.size(); i++){
    d(i) = P(i).cols(); // Number of spline in RE #i
    sum_d += d(i);
  }

  vector<int> beta_fixed_dim(Xf.size());
  for (int i = 0; i < Xf.size(); i++){
    // beta_fixed_dim(0) is for the intercept 
    // beta_fixed_dim(1) is for the fixed effect 1, etc
    beta_fixed_dim(i) = Xf(i).cols(); 
  }

  // Parameter
  PARAMETER_VECTOR(W); // W = c(U1,U2, beta1, beta2, beta_fixed0, beta_fixed1, beta_fixed2)

  vector<vector<Type>> U(d.size());
  for (int i = 0; i < d.size(); i++){
    vector<Type> cur_U(d(i));
    U(i) = cur_U;
  }

  vector<vector<Type>> beta(betadim.size());
  for (int i = 0; i < betadim.size(); i++){
    vector<Type> cur_beta(betadim(i));
    beta(i) = cur_beta;
  }

  vector<vector<Type>> beta_fixed(beta_fixed_dim.size());
  for (int i = 0; i < beta_fixed_dim.size(); i++){
    vector<Type> cur_beta_fixed(beta_fixed_dim(i));
    beta_fixed(i) = cur_beta_fixed;
  }

  int cur_dim_sum  = 0;
  for (int j = 0; j < d.size(); j++){
    if (j == 0){
      for (int i=0;i<d(j);i++) U(j)(i) = W(i);
    }
    else{
      for (int i=0;i<d(j);i++) U(j)(i) = W(i + cur_dim_sum);
    }
    cur_dim_sum += d(j);
  }

  for (int j=0;j<betadim.size();j++){
    if (j == 0){
      for (int i=0;i<betadim(j);i++) beta(j)(i) = W(i + sum_d);
    }
    else{
      int prev_betadim_sum = 0;
      for (int i=0;i<j;i++) prev_betadim_sum += betadim(i);
      for (int i=0;i<betadim(j);i++) beta(j)(i) = W(i + sum_d + prev_betadim_sum);
    }
  }

  for (int j=0;j<beta_fixed_dim.size();j++){
    if (j == 0){
      for (int i=0;i<beta_fixed_dim(j);i++) beta_fixed(j)(i) = W(i + sum_d + sum_betadim);
    }
    else{
      int prev_beta_fixed_dim_sum = 0;
      for (int i=0;i<j;i++) prev_beta_fixed_dim_sum += beta_fixed_dim(i);
      for (int i=0;i<beta_fixed_dim(j);i++) beta_fixed(j)(i) = W(i + sum_d + sum_betadim + prev_beta_fixed_dim_sum);
    }
  }

  // The last element is for the variance of Gaussian family 
  PARAMETER_VECTOR(theta);

  // Transformations
  vector<Type> eta(y.size());
  if (y.size() == offset_sum.size()){
    eta += offset_sum;
  }
  for (int i = 0; i < Xf.size(); i++){
    eta += Xf(i) * beta_fixed(i); // adding each fixed effect
  }

  for (int i = 0; i < X.size(); i++){
    eta += X(i) * beta(i);  // adding each boundary condition for RE
  }

  for (int i = 0; i < B.size(); i++){
    eta += B(i) * U(i); // adding each local spline effect for RE
  }


  vector<Type> sigma(theta.size());
  for (int i = 0; i < theta.size(); i++){
    sigma(i) = exp(-0.5*theta(i)); // transforming each variance parameter (from log precision to standard deviation)
    REPORT(sigma(i));
  }

  // START OF YOUR SETUP
  
  // END OF YOUR SETUP


  // Log likelihood
  Type ll = 0;
  // Family types: Gaussian - 0, Poisson - 1, Binomial - 2,
  // CoxPH - 3, Case-crossover - 4
  if (family_type == 0){
    Type obs_sd = gaussian_sd_known > Type(0.5) ? Type(gaussian_sd_value) : sigma(theta.size() - 1);
    ll = sum(dnorm(y, eta, obs_sd, TRUE));
  } 
  else if (family_type == 1){
    ll = sum(dpois(y, exp(eta), TRUE));
  } 
  else if (family_type == 2){
    DATA_VECTOR(size); // Size vector for binomial data
    ll = sum(dbinom_robust(y, size, eta, TRUE));
  } 
     // family = 'sum(dbinom_robust(y, size, eta, TRUE)'
  else if (family_type == 3){
    DATA_IVECTOR(cens); // cens vector for coxph data
    DATA_IVECTOR(ranks);
    DATA_SPARSE_MATRIX(D); // Differencing matrix to compute delta;
    // Log likelihood
    int n = y.size(); // Sample size
    vector<Type> delta_red = D * eta;
    vector<Type> delta(n);
    delta(0) = 0;
    for (int i=1;i<n;i++) {
      delta(i) = delta_red(i-1);
      }
    for (int i=0;i<n;i++){
      int nn = n-ranks(i)+1;
      vector<Type> delta_vec_i(nn); //rank starts at 1!!!
      for(int j=0;j<nn;j++) {
        delta_vec_i(j) = delta(n - nn + j);
      }
      vector<Type> diffvec(nn);
      for (int j=0;j<nn;j++) {
        diffvec(j) = exp(delta(i) - delta_vec_i(j));
      }
      ll += -cens(i) * log(diffvec.sum());
  }
  } 

  else if (family_type == 4){
    DATA_IMATRIX(stratum_members);
    DATA_VECTOR(stratum_count);
    int n_strata = stratum_members.rows();
    int n_members = stratum_members.cols();

    for (int s = 0; s < n_strata; s++) {
      Type log_denom = 0;
      Type weighted_eta_sum = 0;
      bool has_member = false;

      for (int j = 0; j < n_members; j++) {
        int row_index = stratum_members(s, j);
        if (row_index == 0) continue;

        Type row_eta = eta(row_index - 1);
        weighted_eta_sum += y(row_index - 1) * row_eta;

        if (!has_member) {
          log_denom = row_eta;
          has_member = true;
        } else {
          log_denom = logspace_add(log_denom, row_eta);
        }
      }

      if (has_member) {
        ll += weighted_eta_sum - stratum_count(s) * log_denom;
      }
    }
}

  else if (family_type == -2){
    ll = 0;
  }

  else if (family_type == -1){
    // START OF YOUR LOG_LIKELIHOOD: ll
    ll = 0;
    // END OF YOUR LOG_LIKELIHOOD
  }

  REPORT(ll);

  // Log prior on W
  Type lpW = 0;
  Type log_two_pi = log(Type(2.0 * M_PI));
  // Cross product (for each RE and its boundary, and for fixed effect)
  // For Random Effects:
  int cur_dim_sum_boundary = 0;
  for (int i = 0; i < X.size(); i++){
    for (int j = 0; j < X(i).cols(); j++){
      Type bb = (beta(i)(j) - betamean(cur_dim_sum_boundary + j)) * (beta(i)(j) - betamean(cur_dim_sum_boundary + j));
      lpW += 0.5 * log(betaprec(cur_dim_sum_boundary + j)) - 0.5 * log_two_pi;
      lpW += -0.5 * betaprec(cur_dim_sum_boundary + j) * bb; // Beta part (boundary condition)
    }
    cur_dim_sum_boundary += X(i).cols();
  }

  for (int i = 0; i < P.size(); i++){
    lpW += -0.5 * exp(theta(i)) * ((U(i) * (P(i) * U(i))).sum()); // U part (spline effect)
    // Log determinant
    Type logdet = d(i) * theta(i) + logPdet(i);
    lpW += -0.5 * Type(d(i)) * log_two_pi;
    lpW += 0.5 * logdet; // P part
  }

  // For Fixed Effects;
  for (int i = 0; i < beta_fixed.size(); i++){
    Type bbf = ((beta_fixed(i) - beta_fixed_mean(i)) * (beta_fixed(i) - beta_fixed_mean(i))).sum(); // Fixed effect
    lpW += Type(beta_fixed(i).size()) * (0.5 * log(beta_fixed_prec(i)) - 0.5 * log_two_pi);
    lpW += -0.5 * beta_fixed_prec(i) * bbf; //
  }

  // START OF YOUR LOG_PRIOR: lpT
  // Log prior for theta
  Type lpT = 0;
  // Variance of each random effect (and the family)
  for (int i = 0; i < alpha.size(); i++){
    Type phi = -log(alpha(i)) / u(i);
    lpT += log(0.5 * phi) - phi*exp(-0.5*theta(i)) - 0.5*theta(i);
  }
  // END OF YOUR LOG_PRIOR

  // Final result!
  Type logpost = -1 * (ll + lpW + lpT);

  
  return logpost;
}
