//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {
 real partial_sum_lpmf(array[] int Y, int start, int end, vector eta) {
   return bernoulli_logit_lupmf(Y | eta[start:end]);
   }
}


// Input data
data {
  int<lower=0> N; //sample size
  int<lower=0> J; //number of visits per subject
  int<lower=0> K; //number of observations per visit
  int<lower=0> L; //number of level 1 eigenfunctions
  int<lower=0> M; //number of level 2 eigenfunctions
  int<lower=0> Nobs; // number of total observations N*J*K
  array[Nobs] int<lower=0,upper=1> Y; //observed binary outcome
  matrix[L, K] efuncs_l1; //level 1 eigenfunctions
  matrix[M, K] efuncs_l2; //level 2 eigenfunctions
  int<lower=1> grainsize;
}

// The parameters to be re-estimated
parameters {
  real b0; //mean function
  vector<lower=0>[L] lambda; //level1 eigenvalues
  vector<lower=0>[M] gamma; //level2 eigenvalues
  matrix[N, L] xi; //level 1 random slopes
  matrix[N*J, M] zeta; //level 2 random slopes
}

// Some intermediate values
transformed parameters{
  matrix[N, L] xi_sc; // level 1 random effects
  matrix[N*J, M] zeta_sc; //level 2 random effects
  for(l in 1:L){
    xi_sc[, l] = xi[,l]*sqrt(lambda[l]);
  }
  for(m in 1:M){
    zeta_sc[, m] = zeta[,m]*sqrt(gamma[m]);
  }
}



// model
model {
  int inx; // observation index
  int inx_ij;
  vector[Nobs] eta; // linear predictors
  matrix[N, K] eta_l1;
  matrix[N*J, K] eta_l2; 
  
  // prior
  b0 ~ normal(0, 10);
  lambda ~ inv_gamma(1, 1);
  gamma ~ inv_gamma(1, 1);
  to_vector(xi) ~ normal(0, 1);
  to_vector(zeta) ~ normal(0, 1);
  
  // update linear predictor
  eta_l1 = xi_sc * efuncs_l1;
  eta_l2 = zeta_sc * efuncs_l2;
  inx = 1;
  for(i in 1:N){
    for(j in 1:J){
      eta_l2[inx, ] = eta_l2[inx,]+eta_l1[i, ];
      inx = inx + 1;
    }
  }
  
  eta = to_vector((eta_l2)');
  
  target+=reduce_sum(partial_sum_lpmf, Y, grainsize, eta+b0);  
}

