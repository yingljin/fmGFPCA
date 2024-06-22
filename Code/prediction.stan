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

## This script is the code for out-of-sample prediction of 1 subject

// The input data 
data {
  int<lower=0> J; //number of total visits
  int<lower=0> K; //number of observations along the entire track
  
  int<lower=0> Ju; //number of visits for subject u
  int<lower=0> Ku[Ju]; //number of observations per visit
  int<lower=0> Nobs; // number of total observations for subject u
  array[Nobs] int<lower=0,upper=1> Y; //observed binary outcome
  
  int<lower=0> L; //number of level 1 eigenfunctions
  int<lower=0> M; //number of level 2 eigenfunctions
  
  matrix[K, L] efuncs_l1; //level 1 eigenfunctions
  matrix[K, M] efuncs_l2; //level 2 eigenfunctions
  vector[K] b0; //mean function
  vector<lower=0>[L] lambda; //level1 eigenvalues
  vector<lower=0>[M] gamma; //level2 eigenvalues
}


// The parameters to estiamte
// individual-level random slopes
parameters {
  matrix[L, 1] xi; //level 1 random slopes
  matrix[M, J] zeta; //level 2 random slopes
}


// The model 
model {
  //int kuj; //number of observations at a specific visi
  vector[Nobs] eta; // linear predictors
  //matrix[kuj, 1] eta_l1;
  //matrix[kuj, 1] eta_l2;
  int pre_nobs; //number of observations up until now

  // prior
  to_vector(xi) ~ normal(0, lambda);
  for(j in 1:J){
    zeta[ , j] ~ normal(0, gamma);
  }
  
  // update linear predictor
  pre_nobs = 0;
  for(j in 1:Ju){
    //number of observation in this visit is Ku[j]
    //kuj = Ku[j];
    //eta_l1 = efuncs_l1[1:Ku[j], ] * xi;
    //eta_l2 = efuncs_l2[1:Ku[j], ] * to_matrix(zeta[, j]);
    eta[(pre_nobs+1):(pre_nobs+Ku[j])] = b0[1:Ku[j]] + to_vector(efuncs_l1[1:Ku[j], ] * xi) + to_vector(efuncs_l2[1:Ku[j], ] * to_matrix(zeta[, j]));
    pre_nobs = pre_nobs+Ku[j];
    //eta = 0.5;
  }
  
  //eta = to_vector(eta);
  
  // posterior
  Y ~ bernoulli_logit(eta);
  // target+=reduce_sum(partial_sum_lpmf, Y, grainsize, eta);  
}


