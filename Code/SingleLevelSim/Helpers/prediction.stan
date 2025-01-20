//
// This Stan program estimates the subject-level PC scores
// Used for dynamic prediction


## This script is the code for out-of-sample prediction of 1 subject

// The input data 
data {
  int<lower=0> J; //number of observations along the entire track
  int<lower=0> Ju; //number of observations avaialble
  
  array[Ju] int<lower=0,upper=1> Y; //observed binary outcome
  
  int<lower=0> K; //number of  eigenfunctions
  
  matrix[J, K] efuncs; //leigenfunctions
  vector[J] b0; //mean function
  vector<lower=0>[K] lambda; //eigenvalues
}


// The parameters to estiamte
// individual-level random slopes
parameters {
  matrix[K, 1] xi; 
}


// The model 
model {
  vector[Ju] eta; // linear predictors
  
  // prior
  to_vector(xi) ~ normal(0, lambda);

  // update linear predictor
  eta = b0[1:Ju] + to_vector(efuncs[1:Ju, ] * xi);
  
  // posterior
  Y ~ bernoulli_logit(eta);
}


