# This script writes the function for out-of-sample prediction 
# using outputs from function gm_FPCA

library(tidyverse)
library(lme4)
library(mgcv)
library(refund)
library(here)
library(rstan)

# load("Data/data_N20_J2_K100.RData")
# data <- data_N20J2K100$data


#### Out-of-sample prediction for one subject #### 

# data: dataframe to make prediction on. One subject. Include column id, t, visit
# mu, evals, efuncts: output from gm_FPCA
# J: maximum number of visits in the training data
# K: maximum numbeof of observations per 

mDynPred <- function(data, mu, 
                     evals1, evals2, 
                     efuncs1, efuncs2, J=2, 
                     qt_int = c(0.025, 0.975)){
  
  # Usefule scalars
  t_max <- max(data$t) # end of observation point
  ku_vec <- table(data$visit) # number of observations in each visit
  Ju <- length(unique(data$visit)) # number of visits with observations
  K <- nrow(efuncs1) # maximum number of observation per visit
  L <- ncol(efuncs1) # number of PC
  M <- ncol(efuncs2)
  
  # data
  pred_data <- list(J = J, K = K,
                    Ju = Ju, Ku = ku_vec, 
                    Nobs = nrow(data), Y = data$Y,
                    L = L, M = M,
                    efuncs_l1 = efuncs1,
                    efuncs_l2 = efuncs2,
                    b0 = mu,
                    lambda = evals1,
                    gamma = evals2)
    
    # stan model
    fit <- stan(
      file = here("Code/prediction.stan"),  # Stan program
      data = pred_data,    # named list of data
      chains = 4,             # number of Markov chains
      warmup = 1000,          # number of warmup iterations per chain
      iter = 2000,            # total number of iterations per chain
      cores = 1,              # number of cores (could use one per chain)
      refresh = 0             # no progress shown
    )
    
    # scores
    sum_scores <- summary(fit)$summary
    xi_pred <- sum_scores[1:L, "mean"]
    zeta_pred <- sum_scores[-c(1:L, nrow(sum_scores)), "mean"]
    zeta_pred <- matrix(zeta_pred, M, J, byrow = T) # M by J
    
    # linear predictors
    eta_l1 <- efuncs1 %*% xi_pred
    eta_l2 <- as.vector(efuncs2 %*% zeta_pred)
    eta_pred <- rep(mu, J)+rep(eta_l1, J)+eta_l2
    eta_pred <- data.frame(visit = as.factor(rep(1:J, each = K)),
                           eta_pred) 
    
    # interval
    ## level 1
    samp_scores <- as.matrix(fit)
    eta_l1_samp <- efuncs1 %*% t(samp_scores[, 1:L]) # K by number of samples
    
    ## level 2
    eta_l2_sample <- array(NA, dim= c(nrow(samp_scores), K, J))
    for(i in 1:nrow(samp_scores)){
      this_zeta <-  matrix(samp_scores[i, -c(1:L, nrow(sum_scores))], 
                           M, J, byrow = T)
      eta_l2_sample[i,,] <- efuncs2 %*% this_zeta
      
    }
    
    ## whole
    eta_sample <- array(NA, dim= c(nrow(samp_scores), K, J))
    for(i in 1:nrow(samp_scores)){
      eta_sample[i,,] <- eta_l2_sample[i,,]+matrix(rep(eta_l1_samp[,i], J), ncol = J)
    }
    
    ## quantile
    eta_qt <- apply(eta_sample, 2:3, quantile, probs = qt_int) 
 
    ## clean
    eta_pred$eta_pred_lb <- as.vector(eta_qt[1,,])
    eta_pred$eta_pred_ub <- as.vector(eta_qt[2,,])
                    
    
    # output
    return(list(score1 = xi_pred, score2 = zeta_pred, eta_pred = eta_pred))
    
}


#### test ####




# df_test <- gen_bi(N = 1, J = 2, K = 100, t_vec = seq(0, 1, length.out = 100),
#                   L = 4, M = 4, lambda = 0.5^((1:4)-1), gamma =  0.5^((1:4)-1))
# 
# data <- df_test$data %>% filter(t<=0.5)

  