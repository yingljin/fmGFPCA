# This script implements the full simulation 
# using single level fGFPCA
# treating the multiple visits as one very long series
# as a reference method

#### set up ####

set.seed(120)

# options(mc.cores=parallel::detectCores())

library(here)
library(tidyverse)
library(ggpubr)
library(refund)
library(lme4)
library(ROCR)
library(GLMMadaptive)
library(kableExtra)
library(knitr)
library(mvtnorm)
library(mgcv)
library(splines)
library(rstan)
theme_set(theme_minimal())


#### load simulated data ####

load(here("Data/SimData/SimData.RData")) # data generated from Code/GenData.R

N <- length(unique(data_list[[1]]$id))
J <- length(unique(data_list[[1]]$t))*length(unique(data_list[[1]]$visit))
# since we are treating the multiple visit as one series
# number of repeated measures would be number of visits * number of measures within each visit
t <- unique(data_list[[1]]$t)
M <- length(data_list)
K <- 4 # number of eigenfunctions to use


#### fGFPCA Model fitting ####

# functions
source(here("Code/Functions/SingleLevelLocalGLMM.R")) 
# source(here("Code/OutSampBayes.R"))

# split data into train and test set
N_train <- 100
N_test <- 100

# binning 
bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points on time index
mid <- (brks+bin_w/2)[1:n_bin] # interval mid point time index
# because we are treating multiple visitis as one major visit
# the original time from [0,1] cannot be used
# we need to use a delegate time index from 1:400

# values used for projection in step 4
p <- 3 # order of b splines 
knots <- 35 # number of knots (same from FPCA model)
knots_values <- seq(-p, knots + p, length = knots + 1 + 2 *p)/knots
knots_values <- knots_values * (max(mid) - min(mid)) + min(mid)

# result container
M <- 3 # try on a few datasets first

# M
pred_list_fGFPCA <- list()
# score_list_all <- list()
# sd_list_all <- list()
# converge_state_list <- list()
time_fGFPCA <- rep(NA, M)


pb = txtProgressBar(min = 0, max = M, initial = 0, style = 3) 

# The whole process
for(m in 1:M){
  
  # for every simulated dataset
  df <- data_list[[m]]
  
  ## Step 1: bin observations
  df$sind <- rep(1:J, N)
  df$bin <- cut(df$sind, breaks = brks, include.lowest = T, labels = mid)
  df$bin <- as.numeric(as.character(df$bin))
  
  # split into training and test set
  train_df <- df %>% filter(id %in% 1:N_train)
  test_df <- df %>% filter(!id %in% 1:N_train)
  
  # fit fGFPCA model on the training set
  t1 <- Sys.time()
  ## Step 2: local GLMM and estimate latent function
  ## on training set
  bin_lst <- split(train_df, f = train_df$bin)
  df_est_latent <- lapply(bin_lst, function(x){pred_latent(x, n_node = 0)}) 
  df_est_latent <- bind_rows(df_est_latent) 
  ## Step 3: FPCA
  uni_eta_hat <- df_est_latent %>% filter(bin==sind)
  mat_est_unique <- matrix(uni_eta_hat$eta_hat,
                           nrow=N_train, 
                           ncol=n_bin, byrow = F) 
  fpca_mod <- fpca.face(mat_est_unique, argvals = mid, var=T)
  ## Step 4: project and debias
  ## Projection
  B <- spline.des(knots = knots_values, x = mid, ord = p + 1,
                  outer.ok = TRUE)$design  # evaluate B-splines on binned grid
  Bnew <- spline.des(knots = knots_values, x = 1:J, ord = p + 1,
                     outer.ok = TRUE)$design  # evaluate B-splines on original grid
  df_phi <- matrix(NA, J, K) 
  for(k in 1:K){
    lm_mod <- lm(fpca_mod$efunctions[,k] ~ B-1)
    df_phi[,k] <- Bnew %*% coef(lm_mod)
  }# project binned eigenfunctions onto the original grid
  
  ## debias
  df_phi <- data.frame(sind = 1:J, df_phi)
  colnames(df_phi) <- c("sind", paste0("phi", 1:4))
  train_df <- train_df %>% left_join(df_phi, by = "sind")
  train_df$id <- as.factor(train_df$id)
  debias_glmm <- bam(Y ~ s(sind, bs="cc", k=10)+
                       s(id, by=phi1, bs="re")+
                       s(id, by=phi2, bs="re")+
                       s(id, by=phi3, bs="re")+
                       s(id, by=phi4, bs="re"), 
                     family = binomial, data=train_df, 
                     method = "fREML",
                     discrete = TRUE) # about half second
  new_mu <- predict(debias_glmm, type = "terms")[1:J, 1]+coef(debias_glmm)[1] # extract re-evaluated mean
  ## extract variacne from the debiased model
  new_lambda <- 1/debias_glmm$sp[2:5] # extract re-evaluated lambda
  
  # rescale
  new_phi <- df_phi %>% select(starts_with("phi"))*sqrt(n_bin)
  new_phi <- as.matrix(new_phi)
  new_lambda <- new_lambda/n_bin
  # t2 <- Sys.time()
  # fit_time[m] <- difftime(t2, t1, units = "mins")
  
  # Prediction pf test sample using stan
  # container
  pred_list_m <- list()
  # converge_state_m <- matrix(NA, N_test, 4)
  # score_list_m <- array(NA, dim = c(N_test, K, 4))
  # sd_list_m <- array(NA, dim = c(N_test, K, 4)) # subject by PC by case
  test_id <- unique(test_df$id)
  tmax_vec <- c(50, 150)
  
  # t1 <- Sys.time()
  for(i in seq_along(test_id)){
    
    df_i <- test_df %>% filter(id==test_id[i])
    
    # container for one subject
    # scores <- matrix(NA, K, length(tmax_vec))
    # sd_scores <- matrix(NA, K, length(tmax_vec))
    
    for(tid in seq_along(tmax_vec)){
      
      tmax <- tmax_vec[tid]
      df_it <- df_i %>% filter(sind <= tmax)
      
      # into a list
      stanData <- list(
        J = J, Ju = nrow(df_it), Y = df_it$Y, K = K,
        efuncs = new_phi, 
        b0 = new_mu,
        lambda = new_lambda
      )
      
      fit <- stan(
        file = here("Code/Functions/SingleLevelPred.stan"),  # Stan program
        data = stanData,    # named list of data
        chains = 2,             # number of Markov chains
        warmup = 1000,          # number of warmup iterations per chain
        iter = 2000,            # total number of iterations per chain
        cores = 1,              # number of cores (could use one per chain)
        refresh = 0             # no progress shown
      )
      
      # point prediction
      scores_tmax <- summary(fit)$summary[1:K, "mean"]
      eta_pred_out <- new_mu+new_phi%*%scores_tmax
      df_i[ , paste0("pred", tmax)] <- eta_pred_out[,1]
      
      # prediction interval using sampling quantiles
      # score_draws <- as.matrix(fit)[, 1:4]
      # eta_draws <- new_phi %*% t(score_draws)
      # eta_draws <- apply(eta_draws, 1, quantile, probs = c(0.025, 0.975))
      # quantile interval
      # df_i[ , paste0("pred", tmax, "_lb")] <- as.vector(new_mu+eta_draws[1, ])
      # df_i[ , paste0("pred", tmax, "_ub")] <- as.vector(new_mu+eta_draws[2, ])
      
      
    }
    
   pred_list_m[[i]] <- df_i
    
  }
  t2 <- Sys.time()
  # t_pred <- t2-t1 # About 3.5 minutes
  time_fGFPCA[m] <- difftime(t2, t1, units = "mins")
  
  # save results
  pred_list_fGFPCA[[m]]<- bind_rows(pred_list_m)
  # converge_state_list[[m]] <- converge_state_m
  
  setTxtProgressBar(pb, m)
}

close(pb)


#### Check output ####

# time 
mean(time_fGFPCA)

# prediction
# figure
rand_id <- sample(test_id, 4)

df_exp <- pred_list_fGFPCA[[1]] %>% 
  filter(id == 101) %>%
  mutate_at(vars(starts_with("pred")), 
            function(x){exp(x)/(1+exp(x))})  
df_exp[df_exp$sind<=50, "pred50"] <- NA
df_exp[df_exp$sind<=150, "pred150"] <- NA


df_exp %>%
  ggplot()+
  geom_point(aes(x=t, y=Y), size = 0.2)+
  geom_line(aes(x=t, y=probs, col = "True"))+
  geom_line(aes(x=t, y=pred50, col = "0.2"), na.rm = T)+
  geom_line(aes(x=t, y=pred150, col = "0.4"), na.rm = T)+
  facet_grid(rows = vars(id), cols = vars(visit))+
  guides(alpha = "none", col="none")
  

#### Save results ####

save(pred_list_fGFPCA, time_fGFPCA,
     file = here("Data/SimOutput/SimOutputSingle.RData"))

# save(pred_list_fGFPCA, time_fGFPCA,
#      file = here("Data/SimN500/SimOutput_fGFPCA.RData"))
