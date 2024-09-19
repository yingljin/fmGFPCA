library(here)
library(tidyverse)
library(lme4)
library(mgcv)
library(refund)
library(rstan)
library(instantiate)


source(here("Code/Functions/Simulation.R"))
source(here("Code/Functions/gmFPCA_func.R"))
source(here("Code/Functions/ReEvalEfuncs.R"))
# source(here("Code/Functions/Predict.R"))

set.seed(915)
# devtools::install_github("julia-wrobel/fastGFPCA", build_vignettes = TRUE)

# toy simulation scheme
# N_train = 40,  N_test = 10, J = 3, K = 100
# observations stops at the the midpoint of each visit
# bin every 10 observations
# number of FPC = 4

#### Simulation ####

M <- 5
Ntr <- 100
Nte <- 100

Jmax <- 3 # max number of observed visits
J <- 3 # total number of visits

# containers
fit_time_vec <- rep(NA, M)
pred_time_vec <- rep(NA, M)
pred_list_allM <- list()
data_list_allM <- list()

for(m in 1:M){
  
  # generate training data
  data_all <- gen_bi(N = Ntr+Nte, J = J, K = 100, 
                    t_vec = seq(0, 1, length.out = 100),
                    L = 4, M = 4, 
                    lambda = 0.5^((1:4)-1), 
                    gamma =  0.5^((1:4)-1))
  data_all$data$id <- as.factor(data_all$data$id)
  data_all$data$visit <- as.factor(data_all$data$visit)
  data_list_allM[[m]] <- data_all$data
  
  # data split
  ## training '
  data_tr <- data_all$data %>% filter(id %in% 1:Ntr) %>% 
    mutate(id = droplevels(id))
  ## test
  data_te <- data_all$data %>% filter(!id %in% 1:Ntr) %>%
    filter(!(visit==Jmax & t>0.5)) %>%  # truncate
    mutate(id = droplevels(id), visit=droplevels(visit))
    
  tic <- Sys.time()
  # gmFPCA
  gmfpca_fit <- gm_FPCA(data = data_tr, bin_w=10, L = 4)
  toc <- Sys.time()
  fit_time_vec[m] <- difftime(toc, tic, units = "mins")
  
  
  # prediction
  te_id <- unique(data_te$id) # unique test id
  
  ## container
  pred_df_m <- data_all$data %>% filter(id %in% te_id) %>% 
    mutate(eta_pred = NA, # eta_sd = NA,
           id = droplevels(id))
  tic <- Sys.time()
  for(i in seq_along(te_id)){
    
    df_te_i <- data_te %>% filter(id == te_id[i]) 
    
    pred_data <- list(J = J, K = 100,
                      Ju = Jmax, Ku = table(df_te_i$visit), 
                      Nobs = nrow(df_te_i), Y = df_te_i$Y,
                      L = 4, M = 4,
                      efuncs_l1 = gmfpca_fit$efuncs_l1,
                      efuncs_l2 = gmfpca_fit$efuncs_l2,
                      b0 = gmfpca_fit$mu,
                      lambda = gmfpca_fit$evals1,
                      gamma = gmfpca_fit$evals2)
    
    # stan model
    fit <- stan(
      file = here("Code/prediction.stan"),  # Stan program
      data = pred_data,    # named list of data
      chains = 2,             # number of Markov chains
      warmup = 1000,          # number of warmup iterations per chain
      iter = 2000,            # total number of iterations per chain
      cores = 1,              # number of cores (could use one per chain)
      refresh = 0             # no progress shown
    )
    
    # scores
    sum_scores <- summary(fit)$summary[1:(4*(J+1)), "mean"]
    xi_pred <- sum_scores[1:4]
    zeta_pred <- sum_scores[-c(1:4)]
    zeta_pred <- matrix(zeta_pred, 4, J, byrow = T) # M by J
    
    # linear predictors
    eta_l1 <- gmfpca_fit$efuncs_l1 %*% xi_pred
    eta_l2 <- as.vector(gmfpca_fit$efuncs_l2 %*% zeta_pred)
    eta_pred <- rep(gmfpca_fit$mu, J)+rep(eta_l1, J)+eta_l2
    pred_df_m[pred_df_m$id == te_id[i], "eta_pred"] <- eta_pred
    
    # standard deviation
    # scores
    # sum_sd <- summary(fit)$summary[1:24, "sd"]
    # xi_sd <- sum_sd[1:4]
    # zeta_sd <- sum_sd[-c(1:4)]
    # zeta_sd <- matrix(zeta_sd, 4, 5, byrow = T) # M by J
    # 
    # # linear predictors
    # var_l1 <- gmfpca_fit$efuncs_l1^2 %*% xi_sd^2
    # var_l2 <- as.vector(gmfpca_fit$efuncs_l2^2 %*% zeta_sd^2)
    # var_eta <- rep(var_l1, 5)+var_l2
    # pred_df_m[pred_df_m$id == te_id[i], "eta_sd"] <- sqrt(var_eta)
  }
  toc <- Sys.time()
  # computation time
  pred_time_vec[m] <- difftime(toc, tic, units = "mins")
  
  pred_list_allM[[m]] <- pred_df_m
  
  # print progress
  print(paste0(m, "/", M, " simulation completed"))
}

save(fit_time_vec, pred_time_vec, pred_list_allM, 
     file = here("Data/Sim_J3.RData"))


# There are numeric problems here:
## Tail Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
## estimated Bayesian Fraction of Missing Information was low
## chains have not mixed.
## Singularity of local models
## I am only able to run very, very small simulation because of the aborting issue

#### Figures ####

# data
data_tr$data %>% filter(id %in% 16:19) %>%
  ggplot()+
  geom_line(aes(x=t, y=probs))+
  geom_point(aes(x=t, y=Y), size = 0.5)+
  facet_grid(row = vars(id), col=vars(visit))

pred_list_allM[[1]] %>% 
  filter(id %in% sample(101:200, 4)) %>%
  # filter(id==41) %>%
  mutate_at(vars(starts_with("eta")), function(x)(exp(x)/(1+exp(x)))) %>%
  ggplot()+
  geom_line(aes(x=t, y=eta, col = "True"))+
  geom_line(aes(x=t, y=eta_pred, col = "Pred"))+
  # geom_line(aes(x=t, y=I(eta_pred-0.96*eta_sd), col = "Pred"), linetype = "dashed")+
  # geom_line(aes(x=t, y=I(eta_pred+0.96*eta_sd), col = "Pred"), linetype = "dashed")+
  facet_grid(rows = vars(id), cols = vars(visit))

mean(fit_time_vec) # 0.5 minutes each simulation
mean(pred_time_vec) # 2 minutes each simulation


