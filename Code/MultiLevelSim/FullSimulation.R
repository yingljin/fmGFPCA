
#### Set up ####

library(here)
library(tidyverse)
library(lme4)
library(mgcv)
library(refund)
library(rstan)
# library(instantiate)
library(GLMMadaptive)
library(splines)

source(here("Code/Functions/Simulation.R"))
source(here("Code/Functions/gmFPCA_func.R"))
source(here("Code/Functions/ReEvalEfuncs.R"))
# source(here("Code/Functions/Predict.R"))

set.seed(915)
# devtools::install_github("julia-wrobel/fastGFPCA", build_vignettes = TRUE)

# toy simulation scheme
# N_train = 100,  N_test = 100, J = 3, K = 200
# observations stops at the the midpoint of each visit
# bin every 10 observations
# number of FPC = 4

M  <- 100
Ntr <- 100
Nte <- 100

J <- 2 # total number of visits
K <- 200

#### Generate data ####

# containers
data_list_allM_J200 <- list()

for(m in 1:M){
  # generate training data
  data_all <- gen_bi(N = Ntr+Nte, J = J, K = K, 
                     t_vec = seq(0, 1, length.out = K),
                     L = 4, M = 4, 
                     lambda = 0.5^((1:4)-1), 
                     gamma =  0.5^((1:4)-1))
  data_all$data$id <- as.factor(data_all$data$id)
  data_all$data$visit <- as.factor(data_all$data$visit)
  data_list_allM_J200[[m]] <- data_all$data
}

# check

data_list_allM_J200[[1]] %>% filter(id %in% sample(1:200, 2)) %>% 
  ggplot()+
  geom_line(aes(x=t, y=probs))+
  geom_point(aes(x=t, y=Y), size = 0.5)+
  facet_grid(row = vars(id), col=vars(visit))

# save(data_list_allM_J200, 
#      file = here("Data/SimData/SimData_J200.RData")) # data that will be used for all of simulation


#### gmFPCA ####

# load(here("Data/SimData.RData"))

# load(here("Data/SimData/SimData_J200.RData"))

tmax <- 0.5 # max time of observation

M <- 3

fit_time_vec <- rep(NA, M)
pred_time_vec <- rep(NA, M)
pred_list_allM <- list()

# data
data_list_allM <- data_list_allM_J200
rm(data_list_allM_J200)


# simulation
for(m in 1:M){
  
  data_all <- data_list_allM[[m]]
  
  # data split
  ## training '
  data_tr <- data_all %>% filter(id %in% 1:Ntr) %>% 
    mutate(id = droplevels(id))
  ## test
  data_te <- data_all %>% filter(!id %in% 1:Ntr) %>%
    mutate(id = droplevels(id), visit=droplevels(visit))
    
  tic <- Sys.time()
  # gmFPCA
  gmfpca_fit <- gm_FPCA(data = data_tr, bin_w=10, L = 4)
  toc <- Sys.time()
  fit_time_vec[m] <- difftime(toc, tic, units = "mins")
  
  
  # prediction
  te_id <- unique(data_te$id) # unique test id
  
  ## container
  pred_df_m <- data_all %>% filter(id %in% te_id) %>% 
    mutate(eta_pred_J1 = NA, eta_pred_J2 = NA, eta_pred_J3 = NA,
           id = droplevels(id))
  
  ## prediction
  tic <- Sys.time()
  for(i in seq_along(te_id)){ # for a subject
    
    for(j in 1:J){ # a given observation track
    
      # prediction conditioning on half visits from the jth visit
      # for visits before the jth visit, full track
      # for visits from the jth visit, half track
      df_te_ij <- data_te %>% filter(id == te_id[i]) %>% 
        filter(as.numeric(visit)<j | t<=tmax) 
      
      pred_data <- list(J = J, K = K,
                        Ju = J, Ku = table(df_te_ij$visit), 
                        Nobs = nrow(df_te_ij), Y = df_te_ij$Y,
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
    pred_df_m[pred_df_m$id == te_id[i], 
              paste0("eta_pred_J", j)] <- eta_pred
    
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
  }
  toc <- Sys.time()
  # computation time
  pred_time_vec[m] <- difftime(toc, tic, units = "mins")
  
  pred_list_allM[[m]] <- pred_df_m
  
  # print progress
  print(paste0(m, "/", M, " simulation completed"))
}

# check
pred_list_allM[[1]] %>% #head()
# pred_df_m %>%
  filter(id %in% sample(101:200, 2)) %>% #head()
  # filter(id==41) %>%
  mutate_at(vars(starts_with("eta_pred")), function(x)(exp(x)/(1+exp(x)))) %>%
  ggplot()+
  geom_line(aes(x=t, y=probs, col = "True"))+
  geom_line(aes(x=t, y=eta_pred_J1, col = "J1"), linetype = "dashed")+
  geom_line(aes(x=t, y=eta_pred_J2, col = "J2"), linetype = "dashed")+
  # geom_line(aes(x=t, y=eta_pred_J3, col = "J3"), linetype = "dashed")+
  # geom_line(aes(x=t, y=I(eta_pred-0.96*eta_sd), col = "Pred"), linetype = "dashed")+
  # geom_line(aes(x=t, y=I(eta_pred+0.96*eta_sd), col = "Pred"), linetype = "dashed")+
  facet_grid(rows = vars(id), cols = vars(visit))


save(fit_time_vec, pred_time_vec, pred_list_allM, 
     file = here("Data/SimOutput/SimOutput_J200.RData"))


mean(fit_time_vec) # 0.4 minutes each simulation
mean(pred_time_vec) # 5.4 minutes each simulation


# When K = 100, a few numeric problems happened: 
## Tail Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
## estimated Bayesian Fraction of Missing Information was low
## chains have not mixed.
## Singularity of local models
## Some datasets had very, very long prediction time









#### GLMMadaptive ####
load(here("Data/SimData/SimData_J200.RData"))

tmax <- 0.5 # max time of observation

# M <- 5
# containers
fit_time_vec_ref <- rep(NA, M)
pred_time_vec_ref <- rep(NA, M)
pred_list_allM_ref <- list()

# data
data_list_allM <- data_list_allM_J200
rm(data_list_allM_J200)

m <- 1
# simulation
for(m in 1:M){
  
  data_all <- data_list_allM[[m]] %>% 
    mutate(id_visit = paste0(id, "_", visit))
  
  # data split
  ## training '
  data_tr <- data_all %>% filter(id %in% 1:Ntr) %>% 
     mutate(id = droplevels(id))
  ## test
  data_te <- data_all %>% filter(!id %in% 1:Ntr) %>%
    mutate(id = droplevels(id))
  
  tic <- Sys.time()
  # mixed model with nested grouping effect
  ## start with a nested, intercept-only model
  glmm_fit <- mixed_model(fixed = Y ~ 1,
                          random = ~ 1  | id/visit,
                          data = data_tr,
                    family = binomial())
  toc <- Sys.time()
  fit_time_vec_ref[m] <- difftime(toc, tic, units = "mins")
  
  
  # prediction
  tic <- Sys.time()
  ## predict conditioning three half visits
  pred_J1 <- predict(glmm_fit, 
                     newdata = data_te %>% filter(t <= tmax),
                     newdata2 = data_te,
                     type_pred = "link", type = "subject_specific",
                     return_newdata = T)$newdata2
  
  ## predict conditioning on visit 1 and half of visit 3
  pred_J2 <- predict(glmm_fit, 
                     newdata = data_te %>% filter(visit==1|(visit %in% 2:3 & t<=tmax)),
                     newdata2 = data_te,
                     type_pred = "link", type = "subject_specific",
                     return_newdata = T)$newdata2

  ## predict conditioning on visit 1 and half of visit 2
  pred_J3 <- predict(glmm_fit, 
                     newdata = data_te %>% filter(visit %in% 1:2|(visit==3&t<=tmax)),
                     newdata2 = data_te,
                     type_pred = "link", type = "subject_specific",
                     return_newdata = T)$newdata2

  toc <- Sys.time()
  pred_time_vec_ref[m] <- difftime(toc, tic, units = "mins")

  pred_df_m <- pred_J3 %>% select(id, visit, t, eta, Y, pred) %>% 
    rename(eta_pred_J3=pred) %>% 
    left_join(pred_J2 %>% select(id, visit, t, pred) %>% rename(eta_pred_J2=pred),
              by = c("id", "visit", "t")) %>% 
    left_join(pred_J1 %>% select(id, visit, t, pred) %>% rename(eta_pred_J1=pred),
             by = c("id", "visit", "t")) 

  pred_list_allM_ref[[m]] <- pred_df_m
  
  # print progress
  print(paste0(m, "/", M, " simulation completed"))
}



# GLMM model: 
## If use nested design (id/visit), R cannot handle more than an intercept
## It also had problems with out-of-sample prediction
## fixed = Y ~ t, random = ~ t | id/visit: vector memory exhausted
## fixed = Y ~ 1, random = ~ 1 | id/visit: take 10 minutes

## If use cross design (id_visit), then cannot predict unobserved track
## fixed = Y ~ t+I(t^2), random = ~ t + I(t^2) | id_visit, takes about 3 minutes? 
## ns(df=2) is the same

## it does not seem to handle nested models very well
## can't predict unobserved visits
## We need a better reference method
## Maybe use Alpine? 
mean(fit_time_vec_ref) # 0.6 minutes each simulation
mean(pred_time_vec_ref) # 0.0115 minutes each simulation

# check
pred_list_allM_ref[[100]] %>% #head()
#pred_df_m %>% 
  filter(id %in% sample(101:200, 4)) %>% #head()
  # filter(id==41) %>%
  mutate_at(vars(starts_with("eta")), function(x)(exp(x)/(1+exp(x)))) %>%
  ggplot()+
  geom_line(aes(x=t, y=eta, col = "True"))+
  geom_line(aes(x=t, y=eta_pred_J1, col = "J1"), linetype = "dashed", alpha = 0.5)+
  geom_line(aes(x=t, y=eta_pred_J2, col = "J2"), linetype = "dashed", alpha = 0.5)+
  geom_line(aes(x=t, y=eta_pred_J3, col = "J3"), linetype = "dashed", alpha = 0.5)+
  # geom_line(aes(x=t, y=I(eta_pred-0.96*eta_sd), col = "Pred"), linetype = "dashed")+
  # geom_line(aes(x=t, y=I(eta_pred+0.96*eta_sd), col = "Pred"), linetype = "dashed")+
  facet_grid(rows = vars(id), cols = vars(visit))

sum(is.na(pred_df_m$eta_pred_J3))



pred_df_m %>% filter(visit==3)

save(fit_time_vec_ref, pred_time_vec_ref, pred_list_allM_ref, 
     file = here("Data/SimOutput/SimOutput_J200_ref.RData"))


