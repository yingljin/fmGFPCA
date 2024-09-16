library(here)

source(here("Code/Functions/Simulation.R"))
source(here("Code/Functions/gmFPCA_func.R"))
source(here("Code/Functions/ReEvalEfuncs.R"))
source(here("Code/Functions/Predict.R"))

set.seed(915)


# toy simulation scheme
# N_train = 80,  N_test = 20, J = 2, K = 100
# observations stops at the second visit
# specific stopping time is generated randomly
# bin every 10 observations
# number of FPC = 4

#### Simulation ####

M <- 10
Ntr <- 40
Nte <- 10

# containers
comp_time_vec <- rep(NA, M)
pred_list_allM <- list()

for(m in 1:M){
  
  # generate training data
  data_tr <- gen_bi(N = Ntr, J = 2, K = 100, 
                    t_vec = seq(0, 1, length.out = 100),
                    L = 4, M = 4, 
                    lambda = 0.5^((1:4)-1), 
                    gamma =  0.5^((1:4)-1))
  data_tr$data$id <- as.factor(data_tr$data$id)
  data_tr$data$visit <- as.factor(data_tr$data$visit)
  
  tic <- Sys.time()
  # gmFPCA
  gmfpca_fit <- gm_FPCA(data = data_tr$data, bin_w=10, L = 4)
  
  # prediction
  ## generate test data
  data_te <- gen_bi(N = Nte, J = 2, K = 100, 
                    t_vec = seq(0, 1, length.out = 100),
                    L = 4, M = 4, 
                    lambda = 0.5^((1:4)-1), 
                    gamma =  0.5^((1:4)-1))
  data_te$data$id <- as.factor(data_te$data$id)
  data_te$data$visit <- as.factor(data_te$data$visit)
  
  
  ## truncate observations
  tmax_vec <- runif(Nte, 0, 1)
  
  ## predict subject by subject
  pred_list <- list()
  for(i in seq_along(unique(data_te$data$id))){
    
    df_test_i <- data_te$data %>% filter(id == i) %>%
      filter(visit == 1 | t<=tmax_vec[i] )
    
    pred_data <- mDynPred(data = df_test_i, 
                     mu = gmfpca_fit$mu,
                     evals1 = gmfpca_fit$evals1,
                     evals2 = gmfpca_fit$evals2,
                     efuncs1 =  gmfpca_fit$efuncs_l1,
                     efuncs2 = gmfpca_fit$efuncs_l2,
                     J = 2, 
                     warmup = 1000,
                     iter = 2000)
    
    pred_list[[i]] <- data_te$data %>% filter(id == i) %>% 
      mutate(eta_pred = pred_data$df_pred$eta_pred,
             eta_pred_lb = pred_data$df_pred$eta_pred_lb,
             eta_pred_ub = pred_data$df_pred$eta_pred_ub) 
  
    ## remove observed part
    pred_list[[i]][pred_list[[i]]$t<=tmax_vec[i] | pred_list[[i]]$visit==1, 
                   c("eta_pred", "eta_pred_lb", "eta_pred_ub")] <- NA
  }
  
  toc <- Sys.time()
  # computation time
  comp_time_vec[m] <- difftime(toc, tic, units = "mins")
  
  pred_list_allM[[m]] <- bind_rows(pred_list)
  
  # print progress
  print(paste0(m, "/", M, " simulation completed"))
}


# There are numeric problems here:
## Tail Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
## estimated Bayesian Fraction of Missing Information was low
## chains have not mixed.

#### Figures ####

# data
data_tr$data %>% filter(id %in% 16:19) %>%
  ggplot()+
  geom_line(aes(x=t, y=probs))+
  geom_point(aes(x=t, y=Y), size = 0.5)+
  facet_grid(row = vars(id), col=vars(visit))

pred_list_allM[[2]] %>% #head()
  filter(id %in% sample(1:20, 4)) %>% 
  mutate_at(vars(starts_with("eta")), function(x)(exp(x)/(1+exp(x)))) %>%
  ggplot()+
  geom_line(aes(x=t, y=eta, col = "True"))+
  geom_line(aes(x=t, y=eta_pred, col = "Pred"))+
  geom_line(aes(x=t, y=eta_pred_ub, col = "Pred"), linetype = "dashed")+
  geom_line(aes(x=t, y=eta_pred_lb, col = "Pred"), linetype = "dashed")+
  facet_grid(rows = vars(id), cols = vars(visit))

pred_list_allM[[10]] %>% View()

