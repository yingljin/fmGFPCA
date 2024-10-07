
# This scirpt mimics a missing data scenario

#### Set up ####

library(here)
library(tidyverse)
library(lme4)
library(mgcv)
library(refund)
library(rstan)
library(instantiate)
library(GLMMadaptive)
library(splines)
library(visdat)

source(here("Code/Functions/Simulation.R"))
source(here("Code/Functions/gmFPCA_func.R"))
source(here("Code/Functions/ReEvalEfuncs.R"))
# source(here("Code/Functions/Predict.R"))

set.seed(915)
# devtools::install_github("julia-wrobel/fastGFPCA", build_vignettes = TRUE)

# toy simulation scheme
# N = 100 (no need to data split in this case), J = 3, K = 200
# bin every 10 observations
# number of FPC = 4

M <- 100
N <- 100 # sample size

J <- 3 # total number of visits
K <- 200


# 10% of the subjects have missing outcomes
# each subject have 50% of missing outcomes at random positions
p1 <- 0.1 
p2 <- 0.5

#### Generate data ####

# generate data with missingness

# containers
data_list_allM <- list()

for(m in 1:M){
  # generate training data
  data_all <- gen_bi(N = N, J = J, K = K, 
                     t_vec = seq(0, 1, length.out = K),
                     L = 4, M = 4, 
                     lambda = 0.5^((1:4)-1), 
                     gamma =  0.5^((1:4)-1))
  data_all$data$id <- as.factor(data_all$data$id)
  data_all$data$visit <- as.factor(data_all$data$visit)
  
  # introduce missing
  ## miss subjects
  miss_id <- sample(2:N, size = p1*N, replace = FALSE)
  ## missing observations
  data_all$data$Yobs <- data_all$data$Y
  for(i in seq_along(miss_id)){
    miss_t <- sample(J*K, size = p2*J*K, replace = FALSE)
    data_all$data$Yobs[data_all$data$id==miss_id[i]][miss_t] <- NA
  }
  
  data_list_allM[[m]] <- data_all
  
}

# check

data_list_allM[[1]]$data %>% #head()
  select(id, visit, Yobs) %>% 
  mutate(tid = rep(1:K, J*N)) %>%
  pivot_wider(values_from = Yobs, names_from = tid) %>%
  vis_miss()

data_list_allM[[67]]$data %>% 
  group_by(id) %>% 
  summarise(miss = sum(is.na(Yobs))) %>% 
  filter(miss>0)


data_list_allM[[1]]$data %>% 
  filter(id %in% c(21,24, 1, 2)) %>% 
  ggplot()+
  geom_line(aes(x=t, y=probs))+
  geom_point(aes(x=t, y=Yobs), size = 0.3)+
  facet_grid(row = vars(id), col=vars(visit))

data_list_allM_mis <- data_list_allM

save(data_list_allM_mis, 
     file = here("Data/SimData/SimData_mis.RData")) # data that will be used for all of simulation


#### gmFPCA ####

# load(here("Data/SimData.RData"))

load(here("Data/SimData/SimData_mis.RData"))




comp_time_vec <- rep(NA, M)
comp_list_allM <- list()

# simulation
for(m in 1:M){
  
  data_all <- data_list_allM_mis[[m]]$data
  
  data_tr <- data_all %>% select(-Y) %>% rename(Y=Yobs) 
  
  # assign missing (by row index)
    
  tic <- Sys.time()
  # gmFPCA
  gmfpca_fit <- gm_FPCA(data = data_tr, bin_w=10, L = 4, scores = T)
  
  # fill in missing values
  miss_id <- data_all %>% group_by(id) %>% summarize(nmis = sum(is.na(Yobs))) %>%
    filter(nmis>0) %>% select(id) %>% unlist()
  
  ## container
  comp_df_m <- data_all %>% filter(id %in% miss_id) %>% mutate(eta_comp = NA)
  
  ## computation
  for(i in seq_along(miss_id)){ # for a subject
    
    # we don't even need to create the out-of-sample context
    # because we will be able to compute the scores directly from gmFPCA
    score1 <- gmfpca_fit$scores_l1[miss_id[i],]
    
    for(j in 1:J){ # a given observation track
    
      
      score2<- gmfpca_fit$scores_l2[(as.numeric(miss_id[i])-1)*J+j, ]
      # print(paste(i, j))
      # print(length(gmfpca_fit$mu))
      # print(length(gmfpca_fit$efuncs_l1 %*% score1))
      # print(length(gmfpca_fit$efuncs_l2 %*% score2))

      # eta <-
      #   tryCatch(
      #     {gmfpca_fit$mu + as.vector(gmfpca_fit$efuncs_l1 %*% score1) + as.vector(gmfpca_fit$efuncs_l2 %*% score2)},
      #     warning=function(w) {
      #       stop("converted from warning: ", conditionMessage(w))
      #     }
      #   )
      eta <- gmfpca_fit$mu + as.vector(gmfpca_fit$efuncs_l1 %*% score1) + as.vector(gmfpca_fit$efuncs_l2 %*% score2)
      # put in container
      comp_df_m[comp_df_m$id==miss_id[i]&comp_df_m$visit==j, "eta_comp"] <- eta

    }
}
  toc <- Sys.time()
  # computation time
  comp_time_vec[m] <- difftime(toc, tic, units = "mins")
  
  comp_list_allM[[m]] <- comp_df_m
  
  # print progress
  print(paste0(m, "/", M, " simulation completed"))
}

# check
unique(comp_df_m$id)

# comp_df_m %>% 
comp_list_allM[[100]] %>%
  filter(id %in% sample(unique(comp_list_allM[[100]]$id), 4)) %>% 
  # filter(id==41) %>%
  mutate_at(vars(starts_with("eta")), function(x)(exp(x)/(1+exp(x)))) %>%
  mutate(eta=ifelse(is.na(Yobs), NA, eta),
         eta_comp=ifelse(is.na(Yobs), eta_comp, NA)) %>% 
  ggplot()+
  geom_point(aes(x=t, y=eta, col="True"), size = 0.1)+
  geom_point(aes(x=t, y=eta_comp, col="Computed"), size = 0.1)+
  geom_point(aes(x=t, y=Yobs, col = "Observed"), size = 0.2)+
  # geom_line(aes(x=t, y=I(eta_pred-0.96*eta_sd), col = "Pred"), linetype = "dashed")+
  # geom_line(aes(x=t, y=I(eta_pred+0.96*eta_sd), col = "Pred"), linetype = "dashed")+
  facet_grid(rows = vars(id), cols = vars(visit))


save(comp_list_allM, comp_time_vec, 
     file = here("Data/SimOutput/SimOutput_miss.RData"))


mean(comp_time_vec) # 0.4 minutes each simulation



# When K = 100, a few numeric problems happened: 
## Tail Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
## estimated Bayesian Fraction of Missing Information was low
## chains have not mixed.
## Singularity of local models
## Some datasets had very, very long prediction time





#### GLMMadaptive ####

load(here("Data/SimData/SimData_mis.RData"))




# M <- 5
# containers
comp_time_vec_ref <- rep(NA, M)

comp_list_allM_ref <- list()

 # m <- 1
# simulation
for(m in 1:M){
  
  data_all <- data_list_allM_mis[[m]]$data %>% 
    mutate(id_visit = paste0(id, "_", visit))
  
  # model fit
  tic <- Sys.time()
  # mixed model with nested grouping effect
  ## technically, this is not nesting. But the real nesting one is not fittable
  # with spline basis? 
  glmm_fit <- tryCatch({mixed_model(fixed = Yobs ~ t,
                          random = ~ t | id_visit,
                          data = data_all,
                    family = binomial())},
                    error=function(e){
                      return(list(coefficients=NA))})
  # toc <- Sys.time()
  
  if(!is.na(glmm_fit$coefficients[1])){
  
   # prediction
    miss_id <- data_all %>% group_by(id) %>% 
      summarize(nmis = sum(is.na(Yobs))) %>%
      filter(nmis>0) %>% select(id) %>% unlist()
    
    ## predict conditioning on observed outcomes
    ## container
    comp_df_m <- data_all %>% filter(id %in% miss_id) %>% mutate(eta_comp = NA)
    fixed_coef <- fixef(glmm_fit)
    rand_coef <- ranef(glmm_fit, type = "subject_specific")
    Xmat <- cbind(rep(1, K), unique(data_all$t))
    
    ## computation
    for(i in seq_along(miss_id)){ # for a subject
      
      # we don't even need to create the out-of-sample context
      # because we will be able to compute the scores directly from gmFPCA
      
      for(j in 1:J){ # a given observation track
        
        coef_ij <- rand_coef[paste0(miss_id[i], "_", j), ]
        
      
        eta <- Xmat %*% fixed_coef + Xmat %*% coef_ij
        # put in container
        comp_df_m[comp_df_m$id==miss_id[i]&comp_df_m$visit==j, "eta_comp"] <- eta
        
      }
    }
    
    toc <- Sys.time()
    comp_time_vec_ref[m] <- difftime(toc, tic, units = "mins")
    comp_list_allM_ref[[m]] <- comp_df_m
  }
  
  else(print(paste("Skip", m)))
  
}
  
  

# skiped the 43th dataset
  
  
  

  

  
  
 

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
mean(comp_time_vec_ref, na.rm = T) # 0.6 minutes each simulation
 # 0.0115 minutes each simulation

# check
comp_list_allM_ref[[100]] %>%
# comp_df_m %>% 
  filter(id %in% sample(unique(comp_list_allM_ref[[100]]$id), 4)) %>% #head()
  # filter(id==41) %>%
  mutate_at(vars(starts_with("eta")), function(x)(exp(x)/(1+exp(x)))) %>%
  mutate(eta=ifelse(is.na(Yobs), NA, eta),
       eta_comp=ifelse(is.na(Yobs), eta_comp, NA)) %>% 
  ggplot()+
  geom_point(aes(x=t, y=eta, col="True"), size = 0.1)+
  geom_point(aes(x=t, y=eta_comp, col="Computed"), size = 0.1)+
  geom_point(aes(x=t, y=Yobs, col = "Observed"), size = 0.2)+
  # geom_line(aes(x=t, y=I(eta_pred-0.96*eta_sd), col = "Pred"), linetype = "dashed")+
  # geom_line(aes(x=t, y=I(eta_pred+0.96*eta_sd), col = "Pred"), linetype = "dashed")+
  facet_grid(rows = vars(id), cols = vars(visit))
sum(is.na(pred_df_m$eta_pred_J3))



pred_df_m %>% filter(visit==3)

save(comp_list_allM_ref, comp_time_vec_ref, 
     file = here("Data/SimOutput/SimOutput_miss_ref.RData"))


