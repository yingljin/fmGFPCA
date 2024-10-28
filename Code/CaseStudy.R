
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

set.seed(1022)
# devtools::install_github("julia-wrobel/fastGFPCA", build_vignettes = TRUE)


#### Laod and clean NHANES data ####

# list.files("RawData")
df_nhanes <- read_rds(here("RawData/nhanes_gh_act.rds"))
head(df_nhanes)

N <- length(unique(df_nhanes$SEQN)) # N = 13603

# clean & format
mat_mims<- unclass(df_nhanes$MIMS)
colnames(mat_mims) <- 1:1440
## A data frame with participant ID, weekday, MIMS
df_mpa <- bind_cols(df_nhanes %>% select(SEQN, PAXDAYWM), 
                    as.data.frame(mat_mims))
## dichotomize MIMS
df_mpa <- df_mpa %>% 
  pivot_longer(3:1442, names_to = "Minute", values_to = "MIMS") %>%
  mutate(Minute = as.numeric(Minute), 
         PAXDAYWM =as.numeric(PAXDAYWM)) %>% 
  mutate(Indicator = ifelse(MIMS >= 10.588, 1, 0)) 

head(df_mpa)



#### Data spliting ####

# subsamples: N = 100 for both training and testing
N <- 100

id_tr <- sample(unique(df_mpa$SEQN), size = N)

data_tr <- df_mpa %>% filter(SEQN %in% id_tr) %>% 
  rename(id = SEQN, visit = PAXDAYWM, t=Minute, Y = Indicator)

# level-1 eigenvalues very very small
# maybe we don't need subject-level random effects

# #### gmFPCA #####
# 
# # check data format
# head(data_tr)
# 
# # gmFPCA process
# fit_time <- system.time({
#   gmfpca_fit <- gm_FPCA(data = data_tr, bin_w=10, L = 4)
# })
# 
# 
# fit_time/60/60 # about three hours for 500 training subjects
# 
# gmfpca_fit$evals1 
# Seems like the level 1 eigenvalues are very small
# I wonder if we should just use 1 subject-level PC

#### gmFPCA using PVE ####

# other parameters
t_vec <- unique(data_tr$t) # measurement grid
J <- max(as.numeric(unique(data_tr$visit))) # number of visits per subject
K <- length(t_vec) # number of observations per visit

# Step 1
print("Step 1: bin data")
## create an index frame for bins
bin_w <- 10
df_bin <- data.frame(t = t_vec) %>% 
  mutate(tind = 1:K) %>% 
  mutate(bin = ceiling(tind/bin_w))

## bin midpoints
bin_mid <- df_bin %>% group_by(bin) %>% summarise(bin_mid = median(t)) 

## data
data_tr <- data_tr %>% 
  left_join(df_bin, by = "t")

##### Step 2 #####
print("Step 2: local GLMM")
df_s2 <- data_tr %>% 
  filter(complete.cases(.)) %>% ## make it accommodate missing values
  group_by(bin) %>%
  group_map(~{
    fit_local <- glmer(Y~1+(1|id)+(1|id:visit), data = .x, 
                       family = binomial, nAGQ=0)
    eta_hat <- predict(fit_local, type = "link")
    .x$eta_hat <- eta_hat
    return(data.frame(.x))
  }, .keep = T) %>% bind_rows()

head(df_s2, 20)

##### plot check 
df_s2 %>% filter(id %in% sample(id_tr, 4)) %>% 
  mutate(eta_hat = exp(eta_hat)/(1+exp(eta_hat))) %>%
  ggplot()+
  geom_line(aes(x=t, y=eta_hat, col = "est"))+
  geom_point(aes(x=t, y=Y), size = 0.1)+
  facet_grid(row = vars(id), col = vars(visit))

##### step 3 #####
print("Step 3: mFPCA")

df_s3 <- df_s2 %>% select(id, visit, bin, eta_hat) %>% 
  distinct(.) %>% 
  pivot_wider(., names_from = bin, values_from = eta_hat) 

## number of knots
nknot <- min(nrow(bin_mid)-4, 30)
fit_mfpca <- mfpca.face(as.matrix(df_s3 %>% select(-id, -visit)),
                        id = df_s3$id, visit = df_s3$visit, 
                        argvals = as.vector(bin_mid$bin_mid),
                        pve = 0.8, knots = nknot)
### Used pve instead of npc so that the eigenvalues are not too small

fit_mfpca$evalues

L <- length(fit_mfpca$evalues$level1)
M <- length(fit_mfpca$evalues$level2)

# plot check
plot(1:144, fit_mfpca$mu)








#### step 4 ####
print("Step 4: re-evaluation")

## eigenfunctions:
reeval_efunc_l1 <- reeval_efunctions(
  knots = nrow(bin_mid)/2,
  argvals_bin = bin_mid$bin_mid,
  argvals = t_vec,
  efunctions = fit_mfpca$efunctions$level1,
  npc = ncol(fit_mfpca$efunctions$level1)
)

reeval_efunc_l2 <- reeval_efunctions(
  knots = nrow(bin_mid)/2,
  argvals_bin = bin_mid$bin_mid,
  argvals = t_vec,
  efunctions = fit_mfpca$efunctions$level2,
  npc = ncol(fit_mfpca$efunctions$level2)
)


colnames(reeval_efunc_l1) <- paste0("phi", 1:L)
colnames(reeval_efunc_l2) <- paste0("psi", 1:M)

head(reeval_efunc_l2)

## eigenvalues and mean
dim(data_tr)
data_tr <- data_tr[, 1:7]
data_tr <- data_tr %>% 
  left_join(data.frame(t=t_vec, reeval_efunc_l1), by = "t") %>%
  left_join(data.frame(t=t_vec, reeval_efunc_l2), by = "t") 



### create an ID-Visit variable
data_tr <- data_tr %>% mutate(id_visit = paste0(id, "_", visit))
data_tr$id_visit <- as.factor(data_tr$id_visit)
### global model
system.time({
    global_glmm <- bam(Y ~ s(t, bs="cc", k=10)+
                         s(id, by=phi1, bs="re", k=10)+
                         s(id, by=phi2, bs="re", k=10)+
                         s(id, by=phi3, bs="re", k=10)+
                         s(id, by=phi4, bs="re", k=10)+
                         s(id, by=phi5, bs="re", k=10)+
                         s(id_visit, by=psi1, bs="re", k=10)+
                         s(id_visit, by=psi2, bs="re", k=10)+
                         s(id_visit, by=psi3, bs="re", k=10)+
                         s(id_visit, by=psi4, bs="re", k=10)+
                         s(id_visit, by=psi5, bs="re", k=10)+
                         s(id_visit, by=psi6, bs="re", k=10)+
                         s(id_visit, by=psi7, bs="re", k=10)+
                         s(id_visit, by=psi8, bs="re", k=10)+
                         s(id_visit, by=psi9, bs="re", k=10)+
                         s(id_visit, by=psi10, bs="re", k=10)+
                         s(id_visit, by=psi11, bs="re", k=10)+
                         s(id_visit, by=psi12, bs="re", k=10)+
                         s(id_visit, by=psi13, bs="re", k=10)+
                         s(id_visit, by=psi14, bs="re", k=10)+
                         s(id_visit, by=psi15, bs="re", k=10), 
                       family = binomial, data=data_tr, 
                       method = "fREML",
                       discrete = TRUE)
})

#output
# eigenvalues (scaled)
evals <- 1/global_glmm$sp[-1]
evals1 <- evals[1:L]
evals2 <- evals[-c(1:L)]

# population mean (scaled)
mu <- predict(global_glmm, type = "terms")[1:K,1]
mu <- mu+coef(global_glmm)[1]

#### Prediction/computation ####

# select subjects with complete record as test sample
id_all <- unique(df_nhanes$SEQN)
id_te_all <- id_all[!id_all%in% id_tr]

id_te <- df_mpa %>% filter(SEQN %in% id_te_all) %>% 
  select(SEQN, PAXDAYWM) %>% distinct(.) %>%
  group_by(SEQN) %>% summarise(nday = n()) %>%
  filter(nday == 7) %>% select(SEQN) %>% as.vector()
id_te <- id_te[[1]]

data_te <- df_mpa %>% filter(SEQN %in% id_te) %>% 
  rename(id = SEQN, visit = PAXDAYWM, t=Minute, Y = Indicator)


# select 4 subjects
pred_id <- sample(id_te, 4)

##### Case 1: predict day 3 pm and evening #####

# container
df_pred1 <- data_te %>% filter(id %in% pred_id) %>%
  filter(visit <=3) %>%
  mutate(pred_2.5=NA)
pred_time <- rep(NA, 4)

Jmax <- 3

for(i in seq_along(pred_id)){
  
  ## prediction observed
  df_te_ij <- df_pred1 %>% filter(id == pred_id[i]) %>%
    # filter(visit <= 3) %>%
    filter(!(visit == 3 & t > 720))
  
  
  pred_data <- list(J = 3, K = 1440,
                    Ju = 3, Ku = table(df_te_ij$visit), 
                    Nobs = nrow(df_te_ij), Y = df_te_ij$Y,
                    L = 5, M = 15,
                    efuncs_l1 = reeval_efunc_l1,
                    efuncs_l2 = reeval_efunc_l2,
                    b0 = mu,
                    lambda = evals1,
                    gamma = evals2)
  
  # scores
  time_i <- system.time({
      fit <- stan(
        file = here("Code/prediction.stan"),  # Stan program
        data = pred_data,    # named list of data
        chains = 2,             # number of Markov chains
        warmup = 1000,          # number of warmup iterations per chain
        iter = 2000,            # total number of iterations per chain
        cores = 1,              # number of cores (could use one per chain)
        refresh = 0)             # no progress shown
      })
      
    # scores
    sum_scores <- summary(fit)$summary[, "mean"]
    xi_pred <- sum_scores[1:L]
    zeta_pred <- sum_scores[(L+1):(length(sum_scores)-1)]
    zeta_pred <- matrix(zeta_pred, M, 3, byrow = T) # M (npc) by J (day)
      
    # linear predictors
    eta_l1 <- reeval_efunc_l1 %*% xi_pred
    eta_l2 <- as.vector(reeval_efunc_l2 %*% zeta_pred)
    eta_pred <- rep(mu, Jmax)+rep(eta_l1, Jmax)+eta_l2
    
    # put in container
    df_pred1[df_pred1$id==pred_id[i], "pred_2.5"]  <- eta_pred
    pred_time[i] <- time_i[3]/60
    
    # progress 
    print(paste0("Subjct ", pred_id[i], " completed"))
    
}
  


# about 10 minutes for 1 subject
pred_time



##### Case 2: predict entire day 3 #####

# container
df_pred1 <- df_pred1 %>% mutate(pred_2 = NA)
head(df_pred1)

pred_time2 <- rep(NA, 4)
Jmax <- 3

for(i in seq_along(pred_id)){
  
  # observed
  df_te_ij <- df_pred1 %>% filter(id == pred_id[i]) %>% 
    # filter(visit <= 3) %>% 
    filter(visit < 3)
  
  pred_data <- list(J = 3, K = 1440,
                    Ju = 3, Ku = c(table(df_te_ij$visit), 0), 
                    Nobs = nrow(df_te_ij), Y = df_te_ij$Y,
                    L = 5, M = 15,
                    efuncs_l1 = reeval_efunc_l1,
                    efuncs_l2 = reeval_efunc_l2,
                    b0 = mu,
                    lambda = evals1,
                    gamma = evals2)
  
  # scores
  time_i <- system.time({
    fit <- stan(
      file = here("Code/prediction.stan"),  # Stan program
      data = pred_data,    # named list of data
      chains = 2,             # number of Markov chains
      warmup = 1000,          # number of warmup iterations per chain
      iter = 2000,            # total number of iterations per chain
      cores = 1,              # number of cores (could use one per chain)
      refresh = 0)             # no progress shown
  })
  
  # scores
  sum_scores <- summary(fit)$summary[, "mean"]
  xi_pred <- sum_scores[1:L]
  zeta_pred <- sum_scores[(L+1):(length(sum_scores)-1)]
  zeta_pred <- matrix(zeta_pred, M, 3, byrow = T) # M (npc) by J (day)
  
  # linear predictors
  eta_l1 <- reeval_efunc_l1 %*% xi_pred
  eta_l2 <- as.vector(reeval_efunc_l2 %*% zeta_pred)
  eta_pred <- rep(mu, Jmax)+rep(eta_l1, Jmax)+eta_l2
  
  # put in container
  df_pred1[df_pred1$id==pred_id[i], "pred_2"]  <- eta_pred
  pred_time2[i] <- time_i[3]/60
  
  # progress 
  print(paste0("Subjct ", pred_id[i], " completed"))
  
}



# about 10 minutes for 1 subject
pred_time2


##### Figure #####
library(RColorBrewer)
cols <- c("#66C2A5", "#FC8D62", "#000000")
names(cols) <- c("Half visit", "Full visit", "True")

df_pred1 %>%
  filter(visit==3) %>%
  mutate_at(vars(starts_with("pred_")), function(x){exp(x)/(1+exp(x))}) %>%
  mutate(pred_2.5 = ifelse(t>720, pred_2.5, NA)) %>%
  ggplot()+
  geom_point(aes(x=t, y=Y, col = "True"), size = 0.1)+
  geom_line(aes(x=t, y=pred_2.5, col = "Half visit"))+
  geom_line(aes(x=t, y=pred_2, col = "Full visit"))+
  facet_wrap(~id, nrow = 1)+
  scale_color_manual(values = cols)+
  labs(col = "", x = "", y = "Activity")+
  scale_x_continuous(breaks = seq(0,1440, by = 360), 
                     labels = c("0am", "6am", "12pm", "6pm", "12am"))+
  theme_minimal()
ggsave(file=here("Images/case_study_pred.png"), 
       height = 3, width = 12, bg = "white")

##### Save #####

save(df_pred1, file = here("Data/TempCaseStudyPred.RData"))

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
                        efuncs_l1 = reeval_efunc_l1,
                        efuncs_l2 = gre
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
pred_list_allM[[100]] %>% #head()
# pred_df_m %>%
  filter(id %in% sample(101:200, 4)) %>% #head()
  # filter(id==41) %>%
  mutate_at(vars(starts_with("eta_pred")), function(x)(exp(x)/(1+exp(x)))) %>%
  ggplot()+
  geom_line(aes(x=t, y=probs, col = "True"))+
  geom_line(aes(x=t, y=eta_pred_J1, col = "J1"), linetype = "dashed")+
  geom_line(aes(x=t, y=eta_pred_J2, col = "J2"), linetype = "dashed")+
  geom_line(aes(x=t, y=eta_pred_J3, col = "J3"), linetype = "dashed")+
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
  ## technically, this is not nesting. But the real nesting one is not fittable
  # with spline basis? 
  glmm_fit <- mixed_model(fixed = Y ~ ns(t, df = 2),
                          random = ~ ns(t, df = 2)  | id_visit,
                          data = data_tr,
                    family = binomial())
  # linear model
  # glmm_fit <- mixed_model(fixed = Y ~ t,
  #                         random = ~ t | id_visit,
  #                         data = data_tr,
  #                         family = binomial())
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


