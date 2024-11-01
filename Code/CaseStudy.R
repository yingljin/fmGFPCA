
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


##### Case 2: predict day 3 evening #####

# container
df_pred1 <- df_pred1 %>% mutate(pred_2.75 = NA)
head(df_pred1)

pred_time2 <- rep(NA, 4)
Jmax <- 3

for(i in seq_along(pred_id)){
  
  # observed
  df_te_ij <- df_pred1 %>% filter(id == pred_id[i]) %>% 
    # filter(visit <= 3) %>% 
    filter(!(visit == 3 & t > 1020))
  
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
names(cols) <- c("12pm-12am", "5pm-12am", "True")

df_pred1 %>%
  filter(visit==3) %>%
  mutate_at(vars(starts_with("pred_")), function(x){exp(x)/(1+exp(x))}) %>%
  mutate(pred_2.5 = ifelse(t>720, pred_2.5, NA),
         pred_2 = ifelse(t>1020, pred_2, NA)) %>%
  ggplot()+
  geom_point(aes(x=t, y=Y, col = "True"), size = 0.1)+
  geom_line(aes(x=t, y=pred_2.5, col = "12pm-12am"))+
  geom_line(aes(x=t, y=pred_2, col = "5pm-12am"))+
  facet_wrap(~id, nrow = 1)+
  scale_color_manual(values = cols)+
  labs(col = "Prediction window", x = "", y = "Activity")+
  scale_x_continuous(breaks = seq(0,1440, by = 360), 
                     labels = c("0am", "6am", "12pm", "6pm", "12am"))+
  theme_minimal()+
  theme(legend.position = "bottom")
ggsave(file=here("Images/case_study_pred.png"), 
       height = 4, width = 12, bg = "white")

#### GLMMadpative ####

##### Model fit #####

# fit model
data_tr2 <- data_tr[, 1:7] %>% mutate(id_visit = paste0(id, "_" , visit))

system.time({
  fit_glmmad <- mixed_model(fixed = Y ~ t,
                            random = ~ t| id_visit,
                            data = data_tr2,
                            family = binomial())
})

#####  predict Case 1: predict day 3 pm and evening #####

# container
df_pred_given <- data_te %>% filter(id %in% pred_id) %>%
  filter(visit <=3) %>%
  filter(!(visit == 3 & t > 720)) %>%
  mutate(id_visit = paste0(id, "_" , visit))
  

df_pred_new <- data_te %>% filter(id %in% pred_id) %>%
  filter(visit <=3) %>%
  mutate(id_visit = paste0(id, "_" , visit))

unique(df_pred_given$id)
unique(df_pred_new$id)

# pred_time2 <- rep(NA, 4)

pred_2.5 <- predict(fit_glmmad, newdata = df_pred_given, 
                    newdata2 = df_pred_new,
                    type_pred = "link", type = "subject_specific",
                    return_newdata = T)$newdata2

#####  predict Case 1: predict day 3 all day #####

# container
df_pred_given <- data_te %>% filter(id %in% pred_id) %>%
  filter(visit <=3) %>%
  filter(!(visit == 3 & t > 1020)) %>%
  mutate(id_visit = paste0(id, "_" , visit))


unique(df_pred_given$id)
unique(df_pred_new$id)

pred_2 <- predict(fit_glmmad, newdata = df_pred_given, 
                    newdata2 = df_pred_new,
                    type_pred = "link", type = "subject_specific",
                    return_newdata = T)$newdata2
pred_2

