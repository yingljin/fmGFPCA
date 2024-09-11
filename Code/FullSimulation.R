library(here)

source(here("Code/Functions/Simulation.R"))
source(here("Code/Functions/gmFPCA_func.R"))
source(here("Code/Functions/ReEvalEfuncs.R"))
source(here("Code/Functions/Predict.R"))

#### generate data ####

data <- gen_bi(N = 100, J = 3, K = 100, t_vec = seq(0, 1, length.out = 100),
       L = 4, M = 4, lambda = 0.5^((1:4)-1), gamma =  0.5^((1:4)-1))


df_train <- data$data
head(df_train)

#### GMFPCA ####
gmfpca_fit <- gm_FPCA(df_train, 10, L = 4)
gmfpca_fit$evals1
gmfpca_fit$evals2
gmfpca_fit$mu


#### Prediction ####

df_test <- gen_bi(N = 1, J = 3, K = 100, t_vec = seq(0, 1, length.out = 100),
       L = 4, M = 4, lambda = 0.5^((1:4)-1), gamma =  0.5^((1:4)-1))

df_test <- df_test$data %>% filter(t <= 0.5)

pred <- mDynPred(data = df_test, 
         mu = gmfpca_fit$mu,
         evals1 = gmfpca_fit$evals1,
         evals2 = gmfpca_fit$evals2,
         efuncs1 =  gmfpca_fit$efuncs_l1,
         efuncs2 = gmfpca_fit$efuncs_l2,
         J = 3
         )
pred$eta_pred
pred$score1
pred$score2
