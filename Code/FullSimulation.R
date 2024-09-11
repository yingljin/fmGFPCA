library(here)

source(here("Code/Functions/Simulation.R"))
source(here("Code/Functions/gmFPCA_func.R"))
source(here("Code/Functions/ReEvalEfuncs.R"))
source(here("Code/Functions/Predict.R"))

set.seed(911)

#### generate data ####

data <- gen_bi(N = 100, J = 3, K = 100, t_vec = seq(0, 1, length.out = 100),
       L = 4, M = 4, lambda = 0.5^((1:4)-1), gamma =  0.5^((1:4)-1))

data$l2score
df_train <- data$data %>% mutate(id = factor(id), visit = factor(visit))
head(df_train)

# scores
apply(data$l1score, 2, var)
apply(data$l2score, MARGIN = c(2,3), var)

# data
df_train %>% filter(id %in% 1:4) %>%
  ggplot()+
  geom_line(aes(x=t, y=probs))+
  geom_point(aes(x=t, y=Y), size = 0.5)+
  facet_grid(row = vars(id), col=vars(visit))

#### GMFPCA ####
gmfpca_fit <- gm_FPCA(data = df_train, bin_w=10, L = 4)


gmfpca_fit$evals1
gmfpca_fit$evals2
plot(gmfpca_fit$mu)


#### Prediction ####

df_test <- gen_bi(N = 1, J = 3, K = 100, t_vec = seq(0, 1, length.out = 100),
       L = 4, M = 4, lambda = 0.5^((1:4)-1), gamma =  0.5^((1:4)-1))

# truncate
df_test_tr <- df_test$data %>% filter(t <= 0.5)

pred <- mDynPred(data = df_test_tr, 
         mu = gmfpca_fit$mu,
         evals1 = gmfpca_fit$evals1,
         evals2 = gmfpca_fit$evals2,
         efuncs1 =  gmfpca_fit$efuncs_l1,
         efuncs2 = gmfpca_fit$efuncs_l2,
         J = 3
         )

pred$score1
df_test$l1score

pred$score2
df_test$l2score

head(df_test_tr)
head(pred$eta_pred)


df_pred <- df_test$data %>% mutate(eta_pred = pred$eta_pred$eta_pred,
                                   eta_pred_lb = pred$eta_pred$eta_pred_lb,
                                   eta_pred_ub = pred$eta_pred$eta_pred_ub)
df_pred %>% #head()
  mutate_at(vars(starts_with("eta")), function(x)(exp(x)/(1+exp(x)))) %>%
  ggplot()+
  geom_line(aes(x=t, y=eta, col = "True"))+
  geom_line(aes(x=t, y=eta_pred, col = "Pred"))+
  geom_line(aes(x=t, y=eta_pred_ub, col = "Pred"), linetype = "dashed")+
  geom_line(aes(x=t, y=eta_pred_lb, col = "Pred"), linetype = "dashed")+
  facet_wrap(~visit)



pred$eta_pred
pred$score1
pred$score2
