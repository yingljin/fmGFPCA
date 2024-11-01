
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

# subsamples: N = 100 for model training
N <- 100

id_tr <- sample(unique(df_mpa$SEQN), size = N)

data_tr <- df_mpa %>% filter(SEQN %in% id_tr) %>% 
  rename(id = SEQN, visit = PAXDAYWM, t=Minute, Y = Indicator)
head(data_tr)
data_tr$Yobs <- data_tr$Y

#### Introducing missing ####

# plug in missing values of 4 subjects
# each missing for 3-5 hours

id_mis <- sample(id_tr, 4)
id_mis

df_mis <- data.frame(id = id_mis, visit = NA, start = NA, end = NA)

for(i in seq_along(id_mis)){
  
  # visit with missing
  day_mis <- sample(unique(data_tr$visit[data_tr$id==id_mis[i]]), size = 1)
  
  # time the missing period start
  time_mis <- sample(180:300, size = 1)
  
  # missing period
  start_mis <- sample(1:(1440-time_mis-1), size = 1)
  end_mis <- start_mis + time_mis
  
  data_tr$Yobs[data_tr$id==id_mis[i] & data_tr$visit == day_mis  & data_tr$t > start_mis &  data_tr$t < end_mis ] <- NA
  
  
  df_mis$visit[i] <- day_mis
  df_mis$start[i] <- start_mis
  df_mis$end[i] <- end_mis
}

# check


data_tr %>% 
  select(id, visit, t, Yobs) %>%
  pivot_wider(values_from = Yobs, names_from = t) %>%
  visdat::vis_miss()

df_mis
data_tr %>% 
  select(id, visit, t, Yobs) %>%
  group_by(id, visit) %>%
  summarise(nmis = sum(is.na(Yobs))) %>%
  filter(nmis > 0)

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
    fit_local <- glmer(Yobs~1+(1|id)+(1|id:visit), data = .x, 
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

##### step 4 #####
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
                         s(id_visit, by=psi14, bs="re", k=10), 
                       family = binomial, data=data_tr, 
                       method = "fREML",
                       discrete = TRUE)
})

#### Computation ####
# sanity check
df_mis

data_tr %>%
  group_by(id, visit) %>%
  summarize(nmis = sum(is.na(Yobs))) %>%
  filter(nmis > 0)

# prediction
df_pred <- data_tr %>%
  filter((id %in% df_mis$id[1] & visit %in% df_mis$visit[1])|
         (id %in% df_mis$id[2] & visit %in% df_mis$visit[2])|
         (id %in% df_mis$id[3] & visit %in% df_mis$visit[3])|
         (id %in% df_mis$id[4] & visit %in% df_mis$visit[4]))
pred_eta <- predict(global_glmm, newdata = df_pred)

df_pred <- df_pred %>% 
  select(id, t, visit, Y, Yobs) %>%
  mutate(eta_comp = pred_eta)

df_pred <- df_pred %>% left_join(df_mis)

#### plot ####
library(RColorBrewer)
brewer.pal(3, "Set2")
cols <- c("#FC8D62","#66C2A5", "#000000")
names(cols) <- c("Computed", "True", "Observed")
df_pred %>%
  mutate(eta_comp=exp(eta_comp)/(1+exp(eta_comp))) %>%
  mutate(eta_comp = ifelse(is.na(Yobs), eta_comp, NA)) %>%
  mutate(color = ifelse(is.na(Yobs), "True", "Observed")) %>%
  ggplot()+
  geom_point(aes(x=t, y=Y, col = color), size = 0.2)+
  geom_line(aes(x=t, y=eta_comp, col = "Computed"))+
  facet_wrap(~id, nrow = 1)+
  scale_color_manual(values = cols)+
  scale_x_continuous(breaks = seq(0,1440, by = 360), 
                     labels = c("0am", "6am", "12pm", "6pm", "12am"))+
  theme_minimal()+
  labs(x="Minute", col=" ")
ggsave(filename = "Images/case_study_miss.png", 
       height = 3, width = 12, bg = "white")


data_tr %>% 
  group_by(id, visit) %>%
  summarize(nmis = sum(is.na(Y))) %>%
  filter(nmis > 0)

data_tr %>% 
  select(id, visit, t, Y) %>% 
  pivot_wider(values_from = Y, names_from = t) %>%
  visdat::vis_miss()
  

data_pred <- data_tr %>%
  filter(id %in% id_mis) 

#### Fill missing value with GLMMadaptive ####


# data
data_tr2 <- data_tr[, 1:6]
data_tr2$id_visit <- paste0(data_tr2$id, "_", data_tr2$visit)

# fit model
system.time({
  fit_glmmad <- mixed_model(fixed = Yobs ~ ns(t, df = 2),
                                      random = ~ ns(t, df = 2) | id_visit,
                                      data = data_tr2,
                                      family = binomial())
})


# predict
df_pred2 <- df_pred[, 1:5]
pred_eta2 <- predict(fit_glmmad, newdata = df_pred2)

df_pred$eta_comp2 <- pred_eta2


#### save outcome ####
save(df_pred, global_glmm,
     file = here("Data/CaseStudyMiss.RData"))

df_pred <- df_pred %>% 
  select(id, t, visit, Y, Yobs) %>%
  mutate(eta_comp = pred_eta)

df_pred <- df_pred %>% left_join(df_mis)


#### plot ####

brewer.pal(3, "Set2")
cols <- c(brewer.pal(3, "Set2"), "#000000")
names(cols) <- c("gmFPCA", "GLMMadaptive", "True", "Observed")

df_pred %>%
  mutate(eta_comp=exp(eta_comp)/(1+exp(eta_comp)),
         eta_comp2=exp(eta_comp2)/(1+exp(eta_comp2))) %>%
  mutate(eta_comp = ifelse(is.na(Yobs), eta_comp, NA),
         eta_comp2 = ifelse(is.na(Yobs), eta_comp2, NA)) %>%
  mutate(color = ifelse(is.na(Yobs), "True", "Observed")) %>%
  ggplot()+
  geom_point(aes(x=t, y=Y, col = color), size = 0.2)+
  geom_line(aes(x=t, y=eta_comp, col = "gmFPCA"))+
  geom_line(aes(x=t, y=eta_comp2, col = "GLMMadaptive"))+
  facet_wrap(~id, nrow = 1)+
  scale_color_manual(values = cols)+
  scale_x_continuous(breaks = seq(0,1440, by = 360), 
                     labels = c("0am", "6am", "12pm", "6pm", "12am"))+
  theme_minimal()+
  theme(legend.position="bottom")+
  labs(x="Minute", col=" ")
ggsave(filename = "Images/case_study_miss.png", 
       height = 4, width = 12, bg = "white")
