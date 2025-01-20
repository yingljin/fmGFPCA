# Figures out how to make "corrected" out of sample prediction with 
# GLMMadaptive


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

set.seed(116)

#### Generate data ####

Ntr <- 100
Nte <- 100

J <- 2 # total number of visits
K <- 200

# generation
data_all <- gen_bi(N = Ntr+Nte, J = J, K = K, 
                   t_vec = seq(0, 1, length.out = K),
                   L = 4, M = 4, 
                   lambda = 0.5^((1:4)-1), 
                   gamma =  0.5^((1:4)-1))
data_all$data$id <- as.factor(data_all$data$id)
data_all$data$visit <- as.factor(data_all$data$visit)

# split 
data_tr <- data_all$data %>% filter(id %in% 1:Ntr) %>% 
  mutate(id = droplevels(id))
## test
data_te <- data_all$data %>% filter(!id %in% 1:Ntr) %>%
  mutate(id = droplevels(id), visit=droplevels(visit))

#### GLMMadaptive ####

# fit a nexted model
# about 102 seconds
# can only fit intercept-only nested model
system.time({
  glmm_fit <- mixed_model(fixed = Y ~ 1,
                          random = ~ 1 | id/visit,
                          data = data_tr,
                          family = binomial())
})

# on the test data
# try to get individual-level prediction by term
# first, truncate data at half of second visit
data_te1 <- data_te %>% filter(!(visit == 2 & t > 0.5))
data_te2 <- data_te %>% filter(visit == 2 & t > 0.5)

pred <- predict(glmm_fit, newdata = data_te1, newdata2 = data_te2,
                type = "subject_specific")


#### Try BAM ####
data_tr$id_visit <- paste0(data_tr$id, "_", data_tr$visit)
data_tr$id_visit <- as.factor(data_tr$id_visit)

# fit model
# about 18 seconds
system.time(
  fit_bam <- bam(Y ~ s(t) + s(id, by = t, bs = "re") + s(visit, id, by = t, bs = "re"),
                 data = data_tr,
                 family = binomial)
)

## coefficients
length(fit_bam$coefficients)
fit_bam$coefficients[211:310]

## fitted 

fitted <- predict(fit_bam, type = "link")
head(fitted)

data_tr %>%
  mutate(fitted = fit_bam$fitted) %>% 
  filter(id %in% 1:2) %>% 
  ggplot()+
  geom_line(aes(x=t, y=probs))+
  geom_point(aes(x=t, y=Y), size = 0.5)+
  geom_line(aes(x=t, y = fitted, col = "fitted"), linetype = "dashed")+
  facet_grid(row = vars(id), col=vars(visit))

data_tr %>%
  mutate(fitted = fitted) %>% 
  mutate(fitted = exp(fitted)/(1+exp(fitted))) %>%
  filter(id == 1) %>% 
  select(visit, t, fitted) %>% 
  pivot_wider(names_from = visit, values_from = fitted) %>% 
  rename("fit1" = "1", "fit2" = "2") %>%
  filter(abs(fit1-fit2) > 0.001)

# why are the fitted tracks from different visits within the same subject 
# stay the same?
# Well actually the fitted values are different but the difference is too small (<=0.001)
# after adding a random slope the difference starts to show on the linear scale 
# but on the probability scale the difference is still very small
# The model did not pick up within-subject variation? 

# test set

pred <- predict(fit_bam, newdata = data_te)

