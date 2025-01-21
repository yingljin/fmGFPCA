# This script generates data for simulation case 2
# where missing is induced

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

set.seed(121)

M <- 100 # number of simulations
N <- 100 # sample size
J <- 2 # total number of visits
K <- 200


# 10% of the subjects have missing outcomes
# each subject has a chunk of missing from 1 to 10 hours
p1 <- 0.1 

# containers
data_list_mis <- list()

# generate data
pb <- txtProgressBar(min = 0, max = M, style = 3)

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
  data_m <- data_all$data
  miss_id <- sample(1:N, size = p1*N, replace = FALSE)
  ## add a time index
  data_m <- data_m %>%
    mutate(Yobs = Y, tind = rep(1:(J*K), N))
  ## for each individual with missing, generate a random missing chunk
  for(i in seq_along(miss_id)){
    start_mis <- sample(1:(J*K), 1) # start of missing chunk
    len_mis <- sample(60:300, size = 1) # a random chuck with length 1-5 hours
    end_mis <- min(start_mis+len_mis, J*K)
    data_m[data_m$id==miss_id[i]&data_m$tind>=start_mis&data_m$tind<=end_mis, "Yobs"] <- NA
  }
  
  data_list_mis[[m]] <- data_m
  setTxtProgressBar(pb, m)
  
}

close(pb)

# check
data_list_mis[[1]] %>% #head()
  select(id, Yobs, tind) %>% 
  pivot_wider(values_from = Yobs, names_from = tind) %>%
  vis_miss()

save(data_list_mis, 
     file = here("Data/SimData/SimDataMis.RData"))
