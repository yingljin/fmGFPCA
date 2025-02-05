
# This code script generates multi-level generalized functional data
# that will be used for the simulation portion in the manuscript


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


set.seed(915)
# devtools::install_github("julia-wrobel/fastGFPCA", build_vignettes = TRUE)

# simulation scheme
# N_train = 100,  N_test = 100, J = 2, K = 200
# observations stops at the the midpoint of each visit
# number of FPC = 4

M  <- 100
Ntr <- 100
Nte <- 100

# J <- 2 # total number of visits
J <- 3 # total number of visits

K <- 200

#### Generate data ####

# containers
data_list <- list()


# generation

pb = txtProgressBar(min = 0, max = M, initial = 0, style = 3) 

for(m in 1:M){
  # generate training data
  data_all <- gen_bi(N = Ntr+Nte, J = J, K = K, 
                     t_vec = seq(0, 1, length.out = K),
                     L = 4, M = 4, 
                     lambda = 0.5^((1:4)-1), 
                     gamma =  0.5^((1:4)-1))
  data_all$data$id <- as.factor(data_all$data$id)
  data_all$data$visit <- as.factor(data_all$data$visit)
  data_list[[m]] <- data_all$data
  
  setTxtProgressBar(pb, m)
}

close(pb)

# check

data_list[[1]] %>% filter(id %in% sample(1:200, 2)) %>% 
  ggplot()+
  geom_line(aes(x=t, y=probs))+
  geom_point(aes(x=t, y=Y), size = 0.5)+
  facet_grid(row = vars(id), col=vars(visit))

save(data_list,
     file = here("Data/SimData/SimDataJ3.RData")) # data that will be used for all of simulation
