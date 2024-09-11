
# This code simulate multi-level binary functions



# set.seed(625)
# 
# library(tidyverse)
# library(here)

#### Function ####

# N: sample size
# J: number of visit per subject
# K: number of observations per visit
# t_vec: observation grid
# L: number of level 1 PC
# M: number of level 2 PC
# lambda: level 1 eigenvalues
# gamma: level 2 eigenvalues

gen_bi <- function(N, J, K, t_vec, L, M, lambda, gamma){
  
  # eigenfunctions
  ## Level 1
  phi_mat <- matrix(NA, nrow = K, ncol = L)
  phi_mat[, 1] <- sqrt(2)*sin(2*pi*t_vec)
  phi_mat[, 2] <- sqrt(2)*cos(2*pi*t_vec)
  phi_mat[, 3] <- sqrt(2)*sin(4*pi*t_vec)
  phi_mat[, 4] <- sqrt(2)*cos(4*pi*t_vec)
  colnames(phi_mat) <- paste0("phi", 1:L)
  
  
  # level 2 eigenfunctions:
  psi_mat <- matrix(NA, nrow = K, ncol = M)
  psi_mat[, 1] <- sqrt(2)*sin(6*pi*t_vec)
  psi_mat[, 2] <- sqrt(2)*cos(6*pi*t_vec)
  psi_mat[, 3] <- sqrt(2)*sin(8*pi*t_vec)
  psi_mat[, 4] <- sqrt(2)*cos(8*pi*t_vec)
  colnames(psi_mat) <- paste0("psi", 1:M)
  
  df_full <- expand.grid(id=1:N, visit = 1:J, t=t_vec) %>% arrange(id, visit, t)
  df_full$eta <- NA
  
  # containers for score: 
  xi_true <- matrix(NA, nrow = N, ncol = L) # subject by eigenfunction
  zeta_true <- array(NA, dim = c(N, J, M)) # subject by visit by eigenfunction
  
  # generate latent function
  for(i in 1:N){
    # level 1 score (subject): L by 1
    xi <- rnorm(L, 0, sqrt(lambda))
    xi_true[i, ] <- xi
    # level 2 score (subject-visit level): J by M
    # zeta <- sapply(gamma, function(x){rnorm(J, mean=0, sd = sqrt(x))})
    zeta <- rnorm(J*M, mean = 0, sd = rep(sqrt(gamma), J))
    zeta <- matrix(zeta, nrow = J, ncol=M, byrow = T) # visit by eigenfunction
    zeta_true[i,,] <- zeta
    # latent function
    ## subjectlevel
    eta_l1 <- phi_mat %*% xi
    eta_l1 <- rep(eta_l1, J)
    ## subject-visit level 
    eta_l2 <- psi_mat %*% t(zeta)
    eta_l2 <- as.vector(eta_l2)
    ## overall
    eta_i <- 0 + eta_l1 + eta_l2
    df_full[df_full$id==i, "eta"] <-  eta_i
  }
  
  # generate binary outcome
  df_full <- df_full %>% mutate(probs = exp(eta)/(1+exp(eta))) 
  df_full$Y <- rbinom(nrow(df_full), 1, df_full$probs)
  
  return(list(data = df_full, l1score = xi_true, l2score = zeta_true))
  
}

#### Test ####

# df1 <- gen_bi(N = 20, J = 2, K = 100, t_vec = seq(0, 1, length.out = 100), 
#               L = 4, M = 4, lambda = 0.5^((1:4)-1), gamma =  0.5^((1:4)-1))
# 
# df1$data



# #### Check data ####
# # scores
# apply(df1$l1score, 2, var)
# apply(df1$l2score, MARGIN = c(2,3), var)
# 
# # data
# df1$data %>% filter(id %in% 1:4) %>% 
#   ggplot()+
#   geom_line(aes(x=t, y=probs))+
#   geom_point(aes(x=t, y=Y), size = 0.5)+
#   facet_grid(row = vars(id), col=vars(visit))
# # Since the density of observation is much higher
# # is very hard to see the change of activity with the flunctuations of eita
# 
# df_full %>% group_by(id, visit) %>% select(id, visit, Y) %>%
#   summarise(freq = sum(Y))
# 
# #### save data ####
# save(df_full, xi_true, zeta_true, file = here("Data/simData_N500_K500.RData"))
