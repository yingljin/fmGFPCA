
# This code simulate multi-level binary functions



set.seed(625)

library(tidyverse)
library(here)


#### Simulation set up ####

N <- 500 # sample size 
J <- 5 # number of visit per subject
K <- 500 # number of observations per visit
t_vec <- seq(0, 1, length.out = K) # time grid
L <- 4 # number of level 1 eigenfunctions
M <- 4 # number of level 2 eigenfunctions

#### Eigenfunctions ####

# level 1 eigenfunctions:
phi_mat <- matrix(NA, nrow = K, ncol = L)
phi_mat[, 1] <- sqrt(2)*sin(2*pi*t_vec)
phi_mat[, 2] <- sqrt(2)*cos(2*pi*t_vec)
phi_mat[, 3] <- sqrt(2)*sin(4*pi*t_vec)
phi_mat[, 4] <- sqrt(2)*cos(4*pi*t_vec)

colnames(phi_mat) <- paste0("phi", 1:L)
head(phi_mat)

phi_mat %>% data.frame() %>%
  mutate(t=t_vec) %>%
  pivot_longer(starts_with("phi")) %>%
  ggplot()+
  geom_line(aes(x=t, y=value))+
  facet_wrap(~name, nrow= 1)+
  labs(x="t", y="")

# level 1 eigenvalues
lambda <- 0.5^((1:L)-1)

# level 2 eigenfunctions:
psi_mat <- matrix(NA, nrow = K, ncol = M)
psi_mat[, 1] <- sqrt(2)*sin(6*pi*t_vec)
psi_mat[, 2] <- sqrt(2)*cos(6*pi*t_vec)
psi_mat[, 3] <- sqrt(2)*sin(8*pi*t_vec)
psi_mat[, 4] <- sqrt(2)*cos(8*pi*t_vec)

colnames(psi_mat) <- paste0("psi", 1:M)
head(psi_mat)

psi_mat %>% data.frame() %>%
  mutate(t=t_vec) %>%
  pivot_longer(starts_with("psi")) %>%
  ggplot()+
  geom_line(aes(x=t, y=value))+
  facet_wrap(~name, nrow = 1)+
  labs(x="t", y="")

# level 2 eigenvalues
gamma <- 0.5^((1:M)-1)


#### Generate a full datasets ####

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

#### Check data ####
# scores
apply(xi_true, 2, var)
apply(zeta_true, MARGIN = c(2,3), var)

# data
df_full %>% filter(id %in% 1:4) %>%
  ggplot()+
  geom_line(aes(x=t, y=probs))+
  geom_point(aes(x=t, y=Y), size = 0.5)+
  facet_grid(row = vars(id), col=vars(visit))
# Since the density of observation is much higher
# is very hard to see the change of activity with the flunctuations of eita

df_full %>% group_by(id, visit) %>% select(id, visit, Y) %>%
  summarise(freq = sum(Y))

#### save data ####
save(df_full, xi_true, zeta_true, file = here("Data/simData_N500_K500.RData"))
