# writes function for the model fitting part of gmFPCA algorithm
# but for case study along
# allow for different number of PCs in both levels

# library(tidyverse)
# library(lme4)
# library(mgcv)
# library(refund)
# library(fastGFPCA)

# data: training data frame. must include colums id, t, visit, Y
#       id and visit need to be factor
#       visit need to be integer value
# bin_w: how many observations are included in each bin
# K = number of visits per subject
# L: number of l1, l2 eigenfunctions

# Note: this functions can only deal with regularly measured data
#       I plan to modify this to accomodate irregular data as well

gm_FPCA <- function(data, bin_w=10, L=0, M = 4, scores=F){
  
  
  # other parameters
  t_vec <- unique(data$t) # measurement grid
  J <- max(as.numeric(unique(data$visit))) # number of visits per subject
  K <- length(t_vec) # number of observations per visit
  
  # add a function index
  # data$ind <- rep(1:K, J)
  
  # Step 1
  print("Step 1: bin data")
  ## create an index frame for bins
  df_bin <- data.frame(t = t_vec) %>% 
    mutate(tind = 1:K) %>% 
    mutate(bin = ceiling(tind/bin_w))
  
  ## bin midpoints
  bin_mid <- df_bin %>% group_by(bin) %>% summarise(bin_mid = median(t)) 
  
  ## data
  data <- data %>% 
    left_join(df_bin, by = "t")
  
  # Step 2
  print("Step 2: local GLMM")
  df_s2 <- data %>% 
    filter(complete.cases(.)) %>% ## make it accommodate missing values
    group_by(bin) %>%
    group_map(~{
      fit_local <- glmer(Y~1+(1|id)+(1|id:visit), data = .x, 
                         family = binomial, nAGQ=0)
      eta_hat <- predict(fit_local, type = "link")
      .x$eta_hat <- eta_hat
      return(data.frame(.x))
    }, .keep = T) %>% bind_rows()
  
  # step 3
  print("Step 3: mFPCA")
  
  df_s3 <- df_s2 %>% select(id, visit, bin, eta_hat) %>% 
    distinct(.) %>% 
    pivot_wider(., names_from = bin, values_from = eta_hat) 
  
  ## number of knots
  nknot <- min(nrow(bin_mid)-4, 30)
  fit_mfpca <- mfpca.face(as.matrix(df_s3 %>% select(-id, -visit)),
                          id = df_s3$id, visit = df_s3$visit, 
                          argvals = as.vector(bin_mid$bin_mid),
                          npc = L, knots = nknot)
  
  # step 4
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
  colnames(reeval_efunc_l2) <- paste0("psi", 1:L)
  
  ## eigenvalues and mean
  data <- data %>% 
    left_join(data.frame(t=t_vec, reeval_efunc_l1), by = "t") %>%
    left_join(data.frame(t=t_vec, reeval_efunc_l2), by = "t") 
  ### create an ID-Visit variable
  data <- data %>% mutate(id_visit = paste0(id, "_", visit))
  data$id_visit <- as.factor(data$id_visit)
  ### global model
  global_glmm <- bam(Y ~ s(t, bs="cc", k=nknot)+
                       s(id, by=phi1, bs="re")+
                       s(id, by=phi2, bs="re")+
                       s(id, by=phi3, bs="re")+
                       s(id, by=phi4, bs="re")+
                       s(id_visit, by=psi1, bs="re")+
                       s(id_visit, by=psi2, bs="re")+
                       s(id_visit, by=psi3, bs="re")+
                       s(id_visit, by=psi4, bs="re"), 
                     family = binomial, data=data, 
                     method = "fREML",
                     discrete = TRUE)
  
  #output
  # eigenvalues (scaled)
  evals <- 1/global_glmm$sp[-1]
  evals1 <- evals[1:L]
  evals2 <- evals[-c(1:L)]
  
  # population mean (scaled)
  mu <- predict(global_glmm, type = "terms")[1:K,1]
  mu <- mu+coef(global_glmm)[1]
  
  # output
  if(scores==F){
    outlist <- list(evals1 = evals1, evals2 = evals2,
                    mu = mu, efuncs_l1 = reeval_efunc_l1,
                    efuncs_l2 = reeval_efunc_l2)
  }
  else{
    outlist <- list(evals1 = evals1, evals2 = evals2,
                    mu = mu, efuncs_l1 = reeval_efunc_l1,
                    efuncs_l2 = reeval_efunc_l2,
                    scores_l1 = fit_mfpca$scores$level1,
                    scores_l2 = fit_mfpca$scores$level2)
  }
  
  return(outlist)
  
}

