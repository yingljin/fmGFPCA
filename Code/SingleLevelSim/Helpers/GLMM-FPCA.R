##### Local GLMMs #####

# get predicted latent function
# df must have: outcome Y, subject id

pred_latent <- function(df, n_node){
    this_glm <- glmer(Y ~ 1 + (1|id), data = df, family = binomial, nAGQ=n_node)
    eta_hat <- predict(this_glm, type = "link")
    df$eta_hat <- eta_hat
    return(df)
}
  
# Additionally, get the average estimation across subject 
# by extacting the random intercepts from each bin
# mean_latent <- function(df, n_node){
#     this_glm <- glmer(Y ~ 1 + (1|id), data = df, family = binomial, nAGQ=n_node)
#     beta_hat <- coef(summary(this_glm))[1]
#     return(beta_hat)
# }
#   




