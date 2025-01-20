#' This function uses the output from fast_fgpca to do out-of-sample dynamic prediction of binary outocmes.
#' For a single subject with a partially observed function,
#' this function makes prediction of the latent process of the unobserved part.
#'
#' @author Andrew Leroux \email{andrew.leroux@@cuanschutz.edu},
#' Julia Wrobel \email{julia.wrobel@@cuanschutz.edu}
#' Ying jin \email{ying.jin@cuanschutz.edu}
#' @import dplyr
#' @importFrom rstan stan
#' @importFrom stats coef predict binomial lm median as.formula quantile
#' @importFrom refund fpca.face
#' @importFrom lme4 glmer
#' @importFrom utils txtProgressBar setTxtProgressBar data
#' @importFrom mgcv bam predict.bam
#'
#' @return A list containing:
#' \item{pred}{A dataframe of predicted latent functions. If interval = TRUE, the dataframe also includes prediction interval.}
#' \item{scores}{Estimated FPC scores.}
#' \item{mu}{Estimated population-level mean.}
#' @references Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C. (2016).
#' Fast covariance estimation for high-dimensional functional data.
#' \emph{Statistics and Computing}, 26, 409-421.
#' DOI: 10.1007/s11222-014-9485-x.
#'
#' @examples
#' # simulate training data
#' set.seed(1001)
#'
#' # binomial data, with bins that do not overlap
#' df_gfpca <- sim_gfpca(N = 200, J = 200, case = 2)$df_gfpca
#' gfpca_mod <- fast_gfpca(df_gfpca, overlap = FALSE, binwidth = 10, family = binomial)
#'
#' # simulated test data
#' df_test <- sim_gfpca(N = 1, J = 200, case = 2)$df_gfpca
#'
#' # truncate function track
#' df_test <- df_test %>% filter(index <= 0.75)
#'
#' # prediction
#' df_pred <- DynPred(data_i = df_test, efuncs = gfpca_mod$efunctions,
#'                    evals = gfpca_mod$evalues, mu = gfpca_mod$mu)
#' df_pred$scores
#' df_pred$pred
#'
#' @param data_i Dataframe of partially observed function of one single subject. Should have variables id, value, index.
#' @param efuncs Matrix of estimated FPC basis functions. Can be obtained from fast_gfpca.
#' @param evals Estimated variance of the FPC scores. Can be obtained from fast_gfpca.
#' @param mu Estimated population-level mean function. Can be obtained from fast_gfpca.
#' @param interval Whether or not to return the prediction interval of the latent function.
#' @param int_qt A vector with two elements, including the upper and lower confidence bound of prediction interval.
#'
#' @param ... Additional arguments passed to or from other functions
#'@export



DynPred <- function(data_i, efuncs, evals, mu,
                    interval = T, int_qt = c(0.025, 0.075),
                    ...){

  tmax <- max(data_i$index)
  J <- nrow(efuncs)
  K <- ncol(efuncs)

  # stan data
  stanData <- list(
      J = J, # full measurement grid
      Ju = nrow(data_i), # number of observations
      Y = data_i$value, # outcome
      K = K, # number of PC functions
      efuncs = efuncs, # estimated PC functions
      b0 = mu, # estimated mean
      lambda = evals # estimated eigenvalues
  )

  fit <- stan(
    file = "Code/prediction.stan",  # Stan program
    data = stanData,    # named list of data
    chains = 2,             # number of Markov chains
    warmup = 1000,          # number of warmup iterations per chain
    iter = 2000,            # total number of iterations per chain
    cores = 1,              # number of cores (could use one per chain)
    refresh = 0             # no progress shown
  )

  # point prediction

  scores_tmax <- summary(fit)$summary[1:K, "mean"]
  pred_eta <- mu+efuncs%*%scores_tmax

  if(interval == TRUE){
  # credible prediction interval using sampling quantiles
    score_draws <- as.matrix(fit)[, 1:K]
    eta_draws <- mu+efuncs %*% t(score_draws)
    eta_draws <- apply(eta_draws, 1, quantile, probs = c(0.025, 0.975))
    pred = data.frame(pred=pred_eta,
                      pred_lb = eta_draws[1, ],
                      pred_ub = eta_draws[2, ])
  }
  else{
    pred = data.frame(pred=pred_eta)
  }


  return(list(pred = pred_eta, scores = scores_tmax))

}
