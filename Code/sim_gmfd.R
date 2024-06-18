#' Simulate generalized multilevel functional data
#'
#' @param family Type of data to generate. Currently supports "binomial" and "poisson".
#' @param I Number of subjects.
#' @param Jmean Mean number of functions per subject.
#' @param L Number of grid points on the functional domain.
#' @param orthonormal Logical. Whether eigenfunctions from level one should be orthogonal
#' to those from level two.
#'
#' @export sim_gmfd
#'
#' @return A list containing:
#' \item{data}{A data frame of the simulated generalized multilevel functional
#' data containing variables \code{id}, \code{visit}, and \code{Y}}
#' \item{eta}{}
#' \item{npc}{}
#' \item{efunctions}{}
#' \item{evalues}{}
#' \item{scores}{}

sim_gmfd <- function(
    family = "binomial",
    I  = 100, Jmean = 5, L = 100,
    orthonormal = TRUE,
    npc        = list(level1 = 4,    level2 = 4),
    efunctions = list(level1 = NULL, level2 = NULL),
    evalues    = list(level1 = NULL, level2 = NULL),
    scores     = list(level1 = NULL, level2 = NULL)
){
  times <- seq(0, 1, length.out = L)

  ## generate the number of visits for each subject from the Poisson distribution
  Js <- pmax(rpois(I, Jmean), 1)
  nfun <- sum(Js)

  # generate eigenfunctions and scores
  if(is.null(efunctions$level1)){
    efunctions$level1 <- sapply(
      seq(1, npc$level1, by = 1), function(k) sqrt(2) * sin(k * pi * times)
    )
  }
  if(orthonormal){
    efunctions$level2 <- sapply(
      seq(1, npc$level2, by = 1), function(k) sqrt(2) * cos(k * pi * times)
    )
  } else {
    l2temp <- funData::eFun(times, M = L, type = "Poly")
    efunctions$level2 <- t(l2temp@X)
  }

  if(is.null(evalues$level1)){
    evalues$level1 <- 0.5^seq(from = 0, to = npc$level1 - 1, by = 1)
  }
  if(is.null(evalues$level2)){
    evalues$level2 <- 0.5^seq(from = 0, to = npc$level2 - 1, by = 1)
  }
  if(is.null(scores$level1)){
    scores$level1 <- sapply(
      seq(1, npc$level1, by = 1), function(k) rnorm(I, mean = 0, sd = sqrt(evalues$level1[k]))
    )
  }
  if(is.null(scores$level2)){
    scores$level2 <- sapply(
      seq(1, npc$level2, by = 1), function(k) rnorm(nfun, mean = 0, sd = sqrt(evalues$level2[k]))
    )
  }

  ## generate random effects
  subid <- as.factor(rep(1:I, times = Js))
  eta.design <- model.matrix( ~ 0 + subid)
  eta.l1 <- scores$level1 %*% t(efunctions$level1)
  eta <- eta.design %*% eta.l1
  eta <- eta + scores$level2 %*% t(efunctions$level2)

  ## generate generalized multilevel functional data
  Y <- matrix(NA, nfun, L)
  if(family == "binomial"){
    for(i in 1:nfun){
      Y[i, ] <- rbinom(L, size = 1, prob = plogis(eta[i, ]))
    }
  }else if(family == "poisson"){
    for(i in 1:nfun){
      Y[i, ] <- rpois(n = 1, lambda = exp(eta[i, ]))
    }
  }
  visit <- sequence(Js)
  dat.long <- tidyr::pivot_longer(
    cbind(
      data.frame(id = subid, visit = visit),
      data.frame(Y)
    ),
    cols = -c(id, visit), names_to = "time", values_to = "y"
  )
  dat.long$time <- as.numeric(gsub("X", "", dat.long$time))
  return(
    list(
      data = dat.long,
      eta = eta,
      npc = npc,
      efunctions = efunctions,
      evalues = evalues,
      scores = scores
    )
  )
}




# SimulateGMFD_Easy <- function(
#     family = "gaussian",
#     I,
#     L,
#     Js, # number of visits per subject
#     eta
# ){
#   times <- seq(0, 1, length.out = L)
#   nfun <- sum(Js)
#   subid <- as.factor(rep(1:I, times = Js))
#
#   ## generate generalized multilevel functional data
#   if(family == "binomial"){
#     for(i in 1:nfun){
#       Y[i, ] <- rbinom(L, size = 1, prob = plogis(eta[i, ]))
#     }
#   }else if(family == "poisson"){
#     for(i in 1:nfun){
#       Y[i, ] <- rpois(n = 1, lambda = exp(eta[i, ]))
#     }
#   }
#
#   visit <- sequence(Js)
#   return(
#     list(
#       data = data.frame(
#         id = subid,
#         visit = visit,
#         Y = I(Y)
#       )
#     )
#   )
# }







