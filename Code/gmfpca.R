#' Generalized multilevel functional principal components analysis
#'
#' The GM-FPCA method comprises of four steps:
#' \enumerate{
#' \item Divide the domain into bins;
#' \item Fit a mixed model in each bin and extract the latent linear predictors;
#' \item Apply multilevel FPCA on linear predictors to obtain eigenfunctions;
#' \item Get posterior samples of scores using the Bayesian method.
#' }
#'
#' @param formula Two-sided formula object in lme4 formula syntax.
#' @param data A data frame in long format containing subject identifier, visit identifier, location of observations on the functional domain, and functional observations.
#' @param family GLM family of the response. Currently supports \code{binomial} and \code{poisson}.
#' @param argvals A vector containing locations of observations on the
#' functional domain. If not specified, a regular grid across the range of
#' the domain is assumed.
#' @param bin.pars A named list of parameters related to binning. \code{bin_midpoints} is a vector specifying
#' the midpoints of bins, \code{bin_len} is an integer specifying the bin width, and \code{cyclic} is a logical
#' specifying whether bins should be cyclic.
#' @param var.names A named list of the variable names for id, visit, and domain in the data set.
#' @param mfpca.pars A named list of parameters related to multilevel FPCA. \code{pve} specifies the minimum
#' proportion of variance explained; \code{npc_l1} and \code{npc_l2} specify the number of principal components
#' for level 1 and level 2, respectively. See \code{mfpca.face()} from the package \code{refund} for details.
#' @param stan.pars A named list of parameters related to running STAN. List components include
#' \code{iter_warmup}, \code{iter_sampling}, \code{chains}, and \code{refresh}.
#' See package \code{cmdstanr} for details.
#' @param dir.path A named list of parameters related to saving intermediate results.
#' \code{save_intermediate_res} is logical specifying whether intermediate results should be saved. Default to \code{FALSE}.
#' \code{intermediate_res_path} is a string specifying the path for where intermediate results are saved.
#' @param nknots_min Minimum number of knots in penalized smoothing when projecting eigenfunctions back to the original grid.
#' This parameter will only be used when the number of bins is less than the number of observations on the functional domain.
#' Defaults to \code{35},
#' @param parallel Logical, indicating whether to do parallel computing.
#' Defaults to \code{FALSE}.
#' @param num_cores Number of cores for parallelization. Defaults to 1.
#' @param silent Logical, indicating whether to show descriptions of each step.
#' Defaults to \code{FALSE}.
#'
#' @return A list containing:
#' \item{argvals}{A vector of locations of functional observations.}
#' \item{id_visit}{A data frame containing the subject and visit identifers for each functional observation.}
#' \item{efunctions}{A list containing eigenfunctions for level 1 and level 2.}
#' \item{scores_posterior}{A list containing the full posterior sample of scores for level 1 and level 2.}
#' \item{scores_postmean}{A list containing the posterior mean of scores for level 1 and level 2.}
#'
#' @references
#'
#' @import dplyr
#' @import cmdstanr
#'
#' @export gmfpca
#'
#' @examples
#' set.seed(1)
#' dat <- sim_gmfd(family = "binomial", I  = 100, Jmean = 5, L = 100)
#' fit_gmfpca <- gmfpca(y ~ (1|id) + (1|id:visit), data = dat$data, family = "binomial", argvals = 1:100)

# library(dplyr)
# set.seed(1)
# dat <- sim_gmfd()
# data <- dat$data
# formula <- y ~ 1 + (1|id) + (1|id:visit)
# argvals <- 1:100
# family = "binomial"

gmfpca <- function(
    formula,
    data,
    family,
    argvals = NULL,
    bin.pars = list(
      bin_midpoints = NULL,
      bin_len = 20,
      cyclic  = T
    ),
    var.names = list(
      id_name = "id",
      visit_name = "visit",
      domain_name = "time"
    ),
    mfpca.pars = list(
      pve = 0.99,
      npc_l1 = 4,
      npc_l2 = 4
    ),
    stan.pars = list(
      iter_warmup = 1000,
      iter_sampling = 2000,
      chains = 4,
      refresh = 100
    ),
    dir.path = list(
      save_intermediate_res = FALSE,
      intermediate_res_path = "results/sim/"
    ),
    nknots_min = 35,
    parallel = FALSE,
    num_cores = 8,
    silent = FALSE
    ){

  # If doing parallel computing, set up the number of cores
  if(parallel & !is.integer(num_cores) ) num_cores <- parallel::detectCores() - 1

  # formula <- y ~ (1 | id) + (1 | id:visit)
  # Organize the input from the model formula
  model_formula <- as.character(formula)
  if(!(model_formula[1] == "~" & length(model_formula) == 3)){
    stop(
      paste0(
        "The formula is specified incorrectly. It should look like ",
        "y ~ (1 | id) + (1 | id:visit) with the appropriate variable names."
      )
    )
  }

  # create a folder for saving intermediate results
  if(dir.path$save_intermediate_res){
    res.path <- dir.path$intermediate_res_path
    if(!file.exists(res.path)) dir.create(res.path)
  }

  ##############################################################################
  ## Step 1 Create bins
  ##############################################################################
  if(!silent) print("Step 1: Create bins.")

  # L: dimension of the functional domain
  L <- max(data[[var.names$domain_name]])
  if(is.null(argvals)){
    argvals <- 1:L
  }else{
    # when only using a subset of the functional domain
    if(max(argvals) > L){
      stop(
        paste0(
          "Maximum index specified in argvals is greater than the",
          "maximum number of observations in the functional domain."
        )
      )
    }
    # L <- length(argvals)
  }
  if(L > 400){
    message(
      paste0(
        "Your functional data has ", L, "measurements along the functional domain, ",
        "consider subsampling to cut the computational cost."
      )
    )
  }

  # Create bins
  if(is.null(bin.pars$bin_midpoints)){
    # If bin midpoints are not specified,
    # create a bin at each point along the functional domain.
    bin.pars$bin_midpoints <- argvals
    # c(
    #   seq(argvals[1], argvals[length(argvals)], by = 10),
    #   argvals[length(argvals)]
    # )
  }
  bin.index <- CreateBins(
    argvals = argvals,
    bin_midpoints = bin.pars$bin_midpoints,
    bin_len = bin.pars$bin_len,
    cyclic = bin.pars$cyclic
  )
  # Count the number of visits per subject
  nvisits <- data %>% dplyr::summarise(
    nvisits = n_distinct(.data[[var.names$visit_name]]),
    .by = all_of(var.names$id_name)
  )
  # Create a matrix to store AICs
  # AIC_mat <- matrix(NA, nrow = L, ncol = 2)


  ##############################################################################
  ## Step 2 Fit mixed models
  ##############################################################################
  if(!silent) print("Step 2: Fit mixed models.")

  ### Create a function that fit a mixed model in bin l
  ### Input: l: index of the list bin.index.
  ###           Each element in the list is a vector containing domain points for that bin.
  ### Output: A list containing point estimates, variance estimates, etc.
  get_eta <- function(l){
    dtemp <- data %>% dplyr::filter(.data[[var.names$domain_name]] %in% l)
    fit_mm <- suppressMessages(
      lme4::glmer(
        formula = formula,
        data = dtemp,
        family = family,
        control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 5000))
      )
    )
    coef_df <- coef(fit_mm)
    # re_df <- lme4::ranef(fit_mm) ## random effects
    return(
      list(
        re_names = rownames(coef_df[[paste0(var.names$id_name, ":", var.names$visit_name)]]),
        # eta_vals = rep(re_df[[var.names$id_name]][[1]], times = nvisits$nvisits) +
        #   re_df[[paste0(var.names$id_name, ":", var.names$visit_name)]][[1]],
        eta_vals = rep(coef_df[[var.names$id_name]][[1]], times = nvisits$nvisits) +
          coef_df[[paste0(var.names$id_name, ":", var.names$visit_name)]][[1]],
        isSingular = lme4::isSingular(fit_mm),
        warn = paste(unlist(summary(fit_mm)$optinfo$conv$lme4$messages), collapse = ". ")
      )
    )
  }
  # Fit mixed models
  if(parallel == TRUE){
    mms <- parallel::mclapply(bin.index, get_eta, mc.cores = num_cores)
  }else{
    mms <- lapply(bin.index, get_eta)
  }
  if(dir.path$save_intermediate_res) saveRDS(mms, paste0(res.path, "mixed-model-results.rds"))
  # mms <- read_rds(paste0(res.path, "mixed-model-results.rds"))
  # mms <- read_rds("results/sim/mixed-model-results.rds")
  eta_vals <-  t(do.call(rbind, lapply(mms, '[[', "eta_vals")))
  # Extract warning messages
  mm_summary <- list(
    n_singular = do.call(sum, lapply(mms, '[[', "isSingular")),
    warning = do.call(c, lapply(mms, '[[', "warn"))
  )

  ##############################################################################
  ## Step 3 MFPCA
  ##############################################################################
  if(!silent) print("Step 3: MFPCA.")

  id.visit <- strsplit(mms[[1]][["re_names"]], ":")
  id.visit <- data.frame(
    id = sapply(id.visit, function(x) as.integer(x[1])),
    visit = sapply(id.visit, function(x) as.integer(x[2]))
  )
  mfpca.res <- refund::mfpca.face(
    eta_vals,
    id = id.visit$id,
    visit = id.visit$visit,
    argvals = bin.pars$bin_midpoints,
    pve = mfpca.pars$pve
  )
  if(dir.path$save_intermediate_res) saveRDS(mfpca.res, paste0(res.path, "mfpca-results.rds"))
  # mfpca.res <- read_rds(paste0(res.path, "mfpca-results.rds"))
  # mfpca.res <- read_rds("results/sim/mfpca-results.rds")

  # mfpca.res1 <- mfpca.res
  # par(mfrow = c(1,2))
  # plot(mfpca.res1$mu, main = "use Eta", ylim = c(-0.1, .25), type = "l")
  # plot(mfpca.res$mu, main = "Use RE", ylim = c(-0.1, .25), type = "l")
  # par(mfrow = c(2,4))
  # plot(mfpca.res1$efunctions$level1[, 1], type = "l", main = "use Eta", ylab = "l1e1")#, ylim = c(-0.1, .25))
  # plot(mfpca.res1$efunctions$level1[, 2], type = "l", main = "use Eta", ylab = "l1e2")
  # plot(mfpca.res1$efunctions$level1[, 3], type = "l", main = "use Eta", ylab = "l1e3")
  # plot(mfpca.res1$efunctions$level1[, 4], type = "l", main = "use Eta", ylab = "l1e4")
  # plot(mfpca.res$efunctions$level1[, 1], type = "l", main = "Use RE", ylab = "l1e1")#, ylim = c(-0.1, .25))
  # plot(mfpca.res$efunctions$level1[, 2], type = "l", main = "Use RE", ylab = "l1e2")
  # plot(mfpca.res$efunctions$level1[, 3], type = "l", main = "Use RE", ylab = "l1e3")
  # plot(mfpca.res$efunctions$level1[, 4], type = "l", main = "Use RE", ylab = "l1e4")
  # par(mfrow = c(2,4))
  # plot(mfpca.res1$efunctions$level2[, 1], type = "l", main = "use Eta", ylab = "l2e1")#, ylim = c(-0.1, .25))
  # plot(mfpca.res1$efunctions$level2[, 2], type = "l", main = "use Eta", ylab = "l2e2")
  # plot(mfpca.res1$efunctions$level2[, 3], type = "l", main = "use Eta", ylab = "l2e3")
  # plot(mfpca.res1$efunctions$level2[, 4], type = "l", main = "use Eta", ylab = "l2e4")
  # plot(mfpca.res$efunctions$level2[, 1], type = "l", main = "Use RE", ylab = "l2e1")#, ylim = c(-0.1, .25))
  # plot(mfpca.res$efunctions$level2[, 2], type = "l", main = "Use RE", ylab = "l2e2")
  # plot(mfpca.res$efunctions$level2[, 3], type = "l", main = "Use RE", ylab = "l2e3")
  # plot(mfpca.res$efunctions$level2[, 4], type = "l", main = "Use RE", ylab = "l2e4")


  ##############################################################################
  ## Step 4 Run STAN
  ##############################################################################
  if(!silent) print("Step 4: STAN")
  # deigfun.l1.wide <- mfpca.res$efunctions$level1[, 1:mfpca.pars$npc_l1]
  # deigfun.l2.wide <- mfpca.res$efunctions$level2[, 1:mfpca.pars$npc_l2]

  # Project eigenfunctions back to the full grid
  if(length(bin.pars$bin_midpoints) < length(argvals)){
    nknots <- min(round(length(bin.pars$bin_midpoints)/2), nknots_min)
    deigfun.l1.wide <- reeval_efunctions(
      knots = nknots,
      argvals_bin = bin.pars$bin_midpoints,
      argvals = argvals,
      efunctions = mfpca.res$efunctions$level1[, 1:mfpca.pars$npc_l1],
      npc = mfpca.pars$npc_l1
    )
    deigfun.l2.wide <- reeval_efunctions(
      knots = nknots,
      argvals_bin = bin.pars$bin_midpoints,
      argvals = argvals,
      efunctions = mfpca.res$efunctions$level2[, 1:mfpca.pars$npc_l2],
      npc = mfpca.pars$npc_l2
    )
    b0.mfpca <- reeval_efunctions(
      knots = nknots,
      argvals_bin = bin.pars$bin_midpoints,
      argvals = argvals,
      efunctions = data.frame(mfpca.res$mu),
      npc = 1
    )[, 1]
  } else {
    deigfun.l1.wide <- mfpca.res$efunctions$level1[, 1:mfpca.pars$npc_l1]
    deigfun.l2.wide <- mfpca.res$efunctions$level2[, 1:mfpca.pars$npc_l2]
    b0.mfpca  <- mfpca.res$mu
  }

  standata <- list(
    Nsubs = length(unique(data[[var.names$id_name]])),
    Nfuns = sum(nvisits$nvisits),
    S     = length(argvals),
    L     = mfpca.pars$npc_l1,
    M     = mfpca.pars$npc_l2,
    ylen  = nrow(data),
    y     = data[[model_formula[[2]]]],
    nvisits   = nvisits$nvisits,
    beta_0    = b0.mfpca,
    efuncs_l1 = t(deigfun.l1.wide),
    efuncs_l2 = t(deigfun.l2.wide),
    grainsize = 1
  )
  if(dir.path$save_intermediate_res){
    write_stan_json(
      standata,
      file = paste0(res.path, "standata.json"),
      always_decimal = FALSE
    )
  }

  # standata <- jsonlite::read_json("results/sim/standata.json")
  # set_cmdstan_path(dir.path$cmdstan_installation_path)
  # stanmodel <- system.file(paste0("stanmodel-", family), package = "fli")
  # model_parallel <- cmdstanr::cmdstan_model(exe_file = stanmodel)
  # if(file.exists(stanmodel)){
  # } else {
  #   model_parallel <-  cmdstan_model(
  #     stan_file = paste0(paste0(dir.path$stanfile_path, "stanmodel-", family), ".stan")
  #   )
  # }

  # new way, using {instantiate}
  stanmodel <- instantiate::stan_package_model(
    name = paste0("stanmodel-", family), package = "fli"
  )

  tic <- Sys.time()
  stanfit <- stanmodel$sample(
    # data = paste0(res.path, "standata.json"),
    data = standata,
    iter_warmup = stan.pars$iter_warmup,
    iter_sampling = stan.pars$iter_sampling,
    chains = stan.pars$chains,
    parallel_chains = num_cores,
    refresh = stan.pars$refresh
  )
  toc <- Sys.time()
  toc - tic
  draws.df <- stanfit$draws(format = "df")
  if(dir.path$save_intermediate_res) saveRDS(draws.df, file = paste0(res.path, "mcmcdraws.rds"))
  # draws.df <- read_rds(paste0(res.path, "mcmcdraws.rds"))
  # draws.df <- read_rds(paste0("results/sim/mcmcdraws.rds"))
  # sigma.sq <- draws.df[, grep("sigma_epi_sq", colnames(draws.df))]
  # plot(sigma.sq$sigma_epi_sq, type = "l")


  l1score.posterior <- array(NA, c(standata$Nsubs, stan.pars$iter_sampling, mfpca.pars$npc_l1))
  for(i in 1:standata$Nsubs){
    # print(i)
    for(k in 1:mfpca.pars$npc_l1){
      temp.varname <- paste0("xil1_sc[", i, ",", k, "]")
      l1score.posterior[i, , k] <- draws.df[[temp.varname]]
    }
  }
  l1score.post.mean <- apply(l1score.posterior, c(1,3), mean)


  l2score.posterior <- array(NA, c(standata$Nfuns, stan.pars$iter_sampling, mfpca.pars$npc_l2))
  for(i in 1:standata$Nfuns){
    # print(i)
    for(k in 1:mfpca.pars$npc_l2){
      temp.varname <- paste0("xil2_sc[", i, ",", k, "]")
      l2score.posterior[i, , k] <- draws.df[[temp.varname]]
    }
  }
  l2score.post.mean <- apply(l2score.posterior, c(1,3), mean)

  return(
    list(
      argvals  = argvals,
      id_visit = id.visit,
      efunctions = list(
        level1 = deigfun.l1.wide,
        level2 = deigfun.l2.wide
      ),
      scores_posterior = list(
        level1 = l1score.posterior,
        level2 = l2score.posterior
      ),
      scores_postmean = list(
        level1 = l1score.post.mean,
        level2 = l2score.post.mean
      )
    )
  )
}


# fli.res <- list(
#   argvals  = argvals,
#   id_visit = id.visit,
#   betaHat  = betaHat,
#   betaSE   = betaSE,
#   beta_posterior = beta.posterior,
#   beta_post_mean = beta.post.mean,
#   efunctions = list(
#     level1 = deigfun.l1.wide,
#     level2 = deigfun.l2.wide
#   ),
#   scores_posterior = list(
#     level1 = l1score.posterior,
#     level2 = l2score.posterior
#   ),
#   scores_postmean = list(
#     level1 = l1score.post.mean,
#     level2 = l2score.post.mean
#   ),
#   ci.level = ci.level
# )

# format scores
# colQuantile <- function(A, prob){
#   apply(A, 2, function(x) quantile(x, probs = prob))
# }
# l1score.posterior.m1 <- draws.df[, grep("xil1_sc", colnames(draws.df))]
# l1score.post.mean.m1 <- cbind(
#     data.frame(id = unique(id.visit$id)),
#     data.frame(matrix(colMeans(l1score.posterior.m1), ncol = mfpca.pars$npc_l1, byrow = F))
# )
# l2score.posterior.m1 <- draws.df[, grep("xil2_sc", colnames(draws.df))]
# l2score.post.mean.m1 <- cbind(
#   id.visit,
#   data.frame(
#     matrix(colMeans(l2score.posterior.m1), ncol = mfpca.pars$npc_l2, byrow = F)
#   )
# )
# l1score.sd <- cbind(
#   data.frame(id = unique(id.visit$id)),
#   data.frame(
#     matrix(
#       apply(l1score.posterior, 2, sd),
#       ncol = mfpca.pars$npc_l1, byrow = F
#     )
#   )
# )
# l2score.sd <- cbind(
#   id.visit,
#   data.frame(
#     matrix(
#       apply(l2score.posterior, 2, sd),
#       ncol = mfpca.pars$npc_l2, byrow = F
#     )
#   )
# )
# scores_sd = list(
#   ci_level = ci.level,
#   level1 = l1score.sd,
#   level2 = l2score.sd
# )
# scores_quantile = list(
#   ci_level = ci.level,
#   level1_lower = l1score.lower,
#   level1_upper = l1score.upper,
#   level2_lower = l2score.lower,
#   level2_upper = l2score.upper
# )


# fli.obj$id_visit <- id.visit
# check traceplots
# plot(l1score.posterior$`xil1_sc[1,1]`, type = "l")

# fli.res$efunctions <- list(
#   level1 = deigfun.l1.wide,
#   level2 = deigfun.l2.wide
# )
# fli.res$scores_postmean <- scores_postmean
# fli.res$scores_quantile <- scores_quantile
# str(fli.res)


# fli.res$scores_sd = list(
#     ci_level = ci.level,
#     level1 = l1score.sd,
#     level2 = l2score.sd
#   )

# fli.res$scores_posterior = list(
#   level1 = l1score.posterior,
#   level2 = l2score.posterior
# )
