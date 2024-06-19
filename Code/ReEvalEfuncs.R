
#' Get new spline basis for eigenfunctions
#'
#' Internal function for generating spline basis from fpca.face results to evaluate eigenfunctions
#' on expanded grid.
#'
#' @param knots Number of knots for spline basis
#' @param argvals_bin Index on for binned data.
#' @param efunctions Eigenfunctions from an fpca object.
#' @param argvals Full index on which to generate new basis.
#' @param npc Number of principal components.
#'
#' @importFrom stats lm
#' @importFrom splines spline.des
#'
reeval_efunctions <- function(knots, argvals_bin, argvals, efunctions, npc){
  
  p <- 3 # order of b splines
  
  ## get knot location in the same way that fpca.face() does internally
  ## (see refund:::pspline.setting())
  knots_values <- seq(-p, knots + p, length = knots + 1 + 2 *p)/knots
  knots_values <- knots_values * (max(argvals_bin) - min(argvals_bin)) + min(argvals_bin)
  
  ## evaluate B splines on original grid
  B <- splines::spline.des(knots = knots_values, x = argvals_bin, ord = p + 1,
                           outer.ok = TRUE)$design
  
  Bnew <- splines::spline.des(knots = knots_values, x = argvals, ord = p + 1,
                              outer.ok = TRUE)$design
  
  ## project first eigenfunction onto the B spline basis used in model fitting
  efunctions_new <- matrix(NA, length(argvals), npc)
  for(k in 1:npc){
    lm_mod <- lm(efunctions[,k] ~ B-1)
    efunctions_new[,k] <- Bnew %*% coef(lm_mod)
  }
  
  efunctions_new
  
}
