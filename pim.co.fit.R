#' Fitter function for a probabilistic index model for composite outcomes
#' 
#' #' This is the basic computing engine called by \code{\link{pim.co}}
#' to get the estimates for the coefficients and the variance-
#' covariance matrices. This function currently only spits out
#' these components using the sandwich estimators.
#' 
#' #'
#'
#' @param x a model matrix with as many rows as \code{y}.
#' @param y a vector with the pseudo-responses
#' @param link a character vector with a link function
#' @param start A numeric vector with the starting values for the fitting
#' algorithm, if required.
#' @param method "Broyden" and "Newton" can be chosen in \code{\link[nleqslv]{nleqslv}}
#' @param control conditons for convergence of \code{\link[nleqslv]{nleqslv}}
#' @param vcov.estim a function to determine the variance-covariance matrix.
#' Possibilities are \code{\link{sandwich.vcov}} and \code{link{score.vcov}}.
#' Defaults to \code{sandwich.vcov}
#' @param weights currently not implemented
#' @param penv An environment
#'


#'
#' @return A list with the following elements
#' \describe{
#'  \item{coefficients}{a numeric vector with the coefficients}
#'  \item{vcov}{a numeric matrix with the variance-covarianc matrix for
#'  the coefficients}
#'  \item{fitted}{a numeric vector with the raw fitted values}
#'  \item{flag}{a flag indicating nleqslv is converged or not, its value is equal to termcd of nleqslv}
#'  \item{estim}{a list with two components named \code{coef} and \code{vcov}
#'  containing information on the used estimators for both.}
#' }
#' 
#' 
#' This function is created based on R pacakge "pim" by Joris Meys et al.(2020) (https://cran.r-project.org/web/packages/pim/index.html)


pim.co.fit <- function(x, y, link = "logit", start = rep(0, ncol(x)), 
                       method = c("Broyden", "Newton"), control = list(), 
                       vcov.estim = "sandwich.vcov", weights = NULL, penv){
  fn <- CreateScoreFun(x, y, link=link)
  res <- nleqslv(start, fn, method=method, control = control)
  fits <- x %*% res$x
  flag <- res$termcd    #1-converge; others-not fully converged
  vcov.estimF <- match.fun(vcov.estim)
  dim(fits) <- NULL
  if (!is.list(penv)) penv <- poset(penv, as.list = TRUE)
  vc <- vcov.estimF(fits, x, y, weights, link, penv)
  return(list(coefficients = res$x, vcov = vc, fitted = fits, flag = flag,
              estim = list(coef = "estimator.nleqslv", vcov = vcov.estim)))
}