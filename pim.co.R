#' Fitting a Probabilistic Index Model for composite outcomes
#'
#' @param ID a vector of unique subject-level identifiers.
#' @param time a vector of event times.
#' @param status a vector of event type labels. 0: censoring, 1:death
#' and 2: non-fatal event.
#' @param Z a matrix or a vector of covariates.
#' @param model a single character value with possible values "difference"
#' (the default) and "marginal"
#' @param start start values for fitting estimating equations, default is 0 vector
#' @param link a character vector with a single value that determines the
#' used link function. Possible values are "logit", "probit" and "identity".
#' The default is "logit".
#' @param method "Broyden" and "Newton" can be chosen in \code{\link[nleqslv]{nleqslv}}
#' @param control conditons for convergence of \code{\link[nleqslv]{nleqslv}}
#' @heter when heter="TRUE", data has heterogeneity; when heter="FALSE", then data has no heterogeneity
#' note: heter="TRUE" is only work for continuous covariates 
#' @param weights Currently not implemented.
#' @param keep.data a logical value indicating whether the model
#' matrix should be saved in the object. Defaults to \code{FALSE}.
#' 
#' @return whether model is converged or not and corresponding estimating results
#' This function is created based on R pacakge "pim" by Joris Meys et al.(2020) (https://cran.r-project.org/web/packages/pim/index.html)
#'
#'@examples
#'library(WR)
#'library(pim)
#'library(nleqslv)
#'
#'
#'head(non_ischemic)
#'colnames(non_ischemic)[4:16]=c("Training vs Usual","Age (year)","Male vs Female","Black vs White", 
#  "Other vs White", "BMI","LVEF","Hypertension","COPD","Diabetes",
#  "ACE Inhibitor","Beta Blocker", "Smoker")
#' sample size
#' length(unique(non_ischemic$ID))
#' median(non_ischemic$time[non_ischemic$status<2])/30.5
#' table(non_ischemic$status)
#' get the number of rows and number of covariates.
#' nr <- nrow(non_ischemic)
#' p <- ncol(non_ischemic)-3
#' extract ID, time, status and covariates matrix Z from the data.
#' note that: ID, time and status should be column vector.
#' covariatesZ should be (nr, p) matrix.
#' ID <- non_ischemic[,"ID"]
#' time <- non_ischemic[,"time"]
#' status <- non_ischemic[,"status"]
#' Z <- as.matrix(non_ischemic[,4:(3+p)],nr,p)
#' 
#' using logit link
#' start1 = rep(0, ncol(Z))
#' pimco1 = pim.co(ID, time, status, Z, start=start1)
#' pimco1$flag
#' pimco1$coefficient
#' 
#' using probit link
#' pimco2 = pim_m(ID, time, status, Z, link="probit", model="difference",start=start1)
#' pimco2$coefficient
#' pimco2$flag
#' 
#' another example
#' data("FEVData")
#' n = dim(FEVData)[1]
#' ID = 1:n
#' time = FEVData$FEV
#' status = rep(1,n) ##assumed these subjects are failure time
#' Z <- as.matrix(FEVData[,c(1,4)],n,p)
#' pimco3 = pim.co(ID, time, status, Z, model = "marginal", start = rep(0, ncol(Z)), link = "identity", heter=FALSE)
#' pimco3$coefficient
#' pimco3$flag
#' 



pim.co = function(ID, time, status, Z, model = c("difference", "marginal"), 
                  start = rep(0, ncol(Z)), link = c("logit", "probit", "identity"), 
                  method = c("Broyden", "Newton"), control = list(),
                  heter = FALSE, weights = NULL, keep.data = FALSE){
  ID.save <- ID
  if (is.character(ID) | is.factor(ID)) {
    ID <- as.numeric(factor(ID, levels = unique(ID)))
  }
  data <- cbind(ID, time, status, Z)
  obj <- extract.times(data)
  datap <- obj$datap
  n <- nrow(datap)
  tau <- obj$tau
  p <- ncol(datap) - 4
  
  comij_res <- delta.data.co(datap, tau, model = model, heter = heter)
  outcome <- comij_res[[1]]
  x <- as.matrix(outcome[,-c(1:5,p+6)])
  w <- outcome$w
  t2 <- outcome$t2
  y <- abs(w)*eq(w,1)-w*(t2<tau)
  #make those fair outcome comparison not be used in score function
  x <- as.matrix(abs(w)*x, ncol=p)
  #when model = "difference" or "marginal", penv is also different
  penv = comij_res[[2]]
  res <- pim.co.fit(x, y, link = link, start = start, method = method,
                    control = control, weights = weights, penv = penv)
  names(res$coef) <- colnames(Z)
  if (!keep.data) {
    x <- matrix(nrow = 0, ncol = 0)
    y <- numeric(0)
  }
  flag = res$flag
  Estimate = res$coef
  Std.Error = sqrt(diag(res$vcov))
  z.value = Estimate/Std.Error
  p.value = rep(0,p)
  for(i in 1:p){
    z.valuei = z.value[i]
    if(z.valuei<0) p.value[i] = 2*pnorm(z.valuei) else p.value[i] = 2*(1-pnorm(z.valuei))
  }
  coefficient = data.frame(Estimate,Std.Error,z.value,p.value)
  res = list(coefficient = coefficient, coef = res$coef, vcov = res$vcov, fitted = res$fitted, 
             tau = tau, flag = flag, penv = penv, link = link, estimators = res$estim, 
             model.matrix = x, response = y, keep.data = keep.data, model = model)
  return(res)
}