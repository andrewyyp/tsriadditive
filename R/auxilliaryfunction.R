#' sone function
#' 
#' The s_one function defined in the paper
#' @param survtime the event time
#' @param cause the indicator records the cause. Default to all one. Zero means right censoring. Greater than
#' or equal to two means other cause.
#' @param Z a variable contains all the regressors
#' @param comp the indicator of whether modeling subdistribution hazard
#' @param censorsurv the estimate for the censoring distribution
#' @return s_one defined in the paper
sone <- function(survtime, cause = NULL, Z, comp = FALSE, censorsurv = NULL)
{
  #survtime in increasing order
  s_one <- apply(apply(apply(Z, 2, rev), 2, cumsum), 2, rev)
  if(comp){
    censorest <- apply(Z * (cause >= 2) / censorsurv, 2, cumsum) * censorsurv
    censorest <- rbind(rep(0, ncol(censorest)), censorest[1:(nrow(censorest) - 1), ])
    s_one <- s_one + censorest
  }
  s_one
}

#' szero function
#' 
#' The s_zero function defined in the paper
#' @param survtime the event time
#' @param cause the indicator records the cause. Default to all one. Zero means right censoring. Greater than
#' or equal to two means other cause.
#' @param comp the indicator of whether modeling subdistribution hazard
#' @param censorsurv the estimate for the censoring distribution
#' @return s_zero defined in the paper
szero <- function(survtime, cause = NULL, comp = FALSE, censorsurv = NULL)
{
  #survtime in increasing order
  s_zero <- seq(length(survtime), 1)
  if(comp){
    censorest <- cumsum((cause >= 2) / censorsurv) * censorsurv
    censorest <- c(0, censorest[1:(length(censorest) - 1)])
    s_zero <- s_zero + censorest
  }
  s_zero
}

#' Zbar function
#' 
#' A function for computing Z_bar in the paper
#' @param s_zero s_zero function 
#' @param s_one s_one function
#' @return the Z_bar value in the paper, no difference for two settings
Zbar <- function(s_zero, s_one)
{
  Z_bar <- s_one/s_zero
  Z_bar
}

#' Zint function
#' 
#' A function for computing the integration for Z_bar from zero to some time t
#' @param survtime the event time
#' @param Z_bar the Z_bar value
#' @return the integration of Z_bar from zero to all the event time t
Zint <- function(survtime, Z_bar)
{
  timediff <- c(survtime[1], diff(survtime))
  Z_int <- apply(Z_bar * timediff, 2, cumsum)
  Z_int
}

#' Gint function
#' 
#' This is the integration of our censorsurv from somewhere to infinity
#' @param survtime the event time
#' @param censorsurv the estimate for the censoring distribution
Gint <- function(survtime, censorsurv)
{
  timediff <- c(survtime[1], diff(survtime))
  temp <- cumsum(rev(timediff * censorsurv))
  temp <- c(0, temp[1:(length(temp) - 1)])
  G_int <- rev(temp)
  G_int
}

#' GZ_int function
#' 
#' This is the integration of our censorsurv from somewhere to infinity
#' @param survtime the event time
#' @param Z_bar defined in the paper
#' @param censorsurv the estimate for the censoring distribution
GZint <- function(survtime, Z_bar, censorsurv)
{
  #This is the integration of our censorsurv * Z_bar from somewhere to infinity
  timediff <- c(survtime[1], diff(survtime))
  temp <- apply(apply(timediff * censorsurv * Z_bar, 2, rev), 2, cumsum)
  temp <- rbind(rep(0, ncol(Z_bar)), temp[1:(nrow(temp) - 1), ])
  GZ_int <- apply(temp, 2, rev)
  GZ_int
}

#' A scoreprocess function
#' 
#' This gives us the scoreprocess function defined in the paper
#' @param Z a variable contains all the regressors
#' @param Z_bar defined in the paper
#' @param cause the indicator records the cause. Default to all one. Zero means right censoring. Greater than
#' or equal to two means other cause.
#' @return the score_process
scoreprocess <- function(Z, Z_bar, cause)
{
  score_process <- (Z - Z_bar) * (cause == 1)
  score_process
}

#' Omega inverse function
#' 
#' This is for inverse of omega, which is part of the sandwich estimator
#' @param N the sample size
#' @param survtime the event time
#' @param Z a variable contains all the regressors
#' @param Z_int the integration of Z_bar
#' @param cause the indicator records the cause. Default to all one. Zero means right censoring. Greater than
#' or equal to two means other cause
#' @param comp indicator of whether we are under competing risks setting
#' @param censorsurv the estimate for the censoring distribution
#' @param G_int the integration of G function
#' @param GZ_int the integration of GZ
#' @return the Omega matrix appear in the paper
omegainv <- function(N, survtime, Z, Z_int, cause, comp = FALSE, censorsurv = NULL, G_int = NULL, GZ_int = NULL)
{
  if(!comp) omega <- t(Z) %*% (Z * survtime - Z_int)
  else omega <- t(Z) %*% (Z * survtime - Z_int + (cause >= 2) * (Z * G_int - GZ_int) / censorsurv)
  omega_inv <- solve(omega / N)
  omega_inv
}

#' The finite dimensional coefficients estimator
#' 
#' This function returns us the estimate the finite dimensional part
#' @param N the sample size
#' @param omega_inv the inverse of omega matrix in the paper
#' @param score_process the score process in the paper
#' @return the coefficent estimate from the fit
coefest <- function(N, omega_inv, score_process)
{
  coef_est <- omega_inv %*% apply(score_process / N, 2, sum)
  coef_est
}

#' The baseline hazards function
#' 
#' This gives the estimate of baseline hazards function
#' @param cause the indicator records the cause. Default to all one. Zero means right censoring. Greater than
#' or equal to two means other cause.
#' @param s_zero the S_zero in the paper
#' @param Z_int the integration of Z_bar
#' @param coef_est the coefficent estimate from the fit
#' @return the baseline hazard function estimate
baselineest <- function(cause, s_zero, Z_int, coef_est)
{
  baseline_est <- cumsum((cause == 1) / s_zero) - Z_int %*% coef_est
  #baseline_est <- cummax(pmax(baseline_est, 0))
  baseline_est
}

#' A hazard prediction function
#'
#' the predict function associated with our class
#' @param fit the fitting object after fitting our model
#' @param newtreatment the new treatment value
#' @param newIV new instrumental variable value
#' @param newcovariates new observed covariates
#' @return an object recording the corresponding predicted survival curve and corresponding pointwise
#' confidence interval
#' @importFrom stats predict
hazardpred <- function(fit, newtreatment = NULL, newIV = NULL, newcovariates = NULL)
{
  coef_est <- fit$coef
  baseline_est <- fit$baseline
  useIV <- fit$byprod$useIV
  firstfit <- fit$byprod$firstfit
  survtime <- fit$byprod$survtime
  binary <- fit$byprod$binary
  if(!useIV){
    newobsz <- c(newtreatment, newcovariates) 
    score_pred <- sum(newobsz * coef_est)
    hazard_pred = baseline_est + sum(newobsz * coef_est) * survtime
  }
  else{
    newdata = data.frame(1, newIV, t(newcovariates))
    colnames(newdata) <- names(firstfit$coefficients)
    if(binary) newresidual = newtreatment - predict(firstfit, newdata = newdata, type = "response")
    else newresidual = newtreatment - predict(firstfit, newdata)
    newobsz = c(newtreatment, newresidual, newcovariates) 
    score_pred <- sum(newobsz * coef_est)
    hazard_pred <- baseline_est + score_pred * survtime
  }
  hazard_pred <- cummax(pmax(hazard_pred, 0))
  res <- list(score_pred = score_pred,
              hazard_pred = hazard_pred,
              newobsz = newobsz)
  res
}

#' A survival prediction function
#'
#' The predict function associated with our class
#' @param fit the fitting object after fitting our model
#' @param newtreatment the new treatment value
#' @param newIV new instrumental variable value
#' @param newcovariates new observed covariates
#' @return an object recording the corresponding predicted survival curve and corresponding pointwise
#' confidence interval
survivalpred <- function(fit, newtreatment = NULL, newIV = NULL, newcovariates = NULL)
{
  hazard_pred <- hazardpred(fit, newtreatment, newIV, newcovariates)
  survival_pred <- exp(-hazard_pred$hazard_pred)
  res <- append(list(survival_pred = survival_pred), hazard_pred)
  res
}

