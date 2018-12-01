#' the confidence interval for our beta estimator
#' 
#' @param coef_est the estimate for beta
#' @param vcov the variance covariance matrix for beta estimate
#' @param alpha the prespecified level
#' @return a list containing the (1 - alpha) level confidence interval
#' @importFrom stats qnorm
betaconfint <- function(coef_est, vcov, alpha)
{
  upper <- coef_est + qnorm(1 - alpha/2) * sqrt(diag(vcov))
  lower <- coef_est - qnorm(1 - alpha/2) * sqrt(diag(vcov))
  res <- list(upper = upper,
              lower = lower)
  res
}

#' the pointwise upper confidence interval for the survival curve
#' 
#' @param hazard_pred the predicted hazard function
#' @param hazardpredvar_est the variance of the estimator of the hazard function
#' @param newobsz the new obtained Z value
#' @param alpha the prespecified level
#' @return the upper (1 - alpha) level pointwise confidence interval for the hazard function
#' @importFrom stats qnorm
upperconfint <- function(hazard_pred, hazardpredvar_est, newobsz, alpha)
{
  upperconf_int <- pmin(1, exp(-(hazard_pred - qnorm(1 - alpha/2) * sqrt(hazardpredvar_est))))
  upperconf_int
}

#' the pointwise lower confidence interval for the survival curve
#' 
#' @param hazard_pred the predicted hazard function
#' @param hazardpredvar_est the variance of the estimator of the hazard function
#' @param newobsz the new obtained Z value
#' @param alpha the prespecified level
#' @return the lower (1 - alpha) level pointwise confidence interval for the hazard function
#' @importFrom stats qnorm
lowerconfint <- function(hazard_pred, hazardpredvar_est, newobsz, alpha)
{
  #constant <- 0.05
  #lowerconf_int <- exp(-exp(log(hazard_pred + constant) + qnorm(1 - alpha/2) * sqrt(hazardpredvar_est) / (hazard_pred + constant)))
  #lowerconf_int[hazard_pred == 0] = 1
  lowerconf_int <- exp(-(hazard_pred + qnorm(1 - alpha/2) * sqrt(hazardpredvar_est)))
  lowerconf_int
}

#' the pointwise confidence interval for the survival curve
#' 
#' @param hazard_pred the predicted hazard function
#' @param newobsz the new obtained Z value
#' @param fit the fitting object after fitting our model
#' @param alpha the prespecified level
#' @return a list containing (1 - alpha) level pointwise confidence interval for the hazard function
survprobconfint <- function(hazard_pred, newobsz, fit = NULL, alpha)
{
  fit <- hazardpredvarest(newobsz, fit)
  hazardpredvar_est <- fit$hazardpredvar_est
  upper <- upperconfint(hazard_pred, hazardpredvar_est, newobsz, alpha)
  lower <- lowerconfint(hazard_pred, hazardpredvar_est, newobsz, alpha)
  res <- list(upper = upper,
              lower = lower)
  fit$byprod <- append(fit$byprod, res)
  fit
}
