#' tsriadditive generic
#' @param ... the other arguments
tsriadditive <- function(...) UseMethod("tsriadditive")

#' fit an additive hazards model with two stage residual inclusion method
#' 
#' @param survtime the event time
#' @param cause the indicator records the cause. Default to all one. Zero means right censoring. Greater than
#' or equal to two means other cause.
#' @param treatment the treatment variable, can be null
#' @param IV the instrumental variable
#' @param covariates all the observed confounders
#' @return the fitting result, a class called "tsriadditive"
tsriadditive.default <- function(survtime, cause = NULL, treatment = NULL, IV = NULL, covariates = NULL)
{
  #if there are ties in surtime, in this function we will just break it down randomly
  #cause: 0 as censored, 1 as event at time, >=2 as other cause event
  #treatment can either be binary or continuous, it can also be None, in which case
  #on IV additive model will be fit even if you pass in IV
  #IV can be None, if so, the non IV additive model will be fit
  
  #error message leave for later
  
  #if cause is not provided, we simply assume no censored
  if(is.null(cause)) cause <- rep(1, length(survtime))

  #randomly break ties, order survtime in increasing order
  #reorder all the other variables and save the order into timeorder
  timeorder <- order(survtime, sample(length(survtime)))
  survtime <- survtime[timeorder]
  cause <- cause[timeorder]
  treatment <- treatment[timeorder]
  IV <- IV[timeorder]
  #convert the covariates into matrix object just for later use if it is only dimensional one
  if(is.vector(covariates)) covariates <- matrix(covariates, ncol = 1)
  covariates <- covariates[timeorder, ]
  #comp is the indicator for whether there are competing risks
  comp = !all(cause <= 1)

  #if treatment or IV is None, that means we want to perform regular fit
  if(is.null(treatment) || is.null(IV)) fit <- regadditivefit(survtime, cause, comp, treatment, covariates)
  #so else we have both treatment and IV, we will use 2sri
  else fit <- tsriadditivefit(survtime, cause, comp = comp, treatment, IV, covariates)
  
  fit$byprod$timeorder <- timeorder
  result <- fit
  class(result) <- "tsriadditive"
  result
}

#' the print function associated with our class
#' 
#' @param fit the fitting object after fitting our model
#' @description This function will print our coefficients, the variance covariance matrix of
#' the coeffieients, and the estimate for the baseline hazard function
print.tsriadditive <- function(fit)
{
  cat("\nCoefficients:\n")
  print(fit$coef)
  cat("\ncovariance matrix:\n")
  print(fit$vcov)
  cat("\nbaseline hazard function:\n")
  print(fit$baseline)
}

#' the summary function associated with our class
#' 
#' @param fit the fitting object after fitting our model
#' @description This function will print our coefficients, the variance covariance matrix of
#' the coeffieients, and the corresponding P-values
#' @importFrom stats pnorm
summary.tsriadditive <- function(fit)
{
  result <- data.frame(fit$coef, sqrt(diag(fit$vcov)), fit$coef/sqrt(diag(fit$vcov)), 2 * (1 - pnorm(abs(fit$coef)/sqrt(diag(fit$vcov)))))
  colnames(result) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  result               
}


#' the predict function associated with our class
#' 
#' @param fit the fitting object after fitting our model
#' @param newtreatment the new treatment value
#' @param newIV new instrumental variable value
#' @param newcovariates new observed covariates
#' @return an object recording the corresponding predicted survival curve and corresponding pointwise
#' confidence interval
predict.tsriadditive <- function(fit, newtreatment = NULL, newIV = NULL, newcovariates = NULL)
{
  res = survivalpred(fit, newtreatment, newIV, newcovariates)
  res
}

print.predict.tsriadditive <- function(fit, newtreatment = NULL, newIV = NULL, newcovariates = NULL, ...)
{
  res = predict.tsriadditive(fit, newtreatment, newIV, newcovariates,...)
  cat("The predicted survival function is:\n")
  print(res$survival_pred)
}

#' plot function associated with our class
#' 
#' @param fit the fitting object after fitting our model
#' @param newtreatment the new treatment value
#' @param newIV new instrumental variable value
#' @param newcovariates new observed covariates
#' @param alpha the level for confidence interval
#' @param unit the time unit we focus
#' @param ... the other arguments you want to put in the built-in plot function
#' @description this function will plot the predicted curve and corresponding pointwise confidence
#' interval at level alpha
#' @importFrom stats stepfun 
#' @importFrom graphics plot
plot.tsriadditive <- function(fit, newtreatment = NULL, newIV = NULL, newcovariates = NULL, alpha = 0.05, unit = "", ...)
{
  comp <- fit$byprod$comp
  pred <- survivalpred(fit, newtreatment, newIV, newcovariates)
  survival_pred <- pred$survival_pred
  hazard_pred <- pred$hazard_pred
  newobsz <- pred$newobsz
  survtime <- fit$byprod$survtime
  if(!comp){
    survcurve <- stepfun(fit$byprod$survtime, c(1, survival_pred), right = TRUE)
    plot(survcurve, verticals = TRUE, do.points = FALSE, xlim = c(0, max(survtime)), 
         ylim = c(0, 1), lwd = 1.2, xaxs = 'i', yaxs = 'i',
         xlab = paste0("time (", unit, ")"), ylab = "Survival probability", main = NULL, ...)
    fit <- survprobconfint(hazard_pred, newobsz, fit, alpha)
    survuppercurve <- stepfun(fit$byprod$survtime, c(1, fit$byprod$upper), right = TRUE)
    plot(survuppercurve, verticals = TRUE, do.points = FALSE, add = TRUE, lty = 3,...)
    survlowercurve <- stepfun(fit$byprod$survtime, c(1, fit$byprod$lower), right = TRUE)
    plot(survlowercurve, verticals = TRUE, do.points = FALSE, add = TRUE, lty = 3, ...)
  }
  else{
    distribution_pred <- 1 - survival_pred
    survcurve <- stepfun(fit$byprod$survtime, c(0, distribution_pred), right = TRUE)
    fit <- survprobconfint(hazard_pred, newobsz, fit, alpha)
    survuppercurve <- stepfun(fit$byprod$survtime, c(0, 1 - fit$byprod$upper), right = TRUE)
    survlowercurve <- stepfun(fit$byprod$survtime, c(0, 1 - fit$byprod$lower), right = TRUE)
    plot(survcurve, verticals = TRUE, do.points = FALSE, xlim = c(0, max(survtime)), 
         ylim = c(0, max(1 - fit$byprod$lower)), lwd = 1.2, xaxs = 'i', yaxs = 'i',
         xlab = paste0("time (", unit, ")"), ylab = "Cumulative incidence function", main = NULL, ...) 
    plot(survuppercurve, verticals = TRUE, do.points = FALSE, add = TRUE, lty = 3, ...)
    plot(survlowercurve, verticals = TRUE, do.points = FALSE, add = TRUE, lty = 3, ...)
    fit$byprod <- append(fit$byprod, list(survcurve = survcurve,
                                          survuppercurve = survuppercurve,
                                          survlowercurve = survlowercurve))
  }
  fit
}

