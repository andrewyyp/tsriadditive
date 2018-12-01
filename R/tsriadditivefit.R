#' fit an additive hazard using two stage residual inclusion method
#' 
#' @param survtime the event time
#' @param cause the indicator records the cause. Default to all one. Zero means right censoring. Greater than
#' or equal to two means other cause.
#' @param comp the indicator of whether modeling subdistribution hazard
#' @param treatment the treatment variable, can be null
#' @param IV the instrumental variable
#' @param covariates all the observed confounders
#' @return the fitting result, a list containing the cofficients, the baseline function, the variance covariance
#' @importFrom stats glm residuals lm
tsriadditivefit <- function(survtime, cause, comp = FALSE, treatment = NULL, IV = NULL, covariates = NULL)
{
  #first stage
  binary <- treatment == 0 || treatment == 1
  #if binary is TRUE we treatment as binary treatment o/w continuous
  if(binary){
    firstfit <- glm(treatment ~. , data = data.frame(cbind(treatment, IV, covariates)), family = "binomial")
    residual <- residuals(firstfit, type = "response")
  }
  else{
    firstfit <- lm(treatment ~. , data = data.frame(cbind(treatment, IV, covariates)))
    residual <- residuals(firstfit, type = "response")
  } 
  X <- cbind(1, IV, covariates)
  Z <- cbind(treatment, residual, covariates)
  if(!comp) fit <- tsrisurvadditivefit(survtime, cause, Z)
  else fit <- tsricompadditivefit(survtime, cause, Z)
  firststage <- list(X = X, firstfit = firstfit)
  fit$byprod <- append(fit$byprod, firststage)
  fit$byprod$binary <- binary
  fit <- betavarest(fit)
  fit
}

#' fit an additive hazard using IV method under survival settings
#' 
#' @param survtime the event time
#' @param cause the indicator records the cause. Default to all one. Zero means right censoring. Greater than
#' or equal to two means other cause.
#' @param Z a variable contains all the regressors
#' @return the fitting result, a list containing the cofficients, the baseline function and 
#' the byproduct including some pieces during the computing process
tsrisurvadditivefit <- function(survtime, cause = NULL, Z)
{
  fit <- regsurvadditivefit(survtime, cause, Z)
  fit$byprod$useIV = TRUE
  fit
}

#' fit an additive hazard without using IV method under competing risks settings
#' 
#' @param survtime the event time
#' @param cause the indicator records the cause. Default to all one. Zero means right censoring. Greater than
#' or equal to two means other cause.
#' @param Z a variable contains all the regressors
#' @return the fitting result, a list containing the cofficients, the baseline function and 
#' the byproduct including some pieces during the computing process
tsricompadditivefit <- function(survtime, cause = NULL, Z)
{
  fit <- regcompadditivefit(survtime, cause, Z)
  fit$byprod$useIV <- TRUE
  fit
}