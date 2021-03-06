% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tsriadditive.R
\name{tsriadditive.default}
\alias{tsriadditive.default}
\title{fit an additive hazards model with two stage residual inclusion method}
\usage{
\method{tsriadditive}{default}(survtime, cause = NULL,
  treatment = NULL, IV = NULL, covariates = NULL)
}
\arguments{
\item{survtime}{the event time}

\item{cause}{the indicator records the cause. Default to all one. Zero means right censoring. Greater than
or equal to two means other cause.}

\item{treatment}{the treatment variable, can be null}

\item{IV}{the instrumental variable}

\item{covariates}{all the observed confounders}
}
\value{
the fitting result, a class called "tsriadditive"
}
\description{
fit an additive hazards model with two stage residual inclusion method
}
