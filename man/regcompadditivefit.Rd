% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regadditivefit.R
\name{regcompadditivefit}
\alias{regcompadditivefit}
\title{regcompadditivefit function}
\usage{
regcompadditivefit(survtime, cause, Z = NULL)
}
\arguments{
\item{survtime}{the event time}

\item{cause}{the indicator records the cause. Default to all one. Zero means right censoring. Greater than
or equal to two means other cause.}

\item{Z}{a variable contains all the regressors}
}
\value{
the fitting result, a list containing the cofficients, the baseline function and
}
\description{
fit an additive hazard without using IV method under competing risks settings
}
