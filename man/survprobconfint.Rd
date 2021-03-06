% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confintest.R
\name{survprobconfint}
\alias{survprobconfint}
\title{the pointwise confidence interval for the survival curve}
\usage{
survprobconfint(hazard_pred, newobsz, fit = NULL, alpha)
}
\arguments{
\item{hazard_pred}{the predicted hazard function}

\item{newobsz}{the new obtained Z value}

\item{fit}{the fitting object after fitting our model}

\item{alpha}{the prespecified level}
}
\value{
a list containing (1 - alpha) level pointwise confidence interval for the hazard function
}
\description{
the pointwise confidence interval for the survival curve
}
