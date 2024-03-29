% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{rtaur_fc}
\alias{rtaur_fc}
\title{rtaur_fc}
\usage{
rtaur_fc(d, cinv, a, b)
}
\arguments{
\item{d}{current value of d}

\item{cinv}{inverse of Gaussian process covariance matrix}

\item{a}{inverse-gamma shape parameter}

\item{b}{inverse-gamma scale parameter}
}
\value{
One sample estimate of tau_r from inverse-gamma distribution
}
\description{
Sample the variance in d
}
