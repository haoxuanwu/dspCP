% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{initSV}
\alias{initSV}
\title{Initialize the stochastic volatility parameters}
\usage{
initSV(omega)
}
\arguments{
\item{omega}{\code{T x p} matrix of errors}
}
\value{
List of relevant components: \code{sigma_wt}, the \code{T x p} matrix of standard deviations,
and additional parameters (unconditional mean, AR(1) coefficient, and standard deviation).
}
\description{
Compute initial values for normal stochastic volatility parameters.
The model assumes an AR(1) for the log-volatility.
}
