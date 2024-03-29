% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abco.R
\name{t_sampleR_mh}
\alias{t_sampleR_mh}
\title{Sample the threshold parameter}
\usage{
t_sampleR_mh(
  h_yc,
  h_phi,
  h_phi2,
  h_sigma_eta_t,
  h_sigma_eta_0,
  h_st,
  h_r,
  lower_b,
  upper_b,
  omega,
  D
)
}
\arguments{
\item{h_yc}{the \code{T} vector of centered log-volatilities
(i.e., the log-vols minus the unconditional means \code{dhs_mean})}

\item{h_phi}{the \code{1} vector of previous AR(1) coefficient(s)}

\item{h_phi2}{the \code{1} vector of previous penalty coefficient(s)}

\item{h_sigma_eta_t}{the \code{T} vector of log-vol innovation standard deviations}

\item{h_sigma_eta_0}{the \code{1} vector of initial log-vol innovation standard deviations}

\item{h_st}{the \code{T} vector of indicators on whether each time-step exceed the estimated threshold}

\item{h_r}{\code{1} the previous draw of the threshold parameter}

\item{lower_b}{the lower bound in the uniform prior of the threshold variable}

\item{upper_b}{the upper bound in the uniform prior of the threshold variable}

\item{omega}{\code{T} vector of evolution errors}

\item{D}{the degree of differencing (one or two)}
}
\value{
the sampled threshold value \code{r}
}
\description{
Compute one draw of thee threshold parameter in th TAR(1) model with Gaussian innovations
and time-dependent innovation variances. The sampler utilizes metropolis hasting to draw
from uniform prior.
}
