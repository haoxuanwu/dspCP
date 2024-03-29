% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abco.R
\name{abco}
\alias{abco}
\title{Adaptive Bayesian Changepoint with Outliers}
\usage{
abco(
  y,
  D = 1,
  useAnom = TRUE,
  useObsSV = TRUE,
  nsave = 1000,
  nburn = 1000,
  nskip = 4,
  mcmc_params = list("mu", "omega", "r"),
  verbose = TRUE,
  cp_thres = 0.5
)
}
\arguments{
\item{y}{the \code{T} vector of time series observations}

\item{D}{degree of differencing (D = 1, or D = 2)}

\item{useAnom}{logical; if TRUE, include an anomaly component in the observation equation}

\item{useObsSV}{logical; if TRUE, include a (normal) stochastic volatility model
for the observation error variance}

\item{nsave}{number of MCMC iterations to record}

\item{nburn}{number of MCMC iterations to discard (burnin)}

\item{nskip}{number of MCMC iterations to skip between saving iterations,
i.e., save every (nskip + 1)th draw}

\item{mcmc_params}{named list of parameters for which we store the MCMC output;
must be one or more of:
\itemize{
\item "mu" (conditional mean)
\item "yhat" (posterior predictive distribution)
\item "evol_sigma_t2" (evolution error variance)
\item "obs_sigma_t2" (observation error variance)
\item "zeta_sigma_t2" (outlier error variance)
\item "dhs_phi" (DHS AR(1) coefficient)
\item "dhs_mean" (DHS AR(1) unconditional mean)
}}

\item{verbose}{logical; should R report extra information on progress?}

\item{cp_thres}{Percentage of posterior samples needed to declare a changepoint}
}
\value{
A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
}
\description{
Run the MCMC sampler for ABCO with a penalty on
first (D = 1), or second (D = 2) differences of the conditional expectation.
The penalty utilizes the dynamic horseshoe prior on the evolution errors.
Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
mostly relying on a dynamic linear model representation.
}
