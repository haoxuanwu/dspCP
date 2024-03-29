% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_samplers.R
\name{dsp_cp}
\alias{dsp_cp}
\title{MCMC Sampler for Dynamic Shrinkage with Changepoint}
\usage{
dsp_cp(
  y,
  cp = FALSE,
  evol_error = "DHS",
  D = 1,
  useObsSV = TRUE,
  useAnom = FALSE,
  nsave = 1000,
  nburn = 1000,
  nskip = 4,
  mcmc_params = list("mu", "omega", "r"),
  computeDIC = TRUE,
  verbose = TRUE,
  cp_thres = 0.4,
  return_full_samples = TRUE
)
}
\arguments{
\item{y}{the \code{T x 1} vector of time series observations}

\item{cp}{flag indicator to determine whether to use threshold shrinkage with changepoints.}

\item{evol_error}{the evolution error distribution; must be one of
'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)}

\item{D}{degree of differencing (D = 0, D = 1, or D = 2)}

\item{useObsSV}{logical; if TRUE, include a (normal) stochastic volatility model
for the observation error variance}

\item{useAnom}{logical; if TRUE, include an anomaly component in the observation equation
(only for threshold shrinkage with changepoints.)}

\item{nsave}{number of MCMC iterations to record}

\item{nburn}{number of MCMC iterations to discard (burin-in)}

\item{nskip}{number of MCMC iterations to skip between saving iterations,
i.e., save every (nskip + 1)th draw}

\item{mcmc_params}{named list of parameters for which we store the MCMC output;
must be one or more of:
\itemize{
\item "mu" (conditional mean)
\item "yhat" (posterior predictive distribution)
\item "evol_sigma_t2" (evolution error variance)
\item "obs_sigma_t2" (observation error variance)
\item "dhs_phi" (DHS AR(1) coefficient)
\item "dhs_mean" (DHS AR(1) unconditional mean)
}}

\item{computeDIC}{logical; if TRUE, compute the deviance information criterion \code{DIC}
and the effective number of parameters \code{p_d}}

\item{verbose}{logical; should R report extra information on progress?}

\item{cp_thres}{Percentage of posterior samples needed to declare a changepoint}

\item{return_full_samples}{logical; if TRUE, return full MCMC samples for desired parameters;
if FALSE, return the posterior meean for desired parameters}
}
\value{
A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
if threshold shrinkage with changepoints is used, also return detected changepoint locations
}
\description{
Run the MCMC for Bayesian trend filtering with a penalty on zeroth (D = 0),
first (D = 1), or second (D = 2) differences of the conditional expectation.
The user has the option to utilize threshold shrinkage to identify changepoints
for D = 1 or 2.
The penalty is determined by the prior on the evolution errors, which include:
\itemize{
\item the dynamic horseshoe prior ('DHS');
\item the static horseshoe prior ('HS');
\item the Bayesian lasso ('BL');
\item the normal stochastic volatility model ('SV');
\item the normal-inverse-gamma prior ('NIG').
}
In each case, the evolution error is a scale mixture of Gaussians.
Note that only 'DHS' works for threshold shrinkage with changepoints.
Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
mostly relying on a dynamic linear model representation.
}
\note{
The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
deviation is recommended to avoid numerical issues.
}
\examples{
\dontrun{
# Example 1: Change in mean with stochastic volatility
signal = c(rep(0, 50), rep(10, 50))
noise = rep(1, 100)
noise_var = rep(1, 100)
for (k in 2:100){
  noise_var[k] = exp(0.9*log(noise_var[k-1]) + rnorm(1, 0, 0.5))
  noise[k] = rnorm(1, 0, sqrt(noise_var[k])) }

y = signal + noise
mcmc_output = dsp_cp(y, cp=TRUE, mcmc_params = list('yhat', 'mu', "omega", "r"))
cp = mcmc_output$cp
print(paste0('Changepoint Locations: ', cp))
plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, y_true = signal)

# Example 2: Change in mean with Outliers
signal = c(rep(0, 100), rep(5, 100))
y = c(rep(0, 100), rep(5, 100)) + rnorm(200)
y[50] = 10
y[150] = -10

mcmc_output = dsp_cp(y, cp=TRUE, useAnom = TRUE, mcmc_params = list('yhat', 'mu', "omega", "r"))
cp = mcmc_output$cp
plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, y_true = signal)

# Example 3: Change in linear trend
signal = c(seq(1, 50), seq(51, 2))
y = c(seq(1, 50), seq(51, 2)) + rnorm(100)

mcmc_output = dsp_cp(y, cp=TRUE, D=2, mcmc_params = list('yhat', 'mu', "omega", "r"))
cp = mcmc_output$cp
plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, y_true = signal)
}

}
