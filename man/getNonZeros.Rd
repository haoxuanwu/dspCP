% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{getNonZeros}
\alias{getNonZeros}
\title{Compute Non-Zeros (Signals)}
\usage{
getNonZeros(post_evol_sigma_t2, post_obs_sigma_t2 = NULL)
}
\arguments{
\item{post_evol_sigma_t2}{the \code{Nsims x T} or \code{Nsims x T x p}
matrix/array of posterior draws of the evolution error variances.}

\item{post_obs_sigma_t2}{the \code{Nsims x 1} or \code{Nsims x T} matrix of
posterior draws of the observation error variances.}
}
\value{
A vector (or matrix) of indices identifying the signals according to the
horsehoe-type thresholding rule.
}
\description{
Estimate the location of non-zeros (signals) implied by
horseshoe-type thresholding.
}
\details{
Thresholding is based on \code{kappa[t] > 1/2}, where
\code{kappa = 1/(1 + evol_sigma_t2/obs_sigma_t2)}, \code{evol_sigma_t2} is the
evolution error variance, and \code{obs_sigma_t2} is the observation error variance.
In particular, the decision rule is based on the posterior mean of \code{kappa}.
}
\note{
The thresholding rule depends on whether the prior variance for the state
variable \code{mu} (i.e., \code{evol_sigma_t2}) is scaled by the observation standard
deviation, \code{obs_sigma_t2}. Explicitly, if \code{mu[t]} ~ N(0, \code{evol_sigma_t2[t]})
then the correct thresholding rule is based on \code{kappa = 1/(1 + evol_sigma_t2/obs_sigma_t2)}.
However, if \code{mu[t]} ~ N(0, \code{evol_sigma_t2[t]*obs_sigma_t2[t]})
then the correct thresholding rule is based on \code{kappa = 1/(1 + evol_sigma_t2)}.
The latter case may be implemented by omitting the input for \code{post_obs_sigma_t2}
(or setting it to NULL).
}
\examples{
\dontrun{
# Simulate a function with many changes:
simdata = simUnivariate(signalName = "blocks", T = 128, RSNR = 7, include_plot = TRUE)
y = simdata$y

# Run the MCMC:
out = btf(y, D = 1, evol_error = "HS",
                 mcmc_params = list('mu','evol_sigma_t2', 'obs_sigma_t2'))
# Compute the CPs:
nz = getNonZeros(post_evol_sigma_t2 = out$evol_sigma_t2,
                post_obs_sigma_t2 = out$obs_sigma_t2)
# True CPs:
cp_true = 1 + which(abs(diff(simdata$y_true)) > 0)

# Plot the results:
plot_cp(y, nz)
plot_cp(colMeans(out$mu), nz)
# abline(v = cp_true)

# Regression example:
simdata = simRegression(T = 200, p = 5, p_0 = 2)
y = simdata$y; X = simdata$X
# Run the MCMC:
out = btf_reg(y, X, D = 1, evol_error = 'DHS',
                     mcmc_params = list('mu', 'beta', 'yhat',
                                        'evol_sigma_t2', 'obs_sigma_t2'))
for(j in 1:ncol(X))
 plot_fitted(rep(0, length(y)),
             mu = colMeans(out$beta[,,j]),
             postY = out$beta[,,j],
             y_true = simdata$beta_true[,j])

# Compute the CPs
nz = getNonZeros(post_evol_sigma_t2 = out$evol_sigma_t2,
                post_obs_sigma_t2 = out$obs_sigma_t2)
for(j in 1:ncol(X))
 plot_cp(mu = colMeans(out$beta[,,j]),
         cp_inds = nz[nz[,2]==j,1])
}

}
