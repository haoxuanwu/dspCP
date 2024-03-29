% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/component_samplers.R
\name{sampleLogVolMu}
\alias{sampleLogVolMu}
\title{Sample the AR(1) unconditional means}
\usage{
sampleLogVolMu(h, h_mu, h_phi, h_sigma_eta_t, h_sigma_eta_0, h_log_scale = 0)
}
\arguments{
\item{h}{the \code{T x p} matrix of log-volatilities}

\item{h_mu}{the \code{p x 1} vector of previous means}

\item{h_phi}{the \code{p x 1} vector of AR(1) coefficient(s)}

\item{h_sigma_eta_t}{the \code{T x p} matrix of log-vol innovation standard deviations}

\item{h_sigma_eta_0}{the standard deviations of initial log-vols}

\item{h_log_scale}{prior mean from scale mixture of Gaussian (Polya-Gamma) prior, e.g. log(sigma_e^2) or dhs_mean0}
}
\value{
a list containing
\itemize{
\item the sampled mean(s) \code{dhs_mean} and
\item the sampled precision(s) \code{dhs_mean_prec_j} from the Polya-Gamma parameter expansion
}
}
\description{
Compute one draw of the unconditional means in an AR(1) model with Gaussian innovations
and time-dependent innovation variances. In particular, we use the sampler for the
log-volatility AR(1) process with the parameter-expanded Polya-Gamma sampler. The sampler also applies
to a multivariate case with independent components.
}
