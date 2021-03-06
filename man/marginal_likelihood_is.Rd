% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihood.R
\name{marginal_likelihood_is}
\alias{marginal_likelihood_is}
\title{Marginal likelihood with Importance Sampling}
\usage{
marginal_likelihood_is(V, W, alpha = 1, beta = 1, nsamples = 100)
}
\arguments{
\item{V}{data matrix}

\item{W}{dictionary matrix}

\item{alpha}{parameter for the Negative Multinomial}

\item{beta}{parameter for the Negative Multinomial}

\item{nsamples}{number of samples from the Importance distribution}
}
