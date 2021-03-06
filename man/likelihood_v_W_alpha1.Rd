% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/closeform.R
\name{likelihood_v_W_alpha1}
\alias{likelihood_v_W_alpha1}
\title{likelihood p(v | W)}
\usage{
likelihood_v_W_alpha1(v, W, beta, log = TRUE, vectorizable = TRUE)
}
\arguments{
\item{v}{a column of V where every component is a feature.}
}
\value{
log probability
}
\description{
computes p(v|W) using the simplified expression (if alpha=1)
and F = 1
}
\details{
Computes the likelihood of a column of V. Since each column
is independent, the likelihood of the whole matrix V is the product of the
likelihoods of each column in V.
}
