% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_gibbs.R
\name{sample_gibbs_c}
\alias{sample_gibbs_c}
\title{Conditional sampler of c_fn}
\usage{
sample_gibbs_c(v_fn, w_f, h_n)
}
\arguments{
\item{v_fn}{value of V[f,n] the drawn vector will sum up to this value}

\item{w_f}{f-th row of the W matrix}

\item{h_n}{h-th column of the H matrix}
}
\description{
Samples a multinomial vector c_fn
}
\details{
Deals with cases where v_fn = 0, therefore c_fn = 0,..0
}
