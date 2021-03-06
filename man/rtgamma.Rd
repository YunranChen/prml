% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prml_analysis.R
\name{rtgamma}
\alias{rtgamma}
\title{Generate random variable from truncated gamma}
\usage{
rtgamma(size, shape_, rate_, a, b)
}
\arguments{
\item{size}{A number, size of the random variable.}

\item{shape_}{A number, shape of the truncated gamma.}

\item{rate_}{A number, rate of the truncated gamma.}

\item{a}{lower bound of the range}

\item{b}{upper bound of the range}
}
\value{
A vector, iid samples from \eqn{ Ga_{[a,b]}(shape\_,rate\_) }
}
\description{
Generate random variable from truncated gamma
}
