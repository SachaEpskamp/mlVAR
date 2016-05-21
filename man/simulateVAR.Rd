\name{simulateVAR}
\alias{simulateVAR}

\title{
Simulate data from VAR model
}
\description{
Simulates a timeseries using VAR parameters
}
\usage{
simulateVAR(pars,  means = 0, lags = 1, Nt = 100, init, residuals = 0.1,
                 burnin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pars}{
A square matrix or a list of square matrices indicating the VAR parameters
}
\item{means}{
A vector of means.
}
  \item{lags}{
The lags to which the 'pars' argument parameters correspond. If 'pars' is a list then this argument should be a vector indicating which lags are represented by each element of the 'pars' list.
}
  \item{Nt}{
Number of time points
}
  \item{init}{
Initial setup. Must be a matrix of the first lags with rows corresponding to time points and columns corresponding to variables (e.g., if only two lags are used then the matrix must have two rows indicating the first two times points.)
}
  \item{residuals}{
Standard deviation of the residuals or a residual covariance matrix
}
  \item{burnin}{
Initial simulations not returned. Defaults to \code{ min(round(Nt/2), 100)}.
}
}
\author{
  Sacha Epskamp (mail@sachaepskamp.com), Marie K. Deserno (m.k.deserno@uva.nl) and Laura F. Bringmann (laura.bringmann@ppw.kuleuven.be)
}
