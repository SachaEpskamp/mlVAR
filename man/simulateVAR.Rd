\name{simulateVAR}
\alias{simulateVAR}

\title{
Simulate data from VAR model
}
\description{
Simulates a timeseries using VAR parameters
}
\usage{
simulateVAR(pars, lags = 1, Nt, init, errorSD = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pars}{
A square matrix or a list of square matrices indicating the VAR parameters
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
  \item{errorSD}{
Standard deviation of the error
}
}
\author{
  Sacha Epskamp (mail@sachaepskamp.com), Marie K. Deserno (m.k.deserno@uva.nl) and Laura F. Bringmann (laura.bringmann@ppw.kuleuven.be)
}
