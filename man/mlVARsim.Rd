\name{mlVARsim}
\alias{mlVARsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates an mlVAR model and data
}
\description{
Simulates an mlVAR model and data with a random variance-covariance matrix for the random effects.
}
\usage{
mlVARsim(nPerson = 10, nNode = 5, nTime = 100, lag = 1, thetaVar = rep(1,nNode),
DF_theta = nNode * 2, mu_SD = c(1, 1), init_beta_SD = c(0.1, 1), fixedMuSD = 1, 
shrink_fixed = 0.9, shrink_deviation = 0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nPerson}{
Number of subjects
}
  \item{nNode}{
Number of variables
}
  \item{nTime}{
Number of observations per person
}
  \item{lag}{
The maximum lag to be used
}
\item{thetaVar}{
  Contemporaneous fixed effect variances
}
  \item{DF_theta}{
Degrees of freedom in simulating person-specific contemporaneous covariances (e.g., the individual differences in contemporaneous effects)
}
  \item{mu_SD}{
Range of standard deviation for the means
}
  \item{init_beta_SD}{
Initial range of standard deviations for the temporal effects
}
  \item{fixedMuSD}{
Standard deviation used in sampling the fixed effects
}
  \item{shrink_fixed}{
Shrinkage factor for shrinking the fixed effects if the VAR model is not stationary
}
  \item{shrink_deviation}{
Shrinkage factor for shrinking the random effects variance if the VAR model is not stationary
}
}
\author{
Sacha Epskamp (mail@sachaepskamp.com)
}

