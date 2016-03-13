\name{mlVARsim}
\alias{mlVARsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates an mlVAR model and dataset
}
\description{
Simulates a lag-1 model and dataset.
}
\usage{
mlVARsim(nPerson = 10, nNode = 5, nTime = 100, sparsity = 0, 
      parRange = c(0.22, 0.4), propPositive = 0.5, diagPositive = TRUE, 
      diagIncluded = TRUE, sdRange = c(0.01, 0.2), shrinkFactor = 0.95, 
      residualStyle = c("full", "diag"), residualShared = TRUE, 
      residualSDrange = c(0.05, 0.1), verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nPerson}{
  Number of persons to simulate
}
  \item{nNode}{
  Number of nodes
}
  \item{nTime}{
  Number of measurements per person
}
  \item{sparsity}{
 Proportion of edges always set to zero
}
  \item{parRange}{
    Range of parameters in the lag-1 network
}
  \item{propPositive}{
    Proportion of positive edges in the lag-1 network
}
  \item{diagPositive}{
    Enforce diagonal of lag-1 to be positive?
}
  \item{diagIncluded}{
    Always include diagonal of lag-1?
}
  \item{sdRange}{
    Range of SD's of parameters
}
  \item{shrinkFactor}{
    Factor to shrink the network with if any eigenvalue is outside the unit circle
}
  \item{residualStyle}{
    Residual covariances?
}
  \item{residualShared}{
    Do all persons share the same residual structure?
}
  \item{residualSDrange}{
 Range of residual SD's
}
  \item{verbose}{
  Logical to output progress
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}