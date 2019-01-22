\name{mlVARsample}
\alias{mlVARsample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulator function given an mlVAR object
}
\description{
Simulates data based on an mlVAR object
}
\usage{
mlVARsample(object, nTime = c(25, 50, 100, 200), nReps = 100, nCores = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
mlVAR object
}
  \item{nTime}{
Vector with number of time points to test
}
  \item{nReps}{
Number of repititions for each condition
}
  \item{nCores}{
Number of cores to use
}
  \item{\dots}{
Arguments sent to mlVAR
}
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

