\name{summary.mlVAR}
\alias{summary.mlVAR}
\alias{print.mlVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of mlVAR results
}
\description{
Prints tables with fit indices and parameter estimates.
}
\usage{
\method{summary}{mlVAR}(object, show = c("fit", "temporal", "contemporaneous", "between"), 
      round = 3, ...)
\method{print}{mlVAR}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An \code{mlVAR} object.
}
  \item{show}{
Which tables to show?
}
  \item{round}{
Number of digits.
}
  \item{x}{
An \code{mlVAR} object.
}
  \item{\dots}{
Not used
}
}
\author{
Sacha Epskamp (mail@sachaepskamp.com), Marie K. Deserno (m.k.deserno@uva.nl) and Laura F. Bringmann (laura.bringmann@ppw.kuleuven.be)
}