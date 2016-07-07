\name{plot.mlVAR}
\alias{plot.mlVAR}
\alias{plot.mlVARsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Method for mlVAR
}
\description{
  The function \code{plot.mlVAR} plots estimated model coefficients as networks using \code{qgraph}. These can be three networks: temporal, contemporaneous and between-subjects effects, of which the latter two can be plotted as a correlation or a partial correlation network.
}
\usage{
  \method{plot}{mlVAR}(x, type = c("temporal", "contemporaneous", "between"), lag = 1, 
      partial = TRUE, SD = FALSE, subject, order, nonsig = c("show", "dashed", "hide"), 
      rule = c("or","and"), alpha = 0.05, onlySig = FALSE, layout = "spring", ...)
  \method{plot}{mlVARsim}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An \code{mlVAR} object.
}
  \item{type}{
What network to plot?
}
  \item{lag}{
The lag to use when \code{type = "temporal"}
}
  \item{partial}{
Logical, should partial correlation matrices be plotted instead of correlation methods? Only used if \code{type} is \code{"contemporaneous"} or \code{"between"}. Defaults to \code{TRUE}.
}
  \item{SD}{
Logical. Plot the standard-deviation of random effects instead of the fixed effect estimate?
}
  \item{subject}{
Subject number. If not missing, will plot the network of a specific subject instead.
}
  \item{order}{
An optional character vector used to set the order of nodes in the network.
}
  \item{nonsig}{
How to handle non-significant edges? 
}
\item{rule}{
How to choose significance in node-wise estimated GGMs (contemporaneous and between-subjects). \code{"or"} selects an edge as being significant if one node predicting the other is significant, and \code{"and"} requires both predictions to be significant.
}
  \item{alpha}{
Alpha level to test for significance
}
  \item{onlySig}{
Deprecated argument only used for backward competability.
}
  \item{layout}{
The layout argument used by \code{\link[qgraph]{qgraph}}
}
  \item{\dots}{
Arguments sent to \code{\link[qgraph]{qgraph}}
}
}

\author{
Sacha Epskamp (mail@sachaepskamp.com)
}
