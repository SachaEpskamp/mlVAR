\name{plot.psyvar}
\alias{plot.psyvar}

\title{
  Plot Method for \code{psyvar}
}
\description{
  The function \code{plot.psyvar} plots estimated model coefficients as a network using \code{qgraph}.
}
\usage{
  \method{plot}{qgraph}(object, type = c("fixed","se","random","subject"), lag = 1,subject,...)
}

\arguments{
  \item{object}{
    A \code{psyvar} object obtained through the \code{\link{psyvar}}-function
  }
  \item{type}{
  Indicates whether to plot a network of fixed effects coefficients ("fixed"), the standard deviations of the fixed effects ("se"), the variances of the random effect terms ("random") or an individual subject's random effects network ("subject").
  }
  \item{lag}{
  Vector indicating the lags to include
  }
  \item{subject}{
  If type="subject", vector indicating the ID subject number
  
  }
  \item{\dots}{
    Arguments sent to \code{\link[qgraph]{qgraph}}
  }
}

\details{
 
}

\author{
  Sacha Epskamp (mail@sachaepskamp.com)
}
