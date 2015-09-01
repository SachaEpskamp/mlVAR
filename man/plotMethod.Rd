\name{plot.mlVAR}
\alias{plot.mlVAR}

\title{
  Plot Method for mlVAR
}
\description{
  The function \code{plot.mlVAR} plots estimated model coefficients as a network using \code{qgraph}.
}
\usage{
  \method{plot}{mlVAR}(x, type = c("fixed","se","random","subject"), lag = 1,subject,...)
}

\arguments{
  \item{x}{
    A \code{mlVAR} object obtained through the \code{\link{mlVAR}}-function
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


\author{
Sacha Epskamp (mail@sachaepskamp.com), Marie K. Deserno (m.k.deserno@uva.nl) and Laura F. Bringmann (laura.bringmann@ppw.kuleuven.be)
}
