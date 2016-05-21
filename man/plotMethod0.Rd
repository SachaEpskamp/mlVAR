\name{plot.mlVAR0}
\alias{plot.mlVAR0}
\alias{plot.mlVARsim0}

\title{
  Plot Method for mlVAR0
}
\description{
  The function \code{plot.mlVAR0} plots estimated model coefficients as a network using \code{qgraph}.
  
  FUNCTION IS DEPRECATED AND WILL BE REMOVED SOON.
}
\usage{
  \method{plot}{mlVAR0}(x, type = c("fixed", "SD", "subject"), lag = 1,
                 subject, order, onlySig = FALSE, alpha, ...)
  \method{plot}{mlVARsim0}(x, type = c("fixed", "SD", "subject"), lag = 1,
                 subject, order,  ...)
}

\arguments{
  \item{x}{
    A \code{mlVAR0} object obtained through the \code{\link{mlVAR0}}-function
  }
  \item{type}{
  Indicates whether to plot a network of fixed effects coefficients ("fixed"), the standard deviations of the random effect terms ("SD") or an individual subject's random effects network ("subject").
  }
  \item{lag}{
  Vector indicating the lags to include
  }
  \item{subject}{
  If type="subject", vector indicating the ID subject number
  
  }
    \item{order}{
    Order of nodes
  }
    \item{onlySig}{
    Logical. Set to \code{TRUE} to only plot significant fixed effects.
  }
      \item{alpha}{
    Significance level to test edges at if \code{onlySig == TRUE}. Defaults to Bonferonni corrected alpha level of 0.05 divided by the number of fixed effects.
  }
  \item{\dots}{
    Arguments sent to \code{\link[qgraph]{qgraph}}
  }
}


\author{
Sacha Epskamp (mail@sachaepskamp.com), Marie K. Deserno (m.k.deserno@uva.nl) and Laura F. Bringmann (laura.bringmann@ppw.kuleuven.be)
}
