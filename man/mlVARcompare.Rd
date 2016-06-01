\name{mlVARcompare}
\alias{mlVARcompare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare mlVAR model fit
}
\description{
This function compares the fit of several mlVAR models. Since an mlVAR model is a 
combination of univariate models this function will compare the fits for each univariate model.
}
\usage{
mlVARcompare(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
Any number of objects obtained from \code{\link{mlVAR}}
}
}
\details{
Important to note is that the number of observations must be equal to make models comparable. If the lags are different and \code{compareToLags} was not used in \code{mlVAR} this function will stop with an informative error message.
}

\author{
Sacha Epskamp (mail@sachaepskamp.com)
}
\examples{
\dontrun{
### Small example ###
# Simulate data:
Model <- mlVARsim(nPerson = 50, nNode = 3, nTime = 50, lag=1)

# Estimate using different methods:
fit1 <- mlVAR(Model$Data, vars = Model$vars, idvar = Model$idvar, lags = 1, 
    temporal = "correlated")
fit2 <- mlVAR(Model$Data, vars = Model$vars, idvar = Model$idvar, lags = 1, 
    temporal = "orthogonal")
fit3 <- mlVAR(Model$Data, vars = Model$vars, idvar = Model$idvar, lags = 1, 
    temporal = "fixed")

# Compare models:
mlVARcompare(fit1,fit2,fit3)
}
}
