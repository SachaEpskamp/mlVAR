\name{mlVARsample}
\alias{mlVARsample}
\alias{summary.mlVARsample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulator function given an mlVAR object
}
\description{
Simulates data based on an mlVAR object, estimates the mlVAR network model based on the simulated data and compares the estimated network to the mlVAR object network.  
}
\usage{
mlVARsample(object, nTime = c(25,50,100,200), nSample = 100, pMissing = 0, 
  nReps = 100, nCores = 1, ...)

\method{summary}{mlVARsample}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
mlVAR object, or mlVARsample object in the summary method
}
  \item{nTime}{
Vector with number of time points to test. 
}
  \item{nSample}{
Number of individuals in the dataset. It is possible to decrease the number of individuals compared to the individuals in the mlVAR object. However, it is not possible to have more individuals than there are in the mlVAR object.     
}
  \item{pMissing}{
Percentage of missing data to be simulated.   
}
  \item{nReps}{
Number of repetitions for each condition. 
}
  \item{nCores}{
Number of cores to use. 
}
  \item{\dots}{
Arguments sent to mlVAR. 
}
}

\details{
This function simulates data based on the mlVAR object. The individual networks (random effects) are used to simulate data using the \code{graphicalVARsim} function from the graphicalVAR package (Epskamp, 2020). The individual data is combined into one dataset. This dataset is used to estimate the mlVAR network. 

For every condition, the function returns four values per network comparison measure (correlation, sensitivity, specificity, bias, and precision): one for the fixed temporal effects, one for the fixed contemporaneous effects, the mean comparison value of the random temporal effects, and the mean comparison value of the random contemporaneous effects. 
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\references{
Sacha Epskamp (2020). graphicalVAR: Graphical VAR for Experience Sampling Data. R package version 0.2.3. https://CRAN.R-project.org/package=graphicalVAR
}

\seealso{
\code{\link{mlVARsim}}, \code{\link{mlVAR}} 
}

\examples{
\dontrun{
### Small example ###
# Simulate data:
Model <- mlVARsim(nPerson = 100, nNode = 3, nTime = 50, lag=1)

# Estimate using correlated random effects:
fit <- mlVAR(Model$Data, vars = Model$vars, 
             idvar = Model$idvar, lags = 1, 
             temporal = "correlated")

# Sample from fitted model: 
samples <- mlVARsample(fit, nTime = 50, nSample = 50, pMissing = 0.1,
                       nReps = 5, nCores = 1)

# Summarize results:
summary(samples)
} 
}
