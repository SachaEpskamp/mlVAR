\name{psyvar}
\alias{psyvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multilevel VAR estimation for multiple time series
}
\description{
Psyvar is based on a method introduced by Bringmann et al. (2013). It has been built to extract individual network dynamics by estimating a multilevel vector autoregression model that models the time dynamics of selected variables both within an individual and on group level. Each variable (node in the network) at time point $t$ is regressed to a lagged version of itself at time point $t-1$ and all other variables at time point $t-1$. In psychological research, for example, this analysis can be used to relate the dynamics of symptoms on one day to the dynamics of these symptoms on the consecutive day. 
}
\usage{
psyvar(data, vars, idvar, dayvar, beepvar, periodvar, lags = 1, treatmentvar, covariates, maxiter = 800)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{data}{
Data frame
}
  \item{vars}{
Vectors of variables to include in the analysis
}
  \item{idvar}{
String indicating the subject ID
}
  \item{dayvar}{
String indicating assessment day (if missing, every assessment is set to one day)
}
  \item{beepvar}{
String indicating assessment beep per day (if missing, is added)
}
  \item{periodvar}{
String indicating the period (baseline, treatment period, etc.) of assessment (if missing, every assessment is set to one period)
}
  \item{lags}{
Vector indicating the lags to include
}
  \item{treatmentvar}{
Character vector indicating treatment
}
  \item{covariates}{
Character indicating covariates independent of assessment.
}
  \item{maxiter}{
Parameter for the maximal number of iterations passed to the optimizer
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
