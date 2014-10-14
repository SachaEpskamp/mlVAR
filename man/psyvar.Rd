\name{psyvar}
\alias{psyvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multilevel VAR Estimation for Multiple Time Series
}
\description{
The function \code{psyvar} computes estimates of the multivariate vector autoregression model as introduced by Bringmann et al. (2013) which can be extended through treatment effects, covariates and pre- and post assessment effects. 
}
\usage{
psyvar(data, vars, idvar, dayvar, beepvar, periodvar, lags = 1, treatmentvar, covariates, 
        control = list(optimizer = "bobyqa"))
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
  \item{control}{
 A list of arguments sent to \code{\link[lme4]{lmerControl}}
}

}
\details{
Psyvar has been built to extract individual network dynamics by estimating a multilevel vector autoregression model that models the time dynamics of selected variables both within an individual and on group level. For example, in a lag-1-model each variable at time point t is regressed to a lagged version of itself at time point t-1 and all other variables at time point t-1. In psychological research, for example, this analysis can be used to relate the dynamics of symptoms on one day (as assessed by experience sampling methods) to the dynamics of these symptoms on the consecutive day. }
\value{
Psyvar returns a 'psyvar' object containing
 \item{fixedEffects}{A matrix that contains all fixed effects coefficients with dependent variables as rows and the lagged independent variables as columns.}
 \item{se.fixedEffects}{A matrix that contains all standard errors of the fixed effects.}
 \item{randomEffects}{A list of matrices that contain the random effects coefficients.}
 \item{randomEffectsVariance}{A matrix containing the estimated variances between the random-effects terms}
 \item{pvals}{A matrix that contains p-values for all fixed effects.}
 \item{pseudologlik}{The pseudo log-likelihood.}
 \item{BIC}{Bayesian Information Criterion, i.e. the sum of all univariate models' BICs}
 \item{input}{List containing the names of variables used in the analysis}
}
\references{
Bringmann, L. F., Vissers, N., Wichers, M., Geschwind, N., Kuppens, P., Peeters, F., ... & Tuerlinckx, F. (2013). A network approach to psychopathology: New insights into clinical longitudinal data. PloS one, 8(4), e60188.
}
\author{
Sacha Epskamp (mail@sachaepskamp.com), Marie K. Deserno (m.k.deserno@uva.nl) and Laura F. Bringmann (laura.bringmann@ppw.kuleuven.be)
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
