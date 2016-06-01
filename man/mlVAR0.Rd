\name{mlVAR0}
\alias{mlVAR0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multilevel VAR Estimation for Multiple Time Series
}
\description{
The function \code{mlVAR0} computes estimates of the multivariate vector autoregression model as introduced by Bringmann et al. (2013) which can be extended through treatment effects, covariates and pre- and post assessment effects. 

FUNCTION IS DEPRECATED AND WILL BE REMOVED SOON.
}
\usage{
mlVAR0(data, vars, idvar, lags = 1, dayvar, beepvar,
                 periodvar, treatmentvar, covariates, timevar,
                 maxTimeDiff, control = list(optimizer = "bobyqa"),
                 verbose = TRUE, orthogonal, estimator = c("lmer",
                 "lmmlasso"), method = c("default", "stepwise",
                 "movingWindow"), laginteractions = c("none", "mains",
                 "interactions"), critFun = BIC, lambda = 0,
                 center = c("inSubject","general","none"))
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
  \item{lags}{
Vector indicating the lags to include
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
  \item{treatmentvar}{
Character vector indicating treatment
}
  \item{covariates}{
Character indicating covariates independent of assessment.
}
  \item{timevar}{
  Character indicating the time variable
}
  \item{maxTimeDiff}{
  Maximum time differece to include observation pairs
}
  \item{control}{
 A list of arguments sent to \code{\link[lme4]{lmerControl}}
}

  \item{verbose}{
 Logical to print progress to the console
}

  \item{orthogonal}{
 Logical to indicate if orthogonal estimation (no correlated random effects) should be used. Defaults to \code{FALSE} if the number of nodes is less than 6 and \code{TRUE} otherwise
}

  \item{estimator}{
 Estimator to use. Note: \code{lmmlasso} implementation is very experimental
}

  \item{method}{
 Method to use. Experimental
}

  \item{laginteractions}{
 Experimental, do not use.
}

  \item{critFun}{
 Experimental, do not use.
}

  \item{lambda}{
 lmmlasso lambda parameter
}

  \item{center}{
 Centering to be used. \code{"inSubject"} uses within-person centering, \code{"general"} uses grand-mean centering and \code{"none"} does not use centering. IMPORTANT NOTE: \code{"inSubject"} leads to coefficients to resemble within-person slopes, the other centering option leads to coefficients to be a blend of within and between person slopes.
}
}
\details{
mlVAR0 has been built to extract individual network dynamics by estimating a multilevel vector autoregression model that models the time dynamics of selected variables both within an individual and on group level. For example, in a lag-1-model each variable at time point t is regressed to a lagged version of itself at time point t-1 and all other variables at time point t-1. In psychological research, for example, this analysis can be used to relate the dynamics of symptoms on one day (as assessed by experience sampling methods) to the dynamics of these symptoms on the consecutive day. }
\value{
mlVAR0 returns a 'mlVAR0' object containing
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


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fixedEffects}}, \code{\link{fixedEffects}}
}
\examples{
\dontrun{
### Small network ###
nVar <- 3
nPerson <- 25
nTime <- 25

# Simulate model and data:
Model <- mlVARsim0(nPerson,nVar,nTime,sparsity = 0.5)

# Run mlVAR0:
Res <- mlVAR0(Model)

# Compare true fixed model with significant edges of estimated fixed model:
layout(t(1:2))
plot(Model,"fixed", title = "True model",layout="circle", edge.labels = TRUE)
plot(Res,"fixed", title = "Estimated model", layout = "circle", onlySig = TRUE, 
        alpha = 0.05, edge.labels = TRUE)

# Compare true and estimated individual differences in parameters:
layout(t(1:2))
plot(Model,"fixed", title = "True model",layout="circle", edge.color = "blue", 
        edge.labels = TRUE)
plot(Res,"fixed", title = "Estimated model", layout = "circle", edge.color = "blue", 
        edge.labels = TRUE)

# Compare networks of subject 1:
layout(t(1:2))
plot(Model,"subject",subject = 1, title = "True model",layout="circle", 
        edge.labels = TRUE)
plot(Res,"subject",subject = 1,title = "Estimated model", layout = "circle", 
        edge.labels = TRUE)



### Large network ###
nVar <- 10
nPerson <- 50
nTime <- 50

# Simulate model and data:
Model <- mlVARsim0(nPerson,nVar,nTime, sparsity = 0.5)

# Run orthogonal mlVAR:
Res <- mlVAR0(Model, orthogonal = TRUE)

# Compare true fixed model with significant edges of estimated fixed model:
layout(t(1:2))
plot(Model,"fixed", title = "True model",layout="circle")
plot(Res,"fixed", title = "Estimated model", layout = "circle", onlySig = TRUE, 
        alpha = 0.05)

# Compare true and estimated individual differences in parameters:
layout(t(1:2))
plot(Model,"fixed", title = "True model",layout="circle", edge.color = "blue")
plot(Res,"fixed", title = "Estimated model", layout = "circle", edge.color = "blue")

# Compare networks of subject 1:
layout(t(1:2))
plot(Model,"subject",subject = 1, title = "True model",layout="circle")
plot(Res,"subject",subject = 1,title = "Estimated model", layout = "circle")
}
}
