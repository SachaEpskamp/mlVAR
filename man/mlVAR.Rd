\name{mlVAR}
\alias{mlVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multilevel VAR Estimation for Multiple Time Series
}
\description{
The function \code{mlVAR0} computes estimates of the multivariate vector autoregression model. This model returns three stuctures: temporal effects (e.g., lag-1 regression weights), contemporaneous relationships (correlations or partial correlations) and between-subject effects (correlations and partial correlations). See details.
}
\usage{
mlVAR(data, vars, idvar, lags = 1, dayvar, beepvar,
                 estimator = c("default", "lmer", "lm"),
                 contemporaneous = c("default", "correlated",
                 "orthogonal", "fixed", "unique"), temporal =
                 c("default", "correlated", "orthogonal", "fixed",
                 "unique"), nCores = 1, verbose = TRUE, compareToLags,
                 scale = TRUE, scaleWithin = FALSE, orthogonal)
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
String indicating assessment day. Adding this argument makes sure that the first measurement of a day is not regressed on the last measurement of the previous day. IMPORTANT: only add this if the data has multiple observations per day.
}
  \item{beepvar}{
String indicating assessment beep per day (if missing, is added). Adding this argument will cause non-consecutive beeps to be treated as missing!
}
  \item{estimator}{
The estimator to be used. Currently only \code{"lmer"} is supported.
}
\item{contemporaneous}{
  How should the contemporaneous networks be estimated? These networks are always estimated post-hoc by investigating the residuals of the temporal models. \code{"correlated"} and \code{"orthogonal"} run second multi-level models in which the networks are estimated using node-wise estimation. \code{"fixed"} and \code{"unique"} simply correlate the residuals, either by computing one network for all subjects (fixed) or a single network per per subject.
}
  \item{temporal}{
How should the temporal effects be estimated? \code{"correlated"} estimates correlated random effects, \code{"orthogonal"} estimates non-correlated random effects and \code{"fixed"} estimates a model in which only the intercept is random. Defaults to \code{"correlated"} when the number of variables is less than 6 and \code{"orthogonal"} otherwise. \code{"unique"} uses \code{lm} to estimate an unique model for each subject.
}
\item{nCores}{Number of cores to use in computation}
  \item{verbose}{
Logical indicating if console messages and the progress bar should be shown.
}
\item{scale}{Logical, should variables be standardized before estimation?}
\item{scaleWithin}{Logial, should variables be scaled within-person (set to \code{FALSE} to only center within-person)}
  \item{compareToLags}{
A vector indicating which lags to base the data on. If the model is to be compared with a model with multiple lags using \code{\link{mlVARcompare}}, this argument must be used to make sure the number of observations is the same in both models (e.g., a lag 1 model can model the second observation of a day and a lag-2 model can't, causing different number of observations and incomparable models). It is suggested to not use this argument unless you want to compare models, and always run \code{mlVAR} without using this argument afterwards in the selected model.
}
  \item{orthogonal}{
Deprecated argument only added for backward competability. Ignore.
}
}
\details{
This function estimates the multi-level VAR model to obtain temporal, contemporaneous and between-subject effects using nodewise estimation. Temporal and between-subject effects are obtained directly from the models and contemporaneous effects are estimated post-hoc by correlating the residuals.

For more information on multi-level VAR and centering please see Hamaker and Grasman (2014) and Bringmann et al., (2013). A draft on the paper describing this package and metho can be obtained from the maintainer (mail@sachaepskamp.com).
}
\value{
An \code{mlVAR} object 
}
\references{
Bringmann, L. F., Vissers, N., Wichers, M., Geschwind, N., Kuppens, P., Peeters, F., ... & Tuerlinckx, F. (2013). A network approach to psychopathology: New insights into clinical longitudinal data. PloS one, 8(4), e60188.

Hamaker, E. L., & Grasman, R. P. (2014). To center or not to center? Investigating inertia with a multilevel autoregressive model. Frontiers in psychology, 5.
}
\author{
Sacha Epskamp (mail@sachaepskamp.com)
}

\seealso{
\code{\link{mlVARcompare}}, \code{\link{summary.mlVAR}}, \code{\link{plot.mlVAR}}
}
\examples{
\dontrun{
### Small example ###
# Simulate data:
Model <- mlVARsim(nPerson = 50, nNode = 3, nTime = 50, lag=1)

# Estimate using correlated random effects:
fit1 <- mlVAR(Model$Data, vars = Model$vars, idvar = Model$idvar, lags = 1, temporal = "correlated")

# Print some pointers:
print(fit1)

# Summary of all parameter estimates:
summary(fit1)

# Compare temporal relationships:
layout(t(1:2))
plot(Model, "temporal", title = "True temporal relationships", layout = "circle")
plot(fit1, "temporal", title = "Estimated temporal relationships", layout = "circle")

# Compare contemporaneous partial correlations:
layout(t(1:2))
plot(Model, "contemporaneous", title = "True contemporaneous relationships", 
    layout = "circle")
plot(fit1, "contemporaneous", title = "Estimated contemporaneous relationships", 
    layout = "circle")

# Compare between-subjects partial correlations:
layout(t(1:2))
plot(Model, "between", title = "True between-subjects relationships", layout = "circle")
plot(fit1, "between", title = "Estimated between-subjects relationships", 
    layout = "circle")

# Run same model with non-correlated temporal relationships and fixed-effect model:
fit2 <- mlVAR(Model$Data, vars = Model$vars, idvar = Model$idvar, lags = 1, 
    temporal = "orthogonal")
fit3 <- mlVAR(Model$Data, vars = Model$vars, idvar = Model$idvar, lags = 1, 
    temporal = "fixed")

# Compare models:
mlVARcompare(fit1,fit2,fit3)

# Inspect true parameter correlation matrix:
Model$model$Omega$cor$mean
# Even though correlations are high, orthogonal model works well often!


### Large example ###
Model <- mlVARsim(nPerson = 100, nNode = 10, nTime = 100,lag=1)

# Correlated random effects no longer practical. Use orthogonal or fixed:
fit4 <- mlVAR(Model$Data, vars = Model$vars, idvar = Model$idvar, lags = 1, 
    temporal = "orthogonal")
fit5 <- mlVAR(Model$Data, vars = Model$vars, idvar = Model$idvar, lags = 1, 
    temporal = "fixed")

# Compare models:
mlVARcompare(fit4, fit5)

# Compare temporal relationships:
layout(t(1:2))
plot(Model, "temporal", title = "True temporal relationships", layout = "circle")
plot(fit4, "temporal", title = "Estimated temporal relationships", layout = "circle")

# Compare contemporaneous partial correlations:
layout(t(1:2))
plot(Model, "contemporaneous", title = "True contemporaneous relationships", 
    layout = "circle")
plot(fit4, "contemporaneous", title = "Estimated contemporaneous relationships", 
    layout = "circle")

# Compare between-subjects partial correlations:
layout(t(1:2))
plot(Model, "between", title = "True between-subjects relationships", layout = "circle")
plot(fit4, "between", title = "Estimated between-subjects relationships", 
    layout = "circle")
}
}

