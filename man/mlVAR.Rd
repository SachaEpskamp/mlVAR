\name{mlVAR}
\alias{mlVAR}
\alias{predict.mlVAR}
\alias{residuals.mlVAR}
\alias{resimulate}
\alias{resimulate.mlVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multilevel VAR Estimation for Multiple Time Series
}
\description{
The function \code{mlVAR} computes estimates of the multivariate vector autoregression model. This model returns three stuctures: temporal effects (e.g., lag-1 regression weights), contemporaneous relationships (correlations or partial correlations) and between-subject effects (correlations and partial correlations). See details.

The \code{predict} method computes fitted values and residuals from the temporal model (step 1). It returns predictions aligned to the original input data, with \code{NA} for rows that cannot be predicted (e.g., the first observation per person-day lost to lagging, or rows with missing data). When \code{newdata} is supplied, predictions are computed for the new data using the fitted model(s). For \code{estimator = "lmer"}, new subjects are predicted using fixed effects only (\code{allow.new.levels = TRUE}); for \code{estimator = "lm"}, prediction is only possible for subjects seen during training.

The \code{residuals} method is a convenience wrapper around \code{predict} that returns only the residuals.

The \code{resimulate} method generates simulated data from a fitted \code{mlVAR} model using person-specific estimated parameters (fixed + random effects). For each person, a VAR time series is simulated using their estimated temporal effects (Beta), contemporaneous precision matrix (Theta), and person-specific means (mu). The output has the same dimensions and missingness pattern as the original input data, making it suitable for posterior predictive checks and model diagnostics.
}
\usage{
mlVAR(data, vars, idvar, lags = 1, dayvar, beepvar,
                 estimator = c("default", "lmer", "lm","Mplus"),
                 contemporaneous = c("default", "correlated",
                 "orthogonal", "fixed", "unique"), temporal =
                 c("default", "correlated", "orthogonal", "fixed",
                 "unique"), nCores = 1, verbose = TRUE, compareToLags,
                 scale = TRUE, scaleWithin = FALSE, AR = FALSE,
                 MplusSave = TRUE, MplusName = "mlVAR", iterations = "(2000)",
                 chains = nCores, signs, orthogonal,
                 trueMeans, na.rm = TRUE,
                 full_detrend = FALSE
  )

\method{predict}{mlVAR}(object, newdata, scale_back = TRUE,
        include_ids = TRUE, \dots)

\method{residuals}{mlVAR}(object, scale_back = TRUE,
        include_ids = TRUE, \dots)

resimulate(object, \dots)

\method{resimulate}{mlVAR}(object, scale_back = TRUE,
        include_ids = TRUE, nTime = NULL,
        keep_missing = TRUE, \dots)
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
Optional string indicating assessment beep per day. Adding this argument will cause non-consecutive beeps to be treated as missing!
}
  \item{estimator}{
The estimator to be used. \code{"lmer"} for sequential univariate multi-level estimation, \code{"Mplus"} for multivariate Bayesian estimation (requires Mplus), and \code{"lm"} for fixed effects estimation.
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
\item{AR}{Logical, should an auto-regression only model be fitted?}
\item{MplusSave}{Logical, should the Mplus model file and output be saved?}
\item{MplusName}{Name of the Mplus model file and output (without extensions)}
\item{iterations}{The string used to define the number of iterations in Mplus}
\item{chains}{Number of Mplus chains}
\item{signs}{Optional matrix fixing the signs of contemporaneous correlations. Is estimated by running mlVAR with \code{estimator = "lmer"} if missing.}
  \item{orthogonal}{
Deprecated argument only added for backward competability. Ignore.
}
  \item{trueMeans}{
Optional data frame with true person-wise means to use instead of estimated person means. Must contain the \code{idvar} column and all \code{vars} columns.
}
  \item{na.rm}{
Logical, should missing values be removed when checking for rank deficiency and when removing incomplete cases? Defaults to \code{TRUE}.
}
  \item{full_detrend}{
Logical. If \code{TRUE}, standardize each variable per observation position across all subjects before estimation. For example, all first measurements are standardized together, all second measurements together, etc. This removes systematic occasion effects (e.g., time-of-day trends in ESM data). Requires a balanced design (equal number of observations per subject). Applied before grand-mean standardization (\code{scale}). Defaults to \code{FALSE}.
}
  \item{object}{
An \code{mlVAR} object (for \code{predict}, \code{residuals}, and \code{resimulate} methods).
}
  \item{newdata}{
Optional data frame with the same structure as the training data. Must contain the \code{idvar}, \code{dayvar}, \code{beepvar}, and all \code{vars} columns. If supplied, predictions are computed for this new data using the fitted model. For \code{estimator = "lmer"}, new subjects (not seen during training) receive fixed-effects-only predictions. For \code{estimator = "lm"}, prediction is only available for subjects in the training data.
}
  \item{scale_back}{
Logical. If \code{TRUE} (default) and the model was fitted with \code{scale = TRUE}, values are back-transformed to the original scale of the data. If \code{FALSE}, values are returned on the standardized (model) scale. Applies to \code{predict}, \code{residuals}, and \code{resimulate}.
}
  \item{include_ids}{
Logical. If \code{TRUE} (default), the \code{idvar}, \code{dayvar}, and \code{beepvar} columns are prepended as columns in the returned data frame. Applies to \code{predict}, \code{residuals}, and \code{resimulate}. When \code{nTime} is specified in \code{resimulate}, only the \code{idvar} column is included.
}
  \item{nTime}{
Optional integer specifying the number of time points to simulate per person. When \code{NULL} (default), the simulated data has the same number of rows as the original input data. When specified, each person gets exactly \code{nTime} time points, and the output has \code{nIDs * nTime} rows with the \code{idvar} and variable columns only (no \code{dayvar} or \code{beepvar}).
}
  \item{keep_missing}{
Logical. If \code{TRUE} (default) and \code{nTime = NULL}, the original variable-level missingness pattern is applied to the simulated data: any cell that was \code{NA} in the original data is set to \code{NA} in the output. If \code{FALSE}, all valid (non-structural-NA) positions receive simulated values. Ignored when \code{nTime} is specified.
}
}
\details{
This function estimates the multi-level VAR model to obtain temporal, contemporaneous and between-subject effects using nodewise estimation. Temporal and between-subject effects are obtained directly from the models and contemporaneous effects are estimated post-hoc by correlating the residuals. See arxiv.org/abs/1609.04156 for details.

Setting \code{estimator = "Mplus"} will generate a Mplus model, run the analysis and read the results into R. Mplus 8 is required for this estimation. It is recommended to set \code{contemporaneous = "fixed"}, though not required. For the estimation of contemporaneous random effects, the signs of contemporaneous *correlations * (not partial correlations) need be set (or estimated) via the \code{signs} argument.
}
\value{
\code{mlVAR} returns an \code{mlVAR} object.

\code{predict.mlVAR} returns a data frame of predicted (fitted) values with the same number of rows as the original input data (or \code{newdata}). Rows that could not be predicted (e.g., first observation per person-day lost to lagging, or rows with missing data) contain \code{NA}. If \code{include_ids = TRUE}, the \code{idvar}, \code{dayvar}, and \code{beepvar} columns are prepended.

\code{residuals.mlVAR} returns a data frame of residuals with the same structure as \code{predict.mlVAR}.

\code{resimulate.mlVAR} returns a data frame of simulated values. By default (when \code{nTime = NULL}), it has the same number of rows as the original input data, with NAs in the same positions as the original data (unless \code{keep_missing = FALSE}). If \code{include_ids = TRUE}, the \code{idvar}, \code{dayvar}, and \code{beepvar} columns are prepended. When \code{nTime} is specified, the output has \code{nIDs * nTime} rows with the \code{idvar} column and variable columns only. Person-specific stationary means are computed from the model parameters, and innovation covariances are derived from the estimated contemporaneous precision matrices.
}
\references{
Bringmann, L. F., Vissers, N., Wichers, M., Geschwind, N., Kuppens, P., Peeters, F., ... & Tuerlinckx, F. (2013). A network approach to psychopathology: New insights into clinical longitudinal data. PloS one, 8(4), e60188.

Hamaker, E. L., & Grasman, R. P. (2014). To center or not to center? Investigating inertia with a multilevel autoregressive model. Frontiers in psychology, 5.

Epskamp, S., Waldorp, L. J., Mottus, R., & Borsboom, D. (2017). Discovering Psychological Dynamics: The Gaussian Graphical Model in Cross-sectional and Time-series Data. arxiv.org/abs/1609.04156.
}
\author{
Sacha Epskamp (mail@sachaepskamp.com)
}

\seealso{
\code{\link{mlVARcompare}}, \code{\link{summary.mlVAR}}, \code{\link{plot.mlVAR}}, \code{\link{mlGGM}}
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

# Predictions and residuals (same number of rows as input data):
pred <- predict(fit1)
res <- residuals(fit1)
dim(pred)
dim(res)

# Resimulate data from the fitted model (posterior predictive check):
sim <- resimulate(fit1)
dim(sim)  # same dimensions as original data

# Resimulate with custom time series length (e.g., longer series):
sim_long <- resimulate(fit1, nTime = 500)
dim(sim_long)  # nIDs * 500 rows

# Resimulate without applying original missingness pattern:
sim_complete <- resimulate(fit1, keep_missing = FALSE)
sum(is.na(sim_complete))  # only structural NAs (day/beep boundaries)

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

