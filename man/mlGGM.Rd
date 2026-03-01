\name{mlGGM}
\alias{mlGGM}
\alias{print.mlGGM}
\alias{summary.mlGGM}
\alias{plot.mlGGM}
\title{
Multi-Level Gaussian Graphical Model
}
\description{
Estimates within-cluster and between-cluster partial correlation networks from cross-sectional multilevel data (e.g., people nested in clusters) using a single-step nodewise regression approach. For each variable, a mixed-effects model is fitted with within-cluster centered predictors and cluster-mean predictors. The within-cluster and between-cluster GGMs are then extracted from the regression coefficients.

NOTE: This function is currently experimental and may change in future versions.
}
\usage{
mlGGM(data, vars, idvar = "id",
      estimator = c("lmer"),
      randomeffects = c("default", "correlated", "orthogonal", "fixed"),
      scale = TRUE, na.rm = TRUE, verbose = TRUE)

\method{print}{mlGGM}(x, \dots)
\method{summary}{mlGGM}(object, show = c("fit", "within", "between"),
        round = 3, \dots)
\method{plot}{mlGGM}(x, type = c("within", "between"),
     partial = TRUE, SD = FALSE, subject, order,
     nonsig = c("default", "show", "hide", "dashed"),
     rule = c("or", "and"), alpha = 0.05,
     layout = "spring", verbose = TRUE, \dots)
}
\arguments{
  \item{data}{
A data frame containing the variables and cluster ID.
}
  \item{vars}{
A character vector of variable names to include in the analysis.
}
  \item{idvar}{
A string indicating the cluster/group ID variable name.
}
  \item{estimator}{
The estimator to be used. Currently only \code{"lmer"} is supported.
}
  \item{randomeffects}{
How should random effects be estimated? \code{"correlated"} estimates correlated random slopes for within-cluster predictors, \code{"orthogonal"} estimates uncorrelated random slopes, and \code{"fixed"} estimates only a random intercept (no random slopes). Defaults to \code{"correlated"} when the number of variables is 6 or fewer and \code{"orthogonal"} otherwise.
}
  \item{scale}{
Logical, should variables be grand-mean standardized before estimation?
}
  \item{na.rm}{
Logical, should rows with missing values be removed?
}
  \item{verbose}{
Logical indicating if console messages and progress bars should be shown.
}
  \item{x}{
An \code{mlGGM} object.
}
  \item{object}{
An \code{mlGGM} object.
}
  \item{show}{
Character vector indicating which sections to show in the summary. Options: \code{"fit"}, \code{"within"}, \code{"between"}.
}
  \item{round}{
Number of decimal places to round to in the summary output.
}
  \item{type}{
The type of network to plot: \code{"within"} for the within-cluster network, \code{"between"} for the between-cluster network.
}
  \item{partial}{
Logical. If \code{TRUE} (default), plots partial correlations. If \code{FALSE}, plots zero-order correlations.
}
  \item{SD}{
Logical. If \code{TRUE}, plots random effect standard deviations instead of fixed effects. Only available for the within-cluster network.
}
  \item{subject}{
Integer specifying which cluster's network to plot (within-cluster only).
}
  \item{order}{
Character vector or numeric vector specifying the order of nodes in the plot.
}
  \item{nonsig}{
How to handle non-significant edges: \code{"show"} shows all edges, \code{"hide"} removes non-significant edges, \code{"dashed"} shows them as dashed lines. \code{"default"} hides non-significant edges for partial correlation networks and shows them otherwise.
}
  \item{rule}{
Significance rule for undirected networks: \code{"or"} requires at least one direction to be significant, \code{"and"} requires both directions to be significant.
}
  \item{alpha}{
Significance level for thresholding edges.
}
  \item{layout}{
Layout algorithm for the network plot (passed to \code{\link[qgraph]{qgraph}}).
}
  \item{\dots}{
Additional arguments passed to \code{\link[qgraph]{qgraph}} (for plot) or ignored.
}
}
\details{
For each variable \eqn{y_i}, the following mixed-effects model is fitted:

\deqn{y_i \sim 1 + \text{centered}_j + \text{mean}_j + (1 + \text{centered}_j | \text{cluster})}

where \eqn{\text{centered}_j = y_j - \bar{y}_j^{(\text{cluster})}} are within-cluster centered values and \eqn{\text{mean}_j = \bar{y}_j^{(\text{cluster})}} are cluster means, for all \eqn{j \neq i}.

Two GGMs are extracted using the formula \eqn{K = D(I - \Gamma)}, which is then symmetrized and forced to be positive definite:
\itemize{
  \item \bold{Within-cluster}: \eqn{D = \text{diag}(1/\sigma^2)} from residual variances; \eqn{\Gamma} from centered predictor fixed effects.
  \item \bold{Between-cluster}: \eqn{D = \text{diag}(1/\mu_{SD}^2)} from intercept random effect SDs; \eqn{\Gamma} from mean predictor fixed effects.
}

Use \code{\link{getNet}} with \code{type = "within"} or \code{type = "between"} to extract estimated network matrices.
}
\value{
An object of class \code{"mlGGM"} containing:
  \item{results}{A list with estimated network parameters, including \code{within}, \code{between}, \code{Gamma_within}, \code{Gamma_between}, and \code{mu}.}
  \item{output}{A list of fitted lmer model objects.}
  \item{fit}{A data frame with AIC and BIC per variable.}
  \item{data}{The augmented data used in estimation.}
  \item{model}{The predictor model specification.}
  \item{input}{A list of input arguments (vars, estimator, randomeffects, idvar).}
}
\author{
Sacha Epskamp
}
\examples{
\dontrun{
library(bootnet)

# Simulate data:
net_within <- genGGM(6, propPositive = 0.8)
net_between <- genGGM(6, p = 1, propPositive = 0.5)
within_data <- ggmGenerator()(5000, net_within)
within_data <- as.data.frame(within_data)
vars <- names(within_data)
within_data$cluster <- rep(1:50, each = 100)
between_data <- ggmGenerator()(50, net_between)
between_data <- as.data.frame(between_data)
between_data$cluster <- 1:50
data <- dplyr::left_join(within_data, between_data, by = "cluster",
                         suffix = c("_within", "_between"))
newdata <- data[, grep("_within", names(data))] +
           data[, grep("_between", names(data))]
names(newdata) <- vars
newdata$cluster <- data$cluster

# Fit mlGGM:
fit <- mlGGM(newdata, vars = vars, idvar = "cluster")

# Print and summary:
print(fit)
summary(fit)

# Plot networks:
plot(fit, "within")
plot(fit, "between")

# Extract network matrices:
getNet(fit, "within")
getNet(fit, "between")
}
}
\seealso{
\code{\link{mlVAR}}, \code{\link{getNet}}
}
