\name{mlGGM}
\alias{mlGGM}
\alias{print.mlGGM}
\alias{summary.mlGGM}
\alias{plot.mlGGM}
\title{
Multi-Level Gaussian Graphical Model Estimation
}
\description{
Estimates within-cluster and between-cluster partial correlation networks from cross-sectional multilevel data (e.g., individuals nested in classrooms, patients nested in clinics) using a single-step nodewise regression approach. Unlike \code{\link{mlVAR}}, which uses a two-step procedure (temporal model first, then contemporaneous), \code{mlGGM} fits all effects in a single multilevel model per variable. This makes it suitable for cross-sectional data where no temporal ordering exists.

For each variable \eqn{y_i}, a single mixed-effects model is fitted with within-cluster centered predictors and cluster-mean predictors. Two undirected networks (GGMs) are extracted: one capturing partial correlations at the within-cluster level (individual deviations from cluster means) and one at the between-cluster level (cluster means).

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
A data frame containing the observed variables and a cluster identifier variable.
}
  \item{vars}{
A character vector of variable names to include in the GGM estimation. Must contain at least 2 variables.
}
  \item{idvar}{
A string indicating the name of the cluster/group identifier variable in \code{data}. Defaults to \code{"id"}.
}
  \item{estimator}{
The estimator to be used. Currently only \code{"lmer"} (sequential univariate multilevel estimation via \code{\link[lme4]{lmer}}) is supported.
}
  \item{randomeffects}{
How should random effects be estimated? \code{"correlated"} estimates correlated random slopes for within-cluster predictors (using \code{(1 + predictors | cluster)}), \code{"orthogonal"} estimates uncorrelated random slopes (using \code{(1 + predictors || cluster)}), and \code{"fixed"} estimates only a random intercept with no random slopes. \code{"default"} selects \code{"correlated"} when the number of variables is 6 or fewer and \code{"orthogonal"} otherwise, as correlated random effects become computationally expensive with many variables.
}
  \item{scale}{
Logical. Should variables be grand-mean standardized before estimation? Defaults to \code{TRUE}. Standardization ensures that edge weights are comparable across variables.
}
  \item{na.rm}{
Logical. Should rows with missing values be removed? Defaults to \code{TRUE}.
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
Character vector indicating which sections to show in the summary. Options: \code{"fit"} for information criteria (AIC/BIC), \code{"within"} for within-cluster partial correlations and p-values, \code{"between"} for between-cluster partial correlations and p-values.
}
  \item{round}{
Number of decimal places for rounding in the summary output. Defaults to 3.
}
  \item{type}{
The type of network to plot or extract: \code{"within"} for the within-cluster partial correlation network, \code{"between"} for the between-cluster partial correlation network.
}
  \item{partial}{
Logical. If \code{TRUE} (default), plots partial correlations. If \code{FALSE}, plots zero-order correlations.
}
  \item{SD}{
Logical. If \code{TRUE}, plots random effect standard deviations instead of fixed effects. Only available for the within-cluster network.
}
  \item{subject}{
Integer specifying which cluster's network to plot. Only available for the within-cluster network. When specified, the cluster-specific network (fixed effects + random effects for that cluster) is plotted.
}
  \item{order}{
Character vector or numeric vector specifying the order of nodes in the plot.
}
  \item{nonsig}{
How to handle non-significant edges in the plot: \code{"show"} shows all edges regardless of significance, \code{"hide"} removes non-significant edges, \code{"dashed"} shows non-significant edges as dashed lines. \code{"default"} hides non-significant edges for partial correlation networks and shows them otherwise.
}
  \item{rule}{
Significance rule for the undirected network: \code{"or"} requires at least one regression direction to be significant, \code{"and"} requires both directions to be significant. The \code{"and"} rule is more conservative and tends to improve specificity, especially for the between-cluster network.
}
  \item{alpha}{
Significance level for thresholding edges. Defaults to 0.05.
}
  \item{layout}{
Layout algorithm for the network plot, passed to \code{\link[qgraph]{qgraph}}. Defaults to \code{"spring"}.
}
  \item{\dots}{
Additional arguments passed to \code{\link[qgraph]{qgraph}} (for plot) or ignored.
}
}
\details{
\code{mlGGM} estimates multi-level Gaussian graphical models by fitting one mixed-effects regression model per variable. For each variable \eqn{y_i}, the model is:

\deqn{y_i = \beta_0 + \sum_{j \neq i} \gamma^{(w)}_{ij} (y_j - \bar{y}_j) + \sum_{j \neq i} \gamma^{(b)}_{ij} \bar{y}_j + u_{0} + \sum_{j \neq i} u_j (y_j - \bar{y}_j) + \varepsilon}

where \eqn{(y_j - \bar{y}_j)} are within-cluster centered predictors and \eqn{\bar{y}_j} are cluster-mean predictors. The random effects \eqn{u_0, u_j} are placed on the intercept and within-cluster predictors only (not on between-cluster means).

Two GGMs are extracted using the relationship between nodewise regression coefficients and the precision matrix:

\bold{Within-cluster network:} The precision matrix is computed as \eqn{K_w = D_w (I - \Gamma_w)}, where \eqn{D_w = \mathrm{diag}(1/\sigma^2_i)} contains the inverse residual variances and \eqn{\Gamma_w} contains the fixed effects of the within-cluster centered predictors. The matrix is symmetrized and forced to be positive definite. Partial correlations are derived from the precision matrix.

\bold{Between-cluster network:} The precision matrix is computed as \eqn{K_b = D_b (I - \Gamma_b)}, where \eqn{D_b = \mathrm{diag}(1/\mu^2_{SD,i})} contains the inverse squared random intercept standard deviations and \eqn{\Gamma_b} contains the fixed effects of the cluster-mean predictors.

\bold{Cluster-specific networks:} When \code{randomeffects} is \code{"correlated"} or \code{"orthogonal"}, cluster-specific within-cluster networks are also computed by adding the estimated random effects to the fixed effects for each cluster. These can be accessed via the \code{subject} argument in the plot method.

Network matrices can be extracted using \code{\link{getNet}} with \code{type = "within"} or \code{type = "between"}. The \code{rule} argument controls the significance thresholding rule.
}
\value{
An object of class \code{"mlGGM"} containing:
  \item{results}{A list with:
    \itemize{
      \item \code{within}: Within-cluster network estimates (covariance, correlation, partial correlation matrices, each with \code{$mean} and \code{$subject} components)
      \item \code{between}: Between-cluster network estimates (covariance, correlation, partial correlation matrices, with \code{$mean} component)
      \item \code{Gamma_within}: Within-cluster regression coefficients with \code{$mean}, \code{$SE}, \code{$P}, and \code{$subject} components
      \item \code{Gamma_between}: Between-cluster regression coefficients with \code{$mean}, \code{$SE}, and \code{$P} components
      \item \code{mu}: Intercept estimates with \code{$mean}, \code{$SE}, and \code{$SD} components
    }
  }
  \item{output}{A list of fitted \code{\link[lme4]{lmer}} model objects, one per variable.}
  \item{fit}{A data frame with AIC and BIC for each variable's model.}
  \item{data}{The augmented data frame used in estimation.}
  \item{model}{The predictor model specification data frame.}
  \item{input}{A list of input arguments: \code{vars}, \code{estimator}, \code{randomeffects}, \code{idvar}.}
}
\references{
Epskamp, S., Waldorp, L. J., Mottus, R., & Borsboom, D. (2018). The Gaussian graphical model in cross-sectional and time-series data. \emph{Multivariate Behavioral Research}, 53(4), 453-480.

Epskamp, S., Rhemtulla, M. T., & Borsboom, D. (2017). Generalized network psychometrics: Combining network and latent variable models. \emph{Psychometrika}, 82(4), 904-927.
}
\author{
Sacha Epskamp
}
\examples{
\donttest{
# Generate simulated multilevel data with known GGM structure:
library(bootnet)

nNode <- 4
nCluster <- 30
clusterSize <- 50

# True within-cluster and between-cluster networks:
set.seed(1)
net_within <- genGGM(nNode, propPositive = 0.8)
net_between <- genGGM(nNode, p = 1, propPositive = 0.5)

# Generate within-cluster deviations:
within_data <- ggmGenerator()(nCluster * clusterSize, net_within)
within_data <- as.data.frame(within_data)
vars <- names(within_data)
within_data$cluster <- rep(1:nCluster, each = clusterSize)

# Generate between-cluster means:
between_data <- ggmGenerator()(nCluster, net_between)
between_data <- as.data.frame(between_data)
between_data$cluster <- 1:nCluster

# Combine: observed = within deviation + cluster mean
merged <- dplyr::left_join(within_data, between_data, by = "cluster",
                           suffix = c("_w", "_b"))
data <- merged[, grep("_w$", names(merged))] +
        merged[, grep("_b$", names(merged))]
names(data) <- vars
data$cluster <- merged$cluster

# Fit the multi-level GGM:
fit <- mlGGM(data, vars = vars, idvar = "cluster")

# Inspect results:
summary(fit)

# Plot within-cluster and between-cluster networks:
plot(fit, type = "within")
plot(fit, type = "between")

# Extract network matrices (with AND rule for better specificity):
within_net <- getNet(fit, "within", rule = "and")
between_net <- getNet(fit, "between", rule = "and")
}
}
\seealso{
\code{\link{mlVAR}} for multilevel VAR estimation with temporal data, \code{\link{getNet}} for extracting network matrices.
}
