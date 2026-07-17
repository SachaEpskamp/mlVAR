# Audit item 1: estimator = "lm" crashed for non-numeric (character) subject
# IDs because per-subject residuals were collected via cbind() with the id
# column, coercing the whole residual matrix to character. This broke cov()
# (contemporaneous = "unique") and the residual lmer() step (correlated/
# orthogonal). Fixed by building the residuals as a type-preserving
# data.frame in R/lm_mlVAR.R.
#
# Numeric-ID results being unchanged vs v0.6.1 is pinned by
# tests/test-backcompat.R (block c).

library(mlVAR)

set.seed(123)
sim <- mlVARsim(nPerson = 4, nNode = 3, nTime = 60)

dat_num <- sim$Data
stopifnot(is.numeric(dat_num$ID))

# Character-ID copy of the exact same data:
dat_char <- dat_num
dat_char$ID <- paste0("p", dat_char$ID)
stopifnot(is.character(dat_char$ID))

# --- 1. lm estimator with character IDs, contemporaneous = "unique" --------
fit_char_uni <- suppressWarnings(suppressMessages(
  mlVAR(dat_char, vars = sim$vars, idvar = "ID", lags = 1,
        estimator = "lm", temporal = "unique", contemporaneous = "unique",
        verbose = FALSE)))
stopifnot(inherits(fit_char_uni, "mlVAR"))
cat("OK: lm estimator with character IDs (contemporaneous = 'unique') completes\n")

# --- 2. lm estimator with character IDs, contemporaneous = "correlated" ----
fit_char_cor <- suppressWarnings(suppressMessages(
  mlVAR(dat_char, vars = sim$vars, idvar = "ID", lags = 1,
        estimator = "lm", temporal = "unique",
        contemporaneous = "correlated", verbose = FALSE)))
stopifnot(inherits(fit_char_cor, "mlVAR"))
cat("OK: lm estimator with character IDs (contemporaneous = 'correlated') completes\n")

# --- 3. Equivalence: character-ID results == numeric-ID results ------------
fit_num_uni <- suppressWarnings(suppressMessages(
  mlVAR(dat_num, vars = sim$vars, idvar = "ID", lags = 1,
        estimator = "lm", temporal = "unique", contemporaneous = "unique",
        verbose = FALSE)))
fit_num_cor <- suppressWarnings(suppressMessages(
  mlVAR(dat_num, vars = sim$vars, idvar = "ID", lags = 1,
        estimator = "lm", temporal = "unique",
        contemporaneous = "correlated", verbose = FALSE)))

stopifnot(isTRUE(all.equal(
  getNet(fit_char_uni, "temporal", nonsig = "show"),
  getNet(fit_num_uni, "temporal", nonsig = "show"),
  tolerance = 1e-8, check.attributes = FALSE)))
stopifnot(isTRUE(all.equal(
  getNet(fit_char_uni, "contemporaneous", nonsig = "show"),
  getNet(fit_num_uni, "contemporaneous", nonsig = "show"),
  tolerance = 1e-8, check.attributes = FALSE)))
cat("OK: character-ID == numeric-ID networks (contemporaneous = 'unique')\n")

stopifnot(isTRUE(all.equal(
  getNet(fit_char_cor, "temporal", nonsig = "show"),
  getNet(fit_num_cor, "temporal", nonsig = "show"),
  tolerance = 1e-8, check.attributes = FALSE)))
stopifnot(isTRUE(all.equal(
  getNet(fit_char_cor, "contemporaneous", nonsig = "show"),
  getNet(fit_num_cor, "contemporaneous", nonsig = "show"),
  tolerance = 1e-8, check.attributes = FALSE)))
cat("OK: character-ID == numeric-ID networks (contemporaneous = 'correlated')\n")

cat("All item 1 (lm estimator with character IDs) checks passed.\n")
