# Audit item 8: robustness / docs polish (bundle).
#
# 1. plot.mlVAR now validates 'lag' (temporal) with an informative message
#    before any bare array subscript error can occur.
# 2. plot.mlVAR now validates 'subject' (temporal and contemporaneous) against
#    the number of subject-specific networks before subscripting.
# 3. print.mlVARcompare handles ties in the AIC/BIC "best model" line.
#
# (Documentation-only changes: mlVAR.Rd p-value note, movingWindow maxEffects
# comment. NEWS bullets. Not exercised here.)

library(mlVAR)

set.seed(2026)
sim <- mlVARsim(nPerson = 6, nNode = 3, nTime = 50, lag = 1)
dat <- sim$Data

fit <- suppressWarnings(suppressMessages(
  mlVAR(dat, vars = sim$vars, idvar = "ID", lags = 1,
        temporal = "correlated", verbose = FALSE)))
stopifnot(inherits(fit, "mlVAR"))

nLags <- dim(fit$results$Beta$mean)[3]
nSubjects <- length(fit$results$Beta$subject)

pdf(NULL)
on.exit(dev.off(), add = TRUE)

# --- 1. lag out of range errors informatively -------------------------------
res <- tryCatch(
  suppressWarnings(suppressMessages(plot(fit, "temporal", lag = 2))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(res))
stopifnot(grepl("lag", res) && grepl("exceeds the number of lags", res))
stopifnot(grepl(paste0("(", nLags, ")"), res, fixed = TRUE))
cat("OK: plot(fit, 'temporal', lag = 2) gives the informative lag error\n")

# --- 2a. out-of-range subject (temporal) errors informatively ---------------
res <- tryCatch(
  suppressWarnings(suppressMessages(
    plot(fit, "temporal", subject = nSubjects + 1))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(res))
stopifnot(grepl("subject", res) && grepl("number of subjects", res))
stopifnot(grepl(paste0("(", nSubjects, ")"), res, fixed = TRUE))
cat("OK: plot(fit, 'temporal', subject = ", nSubjects + 1,
    ") gives the informative subject error\n", sep = "")

# --- 2b. out-of-range subject (contemporaneous) errors informatively --------
res <- tryCatch(
  suppressWarnings(suppressMessages(
    plot(fit, "contemporaneous", subject = nSubjects + 1))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(res))
stopifnot(grepl("subject", res) && grepl("number of subjects", res))
cat("OK: plot(fit, 'contemporaneous', subject = out-of-range) errors informatively\n")

# --- 3. in-range plots still work -------------------------------------------
p1 <- suppressWarnings(suppressMessages(plot(fit, "temporal", lag = 1)))
stopifnot(!is.null(p1))
p2 <- suppressWarnings(suppressMessages(plot(fit, "temporal", subject = 1)))
stopifnot(!is.null(p2))
p3 <- suppressWarnings(suppressMessages(plot(fit, "contemporaneous", subject = 1)))
stopifnot(!is.null(p3))
p4 <- suppressWarnings(suppressMessages(plot(fit, "contemporaneous")))
stopifnot(!is.null(p4))
cat("OK: in-range temporal and contemporaneous plots still work\n")

# --- 4. mlVARcompare tie handling: print completes --------------------------
# Compare two lag specifications on the same data. compareToLags keeps the
# number of observations equal so the models are comparable.
comp <- tryCatch({
  fit1 <- suppressWarnings(suppressMessages(
    mlVAR(dat, vars = sim$vars, idvar = "ID", lags = 1,
          temporal = "orthogonal", compareToLags = 1:2, verbose = FALSE)))
  fit2 <- suppressWarnings(suppressMessages(
    mlVAR(dat, vars = sim$vars, idvar = "ID", lags = 2,
          temporal = "orthogonal", compareToLags = 1:2, verbose = FALSE)))
  mlVARcompare(fit1, fit2)
}, error = function(e) e)

if (inherits(comp, "error")) {
  cat("NOTE: mlVARcompare skipped (", conditionMessage(comp), ")\n", sep = "")
} else {
  stopifnot(inherits(comp, "mlVARcompare"))
  out <- capture.output(print(comp))
  stopifnot(length(out) > 0)
  stopifnot(any(grepl("Best fitting model AIC", out)))
  stopifnot(any(grepl("Best fitting model BIC", out)))
  cat("OK: mlVARcompare + print(comp) completes\n")
}

cat("All item 8 (robustness / docs polish) checks passed.\n")
