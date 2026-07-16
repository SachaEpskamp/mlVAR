# Audit item 2: temporal = "fixed" produced a Beta$SD array with a spurious
# extra lag slice. In R/lmer_murmur.R the "fixed" branch built the SD array as
#   array(0, c(nVar, nVar, length(unique(model$lag))))
# but model$lag contains NA for between-subject predictor rows, so unique()
# yielded one lag level too many (e.g. lag-1 fit gave dim(Beta$mean) = 3 3 1 but
# dim(Beta$SD) = 3 3 2). Fixed to use the pre-computed `lags` (NA-stripped).
#
# summary.mlVAR silently recycled c(Beta$SD) against nVar^2 * nLag rows; a
# defensive length check now surfaces any future mismatch instead of recycling.

library(mlVAR)

set.seed(1)
sim <- mlVARsim(nPerson = 6, nNode = 3, nTime = 80, lag = 1)

fit <- suppressWarnings(suppressMessages(
  mlVAR(sim$Data, vars = sim$vars, idvar = "ID", lags = 1,
        temporal = "fixed", verbose = FALSE)))
stopifnot(inherits(fit, "mlVAR"))

# --- 1. Beta$SD must have the same shape as Beta$mean -----------------------
stopifnot(identical(dim(fit$results$Beta$SD), dim(fit$results$Beta$mean)))
stopifnot(identical(dim(fit$results$Beta$SD), c(3L, 3L, 1L)))
cat("OK: dim(Beta$SD) == dim(Beta$mean) for temporal = 'fixed'\n")

# --- 2. summary() completes and its temporal DF has nVar^2 * nLag rows -------
out <- capture.output(s <- summary(fit))
nVar <- length(sim$vars)
nLag <- 1L
stopifnot(nrow(s$temporal) == nVar^2 * nLag)
cat("OK: summary() completes; temporal data frame has", nVar^2 * nLag, "rows\n")

# --- 3. plot(..., SD = TRUE) consumes Beta$SD without error -----------------
pdf(NULL)
p <- plot(fit, "temporal", SD = TRUE, nonsig = "show", verbose = FALSE, DoNotPlot = TRUE)
dev.off()
cat("OK: plot(fit, 'temporal', SD = TRUE) runs without error\n")

cat("All item 2 (temporal = 'fixed' Beta$SD) checks passed.\n")
