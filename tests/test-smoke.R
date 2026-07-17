# Smoke test: fastest possible end-to-end run of the flagship mlVAR pipeline.
# Simulate tiny data -> fit with defaults -> summary() -> plot().

library(mlVAR)

set.seed(123)
sim <- mlVARsim(nPerson = 5, nNode = 3, nTime = 30)
stopifnot(inherits(sim, "mlVARsim"), is.data.frame(sim$Data))
cat("OK: mlVARsim\n")

fit <- suppressWarnings(suppressMessages(
  mlVAR(sim$Data, vars = sim$vars, idvar = sim$idvar, lags = 1,
        verbose = FALSE)))
stopifnot(inherits(fit, "mlVAR"))
cat("OK: mlVAR fit (class mlVAR)\n")

invisible(capture.output(s <- summary(fit)))
stopifnot(is.data.frame(s$temporal), nrow(s$temporal) > 0)
cat("OK: summary()\n")

# plot.mlVAR draws via qgraph, which needs a graphics device:
pdf(NULL)
res <- try(suppressMessages(invisible(capture.output(
  plot(fit, "temporal")))), silent = TRUE)
invisible(dev.off())
stopifnot(!inherits(res, "try-error"))
cat("OK: plot(fit, \"temporal\")\n")

cat("Smoke test passed.\n")
