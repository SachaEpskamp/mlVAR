# Backward-compatibility reference values captured from mlVAR 0.6.1
# (commit e43dab5) on R 4.6.1 — must not change unless a NEWS-documented
# intentional break is made.
#
# All reference literals below were generated with dput(round(..., 8)) from
# the installed, unmodified mlVAR 0.6.1 and verified byte-identical across
# two independent Rscript runs. Deterministic data-generation checks use a
# 1e-6 tolerance; checks involving iterative model fits (lmer/lm) use 1e-3,
# since optimizer termination varies across platforms/BLAS builds by up to
# ~1e-4 in relative terms (observed on win-builder), while a real behavioral
# break would move results by far more.

library(mlVAR)

check <- function(name, actual, expected, tol = 1e-3) {
  actual <- round(unname(actual), 8)
  res <- all.equal(actual, expected, tolerance = tol, check.attributes = FALSE)
  if (!isTRUE(res)) {
    cat("MISMATCH in", name, ":\n")
    print(res)
    cat("Actual:\n"); print(actual)
    cat("Expected:\n"); print(expected)
    stop("Backward-compatibility check failed: ", name)
  }
  cat("OK:", name, "\n")
}

# ---------------------------------------------------------------------------
# Simulated data (generated once, reused by every block)
# ---------------------------------------------------------------------------
set.seed(1)
sim <- mlVARsim(nPerson = 8, nNode = 3, nTime = 60)

# e. mlVARsim model truth + generated-data digest.
#    NOTE: sim$model$mu$SD is deliberately NOT pinned — v0.6.1 stores the
#    wrong quantity there (the mu_SD argument range instead of the realised
#    SDs), and a documented bug fix will change its contents. The Beta array
#    and the generated data itself must stay identical regardless.
check("mlVARsim Beta mean (model truth)", sim$model$Beta$mean,
      structure(c(0.25982487, 0, 0, 0, 0.20877215, -0.34575695,
                  -0.47943937, 0.68017766, 0.04683755), dim = c(3L, 3L, 1L)),
      tol = 1e-6)

stopifnot(identical(dim(sim$Data), c(480L, 4L)))
cat("OK: mlVARsim Data dimensions\n")

check("mlVARsim Data digest (sum of abs)",
      sum(abs(as.matrix(sim$Data[, sim$vars]))), 2265.68404, tol = 1e-6)
check("mlVARsim Data column sums",
      colSums(as.matrix(sim$Data[, sim$vars])),
      c(655.59966922, -453.32070334, 242.8679717), tol = 1e-6)

# ---------------------------------------------------------------------------
# a. Default lmer estimation (flagship path): temporal = contemporaneous =
#    "correlated" for 3 nodes.
# ---------------------------------------------------------------------------
set.seed(42)
fit <- suppressWarnings(suppressMessages(
  mlVAR(sim$Data, vars = sim$vars, idvar = sim$idvar, lags = 1,
        verbose = FALSE)))
stopifnot(inherits(fit, "mlVAR"))

check("default lmer: Beta fixed effects (mean)", fit$results$Beta$mean,
      structure(c(0.19834461, 0.05709925, -0.16587381, -0.18436624,
                  0.18707047, -0.55290386, -0.39984319, 0.71747531,
                  0.02210775), dim = c(3L, 3L, 1L)))

check("default lmer: Beta random-effect SDs", fit$results$Beta$SD,
      structure(c(0.36950579, 0.20379508, 0.13754687, 0.67260263,
                  0.11762314, 0.24060748, 0.12654103, 0.22480877,
                  0.30915927), dim = c(3L, 3L, 1L)))

check("default lmer: temporal network (getNet)",
      getNet(fit, "temporal", nonsig = "show"),
      structure(c(0.19834461, -0.18436624, -0.39984319, 0.05709925,
                  0.18707047, 0.71747531, -0.16587381, -0.55290386,
                  0.02210775), dim = c(3L, 3L)))

check("default lmer: contemporaneous network (getNet)",
      getNet(fit, "contemporaneous", nonsig = "show"),
      structure(c(0, -0.56710865, 0.18022122, -0.56710865, 0,
                  -0.62252982, 0.18022122, -0.62252982, 0),
                dim = c(3L, 3L)))

check("default lmer: between-subjects network (getNet)",
      getNet(fit, "between", nonsig = "show"),
      structure(c(0, 0.81362429, 0.88286285, 0.81362429, 0,
                  -0.515902, 0.88286285, -0.515902, 0),
                dim = c(3L, 3L)))

# ---------------------------------------------------------------------------
# b. Orthogonal estimator path
# ---------------------------------------------------------------------------
set.seed(42)
fit2 <- suppressWarnings(suppressMessages(
  mlVAR(sim$Data, vars = sim$vars, idvar = sim$idvar, lags = 1,
        temporal = "orthogonal", contemporaneous = "orthogonal",
        verbose = FALSE)))

check("orthogonal: temporal network (getNet)",
      getNet(fit2, "temporal", nonsig = "show"),
      structure(c(0.19571391, -0.24888247, -0.44896941, 0.04693608,
                  0.19800571, 0.72517566, -0.20094253, -0.61454947,
                  0.01332294), dim = c(3L, 3L)))

check("orthogonal: contemporaneous network (getNet)",
      getNet(fit2, "contemporaneous", nonsig = "show"),
      structure(c(0, -0.55585573, 0.15771588, -0.55585573, 0,
                  -0.66979703, 0.15771588, -0.66979703, 0),
                dim = c(3L, 3L)))

# ---------------------------------------------------------------------------
# c. lm estimator with numeric IDs (temporal/contemporaneous = "unique").
#    sim$Data$ID is integer; results for numeric IDs must stay identical
#    after any character-ID bug fix.
# ---------------------------------------------------------------------------
stopifnot(is.numeric(sim$Data$ID))
set.seed(42)
fit3 <- suppressWarnings(suppressMessages(
  mlVAR(sim$Data, vars = sim$vars, idvar = sim$idvar, lags = 1,
        estimator = "lm", temporal = "unique", contemporaneous = "unique",
        verbose = FALSE)))

check("lm/unique: mean temporal network (getNet)",
      getNet(fit3, "temporal", nonsig = "show"),
      structure(c(0.22783736, -0.15432406, -0.41340804, 0.01594413,
                  0.14881262, 0.7160717, -0.14986572, -0.47754751,
                  0.04558778), dim = c(3L, 3L)))

check("lm/unique: subject 1 Beta", fit3$results$Beta$subject[[1]],
      structure(c(0.19079773, -0.19164951, 0.18782874, 0.33190971,
                  -0.11709151, -0.25586798, -0.57874105, 0.66043975,
                  0.24172542), dim = c(3L, 3L, 1L)))

# ---------------------------------------------------------------------------
# d. summary() numeric columns of the temporal data frame
#    (summary.mlVAR prints as a side effect; capture and discard)
# ---------------------------------------------------------------------------
invisible(capture.output(s <- summary(fit)))
stopifnot(is.data.frame(s$temporal), nrow(s$temporal) == 9L)

check("summary: temporal fixed effects", round(s$temporal$fixed, 6),
      c(0.198, 0.057, -0.166, -0.184, 0.187, -0.553, -0.4, 0.717, 0.022))
check("summary: temporal SEs", round(s$temporal$SE, 6),
      c(0.147, 0.091, 0.09, 0.248, 0.06, 0.103, 0.074, 0.093, 0.126))
check("summary: temporal random-effect SDs", round(s$temporal$ran_SD, 6),
      c(0.37, 0.204, 0.138, 0.673, 0.118, 0.241, 0.127, 0.225, 0.309))

cat("All backward-compatibility checks passed.\n")
