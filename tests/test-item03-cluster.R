# Audit item 3: parallel-cluster lifecycle in R/lmer_murmur.R (lmer_mlVAR).
# With nCores > 1 a PSOCK cluster was created but stopCluster() only ran inside
# the two-step contemporaneous branch, so contemporaneous = "fixed"/"unique"
# runs (and any error path) leaked worker processes and their connections.
# The cluster is now stopped unconditionally via on.exit(). This test runs the
# previously-leaking path and asserts open connections return to baseline.
#
# Gated behind NOT_CRAN=true: spawning PSOCK clusters in tests is flaky on CRAN.

library(mlVAR)

if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
  cat("SKIP: test-item03-cluster.R (set NOT_CRAN=true to run cluster test)\n")
  quit(save = "no", status = 0)
}

set.seed(123)
sim <- mlVARsim(nPerson = 4, nNode = 3, nTime = 50, lag = 1)

# Helper: current number of open connections.
nConn <- function() nrow(showConnections())

# Helper: settle connections and count (workers close asynchronously).
settledConn <- function() {
  Sys.sleep(1)
  gc()
  nConn()
}

# --- 1. Previously-leaking path: contemporaneous = "fixed", nCores = 2 -------
baseline <- nConn()

fit1 <- suppressWarnings(suppressMessages(
  mlVAR(sim$Data, vars = sim$vars, idvar = "ID", lags = 1,
        nCores = 2, temporal = "orthogonal", contemporaneous = "fixed",
        verbose = FALSE)))
stopifnot(inherits(fit1, "mlVAR"))

after1 <- settledConn()
stopifnot(after1 <= baseline)
cat("OK: contemporaneous = 'fixed', nCores = 2 completes and leaks no connections",
    "(baseline =", baseline, ", after =", after1, ")\n")

# --- 2. Two-step path: contemporaneous = "correlated", nCores = 2 -----------
baseline2 <- nConn()

fit2 <- suppressWarnings(suppressMessages(
  mlVAR(sim$Data, vars = sim$vars, idvar = "ID", lags = 1,
        nCores = 2, temporal = "orthogonal", contemporaneous = "correlated",
        verbose = FALSE)))
stopifnot(inherits(fit2, "mlVAR"))

after2 <- settledConn()
stopifnot(after2 <= baseline2)
cat("OK: contemporaneous = 'correlated', nCores = 2 completes and leaks no connections",
    "(baseline =", baseline2, ", after =", after2, ")\n")

cat("All item-03 cluster-lifecycle checks passed\n")
