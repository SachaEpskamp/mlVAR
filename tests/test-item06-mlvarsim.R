# Item 06: mlVARsim now stores the realised per-node SDs in model$mu$SD
# (previously it stored the mu_SD argument range, e.g. c(1,1)).
library(mlVAR)

# ---------------------------------------------------------------------------
# 1. Basic structure: one realised SD per node, all strictly positive.
# ---------------------------------------------------------------------------
set.seed(2)
sim <- mlVARsim(nPerson = 20, nNode = 4, nTime = 50)
stopifnot(length(sim$model$mu$SD) == 4)
stopifnot(all(sim$model$mu$SD > 0))
cat("OK: model$mu$SD has one positive value per node\n")

# ---------------------------------------------------------------------------
# 2. The stored SD should track the empirical spread of the person-specific
#    means. model$mu$subject is a list (length nPerson) of length-nNode
#    vectors; extract node j across subjects with sapply(..., "[", j).
# ---------------------------------------------------------------------------
set.seed(42)
simBig <- mlVARsim(nPerson = 200, nNode = 4, nTime = 50)
for (j in 1:4) {
  empSD <- sd(sapply(simBig$model$mu$subject, "[", j))
  theSD <- simBig$model$mu$SD[j]
  stopifnot(empSD >= 0.3 * theSD, empSD <= 3 * theSD)
}
cat("OK: realised SDs match empirical spread of person means (0.3-3x)\n")

# ---------------------------------------------------------------------------
# 3. RNG-stream guard: the generated data must be byte-identical to the
#    pre-change installed package. The digest literal is reused verbatim from
#    tests/test-backcompat.R (do NOT modify that file).
# ---------------------------------------------------------------------------
set.seed(1)
s1 <- mlVARsim(8, 3, 60)
stopifnot(isTRUE(all.equal(
  sum(abs(as.matrix(s1$Data[, s1$vars]))), 2265.68404, tolerance = 1e-6)))
cat("OK: generated-data digest unchanged (RNG stream intact)\n")

cat("All item-06 checks passed.\n")
