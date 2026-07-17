# Audit item 4: in R/lmer_murmur.R the "Computing random effects" loop
# (orthogonal/correlated contemporaneous branch) recomputed the per-subject
# posthoc residual covariance list Theta_posthoc on EVERY iteration of
# `for (p in 1:nRandom)` even though only element [[p]] was used: O(N^2)
# covariance computations for N subjects. The lapply() is now hoisted out of
# the loop (computed once, O(N)) and the list is named by subject ID so the
# p <-> subject correspondence is explicit (indexed by the ranef rownames)
# rather than relying on positional alignment. Results must be identical.
#
# This test checks CORRECTNESS of the named indexing: fitting on data whose
# rows arrive in non-sorted subject-ID order must give the same per-subject
# contemporaneous networks as fitting on pre-sorted data.

library(mlVAR)

set.seed(4)
sim <- mlVARsim(nPerson = 8, nNode = 3, nTime = 50, lag = 1)
dat <- sim$Data

# Relabel subjects so ID values are not in sorted order relative to the
# original generation order, then permute the per-subject blocks (keeping
# within-subject row order intact, as row order defines the time series) so
# subject IDs appear unsorted in the raw data:
idmap <- c(5, 3, 8, 1, 7, 2, 6, 4)
dat$ID <- idmap[match(dat$ID, sort(unique(dat$ID)))]

block_order <- c(8, 2, 5, 1, 7, 3, 6, 4)  # order in which ID blocks appear
dat_shuffled <- do.call(rbind, lapply(block_order, function(id) {
  dat[dat$ID == id, , drop = FALSE]
}))
rownames(dat_shuffled) <- NULL
dat_sorted <- do.call(rbind, lapply(sort(unique(dat$ID)), function(id) {
  dat[dat$ID == id, , drop = FALSE]
}))
rownames(dat_sorted) <- NULL

# Subjects must genuinely appear unsorted in the shuffled data:
stopifnot(is.unsorted(dat_shuffled$ID))

fit_shuf <- suppressWarnings(suppressMessages(
  mlVAR(dat_shuffled, vars = sim$vars, idvar = "ID", lags = 1,
        temporal = "orthogonal", contemporaneous = "orthogonal",
        verbose = FALSE)))
fit_sort <- suppressWarnings(suppressMessages(
  mlVAR(dat_sorted, vars = sim$vars, idvar = "ID", lags = 1,
        temporal = "orthogonal", contemporaneous = "orthogonal",
        verbose = FALSE)))

stopifnot(inherits(fit_shuf, "mlVAR"), inherits(fit_sort, "mlVAR"))
stopifnot(identical(sort(as.numeric(fit_shuf$IDs)), sort(unique(dat$ID))))

# --- 1. Per-subject contemporaneous networks: symmetric and finite ----------
for (p in seq_along(fit_shuf$IDs)) {
  for (mat in c("pcor", "cor", "cov", "prec")) {
    M <- fit_shuf$results$Theta[[mat]]$subject[[p]]
    stopifnot(is.matrix(M), all(is.finite(M)),
              isTRUE(all.equal(M, t(M), tolerance = 1e-8,
                               check.attributes = FALSE)))
  }
}
cat("OK: per-subject Theta matrices are finite and symmetric\n")

# --- 2. Same subject ID => same network, shuffled vs pre-sorted input -------
# Note: exact (1e-8) equality across input row orders is not attainable: the
# within-person centering (ave/mean) runs before mlVAR sorts the data, so
# floating-point summation order differs (~1e-15 in the augmented data),
# which lmer amplifies to ~1e-5 in per-subject estimates. A p <-> subject
# mix-up, by contrast, gives differences >= ~0.077 here. So we require
# same-ID agreement within 1e-3 AND that the same ID is by far the best
# match among all subjects.
for (id in unique(dat$ID)) {
  p_shuf <- match(as.character(id), as.character(fit_shuf$IDs))
  p_sort <- match(as.character(id), as.character(fit_sort$IDs))
  stopifnot(!is.na(p_shuf), !is.na(p_sort))
  for (mat in c("pcor", "cov")) {
    M_shuf <- fit_shuf$results$Theta[[mat]]$subject[[p_shuf]]
    diffs <- sapply(seq_along(fit_sort$IDs), function(q) {
      max(abs(M_shuf - fit_sort$results$Theta[[mat]]$subject[[q]]))
    })
    stopifnot(diffs[p_sort] < 1e-3)                     # same ID agrees
    stopifnot(all(diffs[-p_sort] > 10 * diffs[p_sort])) # and is the best match
  }
}
cat("OK: per-subject networks match by subject ID for shuffled vs sorted input\n")

# --- 3. getNet subject extraction works ---------------------------------------
pdf(NULL)
net1 <- getNet(fit_shuf, "contemporaneous", subject = 1, nonsig = "show",
               verbose = FALSE)
dev.off()
# getNet symmetrises and zeroes the diagonal; compare off-diagonals:
ref1 <- fit_shuf$results$Theta$pcor$subject[[1]]
ref1 <- (ref1 + t(ref1)) / 2
off <- row(ref1) != col(ref1)
stopifnot(is.matrix(net1), all(is.finite(net1)),
          isTRUE(all.equal(unname(net1)[off], unname(ref1)[off],
                           tolerance = 1e-8)))
cat("OK: getNet(fit, 'contemporaneous', subject = 1) matches Theta$pcor\n")

cat("All item 4 (Theta_posthoc hoist + named indexing) checks passed.\n")
