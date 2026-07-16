# CRAN comments for mlVAR 0.7

This is a version update of the existing CRAN package mlVAR. The current CRAN
release is 0.6.1. Version 0.7 is a bug-fix release; Adela-Maria Isvoranu has
been added as a contributor in Authors@R.

## Release summary

See NEWS for the full list. Highlights:

* Fixed the rank-deficiency check silently reordering the retained variables:
  when linearly dependent variables were dropped, the survivors were kept in
  QR-pivot order rather than input order, changing the node order in all
  outputs.
* Fixed the parametric bootstrap `exclude` argument in the vendored `parSim()`
  being inverted: matching rows were kept rather than dropped, so `exclude`
  had the opposite of its documented effect.
* Fixed `temporal = "fixed"` producing a `Beta$SD` array with a spurious extra
  lag slice (sized from lags including an NA level).
* Fixed `mlVARsim()` storing the `mu_SD` argument range in `mu$SD` rather than
  the realised per-node standard deviations.
* Fixed a parallel-cluster leak: with `nCores > 1` the PSOCK cluster was only
  stopped on one branch, leaving worker processes open on other paths and on
  error. It is now stopped unconditionally via `on.exit()`.
* `nCores` now uses the requested number of worker processes (previously one
  fewer worker was spawned).
* Fixed `estimator = "lm"` crashing for non-numeric subject IDs.
* Performance: the random-effects stage of two-step contemporaneous estimation
  no longer recomputes all subjects' posthoc residual covariances on every
  iteration (O(N^2) -> O(N)). Results are identical.
* Calling `mlVAR()` without `idvar`, or with `estimator = "JAGS"`, now gives an
  informative error.
* Added a test suite under `tests/`.

## Note on documentation restored in this release

CRAN's 0.6.1 was built from a working copy containing two documentation edits
that were never committed to the package's git repository, so they were absent
from the source this release was prepared from. Both have been restored here:
the `\item{\dots}` entry in `mlVAR.Rd` (required, since `\dots` appears in the
`\usage` section) and the `resimulate(variance = "empirical")` example.

The `mlGGM` example is wrapped in `\donttest{}` rather than the `\dontrun{}`
that shipped in 0.6.1. It has been verified to run under
`R CMD check --as-cran --run-donttest` in well under a second.

## Test environments

* Local: macOS Sonoma 14.2.1 (aarch64-apple-darwin20), R 4.5.3 (2026-03-11)

## R CMD check --as-cran results

0 ERRORs, 0 WARNINGs, 1 NOTE.

The NOTE is a property of the local test machine rather than of the package:

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent enough
HTML Tidy.
```

macOS ships an outdated HTML Tidy, so this check is skipped locally; it is
expected to pass on the CRAN build machines.
