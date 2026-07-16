# Audit item 7: dead code and contradictory checks in R/mlVAR.R (+ one in
# R/lmer_murmur.R).
#
# 1. estimator = "JAGS" now gives its informative "not implemented" message
#    (previously unreachable: match.arg() rejected "JAGS" first with a
#    generic error).
# 2. Missing 'idvar' now gives an informative error (previously
#    stopifnot(!missing(idvar)) was immediately followed by a dead
#    auto-"ID" branch; enabling that branch was tested and fails inside
#    lmer, since random effects require >= 2 groups).
# 3. The rank-deficiency variable drop no longer reorders the retained
#    variables into QR-pivot order: keep <- sort(qrX$pivot[1:rnk]).
# 4. For contemporaneous = "fixed"/"unique" (lmer path), output$contemporaneous
#    is now NULL instead of being overloaded with a modelCov object
#    (nothing reads that field; second-step lmer fits only exist for
#    correlated/orthogonal).

library(mlVAR)

set.seed(2026)
sim <- mlVARsim(nPerson = 6, nNode = 3, nTime = 50)
dat <- sim$Data

# --- 1. JAGS estimator gives informative error ------------------------------
res <- tryCatch(
  mlVAR(dat, vars = sim$vars, idvar = "ID", lags = 1,
        estimator = "JAGS", verbose = FALSE),
  error = function(e) conditionMessage(e))
stopifnot(is.character(res))
stopifnot(grepl("JAGS", res) && grepl("not implemented", res))
cat("OK: estimator = 'JAGS' gives the informative 'not implemented' error\n")

# --- 2. Missing idvar gives informative error -------------------------------
res <- tryCatch(
  mlVAR(dat, vars = sim$vars, lags = 1, verbose = FALSE),
  error = function(e) conditionMessage(e))
stopifnot(is.character(res))
stopifnot(grepl("idvar", res) && grepl("graphicalVAR", res))
cat("OK: missing 'idvar' gives an informative error\n")

# --- 3. Rank-deficiency drop preserves input variable order -----------------
# Insert a linearly dependent variable in the MIDDLE of vars, so that the
# old QR-pivot ordering (which put pivoted columns first) would have
# scrambled the order of the retained variables:
dat2 <- dat
dat2$V1b <- dat2$V1 + dat2$V2  # exact linear dependence
vars_rd <- c("V1", "V2", "V1b", "V3")

warns <- character(0)
fit_rd <- withCallingHandlers(
  mlVAR(dat2, vars = vars_rd, idvar = "ID", lags = 1, verbose = FALSE),
  warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  })

# (a) The rank-deficiency branch fired:
stopifnot(any(grepl("linearly dependent", warns)))
# Exactly one variable was dropped, and it is the dependent one:
stopifnot(length(fit_rd$input$vars) == 3)
dropped <- setdiff(vars_rd, fit_rd$input$vars)
stopifnot(length(dropped) == 1)

# (b) Retained variables preserve the original input order:
stopifnot(identical(fit_rd$input$vars,
                    vars_rd[vars_rd %in% fit_rd$input$vars]))
# Node order in the estimated networks matches too:
net_rd <- getNet(fit_rd, "temporal", nonsig = "show")
stopifnot(identical(colnames(net_rd), fit_rd$input$vars))
cat("OK: rank-deficiency drop fires with a warning and preserves input order",
    "(dropped:", dropped, ")\n")

# --- 4. contemporaneous = "fixed" path works; output field not overloaded ---
fit_fix <- suppressWarnings(suppressMessages(
  mlVAR(dat, vars = sim$vars, idvar = "ID", lags = 1,
        contemporaneous = "fixed", verbose = FALSE)))
stopifnot(inherits(fit_fix, "mlVAR"))
# output$contemporaneous is NULL (no second-step lmer fits exist), not a
# type-inconsistent modelCov object:
stopifnot(is.null(fit_fix$output$contemporaneous))
# End-to-end: contemporaneous network is still available:
net_c <- getNet(fit_fix, "contemporaneous", nonsig = "show")
stopifnot(is.matrix(net_c), nrow(net_c) == 3, ncol(net_c) == 3)
# Temporal lmer fits are untouched:
stopifnot(length(fit_fix$output$temporal) == 3)
cat("OK: contemporaneous = 'fixed' fit completes; output$contemporaneous is NULL\n")

cat("All item 7 (dead code / contradictory checks) checks passed.\n")
