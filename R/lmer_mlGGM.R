# lmer estimation backend for mlGGM

lmer_mlGGM <- function(model, augData, idvar, vars, randomeffects = "correlated", verbose = TRUE) {

  # Use vars to preserve original variable ordering
  # (unique(model$dep) can be wrong after removing self-predictions)
  Outcomes <- vars
  nVar <- length(Outcomes)

  # Output list:
  lmerResults <- list()

  if (verbose) {
    message("Estimating within-cluster and between-cluster effects")
    pb <- txtProgressBar(min = 0, max = length(Outcomes), style = 3)
  }

  # Identify within and between predictor IDs for each outcome:
  for (i in seq_along(Outcomes)) {
    subModel <- model %>% filter(.data[['dep']] == Outcomes[i])

    withinIDs <- subModel$predID[subModel$type == "within"]
    betweenIDs <- subModel$predID[subModel$type == "between"]
    allPredIDs <- c(withinIDs, betweenIDs)

    # Formula construction:
    if (randomeffects == "correlated") {
      mod <- paste0(
        Outcomes[i], " ~ ",
        paste(allPredIDs, collapse = " + "),
        " + (1 + ", paste(withinIDs, collapse = " + "),
        " | ", idvar, ")"
      )
    } else if (randomeffects == "orthogonal") {
      mod <- paste0(
        Outcomes[i], " ~ ",
        paste(allPredIDs, collapse = " + "),
        " + (1 + ", paste(withinIDs, collapse = " + "),
        " || ", idvar, ")"
      )
    } else if (randomeffects == "fixed") {
      mod <- paste0(
        Outcomes[i], " ~ ",
        paste(allPredIDs, collapse = " + "),
        " + (1 | ", idvar, ")"
      )
    }

    formula <- as.formula(mod)
    suppressMessages(suppressWarnings(
      lmerResults[[i]] <- lmer(formula, data = augData, REML = FALSE)
    ))

    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbose) {
    close(pb)
  }

  ### Collect the results ###
  Results <- list()

  # Identify predictor ID mappings:
  withinMod <- model %>%
    filter(.data[['type']] == "within") %>%
    group_by(.data[["pred"]]) %>%
    dplyr::summarize(id = unique(.data[['predID']]), .groups = "drop")
  withinIDs_ordered <- withinMod$id[match(Outcomes, withinMod$pred)]

  betweenMod <- model %>%
    filter(.data[['type']] == "between") %>%
    group_by(.data[["pred"]]) %>%
    dplyr::summarize(id = unique(.data[['predID']]), .groups = "drop")
  betweenIDs_ordered <- betweenMod$id[match(Outcomes, betweenMod$pred)]

  ### Intercepts (mu) ###
  mu_fixed <- unname(sapply(lmerResults, function(x) fixef(x)[1]))
  names(mu_fixed) <- Outcomes

  mu_SD <- sapply(lmerResults, function(x) attr(lme4::VarCorr(x)[[1]], "stddev")[1])
  names(mu_SD) <- Outcomes

  mu_SE <- unname(sapply(lmerResults, function(x) se.fixef(x)[1]))
  names(mu_SE) <- Outcomes

  mu_P <- 2 * (1 - pnorm(abs(mu_fixed / mu_SE)))

  # Cluster-specific intercepts:
  rans <- lapply(lmerResults, function(x) ranef(x)[[idvar]][, 1])
  mu_subject <- lapply(seq_along(rans[[1]]), function(p) {
    mu <- mu_fixed + sapply(rans, '[[', p)
    names(mu) <- Outcomes
    mu
  })

  Results[['mu']] <- modelArray(mean = mu_fixed, SE = mu_SE, SD = mu_SD, subject = mu_subject)

  ### WITHIN-CLUSTER Gamma ###
  Gamma_within <- Gamma_within_SE <- matrix(0, nVar, nVar)
  for (i in seq_along(Outcomes)) {
    Gamma_within[i, -i] <- fixef(lmerResults[[i]])[withinIDs_ordered[-i]]
    Gamma_within_SE[i, -i] <- se.fixef(lmerResults[[i]])[withinIDs_ordered[-i]]
  }
  rownames(Gamma_within) <- colnames(Gamma_within) <- Outcomes
  rownames(Gamma_within_SE) <- colnames(Gamma_within_SE) <- Outcomes

  Gamma_within_P <- 2 * (1 - pnorm(abs(Gamma_within / Gamma_within_SE)))
  diag(Gamma_within_P) <- 0

  Results[["Gamma_within"]] <- modelArray(
    mean = Gamma_within,
    SE = Gamma_within_SE,
    P = Gamma_within_P
  )

  ### Within-cluster precision and covariance ###
  sigma2 <- sapply(lmerResults, sigma)^2
  D_within <- diag(1 / sigma2)

  inv_within <- D_within %*% (diag(nVar) - Gamma_within)
  inv_within <- (inv_within + t(inv_within)) / 2
  inv_within <- forcePositive(inv_within)

  cov_within <- corpcor::pseudoinverse(inv_within)
  colnames(cov_within) <- rownames(cov_within) <- Outcomes
  colnames(inv_within) <- rownames(inv_within) <- Outcomes

  Results[["within"]] <- modelCov(
    cov = modelArray(mean = cov_within),
    prec = modelArray(mean = inv_within)
  )

  ### BETWEEN-CLUSTER Gamma ###
  Gamma_between <- Gamma_between_SE <- matrix(0, nVar, nVar)
  for (i in seq_along(Outcomes)) {
    Gamma_between[i, -i] <- fixef(lmerResults[[i]])[betweenIDs_ordered[-i]]
    Gamma_between_SE[i, -i] <- se.fixef(lmerResults[[i]])[betweenIDs_ordered[-i]]
  }
  rownames(Gamma_between) <- colnames(Gamma_between) <- Outcomes
  rownames(Gamma_between_SE) <- colnames(Gamma_between_SE) <- Outcomes

  Gamma_between_P <- 2 * (1 - pnorm(abs(Gamma_between / Gamma_between_SE)))
  diag(Gamma_between_P) <- 0

  Results[["Gamma_between"]] <- modelArray(
    mean = Gamma_between,
    SE = Gamma_between_SE,
    P = Gamma_between_P
  )

  ### Between-cluster precision and covariance ###
  if (any(mu_SD == 0)) {
    warning("Zero SD found in intercept of following variables: ",
            paste(names(mu_SD[which(mu_SD == 0)]), collapse = ", "),
            " - Between-cluster effects could not be estimated")

    cov_between <- inv_between <- matrix(NA, nVar, nVar)
    colnames(cov_between) <- rownames(cov_between) <- Outcomes
    colnames(inv_between) <- rownames(inv_between) <- Outcomes

    Results[["between"]] <- modelCov(
      cov = modelArray(mean = cov_between),
      prec = modelArray(mean = inv_between)
    )
  } else {
    D_between <- diag(1 / mu_SD^2)
    inv_between <- D_between %*% (diag(nVar) - Gamma_between)
    inv_between <- (inv_between + t(inv_between)) / 2
    inv_between <- forcePositive(inv_between)

    cov_between <- corpcor::pseudoinverse(inv_between)
    colnames(cov_between) <- rownames(cov_between) <- Outcomes
    colnames(inv_between) <- rownames(inv_between) <- Outcomes

    Results[["between"]] <- modelCov(
      cov = modelArray(mean = cov_between),
      prec = modelArray(mean = inv_between)
    )
  }

  ### CLUSTER-SPECIFIC within-networks (if random effects) ###
  if (randomeffects != "fixed") {
    nClusters <- nrow(ranef(lmerResults[[1]])[[idvar]])

    if (verbose) {
      message("Computing cluster-specific within-cluster networks")
      pb <- txtProgressBar(min = 0, max = nClusters, style = 3)
    }

    # Get random effects for within predictors:
    # For "correlated": columns 2:(nVar) of ranef (first column is intercept)
    # For "orthogonal": need to match by variable name
    Gamma_within_subject <- within_subject_prec <- within_subject_cov <- vector("list", nClusters)

    for (p in 1:nClusters) {
      Gamma_p <- Gamma_within

      for (i in seq_along(Outcomes)) {
        re_i <- ranef(lmerResults[[i]])[[idvar]]

        if (randomeffects == "correlated") {
          # Columns after intercept are the within random slopes
          re_cols <- withinIDs_ordered[-i]
          re_vals <- unlist(re_i[p, re_cols, drop = TRUE])
        } else {
          # Orthogonal: extract via as.data.frame(VarCorr())
          # ranef columns may have different names
          re_cols <- withinIDs_ordered[-i]
          re_vals <- unlist(re_i[p, re_cols, drop = TRUE])
        }

        Gamma_p[i, -i] <- Gamma_within[i, -i] + re_vals
      }

      Gamma_within_subject[[p]] <- Gamma_p

      # Compute cluster-specific precision/covariance:
      inv_p <- D_within %*% (diag(nVar) - Gamma_p)
      inv_p <- (inv_p + t(inv_p)) / 2
      within_subject_prec[[p]] <- forcePositive(inv_p)
      within_subject_cov[[p]] <- forcePositive(corpcor::pseudoinverse(within_subject_prec[[p]]))

      # Sanity check: if cov trace is too big, fall back to fixed
      if (sum(diag(within_subject_cov[[p]])) > 10 * sum(diag(cov_within))) {
        within_subject_cov[[p]] <- cov_within
        within_subject_prec[[p]] <- inv_within
      }

      if (verbose) {
        setTxtProgressBar(pb, p)
      }
    }
    if (verbose) {
      close(pb)
    }

    # Update Gamma_within with subject-level info:
    Results[["Gamma_within"]] <- modelArray(
      mean = Gamma_within,
      SE = Gamma_within_SE,
      P = Gamma_within_P,
      subject = Gamma_within_subject
    )

    # Update within modelCov with subject-level:
    Results[["within"]] <- modelCov(
      cov = modelArray(mean = cov_within, subject = within_subject_cov),
      prec = modelArray(mean = inv_within, subject = within_subject_prec)
    )

    # Random effect SDs for within Gamma:
    Gamma_within_SDs <- matrix(0, nVar, nVar)
    for (i in seq_along(Outcomes)) {
      if (randomeffects == "correlated") {
        sds <- attr(lme4::VarCorr(lmerResults[[i]])[[idvar]], "stddev")
        # Match names to within predictor IDs (skip intercept)
        Gamma_within_SDs[i, -i] <- sds[withinIDs_ordered[-i]]
      } else {
        # Orthogonal: extract from VarCorr data frame
        df <- as.data.frame(lme4::VarCorr(lmerResults[[i]]))
        Gamma_within_SDs[i, -i] <- df$sdcor[match(withinIDs_ordered[-i], df$var1)]
      }
    }
    colnames(Gamma_within_SDs) <- rownames(Gamma_within_SDs) <- Outcomes

    Results[["Gamma_within"]][["SD"]] <- Gamma_within_SDs
  }

  ### Goodness of fit ###
  names(lmerResults) <- Outcomes

  fit <- data.frame(
    var = Outcomes,
    aic = sapply(lmerResults, AIC),
    bic = sapply(lmerResults, BIC)
  )
  rownames(fit) <- NULL

  ### Combine results ###
  Res <- list(
    results = Results,
    output = lmerResults,
    fit = fit,
    data = augData,
    model = model
  )

  return(Res)
}
