print.mlVAR <- function(x,...){
  name <- deparse(substitute(x))[[1]]
  if (nchar(name) > 10) name <- "object"
  if (name=="x") name <- "object"
  
  est <- x$input$estimator
  
  cat("\nmlVAR estimation completed. Input was:\n",
      "\t- Variables:",x$input$vars,"\n",
      "\t- Lags:",x$input$lags,"\n",
      "\t- Estimator:",est,"\n")
  if (!all(x$input$lags==0)){
    cat("\t- Temporal:",x$input$temporal)
  }
  
  
  ### Tips
  cat("\n\n",
      paste0("Use summary(",name,") to inspect fit and parameter estimates (see ?summary.mlVAR)"),
      "\n",
      paste0("Use plot(",name,") to plot estimated networks (see ?plot.mlVAR)"),
      "\n",
      paste0("Use mlVARcompare(object1, object2) to compare mlVAR models (see ?mlVARcompare)")
  )
  
}


### summary method:
summary.mlVAR <- function(
  object,
  show  = c("fit","temporal","contemporaneous","between"),
  round = 3,
  ...
){
  x <- object
  est <- x$input$estimator
  
  cat("\nmlVAR estimation completed. Input was:\n",
      "\t- Variables:",x$input$vars,"\n",
      "\t- Lags:",x$input$lags,"\n",
      "\t- Estimator:",est,"\n",
      "\t- Temporal:",x$input$temporal)
  
  nVar <- length(object$input$vars)
  nLag <- length(object$input$lags)
  Lags <- object$input$lags
  vars <- object$input$vars
  
  if ("fit" %in% show){
    cat("\n\nInformation indices:\n")
    print(object$fit, row.names=FALSE)
    
    
  }
  
  # Temporal parameters
  if ("temporal" %in% show){
    
    TemporalDF <- data.frame(
      from  = rep(rep(vars,each=nVar),nLag),
      to =  rep(rep(vars,times=nVar),nLag),
      lag = rep(Lags,each=nVar^2),
      fixed = round(c(object$results$Beta$mean),round),  
      SE = round(c(object$results$Beta$SE),round),
      P = round(c(object$results$Beta$P),round),
      ran_SD = round(c(object$results$Beta$SD),round)
    )
    rownames(TemporalDF) <- NULL
    
    cat("\n\nTemporal effects:\n")
    print(TemporalDF,row.names=FALSE)
    
  } else {
    TemporalDF <- NULL
  }
  
  if ("contemporaneous" %in% show){
    
    cor <- object$results$Theta$cor$mean
    corSD <- object$results$Theta$cor$SD
    pcor <- object$results$Theta$pcor$mean
    pcorSD <- object$results$Theta$pcor$SD
    UT <- upper.tri(cor)
    P <- object$results$Gamma_Theta$P
    if (is.null(P)){
      P <- matrix(NA,nrow(UT),ncol(UT))
    }

    cat("\n\nContemporaneous effects (posthoc estimated):\n")
    ContDF <- data.frame(
      node1 = vars[col(cor)][UT],
      node2 = vars[row(cor)][UT],
      "P 1->2" = round(P[UT],round),
      "P 2->1" = round(t(P)[UT],round),
      pcor = round(pcor[UT],round),
      ran_SD_pcor = round(pcorSD[UT],round),
      cor = round(cor[UT],round),
      ran_SD_cor = round(corSD[UT],round)
    )
    names(ContDF) <- c("v1","v2","P 1->2","P 1<-2","pcor","ran_SD_pcor","cor","ran_SD_cor")
    
    print(ContDF,row.names=FALSE)
    
  } else {
    ContDF <- NULL
  }
  
  
  if ("between" %in% show){
    
    P <- object$results$Gamma_Omega_mu$P
    cor <- object$results$Omega_mu$cor$mean
    pcor <- object$results$Omega_mu$pcor$mean
    UT <- upper.tri(cor)
    
    if (is.null(P)){
      P <- matrix(NA,nrow(cor),ncol(cor))
    }
    
    
    cat("\n\nBetween-subject effects:\n")
    BetDF <- data.frame(
      node1 = vars[col(cor)][UT],
      node2 = vars[row(cor)][UT],
      "P 1->2" = round(P[UT],round),
      "P 2->1" = round(t(P)[UT],round),
      pcor = round(pcor[UT],round),
      cor = round(cor[UT],round)
    )
    names(BetDF) <- c("v1","v2","P 1->2","P 1<-2","pcor","cor")
    
    print(BetDF,row.names=FALSE)
    
    
    
  } else {
    BetDF <- NULL
  }
  
  
  invisible(list(temporal = TemporalDF,contemporaneous = ContDF,between = BetDF ))
}


## Internal helper: augment data for prediction (mirrors the augmentation pipeline in mlVAR.R)
.augment_mlVAR_data <- function(data, model, idvar, dayvar, beepvar, scaleWithin, vars, estimator) {
  augData <- data

  # Fill missing beeps:
  beepsPerDay <- dplyr::summarize(
    augData %>% dplyr::group_by(.data[[idvar]], .data[[dayvar]]),
    first = min(.data[[beepvar]], na.rm = TRUE),
    last = max(.data[[beepvar]], na.rm = TRUE),
    .groups = "drop"
  )

  allBeeps <- expand.grid(
    unique(data[[idvar]]),
    unique(data[[dayvar]]),
    seq(min(data[[beepvar]], na.rm = TRUE), max(data[[beepvar]], na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  names(allBeeps) <- c(idvar, dayvar, beepvar)

  allBeeps <- allBeeps %>%
    dplyr::left_join(beepsPerDay, by = c(idvar, dayvar)) %>%
    dplyr::group_by(.data[[idvar]], .data[[dayvar]]) %>%
    dplyr::filter(.data[[beepvar]] >= .data$first, .data[[beepvar]] <= .data$last) %>%
    dplyr::arrange(.data[[idvar]], .data[[dayvar]], .data[[beepvar]])

  augData <- augData %>%
    dplyr::right_join(allBeeps, by = c(idvar, dayvar, beepvar)) %>%
    dplyr::arrange(.data[[idvar]], .data[[dayvar]], .data[[beepvar]])

  # Create predictors from model:
  UniquePredModel <- model[!duplicated(model[, c("pred", "lag", "type")]), c("pred", "lag", "type", "predID")]

  if (!estimator %in% c("Mplus", "JAGS")) {
    for (i in seq_len(nrow(UniquePredModel))) {
      if (UniquePredModel$type[i] == "between") {
        if (estimator == "lmer") {
          augData[[UniquePredModel$predID[i]]] <- ave(
            augData[[UniquePredModel$pred[i]]],
            augData[[idvar]],
            FUN = aveMean
          )
        }
      } else {
        # Lag:
        augData[[UniquePredModel$predID[i]]] <- ave(
          augData[[UniquePredModel$pred[i]]],
          augData[[idvar]], augData[[dayvar]],
          FUN = function(x) aveLag(x, UniquePredModel$lag[i])
        )
        # Center within person:
        augData[[UniquePredModel$predID[i]]] <- ave(
          augData[[UniquePredModel$predID[i]]],
          augData[[idvar]],
          FUN = function(xx) aveCenter(xx, scale = scaleWithin)
        )
      }
    }
  }

  # Within-person standardize dependent vars if scaleWithin = TRUE:
  if (isTRUE(scaleWithin)) {
    for (v in vars) {
      augData[[v]] <- ave(augData[[v]], augData[[idvar]], FUN = function(xx) aveScaleNoCenter(xx))
    }
  }

  # Remove NAs:
  Vars <- unique(c(model$dep, model$predID, idvar, beepvar, dayvar))
  augData <- na.omit(augData[, Vars])

  return(augData)
}


## Internal helper: predict on newdata
.predict_mlVAR_newdata <- function(object, newdata, scale_back = TRUE) {

  vars <- object$input$vars
  idvar <- object$input$idvar
  dayvar <- object$input$dayvar
  beepvar <- object$input$beepvar
  estimator <- object$input$estimator

  # Auto-create DAY and BEEP if they were auto-generated during fitting:
  if (!dayvar %in% names(newdata)) {
    newdata[[dayvar]] <- 1
  }
  if (!beepvar %in% names(newdata)) {
    newdata[[beepvar]] <- ave(seq_len(nrow(newdata)), newdata[[idvar]], newdata[[dayvar]], FUN = seq_along)
  }

  # Validate newdata columns:
  required_cols <- c(idvar, dayvar, beepvar, vars)
  missing_cols <- required_cols[!required_cols %in% names(newdata)]
  if (length(missing_cols) > 0) {
    stop(paste0("The following columns are missing from newdata: ", paste(missing_cols, collapse = ", ")))
  }

  # Store original newdata for output:
  origNewdata <- newdata[, required_cols, drop = FALSE]

  # Apply full_detrend if model used it:
  if (isTRUE(object$input$full_detrend)) {
    obs_idx <- ave(seq_len(nrow(newdata)), newdata[[idvar]], FUN = seq_along)
    obs_per_cluster <- table(newdata[[idvar]])
    if (length(unique(obs_per_cluster)) != 1) {
      warning("'full_detrend' was used in model fitting but newdata has unequal observations per cluster. Detrending not applied to newdata.")
    } else {
      for (v in vars) {
        newdata[[v]] <- ave(newdata[[v]], obs_idx, FUN = Scale)
      }
    }
  }

  # Apply grand-mean scaling (using training set parameters):
  if (isTRUE(object$input$scaled)) {
    for (v in vars) {
      newdata[[v]] <- (newdata[[v]] - object$input$scale_means[v]) / object$input$scale_sds[v]
    }
  }

  # Augment newdata (fill missing beeps, create lagged predictors, center):
  augData <- .augment_mlVAR_data(
    data = newdata,
    model = object$model,
    idvar = idvar,
    dayvar = dayvar,
    beepvar = beepvar,
    scaleWithin = isTRUE(object$input$scaleWithin),
    vars = vars,
    estimator = estimator
  )

  # Initialize result matrices aligned to original newdata:
  nOrig <- nrow(origNewdata)
  predicted_df <- as.data.frame(matrix(NA_real_, nrow = nOrig, ncol = length(vars)))
  colnames(predicted_df) <- vars
  residuals_df <- as.data.frame(matrix(NA_real_, nrow = nOrig, ncol = length(vars)))
  colnames(residuals_df) <- vars

  # Match augData rows back to original newdata:
  aug_key <- paste(augData[[idvar]], augData[[dayvar]], augData[[beepvar]], sep = "_")
  orig_key <- paste(origNewdata[[idvar]], origNewdata[[dayvar]], origNewdata[[beepvar]], sep = "_")
  idx <- match(orig_key, aug_key)
  matched <- !is.na(idx)

  if (estimator == "lmer") {
    lmerResults <- object$output$temporal

    for (i in seq_along(vars)) {
      pred_vals <- predict(lmerResults[[i]], newdata = augData, allow.new.levels = TRUE)
      obs_vals <- augData[[vars[i]]]
      resid_vals <- obs_vals - pred_vals

      predicted_df[matched, vars[i]] <- pred_vals[idx[matched]]
      residuals_df[matched, vars[i]] <- resid_vals[idx[matched]]
    }

  } else if (estimator == "lm") {
    # lm: can only predict for subjects seen during training
    train_IDs <- unique(object$data[[idvar]])
    new_IDs <- unique(origNewdata[[idvar]])
    unknown_IDs <- new_IDs[!new_IDs %in% train_IDs]
    if (length(unknown_IDs) > 0) {
      warning(paste0("The following subject IDs in newdata were not in the training data ",
                     "and cannot be predicted with estimator = 'lm': ",
                     paste(unknown_IDs, collapse = ", ")))
    }

    for (i in seq_along(vars)) {
      pred_all <- rep(NA_real_, nrow(augData))

      for (p in seq_along(train_IDs)) {
        id <- train_IDs[p]
        sub_rows <- which(augData[[idvar]] == id)
        if (length(sub_rows) == 0) next

        sub <- augData[sub_rows, , drop = FALSE]
        lm_obj <- object$output[[i]][[p]]
        pred_vals <- predict(lm_obj, newdata = sub)
        pred_all[sub_rows] <- as.numeric(pred_vals)
      }

      obs_vals <- augData[[vars[i]]]
      resid_vals <- obs_vals - pred_all

      predicted_df[matched, vars[i]] <- pred_all[idx[matched]]
      residuals_df[matched, vars[i]] <- resid_vals[idx[matched]]
    }

  } else {
    stop(paste0("predict.mlVAR with newdata not implemented for estimator = '", estimator, "'"))
  }

  # Scale back to original metric if requested:
  if (scale_back && isTRUE(object$input$scaled)) {
    for (v in vars) {
      predicted_df[[v]] <- predicted_df[[v]] * object$input$scale_sds[v] + object$input$scale_means[v]
      residuals_df[[v]] <- residuals_df[[v]] * object$input$scale_sds[v]
    }
  }

  list(
    predicted = predicted_df,
    residuals = residuals_df,
    ids = origNewdata[, c(idvar, dayvar, beepvar), drop = FALSE]
  )
}


## Internal workhorse: computes both predictions and residuals on training data
.predict_mlVAR_train <- function(object, scale_back = TRUE) {

  vars <- object$input$vars
  idvar <- object$input$idvar
  dayvar <- object$input$dayvar
  beepvar <- object$input$beepvar
  estimator <- object$input$estimator

  origData <- object$input$originalData
  augData <- object$data

  nOrig <- nrow(origData)
  predicted_df <- as.data.frame(matrix(NA_real_, nrow = nOrig, ncol = length(vars)))
  colnames(predicted_df) <- vars
  residuals_df <- as.data.frame(matrix(NA_real_, nrow = nOrig, ncol = length(vars)))
  colnames(residuals_df) <- vars

  # Match augData rows to originalData via composite key:
  aug_key <- paste(augData[[idvar]], augData[[dayvar]], augData[[beepvar]], sep = "_")
  orig_key <- paste(origData[[idvar]], origData[[dayvar]], origData[[beepvar]], sep = "_")
  idx <- match(orig_key, aug_key)
  matched <- !is.na(idx)

  if (estimator == "lmer") {
    lmerResults <- object$output$temporal

    for (i in seq_along(vars)) {
      fit_vals <- fitted(lmerResults[[i]])
      resid_vals <- residuals(lmerResults[[i]])

      # Align to augData rows (names = rownames of data used in lmer):
      aug_fit <- rep(NA_real_, nrow(augData))
      aug_res <- rep(NA_real_, nrow(augData))
      aug_positions <- match(names(fit_vals), rownames(augData))
      aug_fit[aug_positions] <- as.numeric(fit_vals)
      aug_res[aug_positions] <- as.numeric(resid_vals)

      predicted_df[matched, vars[i]] <- aug_fit[idx[matched]]
      residuals_df[matched, vars[i]] <- aug_res[idx[matched]]
    }

  } else if (estimator == "lm") {
    IDs <- unique(augData[[idvar]])

    for (i in seq_along(vars)) {
      aug_fit <- rep(NA_real_, nrow(augData))
      aug_res <- rep(NA_real_, nrow(augData))

      for (p in seq_along(IDs)) {
        lm_obj <- object$output[[i]][[p]]
        fit_vals <- fitted(lm_obj)
        resid_vals <- stats::residuals(lm_obj)

        # Names = rownames from augData subset for this subject:
        aug_positions <- match(names(fit_vals), rownames(augData))
        aug_fit[aug_positions] <- as.numeric(fit_vals)
        aug_res[aug_positions] <- as.numeric(resid_vals)
      }

      predicted_df[matched, vars[i]] <- aug_fit[idx[matched]]
      residuals_df[matched, vars[i]] <- aug_res[idx[matched]]
    }

  } else {
    stop(paste0("predict/residuals not implemented for estimator = '", estimator, "'"))
  }

  # Scale back to original metric if requested:
  if (scale_back && isTRUE(object$input$scaled)) {
    for (v in vars) {
      predicted_df[[v]] <- predicted_df[[v]] * object$input$scale_sds[v] + object$input$scale_means[v]
      residuals_df[[v]] <- residuals_df[[v]] * object$input$scale_sds[v]
    }
  }

  list(predicted = predicted_df, residuals = residuals_df)
}


predict.mlVAR <- function(object, newdata, scale_back = TRUE, include_ids = TRUE, ...) {

  vars <- object$input$vars
  idvar <- object$input$idvar
  dayvar <- object$input$dayvar
  beepvar <- object$input$beepvar

  # Dispatch to newdata handler:
  if (!missing(newdata)) {
    res <- .predict_mlVAR_newdata(object, newdata, scale_back = scale_back)
    result <- res$predicted
    if (include_ids) {
      result <- cbind(res$ids, result)
    }
    return(result)
  }

  # Training data:
  res <- .predict_mlVAR_train(object, scale_back = scale_back)
  result <- res$predicted

  if (include_ids) {
    ids <- object$input$originalData[, c(idvar, dayvar, beepvar), drop = FALSE]
    result <- cbind(ids, result)
  }

  return(result)
}


residuals.mlVAR <- function(object, scale_back = TRUE, include_ids = TRUE, ...) {

  vars <- object$input$vars
  idvar <- object$input$idvar
  dayvar <- object$input$dayvar
  beepvar <- object$input$beepvar

  res <- .predict_mlVAR_train(object, scale_back = scale_back)
  result <- res$residuals

  if (include_ids) {
    ids <- object$input$originalData[, c(idvar, dayvar, beepvar), drop = FALSE]
    result <- cbind(ids, result)
  }

  return(result)
}


### mlGGM S3 methods ###

print.mlGGM <- function(x, ...) {
  name <- deparse(substitute(x))[[1]]
  if (nchar(name) > 10) name <- "object"
  if (name == "x") name <- "object"

  cat("\nmlGGM estimation completed. Input was:\n",
      "\t- Variables:", x$input$vars, "\n",
      "\t- Estimator:", x$input$estimator, "\n",
      "\t- Random effects:", x$input$randomeffects, "\n")

  cat("\n",
      paste0("Use summary(", name, ") to inspect fit and parameter estimates"),
      "\n",
      paste0("Use plot(", name, ") to plot estimated networks"),
      "\n",
      paste0("Use getNet(", name, ", 'within') or getNet(", name, ", 'between') to extract network matrices"),
      "\n")
}


summary.mlGGM <- function(
  object,
  show = c("fit", "within", "between"),
  round = 3,
  ...
) {
  x <- object

  cat("\nmlGGM estimation completed. Input was:\n",
      "\t- Variables:", x$input$vars, "\n",
      "\t- Estimator:", x$input$estimator, "\n",
      "\t- Random effects:", x$input$randomeffects, "\n")

  nVar <- length(object$input$vars)
  vars <- object$input$vars

  if ("fit" %in% show) {
    cat("\nInformation indices:\n")
    print(object$fit, row.names = FALSE)
  }

  if ("within" %in% show) {
    pcor <- object$results$within$pcor$mean
    cor_mat <- object$results$within$cor$mean
    P <- object$results$Gamma_within$P
    UT <- upper.tri(pcor)

    if (is.null(P)) {
      P <- matrix(NA, nrow(UT), ncol(UT))
    }

    cat("\n\nWithin-cluster effects:\n")
    WithinDF <- data.frame(
      v1 = vars[col(pcor)][UT],
      v2 = vars[row(pcor)][UT],
      "P 1->2" = round(P[UT], round),
      "P 2->1" = round(t(P)[UT], round),
      pcor = round(pcor[UT], round),
      cor = round(cor_mat[UT], round)
    )
    names(WithinDF) <- c("v1", "v2", "P 1->2", "P 1<-2", "pcor", "cor")
    print(WithinDF, row.names = FALSE)
  } else {
    WithinDF <- NULL
  }

  if ("between" %in% show) {
    pcor <- object$results$between$pcor$mean
    cor_mat <- object$results$between$cor$mean
    P <- object$results$Gamma_between$P
    UT <- upper.tri(pcor)

    if (is.null(P)) {
      P <- matrix(NA, nrow(pcor), ncol(pcor))
    }

    cat("\n\nBetween-cluster effects:\n")
    BetDF <- data.frame(
      v1 = vars[col(pcor)][UT],
      v2 = vars[row(pcor)][UT],
      "P 1->2" = round(P[UT], round),
      "P 2->1" = round(t(P)[UT], round),
      pcor = round(pcor[UT], round),
      cor = round(cor_mat[UT], round)
    )
    names(BetDF) <- c("v1", "v2", "P 1->2", "P 1<-2", "pcor", "cor")
    print(BetDF, row.names = FALSE)
  } else {
    BetDF <- NULL
  }

  invisible(list(within = WithinDF, between = BetDF))
}


plot.mlGGM <- function(
  x,
  type = c("within", "between"),
  partial = TRUE,
  SD = FALSE,
  subject,
  order,
  nonsig = c("default", "show", "hide", "dashed"),
  rule = c("or", "and"),
  alpha = 0.05,
  layout = "spring",
  verbose = TRUE,
  ...
) {
  rule <- match.arg(rule)
  type <- match.arg(type)
  nonsig <- match.arg(nonsig)

  if (nonsig == "default") {
    if (!partial || !missing(subject)) {
      nonsig <- "show"
    } else {
      nonsig <- "hide"
    }
    if (verbose) {
      message(paste0("'nonsig' argument set to: '", nonsig, "'"))
    }
  }

  if (missing(order)) {
    order <- x$input$vars
  }

  if (is.character(order)) {
    ord <- match(order, x$input$vars)
  } else {
    ord <- order
  }

  sub <- ifelse(partial, "pcor", "cor")

  if (type == "within") {
    if (SD) {
      NET <- x$results$within[[sub]]$SD
      SIG <- matrix(TRUE, nrow(NET), ncol(NET))
    } else if (!missing(subject)) {
      NET <- x$results$within[[sub]]$subject[[subject]]
      SIG <- matrix(TRUE, nrow(NET), ncol(NET))
      if (nonsig != "show") {
        warning("Can not hide non-significant edges for subject network.")
      }
    } else {
      NET <- x$results$within[[sub]]$mean

      if (nonsig != "show") {
        if (partial && !is.null(x$results$Gamma_within) && !all(is.nan(x$results$Gamma_within$P))) {
          diag(x$results$Gamma_within$P) <- 0
          if (rule == "or") {
            SIG <- x$results$Gamma_within$P < alpha | t(x$results$Gamma_within$P) < alpha
          } else {
            SIG <- x$results$Gamma_within$P < alpha & t(x$results$Gamma_within$P) < alpha
          }
        } else if (!any(is.na(x$results$within[[sub]]$P))) {
          SIG <- x$results$within[[sub]]$P < alpha
        } else {
          SIG <- matrix(TRUE, nrow(NET), ncol(NET))
          if (nonsig != "show") {
            stop("No p-values or CI computed. Can not hide non-significant edges.")
          }
        }
      } else {
        SIG <- matrix(TRUE, nrow(NET), ncol(NET))
      }
    }
    NET <- makeSym(NET)
  }

  if (type == "between") {
    if (!missing(subject)) {
      stop("No subject-specific between network possible")
    }
    if (SD) {
      stop("No SD for between-cluster network.")
    }

    NET <- x$results$between[[sub]]$mean

    if (nonsig != "show") {
      if (partial && !is.null(x$results$Gamma_between) && !all(is.nan(x$results$Gamma_between$P))) {
        diag(x$results$Gamma_between$P) <- 0
        if (rule == "or") {
          SIG <- x$results$Gamma_between$P < alpha | t(x$results$Gamma_between$P) < alpha
        } else {
          SIG <- x$results$Gamma_between$P < alpha & t(x$results$Gamma_between$P) < alpha
        }
      } else if (!any(is.na(x$results$between[[sub]]$P))) {
        SIG <- x$results$between[[sub]]$P < alpha
      } else {
        SIG <- matrix(TRUE, nrow(NET), ncol(NET))
        if (nonsig != "show") {
          stop("No p-values or CI computed. Can not hide non-significant edges.")
        }
      }
    } else {
      SIG <- matrix(TRUE, nrow(NET), ncol(NET))
    }

    NET <- makeSym(NET)
  }

  ### PLOT NETWORK ###
  if (nonsig == "dashed") {
    lty <- ifelse(!SIG, 2, 1)
  } else {
    lty <- 1
  }

  if (nonsig == "hide") {
    NET <- NET * SIG
  }

  if (any(is.na(NET[ord, ord][upper.tri(NET[ord, ord])])) || any(is.na(NET[ord, ord][lower.tri(NET[ord, ord])]))) {
    stop("Network not estimated correctly.")
  }

  qgraph::qgraph(NET[ord, ord], lty = lty, labels = x$input$vars[ord],
                  layout = layout, ..., directed = FALSE)
}


makeSym <- function(x) (x + t(x))/2


getNet <- function(x,...){
  qgraph::getWmat(plot(x,...,DoNotPlot=TRUE))
}

plot.mlVARsim <- function(x,...){
  x$results <- x$model
  x$input <- list(vars = x$vars)
  class(x) <- "mlVAR"
  plot.mlVAR(x,...,nonsig="show")
}

plot.mlVAR <- 
  function(x, # mlVAR object
           type = c("temporal","contemporaneous","between"), # # also allows for partial matching. e.g., temp or t.
           lag = 1, # lag of temporal network
           partial = TRUE, # Show partial correlations?
           SD = FALSE, # Plots SD instead of normal parameters
           subject, # If assigned, show the network of a particulair subject instead
           order, # If assigned, re-order nodes
           nonsig = c("default","show","hide","dashed"), # How to handle nonsignificant edges? In Bayesian estimation, checks if 0 is inside interval.
           rule = c("or", "and"), # GGM sig rule
           alpha = 0.05, # alpha value for significance test
           onlySig = FALSE, # Backward competability argument.
           layout = "spring",
           verbose = TRUE,
           ...  #Arguments sent to qgraph
  ){
    rule <- match.arg(rule)
    
    # First some backward competability:
    if (type[[1]] == "fixed"){
      warning("type = 'fixed' is deprecated. type = 'temporal' instead.")
      type <- "temporal"
    }
    if (type[[1]] == "SD"){
      warning("type = 'SD' is deprecated. type = 'temporal' and SD = TRUE instead.")
      type <- "temporal"
      SD <- TRUE
    }
    if (type[[1]] == "subject"){
      warning("type = 'subject' is deprecated. type = 'temporal' and subject = ... instead")
      type <- "temporal"
      if (missing(subject)) stop("'subject' must be assigned")
    }
    if (onlySig){
      warning("'onlySig' is deprecated. Setting nonsig = 'hide'.")
      nonsig <- "hide"
    }
    
    
    # Now check for arguments:
    type <- match.arg(type)
    nonsig <- match.arg(nonsig)
    if (nonsig == "default"){
      if (!SD){
        nonsig <- "hide"
        if (!partial){
          nonsig <- "show"
        }
        if (!missing(subject)){
          nonsig <- "show"
        }
      } else {
        nonsig <- 'show'
      }

      if (verbose){
        message(paste0("'nonsig' argument set to: '",nonsig,"'"))
      }
    }
    
    if (missing(order)){
      order <- x$input$vars
    }
    
    if (!missing(subject) && SD){
      stop("'SD' not available for subject.")
    }
    
    # If order is character, find ord as number. Else just set ord to order:
    if (is.character(order)){
      ord <- match(order,x$input$vars)
    } else {
      ord <- order
    }
    
    ### OBTAIN NETWORK TO PLOT ###
    
    # Temporal:
    if (type == "temporal"){
      
      if (SD){
        SIG <- matrix(TRUE,length(ord),length(ord))
        NET <- t(x$results$Beta$SD[,,lag])
      } else {
        
        if (missing(subject)){
          # Obtain fixed effects network:
#           if (SD){
#             
#           } else {
            NET <- t(x$results$Beta$mean[,,lag])  
          # }
          
          # Attempt to obtain significance:
          # Via P:
          
          if (any(is.na(x$results$Beta$P))){

            # Via CI:
            if (!any(is.na(x$results$Beta$lower)) && !any(is.na(x$results$Beta$upper))){
              SIG <- t(x$results$Beta$lower[,,lag]) > 0 |  t(x$results$Beta$upper[,,lag]) < 0
            } else {
              # No P or CI:
              SIG <- matrix(TRUE, nrow(NET), ncol(NET))
              if (nonsig != "show"){
                warning("No p-values or CI computed. Can not hide non-significant edges.")
              }
            }
          } else {
            SIG <-  t(x$results$Beta$P[,,lag]) < alpha
          }
          
        } else {
          # Obtain subject network:
          NET <- t(x$results$Beta$subject[[subject]][,,lag])
          SIG <- matrix(TRUE, nrow(NET), ncol(NET))
          if (nonsig != "show"){
            warning("Can not hide non-significant edges for subject network.")
          }
        } 
      }
    }
    
    # Contemporaneous:
    if (type == "contemporaneous"){
      
      sub <- ifelse(partial,"pcor","cor")
      
      if (missing(subject)){
        if (SD){
          NET <- x$results$Theta[[sub]]$SD
        } else {
          # Obtain fixed effects network:
          NET <- x$results$Theta[[sub]]$mean
        }
        
        # Attempt to obtain significance:
        # Via P:
        
        # If nonsig is not show:
        if (nonsig != "show"){
          
          # Stop if SD:
          if (SD){
            stop("No significance for SD")
          }
          
          # Not partial:
          if (!any(is.na(x$results$Theta[[sub]]$lower)) && !any(is.na(x$results$Theta[[sub]]$upper))){
            SIG <- x$results$Theta[[sub]]$lower > 0 |  x$results$Theta[[sub]]$upper < 0
          } else if (!any(is.na(x$results$Theta[[sub]]$P))){
            SIG <-  x$results$Theta[[sub]]$P < alpha
          } else {
            
            # If partial, try via gamma:
            if (partial && !is.null(x$results$Gamma_Theta) && !all(is.nan(x$results$Gamma_Theta$P))){
              diag(x$results$Gamma_Theta$P) <- 0
              
              if (rule == "or"){
                SIG <- x$results$Gamma_Theta$P < alpha | t(x$results$Gamma_Theta$P ) < alpha
              } else {
                SIG <- x$results$Gamma_Theta$P < alpha & t(x$results$Gamma_Theta$P ) < alpha
              }
            } else {
              # No P or CI:
              SIG <- matrix(TRUE, nrow(NET), ncol(NET))
              if (nonsig != "show"){
                stop("No p-values or CI computed. Can not hide non-significant edges.")
              }              
            }
          }
        } else {
          SIG <- matrix(TRUE, nrow(NET), ncol(NET))
        }
        
        #         if (!SD && any(is.na(x$results$Theta[[sub]]$P))){
        #           
        #           # Via CI:
        #           if (!any(is.na(x$results$Theta[[sub]]$lower)) && !any(is.na(x$results$Theta[[sub]]$upper))){
        #             SIG <- x$results$Theta[[sub]]$lower > 0 |  x$results$Theta[[sub]]$upper < 0
        #           } else {
        #             # No P or CI:
        #             SIG <- matrix(TRUE, nrow(NET), ncol(NET))
        #             if (nonsig != "show"){
        #               warning("No p-values or CI computed. Can not hide non-significant edges.")
        #             }
        #           }
        #         } else {
        #           SIG <-  x$results$Theta[[sub]]$P < alpha
        #         }
        
      } else {
        # Obtain subject network:
        NET <- x$results$Theta[[sub]]$subject[[subject]]
        SIG <- matrix(TRUE, nrow(NET), ncol(NET))
        if (nonsig != "show"){
          warning("Can not hide non-significant edges for subject network.")
        }
      }
      NET <- makeSym(NET)
    }
    
    # Between-subjects:
    if (type == "between"){
      
      sub <- ifelse(partial,"pcor","cor")
      
      if (!missing(subject)){
        stop("No subject-specific between network possible")
      }
      if (SD){
        stop("No SD for between-subjects network.")
      }
      
      # Obtain fixed effects network:
      NET <- x$results$Omega_mu[[sub]]$mean
      
      #       # Attempt to obtain significance:
      #       # Via P:
      #       if (any(is.na(x$results$Omega_mu[[sub]]$P))){
      #         
      #         # Via CI:
      #         if (!any(is.na(x$results$Omega_mu[[sub]]$lower)) && !any(is.na(x$results$Omega_mu[[sub]]$upper))){
      #           SIG <- x$results$Omega_mu[[sub]]$lower > 0 |  x$results$Omega_mu[[sub]]$upper < 0
      #         } else {
      #           # No P or CI:
      #           SIG <- matrix(TRUE, nrow(NET), ncol(NET))
      #           if (nonsig != "show"){
      #             warning("No p-values or CI computed. Can not hide non-significant edges.")
      #           }
      #         }
      #       } else {
      #         SIG <-  x$results$Omega_mu[[sub]]$P < alpha
      #       }
      #       
      
      # If nonsig is not show:
      if (nonsig != "show"){
        
        # Stop if SD:
        if (SD){
          stop("No significance for SD")
        }
        
        # Not partial:
        if (!any(is.na(x$results$Omega_mu[[sub]]$lower)) && !any(is.na(x$results$Omega_mu[[sub]]$upper))){
          SIG <- x$results$Omega_mu[[sub]]$lower > 0 |  x$results$Omega_mu[[sub]]$upper < 0
        } else if (!any(is.na(x$results$Omega_mu[[sub]]$P))){
          SIG <-  x$results$Omega_mu[[sub]]$P < alpha
        } else {
          
          # If partial, try via gamma:
          if (partial && !is.null(x$results$Gamma_Omega_mu) && !all(is.nan(x$results$Gamma_Omega_mu$P))){
            diag(x$results$Gamma_Omega_mu$P) <- 0
            
            if (rule == "or"){
              SIG <- x$results$Gamma_Omega_mu$P < alpha | t(x$results$Gamma_Omega_mu$P ) < alpha
            } else {
              SIG <- x$results$Gamma_Omega_mu$P < alpha & t(x$results$Gamma_Omega_mu$P ) < alpha
            }
          } else {
            # No P or CI:
            SIG <- matrix(TRUE, nrow(NET), ncol(NET))
            if (nonsig != "show"){
              stop("No p-values or CI computed. Can not hide non-significant edges.")
            }              
          }
        }
      } else {
        SIG <- matrix(TRUE, nrow(NET), ncol(NET))
      }
      
      NET <- makeSym(NET)
    }
    
    ### PLOT NETWORK ###
    if (nonsig == "dashed"){
      lty <- ifelse(!SIG,2,1)
    } else {
      lty <- 1
    }
    
    if (nonsig == "hide"){
      NET <- NET * SIG
    }
    
    if (any(is.na(NET[ord,ord][upper.tri(NET[ord,ord])])) || any(is.na(NET[ord,ord][lower.tri(NET[ord,ord])]))){
      stop("Network not estimated correctly.")
    }
    qgraph::qgraph(NET[ord,ord],lty = lty,labels = x$input$vars[ord],layout=layout,
                   ..., directed = type == "temporal")
  }




# # # logLik function:
# # logLik.mlVAR <- function(object){
# #   res <- object$pseudologlik
# #   class(res) <- "logLik"
# #   attr(res, "df") <- object$df
# #   return(res)
# # } 
# 
# # Plot function:
# # plot.mlVAR_MW <- function(x, lag = 1, ...){
# #   fixef <- fixedEffects(x)
# #   
# #   # Extract only lagged variables:
# #   sub <- fixef %>% filter(grepl(paste0("^L",lag,"_"), Predictor))
# #   
# #   # make matrix:
# #   Nodes <- as.character(unique(sub$Response))
# #   nNode <- length(Nodes)
# #   Network <- matrix(0, nNode, nNode)
# #   for (i in seq_along(Nodes)){
# #     Network[,i] <- sub$effect[sub$Response==Nodes[i]][match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes)]
# #   }
# #   
# #   Graph <- qgraph(Network, labels=Nodes, ...)
# #   invisible(Graph)
# # }
# 
# ### 
# tab2net <- function(x,lag=1){
#   Nodes <- x$dep
#   x <- x[,grepl(paste0("^L",lag,"_"), names(x))]
#   x[is.na(x)] <- 0
#   x <- t(x)
#   colnames(x) <- rownames(x) <- Nodes
#   x
# }
# 
# 
# #### These functions probably should be the other way around... but this works...
# getNet <- function(x, ...){
#   qgraph:::getWmat(plot(x,...,DoNotPlot=TRUE))
# }
# 
# plot.mlVAR <- function(x, type = c("fixed","SD","subject"), lag = 1,subject,order,onlySig = FALSE,
#                 alpha, # alpha value. if missing, do bonferonni on 0.05.
#                 ...){
#   type <- match.arg(type)
#   
#   if (type[[1]]=="subject" & missing(subject)){
#     stop("'subject' is needed to plot individual network")
#   }
# 
#   # Nodes <- rownames(x$fixedEffects)
#   if (type[[1]]=="fixed"){
# 
#     fixef <- fixedEffects(x)
#     Nodes <- as.character(unique(fixef$Response))
# 
#     # Extract only lagged variables:
#     sub <- fixef %>% filter(grepl(paste0("^L",lag,"_"), Predictor))
#     
#     # make matrix:
#     ### Simple trick to remove nonsigs, edit fixed effects to be 0 if not sig:
#     if (onlySig){
#       if (missing(alpha)){
#         alpha <- 0.05 / nrow(sub)
#       }
#       sub <- sub %>% ungroup %>% mutate(effect = ifelse(p < alpha, effect, 0))
#     }
# 
#     if (!missing(order)){
#       if (!all(sort(Nodes)==sort(order)))stop("'order' must contain exact node labels")
#       Nodes <- Nodes[match(order,Nodes)]
#     }
#     nNode <- length(Nodes)
#     Network <- matrix(0, nNode, nNode)
# 
#     for (i in seq_along(Nodes)){
#       # Network[,i] <- sub$effect[sub$Response==Nodes[i]][match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes)]
#       Network[ match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes),i] <- sub$effect[sub$Response==Nodes[i]]
#     }
#     
#     Graph <- qgraph(Network, labels=Nodes, ...)
# #   } else if (type[[1]]=="se"){
# #     fixef <- fixedEffects(x)
# #     Nodes <- as.character(unique(fixef$Response))
# #     
# #     # Extract only lagged variables:
# #     sub <- fixef %>% filter(grepl(paste0("^L",lag,"_"), Predictor))
# #     
# #     # make matrix:
# # 
# #     if (!missing(order)){
# #       if (!all(sort(Nodes)==sort(order)))stop("'order' must contain exact node labels")
# #       Nodes <- Nodes[match(order,Nodes)]
# #     }
# #     nNode <- length(Nodes)
# #     Network <- matrix(0, nNode, nNode)
# #     for (i in seq_along(Nodes)){
# #       # Network[,i] <- sub$effect[sub$Response==Nodes[i]][match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes)]
# #       Network[ match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes),i] <- sub$se[sub$Response==Nodes[i]]
# #     }
# # 
# #     
# #     Graph <- qgraph(Network, labels=Nodes, ...)
# #     
#   } else if (type[[1]]=="SD"){
#   if (onlySig){
#     warning("'onlySig' argument can only be used with type = 'fixed'")
#   }
#     ranef <- randomEffects(x)
#     Nodes <- as.character(unique(ranef$Response))
#     
#     # Extract only lagged variables:
#     sub <- ranef %>% filter(grepl(paste0("^L",lag,"_"), Predictor))
# 
#     # make matrix:
#     nNode <- length(Nodes)
#     Network <- matrix(0, nNode, nNode)
# 
#     for (i in seq_along(Nodes)){
#     Network[ match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes),i] <- sub$variance[sub$Response==Nodes[i]]
#     }
#     
#     Graph <- qgraph(sqrt(Network), labels=Nodes, ...)
#     
#   }  else if (type[[1]]=="subject"){
#     if (onlySig){
#       warning("'onlySig' argument can only be used with type = 'fixed'")
#     }
#     Net <- tab2net(x$randomEffects[[subject]])
#     fixed <- getNet(x,"fixed",onlySig=FALSE)
#     
#     Graph <- qgraph(fixed+Net, labels=rownames(Net), ...)
#      
#   } else stop("'type' is not supported")
#   
#   invisible(Graph)
# }
# 
# 
# # Print and summary:
# summary.mlVAR <- function(object,...){
#   input <- object$input
#   
#   inputstr <- paste(sapply(seq_along(input),function(i)paste0(names(input)[i],":\t\t",paste(input[[i]],collapse=", "))), collapse = "\n")
#   
#   cat(paste0("==== mlVAR results ====\n",inputstr,"\n\nNumber of random effects:\t\t",length(object$randomEffects),"\n",
#              "pseudo log-likeligood:\t\t",round(object$pseudologlik,2),"\n",
#              "Degrees of Freedom:\t\t",round(object$df,2),"\n",
#              "BIC:\t\t",round(object$BIC,2)             
#              ))
# }
# 
# print.mlVAR <- function(x,...) summary.mlVAR(x,...)
# 
# 
# ### Model plot method:
# 
# plot.mlVARsim0 <- function(x, type = c("fixed","SD","subject"), lag = 1,subject,order,...){
#   if (type[[1]]=="subject" & missing(subject)){
#     stop("'subject' is needed to plot individual network")
#   }
#   
#   # Nodes <- rownames(x$fixedEffects)
#   if (type[[1]]=="fixed"){
# 
# 
#     fixef <- x$fixedEffects
#     Nodes <- x$vars
#     
#     Graph <- qgraph(t(fixef), labels=Nodes, ...)
#     
#   } else if (type[[1]]=="se"){
#     
#     stop("No standard errors in true model")
# 
#     
#   } else if (type[[1]]=="SD"){
#     
#     ranef <- x$randomEffectsSD
#     Nodes <- x$vars
# 
#     Graph <- qgraph(t(ranef), labels=Nodes, ...)
#     
#   }  else if (type[[1]]=="subject"){
#     fixef <- x$fixedEffects
#     Nodes <- x$vars
#     Net <- fixef + x$randomEffects[[subject]]
#     Graph <- qgraph(t(Net), labels=Nodes, ...)
#     
#   } else stop("'type' is not supported")
#   
#   invisible(Graph)
# }
