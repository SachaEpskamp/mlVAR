# mlGGM: Multi-Level Gaussian Graphical Model
# Estimates within-cluster and between-cluster partial correlation networks
# from cross-sectional multilevel data using a single-step nodewise regression.

mlGGM <- function(
  data,        # Data frame
  vars,        # Vector of variable names to include
  idvar = "id", # String indicating the cluster id variable name
  estimator = c("lmer"),
  randomeffects = c("default", "correlated", "orthogonal", "fixed"),
  scale = TRUE,
  na.rm = TRUE,
  verbose = TRUE
) {

  # Experimental warning:
  warning("mlGGM is currently experimental and may change in future versions.")

  # Match args:
  estimator <- match.arg(estimator)
  randomeffects <- match.arg(randomeffects)

  # Default random effects: correlated if <=6 vars, else orthogonal
  if (randomeffects == "default") {
    if (length(vars) <= 6) {
      randomeffects <- "correlated"
    } else {
      randomeffects <- "orthogonal"
    }
    if (verbose) {
      message(paste0("'randomeffects' set to '", randomeffects, "' (default for ", length(vars), " variables)"))
    }
  }

  # Validate vars:
  if (missing(vars) || length(vars) < 2) {
    stop("'vars' must contain at least 2 variable names.")
  }
  if (!all(vars %in% names(data))) {
    missing_vars <- vars[!vars %in% names(data)]
    stop(paste0("The following variables were not found in data: ", paste(missing_vars, collapse = ", ")))
  }

  # Validate idvar:
  if (!is.character(idvar) || length(idvar) != 1 || !idvar %in% names(data)) {
    stop("'idvar' must be a string indicating a column name of the data.")
  }

  # Remove rows with NA in idvar:
  data <- data[!is.na(data[[idvar]]), ]

  ### QR rank-deficiency check ###
  if (na.rm) {
    X <- as.matrix(na.omit(data[, vars]))
  } else {
    X <- as.matrix(data[, vars])
  }

  qrX <- qr(X)
  rnk <- qrX$rank

  if (rnk < length(vars)) {
    keep <- qrX$pivot[1:rnk]
    discard <- vars[!seq_along(vars) %in% keep]

    warning(paste0("The following variables are linearly dependent on other columns, and therefore dropped from the mlGGM analysis:\n",
                   paste("\t-", discard, collapse = "\n")))
    vars <- vars[keep]
  }

  # Grand-mean standardize if scale = TRUE:
  if (scale) {
    for (v in vars) {
      data[[v]] <- Scale(data[[v]])
    }
  }

  ### Build predictor model ###
  # For each (dep, pred) pair where dep != pred, we have two types:
  #   - "within": centered predictor (y_j - cluster_mean(y_j))
  #   - "between": cluster-mean predictor (cluster_mean(y_j))
  PredModel <- rbind(
    expand.grid(
      dep = vars,
      pred = vars,
      type = "within",
      stringsAsFactors = FALSE
    ),
    expand.grid(
      dep = vars,
      pred = vars,
      type = "between",
      stringsAsFactors = FALSE
    )
  )
  # Remove self-predictions:
  PredModel <- PredModel[PredModel$dep != PredModel$pred, ]

  # Unique predictors:
  UniquePredModel <- PredModel[!duplicated(PredModel[, c("pred", "type")]), c("pred", "type")]
  UniquePredModel$predID <- paste0("Predictor__", seq_len(nrow(UniquePredModel)))

  # Left join to total:
  PredModel <- PredModel %>% left_join(UniquePredModel, by = c("pred", "type"))

  ### Build augmented data ###
  augData <- data

  for (i in seq_len(nrow(UniquePredModel))) {
    if (UniquePredModel$type[i] == "between") {
      augData[[UniquePredModel$predID[i]]] <- ave(
        augData[[UniquePredModel$pred[i]]],
        augData[[idvar]],
        FUN = aveMean
      )
    } else {
      # within: center within cluster
      augData[[UniquePredModel$predID[i]]] <- ave(
        augData[[UniquePredModel$pred[i]]],
        augData[[idvar]],
        FUN = aveCenter
      )
    }
  }

  # Subset to needed columns:
  Vars <- unique(c(vars, PredModel$predID, idvar))
  augData <- augData[, Vars]

  # Remove missings:
  if (na.rm) {
    augData <- na.omit(augData)
  }

  # Check for small clusters:
  tab <- table(augData[[idvar]])
  if (any(tab < 5)) {
    warning(sum(tab < 5), " cluster(s) detected with < 5 observations.")
  }

  ### Dispatch to lmer backend ###
  if (estimator == "lmer") {
    Res <- lmer_mlGGM(
      model = PredModel,
      augData = augData,
      idvar = idvar,
      vars = vars,
      randomeffects = randomeffects,
      verbose = verbose
    )
  } else {
    stop(paste0("Estimator '", estimator, "' not yet implemented for mlGGM."))
  }

  # Store input:
  Res[['input']] <- list(
    vars = vars,
    estimator = estimator,
    randomeffects = randomeffects,
    idvar = idvar
  )

  class(Res) <- "mlGGM"
  return(Res)
}
