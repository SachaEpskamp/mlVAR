asterix <- function(x){
  y <- rep("",length(x))
  y[x<0.05] <- "*"
  y[x<0.01] <- "**"
  y[x<0.001] <- "***"
  y
}

## fixedEffects:
fixedEffects <- function(object,digits=5){
  # global dummies:
  . <- NULL
  effect <- NULL
  
  if (is(object,"mlVAR_MW")){
    return(object$fixedEffects)
  }
  
  if (class(object)!="mlVAR0"){
    stop("Only works for mlVAR0 objects.")
  }
  
  warning("Function is deprecated and will be removed soon.")

  # Nodes:
  Nodes <- unique(object$fixedEffects$dep)
  # Predictors:
  Predictor = names(object$fixedEffects)[-1]
  
  coef <- as.matrix(object$fixedEffects[,-1])
  if (!is.null(object$se.fixedEffects)){
    s.coef <- as.matrix(object$se.fixedEffects[,-1])
  } else {
    s.coef <- NA
  }
  if (!is.null(object$pval)){
    pvals <- as.matrix(object$pvals[,-1])
  } else {
    pvals <- NA
  }


  # Data frame of fixed effects:
  df <- data.frame(
    Response = object$fixedEffects$dep,
    Predictor = Predictor[col(coef)],
    effect = round(c(coef),digits),
    se = round(c(s.coef),digits),
    p = round(c(pvals),digits),
    sig = asterix(c(pvals)),
    stringsAsFactors=FALSE
  )
  
  if (any(duplicated(df[c("Response","Predictor")]))){
    message("Duplicate effects found (possibly due to moving window approach), averaging effect, se and p-value.")
    df <- df %>% group_by_("Response","Predictor") %>% summarise_each_(funs(mean(., na.rm=TRUE)), vars = c("effect","se","p"))
  }

  # Remove NA rows:
  df <- df %>% filter(!is.na(effect))
  
  return(df)
}


# Random effects:
randomEffects<- function(object, digits=5){
  if (is(object,"mlVAR_MW")) stop("Cannot estimate random effects with moving window approach")

  # global dummies:
  . <- NULL
  variance <- NULL
  
  
  if (class(object)!="mlVAR0"){
    stop("Only works for mlVAR0 objects.")
  }
  
  warning("Function is deprecated and will be removed soon.")
  
  # Predictors:
  Predictor = colnames(object$randomEffectsVariance)[-1]
  
  # Data frame of fixed effects:
  df <- data.frame(
    Response = object$randomEffectsVariance$dep,
    Predictor = Predictor[col(object$randomEffectsVariance[,-1])],
    variance = round(unlist(object$randomEffectsVariance[,-1]),digits)
  )
  rownames(df) <- NULL
  
  if (any(duplicated(df[c("Response","Predictor")]))){
    message("Duplicate effects found (possibly due to moving window approach), averaging variance.")
    df <- df %>% group_by_("Response","Predictor") %>% summarise_each_(funs(mean(., na.rm=TRUE)), vars = c("variance"))
  }
  
  # Remove NA rows:
  df <- df %>% filter(!is.na(variance))
  
  return(df)
  return(df)
}
