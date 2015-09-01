

## fixedEffects:
fixedEffects <- function(object){
  if (is(object,"mlVAR_MW")){
    return(object$fixedEffects)
  }
  
  # Nodes:
  Nodes <- unique(object$fixedEffects$dep)
  # Predictors:
  Predictor = names(object$fixedEffects)[-1]
  
  coef <- as.matrix(object$fixedEffects[,-1])
  s.coef <- as.matrix(object$se.fixedEffects[,-1])
  pvals <- as.matrix(object$pvals[,-1])
  
  # Data frame of fixed effects:
  df <- data.frame(
    Response = object$fixedEffects$dep,
    Predictor = Predictor[col(coef)],
    effect = c(coef),
    se = c(coef),
    p = c(coef)
  )
  
  if (any(duplicated(df[c("Response","Predictor")]))){
    message("Duplicate effects found (possibly due to moving window approach), averaging effect, se and p-value.")
    df <- df %>% group_by_("Response","Predictor") %>% summarise_each_(funs(mean(., na.rm=TRUE)), vars = c("effect","se","p"))
  }
  
  return(df)
}


# Random effects:
randomEffects<- function(object){
  if (is(object,"mlVAR_MW")) stop("Cannot estimate random effects with moving window approach")
  
  # Predictors:
  Predictor = colnames(object$randomEffectsVariance)[-1]
  
  # Data frame of fixed effects:
  df <- data.frame(
    Response = object$randomEffectsVariance$dep,
    Predictor = Predictor[col(object$randomEffectsVariance[,-1])],
    variance = unlist(object$randomEffectsVariance[,-1])
  )
  rownames(df) <- NULL
  
  if (any(duplicated(df[c("Response","Predictor")]))){
    message("Duplicate effects found (possibly due to moving window approach), averaging variance.")
    df <- df %>% group_by_("Response","Predictor") %>% summarise_each_(funs(mean(., na.rm=TRUE)), vars = c("variance"))
  }
  
  
  return(df)
}