# Two wrappers:
# Stepwise and MovingWindow
Stepwise <- function(
  autoLaggedVars,
  laggedVars,
  critFun = BIC,
  maxEffects = 6,
  progress=FALSE,
  ... # Args sent to NodeWise
){
  curModel <- character(0)
  
  # Start with model with only auto-regressions:
  CurResult <- NodeWise(
    autoLaggedVars = autoLaggedVars,
    laggedVars =   curModel,
    ...
  )
  
  curCrit <- critFun(CurResult$Result)
  propModels <- list()
  propResults <- list()
  propCrits <- numeric(0)
  
  repeat{
    if (progress){
      message(paste("Current model: ",deparse(paste(CurResult$formula,collapse = ""))))
    }
    for (i in seq_along(laggedVars)){ 
      
      if (!laggedVars[i] %in% curModel){
        propModels[[i]] <- c(curModel,laggedVars[i])
      } else {
        propModels[[i]] <- curModel[curModel != laggedVars[i]]
      }
      
      if (length(propModels[[i]])>maxEffects){
        propResults[[i]] <- list()
        propCrits[i] <- Inf
      } else {
        propResults[[i]] <- NodeWise(
          autoLaggedVars = autoLaggedVars,
          laggedVars =   propModels[[i]],
          ...
        )
        propCrits[i] <- critFun( propResults[[i]]$Result)        
      }
      
    }
    if (!any(propCrits < curCrit)){
      break
    } else {
      best <- which.min(propCrits)
      curModel <- propModels[[i]]
      CurResult <- propResults[[best]]
      curCrit <- propCrits[[best]]
    }
  }
  
  return(CurResult)
}

