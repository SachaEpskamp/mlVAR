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

  curCrit <- critFun(CurResult$lmerResult)
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
        propCrits[i] <- critFun( propResults[[i]]$lmerResult)        
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



# Inner function for nodewise estimation:
NodeWise <- function(
  aData, # Augmented dataset from mlVAR
  response, # Response variable
  idvar, # String indicating the subject id variable name 
  autoLaggedVars, # auto regression lagged vars (always include in stepup and moving window)
  laggedVars, # Vector indicating the lags to include
  include, # Which lagged variables to include?
  includeType = c("stepwise","movingwindow"),
  periodvar, # string indicating the period of measurement.
  treatmentvar, # character vector indicating treatment
  covariates, # character indicating covariates independent of measurement.
  control, # "bobyqa"  or "Nelder_Mead"
  timevar
){
  includeType <- match.arg(includeType)
  
  # If include is missing, include all laggedVars:
  if (missing(include)){
    include <- seq_along(laggedVars)
  }
  
  # Predictor string:
  Pred <- character(0)
  
  # Vector of lagged variables for fixed effects:
  if (includeType == "stepwise"){
    laggedVars_fixedEffects <- paste(c(autoLaggedVars,laggedVars[include]),collapse = " + ")
  } else {
    laggedVars_fixedEffects <- paste(c(autoLaggedVars,laggedVars),collapse = " + ")
  }
  
  # Add to predictor set:
  Pred <- c(Pred, paste("(",laggedVars_fixedEffects,")"))
  
  # If period has more than 1 level:
  if (length(unique(aData[[periodvar]])) > 1)
  {
    # Fixed effect for period:
    Pred <- c(Pred, periodvar)
    
    # Interaction lagged vars & period:
    Pred <- c(Pred, paste("(",laggedVars_fixedEffects,") : ", periodvar))
    
    # If not missing, interaction with treatment:
    if (!missing(treatmentvar))
    {
      # Inteaction period : treatment:
      Pred <- c(Pred, paste(periodvar," : ", treatmentvar))
      
      # 3 way interaction:
      Pred <- c(Pred, paste(periodvar,  ": (",laggedVars_fixedEffects,") : ", treatmentvar))
    }
  } 
  
  # Covariate:
  if (!missing(covariates))
  {
    Pred <- c(Pred, paste0(covariates,"*(",laggedVars_fixedEffects,")"))
  }
  
  # Combine:
  Pred <- paste(Pred, collapse = " + ")
  
  # Create model:
  Rands <- c(autoLaggedVars,laggedVars[include]) 
  Pred <- paste0(Pred," + (",paste(Rands,collapse="+")," |",idvar,")")
  ff <- as.formula(paste(response,"~",Pred))
  
  # RUN LMER:
  Results <- lme4::lmer(ff,data=aData, control = do.call('lmerControl',control),REML=FALSE)
  formula <- ff
 
  # Extract fixed effects:
  feCur <- fixef(Results)
  FixEf <- feCur[names(feCur) %in% c("(Intercept)",Rands),drop=FALSE]
  feSECur <- se.fixef(Results)
  FixEf_SE <- feSECur[names(feSECur) %in% c("(Intercept)",Rands),drop=FALSE]
  
  ### Format results:
  Coef <- as.data.frame(t(FixEf))
  se.Coef <- as.data.frame(t(FixEf_SE))
  pvals <- as.data.frame(2*(1-pnorm(abs(as.matrix(Coef)/as.matrix(se.Coef)))))
  
  Coef <- cbind(dep = response, Coef,stringsAsFactors=FALSE)
  se.Coef <- cbind(dep = response, se.Coef,stringsAsFactors=FALSE)
  pvals <- cbind(dep = response, pvals,stringsAsFactors=FALSE)
  
  # Random effects:
  ranEffects <- lme4::ranef(Results)
  ranPerID <- list()
  for (i in seq_len(nrow(ranEffects[[idvar]]))){
    ranPerID[[i]] <- ranEffects[[idvar]][i,,drop=FALSE]
    ranPerID[[i]] <- cbind(dep = response,ranPerID[[i]],stringsAsFactors=FALSE)
  }
  
  # Variance of random effects:
  Variance <- as.data.frame(t(diag(lme4::VarCorr(Results)[[idvar]])))
  Variance <- cbind(dep =  response, Variance,stringsAsFactors=FALSE)
  
  
  return(list(
    lmerResult = Results,
    formula = formula,
    FixEf = as.data.frame(t(FixEf)),
    FixEf_SE = as.data.frame(t(FixEf_SE)),
    Coef=Coef,
    se.Coef=se.Coef,
    pvals=pvals,
    ranEffects=ranEffects,
    ranPerID=ranPerID,
    Variance=Variance
  ))
}