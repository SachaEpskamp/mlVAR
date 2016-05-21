

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
  estimator = c("lmer","lmmlasso"),
  lambda = 0,
  timevar,
  orthogonal = FALSE
){
  includeType <- match.arg(includeType)
  estimator <- match.arg(estimator)
  # pdMat <- match.arg(pdMat)

  # If include is missing, include all laggedVars:
  if (missing(include)){
    include <- seq_along(laggedVars)
  }
  
  if (estimator == "lmer"){
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
    if (orthogonal){
      Pred <- gsub("\\|","||",Pred)      
    }
    ff <- as.formula(paste(response,"~",Pred))
    
    # RUN LMER:
    Results <- suppressWarnings(lme4::lmer(ff,data=aData, control = do.call('lmerControl',control),REML=FALSE))
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
    if(orthogonal){
      ranEf <- as.data.frame(lme4::VarCorr(Results))
      Variance <- as.data.frame(t(ranEf$vcov[grepl(idvar,ranEf$grp)]))      
      names(Variance) <- ranEf$var1[grepl(idvar,ranEf$grp)]
    } else {
      Variance <- as.data.frame(t(diag(lme4::VarCorr(Results)[[idvar]])))
    }

    Variance <- cbind(dep =  response, Variance,stringsAsFactors=FALSE)
    
    
    return(list(
      Result = Results,
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
    
  } else if (estimator == "lmmlasso"){
    stop("'lmmlasso' not supported")
#     if (orthogonal){
#       pdMat <- "pdDiag"
#     } else {
#       pdMat <- "pdSym"
#     }
#     # Removing missing values:
#     aData <- aData[rowSums(is.na(aData))==0,]
# 
#     # Random effects matrix:
#     Rand <- cbind(intercept=1,aData[c(autoLaggedVars,laggedVars[include])])
#     
#     # Fixed effects matrix:
#     Pred <- character(0)
#     
#     # Vector of lagged variables for fixed effects:
#     if (includeType == "stepwise"){
#       laggedVars_fixedEffects <- c(autoLaggedVars,laggedVars[include])
#     } else {
#       laggedVars_fixedEffects <-c(autoLaggedVars,laggedVars)
#     }
#     
#     # Add to predictor set:
#     Pred <- c(Pred,laggedVars_fixedEffects)
#     
#     # If period has more than 1 level:
#     if (length(unique(aData[[periodvar]])) > 1)
#     {
#       # Fixed effect for period:
#       Pred <- c(Pred, periodvar)
#       
#       # Add interaction columns:
#       
#       # Interaction lagged vars & period:
#       int <- aData[laggedVars_fixedEffects] * aData[[periodvar]]
#       names(int) <- paste0(names(int),"_x_",periodvar)
#       aData <- cbind(aData,int)
#       Pred <- c(Pred,  names(int))
#       
#       # If not missing, interaction with treatment:
#       if (!missing(treatmentvar))
#       {
#         Pred <- c(Pred, treatmentvar)
#         
#         # Interaction lagged vars & period:
#         int <- aData[treatmentvar] * aData[[periodvar]]
#         names(int) <- paste0(treatmentvar,"_x_",periodvar)
#         aData <- cbind(aData,int)
#         Pred <- c(Pred,  names(int))
#         
#         # 3 way interaction:
#         int <- aData[treatmentvar] * aData[[periodvar]] * aData[[treatmentvar]]
#         names(int) <- paste0(treatmentvar,"_x_",periodvar,"_x_",treatmentvar)
#         aData <- cbind(aData,int)
#         Pred <- c(Pred,  names(int))
#       }
#     } 
#     
#     # Covariate:
#     if (!missing(covariates))
#     {
#       Pred <- c(Pred, covariates)
#       for (CV in seq_along(covariates)){
#         int <- aData[laggedVars_fixedEffects] * aData[[covariates[CV]]]
#         names(int) <- paste0(laggedVars_fixedEffects,"_x_",covariates[[CV]])
#         aData <- cbind(aData,int)
#         Pred <- c(Pred,  names(int))
#       }
#     }
#     
#     # Fixed effects:
#     Fixed <- cbind(intercept=1,aData[Pred])
#     
#     # Run lmmlasso:
#     noPen <- c(1,which(names(Fixed)%in%autoLaggedVars))
#     Res <- lmmlasso(x=as.matrix(Fixed), y=aData[[response]], z = as.matrix(Rand), grp = aData[[idvar]], lambda = lambda, nonpen = noPen, pdMat=pdMat)
# 
#     Coef <- cbind(dep = response, as.data.frame(t(coef(Res)[-1])))
#     return(list(
#       Result = Res,
# #       formula = formula,
# #       FixEf = as.data.frame(t(FixEf)),
# #       FixEf_SE = as.data.frame(t(FixEf_SE)),
#       Coef=Coef
# #       se.Coef=se.Coef,
# #       pvals=pvals,
# #       ranEffects=ranEffects,
# #       ranPerID=ranPerID,
# #       Variance=Variance
#     ))
    
  }

 
}
