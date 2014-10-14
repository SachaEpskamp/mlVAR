psyvar <- function(
  data, # Data frame
  vars, # Vector of variables to include in analysis
  idvar, # String indicating the subject id variable name
  dayvar, # string indicating the measurement id variable name (if missing, every measurement is set to one day)
  beepvar, # String indicating beep per day (is missing, is added)
  periodvar, # string indicating the period of measurement.
  lags = 1, # Vector indicating the lags to include
  treatmentvar, # character vector indicating treatment
  covariates, # character indicating covariates independent of measurement.
  maxiter = 800
)
{
  # Load packages:
  stopifnot(require("plyr"))
  stopifnot(require("lme4"))
  stopifnot(require("arm"))
  
  # Check input:
  stopifnot(!missing(vars))
  stopifnot(!missing(idvar))
  
  # Add day id if missing:
  if (missing(idvar))
  {
    idvar <- "ID"
    data[[idvar]] <- 1
  }
  
  # Add day var if missing:
  if (missing(dayvar))
  {
    dayvar <- "DAY"
    data[[dayvar]] <- 1
  }
  
  # Add beep var if missing:
  if (missing(beepvar))
  {
    beepvar <- "BEEP"
    data[[beepvar]] <- ave(data[[idvar]],data[[idvar]],data[[dayvar]],FUN = seq_along)
  }
  
  # Add period var if missing:
  if (missing(periodvar))
  {
    periodvar <- "PERIOD"
    data[[periodvar]] <- 1
  }
  
  # Remove NA period, day or beeps:
  data <- data[!is.na(data[[idvar]]) & !is.na(data[[periodvar]])  &  !is.na(data[[dayvar]]) & !is.na(data[[beepvar]]), ]
  
  # Create augmented lagged data:
  augData <- ddply(data, c(idvar,dayvar,periodvar), function(x){
    # Check for duplicate beeps:
    if (any(duplicated(x[[beepvar]])))
    {
      stop("Duplicated beepnumbers found. ")
    }
    # Order by beep:
    x <- x[order(x[[beepvar]]),]
    
    
    # Augment missing:
    beepseq <- seq(1, max(x[[beepvar]]))
    if (!identical(x[[beepvar]], beepseq))
    {
      dummy <- data.frame(ID = unique(x[[idvar]]), PERIOD= unique(x[[periodvar]]), DAY = unique(x[[dayvar]]), BEEP = beepseq)
      names(dummy) <- c(idvar,periodvar,dayvar,beepvar)
      x <- join(dummy,x, by = names(dummy)) 
    }
    
    # Lag variables:
    for (l in lags)
    {
      if (l > nrow(x)) stop("Lag is higher than amount of measurements")
      
      lagDF <- x[-(nrow(x)-(1:l)+1),vars]
      lagDF <- lagDF[c(rep(NA,l),seq_len(nrow(lagDF))),]
      names(lagDF) <- paste0("L",l,"_",names(lagDF))
      x <- cbind(x,lagDF)
    }
    
    return(x)
  })
  
  
  ### MULTILEVEL ANALYSIS ###
  # Predictor string:
  Pred <- character(0)
  
  # Fixed effects for each lag:
  lagVars <- paste(sapply(lags, function(l) paste0("L",l,"_",vars, collapse = " + ")),collapse = " + ")
  Pred <- c(Pred, paste("(",lagVars,")"))
  
  # If period has more than 1 level:
  if (length(unique(augData[[periodvar]])) > 1)
  {
    # Fixed effect for period:
    Pred <- c(Pred, periodvar)
    
    # Interaction lagged vars & period:
    Pred <- c(Pred, paste("(",lagVars,") : ", periodvar))
    
    # If not missing, interaction with treatment:
    if (!missing(treatmentvar))
    {
      # Inteaction period : treatment:
      Pred <- c(Pred, paste(periodvar," : ", treatmentvar))
      
      # 3 way interaction:
      Pred <- c(Pred, paste(periodvar,  ": (",lagVars,") : ", treatmentvar))
    }
    
    # Random effects:
    ### -1????
    Pred <- c(Pred, paste0("( factor(",periodvar,") - 1 + ",lagVars," |",idvar,")"))
  } else {
    # Random effects:
    ### -1????
    Pred <- c(Pred, paste0("(",lagVars," |",idvar,")"))
  }
  
  
  # Covariate:
  if (!missing(covariates))
  {
    Pred <- c(Pred, paste0(covariates,"*(",lagVars,")"))
  }
  
  
  
  # Combine:
  Pred <- paste(Pred, collapse = " + ")
  
  Results <- list() # List to store results in
  
  # Run models:
  pb <- txtProgressBar(min=0,max=length(vars),style=3)
  for (j in seq_along(vars)){
    ff <- as.formula(paste(vars[j],"~",Pred))
    Results[[j]] <- lmer(ff,data=augData, control = lmerControl(optCtrl=list(maxIter=maxiter)),REML=FALSE)
    setTxtProgressBar(pb, j)
  }
  close(pb)
  
  # Extract info:
  BIC <- sum(unlist(lapply(Results,BIC)))
  Coef <- do.call(rbind,lapply(Results,fixef))
  se.Coef <- do.call(rbind,lapply(Results,se.fixef))
  colnames(Coef) <- names(fixef(Results[[1]]))
  colnames(se.Coef) <- names(fixef(Results[[1]]))
  rownames(Coef) <- columns
  rownames(se.Coef) <- columns
  
  # Random effects:
  ranEffects <- lapply(Results, ranef)
  ranPerID <- list()
  for (i in seq_len(nrow(ranEffects[[1]][[idvar]]))){
    ranPerID[[i]] <- do.call(rbind, lapply(ranEffects,function(x)x[[idvar]][i,]))  
    rownames(ranPerID[[i]]) <- vars
  }
  names(ranPerID) <- rownames( ranEffects[[1]][[idvar]] )
  
  # Variance of random effects:
  Variance <- do.call(rbind,lapply(Results,function(x)diag(VarCorr(x)[[idvar]])))
  rownames(Variance) <- vars
  
  # Output:
  out <- list(
    fixedEffects = Coef,
    se.fixedEffects = se.Coef,
    randomEffects =  ranPerID,
    randomEffectsVariance = Variance,
    pvals = 2*(1-pnorm(abs(Coef/se.Coef))),
    BIC = BIC)
  #     results = Results)
  
  class(out) <- "psyvar"
  
  return(out)  
}
# 
# fixef.mlvar <- function(object) object$fixedEffects
# se.fixef.mlvar <- function(object) object$se.fixedEffects
# 
# summary.mlvar <- function(x) x[c('coef','se.coef','pvals')]
# coef.mlvar <- function(x) x$coef
