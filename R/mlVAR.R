mlVAR <- function(
  data, # Data frame
  vars, # Vector of variables to include in analysis
  idvar, # String indicating the subject id variable name 
  lags = 1, # Vector indicating the lags to include
  dayvar, # string indicating the measurement id variable name (if missing, every measurement is set to one day)
  beepvar, # String indicating beep per day (is missing, is added)
  periodvar, # string indicating the period of measurement.
  treatmentvar, # character vector indicating treatment
  covariates, # character indicating covariates independent of measurement.
  control = list(optimizer = "bobyqa"), # "bobyqa"  or "Nelder_Mead"
  windowSize, # Assign to use moving window estimation. If missing, does not use moving window estimation
  progress = TRUE, # Include progress bar?
  timevar,
  laginteractions = c("none","mains","interactions") # Include interactions with lag?
)
{
  laginteractions <- match.arg(laginteractions)
  
  # Check input:
  stopifnot(!missing(vars))
  stopifnot(!missing(idvar))
  
  # inout list (to include in output):
  input <- list(vars = vars, lags = lags)
  
  # Add day id if missing:
  if (missing(idvar))
  {
    idvar <- "ID"
    data[[idvar]] <- 1
  } else input$idvar <- idvar
  
  # Add day var if missing:
  if (missing(dayvar))
  {
    dayvar <- "DAY"
    data[[dayvar]] <- 1
  } else input$dayvar <- dayvar
  
  # Add beep var if missing:
  if (missing(beepvar))
  {
    beepvar <- "BEEP"
    data[[beepvar]] <- ave(data[[idvar]],data[[idvar]],data[[dayvar]],FUN = seq_along)
  } else input$beepvar <- beepvar
  
  # Add period var if missing:
  if (missing(periodvar))
  {
    periodvar <- "PERIOD"
    data[[periodvar]] <- 1
  } else input$periodvar <- periodvar
  
  if (!missing(treatmentvar)) input$treatmentvar <- treatmentvar
  if (!missing(covariates)) input$covariates <- covariates
  
  # Remove NA period, day or beeps:
  data <- data[!is.na(data[[idvar]]) & !is.na(data[[periodvar]])  &  !is.na(data[[dayvar]]) & !is.na(data[[beepvar]]), ]
  
  # Create augmented lagged data:
  augData <- plyr::ddply(data, c(idvar,dayvar,periodvar), function(x){
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
  # All lagged variables:
  AllLagVars <- c(sapply(lags, function(l) paste0("L",l,"_",vars)))
  


  # Lag interactions:
  if (laginteractions != "none"){
    if (missing(timevar)) stop("'timevar' needed to include laginteractions")
    
      call <- substitute(augData %>% group_by_(idvar) %>% mutate(LAGDIFF = c(NA,diff(x))),
                         list(x = as.name(timevar)))
      augData <- eval(call)
      
      
      AllLagVars <- c(AllLagVars,"LAGDIFF")
      
      
      if (laginteractions == "interactions"){
        
        for (it in seq_along(vars)){
          augData[[paste0("lag_x_",vars[it])]] <- augData[[vars[it]]] * augData$LAGDIFF
        }    
        
        AllLagVars <- c( AllLagVars ,paste0("lag_x_",vars))

        
      }


  }

  # Vector of lagged variables for fixed effects:
  lagVars <- paste(AllLagVars,collapse = " + ")
  
  # Add to predictor set:
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
    
    #     # Random effects:
    #     ### -1????
    #     Pred <- c(Pred, paste0("( factor(",periodvar,") - 1 + ",lagVarsRand," |",idvar,")"))
  } #else {
  # Random effects:
  #     ### -1????
  #     Pred <- c(Pred, paste0("(",lagVarsRand," |",idvar,")"))
  #   }
  
  
  # Covariate:
  if (!missing(covariates))
  {
    Pred <- c(Pred, paste0(covariates,"*(",lagVars,")"))
  }
  
  # Combine:
  Pred <- paste(Pred, collapse = " + ")
  
  Results <- list() # List to store results in
  formulas <- list()
  FixEf <- list()
  FixEf_SE <- list()
  
  # Run models:
  if (progress){
    pb <- txtProgressBar(min=0,max=length(vars),style=3)    
  }
  
  # Construct the window matrix:
  Neffect <- length(AllLagVars) - 1
  
  
  if (!missing(windowSize)){
    # Number of possible effects to include is amount of lagged - 1 (autocor is always included):
    
    
    # Window matrix is Neffect * windowSize
    # Randomize:
    if (Neffect>2) Rand <- sample(seq_len(Neffect)) else Rand <- 1
    WindowMatrix <- matrix(Rand[(0:(windowSize-1) + rep(0:(Neffect-1),each=windowSize)) %% Neffect + 1],
                           Neffect, windowSize, byrow=TRUE)
    
  } else WindowMatrix <- matrix(1:Neffect,1)
  
  cur <- 1
  for (j in seq_along(vars)){
    
    for (w in seq_len(nrow(WindowMatrix))){
      PredCur <- Pred
      
      # Always include the autoregression
      whichAuto <- grep(paste0("^L1_",vars[j]), AllLagVars)
      RandsCur <- c(AllLagVars[whichAuto],AllLagVars[-j][WindowMatrix[w,]]) 
      PredCur <- paste0(PredCur," + (",paste(RandsCur,collapse="+")," |",idvar,")")
      
      ff <- as.formula(paste(vars[j],"~",PredCur))
      Results[[cur]] <- lme4::lmer(ff,data=augData, control = do.call('lmerControl',control),REML=FALSE)
      formulas[[cur]] <- ff
      
      # Extract fixed effects:
      feCur <- fixef(Results[[cur]])
      FixEf[[cur]] <- feCur[names(feCur) %in% c("(Intercept)",RandsCur)]
      
      feSECur <- se.fixef(Results[[cur]])
      FixEf_SE[[cur]] <- feSECur[names(feSECur) %in% c("(Intercept)",RandsCur)]
      
      cur <- cur + 1
    }
    
    if (progress){
      setTxtProgressBar(pb, j)
    }
  }
  
  if (progress){
    close(pb)    
  }
  
  
  # Extract info:
  if (!missing(windowSize)){
    logLik <- NA
    df <- NA
    BIC <- NA
  } else {
    logLik <- sum(unlist(lapply(Results,logLik)))
    df <- sum(unlist(lapply(lapply(Results,logLik),attr,'df')))
    BIC <- sum(unlist(lapply(Results,BIC)))    
  }
  
  
  Coef <- rbind_all(lapply(FixEf,function(x)as.data.frame(t(x))))
  se.Coef <- rbind_all(lapply(FixEf_SE,function(x)as.data.frame(t(x))))
  
  pvals <- as.data.frame(2*(1-pnorm(abs(as.matrix(Coef)/as.matrix(se.Coef)))))
  
  
  Coef <- cbind(dep = rep(vars,each=nrow(WindowMatrix)), Coef)
  se.Coef <- cbind(dep = rep(vars,each=nrow(WindowMatrix)), se.Coef)
  pvals <- cbind(dep = rep(vars,each=nrow(WindowMatrix)), pvals)
  
  # Random effects:
  ranEffects <- lapply(Results, lme4::ranef)
  ranPerID <- list()
  for (i in seq_len(nrow(ranEffects[[1]][[idvar]]))){
    ranPerID[[i]] <- rbind_all(lapply(ranEffects,function(x)x[[idvar]][i,]))
    ranPerID[[i]] <- cbind(dep = rep(vars,each=Neffect),ranPerID[[i]])
  }
  # names(ranPerID) <- rownames( ranEffects[[1]][[idvar]] )
  
  # Variance of random effects:
  Variance <- rbind_all(lapply(Results,function(x)as.data.frame(t(diag(lme4::VarCorr(x)[[idvar]])))))
  Variance <- cbind(dep =  rep(vars,each=Neffect), Variance)
  
  # Output:
  out <- list(
    fixedEffects = Coef,
    se.fixedEffects = se.Coef,
    randomEffects =  ranPerID,
    randomEffectsVariance = Variance,
    pvals = pvals,
    pseudologlik = logLik,
    df = df,
    BIC = BIC,
    input = input,
    lmerResults = Results,
    lmerFormulas = formulas)
  #     results = Results)
  
  class(out) <- "mlVAR"
  
  return(out)  
}
# 
# fixef.mlvar <- function(object) object$fixedEffects
# se.fixef.mlvar <- function(object) object$se.fixedEffects
# 
# summary.mlvar <- function(x) x[c('coef','se.coef','pvals')]
# coef.mlvar <- function(x) x$coef
