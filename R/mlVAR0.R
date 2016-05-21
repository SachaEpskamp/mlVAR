mlVAR0 <- function(
  data, # Data frame
  vars, # Vector of variables to include in analysis
  idvar, # String indicating the subject id variable name 
  lags = 1, # Vector indicating the lags to include
  dayvar, # string indicating the measurement id variable name (if missing, every measurement is set to one day)
  beepvar, # String indicating beep per day (is missing, is added)
  periodvar, # string indicating the period of measurement.
  treatmentvar, # character vector indicating treatment
  covariates, # character indicating covariates independent of measurement.
  timevar,
  maxTimeDiff, # If not missing. maximum time difference.
  control = list(optimizer = "bobyqa"), # "bobyqa"  or "Nelder_Mead"
  # windowSize, # Assign to use moving window estimation. If missing, does not use moving window estimation
  verbose = TRUE, # Include progress bar?
  # standardize = c("inSubject","none","general"),
  orthogonal, # by default set to TRUE when nvar > 6 and FALSE otherwise
  estimator = c("lmer","lmmlasso"),
  method = c("default","stepwise","movingWindow"),
  laginteractions = c("none","mains","interactions"), # Include interactions with lag?
  critFun = BIC,
  lambda = 0,
  center = c("inSubject","general","none")
)
{

  
  warning("Function is deprecated. Use mlVAR instead.")
  
  # standardize = c("inSubject","none","general")
  center <- match.arg(center)
  
  if (is(data,"mlVARsim0")){
    
    vars <- data$vars
    idvar <- data$idvar
    data <- data$Data
  }
  

  
  method <- match.arg(method)
  estimator <- match.arg(estimator)
  laginteractions <- match.arg(laginteractions)
  
  if (missing(orthogonal)){
    if (length(vars) > 6 & method != "movingWindow"){
      
      if (verbose) message("More than 6 nodes and method != 'movingWindow', correlations between random effects are set to zero (orthogonal = TRUE)")
      orthogonal <- TRUE
      
    } else {
      orthogonal <- FALSE
    }    
  }

  
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
      x <- plyr::join(dummy,x, by = names(dummy)) 
    }
    
    # Lag variables:
    for (l in lags)
    {
      if (l > nrow(x)) stop("Lag is larger than number of measurements")
      
      lagDF <- x[-(nrow(x)-(1:l)+1),vars]
      lagDF <- lagDF[c(rep(NA,l),seq_len(nrow(lagDF))),]
      names(lagDF) <- paste0("L",l,"_",names(lagDF))
      x <- cbind(x,lagDF)
    }
    
    return(x)
  })
  
  
  ### MULTILEVEL ANALYSIS ###
  # All lagged variables:
  AllLagVars <- c(sapply(lags, function(l) paste0("L",l,"_",vars)))
  
  # Lagged interactions:
  # Lag interactions:
  if (laginteractions != "none" || !missing(maxTimeDiff)){
    if (missing(timevar)) stop("'timevar' needed to include time difference")
    
    call <- substitute(augData %>% group_by_(idvar) %>% mutate(LAGDIFF = as.numeric(c(NA,difftime(x[-1],x[-length(x)],units="hours")))),
                       list(x = as.name(timevar)))
    augData <- eval(call)
    
    
  }
  if (laginteractions != "none"){
    
    
    AllLagVars <- c(AllLagVars,"LAGDIFF")
    
    if (laginteractions == "interactions"){
      
      for (it in seq_along(vars)){
        augData[[paste0("lag_x_",vars[it])]] <- augData[[vars[it]]] * augData$LAGDIFF
      }    
      AllLagVars <- c( AllLagVars ,paste0("lag_x_",vars))
    }
  }
  
  # Center predictors:
  # standardize <- match.arg(standardize)
  Center <- function(x) return(x - mean(x,na.rm=TRUE)) 
  if (center =="inSubject"){
    for(i in unique(data[[idvar]])){
      for (var in AllLagVars){
        augData[augData[[idvar]]==i,var] <- Center(augData[augData[[idvar]]==i,var])
      }
    } 
  } else if (center == "general"){
    augData[,AllLagVars] <- sapply(augData[,AllLagVars] ,Center)
  }
  
  if ( !missing(maxTimeDiff)){
    augData <- augData %>% filter_(~LAGDIFF < maxTimeDiff)
  }
  
  ### RUN NODEWISE ANALYSES ###
  NodeWise_Results <- list()
  Results <- list()
  formulas <- list()
  FixEf <- list()
  FixEf_SE <- list()
  
  # Run models:
  if (verbose){
    pb <- txtProgressBar(min=0,max=length(vars),style=3)    
  }
  
  for (j in seq_along(vars)){
    whichAuto <- grepl(paste0("^L1_",vars[j],"$"), AllLagVars)
    
    if (method == "stepwise"){
      NodeWise_Results[[j]] <- Stepwise(
        aData = augData,
        response = vars[j],
        idvar = idvar, # String indicating the subject id variable name 
        autoLaggedVars = AllLagVars[whichAuto],
        laggedVars = AllLagVars[!whichAuto], # Vector indicating the lags to include
        periodvar = periodvar, # string indicating the period of measurement.
        treatmentvar = treatmentvar, # character vector indicating treatment
        covariates = covariates, # character indicating covariates independent of measurement.
        control = control, # "bobyqa"  or "Nelder_Mead"
        timevar = timevar,
        critFun=critFun,
        progress=verbose,
        estimator = estimator,
        lambda = lambda,
        orthogonal = orthogonal
      )
    } else if (method == "default") {
      NodeWise_Results[[j]] <- NodeWise(
        aData = augData,
        response = vars[j],
        idvar = idvar, # String indicating the subject id variable name 
        autoLaggedVars = AllLagVars[whichAuto],
        laggedVars = AllLagVars[!whichAuto], # Vector indicating the lags to include
        periodvar = periodvar, # string indicating the period of measurement.
        treatmentvar = treatmentvar, # character vector indicating treatment
        covariates = covariates, # character indicating covariates independent of measurement.
        control = control, # "bobyqa"  or "Nelder_Mead"
        timevar = timevar,
        estimator = estimator,
        lambda = lambda,
        orthogonal = orthogonal
      )
    } else if (method=="movingWindow") {
      NodeWise_Results[[j]] <- movingWindow(
        aData = augData,
        response = vars[j],
        idvar = idvar, # String indicating the subject id variable name 
        autoLaggedVars = AllLagVars[whichAuto],
        laggedVars = AllLagVars[!whichAuto], # Vector indicating the lags to include
        periodvar = periodvar, # string indicating the period of measurement.
        treatmentvar = treatmentvar, # character vector indicating treatment
        covariates = covariates, # character indicating covariates independent of measurement.
        control = control, # "bobyqa"  or "Nelder_Mead"
        timevar = timevar,
        estimator = estimator,
        lambda = lambda,
        orthogonal = orthogonal
      )
      
    } else stop("Method not implemented")
    
    
    
    
    if (verbose){
      setTxtProgressBar(pb, j)
    }
  }
  
  if (verbose){
    close(pb)    
  }
  
  
  
  # Fixed effects for each lag:
  
  
  
  #   ### MOVE TO MOVING WINDOW FUNCTION ####
  #   
  #   # Construct the window matrix:
  #   Neffect <- length(AllLagVars) - 1
  #   
  #   
  #   if (!missing(windowSize)){
  #     # Number of possible effects to include is amount of lagged - 1 (autocor is always included):
  #     
  #     
  #     # Window matrix is Neffect * windowSize
  #     # Randomize:
  #     if (Neffect>2) Rand <- sample(seq_len(Neffect)) else Rand <- 1
  #     WindowMatrix <- matrix(Rand[(0:(windowSize-1) + rep(0:(Neffect-1),each=windowSize)) %% Neffect + 1],
  #                            Neffect, windowSize, byrow=TRUE)
  #     
  #   } else WindowMatrix <- matrix(1:Neffect,1)
  #   
  #   cur <- 1
  #   for (j in seq_along(vars)){
  #     
  #     for (w in seq_len(nrow(WindowMatrix))){
  #       PredCur <- Pred
  #       
  #       # Always include the autoregression
  #       whichAuto <- grep(paste0("^L1_",vars[j]), AllLagVars)
  #       RandsCur <- c(AllLagVars[whichAuto],AllLagVars[-j][WindowMatrix[w,]]) 
  #       PredCur <- paste0(PredCur," + (",paste(RandsCur,collapse="+")," |",idvar,")")
  #       
  #       ff <- as.formula(paste(vars[j],"~",PredCur))
  #       Results[[cur]] <- lme4::lmer(ff,data=augData, control = do.call('lmerControl',control),REML=FALSE)
  #       formulas[[cur]] <- ff
  #       
  #       # Extract fixed effects:
  #       feCur <- fixef(Results[[cur]])
  #       FixEf[[cur]] <- feCur[names(feCur) %in% c("(Intercept)",RandsCur)]
  #       
  #       feSECur <- se.fixef(Results[[cur]])
  #       FixEf_SE[[cur]] <- feSECur[names(feSECur) %in% c("(Intercept)",RandsCur)]
  #       
  #       cur <- cur + 1
  #     }
  #     
  #     if (progress){
  #       setTxtProgressBar(pb, j)
  #     }
  #   }
  #   
  #   if (progress){
  #     close(pb)    
  #   }
  #   
  
  if (estimator=="lmmlasso"){
    out <- list(
      fixedEffects = as.data.frame(rbind_all(lapply(NodeWise_Results,"[[","Coef")))
      # se.fixedEffects = as.data.frame(rbind_all(lapply(NodeWise_Results,"[[","se.Coef"))),
      # randomEffects =  ranPerID,
      # randomEffectsVariance = as.data.frame(rbind_all(lapply(NodeWise_Results,"[[","Variance"))),
      # pvals = as.data.frame(rbind_all(lapply(NodeWise_Results,"[[","pvals"))),
      # pseudologlik = logLik,
      # df = df,
      # BIC = BIC,
      # input = input,
      # lmerResults = lapply(NodeWise_Results,"[[","lmerResult"),
      # lmerFormulas = lapply(NodeWise_Results,"[[","formula")
    )
    
    class(out) <- "mlVAR"
    return(out)
  }
  
  # Extract info:
    if (method == "movingWindow"){
      logLik <- NA
      df <- NA
      BIC <- NA
    } else {
  logLik <- sum(unlist(lapply(lapply(NodeWise_Results,"[[","Result"),logLik)))
  df <- sum(unlist(lapply(lapply(lapply(NodeWise_Results,"[[","Result"),logLik),attr,'df')))
  BIC <- sum(unlist(lapply(lapply(NodeWise_Results,"[[","Result"),BIC)))    
  }
  
  
  ranlist <- lapply(NodeWise_Results,"[[","ranPerID")
  ranPerID <- lapply(lapply(seq_along(ranlist[[1]]),function(i)lapply(ranlist,function(x)as.data.frame(x[[i]]))),function(xx)as.data.frame(rbind_all(xx)))

  # Output:
  out <- list(
    fixedEffects = as.data.frame(rbind_all(lapply(NodeWise_Results,"[[","Coef"))),
    se.fixedEffects = as.data.frame(rbind_all(lapply(NodeWise_Results,"[[","se.Coef"))),
    randomEffects =  ranPerID,
    randomEffectsVariance = as.data.frame(rbind_all(lapply(NodeWise_Results,"[[","Variance"))),
    pvals = as.data.frame(rbind_all(lapply(NodeWise_Results,"[[","pvals"))),
#     pseudologlik = logLik,
#     df = df,
#     BIC = BIC,
    input = input,
    lmerResults = lapply(NodeWise_Results,"[[","lmerResult"),
    lmerFormulas = lapply(NodeWise_Results,"[[","formula"))
  
  class(out) <- "mlVAR0"
  
  return(out)  
}

