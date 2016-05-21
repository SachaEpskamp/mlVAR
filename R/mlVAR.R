aveMean <- function(x){
  mean(x,na.rm=TRUE)
}

aveCenter <- function(x){
  x - mean(x,na.rm=TRUE)
}

aveLag <- function(x, lag=1){
  # Then lag:
  if (lag > length(x)){
    return(rep(NA,length(x)))
  } else {
    return(c(rep(NA,lag),head(x,length(x)-lag))) 
  }
}


mlVAR <- function(
  data, # Data frame
  
  # Variable names:
  vars, # Vector of variables to include in analysis
  idvar, # String indicating the subject id variable name 
  lags = 1, # Vector indicating the lags to include. Defaults to 1
  dayvar, # string indicating the measurement id variable name (if missing, every measurement is set to one day). Used to not model scores over night
  beepvar, # String indicating beep per day (is missing, is added)
  # timevar, # Only used for maxtimeDifference
  
  # Estimation options:
  # orthogonal, # TRUE or FALSE for orthogonal edges. Defaults to nvar < 6
  estimator = c("lmer"), # Add more? estimator = "least-squares" IGNORES multi-level
  # contemporaneous = c("default","shared","unique"), # Shared or unique contamporaneous relationships?
  temporal = c("default", "correlated","orthogonal","fixed"), # Unique = multi-level!
  # betweenSubjects = c("default","GGM","posthoc"), # Should covariances between means be estimated posthoc or as a GGM? Only used when method = "univariate"
  
  # Misc:
  # maxTimeDiff, # If not missing. maximum time difference.
  # LMERcontrol = list(optimizer = "bobyqa"), # "bobyqa"  or "Nelder_Mead"
  # JAGSoptions = list(),
  verbose = TRUE, # Include progress bar?
  compareToLags,
  orthogonal # Used for backward competability
  # JAGSexport = FALSE, # Exports jags files
  # n.chain = 3,
  # n.iter = 10000,
  # estOmega = FALSE
)
{
  if (0 %in% lags & length(lags) > 1){
    stop("0 in 'lags' ignored; contemporaneous relationships are estimated by default.")
    lags <- lags[lags!=0]
  }
  # First check the estimation options:
  # method <- match.arg(method)
  estimator <- match.arg(estimator)
  # contemporaneous <- match.arg(contemporaneous)
  # betweenSubjects <- match.arg(betweenSubjects)
  # temporal <- match.arg(temporal)
  betweenSubjects <- "GGM"
  
  # Nuclear option on JAGS:
  if (estimator == "JAGS"){
    stop("'JAGS' estimator is not implemented yet and expected in a later version of mlVAR.")
  }
  
#   # Only estimate omega when estimator is JAGS
#   if (estOmega & estimator != "JAGS"){
#     stop("Cannot estimate Omega when estimator is not JAGS.")
#   }
  
  # Some dummies for later versions:
  # temporal <- "unique"
  temporal <- match.arg(temporal)
  contemporaneous <- "shared"
  
  if (!missing(orthogonal)){
    temporal <- ifelse(orthogonal,"orthogonal","correlated")
    warning(paste0("'orthogonal' argument is deprecated Setting temporal = '",temporal,"'"))
  }
#   
#   if (contemporaneous == "default"){
#     contemporaneous <- ifelse(estimator == "least-squares", "unique", "shared")
#   }
#   if (betweenSubjects != "default"){
#     if (estimator != "lmer"){
#       warning("'betweenSubjects' argument not used in estimator.")
#     }
#   } else {
#     betweenSubjects <- ifelse(estimator == "least-squares","posthoc","GGM")
#   }

  
  # Unimplemented methods:
#   if (method == "multivariate" & estimator == "lmer"){
#     stop("Multivariate estimation using 'lmer' is not implemented.")
#   }
#   if (contemporaneous == "unique" & estimator %in% c("lmer")){
#     stop(paste0("Unique contemporaenous effects estimation not implemented for estimator = '",estimator,"'"))
#   }
  
  # Obtain variables from mlVARsim0:
  if (is(data,"mlVARsim0")){
    vars <- data$vars
    idvar <- data$idvar
    data <- data$Data
  }
  
  if (is(data,"mlVARsim")){
    vars <- data$vars
    idvar <- data$idvar
    lags <- data$lag
    data <- data$Data
  }
  
  if (missing(lags)){
    lags <- 1
  }
  
  # Set temporal to fixed if lag-0:
  if (all(lags==0)){
    temporal <- "fixed"
  }
  
  # Set orthogonal:
  if (temporal == "default"){
    if (length(vars) > 6){
      
      if (verbose) message("More than 6 nodes, correlations between random effects are set to zero (temporal = 'orthogonal')")
      temporal <- "orthogonal"
      
    } else {
      temporal <- "correlated"
    }    
  }
  
  # CompareToLags:
  if (missing(compareToLags)){
    compareToLags <- lags
  }
  if (length(compareToLags) < length(lags)){
    stop("'compareToLags' must be at least as long as 'lags'")
  }
 
  # Check input:
  stopifnot(!missing(vars))
  stopifnot(!missing(idvar))
  

  # input list (to include in output):
  # input <- list(vars = vars, lags = lags, estimator=estimator,temporal = temporal)
  
  # Add day id if missing:
  if (missing(idvar))
  {
    idvar <- "ID"
    data[[idvar]] <- 1
  } # else input$idvar <- idvar
  
  # Add day var if missing:
  if (missing(dayvar))
  {
    dayvar <- "DAY"
    data[[dayvar]] <- 1
  } # else input$dayvar <- dayvar
  
  # Add beep var if missing:
  if (missing(beepvar))
  {
    beepvar <- "BEEP"
    data[[beepvar]] <- ave(data[[idvar]],data[[idvar]],data[[dayvar]],FUN = seq_along)
  } # else input$beepvar <- beepvar
  
  # Remove NA day or beeps:
  data <- data[!is.na(data[[idvar]]) & !is.na(data[[dayvar]]) & !is.na(data[[beepvar]]), ]
  
  ### Codes from mlVAR
  # Create mlVAR-like predictor data-frame:
  # Within-subjects model:
  PredModel <- expand.grid(
    dep = vars,
    pred = vars,
    lag = compareToLags[compareToLags!=0],
    type =  "within",
    stringsAsFactors = FALSE
  )
  
   # Between-subjects model:
  if (betweenSubjects == "GGM" & estimator != "JAGS"){
    between <- expand.grid(dep=vars,pred=vars,lag=NA,type="between",
                          stringsAsFactors = FALSE)
    
    between <- between[between$dep != between$pred,]
    
    PredModel <- rbind(PredModel,between )
  }
  
  # Unique predictors:
  UniquePredModel <- PredModel[!duplicated(PredModel[,c("pred","lag","type")]),c("pred","lag","type")]
  
  # Ad ID:
  UniquePredModel$predID <- paste0("Predictor__",seq_len(nrow(UniquePredModel)))
  
  # Left join to total:
  PredModel <- PredModel %>% left_join(UniquePredModel, by = c("pred","lag","type"))
  
  # Augment the data
  augData <- data
  
  # Add missing rows for missing beeps
  beepsPerDay <-  eval(substitute(dplyr::summarize_(data %>% group_by_(idvar,dayvar), 
                                             first = ~ min(beepvar,na.rm=TRUE),
                                             last = ~ max(beepvar,na.rm=TRUE)), 
                                  list(beepvar = as.name(beepvar))))
  
  # all beeps:
  allBeeps <- expand.grid(unique(data[[idvar]]),unique(data[[dayvar]]),seq(min(data[[beepvar]],na.rm=TRUE),max(data[[beepvar]],na.rm=TRUE))) 
  names(allBeeps) <- c(idvar,dayvar,beepvar)
  
  # Left join the beeps per day:
  allBeeps <- eval(substitute({
    allBeeps %>% left_join(beepsPerDay, by = c(idvar,dayvar)) %>% 
      group_by_(idvar,dayvar) %>% filter_(~BEEP >= first, ~BEEP <= last)%>%
      arrange_(idvar,dayvar,beepvar)
  },  list(BEEP = as.name(beepvar))))
  
  
  # Enter NA's:
  augData <- augData %>% right_join(allBeeps, by = c(idvar,dayvar,beepvar))
  
  # Add the predictors (when estimatior != JAGS):
  if (estimator != "JAGS"){
    for (i in seq_len(nrow(UniquePredModel))){
      # between: add mean variable:
      
      if (UniquePredModel$type[i] == "between"){
        augData[[UniquePredModel$predID[i]]] <- ave(augData[[UniquePredModel$pred[i]]],augData[[idvar]], FUN = aveMean)
      } else {
        # First include:
        augData[[UniquePredModel$predID[i]]] <-  ave(augData[[UniquePredModel$pred[i]]],augData[[idvar]],augData[[dayvar]], FUN = function(x)aveLag(x,UniquePredModel$lag[i]))
        
        # Then center:
        ### CENTERING ONLY NEEDED WHEN ESTIMATOR != JAGS ###
        if (!estimator %in% c("stan","JAGS")){
          augData[[UniquePredModel$predID[i]]] <- ave(augData[[UniquePredModel$predID[i]]],augData[[idvar]], FUN = aveCenter) 
        } 
      }
      
    }
  }

  # Remove missings from augData:
  Vars <- unique(c(PredModel$dep,PredModel$predID,idvar,beepvar,dayvar))
  augData <- na.omit(augData[,Vars])
  PredModel <- PredModel[is.na(PredModel$lag) | (PredModel$lag %in% lags),]
  
  
  #### RUN THE MODEL ###
  if (estimator == "lmer"){
    Res <- lmer_mlVAR(PredModel,augData,idvar,verbose=verbose, betweenSubjects=betweenSubjects,temporal=temporal)
  # } else if (estimator == "least-squares"){
    # Res <- leastSquares_mlVAR(PredModel,augData,idvar,orthogonal = orthogonal, verbose=verbose, betweenSubjects=betweenSubjects)
#   } else if (estimator == "JAGS"){
#     Res <- JAGS_mlVAR(augData, vars, 
#                       idvar,
#                       lags, 
#                       dayvar, 
#                       beepvar,
#                       temporal = temporal,
#                       orthogonal = orthogonal, verbose=verbose, contemporaneous=contemporaneous,
#                       JAGSexport=JAGSexport,
#                       n.chain = n.chain, n.iter=n.iter,estOmega=estOmega)
  } else  stop(paste0("Estimator '",estimator,"' not yet implemented."))
  
  
  # Add input:
  Res[['input']] <- list(
    vars = vars, 
    lags = lags,
    compareToLags=compareToLags,
    estimator = estimator,
    temporal = temporal
  )
  
  return(Res)
  
}

