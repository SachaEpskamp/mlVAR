Mplus_mlVAR <- 
  function(model,augData,idvar,contemporaneous = "orthogonal", verbose=TRUE,temporal="orthogonal",
           MplusSave = TRUE,
           MplusName = "mlVAR",
           iterations = "(2000)",
           nCores = 1,
           chains = nCores,
           signs,
           ...){
    
    if (temporal == "unique"){
      stop("temporal = 'unique' not supported in Mplus estimation.")
    }
    
    if (contemporaneous == "unique"){
      stop("contemporaneous = 'unique' not supported in Mplus estimation.")
    }
    
    # Abbreviate all to 8 characters:
    model$dep <- abbreviate(model$dep,8)
    model$pred <- abbreviate(model$pred,8)
    names(augData) <- abbreviate(names(augData),8)
    
    Outcomes <- unique(model$dep)
    nVar <- length(Outcomes)
    lags <- unique(model$lag[!is.na(model$lag)])
    if (length(lags)==0) lags <- 0
    
    
    if (contemporaneous != "fixed" && missing(signs)){
      if (verbose){
        message("'signs' argument not specified. Now running fixed effects lmer models to estimate...")
        lmerRes <- mlVAR(augData, Outcomes, idvar, lags=lags, contemporaneous = "fixed", temporal= "fixed")
      }
      signs <- sign(lmerRes$results$Theta$cor$mean)
    }
    
    if (!missing(signs)){
      if (!is.matrix(signs) || nrow(signs) != length(Outcomes) || ncol(signs) != length(Outcomes)){
        stop("'signs' is not a square matrix with a row/column for each variable.")
      }
    }
    
    
    
    if (!identical(lags,1)){
      stop("Currently only lags = 1 supported in estimator = 'Mplus'")
    }
    
    
    
    # Mplus file names:
    if (MplusSave){
      MplusDir <- getwd()
    } else {
      MplusDir <- tempdir()
    }
    
    codesFile <- paste0(MplusDir,"/",MplusName,".inp")
    dataFile <- paste0(MplusDir,"/",MplusName,".dat")
    outFile <- paste0(MplusDir,"/",MplusName,".out")
    samplesFile <- paste0(MplusDir,"/",MplusName,".txt")
    mlVARmodelFile <- paste0(MplusDir,"/",MplusName,".mlVARmodel")
    signsFile <- paste0(MplusDir,"/",MplusName,".mlVARsigns")
    
    # Only use columns of variables and ID variable in dataset:
    augData <- augData[,c(Outcomes, idvar)]
    
    # Store Mplus datafile:
    suppressWarnings(MplusAutomation::prepareMplusData(augData, dataFile))
    
    # Store signs file:
    if (contemporaneous != "fixed"){
      write.table(signs, signsFile, col.names=FALSE, row.names=FALSE)
    }
    
    ### PREPARE THE MPLUS FILE ###
    mplusCodes <- MPLUS_TEMPLATE
    
    # Replace datafile:
    mplusCodes <- gsub("@@DATAFILE@@",dataFile,mplusCodes)
    
    # Replace variables:
    mplusCodes <- gsub("@@VARIABLENAMES@@",paste(names(augData), collapse = " "),mplusCodes)
    mplusCodes <- gsub("@@USEVARIABLES@@",paste(c(Outcomes,idvar), collapse = " "),mplusCodes)
    
    # Replace cluster:
    mplusCodes <- gsub("@@CLUSTER@@",idvar,mplusCodes)
    
    # Replace lagged:
    lags <- model %>% dplyr::group_by(.data[["pred"]]) %>% summarize(lag=max(.data[['lag']]))
    mplusCodes <- gsub("@@LAGGED@@",paste0(lags$pred,"(",lags$lag,")",collapse=" "),mplusCodes)
    
    # Replace iterations:
    mplusCodes <- gsub("@@ITERATIONS@@",iterations,mplusCodes)
    
    # Replace chains:
    mplusCodes <- gsub("@@CHAINS@@",chains,mplusCodes)
    
    # Replace cores:
    mplusCodes <- gsub("@@CORES@@",nCores,mplusCodes)
    
    # Replace bayesoutput:
    mplusCodes <- gsub("@@BAYESSAMPLES@@",samplesFile,mplusCodes)
    
    ### CREATE THE MODEL ###
    # Temporal model:
    temporalModel <- 
      data.frame(
        dep = model$dep,
        pred = model$pred,
        lag = model$lag,
        model = "temporal",
        random = temporal,stringsAsFactors=FALSE
      )
    
    # Contemporaneous model:
    # Variances:
    contemporaneousModel <- data.frame(
      dep = Outcomes,
      pred = Outcomes,
      lag = 0,
      model = "cont var",
      random = contemporaneous,stringsAsFactors=FALSE
    )
    
    # Correlations:
    contCorrelation <- as.data.frame(t(combn(Outcomes,2)),stringsAsFactors=FALSE)
    names(contCorrelation) <- c("dep","pred")
    contCorrelation$lag <- 0
    contCorrelation$model <- "cont cov"
    contCorrelation$random <- contemporaneous
    
    # Combine:
    contemporaneousModel <- rbind(contemporaneousModel,contCorrelation,stringsAsFactors=FALSE)
    rm(contCorrelation)
    
    # Full model:
    fullModel <- rbind(temporalModel, contemporaneousModel)
    rm(temporalModel, contemporaneousModel)
    
    # Add parameter number:
    fullModel$par <- paste0("par",seq_len(nrow(fullModel)))
    
    # within model string:
    withinMod <- ''
    
    
    # Write temporal model:
    if (any(fullModel$model == "temporal")){
      submodel <- fullModel[fullModel$model == "temporal",]
      
      # Random or not?
      randHeader <- ifelse(submodel$random != "fixed",
                           paste0(submodel$par," | "),
                           ""
      )
      
      # Write:
      withinMod <- paste0(withinMod,"\n",
                          paste0(randHeader,submodel$dep," ON ",submodel$pred,"&",submodel$lag,";",collapse = "\n"))
    }
    
    # Test: if contemporaneous == "fixed", let Mplus take care of it:
    if (contemporaneous != "fixed"){
      
      # # Write contemporaneous variances:
      if (any(fullModel$model == "cont var")){
        submodel <- fullModel[fullModel$model == "cont var",]
        
        # Random or not?
        randHeader <- ifelse(submodel$random != "fixed",
                             paste0(submodel$par," | "),
                             ""
        )
        
        # Write:
        withinMod <- paste0(withinMod,"\n",
                            paste0(randHeader,submodel$dep,";",collapse = "\n"))
      }
      
      # Write contemporaneous correlations:
      if (any(fullModel$model == "cont var")){
        submodel <- fullModel[fullModel$model == "cont cov",]
        
        # For every cov, create a dummy factors with equal factor loadings:
        dummyFactors <- paste0("Fac",seq_len(nrow(submodel)))
        # dummySeq <- seq_along(dummyFactors)
        
        # # Random or not?
        # randHeader <- ifelse(submodel$random != "fixed",
        #                      paste0(submodel$par," | ",dummyFactors),
        #                      dummyFactors
        # )
        
        
        # # Write the first factor loadings:
        # withinMod <- paste0(withinMod,"\n",
        # paste0(randHeader," BY ",submodel$dep,"* (",dummySeq,");",collapse = "\n"))        
        # 
        # # Write the second factor loadings:
        # withinMod <- paste0(withinMod,"\n",
        # paste0(randHeader," BY ",submodel$pred,"* (",dummySeq,");",collapse = "\n"))        
        
        # Trial: unscaled factor loadings
        # # Write the first factor loadings:
        # withinMod <- paste0(withinMod,"\n",
        #                     paste0(randHeader," BY ",submodel$dep,"*;",collapse = "\n"))        
        # 
        # # Write the second factor loadings:
        # withinMod <- paste0(withinMod,"\n",
        #                     paste0(randHeader," BY ",submodel$pred,"*;",collapse = "\n"))  
        # 
        
        # 1 or -1 for second?
        signDummy <- ifelse(signs[lower.tri(signs,diag=FALSE)]>0,1,-1)
        
        # Trial: Factor loadings fixed to 1
        withinMod <- paste0(withinMod,"\n",
                            paste0(dummyFactors," BY ",submodel$dep,"@1 ",submodel$pred,"@",signDummy,";",collapse = "\n"))
        
        
        # Random or not?
        randHeader <- ifelse(submodel$random != "fixed",
                             paste0(submodel$par," | "),
                             ""
                             
        )
        
        # Free variances:
        withinMod <- paste0(withinMod,"\n",
                            paste0(randHeader,dummyFactors,";",collapse = "\n"))
        
        
        # # Constrain their variances:
        # withinMod <- paste0(withinMod,"\n",
        #                     paste0(dummyFactors,"@1;",collapse = "\n"))
        
        # Make them orthogonal:
        # if (contemporaneous == "fixed"){
        # factorsComb <- t(combn(dummyFactors,2))
        # withinMod <- paste0(withinMod,"\n",
        #                     paste0(factorsComb[,1]," WITH ",factorsComb[,2],"@0;",collapse = "\n"))
        
        withinMod <- paste0(withinMod,"\n GDUMMY@1;\n",
                            paste0("GDUMMY BY ",dummyFactors,"@0;",collapse="\n"))
        # }
        
        
        
      }
    }
    
    # Write within model:
    mplusCodes <- gsub("@@WITHINMODEL@@",withinMod,mplusCodes)
    
    
    ## Between subjects model ##
    betModel <- paste0(paste(Outcomes[1]," - ",Outcomes[length(Outcomes)], " WITH ", Outcomes[1]," - ",Outcomes[length(Outcomes)]),";")
    
    if (!all(fullModel$random == "fixed")){
      #   mplusCodes <- gsub("@@BETWEENMODEL@@","",mplusCodes)
      #   mplusCodes <- gsub("%BETWEEN%","",mplusCodes)
      # } else {
      #   
      # Add random variances:
      submodel <- fullModel[fullModel$random != "fixed",]
      betModel <- paste0(betModel,"\n",
                         paste0(submodel$par," WITH ",submodel$par,";",collapse = "\n"))
      
      # Add correlated:
      if (any(fullModel$random == "correlated")){
        submodel <- fullModel[fullModel$random == "correlated",]
        corPars <- t(combn(submodel$par,2))
        betModel <- paste0(betModel,"\n",
                           paste0(corPars[,1]," WITH ",corPars[,2],";",collapse = "\n"))
      }
      
      
      
    }
    
    # Write between model:
    mplusCodes <- gsub("@@BETWEENMODEL@@",betModel,mplusCodes)
    
    ### WRITEMODEL FILE TO DISK ###
    write(mplusCodes,codesFile)
    write.table(fullModel, mlVARmodelFile,row.names = FALSE)
    
    ### RUN THE MPLUS MODEL FILE ###
    if (verbose){
      message("Running Mplus...")
    }
    MplusAutomation::runModels(codesFile, showOutput=verbose)
    
    
    ### Collect the results:
    Results <- importMplus(outFile)
    
    return(Results)
  }
