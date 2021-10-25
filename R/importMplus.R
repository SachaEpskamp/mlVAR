# A function that takes the output of Mplus and reads it to an mlVAR object. Can be used to make changes to codes.
importMplus <- function(outfile){

  # Read output:
  foo <- capture.output(output <- suppressMessages(suppressWarnings(readModels(outfile,what = "all"))))

  # Read results into R:
  # saveInfo <- suppressWarnings(getSavedata_Fileinfo(outfile))
  saveInfo <- output$savedata_info
  
  # Read samples into R:
  bayesSamples <- output$bparameters$valid_draw
    #suppressWarnings(getSavedata_Bparams(outfile) )
  
  if (!is(bayesSamples,"mcmc")){
    bayesSamples <- do.call(rbind,bayesSamples)
  }
  
  # Read data into R:
  Data <- read.table(gsub("\\.out",".dat",outfile),na.strings = ".")
  names(Data) <- unlist(strsplit(output$input$variable$names,split = " "))
  
  # Read mlVAR model into R:
  mlvarFile <- gsub("\\.out",".mlVARmodel",outfile)
  if (!file.exists(mlvarFile)){
    stop(paste0("mlVAR output not found, expected a file ",mlvarFile))
  }
  mlVARinfo <- read.table(mlvarFile,header=TRUE,stringsAsFactors=FALSE)
  

  
  # Names of headers in samples:
  sampleNames <- saveInfo$bayesVarNames
  sampleNames <- gsub("^Parameter\\.\\d*\\_","",sampleNames)
  
  ### Collect the results ####
  Results <- list()
  
  # Outcomes:
  Outcomes <- unique(mlVARinfo$dep)
  nOut <- length(Outcomes)
  
  # Predictors:
  Predictors <- unique(mlVARinfo$pred)
  nPred <- length(Predictors)
  
  # Number of lags:
  Lags <- sort(unique(mlVARinfo$lag))
  Lags <- Lags[Lags>0]
  nLag <- length(Lags)
  
  # Number of samples:
  nSample <- nrow(bayesSamples)
  
  if (!identical(Outcomes,Predictors)){
    stop("Mplus estimation with different predictors than outcomes is not yet supported.")
  }
  
  ### BETWEEN SUBJECTS ###
  # Always random
  # Gather the variance--covariance matrices:
  Omega_mu_cov_samples <- array(,c(nOut, nOut, nSample))
  for (i in 1:nOut){
    for (j in i:nOut){
      # Name of column:
      if (i == j){
        colName <- paste0("%BETWEEN%:.",toupper(Outcomes[i]))
        colID <- which(sampleNames == colName)
        Omega_mu_cov_samples[i,j,] <- bayesSamples[,colID]
      } else {
        colName <- paste0("%BETWEEN%:.",toupper(Outcomes[j]),".WITH.",toupper(Outcomes[i]))
        colID <- which(sampleNames == colName)
        Omega_mu_cov_samples[i,j,] <- Omega_mu_cov_samples[j,i,] <- bayesSamples[,colID]
      }
    }
  }
  # Store to output:
  Results[["Omega_mu"]] <- covSamples_Mplus(Omega_mu_cov_samples)
  
  ### TEMPORAL NETWORK ###
  tempModel <- unique(mlVARinfo[mlVARinfo$model == "temporal","random"])
  if (length(tempModel) > 1){
    stop("Temporal model is a mixture of fixed/orthogonal/correlated")
  }
  
  subModel <- mlVARinfo %>% dplyr::filter(.data[['model']] == "temporal")
  # Gather the parameter samples:
  BetaFixed <- array(,c(nOut, nPred, nLag,nSample))
  BetaSD <- array(0,c(nOut, nPred, nLag,nSample))
  
  if (tempModel != "fixed"){
    for (i in 1:nOut){
      for (j in 1:nPred){
        for (l in 1:nLag){
          # Name of parameter:
          parName <- mlVARinfo$par[mlVARinfo$dep == Outcomes[i] & mlVARinfo$pred == Predictors[j] & mlVARinfo$lag == Lags[l]]
          
          if (length(parName) == 1){
            # Fixed effect:
            sampleParName <- paste0("%BETWEEN%:.MEAN.",toupper(parName))
            sampleParID <- which(sampleNames == sampleParName)
            BetaFixed[i,j,l,] <- bayesSamples[,sampleParID]
            
            # Standard deviation:
            sampleParName <- paste0("%BETWEEN%:.",toupper(parName))
            sampleParID <- which(sampleNames == sampleParName)
            BetaSD[i,j,l,] <- bayesSamples[,sampleParID]
          }
        }
      }
    }
    
  } else {
    ### FIXED EFFECTS ### 
    
    for (i in 1:nOut){
      for (j in 1:nPred){
        for (l in 1:nLag){
          # Name of column:
          colName <- paste0("%WITHIN%:.",toupper(Outcomes[1]),".ON.",toupper(Predictors[j]),"&",Lags[l])

            # Fixed effect:
            sampleParID <- which(sampleNames == colName)
            BetaFixed[i,j,l,] <- bayesSamples[,sampleParID]
        }
      }
    }
    
  }
  
  # Store to output:
  Results[["Beta"]] <- parSamples_Mplus(BetaFixed,BetaSD)
  
  ### CONTEMPORANEOUS EFFECTS ###
  conModel <- unique(mlVARinfo[grepl("cont",mlVARinfo$model),"random"])
  
  # Factor names:
  submodel <- mlVARinfo[mlVARinfo$model == "cont cov",]
  
  # For every cov, create a dummy factors with equal factor loadings:
  dummyFactors <- paste0("Fac",seq_len(nrow(submodel)))
  
  if (length(conModel) > 1){
    stop("Contemporaneous model is a mixture of fixed/orthogonal/correlated")
  }
  if (conModel != "fixed"){

    # First collect the matrix as is:
    Theta_samples <- array(,c(nOut, nOut, nSample))
    
    for (i in 1:nOut){
      for (j in i:nOut){
        # Name of column:
        # Name of parameter:
        parName <- mlVARinfo$par[mlVARinfo$dep == Outcomes[i] & mlVARinfo$pred == Predictors[j] & mlVARinfo$lag == 0]
        
        if (length(parName) == 1){
          # Fixed effect:
          sampleParName <- paste0("%BETWEEN%:.MEAN.",toupper(parName))              
          sampleParID <- which(sampleNames == sampleParName)
          Theta_samples[i,j,] <- Theta_samples[j,i,] <- exp(bayesSamples[,sampleParID])
        }
      }
    }
    
    # Read signs into R:
    signsFile <- gsub("\\.out",".mlVARsigns",outfile)
    if (!file.exists(signsFile)){
      stop(paste0("mlVAR signs file not found, expected a file ",mlvarFile))
    }
    signs <- as.matrix(read.table(signsFile,header=FALSE,stringsAsFactors=FALSE))
    
    # Now add sum of off-diagonals to each diagonal!
    for (t in 1:dim(Theta_samples)[3]){
      diag(Theta_samples[,,t]) <- colSums(Theta_samples[,,t])
      
      # Multiply with signs:
      Theta_samples[,,t] <- Theta_samples[,,t] * signs
    }
    
  } else {

    # First collect the matrix as is:
    Theta_samples <- array(,c(nOut, nOut, nSample))
    for (i in 1:nOut){
      for (j in i:nOut){
        # Name of column:
        if (i == j){
          colName <- paste0("%WITHIN%:.",toupper(Outcomes[i]))
          colID <- which(sampleNames == colName)
          Theta_samples[i,j,] <- bayesSamples[,colID]
        } else {
          # Which dummy?
          colName <- paste0("%WITHIN%:.",toupper(Outcomes[j]),".WITH.",toupper(Outcomes[i]))
          colID <- which(sampleNames == colName)
          Theta_samples[i,j,] <- Theta_samples[j,i,] <- bayesSamples[,colID]
        }
      }
    }
  }
  
  # # Now add sum of off-diagonals to each diagonal!
  # for (t in 1:dim(Theta_samples)[3]){
  #   diag(Theta_samples[,,t]) <- colSums(Theta_samples[,,t])
  #   
  #   # Multiply with signs:
  #   Theta_samples[,,t] <- Theta_samples[,,t] * signs
  # }

  # Store to output:
  Results[["Theta"]] <- covSamples_Mplus(Theta_samples)
  
  # Obtain DIC:
  if (is.null(output$summaries$DIC)){
    warning("No DIC obtained. Model might not have converged.")
    DIC <- NA
    pD <- NA
  } else {
    DIC <- output$summaries$DIC
    pD <- output$summaries$pD
  }

  # Combine results:
  Results <- list(
    results = Results,
    output = output,
    fit = NULL,
    data = Data,
    model = mlVARinfo,
    DIC = DIC,
    pD = pD,
    input = list(
      vars = Outcomes,
      lags = Lags,
      temporal = tempModel,
      contemporaneous = conModel,
      estimator = sprintf("Mplus (version %s)",output$summaries$Mplus.version)
    )
  )
  
  class(Results) <- "mlVAR"
  
  return(Results)
}
