mlVARsample <- function(
  object,
  nTime = c(25,50,100,200),
  nReps = 100,
  nCores = 1,
  ... # mlVAR options
){
  if (!identical(object$input$lags,1)){
    stop("Only supported for lags = 1")
  }
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args)==0){
    args <- 1
  }
  cor0 <- function(x,y,...){
    if (sum(!is.na(x)) < 2 || sum(!is.na(y)) < 2 || sd(x,na.rm=TRUE)==0 | sd(y,na.rm=TRUE) == 0){
      return(0)
    } else {
      return(cor(x,y,...))
    }
  }
  # Read the parSim function:
  bias <- function(x,y) mean(abs(x-y),na.rm=TRUE)
  CompareNetworks <- function(true,est, directed = TRUE){
    if (is.matrix(true) & is.matrix(est)){
      if (directed){
        real <- c(true)
        est <- c(est)
      } else {
        real <- true[upper.tri(true,diag=FALSE)]
        est <- est[upper.tri(est,diag=FALSE)]
      }        
    } else {
      real <- true
    }
    
    # True positives:
    TruePos <- sum(est != 0 &  real != 0)
    
    # False pos:
    FalsePos <- sum(est != 0 & real == 0)
    
    # True Neg:
    TrueNeg <- sum(est == 0 & real == 0)
    
    # False Neg:
    FalseNeg <- sum(est == 0 & real != 0)
    
    out <- list()
    
    ### Sensitivity:
    out$sensitivity <- TruePos / (TruePos + FalseNeg)
    
    # Specificity:
    out$specificity <- TrueNeg / (TrueNeg + FalsePos)
    
    # Correlation:
    out$Correlation <- cor0(est,real)
    
    out$Bias <- bias(est,real)
    
    return(out)
  }
  
  # Dots:
  dots <- list(...)
  
  
  # Beta matrices:
  Beta <- object$results$Beta$subject
  
  # Check if all are stationary:
  stationary <- sapply(Beta,function(b){
    ev <- eigen(b[,,1])$values
    all(Re(ev)^2 + Im(ev)^2 < 1)
  })
  
  # Inverse thetas:
  invTheta <- lapply(object$results$Theta$cov$subject,solve)
  
  # Check for outling variance in contemoraneous cov matrices:
  contVar <- sapply(invTheta,function(x)sum(diag(solve(x))))
  # Remove all 10* higher than median:
  goodVar <- contVar < 10*median(contVar)
  
  
  # Intercepts:
  Means <- object$results$mu$subject
  
  
  if (any(!goodVar | !stationary)){
    warning(paste0("Model for subject(s) ",paste(which(!goodVar | !stationary),collapse=" "),
                   " are not proper removed from simulation"))
    
    keep <- goodVar | stationary
    invTheta <- invTheta[keep]
    Beta <- Beta[keep]
    Means <- Means[keep]
  }
  
  ### SIMULATION ###
  Results <- parSim(
    # timepoints conditions:
    nTime = nTime,
    
    # Setup:
    name = "mlVARsim",
    write=FALSE,
    nCores = nCores, # Change this to use more or less computer cores
    reps = nReps, # Number of repetitions per condition
    debug=FALSE,
    
    export = c("cor0","CompareNetworks","bias","dots","object","invTheta","Means","Beta"),
    
    # The simulation code:
    expression = {
      
      # number of subjects:
      nSubject <- length(Beta)
 
      
      # Input
      input <- c(object$input,dots)
      input$idvar <- "id"
      
      # Simulate data:
      simData <- mapply(id = 1:nSubject,b = Beta, k = invTheta, m = Means, SIMPLIFY = FALSE, 
                        FUN = function(id, b,k,m){
                          data <- graphicalVAR::graphicalVARsim(nTime, b[,,1], k, m)
                          data <- as.data.frame(data)
                          names(data) <- input$vars
                          data$id <- id
                          data
                        })
      simData <- do.call(rbind, simData)
      
      
      # Fit model:
      Res <- do.call(mlVAR::mlVAR,c(list(data=simData),input))
      
      # # Fixed effects, significant thresholded:
      # # True networks:
      # truePDC <- getNet(object, "temporal")
      # truePCC <- getNet(object, "contemporaneous")
      # 
      # # Estimated networks:
      # estPDC <- getNet(Res, "temporal")
      # estPCC <- getNet(Res, "contemporaneous")
      # 
      # PCCres <- CompareNetworks(truePCC, estPCC, directed = FALSE)
      # PCCres$network <- "contemporaneous_sig"
      # PDCres <- CompareNetworks(truePDC, estPDC, directed = TRUE)
      # PDCres$network <- "temporal_sig"
      # 
      # 
      # Results_sig <- rbind(
      #   as.data.frame(PCCres),
      #   as.data.frame(PDCres)
      # )
      
      # Fixed effects, not thresholded:
      # True networks:
      truePDC <- getNet(object, "temporal", nonsig = "show")
      truePCC <- getNet(object, "contemporaneous", nonsig = "show")
      
      # Estimated networks:
      estPDC <- getNet(Res, "temporal", nonsig = "show")
      estPCC <- getNet(Res, "contemporaneous", nonsig = "show")
      
      PCCres <- CompareNetworks(truePCC, estPCC, directed = FALSE)
      PCCres$network <- "contemporaneous"
      PDCres <- CompareNetworks(truePDC, estPDC, directed = TRUE)
      PDCres$network <- "temporal"
      
      # Per subject:
      temporalTrue <- unlist(lapply(1:nSubject, function(i){
        getNet(object, "temporal", nonsig = "show", subject = i)
      }))
      contemporaneousTrue <-  unlist(lapply(1:nSubject, function(i){
        net <- getNet(object, "contemporaneous", nonsig = "show", subject = i)
        net[upper.tri(net)]
      }))
      
      temporalEst <- unlist(lapply(1:nSubject, function(i){
        getNet(Res, "temporal", nonsig = "show", subject = i)
      }))
      contemporaneousEst <-  unlist(lapply(1:nSubject, function(i){
        net <- getNet(Res, "contemporaneous", nonsig = "show", subject = i)
        net[upper.tri(net)]
      }))
      
      # perSubject <- lapply(1:nSubject, function(i){
      #   
      #   # Fixed effects, not thresholded:
      #   # True networks:
      #   truePDC <- getNet(object, "temporal", nonsig = "show", subject = i)
      #   truePCC <- getNet(object, "contemporaneous", nonsig = "show", subject = i)
      #   
      #   # Estimated networks:
      #   estPDC <- getNet(Res, "temporal", nonsig = "show", subject = i)
      #   estPCC <- getNet(Res, "contemporaneous", nonsig = "show", subject = i)
      #   
      #   PCCres <- CompareNetworks(truePCC, estPCC, directed = FALSE)
      #   PCCres$network <- "contemporaneous"
      #   PDCres <- CompareNetworks(truePDC, estPDC, directed = TRUE)
      #   PDCres$network <- "temporal"
      #   
      #   
      #   Results_subj <- rbind(
      #     as.data.frame(PCCres),
      #     as.data.frame(PDCres)
      #   )
      #   
      # })
      
      Results <- data.frame(
        network = c("temporal_fixed","contemporaneous_fixed","temporal_subject","contemporaneous_subject"),
        correlation = c(PDCres$Correlation,PCCres$Correlation,cor(temporalTrue,temporalEst), 
                        cor(contemporaneousTrue,contemporaneousEst))
      )
      
      return(Results)
    })   
 
  class(Results) <- "mlVARsample"
  return(Results) 
}
