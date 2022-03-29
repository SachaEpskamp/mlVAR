mlVARsample <- function(
  object,
  nTime = c(25,50,100,200),
  nSample = 100, 
  pMissing = 0, 
  # measures = c("sensitivity", "specificity", "bias", "precision"), 
  nReps = 100,
  nCores = 1,
  ... # mlVAR options
){       
  if (any(nSample > length(object$IDs))){ # what if vector 
    stop("Not possible to have a number of individuals greater than the number of individuals in original mlVAR object")  # if one of nSample, can be a vector 
  }

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
    
    # correlation:
    out$correlation <- cor0(est,real)
    
    out$bias <- bias(est,real)
    
    # Precision (1 - FDR):
    out$precision <- TruePos / (FalsePos + TruePos)
    
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
  
  
  if(length(which(!stationary)) >= 1){
    warning(paste("Non-stationary Beta matrix detected for subject(s): "), paste(which(!stationary), collapse = ", "))
  }
  
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
    nSample = nSample, 
    pMissing = pMissing, 
    # Setup:
    name = "mlVARsim",
    write=FALSE,
    nCores = nCores, # Change this to use more or less computer cores
    reps = nReps, 
    debug=FALSE,
    export = c("cor0","CompareNetworks","bias","dots","object","invTheta","Means","Beta"),
    
    # The simulation code:
    expression = {
      # number of subjects:
      nSubject <- length(Beta)
      # Vary number of families 
      # Select random families 
      random_fam <- sort(sample(nSubject, nSample))
      Beta2 <- Beta[random_fam]
      invTheta2 <-  invTheta[random_fam] 
      Means2 <- Means[random_fam]
      
      nSubject2 <- length(random_fam)
      
      # Input - adjust
      input <- c(object$input,dots)
      input$idvar <- "id"
      input$verbose <- FALSE
      
      # Simulate data:
      simData <- mapply(id = 1:nSubject2,b = Beta2, k = invTheta2, m = Means2, SIMPLIFY = FALSE, 
                        FUN = function(id, b,k,m){
                          data <- graphicalVAR::graphicalVARsim(nTime, b[,,1], k, m)
                          data <- as.data.frame(data)
                          names(data) <- input$vars
                          data$id <- id
                          data
                        })
      simData <- do.call(rbind, simData)
      
      # adjust missingness 
      # missingness 
      total_obs <- nrow(simData)
      rows_na <- sample(nrow(simData), pMissing * total_obs)
      simData[rows_na,input$vars] <- NA 
      
      # Fit model 
      input <- input[!(names(input) %in% c("nSample", "nTime", "pMissing"))]

      Res <- do.call(mlVAR::mlVAR,c(list(data=simData),input))
      
      # # Fixed effects, significant thresholded:

        # True networks:
        true_temporal_fixed_thresholded <- getNet(object, "temporal", nonsig = "hide")
        true_contemporaneous_fixed_thresholded <- getNet(object, "contemporaneous", nonsig = "hide")
        true_between_fixed_thresholded <- getNet(object, "between", nonsig = "hide")
        
        # Estimated networks:
        est_temporal_fixed_thresholded <- getNet(Res, "temporal", nonsig = "hide")
        est_contemporaneous_fixed_thresholded <- getNet(Res, "contemporaneous", nonsig = "hide")
        est_between_fixed_thresholded <- getNet(Res, "between", nonsig = "hide")
        
        # True networks:
        true_temporal_fixed_nothrehsold <- getNet(object, "temporal", nonsig = "show")
        true_contemporaneous_fixed_nothrehsold <- getNet(object, "contemporaneous", nonsig = "show")
        true_between_fixed_nothrehsold <- getNet(object, "between", nonsig = "show")
        
        # Estimated networks:
        est_temporal_fixed_nothrehsold <- getNet(Res, "temporal", nonsig = "show")
        est_contemporaneous_fixed_nothrehsold <- getNet(Res, "contemporaneous", nonsig = "show")
        est_between_fixed_nothrehsold <- getNet(Res, "between", nonsig = "show")
      
     
      # Networks per subject:
      true_temporal_subject <- unlist(lapply(random_fam, function(i){
        getNet(object, "temporal", nonsig = "show", subject = i)
      }))
      true_contemporaneous_subject <-  unlist(lapply(random_fam, function(i){
        net <- getNet(object, "contemporaneous", nonsig = "show", subject = i)
        net[upper.tri(net)]
      }))
      
      # Random effects per subject:
      true_temporal_subject_re <- unlist(lapply(random_fam, function(i){
        getNet(object, "temporal", nonsig = "show", subject = i) - getNet(object, "temporal", nonsig = "show")
      }))
      true_contemporaneous_subject_re <-  unlist(lapply(random_fam, function(i){
        net <- getNet(object, "contemporaneous", nonsig = "show", subject = i)- getNet(object, "contemporaneous", nonsig = "show")
        net[upper.tri(net)]
      }))
      
      # Estimated networks:
      est_temporal_subject <- unlist(lapply(seq_len(nSubject2), function(i){
        getNet(Res, "temporal", nonsig = "show", subject = i)
      }))
      est_contemporaneous_subject <-  unlist(lapply(seq_len(nSubject2), function(i){
        net <- getNet(Res, "contemporaneous", nonsig = "show", subject = i)
        net[upper.tri(net)]
      }))
      
      # Random effects per subject:
      est_temporal_subject_re <- unlist(lapply(seq_len(nSubject2), function(i){
        getNet(Res, "temporal", nonsig = "show", subject = i) - getNet(Res, "temporal", nonsig = "show")
      }))
      est_contemporaneous_subject_re <-  unlist(lapply(seq_len(nSubject2), function(i){
        net <- getNet(Res, "contemporaneous", nonsig = "show", subject = i)- getNet(Res, "contemporaneous", nonsig = "show")
        net[upper.tri(net)]
      }))
      
  
      ### All comparisons ###
      res_temporal_fixed_thresholded <- CompareNetworks(true_temporal_fixed_thresholded, est_temporal_fixed_thresholded, directed = FALSE)
      res_temporal_fixed_thresholded$network <- "temporal_thresholded"
      
      res_contemporaneous_fixed_thresholded <- CompareNetworks(true_contemporaneous_fixed_thresholded, est_contemporaneous_fixed_thresholded, directed = FALSE)
      res_contemporaneous_fixed_thresholded$network <- "contemporaneous_thresholded"
      
      res_between_fixed_thresholded <- CompareNetworks(true_between_fixed_thresholded, est_between_fixed_thresholded, directed = FALSE)
      res_between_fixed_thresholded$network <- "between_thresholded"
      
      # No model selection:
      nomodselect <- function(x){
        x$specificity <- NA
        x$precision <- NA
        x$sensitivity <- NA
        x
      }
      
      res_temporal_fixed_nothrehsold <- CompareNetworks(true_temporal_fixed_nothrehsold, est_temporal_fixed_nothrehsold, directed = FALSE) %>% nomodselect
      res_temporal_fixed_nothrehsold$network <- "temporal"

      
      res_contemporaneous_fixed_nothrehsold <- CompareNetworks(true_contemporaneous_fixed_nothrehsold, est_contemporaneous_fixed_nothrehsold, directed = FALSE)  %>% nomodselect
      res_contemporaneous_fixed_nothrehsold$network <- "contemporaneous"
      
      res_between_fixed_nothrehsold <- CompareNetworks(true_between_fixed_nothrehsold, est_between_fixed_nothrehsold, directed = FALSE)  %>% nomodselect
      res_between_fixed_nothrehsold$network <- "between"
      
      res_temporal_subject <- CompareNetworks(true_temporal_subject, est_temporal_subject, directed = FALSE) %>% nomodselect
      res_temporal_subject$network <- "temporal_subject"
      
      res_contemporaneous_subject <- CompareNetworks(true_contemporaneous_subject, est_contemporaneous_subject, directed = FALSE) %>% nomodselect
      res_contemporaneous_subject$network <- "contemporaneous_subject"
      
      res_temporal_re <- CompareNetworks(true_temporal_subject_re, est_temporal_subject_re, directed = FALSE) %>% nomodselect
      res_temporal_re$network <- "temporal_random_effect"
      
      res_contemporaneous_re <- CompareNetworks(true_contemporaneous_subject_re, est_contemporaneous_subject_re, directed = FALSE) %>% nomodselect
      res_contemporaneous_re$network <- "contemporaneous_random_effect"
      
      # 
      # truePDC_sub <- lapply(random_fam, function(i){getNet(object, "temporal", nonsig = "show", subject = i)})
      # truePCC_sub <- lapply(random_fam, function(i){getNet(object, "contemporaneous", nonsig = "show", subject = i)})
      # estPDC_sub <- lapply(1:nSubject2, function(i){getNet(Res, "temporal", nonsig = "show", subject = i)})
      # estPCC_sub <- lapply(1:nSubject2, function(i){getNet(Res, "contemporaneous", nonsig = "show", subject = i)})
      # 
      # perSubject <- lapply(1:nSubject2, function(i){
      #   
      #   PCCres_sub <- CompareNetworks(truePCC_sub[[i]], estPCC_sub[[i]], directed = FALSE)
      #   PCCres_sub$network <- "contemporaneous"
      #   PDCres <- CompareNetworks(truePDC_sub[[i]], estPDC_sub[[i]], directed = TRUE)
      #   PDCres$network <- "temporal"
      #   
      #   
      #   Results_subj <- rbind(
      #     as.data.frame(PCCres),
      #     as.data.frame(PDCres)
      #   )
      #   
      # })
      # 
      # df_perSubject <- do.call(rbind, perSubject)
      # df_perSubject_temp <- df_perSubject[df_perSubject$network == "temporal",]
      # df_perSubject_contemp <- df_perSubject[df_perSubject$network == "contemporaneous",]
      # 
      # mean_subject_temp <- list()
      # mean_subject_temp$sensitivity <- mean(df_perSubject_temp$sensitivity)
      # mean_subject_temp$specificity <- mean(df_perSubject_temp$specificity)
      # mean_subject_temp$bias <- mean(df_perSubject_temp$bias)
      # mean_subject_temp$precision <- mean(df_perSubject_temp$precision)
      # 
      # mean_subject_contemp <- list()
      # mean_subject_contemp$sensitivity <- mean(df_perSubject_contemp$sensitivity)
      # mean_subject_contemp$specificity <- mean(df_perSubject_contemp$specificity)
      # mean_subject_contemp$bias <- mean(df_perSubject_contemp$bias)
      # mean_subject_contemp$precision <- mean(df_perSubject_contemp$precision)
      # 
      # 
      # # perSubject <- lapply(1:nSubject, function(i){
      # #   
      # #   # Fixed effects, not thresholded:
      # #   # True networks:
      # #   truePDC <- getNet(object, "temporal", nonsig = "show", subject = i)
      # #   truePCC <- getNet(object, "contemporaneous", nonsig = "show", subject = i)
      # #   
      # #   # Estimated networks:
      # #   estPDC <- getNet(Res, "temporal", nonsig = "show", subject = i)
      # #   estPCC <- getNet(Res, "contemporaneous", nonsig = "show", subject = i)
      # #   
      # #   PCCres <- CompareNetworks(truePCC, estPCC, directed = FALSE)
      # #   PCCres$network <- "contemporaneous"
      # #   PDCres <- CompareNetworks(truePDC, estPDC, directed = TRUE)
      # #   PDCres$network <- "temporal"
      # #   
      # #   
      # #   Results_subj <- rbind(
      # #     as.data.frame(PCCres),
      # #     as.data.frame(PDCres)
      # #   )
      # #   
      # # })
      # 
      Results <- dplyr::bind_rows(
        res_temporal_fixed_thresholded,
        res_contemporaneous_fixed_thresholded,
        res_between_fixed_thresholded,
        
        res_temporal_fixed_nothrehsold,
        res_contemporaneous_fixed_nothrehsold,
        res_between_fixed_nothrehsold,
        
        res_temporal_subject,
        res_contemporaneous_subject,
        res_temporal_re,
        res_contemporaneous_re
        
        
      )
      # 
      # 
      # Results <- data.frame(
      #   network = c("temporal_fixed","contemporaneous_fixed","temporal_subject","contemporaneous_subject"),
      #   correlation = c(PDCres$correlation,PCCres$correlation,cor(temporalTrue,temporalEst), 
      #                   cor(contemporaneousTrue,contemporaneousEst)), 
      #   sensitivity = c(PDCres$sensitivity, PCCres$sensitivity, mean_subject_temp$sensitivity, mean_subject_contemp$sensitivity), 
      #   specificity = c(PDCres$specificity, PCCres$specificity, mean_subject_temp$specificity, mean_subject_contemp$specificity), 
      #   bias = c(PDCres$bias, PCCres$bias, mean_subject_temp$bias, mean_subject_contemp$bias), 
      #   precision = c(PDCres$precision, PCCres$precision, mean_subject_temp$precision, mean_subject_contemp$precision)
      # )
      # 
      # usedmeasures <- measures 
      # totalmeasures <- c("sensitivity", "specificity", "bias", "precision")
      # 
      # # if colname not in measure, remove from results 
      # notmeasures <- totalmeasures[which(!totalmeasures %in% usedmeasures)]
      # # remove from results 
      # Results <- Results[ , !(names(Results) %in% notmeasures)]
      return(Results)
    })   
  
  class(Results) <- c("mlVARsample","data.frame")
  return(Results) 
}

# Summary method:
summary.mlVARsample <- function(object, ...){
  sumfun <- function(x){
    if (!all(is.finite(x))) return("") else {
      paste0(round(mean(x,na.rm=TRUE),2)," (SE: ",round(sd(x,na.rm=TRUE)/sqrt(length(x)),2),")")
    }
  }
  
  res <- object %>% group_by(.data[["network"]],.data[["nTime"]],.data[["nSample"]],.data[["pMissing"]]) %>%
    summarise(across(c("sensitivity","specificity","precision","correlation","bias"), ~ sumfun(.x))) %>%
    as.data.frame
  
  return(res)
}
