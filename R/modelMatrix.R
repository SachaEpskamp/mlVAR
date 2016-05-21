dimV <- function(x){
  d <- dim(x)
  if (is.null(d)){
    d <- length(x)
  }
  d
}

# Wrapper for samples of fixed and samples of subjects:
parSamples <- function(
  samples, # BUGS samples
  fixed, # Name of fixed
  subject, # Name of subject
  SD  # Name of SD samples
){
  
  # Check if missing:
  if (!missing(fixed)){
    fixedSamples <- samples[[fixed]]
  } else {
    fixedSamples <- NULL
  }
  
  if (!missing(subject)){
    subjectSamples <- samples[[subject]]
  } else {
    subjectSamples <- NULL
  }
  
  if (!missing(SD)){
    SDSamples <- samples[[SD]]
  } else {
    SDSamples <- NULL
  }
  
  # Post-hoc estimate missing samples:
  if (is.null(fixedSamples) && is.null(subjectSamples)){
    stop("'fixed' and 'subject' are both missing")
  }
  
  # If fixed is missing, construct fixed from averaging subjects in each iteration:
  if (is.null(fixedSamples) && !is.null(subjectSamples)){
    dimSubject <- dimV(subjectSamples)
    fixedSamples <- apply(subjectSamples,seq_len(length(dimSubject)-1),mean)
  }
  
  # If SD is missing, compute from subjects:
  if (is.null(SDSamples) && !is.null(subjectSamples)){
    dimSubject <- dimV(subjectSamples)
    SDSamples <- apply(subjectSamples,seq_len(length(dimSubject)-1),sd)
  }
  
  # Subject list:
  if (is.null(subjectSamples)){
    subjectList <- NULL
  } else {
    dimSubject <- dimV(subjectSamples)
    nDim <- length(dimSubject)
    # Structure as list:
    subjectMeans <- apply(subjectSamples, seq_len(nDim)[-1], mean)
    subjectList <- list()
    
    for (p in seq_len(dimSubject[nDim])){
      if (nDim == 3){
        # Vector
        subjectList[[p]] <- subjectMeans[,p,drop=FALSE]
        dim(subjectList[[p]]) <- dimSubject[2]
      } else  if (nDim == 4){
        # Vector
        subjectList[[p]] <- subjectMeans[,,p,drop=FALSE]
        dim(subjectList[[p]]) <- dimSubject[2:3]
      }else  if (nDim == 5){
        # Vector
        subjectList[[p]] <- subjectMeans[,,,p,drop=FALSE]
        dim(subjectList[[p]]) <- dimSubject[2:4]
      } else stop(">3 dimensional arrays not supported")
      
    }
    
  }
  
  aggDims <- 2:length(dimV(fixedSamples))
  mean <- apply(fixedSamples, aggDims, mean)
  lower <- apply(fixedSamples, aggDims, quantile, 0.025)
  upper <- apply(fixedSamples, aggDims, quantile, 0.975)
  if (!is.null(SDSamples)){
    SD <- apply(SDSamples, aggDims, mean)    
  } else {
    SD <- array(NA, dim = dimV(mean))
  }
  
  modelArray(
    mean = mean,
    lower = lower,
    upper = upper,
    SD = SD,
    subject = subjectList
  )
}

# Wrapper function to compute this based on samples of covariance:
covSamples <- function(samples, fixed, subject){
  
  if (missing(subject)){
    subjectSamples <- NULL
  } else {
    subjectSamples <- samples[[subject]]
  }
  
  # If fixed is missing, average from subject:
  if (missing(fixed)){
    if (missing(subject)){
      stop("'subject' and 'fixed' may not both be missing.")
    }
    
    fixedSamples <- apply(subjectSamples,1:3,mean)
  } else {
    fixedSamples <- samples[[fixed]]
  }
  
  # Covariance:
  cov_mean <-  apply(fixedSamples, 2:3, mean)
  cov_lower <-  apply(fixedSamples, 2:3, quantile, 0.025)
  cov_upper <-  apply(fixedSamples, 2:3, quantile, 0.975)
  
  if (!is.null(subjectSamples)){
    cov_subject <- list()
    nP <- dim(subjectSamples)[length( dim(subjectSamples))]
    for (p in 1:nP){
      cov_subject[[p]] <- apply(subjectSamples[,,,p],2:3,mean)
    } 
  } else {
    cov_subject <- NULL
  }
  
  # cor samples:
  corSamples <- fixedSamples
  for (i in 1:dimV(corSamples)[1]){
    corSamples[i,,] <- cov2cor(corSamples[i,,])
  }
  
  # Correlation:
  cor_mean <-  apply(corSamples, 2:3, mean)
  cor_lower <-  apply(corSamples, 2:3, quantile, 0.025)
  cor_upper <-  apply(corSamples, 2:3, quantile, 0.975)
  
  if (!is.null(subjectSamples)){
    cor_subject <- list()
    for (p in 1:nP){
      corSamples <- subjectSamples[,,,p]
      for (i in 1:dimV(corSamples)[1]){
        corSamples[i,,] <- cov2cor(corSamples[i,,])
      }
      cor_subject[[p]] <- apply(corSamples,2:3,mean)
    } 
  } else {
    cor_subject <- NULL
  }
  
  # Precision samples:
  precSamples <- fixedSamples
  for (i in 1:dimV(precSamples)[1]){
    precSamples[i,,] <- solve(precSamples[i,,])
  } 
  
  # Precision:
  prec_mean <-  apply(precSamples, 2:3, mean)
  prec_lower <-  apply(precSamples, 2:3, quantile, 0.025)
  prec_upper <-  apply(precSamples, 2:3, quantile, 0.975)
  
  
  if (!is.null(subjectSamples)){
    prec_subject <- list()
    for (p in 1:nP){
      precSamples <- subjectSamples[,,,p]
      for (i in 1:dimV(precSamples)[1]){
        precSamples[i,,] <- solve(precSamples[i,,])
      }
      prec_subject[[p]] <- apply(precSamples,2:3,mean)
    } 
  } else {
    prec_subject <- NULL
  }
  
  # pcor:  
  pcorSamples <- fixedSamples
  for (i in 1:dimV(pcorSamples)[1]){
    pcorSamples[i,,] <- corpcor::cor2pcor(pcorSamples[i,,])
  } 
  
  pcor_mean <-  apply(pcorSamples, 2:3, mean)
  pcor_lower <-  apply(pcorSamples, 2:3, quantile, 0.025)
  pcor_upper <-  apply(pcorSamples, 2:3, quantile, 0.975)
  
  if (!is.null(subjectSamples)){
    pcor_subject <- list()
    for (p in 1:nP){
      pcorSamples <- subjectSamples[,,,p]
      for (i in 1:dimV(pcorSamples)[1]){
        pcorSamples[i,,] <- corpcor::cor2pcor(pcorSamples[i,,])
      }
      pcor_subject[[p]] <- apply(pcorSamples,2:3,mean)
    } 
  } else {
    pcor_subject <- NULL
  }
  
  
  modelCov(
    cov = modelArray(
      mean = cov_mean,
      lower = cov_lower,
      upper = cov_upper,
      subject = cov_subject
    ),
    cor = modelArray(
      mean = cor_mean,
      lower = cor_lower,
      upper = cor_upper,
      subject = cor_subject
    ),
    prec = modelArray(
      mean = prec_mean,
      lower = prec_lower,
      upper = prec_upper,
      subject = prec_subject
    ),
    pcor = modelArray(
      mean = pcor_mean,
      lower = pcor_lower,
      upper = pcor_upper,
      subject = pcor_subject
    )
  )
}

modelCov <- function(
  cov,
  cor, # Covariance
  prec, # Precision (inverse cov)
  pcor # Partial correlation
){
  
  if (missing(cov)){
    stop("'cov' can not be missing")
  }
  
  if (!missing(cor) && !is(cor,"mlVARarray")){
    stop("'cor' must be missing or an object of class 'mlVARarray'")
  }
  
  if (!missing(prec) && !is(prec,"mlVARarray")){
    stop("'prec' must be missing or an object of class 'mlVARarray'")
  }
  
  if (!missing(pcor) && !is(pcor,"mlVARarray")){
    stop("'pcor' must be missing or an object of class 'mlVARarray'")
  }
  
  cov2corNA <- function(x){
    if (any(is.na(x))){
      x[] <- NA
      return(x)
    } else {
      return(cov2cor(x))
    }
  }
  
  solveNA <- function(x){
    if (any(is.na(x))){
      x[] <- NA
      return(x)
    } else {
      return(solve(x))
    }
  }
  cor2pcorNA <- function(x){
    if (any(is.na(x))){
      x[] <- NA
      return(x)
    } else {
      return(cor2pcor(x))
    }
  }
  
  
  
  # Fill missings (only mean and subject):
  if (missing(cor)){
    cor <- modelArray(
      mean = cov2corNA(cov[['mean']]),
      subject = lapply(cov[['subject']],cov2corNA)
    )
  }
  
  if (missing(prec)){
    prec <- modelArray(
      mean = solveNA(cov[['mean']]),
      subject = lapply(cov[['subject']],solveNA)
    )
  }
  
  if (missing(pcor)){
    pcor <- modelArray(
      mean = cor2pcorNA(cov[['mean']]),
      subject = lapply(cov[['subject']],cor2pcorNA)
    )
  }
  
  
  Results <- list(
    cov=cov,
    cor=cor,
    prec=prec,
    pcor=pcor
  )
  class(Results) <- c("mlVarCov","list")
  return(Results)
}


modelArray <- function(
  mean, # If missing computed from subject
  subject, # LIST containing estimate per subject
  SE, # Can be missing
  P, # If missing computed from mean and SE
  lower, # If missing computed from mean and SE
  upper,
  SD
){
  Results <- list()
  
  if (missing(subject)){
    subject <- NULL
  }
  
  # Either mean or subject must not be missing:
  if (is.null(subject) && missing(mean)){
    stop("Either 'subject' or 'mean' must be used")
  }
  
  # If mean missing but not subject, compute:
  if (missing(mean) && !is.null(subject)){
    N <- length(subject)
    mean <- Reduce("+",subject) / N
  }
  
  # Dimensions:
  dim <- dimV(mean)
  if (is.null(dim)){
    dim <- length(mean)
  }
  
  # If SD missing but not subject, compute:
  if (missing(SD)){
    if (!is.null(subject)){
      N <- length(subject)
      SD <-  Reduce('+',lapply(subject,function(x)(x - mean)^2)) / (N-1)
    } else {
      SD <- array(NA, dim)
    }
  }
  
  # If P is missing but SE is not, compute:
  if (missing(P)){
    if (!missing(SE)){
      P <- 2*(1-pnorm(abs(mean/SE)))
    } else {
      P <- array(NA, dim)
    }
  }
  
  # If SE is missing, NA:
  if (missing(SE)){
    SE <- array(NA, dim)
  }
  
  # If lower is missing but mean and SE not, compute:
  if (missing(lower)){
    if (!missing(SE)){
      lower <- mean - 1.959964 * SE
    } else {
      lower <- array(NA, dim)
    }
  }
  
  # If upper is missing but mean and SE not, compute:
  if (missing(upper)){
    if (!missing(SE)){
      upper <- mean + 1.959964 * SE
    } else {
      upper <- array(NA, dim)
    }
  }
  
  # Subject:
  if (is.null(subject)){
    subject <- NULL
  }
  
  
  # Store mean to list:
  Results[['mean']] <- mean
  Results[['SD']] <- SD
  Results[['lower']] <- lower
  Results[['upper']] <- upper
  Results[['SE']] <- SE
  Results[['P']] <- P
  Results[['subject']] <- subject
  
  class(Results) <- c("mlVARarray","list")
  return(Results)
}


