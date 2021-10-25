# Moving window estimation of fixed effects

movingWindow <- function(
  autoLaggedVars,
  laggedVars,
  maxEffects = 6,
  progress = FALSE,
  estimator = "lmer",
  ... # Arguments sent to mlVAR
){
  # Global dummies:
  . <- NULL
  
  nAuto <- length(autoLaggedVars)
  nLagged <- length(laggedVars)
  nEffects <- nAuto + nLagged
  
  # Create the moving window of autoLaggedVars:
  if (maxEffects > nEffects) stop("Window size is longer than number of effects to estimate. Moving window approach is not needed.")
  
  windowSize <- maxEffects - nAuto
  
 
  Combs <- matrix(,nrow = nLagged, ncol = maxEffects)
  samp <- sample(1:nLagged)
  for (i in 1:nrow(Combs)){
    Combs[i,] <- samp[1 + (i + (seq_len(maxEffects))) %% nLagged]
  }

  # Start estimation:
  Results <- list()
  
  # Setup progress bar:
  if (progress){
    pb <- txtProgressBar(min=0,max=nrow(Combs),style=3)    
  }
  
  # Run main loop:
  for (i in seq_len(nrow(Combs))){
    
    Results[[i]] <- NodeWise(autoLaggedVars=autoLaggedVars,laggedVars=laggedVars,include = Combs[i,],
                             includeType = "movingwindow", ...)
    
    
    # Update progress:
    if (progress){
      setTxtProgressBar(pb, i)
    }
  }
  
  # Close progress:
  if (progress){
    close(pb)    
  }

  # Gather fixed effects:
  if (estimator=="lmmlasso"){
    FullResults <- list(
      Result = lapply(Results,"[[","Result"),
      Coef = Results %>% lapply("[[","Coef") %>% bind_rows() %>% group_by(.data[["dep"]])  %>%
        summarize_each(funs(mean(., na.rm=TRUE)))
    )

  } else {
    
    FullResults <- list(
      Result = lapply(Results,"[[","Result"),
      formula = lapply(Results,"[[","formula"),
      FixEf = Results %>% lapply("[[","FixEf") %>% bind_rows() %>% 
        summarize_each(funs(mean(., na.rm=TRUE))),
      FixEf_SE = Results %>% lapply("[[","FixEf_SE") %>% bind_rows() %>% 
        summarize_each(funs(mean(., na.rm=TRUE))),
      Coef = Results %>% lapply("[[","Coef") %>% bind_rows() %>% group_by(.data[["dep"]])  %>%
        summarize_each(funs(mean(., na.rm=TRUE))),
      se.Coef = Results %>% lapply("[[","se.Coef") %>% bind_rows() %>%group_by(.data[["dep"]])  %>%
        summarize_each(funs(mean(., na.rm=TRUE))),
      pvals = Results %>% lapply("[[","pvals") %>% bind_rows() %>% group_by(.data[["dep"]])  %>%
        summarize_each(funs(mean(., na.rm=TRUE))),
      ranEffects = Results %>% lapply("[[", "ranEffects") %>% lapply(function(x){
        x <- x[[1]]
        x$id <- seq_len(nrow(x))
        x
      }) %>% bind_rows() %>% group_by(id) %>% summarize_each(funs(mean(., na.rm=TRUE))),
      ranPerID = lapply(seq_along(Results[[1]]$ranPerID), function(i){
        Results %>% lapply(function(x)x$ranPerID[[i]]) %>% 
          bind_rows() %>% group_by(.data[["dep"]]) %>% summarize_each(funs(mean(., na.rm=TRUE)))
      }),
      Variance = Results %>% lapply("[[","Variance") %>% bind_rows() %>% group_by(.data[["dep"]]) %>% 
        summarize_each(funs(mean(., na.rm=TRUE)))
    )
  }

  
  return(FullResults)  
}
