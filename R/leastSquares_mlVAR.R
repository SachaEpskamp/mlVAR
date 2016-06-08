leastSquares_mlVAR <- 
  function(model,augData,idvar,contemporaneous = "unique", verbose=TRUE,...){
    
    # Global dummies:
    type <- NULL
    ord <- NULL
    pred <- NULL
    
    Outcomes <- unique(model$dep)
    nVar <- length(Outcomes)

    withinMod <- model %>% filter(type == "within") %>%
      group_by(pred, lag) %>% dplyr::summarise(id = unique(predID)) %>%
      mutate(ord = match(pred,Outcomes)) %>% arrange(ord,lag)
    
    predID <- withinMod$id
    predLab <- paste0(withinMod$pred,"_lag",withinMod$lag)
    
    IDs <- unique(augData[[idvar]])
    
    leastSquares <- lapply(seq_along(IDs),function(i){
      
      id <- IDs[i]
      # subject data:
      subjectData <- na.omit(augData[augData[[idvar]]==id,])
      
      # C is outcomes:
      C <- as.matrix(subjectData[,Outcomes])
      
      # Center:
      # C <- t(t(C) - colMeans(C))
      
      # L is predictors:
      L <- as.matrix(subjectData[,predID])
      
      # Append intercept:
      L <- cbind(1,L)

      # Compute beta:
      Beta <- t(C) %*% L %*% solve(t(L) %*% L)

      Theta <- 1/(nrow(C)) * t(C - L %*% t(Beta)) %*% (C - L %*% t(Beta))
      
      Mu <- Beta[,1]
      Beta <- Beta[,-1]
      
      colnames(Beta) <- predID
      
      rownames(Theta) <- colnames(Theta) <- rownames(Beta) <- names(Mu) <- Outcomes
      
      list(
        Mu = Mu,
        Beta = Beta,
        Theta = Theta
      )
    })
    
    # Store results:
    Results <- list()
    
    
    
    # Obtain means:
    mu_subject <- lapply(leastSquares,'[[','Mu')
    mu_table <- do.call(rbind,mu_subject)
    mu_fixed <- colMeans(mu_table)
    mu_cov <- cov(mu_table, use="pairwise.complete.obs")
    mu_SD <- sqrt(diag(mu_cov))
    
    Results[["mu"]] <- modelArray(mean = mu_fixed, SD = mu_SD)
    Results[["Omega_mu"]] <-  modelCov(
      cov = modelArray(mean=mu_cov)
    )
    
    
    # Theta:
    Theta_LeastSquares <- lapply(leastSquares, '[[', 'Theta')
    Theta_fixed_LeastSquares <- apply(do.call(abind,c(Theta_LeastSquares,along=3)),1:2,mean)
    
    
    Results[["Theta"]] <- modelCov(
      cov = modelArray(mean=Theta_fixed_LeastSquares, subject = Theta_LeastSquares)
    )
    
    # Beta:
    withinMod <- model %>% filter(type == "within") %>%
      group_by(pred, lag) %>% dplyr::summarise(id = unique(predID)) %>%
      mutate(ord = match(pred,Outcomes)) %>% ungroup %>% arrange(ord,lag)
    Outcomes <- unique(model$dep)
    
    Beta_subject <- lapply(leastSquares, '[[', 'Beta')
    Beta_fixed <- apply(do.call(abind,c(Beta_subject,along=3)),1:2,mean)
    
    Results[["Beta"]] <- modelArray(mean = BetatoArray(Beta_fixed,withinMod,Outcomes), subject = lapply(Beta_subject,BetatoArray, mod=withinMod, out=Outcomes ))
    
    
    # Combine results:
    Results <- list(
      results = Results,
      output = leastSquares
    )
    
    class(Results) <- "mlVAR"
    
    return(Results)
    
  }
