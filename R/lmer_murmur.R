BetatoArray <- function(x, mod, out){
  Lags <- sort(unique(mod$lag))
  
  y <- array(NA, c(length(out), length(out), length(Lags)))
  for (l in seq_along(Lags)){
    for (v in seq_along(out)){
      y[,v,l] <- x[,mod$id[mod$pred == out[v] & mod$lag == Lags[l]]]
    }
  }
  dimnames(y) <- list(out,out,Lags)
  y
}

lmer_mlVAR <- 
  function(model,augData,idvar,betweenSubjects = "posthoc", verbose=TRUE,temporal="unique",...){


    
    Outcomes <- unique(model$dep)
    nVar <- length(Outcomes)
    lags <- unique(model$lag[!is.na(model$lag)])
    if (length(lags)==0) lags <- 0
    
    # Output list:
    lmerResults <- list()
    
    if (verbose){
      pb <- txtProgressBar(min = 0, max = length(Outcomes), style = 3)
    }
    # start nodewise loop:
    for (i in seq_along(Outcomes)){
      # submodel:
      subModel <- model %>% filter_(~ dep == Outcomes[i])
      
      # Setup model:
      if (temporal != "fixed"){
        mod <- paste0(
          Outcomes[i],
          " ~ ",
          paste(subModel$predID,collapse="+"),
          " + (",
          paste(subModel$predID[subModel$type == "within"],collapse="+"),
          ifelse(temporal == "orthogonal","||","|"),
          idvar,
          ")"
        )        
      } else {
        mod <- paste0(
          Outcomes[i],
          " ~ ",
          paste(subModel$predID,collapse="+"),
          " + (1 |", 
          idvar,
          ")"
        )
      }
      
      
      # Formula:
      formula <- as.formula(mod)
      
      # Run lmer:
      lmerResults[[i]] <- suppressWarnings(lmer(formula, data = augData,REML=FALSE, ...))
      
      if (verbose){
        setTxtProgressBar(pb, i)
      }
    }
    if (verbose){
      close(pb)
    }
    
    
    ### Collect the results:
    Results <- list()
    
    # Fixed effects:
    # Mu fixed:
    mu_fixed <- unname(sapply(lmerResults,function(x)fixef(x)[1]))
    names(mu_fixed) <- Outcomes
    
    # Standard deviation of means:
    
    mu_SD <- sapply(lmerResults, function(x) attr(lme4::VarCorr(x)[[1]], "stddev")[1]  )
    names(mu_SD) <- Outcomes
    
    # Standard error:
    mu_SE <- unname(sapply(lmerResults,function(x)se.fixef(x)[1]))
    names(mu_SE) <- Outcomes
    
    # Pvalues:
    mu_P <- 2*(1-pnorm(abs(mu_fixed/mu_SE)))
    
    rans <- lapply(lmerResults,function(x)ranef(x)[[1]][,1])
    mu_subject <- lapply(seq_along(rans[[1]]),function(i){
      mu <- mu_fixed + sapply(rans,'[[',i)
      names(mu) <- Outcomes
      mu
    })
    
    ### STORE ###
    Results[['mu']] <- modelArray(mean = mu_fixed,SE = mu_SE,SD = mu_SD,subject = mu_subject)
    
    # Mu covariances (need to be estimated):
    ### First from random effects:
    #     if (betweenSubjects == "posthoc"){
    #       mu_cov <- cov(do.call(rbind,mu_subject))  
    #       mu_prec <- solve(mu_cov)
    #       
    #     } else {
    #       
    ### Second as GGM:
    # Which predictor codes are the variables:
    modSum <- model %>% filter_(~type == "between") %>% group_by_("pred") %>% dplyr::summarize_(id = ~unique(predID))
    IDs <- modSum$id[match(Outcomes, modSum$pred)]
    
    # Construct gamma and D:
    Gamma <- Gamma_SE <- Gamma_P <- matrix(0, length(Outcomes), length(Outcomes))
    D <- diag(1/mu_SD^2)
    for (i in seq_along(Outcomes)){
      Gamma[i,-i] <- fixef(lmerResults[[i]])[IDs[-i]] 
      Gamma_SE[i,-i] <- se.fixef(lmerResults[[i]])[IDs[-i]] 
    }
    rownames(Gamma) <- colnames(Gamma) <- rownames(Gamma_SE) <- colnames(Gamma_SE) <- Outcomes
    Results[["Gamma"]] <- modelArray(mean = Gamma, SE = Gamma_SE)
    
    
    # Inverse estimate:
    inv <- D %*% (diag(length(Outcomes)) - Gamma)
    # Average:
    inv <- (inv + t(inv))/2
    # invert:
    mu_cov <- solve(inv)
    mu_prec <- inv
    # }
    colnames(mu_cov) <- rownames(mu_cov) <- Outcomes
    
    ### Store results:
    Results[["Omega_mu"]] <- modelCov(
      cor = modelArray(mean=cov2cor(mu_cov)),
      cov = modelArray(mean=mu_cov),
      prec = modelArray(mean=mu_prec)
    )
    
    if (!all(lags==0)){
      
      ### Compute Betas ###
      # Obtain the names of predictors in order Outcomes / lag:
      withinMod <- model %>% filter_(~type == "within") %>%
        group_by_("pred", "lag") %>% dplyr::summarise_(id = ~unique(predID)) %>%
        mutate_(ord = ~match(pred,Outcomes)) %>% ungroup %>% arrange_("ord","lag")
      
      predID <- withinMod$id
      # predLab <- paste0(withinMod$pred,"_lag",withinMod$lag)
      
      # Fixed effects:
      Beta_fixed <- do.call(rbind,lapply(lmerResults,function(x)fixef(x)[predID]))
      colnames(Beta_fixed) <- predID
      rownames(Beta_fixed) <- Outcomes
      
      # SE; SD; P; subject
      
      # SE:
      Beta_SE <- do.call(rbind,lapply(lmerResults,function(x)se.fixef(x)[predID]))
      colnames(Beta_SE) <- predID
      rownames(Beta_SE) <- Outcomes
      
      # P:
      Beta_P <- 2*(1-pnorm(abs(Beta_fixed/Beta_SE)))
      colnames(Beta_P) <- predID
      rownames(Beta_P) <- Outcomes
      
      
      #SD:
      if (temporal != "fixed"){
        if (temporal == "correlated"){
          
          Beta_SD <- do.call(rbind,lapply(lmerResults, function(x) attr(lme4::VarCorr(x)[[idvar]],"stddev")[-1]))
          colnames(Beta_SD) <- predID
          rownames(Beta_SD) <- Outcomes
          
        } else {
          #       
          #       VarCorr(lmerResults[[3]])
          #       
          Beta_SD <- do.call(rbind,lapply(lmerResults, function(x){
            df <- as.data.frame(lme4::VarCorr(x))
            df$sdcor[match(predID,df$var1)]
          }))
          colnames(Beta_SD) <- predID
          rownames(Beta_SD) <- Outcomes
          
        }
        #     lme4::VarCorr(lmerResults[[1]])
        #     
        # Subject:
        #     if (!orthogonal){
        #       
        rans <- lapply(lmerResults, function(x)ranef(x)[[idvar]][,predID])  
        Beta_subject <- lapply(seq_len(nrow(rans[[1]])),function(i){
          Beta <- Beta_fixed + do.call(rbind,lapply(rans,function(x)x[i,]))
          colnames(Beta) <- predID
          rownames(Beta) <- Outcomes
          Beta
        })
        
        ### STORE RESULTS ###
        Results[["Beta"]] <- modelArray(
          mean=BetatoArray(Beta_fixed,withinMod,Outcomes),
          SD=BetatoArray(Beta_SD,withinMod,Outcomes),
          SE = BetatoArray(Beta_SE,withinMod,Outcomes), 
          P = BetatoArray(Beta_P,withinMod,Outcomes),
          subject = lapply(Beta_subject,BetatoArray, mod=withinMod, out=Outcomes ))
        
        
      } else {
        Beta_SD <- array(0,c(nVar,nVar,length(unique(model$lag))))
        Beta_subject <- NULL
        
        ### STORE RESULTS ###
        Results[["Beta"]] <- modelArray(
          mean=BetatoArray(Beta_fixed,withinMod,Outcomes),
          SD=Beta_SD,
          SE = BetatoArray(Beta_SE,withinMod,Outcomes), 
          P = BetatoArray(Beta_P,withinMod,Outcomes),
          subject = Beta_subject)
        
      }
    } else {
      Results[["Beta"]] <- NULL 
    }
    
    
    
    #     } else {
    #       browser()
    #     }
    
    # Now reorder everything to array:
    
    
    
    
    
    ### Compute Theta ####
    Theta_obtained <- matrix(NA, nVar, nVar)
    # diag(Theta_obtained)  <- sapply(lmerResults, lme4::sigma)^2
    diag(Theta_obtained)  <- sapply(lmerResults, stats::sigma)^2
    
    # Estimate individual via correlating residuals:
    resids <- lapply(seq_along(Outcomes),
                     function(i){
                       resid <- residuals(lmerResults[[i]])
                       id <- as.numeric(names(resid))
                       df <- data.frame(foo = resid, id = id)
                       left_join(data.frame(id = seq_len(nrow(augData))), df, by = "id")[,2]
                     })
    resid <- as.data.frame(do.call(cbind,resids))
    
    # Compute observed residuals covariances:
    Theta_posthoc <- lapply(unique(augData[[idvar]]),function(id){
      cov(resid[augData[[idvar]] == id,],use = "pairwise.complete.obs")
    })
    
    # abind all and compute means:
    Theta_fixed_posthoc <- apply(do.call(abind,c(Theta_posthoc,along=3)),1:2,mean,na.rm=TRUE)
    
    # Rescale to fit obtained Theta:
    D <- sqrt(diag(diag(Theta_obtained)))
    Theta_fixed <- D %*% cov2cor(Theta_fixed_posthoc) %*% D
    
    Results[["Theta"]] <- modelCov(cov = modelArray(mean = Theta_fixed, subject = Theta_posthoc))
    
    # Using least-squares:
    # IDs of beta:
    betaIDs <- as.numeric(rownames(ranef(lmerResults[[1]])[[idvar]]))
    
    #     Theta_LeastSquares <- lapply(seq_along(betaIDs),function(i){
    #       
    #       id <- betaIDs[i]
    #       # subject data:
    #       subjectData <- na.omit(augData[augData[[idvar]]==id,])
    #       
    #       # C is outcomes:
    #       C <- as.matrix(subjectData[,Outcomes])
    #       
    #       # Center:
    #       C <- t(t(C) - colMeans(C))
    #       
    #       # L is predictors:
    #       L <- as.matrix(subjectData[,predID])
    #       
    #       #       t(C) %*% L %*% solve(t(L) %*% L)
    #       #       Beta_subject[[i]]
    #       
    #       Theta <- 1/(nrow(C)) * t(C - L %*% t(Beta_subject[[i]])) %*% (C - L %*% t(Beta_subject[[i]]))
    #       rownames(Theta) <- colnames(Theta) <- Outcomes
    #       Theta
    #     })
    #     
    #     # abind all and compute means:
    #     Theta_fixed_LeastSquares <- apply(do.call(abind,c(Theta_LeastSquares,along=3)),1:2,mean)
    #     
    #  
    
    #  if (theta == "leastSquares"){
    #    Theta <- 
    #  }
    
    # Goodness of fit
    names(lmerResults) <- Outcomes
    
    fit <- data.frame(
      var = Outcomes,
      aic = sapply(lmerResults,AIC),
      bic = sapply(lmerResults,BIC)
    )
    rownames(fit) <- NULL
    
    
    # Combine results:
    Results <- list(
      results = Results,
      output = lmerResults,
      fit = fit,
      data = augData,
      model = model
      # DIC = samples$BUGSoutput$DIC,
      # DIC_pD = samples$BUGSoutput$pD
    )
    
    class(Results) <- "mlVAR"
    
    return(Results)
  }
