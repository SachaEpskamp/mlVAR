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

forcePositive <- function(x){
  x <- (x + t(x))/2
  if (any(eigen(x)$values < 0)){
    return(x - (diag(nrow(x)) * min(eigen(x)$values)-0.001))    
  } else {
    return(x)
  }
  
}

lmer_mlVAR <- 
  function(model,augData,idvar,contemporaneous = "orthogonal", verbose=TRUE,temporal="orthogonal",
           nCores = 1, AR = FALSE, ...){
    
    
    
    Outcomes <- unique(model$dep)
    nVar <- length(Outcomes)
    lags <- unique(model$lag[!is.na(model$lag)])
    if (length(lags)==0) lags <- 0
    
    # Output list:
    lmerResults <- list()
    
    if (verbose){
      message("Estimating temporal and between-subjects effects")
      # pb <- txtProgressBar(min = 0, max = length(Outcomes), style = 3)
    }
    # start nodewise loop:
    
    
    if (nCores > 1){
      # Make clusters:
      nClust <- nCores - 1
      cl <- makePSOCKcluster(nClust)
      
      # Export to cluster:
      clusterExport(cl, c("Outcomes", "model", "temporal", "contemporaneous", "idvar", "augData"), envir = environment())
      
      # Run loop:
      lmerResults <-  parLapply(cl, seq_along(Outcomes), function(i){
        # submodel:
        subModel <-  dplyr::filter_(model, ~ dep == Outcomes[i])
        if (AR){
          subModel <- subModel %>% filter_(~dep == pred | type == "between")
        }
        
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
        return(suppressWarnings(lmer(formula, data = augData,REML=FALSE, ...)))
      })

      
    } else {
      if (verbose){
        pb <- txtProgressBar(min = 0, max = length(Outcomes), style = 3)
      }
      
      for (i in seq_along(Outcomes)){
        # submodel:
        subModel <- model %>% filter_(~ dep == Outcomes[i])
        
        # Remove cross-lagged if AR = TRUE:
        if (AR){
          subModel <- subModel %>% filter_(~dep == pred | type == "between")
        }
        
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
    #       mu_prec <- corpcor::pseudoinverse(mu_cov)
    #       
    #     } else {
    #       
    ### Second as GGM:
    # Which predictor codes are the variables:
    modSum <- model %>% filter_(~type == "between") %>% group_by_("pred") %>% dplyr::summarize_(id = ~unique(predID))
    IDs <- modSum$id[match(Outcomes, modSum$pred)]
    
    # Construct gamma and D:
    Gamma_Omega_mu <- Gamma_Omega_mu_SE <- Gamma_Omega_mu_P <- matrix(0, length(Outcomes), length(Outcomes))
    
    for (i in seq_along(Outcomes)){
      Gamma_Omega_mu[i,-i] <- fixef(lmerResults[[i]])[IDs[-i]] 
      Gamma_Omega_mu_SE[i,-i] <- se.fixef(lmerResults[[i]])[IDs[-i]] 
    }
    rownames(Gamma_Omega_mu) <- colnames(Gamma_Omega_mu) <- rownames(Gamma_Omega_mu_SE) <- colnames(Gamma_Omega_mu_SE) <- Outcomes
    Results[["Gamma_Omega_mu"]] <- modelArray(mean = Gamma_Omega_mu, SE = Gamma_Omega_mu_SE)
    
    
    # Inverse estimate:
    if (any(mu_SD ==0)){
      warning("Zero SD found in mean of following variables: ",paste(names(mu_SD[which(mu_SD==0)]), collapse = ", ")," - Between-subject effects could not be estimated")
      
      mu_cov <- mu_prec <- inv <- matrix(NA, length(Outcomes), length(Outcomes))
      
      colnames(mu_cov) <- rownames(mu_cov) <- Outcomes
      
      ### Store results:
      Results[["Omega_mu"]] <- modelCov(
        cor = modelArray(mean=mu_cov),
        cov = modelArray(mean=mu_cov),
        prec = modelArray(mean=mu_prec)
      )
      
    } else {
      D <- diag(1/mu_SD^2)
      inv <- D %*% (diag(length(Outcomes)) - Gamma_Omega_mu)
      # Average:
      inv <- (inv + t(inv))/2
      # Force positive/definite:
      inv <- forcePositive(inv) # - diag(nrow(inv)) * min(min(eigen(inv)$values),0)
      
      # invert:
      mu_cov <- corpcor::pseudoinverse(inv)
      mu_prec <- inv      
      
      colnames(mu_cov) <- rownames(mu_cov) <- Outcomes
      
      ### Store results:
      Results[["Omega_mu"]] <- modelCov(
        cor = modelArray(mean=cov2cor(mu_cov)),
        cov = modelArray(mean=mu_cov),
        prec = modelArray(mean=mu_prec)
      )
    }

    # }
    
    
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
      
      # AR:
      if (AR){
        Beta_fixed[is.na(Beta_fixed)] <- 0
        Beta_SE[is.na(Beta_SE)] <- 0
        Beta_P[is.na(Beta_P)] <- 0
      }
      
      #SD:
      if (temporal != "fixed"){
        if (temporal == "correlated"){
          
          Beta_SD <- do.call(rbind,lapply(lmerResults, function(x) attr(lme4::VarCorr(x)[[idvar]],"stddev")[-1]))
          if (AR){
            Beta_SD <- diag(c(unlist(Beta_SD)))
          }
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
          if (AR){
            Beta_SD <- diag(c(unlist(Beta_SD)))
          }
          colnames(Beta_SD) <- predID
          rownames(Beta_SD) <- Outcomes
          
        }
        #     lme4::VarCorr(lmerResults[[1]])
        #     
        # Subject:
        #     if (!orthogonal){
        #       
        if (!AR){
          rans <- lapply(lmerResults, function(x)ranef(x)[[idvar]][,predID])  
          Beta_subject <- lapply(seq_len(nrow(rans[[1]])),function(i){
            Beta <- Beta_fixed + do.call(rbind,lapply(rans,function(x)x[i,]))
            colnames(Beta) <- predID
            rownames(Beta) <- Outcomes
            Beta
          })
        } else {
          rans <- lapply(lmerResults, function(x)ranef(x)[[idvar]][,2])  
          Beta_subject <- lapply(seq_len(length(rans[[1]])),function(i){
            Beta <- Beta_fixed + diag(sapply(rans,function(x)x[i]))
            colnames(Beta) <- predID
            rownames(Beta) <- Outcomes
            Beta
          })
        }

        
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
    
    
    
    # Now reorder everything to array:
    
    
    
    
    # Using least-squares:
    #     betaIDs <- as.numeric(rownames(ranef(lmerResults[[1]])[[idvar]]))
    #     
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
    #       #       t(C) %*% L %*% corpcor::pseudoinverse(t(L) %*% L)
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
    
    #### OBTAINING THETA #####
    # Estimate individual via correlating residuals:
    resids <- lapply(seq_along(Outcomes),
                     function(i){
                       resid <- residuals(lmerResults[[i]])
                       id <- names(resid)
                       df <- data.frame(foo = resid, id = id)
                       left_join(data.frame(id =  rownames(augData)), df, by = "id")[,2]
                     })
    resid <- as.data.frame(do.call(cbind,resids))
    names(resid) <- Outcomes
    
    if (contemporaneous %in% c("fixed","unique")){
      # If unique, via posthoc estimation:
      
      ### Compute Theta ####
      Theta_obtained <- matrix(NA, nVar, nVar)
      # diag(Theta_obtained)  <- sapply(lmerResults, lme4::sigma)^2
      diag(Theta_obtained)  <- sapply(lmerResults, sigma)^2
      
      if (contemporaneous == "unique"){
        # Compute observed residuals covariances:
        Theta_posthoc <- lapply(unique(augData[[idvar]]),function(id){
          cov(resid[augData[[idvar]] == id,Outcomes],use = "pairwise.complete.obs")
        })
        
        # abind all and compute means:
        Theta_fixed_posthoc <- apply(do.call(abind,c(Theta_posthoc,along=3)),1:2,mean,na.rm=TRUE)
        
        # Rescale to fit obtained Theta:
        D <- sqrt(diag(diag(Theta_obtained)))
        Theta_fixed <- D %*% cov2cor(Theta_fixed_posthoc) %*% D
        
        Results[["Theta"]] <- modelCov(cov = modelArray(mean = Theta_fixed, subject = Theta_posthoc))
      } else {
        
        Theta_fixed_posthoc <- cov(resid, use = "pairwise.complete.obs")
        D <- sqrt(diag(diag(Theta_obtained)))
        Theta_fixed <- D %*% cov2cor(Theta_fixed_posthoc) %*% D
        
        Theta_posthoc <- lapply(unique(augData[[idvar]]),function(id){
          Theta_fixed
        })
        
        Results[["Theta"]] <- modelCov(cov = modelArray(mean = Theta_fixed, subject = Theta_posthoc))
      }
      
    } else {
      ### TWO STEP METHOD ###
      resid[[idvar]] <- augData[[idvar]]
      
      # Output list:
      lmerResults2 <- list()
      
      if (verbose){
        message("Estimating contemporaneous effects")
      }
      # start nodewise loop:
      if (nCores > 1){
        clusterExport(cl, c("Outcomes", "model", "temporal", "contemporaneous", "idvar", "resid"), envir = environment())
        
        # Run loop:
        lmerResults2 <-  parLapply(cl, seq_along(Outcomes), function(i){
          
          # Setup model:
          mod <- paste0(
            Outcomes[i],
            " ~ 0 + ",
            paste(Outcomes[-i],collapse="+"),
            " + ( 0 + ",
            paste(Outcomes[-i],collapse="+"),
            ifelse(contemporaneous  == "orthogonal","||","|"),
            idvar,
            ")"
          )        
          
          
          # Formula:
          formula <- as.formula(mod)
          
          # Run lmer:
          return(suppressWarnings(lmer(formula, data = resid,REML=FALSE, ...)))
        })
        
        # Stop the cluster:
        stopCluster(cl)
        
      } else {
        if (verbose){
          pb <- txtProgressBar(min = 0, max = length(Outcomes), style = 3)
        }
        for (i in seq_along(Outcomes)){
          
          # Setup model:
          mod <- paste0(
            Outcomes[i],
            " ~ 0 + ",
            paste(Outcomes[-i],collapse="+"),
            " + ( 0 + ",
            paste(Outcomes[-i],collapse="+"),
            ifelse(contemporaneous  == "orthogonal","||","|"),
            idvar,
            ")"
          )        
          
          
          # Formula:
          formula <- as.formula(mod)
          
          # Run lmer:
          lmerResults2[[i]] <- suppressWarnings(lmer(formula, data = resid,REML=FALSE, ...))
          
          if (verbose){
            setTxtProgressBar(pb, i)
          }
        }
        if (verbose){
          close(pb)
        }
      }
      
      
      
      
      
      ### Gamma_Theta is the least squares regression matrix:
      Gamma_Theta_fixed <- matrix(0, nVar, nVar)
      for (i in 1:nVar){
        Gamma_Theta_fixed[i,-i] <- lme4::fixef(lmerResults2[[i]])[Outcomes[-i]]
      }
      colnames(Gamma_Theta_fixed) <- Outcomes
      rownames(Gamma_Theta_fixed) <- Outcomes
      
      # SE; SD; P; subject
      
      # SE:
      Gamma_Theta_SE <- matrix(0, nVar, nVar)
      for (i in 1:nVar){
        Gamma_Theta_SE[i,-i] <- arm::se.fixef(lmerResults2[[i]])[Outcomes[-i]]
      }
      colnames(Gamma_Theta_SE) <- Outcomes
      rownames(Gamma_Theta_SE) <- Outcomes
      
      # P:
      Gamma_Theta_P <- 2*(1-pnorm(abs(Gamma_Theta_fixed/Gamma_Theta_SE)))
      colnames(Gamma_Theta_P) <- Outcomes
      rownames(Gamma_Theta_P) <- Outcomes
      
      # Error variances:
      D <- diag(1/sapply(lmerResults2,sigma)^2)
      
      # Inverse estimate:
      inv <- D %*% (diag(length(Outcomes)) - Gamma_Theta_fixed)
      
      # Average:
      inv <- (inv + t(inv))/2
      inv <- forcePositive(inv)
      
      # invert:
      Theta_fixed_cov <- corpcor::pseudoinverse(inv)
      Theta_fixed_prec <- inv
      
      colnames(Theta_fixed_cov) <- rownames(Theta_fixed_cov) <- 
        colnames(Theta_fixed_prec) <- rownames(Theta_fixed_prec)  <- Outcomes
      
      # Compute random effects:
      Gamma_Theta_subject <- Theta_subject_prec <- Theta_subject_cov <- Theta_subject_cor <- list()
      
      for (p in 1:nrow(ranef(lmerResults2[[1]])[[idvar]])){
        Gamma_Theta_subject[[p]] <- Gamma_Theta_fixed +  do.call(rbind,lapply(seq_along(lmerResults2),function(i){
          res <- rep(0,nVar)
          res[-i] <- unlist(ranef(lmerResults2[[i]])[[idvar]][p,])
          res
        }))
        
        
        # Posthoc thetas for abnormal variances:
        Theta_posthoc <- lapply(unique(augData[[idvar]]),function(id){
          cov(resid[augData[[idvar]] == id,Outcomes],use = "pairwise.complete.obs")
        })
        
        
        Theta_subject_prec[[p]] <- forcePositive(D %*% (diag(length(Outcomes)) - Gamma_Theta_subject[[p]]))
        Theta_subject_cov[[p]] <- forcePositive(corpcor::pseudoinverse(Theta_subject_prec[[p]]))
        if (sum(diag(Theta_subject_cov[[p]])) > 10*sum(diag(Theta_fixed_cov))){
          # if cov is too big, replace with sample cov:
          Theta_subject_cov[[p]]  <- Theta_posthoc[[p]]
        }
        Theta_subject_cor[[p]] <- forcePositive(cov2cor(forcePositive(Theta_subject_cov[[p]])))
      }
      
      ### Store results:
      Results[["Theta"]] <- modelCov(
        cor = modelArray(mean=cov2cor(Theta_fixed_cov),subject = Theta_subject_cor),
        cov = modelArray(mean=Theta_fixed_cov,subject = Theta_subject_cov),
        prec = modelArray(mean=Theta_fixed_prec,subject = Theta_subject_prec)
      )
      
      Results[["Gamma_Theta"]] <- modelArray(
        mean = Gamma_Theta_fixed,
        subject = Gamma_Theta_subject,
        SE = Gamma_Theta_SE,
        P = Gamma_Theta_P
      )
    }
    
    
    
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
