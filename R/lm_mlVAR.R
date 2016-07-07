lm_mlVAR <- 
  function(model,augData,idvar,temporal = "unique",contemporaneous = "unique", verbose=TRUE,...){
    
    stopifnot(temporal=="unique")
    
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
    
    # For every outcome compute model for every ID:
    lmResults <- list()
    
    for (i in seq_along(Outcomes)){
      
      # Unique:
      # if (temporal == "unique"){
        lmResults[[i]] <- list()
        
        # loop over subjects:
        for (p in seq_along(IDs)){
          
          sub <- augData[augData[[idvar]] == IDs[p],]
          ff <- as.formula(paste(Outcomes[i],"~",paste(predID,collapse = " + ")))
          
          lmResults[[i]][[p]] <- stats::lm(ff, sub)
        }
#         
#       } else {
# 
#         Y <- augData[,Outcomes[i]]
#         X <- augData[,Outcomes[i]]          
#         
#         lmResults[[i]] <- lm(Y~X)
#         
#       }
    }
    
    # Store results:
    Results <- list()
    
    
    # Obtain coefficients::
    Coefs <- lapply(lmResults,lapply,stats::coef)
    
    # Reorder coefs to per subject:
    Coefs <- lapply(seq_along(IDs),function(i){
      do.call(rbind,lapply(Coefs,'[[',i))
    })
    
    # Means:
    Mus <- lapply(Coefs,function(x)x[,1])
    
    # Store:
    Results[["mu"]] <- modelArray(subject = Mus)
    Results[["Omega_mu"]] <-  modelCov(
      cov = modelArray(mean=forcePositive(cov(do.call(rbind, Mus), use = "pairwise.complete.obs")))
    )
    
    # Beta:
    Betas <- lapply(Coefs,function(x)x[,-1,drop=FALSE])
    
    # Store:
    Results[["Beta"]] <- modelArray(subject = lapply(Betas,BetatoArray, mod=withinMod, out=Outcomes ))
    
    # Compute THETA
    
    
    # Make residuals data frame:
    
    Resids <- do.call(rbind,lapply(seq_along(IDs),function(id){
      mat <- cbind(
        do.call(cbind,lapply(seq_along(lmResults),function(i){
          stats::resid(lmResults[[i]][[id]])
        })),
        IDs[[id]]
      )
      colnames(mat) <- c(Outcomes,idvar)
      mat
    }))
    
    
    if (contemporaneous %in% c("fixed","unique")){
      # If unique, via posthoc estimation:
      
      if (contemporaneous == "unique"){
        # Compute observed residuals covariances:
        Theta_subject <- lapply(seq_along(IDs),function(i){
          Theta <- forcePositive(cov(Resids[Resids[,idvar] == IDs[i], Outcomes],use="pairwise.complete.obs"))
          colnames(Theta) <- rownames(Theta) <- Outcomes
          Theta
        })
        
        Results[["Theta"]] <- modelCov(cov = modelArray(subject = Theta_subject))
      } else {
        Theta <- forcePositive(cov(Resids[, Outcomes],use="pairwise.complete.obs"))
        colnames(Theta) <- rownames(Theta) <- Outcomes
        Results[["Theta"]] <- modelCov(cov = modelArray(mean = Theta,subject=lapply(seq_along(IDs),function(x)Theta)))
      }
      
    } else {
      ### TWO STEP METHOD ###
      # Output list:
      lmerResults2 <- list()
      
      if (verbose){
        message("Estimating contemporaneous effects")
        pb <- txtProgressBar(min = 0, max = length(Outcomes), style = 3)
      }
      # start nodewise loop:
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
        lmerResults2[[i]] <- suppressWarnings(lmer(formula, data = as.data.frame(Resids),REML=FALSE, ...))
        
        if (verbose){
          setTxtProgressBar(pb, i)
        }
      }
      if (verbose){
        close(pb)
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
        
        Theta_subject_prec[[p]] <- forcePositive(D %*% (diag(length(Outcomes)) - Gamma_Theta_subject[[p]]))
        Theta_subject_cov[[p]] <- forcePositive(corpcor::pseudoinverse(Theta_subject_prec[[p]]))
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
    names(lmResults) <- Outcomes
    
    fit <- data.frame(
      var = Outcomes,
      aic = sapply(lmResults,function(x)sum(sapply(x,AIC))),
      bic = sapply(lmResults,function(x)sum(sapply(x,BIC)))
    )
    rownames(fit) <- NULL
    
    
    # Combine results:
    Results <- list(
      results = Results,
      output = lmResults,
      fit = fit,
      data = augData,
      model = model
      # DIC = samples$BUGSoutput$DIC,
      # DIC_pD = samples$BUGSoutput$pD
    )
    
    class(Results) <- "mlVAR"
    
    return(Results)
    
  }
