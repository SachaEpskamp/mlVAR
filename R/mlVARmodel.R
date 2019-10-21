
# Graph simulator:
simGraph <- function(
  Nvar,
  sparsity = 0.5,
  parRange = c(0.5,1),
  constant = 1.1,
  propPositive = 0.5
){
  ## Approach from 
  # Yin, J., & Li, H. (2011). A sparse conditional gaussian graphical model for analysis of genetical genomics data. The annals of applied statistics, 5(4), 2630.
  
  # Empty matrix:
  trueKappa <- matrix(0,Nvar,Nvar)
  
  # Total edges:
  totEdges <- sum(upper.tri(trueKappa))
  
  # Included edges:
  nEdges <- round((1-sparsity)*totEdges)
  
  # Sample the edges:
  inclEdges <- sample(seq_len(totEdges),nEdges)
  
  # Make edges:
  trueKappa[upper.tri(trueKappa)][inclEdges] <- 1
  
  # Make edges negative and add weights:
  trueKappa[upper.tri(trueKappa)] <- trueKappa[upper.tri(trueKappa)] * sample(c(-1,1),sum(upper.tri(trueKappa)),TRUE,prob=c(propPositive,1-propPositive)) * 
    runif(sum(upper.tri(trueKappa)), min(parRange ),max(parRange ))
  
  # Symmetrize:
  trueKappa[lower.tri(trueKappa)] <- t(trueKappa)[lower.tri(trueKappa)]  
  
  # Make pos def:
  diag(trueKappa) <- constant * rowSums(abs(trueKappa))
  diag(trueKappa) <- ifelse(diag(trueKappa)==0,1,diag(trueKappa))
  trueKappa <- trueKappa/diag(trueKappa)[row(trueKappa)]
  trueKappa <- (trueKappa + t(trueKappa)) / 2
  
  return(as.matrix(qgraph::wi2net(trueKappa)))
}


# 1. fixed effects (means) for every parameter and a variance-covariance matrix.
# 2. Generate Beta's 
# 3. Shrink all beta's untill all eigens are in unit circle.


mlVARsim <- function(
  # Simulation setup:
  nPerson = 10, # # of persons.
  nNode = 5, # # of nodes 
  nTime = 100, # Or vector with time points per person. 
  lag = 1,
  
  # Arguments for parameters:
  # contemporaneous = c("wishart","fixedGGM","randomGGM"),
  # between = c("wishart","GGM"),
  thetaVar = rep(1,nNode),
  DF_theta = nNode*2,
  mu_SD = c(1,1),
  init_beta_SD = c(0.1,1),
  fixedMuSD = 1,
  shrink_fixed = 0.9,
  shrink_deviation = 0.9
){
  contemporaneous <- "wishart"
  # contemporaneous <- match.arg(contemporaneous)
  GGMsparsity = 0.5
  
  if (length(nTime)==1){
    nTime <- rep(nTime,nPerson)
  }
  # Number of temporal effects:
  nTemporal <- nNode^2 * lag
  
  # 1. Generate structures:
  # Generate Omega:

  # Simulate mu means:
  Omega_mu <- cov2cor(solve(diag(nNode)-simGraph(nNode)))
  # Omega_mu <- genPositiveDefMat(nNode, "onion", rangeVar = c(1,1))$Sigma
  # Simulate temporal:
  # Omega_Beta <- genPositiveDefMat(nTemporal, "onion", rangeVar = c(1,1))$Sigma
  Omega_Beta <- cov2cor(solve(diag(nTemporal)-simGraph(nTemporal)))
  
#   mat <- Omega_mu
#   diag(mat) <- 0
#   rowSums(mat)
#   
  Omega <- rbind(
    cbind(Omega_mu, matrix(0,nNode,nTemporal)),
    cbind(matrix(0,nTemporal,nNode), Omega_Beta)
  )

  # Omega <- genPositiveDefMat(nNode + nTemporal, "onion", rangeVar = c(1,1))$Sigma
  
  # Generate SD and scale:
  SD <- runif(nNode + nTemporal, c(rep(mu_SD[1],nNode),rep(init_beta_SD[1],nNode)), c(rep(mu_SD[2],nNode),rep(init_beta_SD[2],nNode)))
  Omega <- diag(SD) %*%Omega %*% diag(SD)

  # Generate fixed contemporaneous:
  if (contemporaneous=="wishart"){
    # Theta_fixed <- genPositiveDefMat(nNode, "onion", rangeVar = c(1,1))$Sigma
    Theta_fixed <- cov2cor(solve(diag(nNode)-simGraph(nNode)))
    Theta_fixed <- diag(sqrt(thetaVar)) %*% Theta_fixed %*% diag(sqrt(thetaVar))
    
    # 2. Generate residual covariance matrices:
    Theta <- rWishart(nPerson, DF_theta, Theta_fixed/DF_theta)
  } else {

    if (contemporaneous == "randomGGM"){
      Theta <- lapply(1:nPerson,function(x){
        net <- simGraph(nNode,GGMsparsity)
        cov2cor(solve(diag(nNode) - net))
      })
      Theta <- do.call(abind,c(Theta,along=3))
     
      Theta_fixed <- apply(Theta,1:2,mean)
    } else {
      net <- simGraph(nNode,GGMsparsity)
      Theta_fixed <- cov2cor(solve(diag(nNode) - net))
      Theta <- lapply(1:nPerson,function(x)Theta_fixed)
      Theta <- do.call(abind,c(Theta,along=3))
    }
  }
 
  
  
  # Generate fixed means:
  mu_fixed <- rnorm(nNode,0,fixedMuSD)
  # Generate fixed betas:
  beta_fixed <- rnorm(nTemporal,0)
  # set weakest 50% to zero:
  beta_fixed[order(abs(beta_fixed))[1:round(nTemporal/2)]] <- 0
  # Include auto-regressions:
  mat <- matrix(0,nNode,nNode*lag)
  diag(mat) <- 1
  beta_fixed[c(mat)==1] <- runif(sum(c(mat)==1),0,1)
  
  # 3. Generate random parameter sets:
  if (lag > 0){
    repeat{
      Pars <- rmvnorm(nPerson, c(mu_fixed,beta_fixed), sigma = Omega)
      Mus <- Pars[,1:nNode]
      
      Betas <- array(c(t(Pars[,-(1:nNode)])), c(nNode,nNode*lag,nPerson))
      
      # 4. Construct the matrices:
      if (lag>1){
        under <- cbind(diag(nNode*(lag-1)),matrix(0,nNode*(lag-1),nNode))
        
        ev <- sapply(seq_len(nPerson),function(i){
          mat <- rbind(Betas[,,i],under)
          eigen(mat)$values
        })
        
      } else {
        
        ev <- sapply(seq_len(nPerson),function(i){
          eigen(Betas[,,i])$values
        })
        
      }
      # 5. Store results:
      allEV <- c(ev)
      
      # 6. Break if all Re(ev)^2 + Im(ev)^2 < 1
      if (all(Re(ev)^2 + Im(ev)^2 < 1)){
        
        # simulate VAR for every person:
        DataList <- lapply(1:nPerson,function(p){
          
          pars <- lapply(seq_len(lag),function(l)array(c(Betas[,,p]),c(nNode,nNode,lag))[,,l])
          # If lag > 0 simulate VAR:
          if (lag > 0){
            res <- simulateVAR(pars, means = Mus[p,], lags = seq_len(lag), Nt = nTime[p],init = Mus[p,],burnin = 100,residuals = Theta[,,p]) 
          } else {
            res <- rmvnorm(nTime[p],Mus[p,],Theta[,,p])
          }
          colnames(res) <- paste0("V",1:nNode)
          res$ID <- p
          res
        })
        
        # Rbind data:
        Data <- do.call(rbind,DataList)
        
        # 10. If any absolute > 10, go to 6a
        if (!any(abs(Data[,1:nNode]) > 100)){
          break
        }
      } 
      
      # Else shrink:
      beta_fixed <- beta_fixed * shrink_fixed
      D <- diag(sqrt(diag(Omega)))
      D[-(1:nNode),-(1:nNode)] <- shrink_deviation * D[-(1:nNode),-(1:nNode)] 
      Omega <- D %*% cov2cor(Omega) %*% D
    }
    
  } else {
    Pars <- rmvnorm(nPerson, mu_fixed, sigma = Omega)
    Mus <- Pars[,1:nNode]
    Betas <- array(dim = c(0,0,nPerson))
    
    # simulate VAR for every person:
    DataList <- lapply(1:nPerson,function(p){
  

        res <- as.data.frame(rmvnorm(nTime[p],Mus[p,],Theta[,,p]))
      colnames(res) <- paste0("V",1:nNode)
      res$ID <- p
      res
    })
    
    # Rbind data:
    Data <- do.call(rbind,DataList)
    
  }
  
  
  # Create the list:
  model <- list(
    mu = modelArray(mean = mu_fixed, SD = mu_SD, subject = lapply(1:nrow(Mus),function(i)Mus[i,])),
    Beta = modelArray(mean = array(beta_fixed,c(nNode,nNode,lag)), SD = array(sqrt(diag(Omega[-(1:nNode),-(1:nNode)])),c(nNode,nNode,lag)), 
                      subject = lapply(1:nPerson, function(p)array(Betas[,,p],c(nNode,nNode,lag)))),
    Omega_mu = modelCov(
      cov = modelArray(mean = Omega[1:nNode,1:nNode])
    ),
    Theta = modelCov(
      cov = modelArray(mean = Theta_fixed, subject = lapply(1:nPerson,function(p)Theta[,,p]))
    ),
    Omega = modelCov(
      cov = modelArray(mean = Omega)
    )
    
  )
  
  
  # Data generated! Now return in sensible list:
  Results <- list(
    Data = Data,
#     beta_fixed = array(beta_fixed,c(nNode,nNode,lag)),
#     beta_SD = array(sqrt(diag(Omega[-(1:nNode),-(1:nNode)])),c(nNode,nNode,lag)),
#     mu_fixed = mu_fixed,
#     mu_SD = sqrt(diag(Omega[1:nNode,1:nNode])),
#     Theta_fixed = Theta_fixed,
#     DF_theta = DF_theta,
#     Omega = Omega,
#     Omega_mu = Omega[1:nNode,1:nNode],
#     Omega_beta = Omega[-(1:nNode),-(1:nNode)],
#     Theta = Theta,
#     Mus = Mus,
#     Betas = Betas,
    vars = paste0("V",1:nNode),
    idvar = "ID",
    lag=lag,
model=model
  )
  
  class(Results) <- "mlVARsim"
  
  return(Results)
}


### OLD ###
mlVARsim0 <- function(
  nPerson = 10, # # of persons
  nNode = 5, # # of nodes 
  nTime = 100,
  sparsity = 0, # Sparsity of the model
  parRange = c(0.22,0.4), # Parameter range of fixed effects
  propPositive = 0.5, # Proportion of positive edges excluding diagonal
  diagPositive = TRUE, # Set diagonal to positive (fixed effects)
  diagIncluded = TRUE, # Include diag no matter what
  sdRange = c(0.01,0.2), # Range of SD
  shrinkFactor = 0.95,
  residualStyle = c("full","diag"),
  residualShared = TRUE, # Shared is residual matrix for all, unique is residual matrix for each person.
  residualSDrange = c(0.05,0.1),
  verbose = TRUE
){
  
  warning("Function is deprecated and will be removed soon. Use mlVARsim instead.")
  
  residualStyle <- match.arg(residualStyle)
  
  # Number of parameters:
  nPar <- nNode^2
  
  # Parameter data frame:
  Pars <- data.frame(
    from = rep(seq_len(nNode),times=nNode),
    to = rep(seq_len(nNode),each=nNode),
    type = c(ifelse(diag(nNode)==1,"auto","edge"))
  )
  
  # Determine included edges:
  Pars$included <- runif(nPar) < ifelse(Pars$type == "auto" & diagIncluded, 1, 1-sparsity)
  
  # Determine fixed effects:
  Pars$fixed <- Pars$included * ifelse(runif(nPar) < propPositive | (diagIncluded & Pars$type == "auto"), 1, -1) * 
    runif(nPar,parRange[1],parRange[2])
  
  # Determine SDs:
  Pars$SD <- Pars$included * 
    runif(nPar,sdRange[1],sdRange[2])
  
  # Generate the correlation matrix of non-zero parameters:
  parCor <- cov2cor(genPositiveDefMat(sum(Pars$included), covMethod = "onion")$Sigma)
  
  ### Generate Betas
  repeat{
    D <- diag(Pars$SD[Pars$included])
    parCov <- D %*% parCor %*% D
    # Generate nPerson beta's:
    RandomEffects <- lapply(seq_len(nPerson),function(x){
      Beta <- matrix(0,nNode,nNode)
      Beta[Pars$included] <- rmvnorm(1,rep(0,sum(Pars$included)),parCov)
      return(Beta)
    })
    Betas <- lapply(RandomEffects,function(x){
      
      xx <- x + Pars$fixed
      xx[!Pars$included] <- 0
      return(xx)
    })
    
    # Test if eigenvalues are in unit circle:
    inUnit <- sapply(Betas,function(x){
      ev <- eigen(x,only.values = TRUE)$values
      all((Im(ev)^2 + Re(ev)^2) < 1)
    })
    
    if (all(inUnit)){
      break
    } else {
      Pars$fixed <- shrinkFactor * Pars$fixed 
      Pars$SD <- shrinkFactor * Pars$SD
      if (verbose) message("EV outside unit circle found; shrinking parameters")
    }
  }
  
  ### Generate residuals (sigmas):
  
  if (residualShared){
    residSDs <- diag(runif(nNode,residualSDrange[1],residualSDrange[2]))
    if (residualStyle == "full"){
      residCor <- cov2cor(genPositiveDefMat(nNode, covMethod = "onion")$Sigma)
      Resids <- lapply(seq_len(nPerson),function(x)residSDs %*% residCor %*% residSDs)
    } else {
      Resids <- lapply(seq_len(nPerson),function(x)residSDs)
    }
  } else {
    Resids <- lapply(seq_len(nPerson),function(x){
      residSDs <- diag(runif(nNode,residualSDrange[1],residualSDrange[2]))
      if (residualStyle == "full"){
        residCor <- cov2cor(genPositiveDefMat(nNode, covMethod = "onion")$Sigma)
        residSDs <- residSDs %*% residCor %*% residSDs
      }
      return(residSDs)
    })
  }

  ### Simulate data:
  Data <- dplyr::bind_rows(lapply(seq_len(nPerson),function(i){
    Data <- simulateVAR(Betas[[i]],lags=1,Nt=nTime,residuals = Resids[[i]])
    Data$ID <- i
    return(Data)
  }))
  
  
  
  fullCov <- matrix(0,nrow(Pars),nrow(Pars))
  fullCov[Pars$included,Pars$included] <- parCov
  
  Res <- list(
    parameters = Pars,
    fixedEffects = matrix(Pars$fixed, nNode, nNode),
    randomEffects = RandomEffects,
    randomEffectsSD = matrix(Pars$SD, nNode, nNode),
    Betas = Betas,
    residuals = Resids,
    parameterCovariances = fullCov, # Only of nonzero ones
    Data = Data,
    idvar = "ID",
    vars = paste0("V",seq_len(nNode))
  )
  
  class(Res) <- "mlVARsim0"
  return(Res)
}
