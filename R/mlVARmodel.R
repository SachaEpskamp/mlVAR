# 1. fixed effects (means) for every parameter and a variance-covariance matrix.
# 2. Generate Beta's 
# 3. Shrink all beta's untill all eigens are in unit circle.


mlVARsim <- function(
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
  Data <- dplyr::rbind_all(lapply(seq_len(nPerson),function(i){
    Data <- simulateVAR(Betas[[i]],1,nTime,residuals = Resids[[i]])
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
  
  class(Res) <- "mlVARsim"
  return(Res)
}

# Function to sample a VAR sequence
library("mvtnorm")

simulateVAR <- function(
  pars, # a single VAR matrix or a list of VAR matrices for each lag
  lags = 1, # A sequence of lags that corresponds to the lags of each matrix in the 'pars' argument
  Nt= 100, # Number of time points to sample
  init, # Initial state, defaults to zeros
  residuals = 0.1, # Can also be a matrix of residual covariances!
  burnin = min(round(Nt / 2),100)
)
{
  if (is.matrix(pars)) pars <- list(pars)
  if (any(sapply(pars,function(x) length(unique(dim(x))) > 1 ))) stop ("non-square graph detected.")
  
  Ni <- ncol(pars[[1]])
  maxLag <- max(lags)
  
  # Set residuals:
  if (length(residuals)==1){
    residuals <- diag(residuals,Ni)
  } else if (length(residuals) == Ni){
    residuals <- diag(residuals)
  } 
  
  if (!is.matrix(residuals) && ncol(residuals) != Ni && nrow(residuals) != Ni){
    stop("'residuals' is not a square matrix")
  }
  
  totTime <- Nt + burnin
  
  if (missing(init)) 
  {
    init <- matrix(0,maxLag,Ni)
  }
  
  
  Res <- matrix(,totTime,Ni)
  Res[1:maxLag,] <- init
  for (t in (maxLag+1):(totTime))
  {
    Res[t,] <- rowSums(do.call(cbind,lapply(seq_along(lags),function(i)pars[[i]] %*% Res[t-lags[i],]))) + rmvnorm(1,rep(0,Ni), residuals)
  }    
  
  return(as.data.frame(Res[-(1:burnin),]))
}