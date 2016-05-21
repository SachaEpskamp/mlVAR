
# Function to sample a VAR sequence

simulateVAR <- function(
  pars, # a single VAR matrix or a list of VAR matrices for each lag
  means = 0,
  lags = 1, # A sequence of lags that corresponds to the lags of each matrix in the 'pars' argument
  Nt= 100, # Number of time points to sample
  init, # Initial state, defaults to zeros
  residuals = 0.1, # Can also be a matrix of residual covariances!
  burnin
)
{
  if (is.matrix(pars)) pars <- list(pars)
  if (any(sapply(pars,function(x) length(unique(dim(x))) > 1 ))) stop ("non-square graph detected.")
  if (missing(burnin)){
    burnin <- min(round(Nt / 2),100)
  }
  Ni <- ncol(pars[[1]])
  if (length(means)==1) means <- rep(means,Ni)
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
    Res[t,] <- means + rowSums(do.call(cbind,lapply(seq_along(lags),function(i)pars[[i]] %*% (Res[t-lags[i],]-means)))) + mvtnorm::rmvnorm(1,rep(0,Ni), residuals)
  }    
  
  return(as.data.frame(Res[-(1:burnin),]))
}


# # Function to sample a VAR sequence
# library("mvtnorm")
# 
# simulateVAR <- function(
#   pars, # a single VAR matrix or a list of VAR matrices for each lag
#   lags = 1, # A sequence of lags that corresponds to the lags of each matrix in the 'pars' argument
#   Nt= 100, # Number of time points to sample
#   init, # Initial state, defaults to zeros
#   residuals = 0.1, # Can also be a matrix of residual covariances!
#   burnin = min(round(Nt / 2),100)
# )
# {
#   if (is.matrix(pars)) pars <- list(pars)
#   if (any(sapply(pars,function(x) length(unique(dim(x))) > 1 ))) stop ("non-square graph detected.")
#   
#   Ni <- ncol(pars[[1]])
#   maxLag <- max(lags)
#   
#   # Set residuals:
#   if (length(residuals)==1){
#     residuals <- diag(residuals,Ni)
#   } else if (length(residuals) == Ni){
#     residuals <- diag(residuals)
#   } 
#   
#   if (!is.matrix(residuals) && ncol(residuals) != Ni && nrow(residuals) != Ni){
#     stop("'residuals' is not a square matrix")
#   }
#   
#   totTime <- Nt + burnin
#   
#   if (missing(init)) 
#   {
#     init <- matrix(0,maxLag,Ni)
#   }
#   
#   
#   Res <- matrix(,totTime,Ni)
#   Res[1:maxLag,] <- init
#   for (t in (maxLag+1):(totTime))
#   {
#     Res[t,] <- rowSums(do.call(cbind,lapply(seq_along(lags),function(i)pars[[i]] %*% Res[t-lags[i],]))) + MASS:: mvrnorm(1,rep(0,Ni), residuals)
#   }    
#   
#   return(as.data.frame(Res[-(1:burnin),]))
# }
