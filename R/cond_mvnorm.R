# stepwise multivariate normal sampler:
stepMVnorm <- function(
  n,
  mean = rep(0,ncol(cov)), # Unconditonal mean vector
  cov = diag(length(mean)), # Unconditional covariance matrix
  steps
){
  nVar <- ncol(cov)
  if (missing(steps)){
    steps <- list(seq_len(nVar))
  }
  Res <- matrix(NA,n,nVar)
  curVars <- numeric()
  for (i in seq_along(steps)){
    curVars <- steps[[i]]
    Res[,steps[[i]]] <- cond_mvnorm(Res[,curVars],mean[curVars],cov[curVars,curVars])
  }
  Res
}

# Sampler from conditional multivariate normal per row.
cond_mvnorm <- function(
  x, # matrix to inpute values in. With NA indicating missings
  mean = rep(0,ncol(cov)), # Unconditonal mean vector
  cov = diag(length(mean)) # Unconditional covariance matrix
){
  nP <- nrow(x)
  nV <- ncol(x)
  
  res <- x
  
  for (i in 1:nP){
    # Observed:
    obs <- which(!is.na(x[i,]))
    mis <- which(is.na(x[i,]))

    if (length(obs)==0){
      res[i,] <- c(rmvnorm(1,mean,cov))
    } else {
      # conditional covariance:
      # Following https://gbhqed.wordpress.com/2010/02/21/conditional-and-marginal-distributions-of-a-multivariate-gaussian/
      A <- cov[mis,mis]
      C <- cov[obs,mis]
      B <- cov[obs,obs]
      
      meanMis <- c(mean[mis] +  t(C) %*% solve(B) %*% (x[i,obs] - mean[obs]))
      covMis <- A - t(C) %*% solve(B) %*% C
      res[i,mis] <- c(rmvnorm(1,meanMis,covMis))
    }

  }
  
  return(res)
}