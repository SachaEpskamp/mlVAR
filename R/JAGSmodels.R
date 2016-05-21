# Model description:
# JAGSModel_:
# M: multivariate // U: Univariate
# S: Shared contemporaenous // U: Unique contemporaneous
# F: Only fixed temporal // R: Fixed + random temporal
# O: Orthogonal // C: Correlated

# No this is nonsense... One Model to rule them all, One Model to find them,
# One Model to bring them all and in the darkness bind them.

# @MACRO@ <- replace in R with something.
JAGSmodel <- '
  model{
    #   Data = The dataset
    #   nVar = Number of variabless
    #   nPerson = Number of persons
    #   R_Tau_Theta = Prior scale matrix for Theta
    #   REprior = Prior scale matrix for random effects
    #   nEffect = Number of random effects
    #   randMean = Vector of zeroes
    #   nLag = Number of lags
    #   nDaysPerPerson = Vector with number of days per person
    #   nMeasurementPerDay = Matrix with number of measurements per day
    
    
    ### Fixed effects ###
    # Just a flat prior
    for (i in 1:nEffect){
      Fixed[i] ~ dnorm(fixedPriorMean[i],1/fixedPriorVar[i])
    }
    
    for (i in 1:nVar){
      mu_fixed[i] <- Fixed[i]
      
      for (j in 1:nVar){
        for (l in 1:nLag){
          # Beta_ix[i,j,l] <- nVar + (l-1)*nVar^2 + (j-1)*nVar + i
          Beta_fixed[i,j,l] <- Fixed[Beta_ix[i,j,l]]

          # IF ORTHOGONAL & temporal != "shared": SD Distribution (if not this is done later)
          @@ORTHSDDIST@@
          # sd_Beta_inverse[i,j,l] ~ dgamma(0.01,0.01)
          # sd_Beta[i,j,l] <- 1/sd_Beta_inverse[i,j,l] 

          # If temporal = "shared", just put this to zeroes:
          @@TEMPORALSD@@
        }
      }
    }
    
    ### Random effects ###
    # Vector of random effects
    # If orthogonal = TRUE OR temporal = "shared", Omega only contains block for means

    Omega_inverse[1:nRandom, 1:nRandom] ~ dwish(nRandom*REprior[1:nRandom, 1:nRandom], nRandom) 
    Omega[1:nRandom, 1:nRandom] <- inverse(Omega_inverse[1:nRandom, 1:nRandom])  
    
    ### Theta structure ###
    @@THETA_SHARED@@
      # If Theta is shared:
      # Theta_inverse[1:nVar,1:nVar] ~ dwish(R_Tau_Theta[1:nVar,1:nVar], nVar)
      # Theta[1:nVar,1:nVar] <- inverse(Theta_inverse[1:nVar,1:nVar])
      
      for (p in 1:nPerson){
        
        @@THETA_UNIQUE@@
          # If Theta is unique:
          # Theta_inverse[1:nVar,1:nVar,p] ~ dwish(R_Tau_Theta[1:nVar,1:nVar], nVar)
          # Theta[1:nVar,1:nVar,p] <- inverse(Theta_inverse[1:nVar,1:nVar,p])
          
          # Distribution of random parameters:
         Pars[1:nRandom,p] ~ dmnorm(Fixed[1:nRandom], Omega_inverse[1:nRandom, 1:nRandom])

        for (i in 1:nVar){
          mu[i,p] <- Pars[i,p]
          
          for (j in 1:nVar){
            for (l in 1:nLag){
              # If orthogonal, assign distirbution to relevant PAR
              @@ORTHDIST@@
              # If orthogonal:
              # Pars[Beta_ix[i,j,l],p] ~ dnorm(Fixed[Beta_ix[i,j,l]],sd_Beta_inverse[i,j,l])
  
              Beta[i,j,l,p] <- @@TEMPORALPARS@@ 
                # if temporal == "unique": Pars[Beta_ix[i,j,l],p]
                # if temporal == "shared": Beta_fixed[i,j,l]
            }
          }
        }
      }
    
    
    ### LIKELIHOOD ###
    # For every person:
    for (p in 1:nPerson){
      # For every day of person p:
      for (d in 1:nDaysPerPerson[p]){
        for (t in 1:nLag){
            for (i in 1:nVar){
              # Flat prior:
              Data[t,i,d,p] ~ dnorm(mu[i,p],0.001)
            }
        }

        # for every time point in that day:
        for (t in (nLag+1):nMeasurementPerDay[d,p]){
          
          # Impute correct model based on lag and theta matrix:
          Data[t,1:nVar,d,p] ~  dmnorm(@@MODELMEAN@@ ,@@MODELVAR@@ )
          
        }
      }
    }
    
    ### Matrices to extract ###
    
    ### Random effect variances/SDs and covariances ###
    for (i in 1:nVar){
      sd_mu[i] <- sqrt(Omega[i,i])
      for (j in 1:nVar){
        Omega_mu[i,j] <- Omega[i,j]
        
        # If !orthogonal, obtain sd from varcor mat:
        @@CORSD@@
        # for (l in 1:nLag){
        #   sd_Beta[i,j,l] <-  sqrt(Omega[Beta_ix[i,j,l],Beta_ix[i,j,l]])
        # }
      }
    }
    
    ### Also save inverse of Omega_mu ###
    Omega_mu_inverse[1:nVar,1:nVar] <- inverse(Omega_mu[1:nVar,1:nVar])
  }
'


