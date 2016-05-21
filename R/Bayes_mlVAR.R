# Commented out... not yet implemented

# JAGS_mlVAR <- function(augData, vars, 
#                        idvar,
#                        lags, 
#                        dayvar, 
#                        beepvar, 
#                        orthogonal = FALSE, 
#                        method = c("multivariate","univariate"), 
#                        temporal = c("unique","shared"),
#                        contemporaneous = c("shared","unique"),
#                        verbose=TRUE,
#                        # prior = c("identity","lmer_orth","lmer_seq"),
#                        JAGSexport = FALSE,
#                        n.chain = 3,
#                        n.iter = 10000,
#                        estOmega = FALSE,
#                        ...){
#   
#   ### JAGS MODEL ###
#   
#   
#   #   augData <- augData[rowSums(is.na(augData))==0,]
#   #   
#   # prior <- match.arg(prior)
#   method <- match.arg(method)
#   contemporaneous <- match.arg(contemporaneous)
#   
#   nVar <- length(vars)
#   
# #   # Random effects prior:
# #   if (prior == "identity"){
# #     
# #     REprior <- diag(nVar + nVar^2)
# #   } else 
# #   {
# #     if (verbose)
# #     {
# #       message("computing prior")
# #     }
# #     #     if (prior == "lmer_seq"){
# #     #       
# #     #       lmerRes <- lmer_murmur(Y=Y,X=X,ID=ID,data=data,orthogonal = FALSE,verbose=verbose)
# #     #       REprior <- diag(diag(lmerRes$Sigma_randomEffects))
# #     #     } else if (prior == "lmer_orth"){
# #     #       lmerRes <- lmer_murmur(Y=Y,X=X,ID=ID,data=data,orthogonal = FALSE,verbose=verbose)
# #     #       REprior <- diag(diag(lmerRes$Sigma_randomEffects))
# #     #     } else stop("Prior method not implemented yet.")
# #     #     diag(REprior) <- pmax(diag(REprior),0.0001)
# #     
# #     
# #     # REprior <- solve(REprior)
# #   }
#   
#   
#   #### RESTRUCTURE THE DATA ####
#   # Data should be array (maxBeep x nVar x maxDay x nPerson)
#   # Observations per day per person:
#   nBeep_DayID <- augData %>% group_by_(idvar,dayvar) %>% tally
#   
#   # Maximum number of beeps:
#   maxBeep <- max( nBeep_DayID$n,na.rm=TRUE)
#   
#   # Maximum number of days:
#   nDay_ID <-  eval(substitute(dplyr::summarize_(augData %>% group_by_(idvar), 
#                                                 n = ~length(unique(dayvar))),
#                               list(dayvar = as.name(dayvar))))
#   
#   # Maximum number of days:
#   maxDay <- max(nDay_ID$n,na.rm=TRUE)
#   
#   # Number of lags:
#   nLag <- length(lags)
#   
#   # All subject ids:
#   subjects <- unique(augData[[idvar]])
#   
#   # Number of subjects:
#   nSubjects <- length(subjects)
#   
#   # Construct the data:
#   Data <- array(NA, dim = c(maxBeep, nVar,  maxDay, nSubjects))
#   
#   # Additional data needed:
#   # Vector containing per person the number of days:
#   data_nDaysPerPerson <- rep(NA,nSubjects)
#   
#   # Matrix containing
#   data_nMeasurementPerDay <- matrix(NA, maxDay, nSubjects)
#   
#   # Fill the data:
#   for (p in seq_len(nSubjects)){
#     # Subject ID:
#     id <- subjects[p]
#     
#     # How many days?
#     nDay <- nDay_ID$n[nDay_ID[[idvar]] == id]
#     
#     # Days vector:
#     Days <- unique(augData[[dayvar]][augData[[idvar]]==id])
#     
#     # Add number of days:
#     data_nDaysPerPerson[p] <- length(Days)
#     
#     # Loop over days:
#     for (d in seq_len(nDay)){
#       day <- Days[d]
#       
#       # Obtain the subset:
#       Mat <- as.matrix(augData[augData[[idvar]] == id & augData[[dayvar]] == day, vars])
#       
# #       # Make rows with partial missing all missing:
#       Mat[rowSums(is.na(Mat)) > 0,] <- NA
#       
#       # Implement in data:
#       Data[seq_len(nrow(Mat)),,d,p] <- Mat
#       
#       # Add number of observations:
#       data_nMeasurementPerDay[d,p] <- nrow(Mat)
#     }
#   }
#   
#   #   
#   #   # Lags in the model:
#   #   lags <- unique(model$lag[!is.na(model$lag)])
#   #   
#   #   Y <- as.matrix(augData[,Outcomes])
#   #   X <- array(NA, dim = c(nrow(Y), ncol(Y), length(lags)))
#   #   
#   #   for (l in seq_along(lags)){
#   #     lagMod <- UniquePredModel[UniquePredModel$lag==l,]
#   #     lagPredictors <- lagMod$predID[match(lagMod$pred,Outcomes)]
#   #     
#   #     X[,,l] <- as.matrix(augData[,lagPredictors])
#   #   }
#   #   
#   # ID <- augData[[idvar]]
#   
#   #### CHANGE THE MODEL FILE ####
#   # Contemporaenous shared or unique?
#   if (contemporaneous == "shared"){
#     # SHARED THETA
#     JAGSmodel <- gsub("@@THETA_SHARED@@",
#                       "
#     Theta_inverse[1:nVar,1:nVar] ~ dwish(R_Tau_Theta[1:nVar,1:nVar], nVar)
#     Theta[1:nVar,1:nVar] <- inverse(Theta_inverse[1:nVar,1:nVar])
# ", JAGSmodel)
#     
#     JAGSmodel <- gsub("@@THETA_UNIQUE@@","", JAGSmodel)
#     
#     JAGSmodel <- gsub("@@MODELVAR@@","Theta_inverse[1:nVar,1:nVar]", JAGSmodel)
#   } else {
#     # UNIQUE THETA PER SUBJECT 
#     JAGSmodel <- gsub("@@THETA_UNIQUE@@",
#                       "
#       Theta_inverse[1:nVar,1:nVar,p] ~ dwish(R_Tau_Theta[1:nVar,1:nVar], nVar)
#       Theta[1:nVar,1:nVar,p] <- inverse(Theta_inverse[1:nVar,1:nVar,p])
# ", JAGSmodel)
#     
#     JAGSmodel <- gsub("@@THETA_SHARED@@","", JAGSmodel)
#     
#     JAGSmodel <- gsub("@@MODELVAR@@","Theta_inverse[1:nVar,1:nVar,p]", JAGSmodel)
#   }
#   
#   ### ORTHOGONAL SETUP ###
#   if (orthogonal){
#     # Distribution for SD of beta:
#     if (temporal == "unique"){
#       JAGSmodel <- gsub("@@ORTHSDDIST@@","
#       sd_Beta_inverse[i,j,l] ~ dgamma(0.01,0.01)
#       sd_Beta[i,j,l] <- pow(sd_Beta_inverse[i,j,l],-0.5)", JAGSmodel)
#     } else {
#       JAGSmodel <- gsub("@@ORTHSDDIST@@","",JAGSmodel)
#     }
#     
#     # Orthogonal distribution: 
#     if (temporal == "unique"){
#       JAGSmodel <- gsub("@@ORTHDIST@@","
#           Pars[Beta_ix[i,j,l],p] ~ dnorm(Fixed[Beta_ix[i,j,l]],sd_Beta_inverse[i,j,l])
#                       ", JAGSmodel)       
#     } else {
#       JAGSmodel <- gsub("@@ORTHDIST@@","", JAGSmodel) 
#     }
# 
# 
#             
#   } else {
#     JAGSmodel <- gsub("@@ORTHSDDIST@@","",JAGSmodel)
#     JAGSmodel <- gsub("@@ORTHDIST@@","", JAGSmodel) 
#     
#   }
#   
#   # If orthogonal or shared temporal, ignore part of the var-cor mat:
#   if (orthogonal | temporal == "shared"){
#     # Ignore obtaining SD from varcor:
#     JAGSmodel <- gsub("@@CORSD@@","",JAGSmodel)
#   } else {
#     # Obtain SD from varcor:
#     JAGSmodel <- gsub("@@CORSD@@","
#         for (l in 1:nLag){
#            sd_Beta[i,j,l] <-  sqrt(Omega[Beta_ix[i,j,l],Beta_ix[i,j,l]])
#         }
#         ",JAGSmodel)
#   }
#   
#   
#   #### TEMPORAL MULTI-LEVEL OR FIXED ONLY ###
#   if (temporal == "unique"){
#     JAGSmodel <- gsub("@@TEMPORALPARS@@","Pars[Beta_ix[i,j,l],p]",JAGSmodel)
#     JAGSmodel <- gsub("@@TEMPORALSD@@","",JAGSmodel)
#   } else {
#     JAGSmodel <- gsub("@@TEMPORALPARS@@","Beta_fixed[i,j,l]",JAGSmodel)    
#     JAGSmodel <- gsub("@@TEMPORALSD@@","sd_Beta[i,j,l] <- 0",JAGSmodel)
#   }
#   
#   ### Construct the model mean ###
#   Mean <- "mu[1:nVar,p]"
#   lags <- lags[lags!=0]
#   add <- sapply(lags, function(l){
#     paste0("Beta[1:nVar,1:nVar,",l,",p] %*% (Data[t-",l,",1:nVar,d,p] - mu[1:nVar,p])")
#   })
#   Mean <- paste(c(Mean,add),collapse = " + ")
#   JAGSmodel <- gsub("@@MODELMEAN@@",Mean,JAGSmodel)
#   
#   # Construct Beta_ix:
#   Beta_ix <- array(nVar + seq_len(nVar^2*nLag) , dim = c(nVar, nVar, nLag))
#   
#   nRandom <- ifelse(temporal == "shared" | orthogonal, nVar, nVar + (nVar^2 * nLag))
#   ### OBTAIN INITS FROM ORTHOGONAL OR LEAST SQUARES ###
#   if (verbose){
#     message("Running LMER estimator to obtain prior and initial values.")
#   }
#     lmerRes <- mlVAR(augData,vars,idvar,lags,dayvar,beepvar,orthogonal = orthogonal)
#     
#     REprior <- c(c(lmerRes$results$mu$SD),c(lmerRes$results$Beta$SD))^2
#     # Oracle prior:
#     # REprior <- diag(Model$Omega)
# 
#     REprior[REprior < 0.001] <- min(REprior[!REprior<0.001])
#     REprior <- diag(REprior)
#     # REprior is prior guess for covariance matrix.
#     
# #     diag(solve(apply(rWishart(100, nRandom, REprior),1:2,mean)))
# #     diag(Model$Omega)
# 
#     # Mean, set to lmer:
#     fixedPriorMean <- c(c(lmerRes$results$mu$mean),c(lmerRes$results$Beta$mean))
#     fixedPriorMean <- c(c(lmerRes$results$mu$mean),c(lmerRes$results$Beta$mean))
#     fixedPriorVar <- c(c(lmerRes$results$mu$SE),c(lmerRes$results$Beta$SE))
#     
#     # Fix this prior!
#     browser()
# 
# 
#     # Starting values of fixed pars:
#     Start <- c(c(lmerRes$results$mu$mean),c(lmerRes$results$Beta$mean))
#     Init <- lapply(seq_len(n.chain),function(x){
#       list(Fixed = Start + rnorm(length(Start),0, c(c(lmerRes$results$mu$SD),c(lmerRes$results$Beta$SD))/10))
#       # list(Fixed = Start)
#     })
#   
#   # obtain data:
#   jagsData <- list(
#     Data = Data, # The dataset
#     nDaysPerPerson = data_nDaysPerPerson,
#     nMeasurementPerDay = data_nMeasurementPerDay,
#     # maxBeep = maxBeep, 
#     nVar = nVar,  
#     # maxDay = maxDay, 
#     nPerson = nSubjects,
#     R_Tau_Theta = diag(nVar),
#     REprior = REprior,
#     nEffect = nVar + nLag * nVar^2,
#     # randMean = rep(0, (nVar + (nVar^2*nLag))),
#     nLag = length(lags),
#     nRandom = nRandom,
#     Beta_ix=Beta_ix,
#     fixedPriorMean=fixedPriorMean,
#     fixedPriorVar=fixedPriorVar
#   )
#   
#   ### RUN JAGS ###
#   #     if (verbose){
#   #       message("Estimating model")
#   #     }
#   # 
#   # 
#   
#   if (JAGSexport){
#     if (verbose){
#       message("Exporting model file to 'mlVARexport' folder")
#       if (!dir.exists("mlVARexport")){
#         dir.create("mlVARexport")
#       }
#       write(JAGSmodel,"mlVARexport/JAGSmodel.txt")
#       saveRDS(jagsData,"mlVARexport/JAGSdata.rds" )
#     }
#   }
# 
#   if (verbose){
#     message("Running JAGS")
#   }
#   
#   parsToSave <- c("mu_fixed","Beta_fixed",
#                   "mu", "Beta",
#                   "sd_mu",
#                   "sd_Beta",
#                   "Omega_mu",
#                   "Theta","Theta_inverse",
#                   "Omega_mu_inverse"
#   )
#     
#     if (estOmega){
#       parsToSave <- c(parsToSave,"Omega")
#     }
#     
#   samples <- jags(jagsData, parameters.to.save = parsToSave, n.chain = n.chain, n.iter = n.iter,  inits =  Init, 
#                   model.file = textConnection(JAGSmodel),...)
# 
# 
#   
#   ### GATHER RESULTS ###
#   results <- list()
#   
#   # Mu:
#   results[['mu']] <- parSamples(samples$BUGSoutput$sims.list,
#                                 "mu_fixed","mu","sd_mu")
#   
#   
#   # Beta:
#   results[['Beta']] <- parSamples(samples$BUGSoutput$sims.list,
#                                   "Beta_fixed","Beta","sd_Beta")
#   
#   if (contemporaneous == "unique"){
#     results[['Theta']] <- covSamples(samples$BUGSoutput$sims.list, subject = "Theta")    
#   } else {
#     results[['Theta']] <- covSamples(samples$BUGSoutput$sims.list, fixed = "Theta")      
#   }
#   
#   #### Omega mu ###
#   results[['Omega_mu']] <- covSamples(samples$BUGSoutput$sims.list, fixed = "Omega_mu")
#   
#   if (estOmega){
#     
#     results[['Omega']] <- covSamples(samples$BUGSoutput$sims.list, fixed = "Omega")
#     
#   }
#   
#   # Combine results:
#   Results <- list(
#     results = results,
#     output = samples,
#     DIC = samples$BUGSoutput$DIC,
#     DIC_pD = samples$BUGSoutput$pD
#   )
#   
#   class(Results) <- "mlVAR"
#   
#   return(Results)
#   
# }