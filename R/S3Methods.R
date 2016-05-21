print.mlVAR <- function(x,...){
  name <- deparse(substitute(x))[[1]]
  if (nchar(name) > 10) name <- "object"
  if (name=="x") name <- "object"
  
  est <- x$input$estimator
  
  cat("\nmlVAR estimation completed. Input was:\n",
      "\t- Variables:",x$input$vars,"\n",
      "\t- Lags:",x$input$lags,"\n",
      "\t- Estimator:",est,"\n")
  if (!all(x$input$lags==0)){
    cat("\t- Temporal:",x$input$temporal)
  }
  
  
  ### Tips
  cat("\n\n",
      paste0("Use summary(",name,") to inspect fit and parameter estimates (see ?summary.mlVAR)"),
      "\n",
      paste0("Use plot(",name,") to plot estimated networks (see ?plot.mlVAR)"),
      "\n",
      paste0("Use mlVARcompare(object1, object2) to compare mlVAR models (see ?mlVARcompare)")
  )
  
}


### summary method:
summary.mlVAR <- function(
  object,
  show  = c("fit","temporal","contemporaneous","between"),
  round = 3,
  ...
){
  x <- object
  est <- x$input$estimator
  
  cat("\nmlVAR estimation completed. Input was:\n",
      "\t- Variables:",x$input$vars,"\n",
      "\t- Lags:",x$input$lags,"\n",
      "\t- Estimator:",est,"\n",
      "\t- Temporal:",x$input$temporal)
  
  nVar <- length(object$input$vars)
  nLag <- length(object$input$lags)
  Lags <- object$input$lags
  vars <- object$input$vars
  
  if ("fit" %in% show){
    cat("\n\nInformation indices:\n")
    print(object$fit, row.names=FALSE)
    
    
  }
  
  # Temporal parameters
  if ("temporal" %in% show){
    
    TemporalDF <- data.frame(
      from  = rep(rep(vars,each=nVar),nLag),
      to =  rep(rep(vars,times=nVar),nLag),
      lag = rep(Lags,each=nVar^2),
      fixed = round(c(object$results$Beta$mean),round),  
      SE = round(c(object$results$Beta$SE),round),
      P = round(c(object$results$Beta$P),round),
      ran_SD = round(c(object$results$Beta$SD),round)
    )
    rownames(TemporalDF) <- NULL
    
    cat("\n\nTemporal effects:\n")
    print(TemporalDF,row.names=FALSE)
    
  } else {
    TemporalDF <- NULL
  }
  
  if ("contemporaneous" %in% show){
    
    cor <- object$results$Theta$cor$mean
    corSD <- object$results$Theta$cor$SD
    pcor <- object$results$Theta$pcor$mean
    pcorSD <- object$results$Theta$pcor$SD
    UT <- upper.tri(cor)
    
    
    cat("\n\nContemporaneous effects (posthoc estimated):\n")
    ContDF <- data.frame(
      node1 = vars[col(cor)][UT],
      node2 = vars[row(cor)][UT],
      pcor = round(pcor[UT],round),
      ran_SD_pcor = round(pcorSD[UT],round),
      cor = round(cor[UT],round),
      ran_SD_cor = round(corSD[UT],round)
    )
    
    print(ContDF,row.names=FALSE)
    
  } else {
    ContDF <- NULL
  }
  
  
  if ("between" %in% show){
    
    P <- object$results$Gamma$P
    cor <- object$results$Omega_mu$cor$mean
    pcor <- object$results$Omega_mu$pcor$mean
    UT <- upper.tri(cor)
    
    
    cat("\n\nBetween-subject effects:\n")
    BetDF <- data.frame(
      node1 = vars[col(cor)][UT],
      node2 = vars[row(cor)][UT],
      "P 1->2" = round(P[UT],round),
      "P 2->1" = round(t(P)[UT],round),
      pcor = round(pcor[UT],round),
      cor = round(cor[UT],round)
    )
    names(BetDF) <- c("v1","v2","P 1->2","P 1<-2","pcor","cor")
    
    print(BetDF,row.names=FALSE)
    
    
    
  } else {
    BetDF <- NULL
  }
  
  
  invisible(list(temporal = TemporalDF,contemporaneous = ContDF,between = BetDF ))
}





getNet <- function(x,...){
  qgraph::getWmat(plot(x,...,DoNotPlot=TRUE))
}

plot.mlVARsim <- function(x,...){
  x$results <- x$model
  x$input <- list(vars = x$vars)
  class(x) <- "mlVAR"
  plot.mlVAR(x,...)
}

plot.mlVAR <- 
  function(x, # mlVAR object
           type = c("temporal","contemporaneous","between"), # # also allows for partial matching. e.g., temp or t.
           lag = 1, # lag of temporal network
           partial = FALSE, # Show partial correlations?
           SD = FALSE, # Plots SD instead of normal parameters
           subject, # If assigned, show the network of a particulair subject instead
           order, # If assigned, re-order nodes
           nonsig = c("show","dashed","hide"), # How to handle nonsignificant edges? In Bayesian estimation, checks if 0 is inside interval.
           alpha = 0.05, # alpha value for significance test
           onlySig = FALSE, # Backward competability argument.
           layout = "spring",
           ...  #Arguments sent to qgraph
  ){
    
    # First some backward competability:
    if (type[[1]] == "fixed"){
      warning("type = 'fixed' is deprecated. type = 'temporal' instead.")
      type <- "temporal"
    }
    if (type[[1]] == "SD"){
      warning("type = 'SD' is deprecated. type = 'temporal' and SD = TRUE instead.")
      type <- "temporal"
      SD <- TRUE
    }
    if (type[[1]] == "subject"){
      warning("type = 'subject' is deprecated. type = 'temporal' and subject = ... instead")
      type <- "temporal"
      if (missing(subject)) stop("'subject' must be assigned")
    }
    if (onlySig){
      warning("'onlySig' is deprecated. Setting nonsig = 'hide'.")
      nonsig <- "hide"
    }
    
    
    # Now check for arguments:
    type <- match.arg(type)
    nonsig <- match.arg(nonsig)
    if (missing(order)){
      order <- x$input$vars
    }
    
    if (!missing(subject) && SD){
      stop("'SD' not available for subject.")
    }
    
    # If order is character, find ord as number. Else just set ord to order:
    if (is.character(order)){
      ord <- match(order,x$input$vars)
    } else {
      ord <- order
    }
    
    ### OBTAIN NETWORK TO PLOT ###
    
    # Temporal:
    if (type == "temporal"){
      if (missing(subject)){
        # Obtain fixed effects network:
        if (SD){
          NET <- t(x$results$Beta$SD[,,lag])
        } else {
          NET <- t(x$results$Beta$mean[,,lag])  
        }
        
        
        # Attempt to obtain significance:
        # Via P:
        if (!SD && any(is.na(x$results$Beta$P))){
          
          # Via CI:
          if (!any(is.na(x$results$Beta$lower)) && !any(is.na(x$results$Beta$upper))){
            SIG <- t(x$results$Beta$lower[,,lag]) > 0 |  t(x$results$Beta$upper[,,lag]) < 0
          } else {
            # No P or CI:
            SIG <- matrix(TRUE, nrow(NET), ncol(NET))
            if (nonsig != "show"){
              warning("No p-values or CI computed. Can not hide non-significant edges.")
            }
          }
        } else {
          SIG <-  t(x$results$Beta$P[,,lag]) < alpha
        }
        
      } else {
        # Obtain subject network:
        NET <- t(x$results$Beta$subject[[subject]][,,lag])
        SIG <- matrix(TRUE, nrow(NET), ncol(NET))
        if (nonsig != "show"){
          warning("Can not hide non-significant edges for subject network.")
        }
      }
    }
    
    # Contemporaneous:
    if (type == "contemporaneous"){
      
      sub <- ifelse(partial,"pcor","cor")
      
      if (missing(subject)){
        if (SD){
          NET <- x$results$Theta[[sub]]$SD
        } else {
          # Obtain fixed effects network:
          NET <- x$results$Theta[[sub]]$mean
        }
        
        # Attempt to obtain significance:
        # Via P:
        if (!SD && any(is.na(x$results$Theta[[sub]]$P))){
          
          # Via CI:
          if (!any(is.na(x$results$Theta[[sub]]$lower)) && !any(is.na(x$results$Theta[[sub]]$upper))){
            SIG <- x$results$Theta[[sub]]$lower > 0 |  x$results$Theta[[sub]]$upper < 0
          } else {
            # No P or CI:
            SIG <- matrix(TRUE, nrow(NET), ncol(NET))
            if (nonsig != "show"){
              warning("No p-values or CI computed. Can not hide non-significant edges.")
            }
          }
        } else {
          SIG <-  x$results$Theta[[sub]]$P < alpha
        }
        
      } else {
        # Obtain subject network:
        NET <- x$results$Theta[[sub]]$subject[[subject]]
        SIG <- matrix(TRUE, nrow(NET), ncol(NET))
        if (nonsig != "show"){
          warning("Can not hide non-significant edges for subject network.")
        }
      }
    }
    
    # Contemporaneous:
    if (type == "between"){
      
      sub <- ifelse(partial,"pcor","cor")
      
      if (!missing(subject)){
        stop("No subject-specific between network possible")
      }
      if (SD){
        stop("No SD for between-subjects network.")
      }
      
      # Obtain fixed effects network:
      NET <- x$results$Omega_mu[[sub]]$mean
      
      # Attempt to obtain significance:
      # Via P:
      if (any(is.na(x$results$Omega_mu[[sub]]$P))){
        
        # Via CI:
        if (!any(is.na(x$results$Omega_mu[[sub]]$lower)) && !any(is.na(x$results$Omega_mu[[sub]]$upper))){
          SIG <- x$results$Omega_mu[[sub]]$lower > 0 |  x$results$Omega_mu[[sub]]$upper < 0
        } else {
          # No P or CI:
          SIG <- matrix(TRUE, nrow(NET), ncol(NET))
          if (nonsig != "show"){
            warning("No p-values or CI computed. Can not hide non-significant edges.")
          }
        }
      } else {
        SIG <-  x$results$Omega_mu[[sub]]$P < alpha
      }
      
      
    }
    
    ### PLOT NETWORK ###
    if (nonsig == "dashed"){
      lty <- ifelse(!SIG,2,1)
    } else {
      lty <- 1
    }
    
    if (nonsig == "hide"){
      NET <- NET * SIG
    }
    
    qgraph::qgraph(NET[ord,ord],lty = lty,labels = x$input$vars[ord],layout=layout,
                   ...)
  }




# # # logLik function:
# # logLik.mlVAR <- function(object){
# #   res <- object$pseudologlik
# #   class(res) <- "logLik"
# #   attr(res, "df") <- object$df
# #   return(res)
# # } 
# 
# # Plot function:
# # plot.mlVAR_MW <- function(x, lag = 1, ...){
# #   fixef <- fixedEffects(x)
# #   
# #   # Extract only lagged variables:
# #   sub <- fixef %>% filter(grepl(paste0("^L",lag,"_"), Predictor))
# #   
# #   # make matrix:
# #   Nodes <- as.character(unique(sub$Response))
# #   nNode <- length(Nodes)
# #   Network <- matrix(0, nNode, nNode)
# #   for (i in seq_along(Nodes)){
# #     Network[,i] <- sub$effect[sub$Response==Nodes[i]][match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes)]
# #   }
# #   
# #   Graph <- qgraph(Network, labels=Nodes, ...)
# #   invisible(Graph)
# # }
# 
# ### 
# tab2net <- function(x,lag=1){
#   Nodes <- x$dep
#   x <- x[,grepl(paste0("^L",lag,"_"), names(x))]
#   x[is.na(x)] <- 0
#   x <- t(x)
#   colnames(x) <- rownames(x) <- Nodes
#   x
# }
# 
# 
# #### These functions probably should be the other way around... but this works...
# getNet <- function(x, ...){
#   qgraph:::getWmat(plot(x,...,DoNotPlot=TRUE))
# }
# 
# plot.mlVAR <- function(x, type = c("fixed","SD","subject"), lag = 1,subject,order,onlySig = FALSE,
#                 alpha, # alpha value. if missing, do bonferonni on 0.05.
#                 ...){
#   type <- match.arg(type)
#   
#   if (type[[1]]=="subject" & missing(subject)){
#     stop("'subject' is needed to plot individual network")
#   }
# 
#   # Nodes <- rownames(x$fixedEffects)
#   if (type[[1]]=="fixed"){
# 
#     fixef <- fixedEffects(x)
#     Nodes <- as.character(unique(fixef$Response))
# 
#     # Extract only lagged variables:
#     sub <- fixef %>% filter(grepl(paste0("^L",lag,"_"), Predictor))
#     
#     # make matrix:
#     ### Simple trick to remove nonsigs, edit fixed effects to be 0 if not sig:
#     if (onlySig){
#       if (missing(alpha)){
#         alpha <- 0.05 / nrow(sub)
#       }
#       sub <- sub %>% ungroup %>% mutate(effect = ifelse(p < alpha, effect, 0))
#     }
# 
#     if (!missing(order)){
#       if (!all(sort(Nodes)==sort(order)))stop("'order' must contain exact node labels")
#       Nodes <- Nodes[match(order,Nodes)]
#     }
#     nNode <- length(Nodes)
#     Network <- matrix(0, nNode, nNode)
# 
#     for (i in seq_along(Nodes)){
#       # Network[,i] <- sub$effect[sub$Response==Nodes[i]][match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes)]
#       Network[ match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes),i] <- sub$effect[sub$Response==Nodes[i]]
#     }
#     
#     Graph <- qgraph(Network, labels=Nodes, ...)
# #   } else if (type[[1]]=="se"){
# #     fixef <- fixedEffects(x)
# #     Nodes <- as.character(unique(fixef$Response))
# #     
# #     # Extract only lagged variables:
# #     sub <- fixef %>% filter(grepl(paste0("^L",lag,"_"), Predictor))
# #     
# #     # make matrix:
# # 
# #     if (!missing(order)){
# #       if (!all(sort(Nodes)==sort(order)))stop("'order' must contain exact node labels")
# #       Nodes <- Nodes[match(order,Nodes)]
# #     }
# #     nNode <- length(Nodes)
# #     Network <- matrix(0, nNode, nNode)
# #     for (i in seq_along(Nodes)){
# #       # Network[,i] <- sub$effect[sub$Response==Nodes[i]][match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes)]
# #       Network[ match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes),i] <- sub$se[sub$Response==Nodes[i]]
# #     }
# # 
# #     
# #     Graph <- qgraph(Network, labels=Nodes, ...)
# #     
#   } else if (type[[1]]=="SD"){
#   if (onlySig){
#     warning("'onlySig' argument can only be used with type = 'fixed'")
#   }
#     ranef <- randomEffects(x)
#     Nodes <- as.character(unique(ranef$Response))
#     
#     # Extract only lagged variables:
#     sub <- ranef %>% filter(grepl(paste0("^L",lag,"_"), Predictor))
# 
#     # make matrix:
#     nNode <- length(Nodes)
#     Network <- matrix(0, nNode, nNode)
# 
#     for (i in seq_along(Nodes)){
#     Network[ match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes),i] <- sub$variance[sub$Response==Nodes[i]]
#     }
#     
#     Graph <- qgraph(sqrt(Network), labels=Nodes, ...)
#     
#   }  else if (type[[1]]=="subject"){
#     if (onlySig){
#       warning("'onlySig' argument can only be used with type = 'fixed'")
#     }
#     Net <- tab2net(x$randomEffects[[subject]])
#     fixed <- getNet(x,"fixed",onlySig=FALSE)
#     
#     Graph <- qgraph(fixed+Net, labels=rownames(Net), ...)
#      
#   } else stop("'type' is not supported")
#   
#   invisible(Graph)
# }
# 
# 
# # Print and summary:
# summary.mlVAR <- function(object,...){
#   input <- object$input
#   
#   inputstr <- paste(sapply(seq_along(input),function(i)paste0(names(input)[i],":\t\t",paste(input[[i]],collapse=", "))), collapse = "\n")
#   
#   cat(paste0("==== mlVAR results ====\n",inputstr,"\n\nNumber of random effects:\t\t",length(object$randomEffects),"\n",
#              "pseudo log-likeligood:\t\t",round(object$pseudologlik,2),"\n",
#              "Degrees of Freedom:\t\t",round(object$df,2),"\n",
#              "BIC:\t\t",round(object$BIC,2)             
#              ))
# }
# 
# print.mlVAR <- function(x,...) summary.mlVAR(x,...)
# 
# 
# ### Model plot method:
# 
# plot.mlVARsim0 <- function(x, type = c("fixed","SD","subject"), lag = 1,subject,order,...){
#   if (type[[1]]=="subject" & missing(subject)){
#     stop("'subject' is needed to plot individual network")
#   }
#   
#   # Nodes <- rownames(x$fixedEffects)
#   if (type[[1]]=="fixed"){
# 
# 
#     fixef <- x$fixedEffects
#     Nodes <- x$vars
#     
#     Graph <- qgraph(t(fixef), labels=Nodes, ...)
#     
#   } else if (type[[1]]=="se"){
#     
#     stop("No standard errors in true model")
# 
#     
#   } else if (type[[1]]=="SD"){
#     
#     ranef <- x$randomEffectsSD
#     Nodes <- x$vars
# 
#     Graph <- qgraph(t(ranef), labels=Nodes, ...)
#     
#   }  else if (type[[1]]=="subject"){
#     fixef <- x$fixedEffects
#     Nodes <- x$vars
#     Net <- fixef + x$randomEffects[[subject]]
#     Graph <- qgraph(t(Net), labels=Nodes, ...)
#     
#   } else stop("'type' is not supported")
#   
#   invisible(Graph)
# }
