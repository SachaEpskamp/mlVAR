
### 
tab2net <- function(x,lag=1){
  Nodes <- x$dep
  x <- x[,grepl(paste0("^L",lag,"_"), names(x))]
  x[is.na(x)] <- 0
  x <- t(x)
  colnames(x) <- rownames(x) <- Nodes
  x
}


#### These functions probably should be the other way around... but this works...
getNet <- function(x, ...){
  qgraph::getWmat(plot(x,...,DoNotPlot=TRUE))
}

plot.mlVAR0 <- function(x, type = c("fixed","SD","subject"), lag = 1,subject,order,onlySig = FALSE,
                alpha, # alpha value. if missing, do bonferonni on 0.05.
                ...){
  
  # global dummies:
  Predictor <- NULL
  p <- NULL
  effect <- NULL
  
  type <- match.arg(type)
  
  if (type[[1]]=="subject" & missing(subject)){
    stop("'subject' is needed to plot individual network")
  }

  # Nodes <- rownames(x$fixedEffects)
  if (type[[1]]=="fixed"){

    fixef <- fixedEffects(x)
    Nodes <- as.character(unique(fixef$Response))

    # Extract only lagged variables:
    sub <- fixef %>% filter(grepl(paste0("^L",lag,"_"), Predictor))
    
    # make matrix:
    ### Simple trick to remove nonsigs, edit fixed effects to be 0 if not sig:
    if (onlySig){
      if (missing(alpha)){
        alpha <- 0.05 / nrow(sub)
      }
      sub <- sub %>% ungroup %>% mutate(effect = ifelse(p < alpha, effect, 0))
    }

    if (!missing(order)){
      if (!all(sort(Nodes)==sort(order)))stop("'order' must contain exact node labels")
      Nodes <- Nodes[match(order,Nodes)]
    }
    nNode <- length(Nodes)
    Network <- matrix(0, nNode, nNode)

    for (i in seq_along(Nodes)){
      # Network[,i] <- sub$effect[sub$Response==Nodes[i]][match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes)]
      Network[ match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes),i] <- sub$effect[sub$Response==Nodes[i]]
    }
    
    Graph <- qgraph(Network, labels=Nodes, ...)
#   } else if (type[[1]]=="se"){
#     fixef <- fixedEffects(x)
#     Nodes <- as.character(unique(fixef$Response))
#     
#     # Extract only lagged variables:
#     sub <- fixef %>% filter(grepl(paste0("^L",lag,"_"), Predictor))
#     
#     # make matrix:
# 
#     if (!missing(order)){
#       if (!all(sort(Nodes)==sort(order)))stop("'order' must contain exact node labels")
#       Nodes <- Nodes[match(order,Nodes)]
#     }
#     nNode <- length(Nodes)
#     Network <- matrix(0, nNode, nNode)
#     for (i in seq_along(Nodes)){
#       # Network[,i] <- sub$effect[sub$Response==Nodes[i]][match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes)]
#       Network[ match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes),i] <- sub$se[sub$Response==Nodes[i]]
#     }
# 
#     
#     Graph <- qgraph(Network, labels=Nodes, ...)
#     
  } else if (type[[1]]=="SD"){
  if (onlySig){
    warning("'onlySig' argument can only be used with type = 'fixed'")
  }
    ranef <- randomEffects(x)
    Nodes <- as.character(unique(ranef$Response))
    
    # Extract only lagged variables:
    sub <- ranef %>% filter(grepl(paste0("^L",lag,"_"), Predictor))

    # make matrix:
    nNode <- length(Nodes)
    Network <- matrix(0, nNode, nNode)

    for (i in seq_along(Nodes)){
    Network[ match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes),i] <- sub$variance[sub$Response==Nodes[i]]
    }
    
    Graph <- qgraph(sqrt(Network), labels=Nodes, ...)
    
  }  else if (type[[1]]=="subject"){
    if (onlySig){
      warning("'onlySig' argument can only be used with type = 'fixed'")
    }
    Net <- tab2net(x$randomEffects[[subject]])
    fixed <- getNet(x,"fixed",onlySig=FALSE)
    
    Graph <- qgraph(fixed+Net, labels=rownames(Net), ...)
     
  } else stop("'type' is not supported")
  
  invisible(Graph)
}


# Print and summary:
summary.mlVAR0 <- function(object,...){
  input <- object$input
  
  inputstr <- paste(sapply(seq_along(input),function(i)paste0(names(input)[i],":\t\t",paste(input[[i]],collapse=", "))), collapse = "\n")
  
  cat(paste0("==== mlVAR results ====\n",inputstr,"\n\nNumber of random effects:\t\t",length(object$randomEffects),"\n"
             # "pseudo log-likeligood:\t\t",round(object$pseudologlik,2),"\n",
             # "Degrees of Freedom:\t\t",round(object$df,2),"\n",
             # "BIC:\t\t",round(object$BIC,2)             
             ))
}

print.mlVAR0 <- function(x,...) summary.mlVAR0(x,...)


### Model plot method:
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
