# # logLik function:
# logLik.mlVAR <- function(object){
#   res <- object$pseudologlik
#   class(res) <- "logLik"
#   attr(res, "df") <- object$df
#   return(res)
# } 

# Plot function:
# plot.mlVAR_MW <- function(x, lag = 1, ...){
#   fixef <- fixedEffects(x)
#   
#   # Extract only lagged variables:
#   sub <- fixef %>% filter(grepl(paste0("^L",lag,"_"), Predictor))
#   
#   # make matrix:
#   Nodes <- as.character(unique(sub$Response))
#   nNode <- length(Nodes)
#   Network <- matrix(0, nNode, nNode)
#   for (i in seq_along(Nodes)){
#     Network[,i] <- sub$effect[sub$Response==Nodes[i]][match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes)]
#   }
#   
#   Graph <- qgraph(Network, labels=Nodes, ...)
#   invisible(Graph)
# }


plot.mlVAR <- function(x, type = c("fixed","se","random","subject"), lag = 1,subject,order,...){
  if (type[[1]]=="subject" & missing(subject)){
    stop("'subject' is needed to plot individual network")
  }

  # Nodes <- rownames(x$fixedEffects)
  if (type[[1]]=="fixed"){

    fixef <- fixedEffects(x)
    
    # Extract only lagged variables:
    sub <- fixef %>% filter(grepl(paste0("^L",lag,"_"), Predictor))
    
    # make matrix:
    Nodes <- as.character(unique(sub$Response))
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
  } else if (type[[1]]=="se"){
    fixef <- fixedEffects(x)
    
    # Extract only lagged variables:
    sub <- fixef %>% filter(grepl(paste0("^L",lag,"_"), Predictor))
    
    # make matrix:
    Nodes <- as.character(unique(sub$Response))
    if (!missing(order)){
      if (!all(sort(Nodes)==sort(order)))stop("'order' must contain exact node labels")
      Nodes <- Nodes[match(order,Nodes)]
    }
    nNode <- length(Nodes)
    Network <- matrix(0, nNode, nNode)
    for (i in seq_along(Nodes)){
      # Network[,i] <- sub$effect[sub$Response==Nodes[i]][match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes)]
      Network[ match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes),i] <- sub$se[sub$Response==Nodes[i]]
    }

    
    Graph <- qgraph(Network, labels=Nodes, ...)
    
  } else if (type[[1]]=="random"){
    
    ranef <- randomEffects(x)
    
    # Extract only lagged variables:
    sub <- ranef %>% filter(grepl(paste0("^L",lag,"_"), Predictor))
    
    # make matrix:
    Nodes <- as.character(unique(sub$Response))
    nNode <- length(Nodes)
    Network <- matrix(0, nNode, nNode)
    for (i in seq_along(Nodes)){
      Network[,i] <- sub$variance[sub$Response==Nodes[i]][match(gsub("^L\\d+_","",sub$Predictor)[sub$Response==Nodes[i]], Nodes)]
    }
    
    Graph <- qgraph(Network, labels=Nodes, ...)
    
  }  else if (type[[1]]=="subject"){
      
      stop("Currently not supported")    
#     Net <- x$randomEffects[[subject]][,grepl(paste0("^L",lag,"_(",paste(Nodes,collapse="|"),")$"),colnames(x$randomEffects[[subject]]))]
#     Graph <- qgraph(t(Net), labels=rownames(Net), ...)
#     
  } else stop("'type' is not supported")
  
  invisible(Graph)
}


# Print and summary:
summary.mlVAR <- function(object,...){
  input <- object$input
  
  inputstr <- paste(sapply(seq_along(input),function(i)paste0(names(input)[i],":\t\t",paste(input[[i]],collapse=", "))), collapse = "\n")
  
  cat(paste0("==== mlVAR results ====\n",inputstr,"\n\nNumber of random effects:\t\t",length(object$randomEffects),"\n",
             "pseudo log-likeligood:\t\t",round(object$pseudologlik,2),"\n",
             "Degrees of Freedom:\t\t",round(object$df,2),"\n",
             "BIC:\t\t",round(object$BIC,2)             
             ))
}

print.mlVAR <- function(x,...) summary.mlVAR(x,...)
