
fixef.psyvar <- function(object) object$fixedEffects
se.fixef.psyvar <- function(object) object$se.fixedEffects

summary.psyvar <- function(x) x[c('coef','se.coef','pvals')]
coef.psyvar <- function(x) x$coef

# logLik function:
logLik.psyvar <- function(object){
  res <- object$pseudologlik
  class(res) <- "logLik"
  attr(res, "df") <- object$df
  return(res)
} 

# Plot function:
plot.psyvar <- function(object, type = c("fixed","se","random","subject"), lag = 1,subject,...){
  if (type[[1]]=="subject" & missing(subject)){
    stop("'subject' is needed to plot individual network")
  }
  
  Nodes <- rownames(res$fixedEffects)
  
  if (type[[1]]=="fixed"){
    
    Net <- object$fixedEffects[,grepl(paste0("^L",lag,"_(",paste(Nodes,collapse="|"),")$"),colnames(object$fixedEffects))]
    Graph <- qgraph(t(Net), labels=rownames(Net), ...)
    
  } else if (type[[1]]=="se"){
    
    Net <- object$se.fixedEffects[,grepl(paste0("^L",lag,"_(",paste(Nodes,collapse="|"),")$"),colnames(object$se.fixedEffects))]
    Graph <- qgraph(t(Net), labels=rownames(Net), ...)
    
  } else if (type[[1]]=="random"){
    
    Net <- object$randomEffectsVariance[,grepl(paste0("^L",lag,"_(",paste(Nodes,collapse="|"),")$"),colnames(object$randomEffectsVariance))]
    Graph <- qgraph(t(Net), labels=rownames(Net), ...)
    
    
  }  else if (type[[1]]=="subject"){
    
    Net <- object$randomEffects[[subject]][,grepl(paste0("^L",lag,"_(",paste(Nodes,collapse="|"),")$"),colnames(object$randomEffects[[subject]]))]
    Graph <- qgraph(t(Net), labels=rownames(Net), ...)
    
  } else stop("'type' is not supported")
  
  invisible(Graph)
}


