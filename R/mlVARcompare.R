mlVARcompare <- function(...){
  dots <- list(...)
  if (any(!sapply(dots,is,"mlVAR"))){
    stop("Object used in input that is not of class 'mlVAR'")
  }
  
 # Check if dots have different lags but same number of obersvations:
  nLags <- sapply(dots,function(x)length(x$input$lags))
  if (length(unique(nLags)) > 1){
    # Check if number of rows in data is same:
    if (length(unique(sapply(dots,function(x)nrow(x$data))))>1){
      maxLags <- deparse(dput(dots[[which.max(nLags)]]$input$lags))
      stop(paste0("Models with different lags can not be compared unless 'compareToLags' is used.\n\nRerun analyses using compareToLags = ",maxLags))
    }
  }
  
  
  # Check if vars are same in each object:
  Vars <- lapply(dots,function(x)x$input$vars)
  unVars <- unique(unlist(Vars))
  if (!all(sapply(Vars,length) == length(unVars)) & !all(sapply(Vars,function(x)all(unVars%in%x)))){
    stop("Models do not have the same variables")
  }

  # Construct data frame with fit for every variable:
  Fits <- lapply(dots,'[[','fit')
  for (i in seq_along(Fits)){
    Fits[[i]]$lags <- paste(dots[[i]]$input$lags, collapse = "; ")
    Fits[[i]]$temporalModel <- ifelse(dots[[i]]$input$AR,"AR","VAR")
    if (all(dots[[i]]$input$lags == 0)){
      Fits[[i]]$temporal <- ""     
    } else {
      Fits[[i]]$temporal <- dots[[i]]$input$temporal
    }

  }
  
  Tables <- lapply(unVars,function(v){
    tab <- do.call(rbind,lapply(Fits,function(f)f[f$var==v,]))
    rownames(tab) <- NULL
    tab
  })
  
  names(Tables) <- unVars
  class(Tables) <- c("mlVARcompare","list")
  
  return(Tables)
}


print.mlVARcompare <- function(x,...){
  
  # cat("mlVAR model comparison\n\n")
  for (i in seq_along(x)){
    cat("\nModels for variable",names(x)[i],"\n")
    print(x[[i]],row.names=FALSE)
    
    cat("Best fitting model AIC:",paste0("lags: ",x[[i]]$lags," & temporal: ",x[[i]]$temporal," & temporal model: ",x[[i]]$temporalModel)[which.min(x[[i]]$aic)],"\n")
    cat("Best fitting model BIC:",paste0("lags: ",x[[i]]$lags," & temporal: ",x[[i]]$temporal," & temporal model: ",x[[i]]$temporalModel)[which.min(x[[i]]$bic)],"\n")    
    cat("\n")
  }
  
  # Overal fit:
  cat("\nCounts for best model fit over all variables\n")
  overalTab <- do.call(rbind,x) %>% group_by_("var") %>%
    mutate_(bestAIC = ~aic == min(aic), bestBIC = ~bic == min(bic)) %>%
    ungroup %>% 
    group_by_("lags","temporal","temporalModel") %>% summarize_(nAIC = ~sum(bestAIC), nBIC = ~sum(bestBIC)) %>%
    as.data.frame
  print(overalTab,row.names = FALSE)
  
}
