# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the function of the parametric bootstrap on classic S-estimator
# ---------------------------------------------------------------------

parametric_S <- function(model, Time, B, level){
  
  ls(globalenv())
  if(length(model$K) > 1){
    P 	= 2
  }else{
    P 	= 1
  }
  
  my.DATA 	  = model$model
  myvar 	    = all.vars(model$terms)
  n_obs 	    = model$nobs
  matfix 	    = as.matrix(matrix(model$X, n_obs, length(model$fixef))) 
  
  if(length(model$K) > 1){
    Zmatrix 				    = matrix(c(rep(1, n_obs), my.DATA$model$temps), nrow = n_obs, ncol = 2)
    
    if(length(model$K) > 2){
      ranefmatrix 		= matrix(c(model$eta[1], model$eta[3], model$eta[3], model$eta[2]), 2, 2)
      
    }else{
      ranefmatrix 		    = matrix(c(model$eta[1], 0, 0, model$eta[2]), 2, 2)}
    
  }else{
    ranefmatrix 			  = model$eta[1]
    Zmatrix 				    = matrix(c(rep(1, n_obs)), nrow = n_obs, ncol = 1)
  }
  
  bet 		            = unname(model$fixef)
  result 	            = NULL
  resultr 	          = NULL
  
  for(b in 1:B){
    OK = FALSE
    
    while(!OK){
      
      resr 			  = rnorm(n_obs, 0, sqrt(model$eta0)) 
      ranefr 		  = mvrnorm(n, rep(0, P), ranefmatrix)
      randomeff 	= data.frame(ranefr)
      randata 		= as.data.frame(randomeff[sort(model$model$`(groups)`[,2]),])
      
      if (length(model$K) > 1){
        br 		    = t(randata)
        
      } else {
        br 		    = randata
      }
      
      ybootr 		= rep(0, n_obs)
      
      if (length(model$K) > 1){
        for (l in 1:n_obs){
          ybootr[l] = matfix[l,]%*%bet +  Zmatrix[l,]%*%br[,l] + resr[l]
        }
      } else {
        for (l in 1:n_obs){
          ybootr[l] = matfix[l,]%*%bet  +  Zmatrix[l,]*br[l] + resr[l]
        }
      }
      
      my.DATA$ybootr 	= ybootr
      formulrest 		  = as.character(formula(model))[3]
      formulboot 		  = paste("ybootr ~", formulrest)
      model.bootr 		= try(varComprob(formulboot, groups = groups, data = my.DATA, varcov = K, control = varComprob.control(lower = c(0, 0, -Inf),method = "S", psi="rocke")),silent = TRUE)
      
      if(class(model.bootr) != "try-error"){OK = TRUE}
    }
    
    result 				  = cbind(t(unname(model.bootr$fixef)), model.bootr$eta0, t(unname(model.bootr$eta)))
    resultr 			  = rbind(resultr, result)
    
    if (length(model$K) > 1) {
      pboot 			      = resultr
      colnames(pboot) 	= c(names(model.bootr$fixef), "sigma2", "sigma2_intercept", "sigma2_time", "covariance")
      
    } else{
      pboot 			      = resultr
      colnames(pboot) 	= c(names(model.bootr$fixef), "sigma2", "sigma2_intercept")   
    }
  }  
  
  estim 				        = cbind(t(unname(model$fixef)), model$eta0, t(unname(model$eta)))
  estimatesS 			      = NULL
  
  for(i in 1:dim(pboot)[2]){
    estimatesS[[i]] 		    = c(unname(quantile(pboot[,i], c(((1-level)/2), 1-((1-level)/2)), na.rm = T)))
    names(estimatesS[[i]]) 	= c("lower bound", "upper bound")
  }
  
  names(estimatesS) 	      = colnames(pboot)
  return(list(estimation=pboot,t(as.data.frame(estimatesS))))
}
