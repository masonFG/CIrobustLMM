# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the function of the parametric bootstrap on classic S-estimator
# ---------------------------------------------------------------------

parametric_S <- function(model, model0, Time, B, level){
  
    ls(globalenv())
  if(length(model$K) > 1){
    P 	= 2
  }else{
    P 	= 1
  }
  
  my.DATA 	  = model$model
  n_obs 	    = model$nobs
  index_fixef = which(names(model$fixef)!=names(model0$fixef))
  matfix 	    = as.matrix(matrix(model$X, n_obs, length(model$fixef))) 
  matfixSMALL = as.matrix(matrix(model0$X, n_obs, length(model0$fixef)))
  b_tested    = model$fixef[index_fixef]
  
  if(length(model0$K) > 1){
    Zmatrix 				    = matrix(c(rep(1, n_obs), Time), nrow = n_obs, ncol = 2)
    
    if(length(model0$K) > 2){
      ranefmatrix 		= matrix(c(model0$eta[1], model0$eta[3], model0$eta[3], model0$eta[2]), 2, 2)
      
    }else{
      ranefmatrix 		    = matrix(c(model0$eta[1], 0, 0, model0$eta[2]), 2, 2)}
    
  }else{
    ranefmatrix 			  = model0$eta[1]
    Zmatrix 				    = matrix(c(rep(1, n_obs)), nrow = n_obs, ncol = 1)
  }
  
  bet 		            = unname(model$fixef)
  betSMALL            = unname(model0$fixef)
  result 	            = NULL
  resultr 	          = NULL
  
  for(b in 1:B){
    OK = FALSE
    
    while(!OK){
      
      resr 			  = rnorm(n_obs, 0, sqrt(model0$eta0)) 
      ranefr 		  = mvrnorm(n, rep(0, P), ranefmatrix)
      randomeff 	= data.frame(ranefr)
      randata 		= as.data.frame(randomeff[sort(model$model$`(groups)`[,2]),])
      
      if (length(model0$K) > 1){
        br 		    = t(randata)
        
      } else {
        br 		    = randata
      }
      
      ybootr 		= rep(0, n_obs)
      
      if (length(model0$K) > 1){
        for (l in 1:n_obs){
          ybootr[l] = matfixSMALL[l,]%*%betSMALL +  Zmatrix[l,]%*%br[,l] + resr[l]
        }
      } else {
        for (l in 1:n_obs){
          ybootr[l] = matfixSMALL[l,]%*%betSMALL  +  Zmatrix[l,]*br[l] + resr[l]
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
  
  estim                 = pboot
  effet 				        = estim[,index_fixef]
  
  if(length(index_fixef)>1){
    p_val      			      = colSums(effet > b_tested)/B
  }else{
    p_val      			      = sum(effet > b_tested)/B
    colnames(p_val)       = names(bet)[index_fixef] 
  }
  
  FinalResults               = list(estimation=estim, p_val)
  
  
  
  return(FinalResults)
}
