# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the function of the parametric bootstrap on composite-TAU estimator
# ---------------------------------------------------------------------

parametric_cTAU <- function(model, model0, id, Time, B, level){
  
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
      
      yboot 		= rep(0, n_obs)
      
      if (length(model0$K) > 1){
        for (l in 1:n_obs){
          yboot[l] = matfixSMALL[l,]%*%betSMALL +  Zmatrix[l,]%*%br[,l] + resr[l]
        }
      } else {
        for (l in 1:n_obs){
          yboot[l] = matfixSMALL[l,]%*%betSMALL  +  Zmatrix[l,]*br[l] + resr[l]
        }
      }
      
      my.DATA$yboot 	= yboot
	  my.DATA$id    = id
      formulrest 		  = as.character(formula(model))[3]
	  formulboot 		= paste("yboot ~", formulrest,"+(time|id)")
      model.bootr 		= lmer(formulboot, data = my.DATA)
      
      if(length(model.bootr@optinfo$conv$lme4$messages) == 0){OK = TRUE}
    }
    
    # Estimates on bootsamples
    summb 					    = summary(model.bootr)
    bet_boot 				    = fixef(model.bootr)
    sigma2_boot 			  = sigma(model.bootr)^2
    sigma2_u0_boot 			= as.matrix(summb$varcor[[1]])[1]
    if (length(model$K) > 1) {
      sigma2_u1_boot 		= as.matrix(summb$varcor[[1]])[4]
      covariance_boot 	= as.matrix(summb$varcor[[1]])[2]
      result 				    = c(bet_boot,sigma2_boot,sigma2_u0_boot,sigma2_u1_boot,covariance_boot)
    } else {
      result 				    = c(bet_boot,sigma2_boot,sigma2_u0_boot)  
    }
    
    resultr<-rbind(resultr,result)
  }
  
  
  estim                 = resultr
  if(length(model$K) > 1){
    colnames(estim)       = c(names(fixef(model)), "sigma2", "sigma2_intercept", "sigma2_time", "covariance")
  }else{
    colnames(estim)       = c(names(fixef(model)), "sigma2", "sigma2_intercept")

  }

  effet 				        = abs(estim[,index_fixef])
  
  if(length(index_fixef)>1){
    p_val      			      = colSums(effet > abs(b_tested))/B
  }else{
    p_val      			      = sum(effet > abs(b_tested))/B
    colnames(p_val)       = names(bet)[index_fixef] 
  }
  
  FinalResults               = list(estimation=estim, p_val)
  
  
  return(FinalResults)
}
