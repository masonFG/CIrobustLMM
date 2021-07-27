# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the function of the parametric bootstrap on classic S-estimator
# ---------------------------------------------------------------------

parametric_MLvarCompRob <- function(model, id, B, level){
  
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
    Zmatrix 				    = matrix(c(rep(1, n_obs), WITHIN), nrow = n_obs, ncol = 2)
    
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
	  my.DATA$id    = id
      formulrest 		  = as.character(formula(model))[3]
      formulboot 		= paste("ybootr ~", formulrest,"+(WITHIN|id)")
      model.bootr 		= lmer(formulboot, data = my.DATA,REML=F)
      
      if(length(model.bootr@optinfo$conv$lme4$messages) == 0){OK = TRUE}
    }
    
    # Estimates on bootsamples
    summb 					= summary(model.bootr)
    bet_boot 				= fixef(model.bootr)
    sigma2_boot 		= sigma(model.bootr)^2
    sigma2_u0_boot 	= as.matrix(summb$varcor[[1]])[1]
    
    if (length(model$K) == 3) {
      sigma2_u1_boot 		= as.matrix(summb$varcor[[1]])[4]
      covariance_boot 	= as.matrix(summb$varcor[[1]])[2]
      result 				    = c(bet_boot, sigma2_boot, sigma2_u0_boot, sigma2_u1_boot, covariance_boot)
    } 
    if (length(model$K) == 2) {
      result 		        = c(bet_boot, sigma2_boot, sigma2_u0_boot, sigma2_u1_boot)  
    }
    
    if (length(model$K) == 1) {
      result 		        = c(bet_boot, sigma2_boot, sigma2_u0_boot)  
    }
    resultr<-rbind(resultr,result)
  }
    
  # Constructing Percentile Confidence Intervals
  results 					= list(resultr)
  names(results) 		= c("effects")
  J 						    = dim(results$effects)[2]
  estim 					  = NULL
  
  for(j in 1:J){
    estim 					= c(estim, unname(quantile(results$effects[,j], (1-level)/2, na.rm = T)), unname(quantile(results$effects[,j], 1-(1-level)/2, na.rm = T)))
  }
  
  CI 						    = t(matrix(estim,2,J))
  
  
    row.names(CI) = c(names(model$fixef), "sigma2",  names(model$eta))
    
  
  colnames(CI) 	  = c("lower bound", "upper bound")
  return(list(estimation=results$effects,CI))
}
