# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the function of the parametric bootstrap on REML
# ---------------------------------------------------------------------

param_lmer <- function(model, B, level){
  
  # Number of random effects
  if (length(model@theta) > 1){                              
    P 	  = 2
  } else{
    P 	  = 1
  }
  
  #interaction
  inter   = grep(pattern = ":" , names(fixef(model)), value = TRUE, fixed = TRUE)
  
  if(!is.null(inter)){
    var.inter = strsplit(inter, split = ":" , fixed = TRUE)
  }
  
  # Dataset informations
  bdd 			  = model@frame
  y 			    = model@resp$y
  effetsfix 	= attr(attr(model@frame,"terms"),"varnames.fixed")[-1]
  matfix 		  = unname(cbind(rep(1, dim(bdd)[1]), bdd[effetsfix]))
  
  if(length(var.inter) > 0){
    for (j in length(var.inter)){
      matfix  = unname(cbind(matfix,bdd[var.inter[[j]][1]]*bdd[var.inter[[j]][2]]))
    }
  }
  
  matfix 	  	= as.matrix(matfix)
  nomID 		  = names(summary(model)$ngrps)[1]
  id 			    = bdd[nomID]
  n_subj 		  = dim(unique(id))[1]
  n 			    = nrow(bdd)
  
  # Estimates on original sample
  summ 			  = summary(model)
  bet 			  = as.vector(unname(fixef(model)))
  sigma2 		  = sigma(model)^2
  sigma2_u0 	= as.matrix(summ$varcor[[1]])[1]
  
  if (length(model@theta) > 1) {
    labtime 		  = names((ranef(model))[[1]])[2]
    matriceZ 		  = as.matrix(cbind(rep(1, n), unname(bdd[labtime])))
    sigma2_u1 		= as.matrix(summ$varcor[[1]])[4]
    covariance 		= as.matrix(summ$varcor[[1]])[2]
    est 			    = c(bet, sigma2, sigma2_u0, sigma2_u1, covariance)
    names(est) 		= c("intercept", "time", "sigma2", "sigma2_u0", "sigma2_u1", "covariance")
  } else {
    matriceZ 		  = matrix(c(rep(1, n)),nrow = n, ncol = 1)
    matriceZ 		  = as.matrix(matriceZ)
    est 			    = c(bet, sigma2, sigma2_u0)
    names(est) 		= c("intercept", "time", "sigma2", "sigma2_u0")  
  }
  
  # Bootstrap scheme
  #matrix for results
  result<-NULL
  resultr<-NULL
  
  #bootstrap
  for (b in 1:B){
    OK                   <- FALSE
    
    while(!OK){
      
      # Generating the bootsample observations (y*)
      resr 					      = rnorm(n, 0, sigma(model))
      effalr 					    = mvrnorm(n_subj, rep(0, P), summ$varcor[[1]])
      effaleatoirer 			= data.frame(effalr)
      effaleatoirer$id 		= nomID[[1]]
      baser 					    = as.data.frame(effaleatoirer[rep(1:n_subj, table(id)),])
      
      if (length(model@theta) > 1){
        baser 				= baser[,-3]
        br 					  = as.matrix(t(baser))
      } else {
        br 					  = as.matrix(baser[,-2])
      }
      
      yboot 				  = rep(0, n)
      
      if (length(model@theta) > 1){
        for (l in 1:n){
          yboot[l] = matfix[l,]%*%bet +  matriceZ[l,]%*%br[,l] + resr[l]
        }
      } else {
        for (l in 1:n){
          yboot[l] = matfix[l,]%*%bet  +  matriceZ[l,]*br[l] + resr[l]
        }
      }
      
      # Fitting the bootsample
      bdd$yboot 			= yboot
      formulrest 			= as.character(formula(model))[3]
      formulboot 			= paste("yboot ~", formulrest)
      model.bootr 		= lmer(formulboot, data = bdd,REML=F) 
      
      if(length(model.bootr@optinfo$conv$lme4$messages) == 0){OK = TRUE}
    }
    
    # Estimates on bootsamples
    summb 					    = summary(model.bootr)
    bet_boot 				    = fixef(model.bootr)
    sigma2_boot 			  = sigma(model.bootr)^2
    sigma2_u0_boot 			= as.matrix(summb$varcor[[1]])[1]
    if (length(model@theta) > 1) {
      sigma2_u1_boot 		= as.matrix(summb$varcor[[1]])[4]
      covariance_boot 	= as.matrix(summb$varcor[[1]])[2]
      result 				    = c(bet_boot,sigma2_boot,sigma2_u0_boot,sigma2_u1_boot,covariance_boot)
    } else {
      result 				    = c(bet_boot,sigma2_boot,sigma2_u0_boot)  
    }
    
    resultr<-rbind(resultr,result)
  }
  
  
  # Constructing Percentile Confidence Intervals
  results 					  = list(resultr)
  names(results) 			= c("effects")
  J 						      = dim(results$effects)[2]
  estim 					    = NULL
  
  for(j in 1:J){
    estim 					  = c(estim, unname(quantile(results$effects[,j], (1-level)/2, na.rm = T)), unname(quantile(results$effects[,j], 1-(1-level)/2, na.rm = T)))
  }
  
  CI 						      = t(matrix(estim,2,J))
  
  if(length(model@theta) > 1){
    row.names(CI) 			= c(names(fixef(model)), "sigma2", "sigma2_u0", "sigma2_u1", "covariance")
  }else{
    row.names(CI) 			= c(names(fixef(model)), "sigma2", "sigma2_u0") 
  }
  
  colnames(CI) 				= c("lower bound", "upper bound")
  return(list(estimation=results$effects,CI))
}
