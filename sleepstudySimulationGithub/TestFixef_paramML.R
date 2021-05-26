# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the function of the parametric bootstrap on with lmer or rlmer
# ---------------------------------------------------------------------

param_lmer <- function(model, model0, B, level){
  
  # Number of random effects
  if (length(model@theta) > 1){                              
    P 	  = 2
  } else{
    P 	  = 1
  }
  
  
  # Dataset informations
  bdd 			  = model@frame
  y 			    = model@resp$y
  matfixSMALL 		  = model.matrix(model0)
  matfix       		  = model.matrix(model)
  index_fixef = which(names(matfix)!=names(matfixSMALL))
  matfix 	  	= as.matrix(matfix)
  matfixSMALL = as.matrix(matfixSMALL)
  nomID 		  = names(summary(model)$ngrps)[1]
  id 			    = bdd[nomID]
  n_subj 		  = dim(unique(id))[1]
  n 			    = nrow(bdd)
  b_tested    = fixef(model)[index_fixef]
  
  # Estimates on original sample
  summ 			  = summary(model0)
  bet 			  = as.vector(unname(fixef(model)))
  betSMALL    = as.vector(unname(fixef(model0)))

  
  if (length(model0@theta) > 1) {
    labtime 		  = names((ranef(model0))[[1]])[2]
    matriceZ 		  = as.matrix(cbind(rep(1, n), unname(bdd[labtime])))
  } else {
    matriceZ 		  = matrix(c(rep(1, n)),nrow = n, ncol = 1)
    matriceZ 		  = as.matrix(matriceZ)
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
      resr 					      = rnorm(n, 0, sigma(model0))
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
      
      if (length(model0@theta) > 1){
        for (l in 1:n){
          yboot[l] = matfixSMALL[l,]%*%betSMALL +  matriceZ[l,]%*%br[,l] + resr[l]
        }
      } else {
        for (l in 1:n){
          yboot[l] = matfixSMALL[l,]%*%betSMALL  +  matriceZ[l,]*br[l] + resr[l]
        }
      }
      
      # Fitting the bootsample
      bdd$yboot 			= yboot
      formulrest 			= as.character(formula(model))[3]
      formulboot 			= paste("yboot ~", formulrest)
      model.bootr 		= lmer(formulboot, data = bdd) 
      
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
  
  
  estim                 = resultr
  colnames(estim)       = c(names(fixef(model)), "sigma2", "sigma2_intercept", "sigma2_time", "covariance")
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
