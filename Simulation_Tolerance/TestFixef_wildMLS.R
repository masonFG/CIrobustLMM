# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the function of the wild bootstrap on classic S-estimator
# ---------------------------------------------------------------------

wild_REMLvarComp <- function(model, model0, id, Time, B, level){
  
  # Number of random effects
  if(length(model$K) > 1){
    P 	= 2
  } else{
    P 	= 1
  }
  
  # Dataset informations
  my.DATA 	= model$model
  myvar 		= all.vars(model$terms)
  y 			  = as.matrix(my.DATA[myvar[1]])
  n_obs     = model$nobs
  index_fixef = which(names(model$fixef)!=names(model0$fixef))
  n_fixef     = length(model$fixef)
  n_fixef0    = length(model0$fixef)
  matfix 	  = as.matrix(matrix(model$X, n_obs, n_fixef)) 
  n 			    = nrow(my.DATA)
  matfixSMALL = as.matrix(matrix(model0$X, n_obs, n_fixef0))
  b_tested    = model$fixef[index_fixef]
  
  # Estimates on original sample
  summ 			  = summary(model)
  bet 			  = as.vector(unname(model$fixef))
  betSMALL    = as.vector(unname(model0$fixef))
  
  
  # Define the probability of weights
  v1 			= -(sqrt(5) - 1) /2 
  v2 			= (sqrt(5) + 1) /2 
  p1 			= (sqrt(5) + 1) /(2 * sqrt(5)) 
  p2 			= 1 - p1
  
  # Define the disturbances vector for each participant
  rand 		= table(id)
  TT 			= length(rand)
  nt 			= unname(rand)
  n_obs 	= model$nobs
  X 			= as.matrix(matfixSMALL)
  XX 			= X                                                 
  Xt 			= as.matrix(unname(vector("list", TT)) )                                
  yt 			= as.matrix(vector("list", TT))
  rt_hat 	= as.matrix(vector("list", TT))
  Pt 			= as.matrix(vector("list", TT))
  Int 		= as.matrix(vector("list", TT))
  tXX 		= as.matrix(unname(solve(t(X)%*%X)))  
  
  for (i in 1:TT) {
    Xt[[i]] 		= XX[1:nt[i],]                
    XX 					= XX[-(1:nt[i]),]             
    yt[[i]] 		= y[1:nt[i]]                  
    y 					= y[-(1:nt[i])]               
    rt_hat[[i]] = yt[[i]] - Xt[[i]]%*%betSMALL   
    Pt[[i]] 		= Xt[[i]]%*%tXX%*%t(Xt[[i]])  
    Int[[i]] 		= diag(1,nt[i],nt[i])         
  }
  
  # Matrix for results
  rt_b 						        = vector("list", TT)
  yt_b 						        = vector("list", TT)
  
  result 	            = NULL
  resultr 	          = NULL
  
  for(b in 1:B){
    OK = FALSE
    
    while(!OK){
      
      # Generating the bootsamples observations (y*)
      wt 				= sample(c(v1, v2), TT, replace = TRUE, prob = c(p1, p2)) 
      
      for (tt in 1:TT) {
        rt_b[[tt]] = sqrt(diag(ginv(Int[[tt]]-Pt[[tt]])))*rt_hat[[tt]]*wt[tt] 
        yt_b[[tt]] = Xt[[tt]]%*%betSMALL + rt_b[[tt]]
      }
      
      # Fitting the bootsample
      yboot 			  = unlist(yt_b)
      my.DATA$id    = id
      my.DATA$yboot = yboot
      formulrest 		= as.character(formula(model))[3]
      formulboot 		= paste("yboot ~", formulrest,"+(time|id)")
      model.bootr 		= lmer(formulboot, data = my.DATA)
      
      if(length(model.bootr@optinfo$conv$lme4$messages) == 0){OK = TRUE}
    }
    
    # Estimates on bootsamples
    summb 					= summary(model.bootr)
    bet_boot 				= fixef(model.bootr)
    sigma2_boot 		= sigma(model.bootr)^2
    sigma2_u0_boot 	= as.matrix(summb$varcor[[1]])[1]
    
    if (length(model$K) > 1) {
      sigma2_u1_boot 		= as.matrix(summb$varcor[[1]])[4]
      covariance_boot 	= as.matrix(summb$varcor[[1]])[2]
      result 				    = c(bet_boot, sigma2_boot, sigma2_u0_boot, sigma2_u1_boot, covariance_boot)
    } else {
      result 		        = c(bet_boot, sigma2_boot, sigma2_u0_boot)  
    }
    
    resultr<-rbind(resultr,result)
  }


  estim                 = resultr
  if(length(model$K) > 1){
    colnames(estim)       = c(names(model$fixef), "sigma2", "sigma2_intercept", "sigma2_time", "covariance")
  }else{
    colnames(estim)       = c(names(model$fixef), "sigma2", "sigma2_intercept")

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
