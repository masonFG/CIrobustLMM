# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the function of the wild bootstrap on composite-TAU estimator
# ---------------------------------------------------------------------

wild_cTAU <- function(model, id, Time, B, level){

  # Number of random effects
  if(length(model$K) > 1){
    P 	= 2
  } else{
    P 	= 1
  }
  
  # Dataset informations
  my.DATA 		= model$model
  myvar 		= all.vars(model$terms)
  y 			= as.matrix(my.DATA[myvar[1]])
  effetsfix 	= names(model$fixef)[-1]
  matfix 		= unname(cbind(rep(1,dim(my.DATA)[1]),my.DATA[effetsfix]))
  n 			= nrow(my.DATA)
  
  # Estimates on original sample
  summ 			= summary(model)
  bet 			= as.vector(unname(model$fixef))
  
  if(length(model$K) > 1){
    if(length(model$K) > 2){
      ranefmatrix 		= matrix(c(model$eta[1], model$eta[3], model$eta[3], model$eta[2]), 2, 2)
      est 				= c(bet, model$eta0, model$eta[1], model$eta[2], model$eta[3])
      names(est) 		= c("intercept", "time", "sigma2", "sigma2_u0", "sigma2_u1", "covariance")
    }else{
	ranefmatrix 		= matrix(c(model$eta[1], 0, 0, model$eta[2]), 2, 2)
    est 				= c(bet, model$eta0, model$eta[1], model$eta[2])
    names(est) 			= c("intercept", "time", "sigma2", "sigma2_u0", "sigma2_u1")}
  }else{
  ranefmatrix 			= model$eta[1]
  est 					= c(bet, model$eta0, model$eta[1])
  names(est) 			= c("intercept", "time", "sigma2", "sigma2_u0")
  }
  
  # Define the probability of weights
  v1 			= -(sqrt(5) - 1) /2 
  v2 			= (sqrt(5) + 1) /2 
  p1 			= (sqrt(5) + 1) /(2 * sqrt(5)) 
  p2 			= 1 - p1

  # Define the disturbances vector for each participant
  rand 			= table(id)
  TT 			= length(rand)
  nt 			= unname(rand)
  n_obs 		= model$nobs
  matfix 		= matrix(model$X, n_obs, length(model$K)-1) 
  X 			= as.matrix(matfix)
  XX 			= X                                                 
  Xt 			= as.matrix(unname(vector("list", TT)) )                                
  yt 			= as.matrix(vector("list", TT))
  rt_hat 		= as.matrix(vector("list", TT))
  Pt 			= as.matrix(vector("list", TT))
  Int 			= as.matrix(vector("list", TT))
  tXX 			= as.matrix(unname(solve(t(X)%*%X)))  
  
  for (i in 1:TT) {
    Xt[[i]] 			= XX[1:nt[i],]                
    XX 					= XX[-(1:nt[i]),]             
    yt[[i]] 			= y[1:nt[i]]                  
    y 					= y[-(1:nt[i])]               
    rt_hat[[i]] 		= yt[[i]] - Xt[[i]]%*%bet   
    Pt[[i]] 			= Xt[[i]]%*%tXX%*%t(Xt[[i]])  
    Int[[i]] 			= diag(1,nt[i],nt[i])         
  }
  
  # Matrix for results
  tabella_boot 				= matrix(0, length(est), B)
  rownames(tabella_boot) 	= names(est)
  rt_b 						= vector("list", TT)
  yt_b 						= vector("list", TT)
  
  # Bootstrap scheme
  for(b in 1:B){
    OK 					= FALSE
	
    while(!OK){
      
	  # Generating the bootsamples observations (y*)
	  wt 				= sample(c(v1, v2), TT, replace = TRUE, prob = c(p1, p2)) 
      
      for (tt in 1:TT) {
        rt_b[[tt]] = sqrt(diag(ginv(Int[[tt]]-Pt[[tt]])))*rt_hat[[tt]]*wt[tt] 
        yt_b[[tt]] = Xt[[tt]]%*%bet + rt_b[[tt]]
      }
      
	  # Fitting the bootsample
	  yboot 			= unlist(yt_b)
      my.DATA$yboot 	= yboot
	  formulrest 		= as.character(formula(model))[3]
      formulboot 		= paste("yboot ~", formulrest)
      model.bootr 		= try(varComprob(formulboot, groups = groups, data = my.DATA, varcov = K, control = varComprob.control(lower = c(0, 0, -Inf))), silent = TRUE)
	  
      if(class(model.bootr) != "try-error"){OK = TRUE}
    }
	
	# Estimates on bootsamples
    summb 				= summary(model.bootr)
    bet_boot 			= model.bootr$fixef
    sigma2_boot 		= model.bootr$eta0
    sigma2_u0_boot 		= model.bootr$eta[1]
	
    if(length(model$K) > 1){
      if(length(model$K) > 2){
      sigma2_u1_boot 	= model.bootr$eta[2]
      covariance_boot 	= model.bootr$eta[3]
      tabella_boot[,b] 	= c(bet_boot, sigma2_boot, sigma2_u0_boot, sigma2_u1_boot, covariance_boot)
    } else {
      sigma2_u1_boot 	= model.bootr$eta[2]
      tabella_boot[,b] 	= c(bet_boot, sigma2_boot, sigma2_u0_boot, sigma2_u1_boot)  
    }
    }else{
      tabella_boot[,b] 	= c(bet_boot, sigma2_boot, sigma2_u0_boot)  
  }
  }
  
  # Constructing Percentile Confidence Intervals
  results 				= list(tabella_boot)
  names(results) 		= c("effects")
  J 					= dim(results$effects)[1]
  estim 				= NULL
  
  for(j in 1:J){
    estim 				= c(estim, unname(quantile(t(results$effects)[,j], (1-level)/2, na.rm = T)), unname(quantile(t(results$effects)[,j], 1-(1-level)/2, na.rm = T)))
  }
  
  CI 					= t(matrix(estim, 2, J))
  
  if(length(model$K) > 1){
    if(length(model$K) > 2){
    row.names(CI) 		= c(names(model$fixef), "sigma2", "sigma2_intercept", "sigma2_time", "covariance")
    } else{
    row.names(CI) 		= c(names(model$fixef), "sigma2", "sigma2_intercept", "sigma2_time") 
    }
  }else{
    row.names(CI) 		= c(names(model$fixef), "sigma2", "sigma2_intercept")
    } 
  
  colnames(CI) 			= c("lower bound", "upper bound")
  return(CI)
}