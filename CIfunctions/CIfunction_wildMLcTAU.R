# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: March, 2021
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the function of the wild bootstrap on composite TAU-estimator
# ---------------------------------------------------------------------

wild_varComp <- function(model, id, Time, B, level){
  
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
  effetsfix = names(model$fixef)[-1] 
  n_obs     = model$nobs
  matfix 	  = as.matrix(matrix(model$X, n_obs, length(model$fixef))) 
  n 			    = nrow(my.DATA)
  
  # Estimates on original sample
  summ 			  = summary(model)
  bet 			  = as.vector(unname(model$fixef))
  
  if(length(model$K) > 1){
    if(length(model$K) > 2){
      ranefmatrix 	= matrix(c(model$eta[1], model$eta[3], model$eta[3], model$eta[2]), 2, 2)
      est 				  = c(bet, model$eta0, model$eta[1], model$eta[2], model$eta[3])
      names(est) 		= c("intercept", "time", "sigma2", "sigma2_u0", "sigma2_u1", "covariance")
    }else{
      ranefmatrix 	= matrix(c(model$eta[1], 0, 0, model$eta[2]), 2, 2)
      est 				  = c(bet, model$eta0, model$eta[1], model$eta[2])
      names(est) 		= c("intercept", "time", "sigma2", "sigma2_u0", "sigma2_u1")}
  }else{
    ranefmatrix 	= model$eta[1]
    est 					= c(bet, model$eta0, model$eta[1])
    names(est) 		= c("intercept", "time", "sigma2", "sigma2_u0")
  }
  
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
  X 			= as.matrix(matfix)
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
    rt_hat[[i]] = yt[[i]] - Xt[[i]]%*%bet   
    Pt[[i]] 		= Xt[[i]]%*%tXX%*%t(Xt[[i]])  
    Int[[i]] 		= diag(1,nt[i],nt[i])         
  }
  
  # Matrix for results
  rt_b 						        = vector("list", TT)
  yt_b 						        = vector("list", TT)
  
  # Bootstrap scheme
  # Calculate the number of cores
  no_cores 		  = detectCores() - 1
  
  # Initiate cluster
  cl 			      = makeCluster(no_cores)
  registerDoParallel(cl)
  
  tAB = foreach (b = 1:B,
                 .combine = "rbind",
                 .export = ls(globalenv()),
                 .packages = c("lme4", "MASS")) %dopar% {
                   OK = FALSE
    
    while(!OK){
      
      # Generating the bootsamples observations (y*)
      wt 				= sample(c(v1, v2), TT, replace = TRUE, prob = c(p1, p2)) 
      
      for (tt in 1:TT) {
        rt_b[[tt]] = sqrt(diag(ginv(Int[[tt]]-Pt[[tt]])))*rt_hat[[tt]]*wt[tt] 
        yt_b[[tt]] = Xt[[tt]]%*%bet + rt_b[[tt]]
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
    
    result
  }
  stopImplicitCluster()
  stopCluster(cl) # shut down the cluster

  # Constructing Percentile Confidence Intervals
  results 					= list(tAB)
  names(results) 		= c("effects")
  J 						    = dim(results$effects)[2]
  estim 					  = NULL
  
  for(j in 1:J){
    estim 					= c(estim, unname(quantile(results$effects[,j], (1-level)/2, na.rm = T)), unname(quantile(results$effects[,j], 1-(1-level)/2, na.rm = T)))
  }
  
  CI 						    = t(matrix(estim,2,J))
  
  if(length(model$K) > 1){
    if(length(model$K) > 2){
    row.names(CI) = c(names(model$fixef), "sigma2", "sigma2_intercept", "sigma2_time", "covariance")
    } else{
    row.names(CI) = c(names(model$fixef), "sigma2", "sigma2_intercept", "sigma2_time") 
    }
  }else{
    row.names(CI) = c(names(model$fixef), "sigma2", "sigma2_intercept")
    } 
  
  colnames(CI) 	  = c("lower bound", "upper bound")
  return(CI)
}
