# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: March, 2021
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the function of the wild bootstrap on Student-ML
# ---------------------------------------------------------------------

wild_heavyLme <- function(model, Data, B, level){
  
  # Number of random effects
  if (length(model$theta) > 1){                              
    P 	= 2
  } else{
    P 	= 1
  }
  
  #interaction
  inter = grep(pattern = ":" , names(coefficients(model)), value = TRUE, fixed = TRUE)
  
  if(!is.null(inter)){
    var.inter = strsplit(inter, split = ":" , fixed = TRUE)
  }
  
  # Dataset informations
  varnames 		= all.vars(model$call)
  bdd 			= Data
  y 			= unlist(unname(bdd[varnames[1]]))
  effetsfix 	= varnames[-c(1, length(varnames)-1, length(varnames))]
  n_fixef = length(coefficients(model))
  matfixSMALL 		  = unname(cbind(rep(1, dim(bdd)[1]), bdd[effetsfix]))
  
  if(length(var.inter) > 0){
    for (j in length(var.inter)){
      matfix  = unname(cbind(matfixSMALL,bdd[var.inter[[j]][1]]*bdd[var.inter[[j]][2]]))
    }
  }
  
  matfix 	  	= as.matrix(matfix)
  matfixSMALL = as.matrix(matfixSMALL)
  bdd$id 		= model$lmeData$grp
  n 			= nrow(bdd)
  b_tested    = coefficients(model)[n_fixef]
  
  # Estimates on original sample
  summ 			= summary(model)
  bet 			= as.vector(unname(model$coefficients))
  betSMALL  = bet[-n_fixef]
  sigma2 		= summ$scale*(model$settings[3])/(model$settings[3]-2)
  sigma2_u0 	= summ$theta[1,1]*(model$settings[3])/(model$settings[3]-2)
  
  if (length(model$theta) > 1) {
    sigma2_u1 			= summ$theta[2,2]*(model$settings[3])/(model$settings[3]-2)
    covariance 			= summ$theta[1,2]*(model$settings[3])/(model$settings[3]-2)
    est 				= c(bet, sigma2, sigma2_u0, sigma2_u1, covariance)
    names(est) 			= c("intercept", "time", "sigma2", "sigma2_intercept", "sigma2_slope", "covariance")
  } else {
    est 				= c(bet, sigma2, sigma2_u0)
    names(est) 			= c("intercept", "time", "sigma2", "sigma2_intercept")  
  }
  
  # Define the probability of weights
  v1 			= -(sqrt(5) - 1) /2 
  v2 			= (sqrt(5) + 1) /2 
  p1 			= (sqrt(5) + 1) /(2 * sqrt(5)) 
  p2 			= 1 - p1 
  
  # Define the disturbances vector for each participant
  rand 			= table(model$lmeData$grp)
  TT 			= length(rand)
  nt 			= unname(rand)
  X 			= as.matrix(matfix)
  XX 			= as.matrix(unname(X))                                                 
  Xt 			= as.matrix(unname(vector("list", TT)) )                                
  yt 			= as.matrix(vector("list", TT))
  rt_hat 		= as.matrix(vector("list", TT))
  Pt 			= as.matrix(vector("list", TT))
  Int 			= as.matrix(vector("list", TT))
  tXX 			= as.matrix(unname(solve(t(X)%*%X)))          
  
  for (i in 1:TT) {
    Xt[[i]] 				= XX[1:nt[i],]                
    XX 						= XX[-(1:nt[i]),]             
    yt[[i]] 				= y[1:nt[i]]                  
    y 						= y[-(1:nt[i])]               
    rt_hat[[i]] 			= yt[[i]] - Xt[[i]]%*%bet   
    Pt[[i]] 				= Xt[[i]]%*%tXX%*%t(Xt[[i]])  
    Int[[i]] 				= diag(1, nt[i], nt[i])         
  }
  
  # Matrix for results
  rt_b 						= vector("list", TT)
  yt_b 						= vector("list", TT)
  
  # Bootstrap scheme
  #matrix for results
  result<-NULL
  resultr<-NULL
  
  #bootstrap
  for (b in 1:B){
    OK                   <- FALSE
	
    while(!OK){
	
      # Generating the bootsamples observations (y*)
	  wt 				= sample(c(v1,v2), TT, replace = TRUE, prob = c(p1, p2))  
      
      for (tt in 1:TT) {
        rt_b[[tt]] = sqrt(diag(ginv(Int[[tt]]-Pt[[tt]])))*rt_hat[[tt]]*wt[tt] 
        yt_b[[tt]] = Xt[[tt]]%*%bet + rt_b[[tt]]
      }
      
	  # Fitting the bootsample
      yboot 			= unlist(yt_b)
	    bdd$yboot 		= yboot
	    formulrest 		= as.character(model$call)[2]
	    formulfixboot 	= str_replace_all(formulrest, pattern = varnames[1], replacement = "yboot")
      formulboot 		= paste(formulfixboot,"+(temps|id)")
      model.bootr 		= lmer(formulboot, data = bdd,REML=F)
      
	
    if(length(model.bootr@optinfo$conv$lme4$messages) == 0){OK = TRUE}
  }
  
    # Estimates on bootsamples
    summb 					= summary(model.bootr)
    bet_boot 				= fixef(model.bootr)
    sigma2_boot 		= sigma(model.bootr)^2
    sigma2_u0_boot 	= as.matrix(summb$varcor[[1]])[1]
  
    if (length(model$theta) > 1) {
    sigma2_u1_boot 		= as.matrix(summb$varcor[[1]])[4]
    covariance_boot 	= as.matrix(summb$varcor[[1]])[2]
    result 				    = c(bet_boot, sigma2_boot, sigma2_u0_boot, sigma2_u1_boot, covariance_boot)
  } else {
    result 		        = c(bet_boot, sigma2_boot, sigma2_u0_boot)  
  }
  
  resultr<-rbind(resultr,result)
}
 
  results 					= list(resultr)
  names(results) 		= c("effects")
  J 						    = dim(results$effects)[2]
  estim 					  = NULL
  
  for(j in 1:J){
    estim 					= c(estim, unname(quantile(results$effects[,j], (1-level)/2, na.rm = T)), unname(quantile(results$effects[,j], 1-(1-level)/2, na.rm = T)))
  }
  
  CI 						    = t(matrix(estim,2,J))
  
  if(length(model$theta) > 1){
    row.names(CI) 			= c(names(coefficients(model)), "sigma2", "sigma2_intercept", "sigma2_time", "covariance")
  }else{
    row.names(CI) 			= c(names(coefficients(model)), "sigma2", "sigma2intercept") 
  }
  
  colnames(CI) 	  = c("lower bound", "upper bound")
  return(list(estimation=results$effects,CI))
}