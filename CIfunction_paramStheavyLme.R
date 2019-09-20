# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the function of the Student-parametric bootstrap on Student-ML
# ---------------------------------------------------------------------

paramSt_heavyLme <- function(model, Data, B, level){
  
  # Number of random effects
  if (length(model$theta) > 1){                              
    P 	= 2
  } else{
    P 	= 1
  }
  
  # Dataset informations
  varnames 		= all.vars(model$call)
  bdd 			= Data
  y 			= unlist(unname(bdd[varnames[1]]))
  effetsfix 	= varnames[-c(1, length(varnames)-1, length(varnames))]
  matfix 		= as.matrix(unname(cbind(rep(1, dim(bdd)[1]), bdd[effetsfix])))
  bdd$id 		= model$lmeData$grp
  n 			= nrow(bdd)
  n_subj 		= length(unique(bdd$id))
  
  # Estimates on original sample
  summ 			= summary(model)
  bet 			= as.vector(unname(model$coefficients))
  sigma2 		= summ$scale*(model$settings[3])/(model$settings[3]-2)
  sigma2_u0 	= summ$theta[1,1]*(model$settings[3])/(model$settings[3]-2)
  
  if (length(model$theta) > 1) {
    labtime 		= colnames((model$ranef))[2]
    matriceZ 		= as.matrix(cbind(rep(1, n), unname(bdd[labtime])))
    sigma2_u1 		= summ$theta[2,2]*(model$settings[3])/(model$settings[3]-2)
    covariance 		= summ$theta[1,2]*(model$settings[3])/(model$settings[3]-2)
    est 			= c(bet, sigma2, sigma2_u0, sigma2_u1, covariance)
    names(est) 		= c("intercept", "time", "sigma2", "sigma2_intercept", "sigma2_slope", "covariance")
  } else {
    matriceZ 		= matrix(c(rep(1, n)), nrow = n, ncol = 1)
    est 			= c(bet, sigma2, sigma2_u0)
    names(est) 		= c("intercept", "time", "sigma2", "sigma2_intercept")  
  }
  
  # Matrix for results
  resulth 					= NULL
  result 					= NULL
  tabella_boot 				= matrix(0, length(est), B)
  rownames(tabella_boot) 	= names(est)

  # Bootstrap scheme
  for (b in 1:B){
    OK 						= FALSE
	
    while(!OK){
	
	  # Generating the bootsamples observations (y*)
      resh 					= rmt(n = n, c(0),(model$settings[3]-2)/(model$settings[3])*diag(model$scale,1,1), df = model$settings[3])
      effalh 				= rmt(n = n_subj, c(0, 0), Scatter = (model$settings[3]-2)/(model$settings[3])*model$theta, df = model$settings[3])
      effaleatoireh 		= data.frame(effalh)
      nomID 				= (labels(table(bdd$id)))
      effaleatoireh$id 		= nomID[[1]]
	  
      if (length(model$theta) > 1){
        baseh 	= as.data.frame(effaleatoireh[rep(1:nrow(effaleatoireh), table(bdd$id)),])
        baseh 	= baseh[,-3]
        bh 		= t(baseh)
      } else {
        baseh 	= as.data.frame(effaleatoireh[rep(1:nrow(effaleatoireh), table(bdd$id)),])
        bh 		= baseh[,-2]
      }
      
	  ybooth 	= rep(0,n)
      
      if (length(model$theta) > 1) {
        for (l in 1:n){
          ybooth[l] = matfix[l,]%*%bet +  matriceZ[l,]%*%bh[,l] + resh[l]
        }
      } else {
        for (l in 1:n){
          ybooth[l] = matfix[l,]%*%bet +  matriceZ[l,]*bh[l] + resh[l]
        }
      }
	  
	  # Fitting the bootsample
	  bdd$ybooth 		= ybooth
      formulrest 		= as.character(model$call)[2]
      formulfixboot 	= as.formula(str_replace_all(formulrest, pattern = varnames[1], replacement = "ybooth"))                  
      formulrandboot 	= as.formula(as.character(model$call)[3])
      model.bootr 		= heavyLme(fixed = formulfixboot, random = formulrandboot, groups = ~ id, data = bdd)
      
      if(model.bootr$converged == TRUE){OK = TRUE}
    }
	
	# Estimates on bootsamples
    summb 				= summary(model.bootr)
    bet_boot 			= unname(coefficients(model.bootr))
    sigma2_boot 		= summb$scale*(model$settings[3])/(model$settings[3]-2)
    sigma2_u0_boot 		= summb$theta[1,1]*(model$settings[3])/(model$settings[3]-2)
	
    if (length(model$theta) > 1) {
      sigma2_u1_boot 	= summb$theta[2,2]*(model$settings[3])/(model$settings[3]-2)
      covariance_boot 	= summb$theta[2,1]*(model$settings[3])/(model$settings[3]-2)
      tabella_boot[,b] 	= c(bet_boot, sigma2_boot, sigma2_u0_boot, sigma2_u1_boot, covariance_boot)
    } else {
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
  
  if(length(model$theta) > 1){
    row.names(CI) 		= c(names(model$coefficients), "sigma2", "sigma2_intercept", "sigma2_time", "covariance")
  }else{
    row.names(CI) 		= c(names(model$coefficients), "sigma2", "sigma2_intercept") 
  }
  
  colnames(CI) 			= c("lower bound", "upper bound")
  return(CI)
}