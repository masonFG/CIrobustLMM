# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: March, 2021
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the function of the wild bootstrap on REML
# ---------------------------------------------------------------------

wild_lmer <- function(model, B, level){
  
  if (length(model@theta) > 1){                              
    P 	= 2
  } else{
    P 	= 1
  }
  
  #interaction
  inter = grep(pattern = ":" , names(fixef(model)), value = TRUE, fixed = TRUE)
  if(!is.null(inter)){
    var.inter = strsplit(inter, split = ":" , fixed = TRUE)
  }
  
  bdd 			  = model@frame
  y 			    = model@resp$y
  effetsfix 	= attr(attr(model@frame,"terms"),"varnames.fixed")[-1] 
  matfix 		  = unname(cbind(rep(1,dim(bdd)[1]),bdd[effetsfix]))
  
  if(length(var.inter)>0){
    for (j in length(var.inter)){
      matfix = unname(cbind(matfix,bdd[var.inter[[j]][1]]*bdd[var.inter[[j]][2]]))
    }
  }
  
  nomID 		  = names(summary(model)$ngrps)[1]
  id 			    = bdd[nomID]
  summ 			  = summary(model)
  n 			    = nrow(bdd)
  bet 			  = as.vector(unname(fixef(model)))
  sigma2 		  = sigma(model)^2
  sigma2_u0 	= as.matrix(summ$varcor[[1]])[1]
  
  if (length(model@theta) > 1) {
    sigma2_u1 	= as.matrix(summ$varcor[[1]])[4]
    covariance 	= as.matrix(summ$varcor[[1]])[2]
    est 		    = c(bet, sigma2, sigma2_u0, sigma2_u1, covariance)
    names(est) 	= c("intercept", "time", "sigma2", "sigma2_u0", "sigma2_u1", "covariance")
  } else {
    est 		    = c(bet, sigma2, sigma2_u0)
    names(est) 	= c("intercept", "time", "sigma2", "sigma2_u0")  
  }
  
  v1 			= -(sqrt(5) - 1) /2 
  v2 			= (sqrt(5) + 1) /2 
  p1 			= (sqrt(5) + 1) /(2 * sqrt(5)) 
  p2 			= 1 - p1 
  rand 		= table(id)
  TT 			= length(rand)
  nt 			= unname(rand)
  X 			= matrix(1,n,1)
  reg 		= attr(terms(formula(model)), "term.labels")
  X 			= model.matrix(model, data = as.data.frame(bdd[,reg]))
  XX 			= as.matrix(unname(X))                                                 
  Xt 			= as.matrix(unname(vector("list", TT)) )                                
  yt 			= as.matrix(vector("list", TT))
  rt_hat 	= as.matrix(vector("list", TT))
  Pt 			= as.matrix(vector("list", TT))
  Int 		= as.matrix(vector("list", TT))
  tXX 		= as.matrix(unname(solve(t(X)%*%X)))          
  
  for (i in 1:TT) {
    Xt[[i]] 	  = XX[1:nt[i],]                
    XX 			    = XX[-(1:nt[i]),]             
    yt[[i]] 	  = y[1:nt[i]]                  
    y 			    = y[-(1:nt[i])]               
    rt_hat[[i]] = yt[[i]] - Xt[[i]]%*%bet   
    Pt[[i]] 	  = Xt[[i]]%*%tXX%*%t(Xt[[i]])  
    Int[[i]] 	  = diag(1,nt[i],nt[i])         
  }
  
  rt_b 			    = vector("list", TT)
  yt_b 			    = vector("list", TT)
  
  # Bootstrap scheme
  # Calculate the number of cores
  no_cores 		  = detectCores() - 1
  
  # Initiate cluster
  cl 			      = makeCluster(no_cores)
  registerDoParallel(cl)
  
  tAB = foreach (b = 1:B,
                 .combine = "rbind",
                 .packages = c("lme4", "MASS")) %dopar% {
                   OK = FALSE
                   
                   while(!OK){
                     wt = sample(c(v1, v2), TT, replace = TRUE, prob = c(p1, p2))  
                     
                     for (tt in 1:TT) {
                       rt_b[[tt]] = sqrt(diag(ginv(Int[[tt]]-Pt[[tt]])))*rt_hat[[tt]]*wt[tt] 
                       yt_b[[tt]] = Xt[[tt]]%*%bet + rt_b[[tt]]
                     }
                     
                     yboot 				= unlist(yt_b)
                     bdd$yboot 		= yboot
                     formulrest 		= as.character(formula(model))[3]
                     formulboot 		= paste("yboot ~", formulrest)
                     model.bootr 	= lmer(formulboot, data = bdd)                  
                     
                     if(length(model.bootr@optinfo$conv$lme4$messages) == 0){OK = TRUE}
                   }
                   
                   summb 					= summary(model.bootr)
                   bet_boot 				= fixef(model.bootr)
                   sigma2_boot 		= sigma(model.bootr)^2
                   sigma2_u0_boot 	= as.matrix(summb$varcor[[1]])[1]
                   
                   if (length(model@theta) > 1) {
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
  
  if(length(model@theta) > 1){
    row.names(CI) 	= c(names(fixef(model)), "sigma2", "sigma2_intercept", "sigma2_time", "covariance")
  }else{
    row.names(CI) 	= c(names(fixef(model)), "sigma2", "sigma2intercept") 
  }
  
  colnames(CI) 		  = c("lower bound", "upper bound")
  return(CI)
}
