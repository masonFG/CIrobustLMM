# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the function of the wild bootstrap on lmer and rlmer
# ---------------------------------------------------------------------

wild_lmer <- function(model, model0, B, level){
  
  if (length(model@theta) > 1){                              
    P 	= 2
  } else{
    P 	= 1
  }
  
  bdd 			  = model@frame
  y 			    = model@resp$y
  matfixSMALL 		  = model.matrix(model0)
  matfix       		  = model.matrix(model)
  index_fixef = which(colnames(matfix)!=colnames(matfixSMALL))
  matfix 	  	= as.matrix(matfix)
  matfixSMALL = as.matrix(matfixSMALL)
  
  nomID 		  = names(summary(model)$ngrps)[1]
  id 			    = bdd[nomID]
  n 			    = nrow(bdd)
  bet 			  = as.vector(unname(fixef(model)))
  betSMALL    = as.vector(unname(fixef(model0)))
  b_tested    = fixef(model)[index_fixef]
  
  v1 			= -(sqrt(5) - 1) /2 
  v2 			= (sqrt(5) + 1) /2 
  p1 			= (sqrt(5) + 1) /(2 * sqrt(5)) 
  p2 			= 1 - p1 
  rand 		= table(id)
  TT 			= length(rand)
  nt 			= unname(rand)
  X 			= as.matrix(matfixSMALL)
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
    rt_hat[[i]] = yt[[i]] - Xt[[i]]%*%betSMALL   
    Pt[[i]] 	  = Xt[[i]]%*%tXX%*%t(Xt[[i]])  
    Int[[i]] 	  = diag(1,nt[i],nt[i])         
  }
  
  rt_b 			    = vector("list", TT)
  yt_b 			    = vector("list", TT)
  
  result 	            = NULL
  resultr 	          = NULL
  
  for(b in 1:B){
    OK = FALSE
                   
                   while(!OK){
                     wt = sample(c(v1, v2), TT, replace = TRUE, prob = c(p1, p2))  
                     
                     for (tt in 1:TT) {
                       rt_b[[tt]] = sqrt(diag(ginv(Int[[tt]]-Pt[[tt]])))*rt_hat[[tt]]*wt[tt] 
                       yt_b[[tt]] = Xt[[tt]]%*%betSMALL + rt_b[[tt]]
                     }
                     
                     yboot 				= unlist(yt_b)
                     bdd$yboot 		= yboot
                     formulrest 		= as.character(formula(model))[3]
                     formulboot 		= paste("yboot ~", formulrest)
                     model.bootr 	= lmer(formulboot, data = bdd,REML=F)                  
                     
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
                   
                   resultr<-rbind(resultr,result)
                 }
  
 
  
  estim                 = resultr
  colnames(estim)       = c(names(fixef(model)), "sigma2", "sigma2_intercept", "sigma2_time", "covariance")
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
