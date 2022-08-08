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
  #matfixSMALL 		  = unname(cbind(rep(1,dim(bdd)[1]),bdd[effetsfix]))
  
  #if(length(var.inter)>0){
  #for (j in length(var.inter)){
  #   matfix = unname(cbind(matfixSMALL,bdd[var.inter[[j]][1]]*bdd[var.inter[[j]][2]]))
  # }
  #}
  
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
  
  result 	            = NULL
  resultr 	          = NULL
  
  for(b in 1:B){
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
                     model.bootr 	= lmer(formulboot, data = bdd,REML=F)                  
                     
                     if(length(model.bootr@optinfo$conv$lme4$messages) == 0){OK = TRUE}
                   }
                   
                   summb 					= summary(model.bootr)
                   bet_boot 				= fixef(model.bootr)
                   sigma_boot 		= as.data.frame(VarCorr(model.bootr))[,4]
                   
                   
                   
                   
                     result 		        = c(bet_boot, sigma_boot) 
                   
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
  


    row.names(CI) <- c(names(fixef(model)),as.data.frame(VarCorr(model.bootr))[,1])
 
  
  colnames(CI) 		  = c("lower bound", "upper bound")
  return(list(estimation=results$effects,CI))
}