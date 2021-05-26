clusjackVarCompRob <- function (model, data, clusterid) {
ls(globalenv())
 
  summ <- summary(model)
  n_obs <- nrow(data) 
  betas <- model$fixef

  randoms <- model$eta
  
  
  res.or.coef <- c(model$fixef,model$eta0,randoms)
  
  
  names(res.or.coef) <- c(names(model$fixef), "sigma2", names(randoms))
  
  p <- length(betas) + length(randoms) + 1
  coefs <- matrix(NA, nrow = length(unique(clusterid)), ncol = p)
  cluster <- as.character(clusterid)
  clusters <- unique(cluster)
  nc <- length(clusters)
  
  Obsno <- split(1:n_obs, cluster)
  
  for (ii in 1:nc) {
    print(ii)
    obs <- unlist(Obsno[-ii])
    print(obs)
    i <- which(model$model$`(groups)`[,2]==ii)
    print(i)
    groups_i <<- groups[-i,]
    print(groups)
    print(groups_i)
	data_i <<- data[obs,]
	print(data_i)
    modeljack <- update(model,data = data_i,groups = groups_i)
    summjack <- summary(modeljack)
    jackcoef <- c(modeljack$fixef,modeljack$eta0,modeljack$eta)
    coefs[ii,] <- as.vector(jackcoef)
  }
  colnames(coefs) <- names(res.or.coef)
  uu <- -sweep(coefs,2,colMeans(coefs,na.rm = T), FUN="-")
  acc<-rep(NA,p)
  for(i in 1:p){
    acc[i] <- sum(uu[,i] * uu[,i]* uu[,i],na.rm = T)/(6 * (sum(uu[,i] * uu[,i],na.rm=T))^1.5)
  }
  return(acc)
}





confint_BCaVarComptob <- function(B, model, data, clusterid, coefs, res.or.coef, p, confint.Zboundaries){
  B_alt <- rep(B,p)
  acc <- clusjackVarCompRob(model,data,clusterid)
  biascorr <- qnorm(colSums(sweep(coefs,2,res.or.coef,"-")<0,na.rm = T)/B_alt) #remplacé B_alt par length(unique(clusterid))?
  tt <- ci_BCa <- matrix(NA, nrow=p, ncol=2)
  ooo <- NA
  for (i in 1:p){
    tt[i,] <- as.vector(pnorm(biascorr[i] + (biascorr[i] + confint.Zboundaries)/(1 - acc[i] * (biascorr[i] + confint.Zboundaries))))
    ooo <- trunc(tt[i,]*B_alt[i])
    tryCatch(ci_BCa[i,]<-sort(coefs[,i])[ooo],error = function(e) { ci_BCa[i,]<-c(NA,NA)})
  }
  return(list(ci_BCa,biascorr,acc))
}





BCaboot_VarCompRob <- function(model, data, clusterid, methodCI, Time = time, B=B, confint.level=.95){
  
  res.or <- model
  summ <- summary(model)
  confint.pboundaries = c((1-confint.level)/2,1-(1-confint.level)/2)
  confint.Zboundaries = qnorm(confint.pboundaries)
  n <- nrow(data) 
  betas <- model$fixef
  

    
      randoms <- model$eta
    
  
  res.or.coef <- c(model$fixef,model$eta0,randoms)
  
  
    names(res.or.coef) <- c(names(model$fixef), "sigma2", names(randoms))
  
  
  
  p <- length(betas) + length(randoms) + 1
  
  resampling <- confint.LMM(model = model, Data = data, id = clusterid, Time = time, method = methodCI, B = B, level = .95)
  
  coefs <- resampling$estimation
  colnames(coefs) <- names(res.or.coef) 
  
  ci_BCa <- confint_BCaVarComptob(B, model, data, clusterid, coefs, res.or.coef, p, confint.Zboundaries)
  
  result <- list(BCa.interval = ci_BCa)
  if(length(dim(result$BCa.interval))>0){
  colnames(result$BCa.interval) <- c("lower bound","upper bound")
  row.names(result$BCa.interval) <- colnames(coefs)
  }
  return(list(result,Percentile=resampling[[2]],bootestimates=coefs,biasBCa=ci_BCa[[2]],acc=ci_BCa[[3]]))
}


