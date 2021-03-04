clusjack <- function (model, data, clusterid) {
  res.or <- model
  summ <- summary(model)
  n <- nrow(data) 
  betas <- fixef(model)
  
  if(length(summ$varcor)<2){
    if(length(summ$varcor$id)==1){
      randoms <- summ$varcor$id[1]
    }else{
      randoms <- as.matrix(summ$varcor[[1]])[c(1,2,4)] 
    }
  }else{
    randoms <- c(summ$varcor$id[1],summ$varcor$id.1[1])
  }
  
  
  res.or.coef <- c(fixef(model),sigma(model)^2,randoms)
  
  if(length(randoms)==1){
    names(res.or.coef) <- c(names(fixef(model)), "sigma2", "sigma2_intercept")
  }
  
  if(length(randoms)==2){
    names(res.or.coef) <- c(names(fixef(model)), "sigma2", "sigma2_intercept", "sigma2_time")
  }
  
  if(length(randoms)==3){
    names(res.or.coef) <- c(names(fixef(model)), "sigma2", "sigma2_intercept", "sigma2_time", "covariance")
  }
  
  
  p <- length(betas) + length(randoms) + 1
  coefs <- matrix(NA, nrow = length(unique(clusterid)), ncol = p)
  cluster <- as.character(clusterid)
  clusters <- unique(cluster)
  nc <- length(clusters)
  nc=2
  Obsno <- split(1:n, cluster)
  for (i in 1:nc) {
    obs <- unlist(Obsno[-i])
    modeljack <- update(model, data = data[obs,])
    summjack <- summary(modeljack)
    if(length(randoms)==1){
    jackcoef <- c(fixef(modeljack),sigma(modeljack)^2,summjack$varcor$id[1])
    }
    if(length(randoms)==2){
      jackcoef <- c(fixef(modeljack),sigma(modeljack)^2,summjack$varcor$id[1],summjack$varcor$id.1[1])
    }
    if(length(randoms)==3){
      jackcoef <- c(fixef(modeljack),sigma(modeljack)^2,as.matrix(summjack$varcor[[1]])[c(1,2,4)])
    }
    coefs[i,] <- as.vector(jackcoef)
  }
  colnames(coefs) <- names(res.or.coef)
  uu <- -sweep(coefs,2,colMeans(coefs,na.rm = T), FUN="-")
  acc<-rep(NA,p)
  for(i in 1:p){
    acc[i] <- sum(uu[,i] * uu[,i]* uu[,i],na.rm = T)/(6 * (sum(uu[,i] * uu[,i],na.rm=T))^1.5)
  }
  return(acc)
}




confint_BCa <- function(B, model, data, clusterid, coefs, res.or.coef, p, confint.Zboundaries){
  B_alt <- rep(B,p)
  acc <- clusjack(model,data,clusterid)
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





BCaboot <- function(model, data, clusterid, methodCI, Time = time, B=B, confint.level=.95){
  
  res.or <- model
  summ <- summary(model)
  confint.pboundaries = c((1-confint.level)/2,1-(1-confint.level)/2)
  confint.Zboundaries = qnorm(confint.pboundaries)
  n <- nrow(data) 
  betas <- fixef(model)
  
  if(length(summ$varcor)<2){
    if(length(summ$varcor$id)==1){
      randoms <- summ$varcor$id[1]
    }else{
      randoms <- as.matrix(summ$varcor[[1]])[c(1,2,4)] 
    }
  }else{
    randoms <- c(summ$varcor$id[1],summ$varcor$id.1[1])
  }
    
  
  res.or.coef <- c(fixef(model),sigma(model)^2,randoms)
  
  if(length(randoms)==1){
    names(res.or.coef) <- c(names(fixef(model)), "sigma2", "sigma2_intercept")
  }
  
  if(length(randoms)==2){
    names(res.or.coef) <- c(names(fixef(model)), "sigma2", "sigma2_intercept", "sigma2_time")
  }
  
  if(length(randoms)==3){
    names(res.or.coef) <- c(names(fixef(model)), "sigma2", "sigma2_intercept", "sigma2_time", "covariance")
  }
  
  
  p <- length(betas) + length(randoms) + 1
  
  resampling <- confint.LMM(model = model, Data = data, id = clusterid, Time = time, method = methodCI, B = B, level = .95)
  
  coefs <- resampling$estimation
  colnames(coefs) <- names(res.or.coef) 
  
  ci_BCa <- confint_BCa(B, model, data, clusterid, coefs, res.or.coef, p, confint.Zboundaries)
  
  result <- list(BCa.interval = ci_BCa)
   if(length(dim(result$BCa.interval))>0){
  colnames(result$BCa.interval) <- c("lower bound","upper bound")
  row.names(result$BCa.interval) <- colnames(coefs)
  }
  return(list(result,Percentile=resampling[[2]],bootestimates=coefs,biasBCa=ci_BCa[[2]],acc=ci_BCa[[3]]))
}

