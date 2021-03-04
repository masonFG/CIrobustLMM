

Simul2 = function(NTot,
                  extp1,
                  extRANI1,
                  extRANP1,
                  J,
                  estimator,
                  bootstrap,
                  rseed){
  rs=1000
  if(extp1>0){rs=10000}
  if(extRANI1>0){rs=100000}
  if(extRANP1>0){rs=1000000}
 

#ECHANTILLON

n_suj=NTot #nb de sujet
n_mes=J #nb de mesure max par sujet
set.seed(rseed+rs*J*n_suj)
time<-as.numeric(rep(c(0:(n_mes-1)),n_suj))
id<-sort(as.numeric(rep(1:n_suj,n_mes)))
n_obs = n_suj*n_mes #nb de lignes
#extp1=0 # pourcentage de valeurs extremes résidus
norp=1-extp1 #pourcentage de valeurs non extremes residus
#extRANI1=0 #pourcentage de valeurs extremes ranef intercept
norpRANI=1-extRANI1 #pourcentage de valeurs non extremes ranef intercept
#extRANP1=0 #pourcentage de valeurs extremes ranef pente
norpRANP=1-extRANP1 #pourcentage de valeurs non extremes ranef pente
al = 2 #nb d'effet aleatoire par sujet
varRES_n = 400#variance erreurs normale
varRES_e = 0.25#variance erreurs extremes 
moyRES_e = -4*sqrt(varRES_n)#moyenne des erreurs etremes
matran_n = matrix(c(790,-8.5,-8.5,40),2,2)#matrice variance covariance normale
matran_e = matrix(c(7.9,-.085,-.085,.4),2,2)#matrice variance covariance extreme
moyran_Ie = c(-4*sqrt(790),0)#vecteur moyennes des random intercept
moyran_Pe = c(0,-4*sqrt(40))#vecteur moyennes des random pente

#BDD
pourcentN=min(norpRANI,norpRANP) #pourcentage d'individu non contamines






#set.seed(rseed)
reg<-mvrnorm(n=round(norp*n_obs),mu=0,Sigma=varRES_n)
indexRES<-rep(0,length(reg))
res<-data.frame(cbind(reg,indexRES))
colnames(res)<-c("r","indexRES")

if (extp1>0){
  #set.seed(rseed*10000)
  ext<-mvrnorm(n=round(extp1*n_obs),mu=moyRES_e,Sigma=varRES_e)  #residus extremes, normale translate
  indexRES<-rep(1,length(ext))
  ext<-data.frame(cbind(ext,indexRES))
  colnames(ext)<-c("r","indexRES")
  r_inter<-rbind(res,ext)
  
  #set.seed(rseed*20000)
  indexRES<-sample(c(1:n_obs),n_obs,replace=F)
  resi<-NULL
  R<-NULL
  for (i in indexRES){
    res<-r_inter[i,] 
    R<-rbind(R,res)
  }
}else{
  R<-cbind(res,id)
  colnames(R)<-c("r","indexRES","id")
}



##### EFFETS ALEATOIRES #####


#AVEC CORRELATION ENTRE INTERCEPT ET PENTE
#normale
#set.seed(rseed+10000)
u_n<-as.data.frame(mvrnorm(round(pourcentN*n_suj),mu=c(0,0),Sigma=matran_n))
u_n$indexI<-0
u_n$indexP<-0

if(extRANI1+extRANP1==0){
  U<-cbind(u_n,unique(id))
  colnames(U)<-c("u0","u1","indexRANI","indexRANP","id") 
}else{
  if(extRANI1>0){
    #extremeI
    #set.seed(rseed+20000)
    u_Ie<-as.data.frame(mvrnorm(round(extRANI1*n_suj),mu=moyran_Ie,Sigma=matran_e))
    u_Ie$indexI<-1
    u_Ie$indexP<-0
    u_inter<-data.frame(rbind(u_n,u_Ie))
    
    #set.seed(rseed+30000)
    indexRAN<-sample(c(1:n_suj),n_suj,replace=F)
    UU<-NULL
    U<-NULL
    for (i in indexRAN){
      UU<-u_inter[i,] 
      U<-rbind(UU,U)
    }
    U<-cbind(U,unique(id))
    colnames(U)<-c("u0","u1","indexRANI","indexRANP","id") 
  }else{
    #extremeP
    #set.seed(rseed+20000)
    u_Pe<-as.data.frame(mvrnorm(round(extRANP1*n_suj),mu=moyran_Pe,Sigma=matran_e))
    u_Pe$indexI<-0
    u_Pe$indexP<-1
    u_inter<-data.frame(rbind(u_n,u_Pe))
    
    #set.seed(rseed+30000)
    indexRAN<-sample(c(1:n_suj),n_suj,replace=F)
    UU<-NULL
    U<-NULL
    for (i in indexRAN){
      UU<-u_inter[i,] 
      U<-rbind(UU,U)
    }
    U<-cbind(U,unique(id))
    colnames(U)<-c("u0","u1","indexRANI","indexRANP","id")
  }
}

bdd<-cbind(id,time)
bdd<-merge(bdd,U)
bdd<-cbind(bdd,R)
bdd1<-bdd[,-9]
row.names(bdd1)<-c(1:(n_mes*n_suj))





#BDD FINALE
bdd <- bdd1



bdd$y=250 + 10 * bdd$time + bdd$u0 + bdd$u1 * bdd$time + bdd$r
bdd$obs<-c(1:dim(bdd)[1])


y = bdd$y #VD
time = bdd$time # pente
id = bdd$id # grouping variable
Testeffetfix =F         #autre effet fixe T si oui sinon F
if (Testeffetfix==T){
  effetfixe<-bdd[,c(2,4)]        #selectionner les variables pour les effets fixes
  bddMC = as.data.frame(cbind(y,id,time,effetfixe))
} else{
  bddMC = as.data.frame(cbind(y,id,time)) #base de donnée
} 











# Identify dataset, time and participant variables
Dataset = bdd
time <- bdd$time
Time <<- bdd$time
participant = bdd$id

###############
# LMM estimations

# 4a) Estimation with varComprob() (see corresponding helpfile for more details)

# Build the argument "groups" of the varComprob() function
n <<- length(unique(participant)) # the number of participants
J <<- length(unique(time)) # the number of repeated observations per participant
groups <<- cbind(rep(1:J, each=n),rep((1:n), J)) # a numeric matrix with two columns used to group the observations according to participant.

# Build the argument "varcov" of the varComprob() function
z1 = rep(1, J) #Value for intercept (=1) for the J observations by clusters
z2 = unique(time) # Value for the time variable

K <- list() # the "varcov" object
K[[1]] <- tcrossprod(z1,z1) # Matrix for intercept
K[[2]] <- tcrossprod(z2,z2) # Matrix for time variable
K[[3]] <- tcrossprod(z1,z2) + tcrossprod(z2,z1) # Matrix of interaction Intercept by time variable
names(K) = c("sigma2_Intercept", "sigma2_Time", "Covariance")
K <<- K

# Define the formula of the two nested models
model.formula <<- y ~ 1 + time
model.formula0 = y ~ 1 +time 

# Estimation with S-estimator
if(estimator == "S"){
model <<- varComprob(model.formula, groups = groups, data = Dataset, varcov = K, control = varComprob.control(lower = c(0,0,-Inf), method = "S", psi = "rocke")) 

if(bootstrap == "wild"){
  if(test == TRUE){
model.S0 = varComprob(model.formula0, groups = groups, data = Dataset, varcov = K, control = varComprob.control(lower = c(0,0,-Inf), method = "S", psi = "rocke")) 

resS<-summary(model)
fixedEffects_S <- model$fixef
randomEffects_S <- model$eta
sigma2_S <- model$eta0
estimates <- c(fixedEffects_S,sigma2_S,randomEffects_S)
inf = model$eta - sqrt(diag(model$vcov.eta))*qnorm(.95+0.05/2)
sup = model$eta + sqrt(diag(model$vcov.eta))*qnorm(.95+0.05/2)

ranTestWaldns_S = inf < 0 & 0 < sup 
pWald_S <- c(resS$zTable[,4])
TestWald_S <- pWald_S>.05
TestWald_S <- c(TestWald_S,NA,ranTestWaldns_S)

wildS <- TestFixef(model = model, model0 = model.S0, Data = Dataset, id = participant, Time = time, method = "wild", B = 5000, level = .95)
SUMMARY <- list(model,model.S0,estimates,cbind(inf,sup),TestWald_S,ranTestWaldns_S,wildS) 
names(SUMMARY)<-c("model","model0","estimates","CIWald","pWaldFixed","pWaldRandom","pWILD_interaction")
 }else{
CIwildS <- BCaboot_VarCompRob(model = model, data = Dataset, clusterid = participant, methodCI = "wild", B = 5000, confint.level = .95)
SUMMARY<-list(model,CIwildS)
names(SUMMARY)<-c("model","CIwild")
}
}

#Testwildns_S = wildS[,1] < 0 & 0 < wildS[,2]
if(bootstrap == "parametric"){
  if(test == TRUE){
    model.S0 = varComprob(model.formula0, groups = groups, data = Dataset, varcov = K, control = varComprob.control(lower = c(0,0,-Inf), method = "S", psi = "rocke")) 
paramS <- TestFixef(model = model, model0 = model.S0, Data = Dataset, id = participant, Time = time, method = "parametric", B = 5000, level = .95)
SUMMARY <-  list(model,model.S0,paramS)
names(SUMMARY)<-c("model","model0","pPARAM_interaction")
}else{
  CIparamS <- BCaboot_VarCompRob(model = model, data = Dataset, clusterid = participant, methodCI = "parametric", B = 5000, confint.level = .95)
  SUMMARY <-  list(model,CIparamS)
  names(SUMMARY)<-c("model","CIparam")
  #Testparamns_S = paramS[,1] < 0 & 0 < paramS[,2]
  }
}
}
# Estimation with composite-TAU estimator
if(estimator == "cTAU"){
model  <<- varComprob(model.formula, groups = groups, data = Dataset, varcov = K, control = varComprob.control(lower = c(0, 0, -Inf))) 

if(bootstrap == "wild"){
  
  if(test == TRUE){
model.cTAU0 = varComprob(model.formula0, groups = groups, data = Dataset, varcov = K, control = varComprob.control(lower = c(0, 0, -Inf))) 

rescTAU<-summary(model)
fixedEffects_cTAU <- model$fixef
randomEffects_cTAU <- model$eta
sigma2_cTAU <- model$eta0
estimates <- c(fixedEffects_cTAU,sigma2_cTAU,randomEffects_cTAU)

inf = model$eta - sqrt(diag(model$vcov.eta))*qnorm(.95+0.05/2)
sup = model$eta + sqrt(diag(model$vcov.eta))*qnorm(.95+0.05/2)

ranTestWaldns_cTAU = inf < 0 & 0 < sup 
pWald_cTAU <- c(rescTAU$zTable[,4])
TestWald_cTAU <- pWald_cTAU>.05
TestWald_cTAU <- c(TestWald_cTAU,NA,ranTestWaldns_cTAU)

wildcTAU <- TestFixef(model = model, model0 = model.cTAU0, Data = Dataset, id = participant, Time = time, method = "wild", B = 5000, level = .95)
SUMMARY <- list(model,model.cTAU0,estimates,cbind(inf,sup),TestWald_cTAU,ranTestWaldns_cTAU,wildcTAU) 
names(SUMMARY)<-c("model","model0","estimates","CIWald","pWaldFixed","pWaldRandom","pWILD_interaction")
}else{
    CIwildcTAU <- BCaboot_VarCompRob(model = model, data = Dataset, clusterid = participant, methodCI = "wild", B = 5000, confint.level = .95)
#Testwildns_cTAU = wildcTAU[,1] < 0 & 0 < wildcTAU[,2]
  SUMMARY <- list(model,CIwildcTAU)
  names(SUMMARY)<-c("model","CIwild")
    }

  }

if(bootstrap == "parametric"){
  
  if(test == TRUE){
    model.cTAU0 = varComprob(model.formula0, groups = groups, data = Dataset, varcov = K, control = varComprob.control(lower = c(0, 0, -Inf))) 
paramcTAU <- TestFixef(model = model, model0 = model.cTAU0, Data = Dataset, id = participant, Time = time, method = "parametric", B = 5000, level = .95)
SUMMARY <- list(model,model.cTAU0,paramcTAU)  
names(SUMMARY)<-c("model","model0","pPARAM_interaction")
}else{
CIparamcTAU <- BCaboot_VarCompRob(model = model, data = Dataset, clusterid = participant, methodCI = "parametric", B = 5000, confint.level = .95)
SUMMARY <- list(model,CIparamcTAU)
names(SUMMARY)<-c("model","CIparam")
#Testparamns_cTAU = paramcTAU[,1] < 0 & 0 < paramcTAU[,2]
  }
  
}

}
# 4b) Estimation with rlmer() (see corresponding helpfile for more details)

# Estimation with SMDM
if(estimator == "SMDM"){
model.SMDM = rlmer(y ~ 1 +  time + (time|id), data = Dataset, rho.sigma.e = psi2propII(smoothPsi, k = 2.28), rho.sigma.b = chgDefaults(smoothPsi, k = 5.11, s = 10))                   

if(bootstrap == "wild"){
  if(test == TRUE){
model.SMDM0 = rlmer(y ~ 1 + time + (time|id), data = Dataset, rho.sigma.e = psi2propII(smoothPsi, k = 2.28), rho.sigma.b = chgDefaults(smoothPsi, k = 5.11, s = 10))                   

summ=summary(model.SMDM)
fixedEffects_SMDM <- fixef(model.SMDM)
randomEffects_SMDM <- c(VarCorr(model.SMDM)$id[1,1], VarCorr(model.SMDM)$id[2,2],VarCorr(model.SMDM)$id[2,1])
sigma2_SMDM <- sigma(model.SMDM)^2
estimates <- c(fixedEffects_SMDM,sigma2_SMDM,randomEffects_SMDM)

inf = c(fixef(model.SMDM) - summ$coefficients[,2]*qnorm(.95+.05/2))
sup = c(fixef(model.SMDM) + summ$coefficients[,2]*qnorm(.95+.05/2))
TestWaldns_SMDM = inf < 0 & 0 < sup
TestWaldns_SMDM = c(TestWaldns_SMDM,NA,NA,NA,NA)

wildSMDM <- TestFixef(model = model.SMDM, model0 = model.SMDM0, Data = Dataset, id = participant, Time = time, method = "wild", B = 5000, level = .95)
SUMMARY <- list(model.SMDM,model.SMDM0,estimates,cbind(inf,sup),TestWaldns_SMDM,wildSMDM)
names(SUMMARY)<-c("model","model0","estimates","CIWald","pWaldFixed","pWILD_interaction")
}else{
CIwildSMDM <- BCaboot(model = model.SMDM, data = Dataset, clusterid = participant, methodCI = "wild", Time = time, B = 5000, confint.level = .95)
SUMMARY <- list(model.SMDM,CIwildSMDM)
names(SUMMARY)<-c("model","CIwild")
#Testwildns_SMDM = wildSMDM[,1] < 0 & 0 < wildSMDM[,2]
}
}

if(bootstrap == "parametric"){
  if(test == TRUE){
    model.SMDM0 = rlmer(y ~ 1 + time + (time|id), data = Dataset, rho.sigma.e = psi2propII(smoothPsi, k = 2.28), rho.sigma.b = chgDefaults(smoothPsi, k = 5.11, s = 10))                   
    paramSMDM <- TestFixef(model = model.SMDM, model0 = model.SMDM0, Data = Dataset, id = participant, Time = time, method = "parametric", B = 5000, level = .95)
    SUMMARY <- list(model.SMDM,model.SMDM0,paramSMDM) 
    names(SUMMARY)<-c("model","model0","pPARAM_interaction")
}else{
CIparamSMDM <- BCaboot(model = model.SMDM, data = Dataset, clusterid = participant, methodCI = "parametric", Time = time, B = 5000, confint.level = .95)
SUMMARY <- list(model.SMDM,CIparamSMDM)  
names(SUMMARY)<-c("model","CIparam")
#Testparamns_SMDM = paramcSMDM[,1] < 0 & 0 < paramcSMDM[,2]
}
}
}
# 4c) Estimation with heavyLme() (see corresponding helpfile for more details)

# Estimation with the multivarite t-ML
if(estimator == "tML"){
model.tML = heavyLme(y ~ 1 +  time, random = ~ 1 + time, groups = ~ id, data = Dataset) 

if(bootstrap == "wild"){
  if(test == TRUE){
model.tML0 = heavyLme(y ~ 1  + time, random = ~ 1 + time, groups = ~ id, data = Dataset) 

summ = summary(model.tML)
fixedEffects_tML 			= model.tML$coefficients
randomEffects_tML 	= c(summ$theta[1,1]*(model.tML$settings[3])/(model.tML$settings[3]-2),
                       summ$theta[2,2]*(model.tML$settings[3])/(model.tML$settings[3]-2),
                       summ$theta[1,2]*(model.tML$settings[3])/(model.tML$settings[3]-2))
sigma2_tML 		= summ$scale*(model.tML$settings[3])/(model.tML$settings[3]-2)
estimates <- c(fixedEffects_tML,sigma2_tML,randomEffects_tML)
inf = c(coefficients(model.tML) - summ$coefficients[,2]*qnorm(.975))
sup = c(coefficients(model.tML) + summ$coefficients[,2]*qnorm(.975))
TestWaldns_tML = inf < 0 & 0 < sup
TestWaldns_tML = c(TestWaldns_tML,NA,NA,NA,NA)

wildtML <- TestFixef(model = model.tML, model0 = model.tML0, Data = Dataset, id = participant, Time = time, method = "wild", B = 5000, level = .95)
SUMMARY <- list(model.tML,model.tML0,estimates,cbind(inf,sup),TestWaldns_tML,wildtML) 
names(SUMMARY)<-c("model","model0","estimates","CIWald","pWaldFixed","pWILD_interaction")
}else{
CIwildtML <- BCaboot_Heavy(model = model.tML, data = Dataset, clusterid = participant, methodCI = "wild", Time = time, B = 5000, confint.level = .95)
SUMMARY <- list(model.tML,CIwildtML) 
names(SUMMARY)<-c("model","CIwild")
#Testwildns_tML = wildtML[,1] < 0 & 0 < wildtML[,2]
  }
}

if(bootstrap == "parametric"){
  if(test == TRUE){
    model.tML0 = heavyLme(y ~ 1 + time, random = ~ 1 + time, groups = ~ id, data = Dataset) 
paramtML <- TestFixef(model = model.tML, model0 = model.tML0, Data = Dataset, id = participant, Time = time, method = "parametric", B = 5000, level = .95)
SUMMARY <- list(model.tML,model.tML0,paramtML)
names(SUMMARY)<-c("model","model0","pPARAM_interaction")
}else{
CIparamtML <- BCaboot_Heavy(model = model.tML, data = Dataset, clusterid = participant, methodCI = "parametric", Time = time, B = 5000, confint.level = .95)
SUMMARY <- list(model.tML,CIparamtML)
names(SUMMARY)<-c("model","CIparam")
#Testparamns_tML = paramtML[,1] < 0 & 0 < paramtML[,2]
}
}
}
# 4d) Estimation with lmer() (see corresponding helpfile for more details)

# Estimation with ML
if(estimator == "ML"){
model.ML = lmer(y ~ 1 +  time + (time|id), data = Dataset, REML = F)

if(bootstrap == "wild"){
  if(test == TRUE){
model.ML0 = lmer(y ~ 1 + time + (time|id), data = Dataset, REML = F)

summ=summary(model.ML)
fixedEffects_ML <- fixef(model.ML)
randomEffects_ML <- c(VarCorr(model.ML)$id[1,1], VarCorr(model.ML)$id[2,2],VarCorr(model.ML)$id[2,1])
sigma2_ML <- sigma(model.ML)^2
estimates <- c(fixedEffects_ML,sigma2_ML,randomEffects_ML)

#inf = c(fixef(model.ML) - summ$coefficients[,2]*qnorm(.95+.05/2))
#sup = c(fixef(model.ML) + summ$coefficients[,2]*qnorm(.95+.05/2))
#TestWaldns_ML = inf < 0 & 0 < sup
#TestWaldns_ML = c(TestWaldns_ML,NA,NA,NA,NA)

#pSatterth_ML <- c(summ$coefficients[,5])

#Profile <- confint(model.ML)
#TestProfilens_ML = Profile[,1] < 0 & 0 < Profile[,2]
#TestProfilens_ML = TestProfilens_ML[c(5,6,7,8,4,1,3,2)]


#model.MLgroup = lmer(y ~ 1 + time+group:time + (time|id), data = Dataset, REML = F)
#model.MLtime = lmer(y ~ 1 + group+group:time + (time|id), data = Dataset, REML = F)
#model.MLinter = lmer(y ~ 1 + group+time + (time|id), data = Dataset, REML = F)
#model.MLI = lmer(y ~ 0 + group*time + (time|id), data = Dataset, REML = F)
#model.MLranS = lmer(y ~ 1 + group*time + (1|id), data = Dataset, REML = F)
#model.MLcov = lmer(y ~ 1 + group*time + (time||id), data = Dataset, REML = F)


#LRT_intercept<-anova(model.ML,model.MLI)
#LRT_group<-anova(model.ML,model.MLgroup)
#LRT_time<-anova(model.ML,model.MLtime)
#LRT_interaction<-anova(model.ML,model.MLinter)
#LRT_ranI<-ranova(model.MLranS)
#LRT_ranS<-anova(model.MLranS,model.MLcov)
#LRT_covariance<-anova(model.ML,model.MLcov)

#pLRT <- c(LRT_intercept[2,8],LRT_group[2,8],LRT_time[2,8],LRT_interaction[2,8],NA,LRT_ranI[2,6],LRT_ranS[2,8],LRT_covariance[2,8])

wildML <- TestFixef(model = model.ML, model0 = model.ML0, Data = Dataset, id = participant, Time = time, method = "wild", B = 5000, level = .95)
SUMMARY <- list(model.ML,model.ML0,estimates,wildML) 
names(SUMMARY)<-c("model","model0","estimates","pWILD_interaction")
 }else{
CIwildML <- BCaboot(model = model.ML, data = Dataset, clusterid = participant, methodCI = "wild", Time = time, B = 5000, confint.level = .95)
SUMMARY <- list(model.ML,CIwildML)
names(SUMMARY)<-c("model","CIwild")
#Testwildns_ML = wildML[,1] < 0 & 0 < wildML[,2]
  }
}

if(bootstrap == "parametric"){
  if(test == TRUE){
    model.ML0 = lmer(y ~ 1 + time + (time|id), data = Dataset, REML = F)
paramML <- TestFixef(model = model.ML, model0 = model.ML0, Data = Dataset, id = participant, Time = time, method = "parametric", B = 5000, level = .95)
SUMMARY <- list(model.ML,model.ML0,paramML) 
names(SUMMARY)<-c("model","model0","pPARAM_interaction")
}else{
CIparamML <- BCaboot(model = model.ML, data = Dataset, clusterid = participant, methodCI = "parametric", Time = time, B = 5000, confint.level = .95)
SUMMARY <- list(model.ML,CIparamML) 
names(SUMMARY)<-c("model","CIparam")

#Testparamns_ML = paramcML[,1] < 0 & 0 < paramcML[,2]
}
}
}
return(SUMMARY)
}




