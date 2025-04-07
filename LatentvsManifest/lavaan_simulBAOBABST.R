 library(lavaan)
 library(MASS)
 library(semTools)
 library(semPlot)
 library(dplyr)
 library(tidyverse)
 library(ggplot2)
 library(pipeR)
 library(psych)
 library(lmerTest)

#Example with 25 samples of 100 participants, 12 items with mean and variance change respectively of 0.3 and 0.01
N=100; items=12; mD0=0; mD1=0.3; varD=0.01;covar=0.25*sqrt(varD)*1;B=25;

SimulmeanScore= function(N,B,items,mD0,mD1,varD,covar){
n=N/2 #sample size
#LOADINGS DEFINITION
if(items<7){
  increment = 0.2/(items-1) #increment in communality value for models "low" and "high"
  comm_low = seq(0.2,0.4,by=increment) #communality for model low
  comm_high = seq(0.6,0.8,by=increment) #communality for congeneric high
}else{
  comm_low = rep(seq(0.2,0.4,by=0.04),each=2) #communality for model low
  comm_high = rep(seq(0.6,0.8,by=0.04),each=2) #communality for congeneric high
}

increment_wide = 0.6/(items-1) #increment in communality value for model wide
comm_wide = seq(0.2,0.8,by=increment_wide) #communality for model low

loads_low <- sqrt(comm_low)
loads_high <- sqrt(comm_high)
loads_wide <- sqrt(comm_wide)

change_mean = mD1 #latent change mean
change_var = varD #latent change var

verbal_mean = 0 #factor mean at time 0
verbal_var = 1 #factor var at time 0
cov_changeVerbal = covar #covariance between change en factor scores at time 0
MatrixChange <- rbind(c(change_var,cov_changeVerbal),c(cov_changeVerbal,verbal_var)) #as matrix
colnames(MatrixChange) <- c("change","factor")
rownames(MatrixChange) <- c("change","factor")


#residual variance definition
var_resid_low = 1-loads_low^2 #residuals variance for model low
var_resid_high = 1-loads_high^2 #residuals variance for model high
var_resid_wide = 1-loads_wide^2 #residuals variance for model wide



#residual covariance 
rescov_low = rep(var_resid_low,items)
rescov_high = rep(var_resid_high,items)
rescov_wide = rep(var_resid_wide,items)

#residuals matrix
sigma_res_low <-diag(c(rescov_low,rescov_low),nrow=2*items,ncol=2*items)
sigma_res_high <-diag(c(rescov_high,rescov_high),nrow=2*items,ncol=2*items)
sigma_res_wide <-diag(c(rescov_wide,rescov_wide),nrow=2*items,ncol=2*items)

#items names
manifeste_t0 = paste("x",1:items,"_t0",sep = "")
manifeste_t1 = paste("x",1:items,"_t1",sep = "")
colnames(sigma_res_low)<-c(manifeste_t0,manifeste_t1)
rownames(sigma_res_low)<-c(manifeste_t0,manifeste_t1)
colnames(sigma_res_high)<-c(manifeste_t0,manifeste_t1)
rownames(sigma_res_high)<-c(manifeste_t0,manifeste_t1)
colnames(sigma_res_wide)<-c(manifeste_t0,manifeste_t1)
rownames(sigma_res_wide)<-c(manifeste_t0,manifeste_t1)

#intercept of manifest variables
intercepts_xi <- rep(0,2*items)

#model to be tested
loads_test = paste("a",2:(items),sep = "")
loads_testL = c("0.447213595499958",loads_test)
loads_testH = c("0.774596669241483",loads_test)

regt0_testL = "factor_t0 =~ "
regt1_testL = "factor_t1 =~ "
regt0_testH = "factor_t0 =~ "
regt1_testH = "factor_t1 =~ "

for(i in 1:items){
  if(i<items){
    regt0_testL = paste0(regt0_testL,loads_testL[i],"*",manifeste_t0[i]," + ")
    regt1_testL = paste0(regt1_testL,loads_testL[i],"*",manifeste_t1[i]," + ")
    regt0_testH = paste0(regt0_testH,loads_testH[i],"*",manifeste_t0[i]," + ")
    regt1_testH = paste0(regt1_testH,loads_testH[i],"*",manifeste_t1[i]," + ")
  }else{
    regt0_testL = paste0(regt0_testL,loads_testL[i],"*",manifeste_t0[i])
    regt1_testL = paste0(regt1_testL,loads_testL[i],"*",manifeste_t1[i])
    regt0_testH = paste0(regt0_testH,loads_testH[i],"*",manifeste_t0[i])
    regt1_testH = paste0(regt1_testH,loads_testH[i],"*",manifeste_t1[i])
  }
}


reg_testL = paste0(regt0_testL,"\n",regt1_testL,sep=" ")
reg_testH = paste0(regt0_testH,"\n",regt1_testH,sep=" ")

script_constant_test <-'
Delta =~ 1*factor_t1
         #Regression
         factor_t1 ~ 1*factor_t0
         #Moyennes
         factor_t0 ~ 1
        factor_t1 ~ 0*1
        Delta ~ c(mu_cont,mu_treat)*1
        #Variances facteur
        factor_t0 ~~ factor_t0
        Delta ~~ Delta
        factor_t0 ~~ Delta
        #interaction
        interaction:= (mu_treat-mu_cont)
'
script_constant_testMOY <-'
Delta =~ 1*factor_t1
         #Regression
         factor_t1 ~ 1*factor_t0
         #Moyennes
         factor_t0 ~ 1
        factor_t1 ~ 0*1
        Delta ~ c(mD,mD)*1
        #Variances facteur
        factor_t0 ~~ factor_t0
        Delta ~~ Delta
        factor_t0 ~~ Delta
'
#covariance error script
cov_err_test = paste0(manifeste_t0,"~~",manifeste_t1)
cov_err_test_char = paste0(cov_err_test,collapse = "\n")
#residual variance script
var_test_t0 = paste0(manifeste_t0,"~~",manifeste_t0)
var_test_t0_char = paste0(var_test_t0,collapse = "\n")
var_test_t1 = paste0(manifeste_t1,"~~",manifeste_t1)
var_test_t1_char = paste0(var_test_t1,collapse = "\n")

#SCRIPT TO DEFINE INTERCEPT = 0
intercept_t0 = paste0(manifeste_t0,"~0*1")
intercept_t0_char = paste0(intercept_t0,collapse = "\n")
intercept_t1 = paste0(manifeste_t1,"~0*1")
intercept_t1_char = paste0(intercept_t1,collapse = "\n")

#THE COMPLETE SCRIPT FOR THE MODEL "TEST LOW"
script_testL <- paste0(reg_testL,
                       script_constant_test,
                       cov_err_test_char,"\n",
                       var_test_t0_char,"\n",
                       var_test_t1_char,"\n",
                       intercept_t0_char,"\n",
                       intercept_t1_char,sep=" ")

#THE COMPLETE SCRIPT FOR THE MODEL "TEST HIGH"
script_testH <- paste0(reg_testH,
                       script_constant_test,
                       cov_err_test_char,"\n",
                       var_test_t0_char,"\n",
                       var_test_t1_char,"\n",
                       intercept_t0_char,"\n",
                       intercept_t1_char,sep=" ")

#THE COMPLETE SCRIPT FOR THE MODEL "TEST CHANGE LOW"
script_testLMOY <- paste0(reg_testL,
                          script_constant_testMOY,
                          cov_err_test_char,"\n",
                          var_test_t0_char,"\n",
                          var_test_t1_char,"\n",
                          intercept_t0_char,"\n",
                          intercept_t1_char,sep=" ")

#THE COMPLETE SCRIPT FOR THE MODEL "TEST CHANGE HIGH"
script_testHMOY <- paste0(reg_testH,
                          script_constant_testMOY,
                          cov_err_test_char,"\n",
                          var_test_t0_char,"\n",
                          var_test_t1_char,"\n",
                          intercept_t0_char,"\n",
                          intercept_t1_char,sep=" ")

result=list()
for (i in 1:B){
  set.seed(i)
  scores_latents <- as.data.frame(mvrnorm(n,mu=c(change_mean,verbal_mean),MatrixChange)) #n latent change and factor scores at time 0, treatment
  scores_latents_controle <- as.data.frame(mvrnorm(n,mu=c(0,verbal_mean),MatrixChange)) #n latent change and factor scores at time 0, control
  resids_low <- as.data.frame(mvrnorm(n*2,mu=rep(0,2*items),Sigma=sigma_res_low)) # residuals with low reliability
  resids_high <- as.data.frame(mvrnorm(n*2,mu=rep(0,2*items),Sigma=sigma_res_high)) #  residuals with high reliability
  resids_wide <- as.data.frame(mvrnorm(n*2,mu=rep(0,2*items),Sigma=sigma_res_wide)) #  residuals with wide reliability
  
  
  #data generating
  dataframe_low <- data.frame(matrix(NA,ncol = 2*items, nrow = N))
  dataframe_high <- data.frame(matrix(NA,ncol = 2*items, nrow = N))
  dataframe_wide <- data.frame(matrix(NA,ncol = 2*items, nrow = N))
  
  colnames(dataframe_low)<-c(manifeste_t0,manifeste_t1)
  colnames(dataframe_high)<-c(manifeste_t0,manifeste_t1)
  colnames(dataframe_wide)<-c(manifeste_t0,manifeste_t1)
  
  dataframe_low$group <- rep(c(0,1),each=n)
  dataframe_high$group <- rep(c(0,1),each=n)
  dataframe_wide$group <- rep(c(0,1),each=n)
  
  for(j in 1:items){
 #manifest at time 0 is obtained multiplying the factor score and loadings and adding the resids. Control are the first n rows
   manifeste_t0_low_control = scores_latents_controle$factor*loads_low[j] + resids_low[1:n,j]
    manifeste_t0_high_control = scores_latents_controle$factor*loads_high[j] + resids_high[1:n,j]
    manifeste_t0_wide_control = scores_latents_controle$factor*loads_wide[j] + resids_wide[1:n,j]

    #manifest at time 1 is obtained by adding the factor and change scores, multiplied by loadings and adding the resids. Control are the first n rows
    manifeste_t1_low_control =  (scores_latents_controle$factor + scores_latents_controle$change)*loads_low[j] + resids_low[1:n,(j+items)]
    manifeste_t1_high_control =  (scores_latents_controle$factor + scores_latents_controle$change)*loads_high[j] + resids_high[1:n,(j+items)]
    manifeste_t1_wide_control =  (scores_latents_controle$factor + scores_latents_controle$change)*loads_wide[j] + resids_wide[1:n,(j+items)]

   #same for treatment
    manifeste_t0_low_treatment = scores_latents$factor*loads_low[j] + resids_low[(n+1):N,j]
    manifeste_t0_high_treatment = scores_latents$factor*loads_high[j] + resids_high[(n+1):N,j]
    manifeste_t0_wide_treatment = scores_latents$factor*loads_wide[j] + resids_wide[(n+1):N,j]
    
    manifeste_t1_low_treatment =  (scores_latents$factor + scores_latents$change)*loads_low[j] + resids_low[(n+1):N,(j+items)]
    manifeste_t1_high_treatment =  (scores_latents$factor + scores_latents$change)*loads_high[j] + resids_high[(n+1):N,(j+items)]
    manifeste_t1_wide_treatment =  (scores_latents$factor + scores_latents$change)*loads_wide[j] + resids_wide[(n+1):N,(j+items)]
    
    dataframe_low[,j] = c(manifeste_t0_low_control,manifeste_t0_low_treatment)
    dataframe_low[,j+items] = c(manifeste_t1_low_control,manifeste_t1_low_treatment)
    dataframe_high[,j] = c(manifeste_t0_high_control,manifeste_t0_high_treatment)
    dataframe_high[,j+items] = c(manifeste_t1_high_control,manifeste_t1_high_treatment)
    dataframe_wide[,j] = c(manifeste_t0_wide_control,manifeste_t0_wide_treatment)
    dataframe_wide[,j+items] = c(manifeste_t1_wide_control,manifeste_t1_wide_treatment)
    
  }
  
  
  
  #Reliability estimates
  relia_low_control1<-omega(dataframe_low[dataframe_low$group==0,c(1:items)],plot = F)
  relia_low_control2<-omega(dataframe_low[dataframe_low$group==0,c((items+1):(2*items))],plot = F)
  relia_low_treat1<-omega(dataframe_low[dataframe_low$group==1,c(1:items)],plot = F)
  relia_low_treat2<-omega(dataframe_low[dataframe_low$group==1,c((items+1):(2*items))],plot = F)
  relia_high_control1<-omega(dataframe_high[dataframe_high$group==0,c(1:items)],plot = F)
  relia_high_control2<-omega(dataframe_high[dataframe_high$group==0,c((items+1):(2*items))],plot = F)
  relia_high_treat1<-omega(dataframe_high[dataframe_high$group==1,c(1:items)],plot = F)
  relia_high_treat2<-omega(dataframe_high[dataframe_high$group==1,c((items+1):(2*items))],plot = F)
  relia_wide_control1<-omega(dataframe_wide[dataframe_wide$group==0,c(1:items)],plot = F)
  relia_wide_control2<-omega(dataframe_wide[dataframe_wide$group==0,c((items+1):(2*items))],plot = F)
  relia_wide_treat1<-omega(dataframe_wide[dataframe_wide$group==1,c(1:items)],plot = F)
  relia_wide_treat2<-omega(dataframe_wide[dataframe_wide$group==1,c((items+1):(2*items))],plot = F)
  
  #Cronbach alpha
  alpha_low_control1<-relia_low_control1$alpha
  alpha_low_control2<-relia_low_control2$alpha
  alpha_low_treat1<-relia_low_treat1$alpha
  alpha_low_treat2<-relia_low_treat2$alpha
  alpha_high_control1<-relia_high_control1$alpha
  alpha_high_control2<-relia_high_control2$alpha
  alpha_high_treat1<-relia_high_treat1$alpha
  alpha_high_treat2<-relia_high_treat2$alpha
  alpha_wide_control1<-relia_wide_control1$alpha
  alpha_wide_control2<-relia_wide_control2$alpha
  alpha_wide_treat1<-relia_wide_treat1$alpha
  alpha_wide_treat2<-relia_wide_treat2$alpha
  
  
  #hierarchical omega
  omega_h_low_control1<-relia_low_control1$omega_h
  omega_h_low_control2<-relia_low_control2$omega_h
  omega_h_low_treat1<-relia_low_treat1$omega_h
  omega_h_low_treat2<-relia_low_treat2$omega_h
  omega_h_high_control1<-relia_high_control1$omega_h
  omega_h_high_control2<-relia_high_control2$omega_h
  omega_h_high_treat1<-relia_high_treat1$omega_h
  omega_h_high_treat2<-relia_high_treat2$omega_h
  omega_h_wide_control1<-relia_wide_control1$omega_h
  omega_h_wide_control2<-relia_wide_control2$omega_h
  omega_h_wide_treat1<-relia_wide_treat1$omega_h
  omega_h_wide_treat2<-relia_wide_treat2$omega_h
  
  #total omega
  omega.tot_low_control1<-relia_low_control1$omega.tot
  omega.tot_low_control2<-relia_low_control2$omega.tot
  omega.tot_low_treat1<-relia_low_treat1$omega.tot
  omega.tot_low_treat2<-relia_low_treat2$omega.tot
  omega.tot_high_control1<-relia_high_control1$omega.tot
  omega.tot_high_control2<-relia_high_control2$omega.tot
  omega.tot_high_treat1<-relia_high_treat1$omega.tot
  omega.tot_high_treat2<-relia_high_treat2$omega.tot
  omega.tot_wide_control1<-relia_wide_control1$omega.tot
  omega.tot_wide_control2<-relia_wide_control2$omega.tot
  omega.tot_wide_treat1<-relia_wide_treat1$omega.tot
  omega.tot_wide_treat2<-relia_wide_treat2$omega.tot
  
  
  #mean scores
  dataframe_low$meanScoreT0 <- rowSums(dataframe_low[,c(1:items)])/items
  dataframe_low$meanScoreT1 <- rowSums(dataframe_low[,c((items+1):(2*items))])/items
  dataframe_low$id <- 1:dim(dataframe_low)[1]
  dataframe_high$meanScoreT0 <- rowSums(dataframe_high[,c(1:items)])/items
  dataframe_high$meanScoreT1 <- rowSums(dataframe_high[,c((items+1):(2*items))])/items
  dataframe_high$id <- 1:dim(dataframe_high)[1]
  dataframe_wide$meanScoreT0 <- rowSums(dataframe_wide[,c(1:items)])/items
  dataframe_wide$meanScoreT1 <- rowSums(dataframe_wide[,c((items+1):(2*items))])/items
  dataframe_wide$id <- 1:dim(dataframe_wide)[1]
  
  #diff scores
  dataframe_low$diffScore <- dataframe_low$meanScoreT1 - dataframe_low$meanScoreT0
  dataframe_high$diffScore <- dataframe_high$meanScoreT1 - dataframe_high$meanScoreT0
  dataframe_wide$diffScore <- dataframe_wide$meanScoreT1 - dataframe_wide$meanScoreT0
  
  #true scores
  dataframe_low$trueFactor <- c(scores_latents_controle$factor,scores_latents$factor)
  dataframe_high$trueFactor <- c(scores_latents_controle$factor,scores_latents$factor)
  dataframe_wide$trueFactor <- c(scores_latents_controle$factor,scores_latents$factor)
  dataframe_low$trueChange <- c(scores_latents_controle$change,scores_latents$change)
  dataframe_high$trueChange <- c(scores_latents_controle$change,scores_latents$change)
  dataframe_wide$trueChange <- c(scores_latents_controle$change,scores_latents$change)
  
  # plot(dataframe_low$diffScore,dataframe_low$trueChange)
  # abline(0,1)
  # plot(dataframe_high$diffScore,dataframe_high$trueChange)
  # abline(0,1)
  # plot(dataframe_wide$diffScore,dataframe_wide$trueChange)
  # abline(0,1)
  
  ###Estimates and tests
  
  #### TESTS
  # WIDE
  if(items==3){estimates <- c(11,13,41,43,61); nb_est = 61  }
  if(items==6){estimates <- c(17,19,68,70,103); nb_est = 103    }
  if(items==12){estimates <- c(29,31,122,124,187); nb_est = 187    }
  # lavaan estimates
  fit_wide <- try(sem(script_testL,
                      data = dataframe_wide,
                      group = "group",
                      group.equal="loadings",
                      meanstructure = T))
  
  chisq_ratiogroup <- 99
  wide_latent_Score0_contr <- NA
  wide_latent_Score0_treat <- NA
  wide_latent_Change_contr <- NA
  wide_latent_Change_treat <- NA
  
  if(class(fit_wide)=="try-error"){
    lavaan_wide<-as.data.frame(matrix(NA,nrow = nb_est,ncol = 6))
    colnames(lavaan_wide)<-c("lhs","op","rhs","group","est","pvalue")
    lavaan_wide$converged <- FALSE
    
    lavaan_wide[estimates,1] <- "Delta"
    lavaan_wide[estimates,2] <- c("~1","~~","~1","~~")
    lavaan_wide[estimates,3] <- c("","Delta","","Delta")
    lavaan_wide[estimates,4] <- c(1,1,2,2)
    
    lavaan_wide$TrueValue <- c(loads_wide,loads_wide,1,1,0,0,mD0,1,varD,rep(0,items),var_resid_wide,var_resid_wide,rep(0,2*items),0,0.25*sqrt(varD),
                               loads_wide,loads_wide,1,1,0,0,mD1,1,varD,rep(0,items),var_resid_wide,var_resid_wide,rep(0,2*items),0,0.25*sqrt(varD),mD1-mD0)
    lavaan_wide$Test <- "factorScore"
    lavaan_wide$chisq <- NA
    lavaan_wide$chisq_df <- NA
    lavaan_wide$chisq_pval <- NA
    lavaan_wide$cfi <- NA
    lavaan_wide$rmsea <- NA
    lavaan_wide$srmr <- NA
    lavaan_wide$chisq_diff <- NA
    
    result_LRT_wide <- data.frame(matrix(NA,3,dim(lavaan_wide)[2]))
    colnames(result_LRT_wide) <- colnames(lavaan_wide)
    result_LRT_wide$TrueValue <-(mD1-mD0)
    result_LRT_wide$Test <-NA
    result_LRT_wide$lhs <- "group by time"
  }else{
    fit_wide <- sem(script_testL,
                    data = dataframe_wide,
                    group = "group",
                    group.equal="loadings",
                    meanstructure = T)
    
    #fitwide_group0 <- sem(script_testL, data = subset(dataframe_wide, group == 0))
    #aa=summary(fit_wide, rsquare=TRUE, fit.measures=T, standardized=TRUE)
    lavaan_wide<-as.data.frame(parameterEstimates(fit_wide)[,c(1:3,5,7,10)])
    lavaan_wide$converged <- lavInspect(fit_wide,"converged")
    
    lavaan_wide$TrueValue <- c(loads_wide,loads_wide,1,1,0,0,mD0,1,varD,rep(0,items),var_resid_wide,var_resid_wide,rep(0,2*items),0,0.25*sqrt(varD),
                               loads_wide,loads_wide,1,1,0,0,mD1,1,varD,rep(0,items),var_resid_wide,var_resid_wide,rep(0,2*items),0,0.25*sqrt(varD),mD1-mD0)
    
    lavaan_wide$chisq <- NA
    lavaan_wide$chisq_df <- NA
    lavaan_wide$chisq_pval <- NA
    lavaan_wide$cfi <- NA
    lavaan_wide$rmsea <- NA
    lavaan_wide$srmr <- NA
    lavaan_wide$chisq_diff <- NA
    
    result_LRT_wide <- data.frame(matrix(NA,3,dim(lavaan_wide)[2]))
    colnames(result_LRT_wide) <- colnames(lavaan_wide)
    result_LRT_wide$TrueValue <-(mD1-mD0)
    result_LRT_wide$Test <-NA
    result_LRT_wide$lhs <- "group by time"
    result_LRT_wide$chisq <- NA
    result_LRT_wide$chisq_df <- NA
    result_LRT_wide$chisq_pval <- NA
    result_LRT_wide$cfi <- NA
    result_LRT_wide$rmsea <- NA
    result_LRT_wide$srmr <- NA
    wide_latent_Score0_contr <- mean(predict(fit_wide)$`0`[,1])
    wide_latent_Score0_treat <- mean(predict(fit_wide)$`1`[,1])
    wide_latent_Change_contr <- mean(predict(fit_wide)$`0`[,3])
    wide_latent_Change_treat <- mean(predict(fit_wide)$`1`[,3])
    
    if(lavInspect(fit_wide,"converged")==TRUE ){#& lavInspect(fitwide_group0,"converged")==TRUE
      lavaan_wide$converged = TRUE
      lavaan_wide$chisq <- tryCatch(fitmeasures(fit_wide)[3])
      lavaan_wide$chisq_df <- tryCatch(fitmeasures(fit_wide)[4])
      lavaan_wide$chisq_pval <- tryCatch(fitmeasures(fit_wide)[5])
      lavaan_wide$cfi <- tryCatch(fitmeasures(fit_wide)[9])
      lavaan_wide$rmsea <- tryCatch(fitmeasures(fit_wide)[23])
      lavaan_wide$srmr <- tryCatch(fitmeasures(fit_wide)[33])
      #chisq_ratiogroup <- tryCatch(fitmeasures(fitwide_group0)[3])/tryCatch(fitmeasures(fit_wide)[3])
      
      fit_wide_m <- try(sem(script_testLMOY,
                            data = dataframe_wide,
                            group = "group",
                            group.equal="loadings",
                            meanstructure = T))
      
      if(class(fit_wide_m)!="try-error"){
        fit_wide_m <- sem(script_testLMOY,
                          data = dataframe_wide,
                          group = "group",
                          group.equal="loadings",
                          meanstructure = T)
        LRT_wide_result <- anova(fit_wide,fit_wide_m)
        result_LRT_wide[1,]$est <- lavaan_wide[estimates[3],"est"]-lavaan_wide[estimates[1],"est"]
        result_LRT_wide[1,]$pvalue <- LRT_wide_result$`Pr(>Chisq)`[2]
        result_LRT_wide$chisq_diff <- LRT_wide_result$`Chisq diff`[2]
      }
    }
  }
  
  
  lavaan_wide$Test <- "factorScore"
  #lmer estimates
  long_wide<-gather(dataframe_wide[,(2*items+1):(2*items+4)],Time, Score,meanScoreT0:meanScoreT1, factor_key = T )
  long_wide$group <- factor(long_wide$group)
  lmer_wide<-lmer(Score~Time*group +(1|id),long_wide,REML=F)
  lmer_wide_result <- summary(lmer_wide)
  
  long_wide$groupB <- factor(long_wide$group, levels = c(1,0))
  lmer_wide2<-lmer(Score~Time*groupB +(1|id),long_wide,REML=F)
  lmer_wide2_result <- summary(lmer_wide2)
  
  result_lmer_wide <- data.frame(matrix(NA,4,dim(lavaan_wide)[2]))
  colnames(result_lmer_wide) <- colnames(lavaan_wide)
  result_lmer_wide[,1] <- "Delta"
  result_lmer_wide[,2] <- c("~1","~~","~1","~~")
  result_lmer_wide[,3] <- c("","Delta","","Delta")
  result_lmer_wide[,4] <- c(1,1,2,2)
  
  result_lmer_wide$est <- c(fixef(lmer_wide)[2],NA,
                            fixef(lmer_wide2)[2],NA)
  result_lmer_wide$pvalue <- c(lmer_wide_result$coefficients[2,5],NA,lmer_wide2_result$coefficients[2,5],NA)
  result_lmer_wide$TrueValue <- c(mD0,varD,mD1,varD)
  result_lmer_wide$chisq <- NA
  result_lmer_wide$chisq_df <- NA
  result_lmer_wide$chisq_pval <- NA
  result_lmer_wide$cfi <- NA
  result_lmer_wide$rmsea <- NA
  result_lmer_wide$srmr <- NA
  result_lmer_wide$Test <- "meanScore"
  result_lmer_wide$chisq_diff <- NA
  #t test estimates
  ttestdata_wide<-data.frame(diffCont=dataframe_wide[dataframe_wide$group==0,]$diffScore,diffTreat=dataframe_wide[dataframe_wide$group==1,]$diffScore)
  ttest_wide_result<-t.test(ttestdata_wide$diffTreat,ttestdata_wide$diffCont,var.equal=F)
  ttest_wide_Cont<-t.test(ttestdata_wide$diffCont,var.equal=F)
  ttest_wide_Treat<-t.test(ttestdata_wide$diffTreat,var.equal=F)
  
  result_ttest_wide <- data.frame(matrix(NA,4,dim(lavaan_wide)[2]))
  colnames(result_ttest_wide) <- colnames(lavaan_wide)
  result_ttest_wide[,1] <- "Delta"
  result_ttest_wide[,2] <- c("~1","~~","~1","~~")
  result_ttest_wide[,3] <- c("","Delta","","Delta")
  result_ttest_wide[,4] <- c(1,1,2,2)
  result_ttest_wide$est <- c(ttest_wide_result$estimate[2],sd(ttestdata_wide$diffCont)^2,
                             ttest_wide_result$estimate[1],sd(ttestdata_wide$diffTreat)^2)
  result_ttest_wide$pvalue <- c(ttest_wide_Cont$p.value,NA,ttest_wide_Treat$p.value,NA)
  result_ttest_wide$TrueValue <- c(mD0,varD,mD1,varD)
  result_ttest_wide$chisq <- NA
  result_ttest_wide$chisq_df <- NA
  result_ttest_wide$chisq_pval <- NA
  result_ttest_wide$cfi <- NA
  result_ttest_wide$rmsea <- NA
  result_ttest_wide$srmr <- NA
  result_ttest_wide$Test <- "diffScoreST"
  result_ttest_wide$chisq_diff <- NA
  
  result_LRT_wide[1,]$Test <- "factorScore"
  result_LRT_wide[2,]$est <- ttest_wide_result$estimate[1]-ttest_wide_result$estimate[2]
  result_LRT_wide[2,]$pvalue <- ttest_wide_result$p.value
  result_LRT_wide[2,]$Test <- "diffScoreST"
  result_LRT_wide[3,]$est <- fixef(lmer_wide)[4]
  result_LRT_wide[3,]$pvalue <- lmer_wide_result$coefficients[4,5]
  result_ttest_wide$chisq <- NA
  result_LRT_wide[3,]$Test <- "meanScore"
  
  results_wide<-rbind(lavaan_wide,result_lmer_wide,result_ttest_wide,result_LRT_wide)
  results_wide$alpha_control1 <- alpha_wide_control1
  results_wide$alpha_treat1 <- alpha_wide_treat1
  results_wide$alpha_control2 <- alpha_wide_control2
  results_wide$alpha_treat2 <- alpha_wide_treat2
  results_wide$omega_h_control1 <- omega_h_wide_control1
  results_wide$omega_h_treat1 <- omega_h_wide_treat1
  results_wide$omega_h_control2 <- omega_h_wide_control2
  results_wide$omega_h_treat2 <- omega_h_wide_treat2
  results_wide$omega.tot_control1 <- omega.tot_wide_control1
  results_wide$omega.tot_treat1 <- omega.tot_wide_treat1
  results_wide$omega.tot_control2 <- omega.tot_wide_control2
  results_wide$omega.tot_treat2 <- omega.tot_wide_treat2
  results_wide$mean_TrueChange_contr <- mean(dataframe_wide[dataframe_wide$group==0,]$trueChange)
  results_wide$mean_TrueFactor_contr <- mean(dataframe_wide[dataframe_wide$group==0,]$trueFactor)
  results_wide$mean_diffScore_contr <- mean(dataframe_wide[dataframe_wide$group==0,]$diffScore)
  results_wide$mean_ScoreT0_contr <- mean(dataframe_wide[dataframe_wide$group==0,]$meanScoreT0)
  results_wide$mean_TrueChange_treat <- mean(dataframe_wide[dataframe_wide$group==1,]$trueChange)
  results_wide$mean_TrueFactor_treat <- mean(dataframe_wide[dataframe_wide$group==1,]$trueFactor)
  results_wide$mean_diffScore_treat <- mean(dataframe_wide[dataframe_wide$group==1,]$diffScore)
  results_wide$mean_ScoreT0_treat <- mean(dataframe_wide[dataframe_wide$group==1,]$meanScoreT0)
  results_wide$bias_TrueChange_contr <- results_wide$mean_TrueChange_contr - results_wide$mean_diffScore_contr
  results_wide$bias_TrueFactor_contr <- results_wide$mean_TrueFactor_contr - results_wide$mean_ScoreT0_contr
  results_wide$bias_TrueChange_treat <- results_wide$mean_TrueChange_treat - results_wide$mean_diffScore_treat
  results_wide$bias_TrueFactor_treat <- results_wide$mean_TrueFactor_treat - results_wide$mean_ScoreT0_treat
  #results_wide$chisq_ratiogroup <- chisq_ratiogroup
  results_wide$latent_Score0_contr <- wide_latent_Score0_contr
  results_wide$latent_Score0_treat <- wide_latent_Score0_treat
  results_wide$latent_Change_contr <- wide_latent_Change_contr
  results_wide$latent_Change_treat <- wide_latent_Change_treat
  
  results_wide$Condition <- "wide"
  
  # HIGH
  fit_high <- try(sem(script_testH,
                      data = dataframe_high,
                      group = "group",
                      group.equal="loadings",
                      meanstructure = T))
  
  
  chisq_ratiogroup <- 99
  high_latent_Score0_contr <- NA
  high_latent_Score0_treat <- NA
  high_latent_Change_contr <- NA
  high_latent_Change_treat <- NA
  
  if(class(fit_high)=="try-error"){
    lavaan_high<-as.data.frame(matrix(NA,nrow = nb_est,ncol = 6))
    colnames(lavaan_high)<-c("lhs","op","rhs","group","est","pvalue")
    lavaan_high$converged <- FALSE
    
    lavaan_high[estimates,1] <- "Delta"
    lavaan_high[estimates,2] <- c("~1","~~","~1","~~")
    lavaan_high[estimates,3] <- c("","Delta","","Delta")
    lavaan_high[estimates,4] <- c(1,1,2,2)
    
    lavaan_high$TrueValue <- c(loads_high,loads_high,1,1,0,0,mD0,1,varD,rep(0,items),var_resid_high,var_resid_high,rep(0,2*items),0,0.25*sqrt(varD),
                               loads_high,loads_high,1,1,0,0,mD1,1,varD,rep(0,items),var_resid_high,var_resid_high,rep(0,2*items),0,0.25*sqrt(varD),mD1-mD0)
    lavaan_high$Test <- "factorScore"
    lavaan_high$chisq <- NA
    lavaan_high$chisq_df <- NA
    lavaan_high$chisq_pval <- NA
    lavaan_high$cfi <- NA
    lavaan_high$rmsea <- NA
    lavaan_high$srmr <- NA
    lavaan_high$chisq_diff <- NA
    
    result_LRT_high <- data.frame(matrix(NA,3,dim(lavaan_high)[2]))
    colnames(result_LRT_high) <- colnames(lavaan_high)
    result_LRT_high$TrueValue <-(mD1-mD0)
    result_LRT_high$Test <-NA
    result_LRT_high$lhs <- "group by time"
    result_LRT_high$chisq_diff <- NA
  }else{
    fit_high <- sem(script_testH,
                    data = dataframe_high,
                    group = "group",
                    group.equal="loadings",
                    meanstructure = T)
    #fithigh_group0 <- sem(script_testH, data = subset(dataframe_high, group == 0))
    
    #aa=summary(fit_high, rsquare=TRUE, fit.measures=T, standardized=TRUE)
    lavaan_high<-as.data.frame(parameterEstimates(fit_high)[,c(1:3,5,7,10)])
    lavaan_high$converged <- lavInspect(fit_high,"converged")
    
    lavaan_high$TrueValue <- c(loads_high,loads_high,1,1,0,0,mD0,1,varD,rep(0,items),var_resid_high,var_resid_high,rep(0,2*items),0,0.25*sqrt(varD),
                               loads_high,loads_high,1,1,0,0,mD1,1,varD,rep(0,items),var_resid_high,var_resid_high,rep(0,2*items),0,0.25*sqrt(varD),mD1-mD0)
    
    lavaan_high$chisq <- NA
    lavaan_high$chisq_df <- NA
    lavaan_high$chisq_pval <- NA
    lavaan_high$cfi <- NA
    lavaan_high$rmsea <- NA
    lavaan_high$srmr <- NA
    lavaan_high$chisq_diff <- NA
    
    result_LRT_high <- data.frame(matrix(NA,3,dim(lavaan_high)[2]))
    colnames(result_LRT_high) <- colnames(lavaan_high)
    result_LRT_high$TrueValue <-(mD1-mD0)
    result_LRT_high$Test <-NA
    result_LRT_high$lhs <- "group by time"
    result_LRT_high$chisq <- NA
    result_LRT_high$chisq_df <- NA
    result_LRT_high$chisq_pval <- NA
    result_LRT_high$cfi <- NA
    result_LRT_high$rmsea <- NA
    result_LRT_high$srmr <- NA
    result_LRT_high$chisq_diff <- NA
    high_latent_Score0_contr <- mean(predict(fit_high)$`0`[,1])
    high_latent_Score0_treat <- mean(predict(fit_high)$`1`[,1])
    high_latent_Change_contr <- mean(predict(fit_high)$`0`[,3])
    high_latent_Change_treat <- mean(predict(fit_high)$`1`[,3])
    
    if(lavInspect(fit_high,"converged")==TRUE){# & lavInspect(fithigh_group0,"converged")==TRUE
      lavaan_high$converged = TRUE
      lavaan_high$chisq <- tryCatch(fitmeasures(fit_high)[3])
      lavaan_high$chisq_df <- tryCatch(fitmeasures(fit_high)[4])
      lavaan_high$chisq_pval <- tryCatch(fitmeasures(fit_high)[5])
      lavaan_high$cfi <- tryCatch(fitmeasures(fit_high)[9])
      lavaan_high$rmsea <- tryCatch(fitmeasures(fit_high)[23])
      lavaan_high$srmr <- tryCatch(fitmeasures(fit_high)[33])
      #chisq_ratiogroup <- tryCatch(fitmeasures(fithigh_group0)[3])/tryCatch(fitmeasures(fit_high)[3])
      
      
      fit_high_m <- try(sem(script_testHMOY,
                            data = dataframe_high,
                            group = "group",
                            group.equal="loadings",
                            meanstructure = T))
      
      if(class(fit_high_m)!="try-error"){
        fit_high_m <- sem(script_testHMOY,
                          data = dataframe_high,
                          group = "group",
                          group.equal="loadings",
                          meanstructure = T)
        LRT_high_result <- anova(fit_high,fit_high_m)
        result_LRT_high[1,]$est <- lavaan_high[estimates[3],"est"]-lavaan_high[estimates[1],"est"]
        result_LRT_high[1,]$pvalue <- LRT_high_result$`Pr(>Chisq)`[2]
        result_LRT_high$chisq_diff <- LRT_high_result$`Chisq diff`[2]
      }
    }
  }
  
  lavaan_high$Test <- "factorScore"
  
  long_high<-gather(dataframe_high[,(2*items+1):(2*items+4)],Time, Score,meanScoreT0:meanScoreT1, factor_key = T )
  long_high$group <- factor(long_high$group)
  lmer_high<-lmer(Score~Time*group +(1|id),long_high,REML=F)
  lmer_high_result <- summary(lmer_high)
  
  long_high$groupB <- factor(long_high$group, levels = c(1,0))
  lmer_high2<-lmer(Score~Time*groupB +(1|id),long_high,REML=F)
  lmer_high2_result <- summary(lmer_high2)
  
  result_lmer_high <- data.frame(matrix(NA,4,dim(lavaan_high)[2]))
  colnames(result_lmer_high) <- colnames(lavaan_high)
  result_lmer_high[,1] <- "Delta"
  result_lmer_high[,2] <- c("~1","~~","~1","~~")
  result_lmer_high[,3] <- c("","Delta","","Delta")
  result_lmer_high[,4] <- c(1,1,2,2)
  
  result_lmer_high$est <- c(fixef(lmer_high)[2],NA,
                            fixef(lmer_high2)[2],NA)
  result_lmer_high$pvalue <- c(lmer_high_result$coefficients[2,5],NA,lmer_high2_result$coefficients[2,5],NA)
  result_lmer_high$TrueValue <- c(mD0,varD,mD1,varD)
  result_lmer_high$chisq <- NA
  result_lmer_high$chisq_df <- NA
  result_lmer_high$chisq_pval <- NA
  result_lmer_high$cfi <- NA
  result_lmer_high$rmsea <- NA
  result_lmer_high$srmr <- NA
  result_lmer_high$Test <- "meanScore"
  result_lmer_high$chisq_diff <- NA
  
  
  ttestdata_high<-data.frame(diffCont=dataframe_high[dataframe_high$group==0,]$diffScore,diffTreat=dataframe_high[dataframe_high$group==1,]$diffScore)
  ttest_high_result<-t.test(ttestdata_high$diffTreat,ttestdata_high$diffCont,var.equal=F)
  ttest_high_Cont<-t.test(ttestdata_high$diffCont,var.equal=F)
  ttest_high_Treat<-t.test(ttestdata_high$diffTreat,var.equal=F)
  
  result_ttest_high <- data.frame(matrix(NA,4,dim(lavaan_high)[2]))
  colnames(result_ttest_high) <- colnames(lavaan_high)
  result_ttest_high[,1] <- "Delta"
  result_ttest_high[,2] <- c("~1","~~","~1","~~")
  result_ttest_high[,3] <- c("","Delta","","Delta")
  result_ttest_high[,4] <- c(1,1,2,2)
  result_ttest_high$est <- c(ttest_high_result$estimate[2],sd(ttestdata_high$diffCont)^2,
                             ttest_high_result$estimate[1],sd(ttestdata_high$diffTreat)^2)
  result_ttest_high$pvalue <- c(ttest_high_Cont$p.value,NA,ttest_high_Treat$p.value,NA)
  result_ttest_high$TrueValue <- c(mD0,varD,mD1,varD)
  result_ttest_high$chisq <- NA
  result_ttest_high$chisq_df <- NA
  result_ttest_high$chisq_pval <- NA
  result_ttest_high$cfi <- NA
  result_ttest_high$rmsea <- NA
  result_ttest_high$srmr <- NA
  result_ttest_high$Test <- "diffScoreST"
  result_ttest_high$chisq_diff <- NA
  
  
  
  result_LRT_high[1,]$Test <- "factorScore"
  result_LRT_high[2,]$est <- ttest_high_result$estimate[1]-ttest_high_result$estimate[2]
  result_LRT_high[2,]$pvalue <- ttest_high_result$p.value
  result_LRT_high[2,]$Test <- "diffScoreST"
  result_LRT_high[3,]$est <- fixef(lmer_high)[4]
  result_LRT_high[3,]$pvalue <- lmer_high_result$coefficients[4,5]
  result_ttest_high$chisq <- NA
  result_LRT_high[3,]$Test <- "meanScore"
  
  results_high<-rbind(lavaan_high,result_lmer_high,result_ttest_high,result_LRT_high)
  results_high$alpha_control1 <- alpha_high_control1
  results_high$alpha_treat1 <- alpha_high_treat1
  results_high$alpha_control2 <- alpha_high_control2
  results_high$alpha_treat2 <- alpha_high_treat2
  results_high$omega_h_control1 <- omega_h_high_control1
  results_high$omega_h_treat1 <- omega_h_high_treat1
  results_high$omega_h_control2 <- omega_h_high_control2
  results_high$omega_h_treat2 <- omega_h_high_treat2
  results_high$omega.tot_control1 <- omega.tot_high_control1
  results_high$omega.tot_treat1 <- omega.tot_high_treat1
  results_high$omega.tot_control2 <- omega.tot_high_control2
  results_high$omega.tot_treat2 <- omega.tot_high_treat2
  results_high$mean_TrueChange_contr <- mean(dataframe_high[dataframe_high$group==0,]$trueChange)
  results_high$mean_TrueFactor_contr <- mean(dataframe_high[dataframe_high$group==0,]$trueFactor)
  results_high$mean_diffScore_contr <- mean(dataframe_high[dataframe_high$group==0,]$diffScore)
  results_high$mean_ScoreT0_contr <- mean(dataframe_high[dataframe_high$group==0,]$meanScoreT0)
  results_high$mean_TrueChange_treat <- mean(dataframe_high[dataframe_high$group==1,]$trueChange)
  results_high$mean_TrueFactor_treat <- mean(dataframe_high[dataframe_high$group==1,]$trueFactor)
  results_high$mean_diffScore_treat <- mean(dataframe_high[dataframe_high$group==1,]$diffScore)
  results_high$mean_ScoreT0_treat <- mean(dataframe_high[dataframe_high$group==1,]$meanScoreT0)
  results_high$bias_TrueChange_contr <- results_high$mean_TrueChange_contr - results_high$mean_diffScore_contr
  results_high$bias_TrueFactor_contr <- results_high$mean_TrueFactor_contr - results_high$mean_ScoreT0_contr
  results_high$bias_TrueChange_treat <- results_high$mean_TrueChange_treat - results_high$mean_diffScore_treat
  results_high$bias_TrueFactor_treat <- results_high$mean_TrueFactor_treat - results_high$mean_ScoreT0_treat
 # results_high$chisq_ratiogroup <- chisq_ratiogroup
  results_high$latent_Score0_contr <- high_latent_Score0_contr
  results_high$latent_Score0_treat <- high_latent_Score0_treat
  results_high$latent_Change_contr <- high_latent_Change_contr
  results_high$latent_Change_treat <- high_latent_Change_treat
  
  results_high$Condition <- "high"
  
  
  # LOW
  fit_low <- try(sem(script_testL,
                     data = dataframe_low,
                     group = "group",
                     group.equal="loadings",
                     meanstructure = T))
  
  chisq_ratiogroup <- 99
  low_latent_Score0_contr <- NA
  low_latent_Score0_treat <- NA
  low_latent_Change_contr <- NA
  low_latent_Change_treat <- NA
  
  if(class(fit_low)=="try-error"){
    lavaan_low<-as.data.frame(matrix(NA,nrow = nb_est,ncol = 6))
    colnames(lavaan_low)<-c("lhs","op","rhs","group","est","pvalue")
    lavaan_low$converged <- FALSE
    
    lavaan_low[estimates,1] <- "Delta"
    lavaan_low[estimates,2] <- c("~1","~~","~1","~~")
    lavaan_low[estimates,3] <- c("","Delta","","Delta")
    lavaan_low[estimates,4] <- c(1,1,2,2)
    
    lavaan_low$TrueValue <- c(loads_low,loads_low,1,1,0,0,mD0,1,varD,rep(0,items),var_resid_low,var_resid_low,rep(0,2*items),0,0.25*sqrt(varD),
                              loads_low,loads_low,1,1,0,0,mD1,1,varD,rep(0,items),var_resid_low,var_resid_low,rep(0,2*items),0,0.25*sqrt(varD),mD1-mD0)
    lavaan_low$Test <- "factorScore"
    lavaan_low$chisq <- NA
    lavaan_low$chisq_df <- NA
    lavaan_low$chisq_pval <- NA
    lavaan_low$cfi <- NA
    lavaan_low$rmsea <- NA
    lavaan_low$srmr <- NA
    lavaan_low$chisq_diff <- NA
    
    result_LRT_low <- data.frame(matrix(NA,3,dim(lavaan_low)[2]))
    colnames(result_LRT_low) <- colnames(lavaan_low)
    result_LRT_low$TrueValue <-(mD1-mD0)
    result_LRT_low$Test <-NA
    result_LRT_low$lhs <- "group by time"
    result_LRT_low$chisq_diff <- NA
    
  }else{
    fit_low <- sem(script_testL,
                   data = dataframe_low,
                   group = "group",
                   group.equal="loadings",
                   meanstructure = T)
    #fitlow_group0 <- sem(script_testL, data = subset(dataframe_low, group == 0))
    
    #aa=summary(fit_low, rsquare=TRUE, fit.measures=T, standardized=TRUE)
    lavaan_low<-as.data.frame(parameterEstimates(fit_low)[,c(1:3,5,7,10)])
    lavaan_low$converged <- lavInspect(fit_low,"converged")
    
    lavaan_low$TrueValue <- c(loads_low,loads_low,1,1,0,0,mD0,1,varD,rep(0,items),var_resid_low,var_resid_low,rep(0,2*items),0,0.25*sqrt(varD),
                              loads_low,loads_low,1,1,0,0,mD1,1,varD,rep(0,items),var_resid_low,var_resid_low,rep(0,2*items),0,0.25*sqrt(varD),mD1-mD0)
    lavaan_low$chisq <- NA
    lavaan_low$chisq_df <- NA
    lavaan_low$chisq_pval <- NA
    lavaan_low$cfi <- NA
    lavaan_low$rmsea <- NA
    lavaan_low$srmr <- NA
    lavaan_low$chisq_diff <- NA
    
    result_LRT_low <- data.frame(matrix(NA,3,dim(lavaan_low)[2]))
    colnames(result_LRT_low) <- colnames(lavaan_low)
    result_LRT_low$TrueValue <-(mD1-mD0)
    result_LRT_low$Test <-NA
    result_LRT_low$lhs <- "group by time"
    result_LRT_low$chisq <- NA
    result_LRT_low$chisq_df <- NA
    result_LRT_low$chisq_pval <- NA
    result_LRT_low$cfi <- NA
    result_LRT_low$rmsea <- NA
    result_LRT_low$srmr <- NA
    result_LRT_low$chisq_diff <- NA  
    low_latent_Score0_contr <- mean(predict(fit_low)$`0`[,1])
    low_latent_Score0_treat <- mean(predict(fit_low)$`1`[,1])
    low_latent_Change_contr <- mean(predict(fit_low)$`0`[,3])
    low_latent_Change_treat <- mean(predict(fit_low)$`1`[,3])
    
    
    
    if(lavInspect(fit_low,"converged")==TRUE){# & lavInspect(fitlow_group0,"converged")==TRUE
      lavaan_low$converged = TRUE
      lavaan_low$chisq <- tryCatch(fitmeasures(fit_low)[3])
      lavaan_low$chisq_df <- tryCatch(fitmeasures(fit_low)[4])
      lavaan_low$chisq_pval <- tryCatch(fitmeasures(fit_low)[5])
      lavaan_low$cfi <- tryCatch(fitmeasures(fit_low)[9])
      lavaan_low$rmsea <- tryCatch(fitmeasures(fit_low)[23])
      lavaan_low$srmr <- tryCatch(fitmeasures(fit_low)[33])
      #chisq_ratiogroup <- tryCatch(fitmeasures(fitlow_group0)[3])/tryCatch(fitmeasures(fit_low)[3])
      
      
      fit_low_m <- try(sem(script_testLMOY,
                           data = dataframe_low,
                           group = "group",
                           group.equal="loadings",
                           meanstructure = T))
      
      if(class(fit_low_m)!="try-error"){
        fit_low_m <- sem(script_testLMOY,
                         data = dataframe_low,
                         group = "group",
                         group.equal="loadings",
                         meanstructure = T)
        LRT_low_result <- anova(fit_low,fit_low_m)
        result_LRT_low[1,]$est <- lavaan_low[estimates[3],"est"]-lavaan_low[estimates[1],"est"]
        result_LRT_low[1,]$pvalue <- LRT_low_result$`Pr(>Chisq)`[2]
        result_LRT_low$chisq_diff <- LRT_low_result$`Chisq diff`[2]
      }
    }
  }
  
  lavaan_low$Test <- "factorScore"
  
  long_low<-gather(dataframe_low[,(2*items+1):(2*items+4)],Time, Score,meanScoreT0:meanScoreT1, factor_key = T )
  long_low$group <- factor(long_low$group)
  lmer_low<-lmer(Score~Time*group +(1|id),long_low,REML=F)
  lmer_low_result <- summary(lmer_low)
  
  long_low$groupB <- factor(long_low$group, levels = c(1,0))
  lmer_low2<-lmer(Score~Time*groupB +(1|id),long_low,REML=F)
  lmer_low2_result <- summary(lmer_low2)
  
  result_lmer_low <- data.frame(matrix(NA,4,dim(lavaan_low)[2]))
  colnames(result_lmer_low) <- colnames(lavaan_low)
  result_lmer_low[,1] <- "Delta"
  result_lmer_low[,2] <- c("~1","~~","~1","~~")
  result_lmer_low[,3] <- c("","Delta","","Delta")
  result_lmer_low[,4] <- c(1,1,2,2)
  
  result_lmer_low$est <- c(fixef(lmer_low)[2],NA,
                           fixef(lmer_low2)[2],NA)
  result_lmer_low$pvalue <- c(lmer_low_result$coefficients[2,5],NA,lmer_low2_result$coefficients[2,5],NA)
  result_lmer_low$TrueValue <- c(mD0,varD,mD1,varD)
  result_lmer_low$chisq <- NA
  result_lmer_low$chisq_df <- NA
  result_lmer_low$chisq_pval <- NA
  result_lmer_low$cfi <- NA
  result_lmer_low$rmsea <- NA
  result_lmer_low$srmr <- NA
  result_lmer_low$Test <- "meanScore"
  result_lmer_low$chisq_diff <- NA
  
  
  ttestdata_low<-data.frame(diffCont=dataframe_low[dataframe_low$group==0,]$diffScore,diffTreat=dataframe_low[dataframe_low$group==1,]$diffScore)
  ttest_low_result<-t.test(ttestdata_low$diffTreat,ttestdata_low$diffCont,var.equal=F)
  ttest_low_Cont<-t.test(ttestdata_low$diffCont,var.equal=F)
  ttest_low_Treat<-t.test(ttestdata_low$diffTreat,var.equal=F)
  
  result_ttest_low <- data.frame(matrix(NA,4,dim(lavaan_low)[2]))
  colnames(result_ttest_low) <- colnames(lavaan_low)
  result_ttest_low[,1] <- "Delta"
  result_ttest_low[,2] <- c("~1","~~","~1","~~")
  result_ttest_low[,3] <- c("","Delta","","Delta")
  result_ttest_low[,4] <- c(1,1,2,2)
  result_ttest_low$est <- c(ttest_low_result$estimate[2],sd(ttestdata_low$diffCont)^2,
                            ttest_low_result$estimate[1],sd(ttestdata_low$diffTreat)^2)
  result_ttest_low$pvalue <- c(ttest_low_Cont$p.value,NA,ttest_low_Treat$p.value,NA)
  result_ttest_low$TrueValue <- c(mD0,varD,mD1,varD)
  result_ttest_low$chisq <- NA
  result_ttest_low$chisq_df <- NA
  result_ttest_low$chisq_pval <- NA
  result_ttest_low$cfi <- NA
  result_ttest_low$rmsea <- NA
  result_ttest_low$srmr <- NA
  result_ttest_low$Test <- "diffScoreST"
  result_ttest_low$chisq_diff <- NA
  
  
  result_LRT_low[1,]$Test <- "factorScore"
  result_LRT_low[2,]$est <- ttest_low_result$estimate[1]-ttest_low_result$estimate[2]
  result_LRT_low[2,]$pvalue <- ttest_low_result$p.value
  result_LRT_low[2,]$Test <- "diffScoreST"
  result_LRT_low[3,]$est <- fixef(lmer_low)[4]
  result_LRT_low[3,]$pvalue <- lmer_low_result$coefficients[4,5]
  result_ttest_low$chisq <- NA
  result_LRT_low[3,]$Test <- "meanScore"
  
  results_low<-rbind(lavaan_low,result_lmer_low,result_ttest_low,result_LRT_low)
  results_low$alpha_control1 <- alpha_low_control1
  results_low$alpha_treat1 <- alpha_low_treat1
  results_low$alpha_control2 <- alpha_low_control2
  results_low$alpha_treat2 <- alpha_low_treat2
  results_low$omega_h_control1 <- omega_h_low_control1
  results_low$omega_h_treat1 <- omega_h_low_treat1
  results_low$omega_h_control2 <- omega_h_low_control2
  results_low$omega_h_treat2 <- omega_h_low_treat2
  results_low$omega.tot_control1 <- omega.tot_low_control1
  results_low$omega.tot_treat1 <- omega.tot_low_treat1
  results_low$omega.tot_control2 <- omega.tot_low_control2
  results_low$omega.tot_treat2 <- omega.tot_low_treat2
  results_low$mean_TrueChange_contr <- mean(dataframe_low[dataframe_low$group==0,]$trueChange)
  results_low$mean_TrueFactor_contr <- mean(dataframe_low[dataframe_low$group==0,]$trueFactor)
  results_low$mean_diffScore_contr <- mean(dataframe_low[dataframe_low$group==0,]$diffScore)
  results_low$mean_ScoreT0_contr <- mean(dataframe_low[dataframe_low$group==0,]$meanScoreT0)
  results_low$mean_TrueChange_treat <- mean(dataframe_low[dataframe_low$group==1,]$trueChange)
  results_low$mean_TrueFactor_treat <- mean(dataframe_low[dataframe_low$group==1,]$trueFactor)
  results_low$mean_diffScore_treat <- mean(dataframe_low[dataframe_low$group==1,]$diffScore)
  results_low$mean_ScoreT0_treat <- mean(dataframe_low[dataframe_low$group==1,]$meanScoreT0)
  results_low$bias_TrueChange_contr <- results_low$mean_TrueChange_contr - results_low$mean_diffScore_contr
  results_low$bias_TrueFactor_contr <- results_low$mean_TrueFactor_contr - results_low$mean_ScoreT0_contr
  results_low$bias_TrueChange_treat <- results_low$mean_TrueChange_treat - results_low$mean_diffScore_treat
  results_low$bias_TrueFactor_treat <- results_low$mean_TrueFactor_treat - results_low$mean_ScoreT0_treat
 # results_low$chisq_ratiogroup <- chisq_ratiogroup
  
  results_low$latent_Score0_contr <- low_latent_Score0_contr
  results_low$latent_Score0_treat <- low_latent_Score0_treat
  results_low$latent_Change_contr <- low_latent_Change_contr
  results_low$latent_Change_treat <- low_latent_Change_treat
  
  results_low$Condition <- "low"
  
  #Satterthwaite with lavaan
  fit_wideSATT <- sem(script_testL,
                      data = dataframe_wide,
                      estimator="MLMVS",
                      group = "group",
                      group.equal="loadings",
                      meanstructure = T)
  lavaan_wideSATT<-as.data.frame(parameterEstimates(fit_wideSATT)[nb_est,c(1:3,5,7,10)])
  lavaan_wideSATT$TrueValue <- mD1-mD0
  lavaan_wideSATT$Test <- "Satterthwaite"
  lavaan_wideSATT$Condition <- "wide"
  fit_highSATT <- sem(script_testH,
                      data = dataframe_high,
                      estimator="MLMVS",
                      group = "group",
                      group.equal="loadings",
                      meanstructure = T)
  lavaan_highSATT<-as.data.frame(parameterEstimates(fit_highSATT)[nb_est,c(1:3,5,7,10)])
  lavaan_highSATT$TrueValue <- mD1-mD0
  lavaan_highSATT$Test <- "Satterthwaite"
  lavaan_highSATT$Condition <- "high"
  
  fit_lowSATT <- sem(script_testL,
                     data = dataframe_low,
                     estimator="MLMVS",
                     group = "group",
                     group.equal="loadings",
                     meanstructure = T)
  lavaan_lowSATT<-as.data.frame(parameterEstimates(fit_lowSATT)[nb_est,c(1:3,5,7,10)])
  lavaan_lowSATT$TrueValue <- mD1-mD0
  lavaan_lowSATT$Test <- "Satterthwaite"
  lavaan_lowSATT$Condition <- "low"
  resultSATT_i <- rbind(lavaan_wideSATT,lavaan_highSATT,lavaan_lowSATT)
  #results
  
  result_i<-rbind(results_wide,results_high,results_low)
  result_i <- bind_rows(result_i,resultSATT_i)
  result_i$seed <- i
  result_i$items <- items
  result_i$effectsize <- mD1
  result_i$changevar <- varD
  result_i$N <- N
  
  result[[i]]=result_i
  print(i)
}

result_matrix <- bind_rows(result)
return(result_matrix)

}

results=SimulmeanScore(N=N,B=B,items=items,mD0=mD0,mD1=mD1,varD=varD,covar=covar)
