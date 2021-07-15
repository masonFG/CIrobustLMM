# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: August, 2020
# R version: 3.6.0
#  R code for analyses presented in

#Mason, F., Cantoni, E., & Ghisletta, P. (submitted). Robust estimation and confidence intervals in linear mixed models.

#This script contains the code to reproduce the sleepstudy example and can be adapted to your own balanced dataset.
# ---------------------------------------------------------------------

###############
# preliminary setup

# 1) Set working directory (select the folder "CIfunctions")
#setwd("..../CIfunctions")
 ML = "ML"
 tML = "tML"
 DAStau = "DAStau"
 S = "S"
 cTAU = "cTAU"
 parametric = "parametric"
 wild="wild"

eval(parse(text=paste(commandArgs(trailingOnly = TRUE), collapse="")))

simname = paste0("res", paste(commandArgs(trailingOnly = TRUE), collapse=""), ".RData")
print(simname)

# 2) Load the packages
library(MASS) # version 7.3-51.1
library(robustvarComp) # version 0.1-2
library(robustlmm) # version 2.3
#library(heavy) # version 0.38.19
library(lme4) # version 1.1-20
library(lmerTest) # version 3.1-2
#library(doParallel) # version 1.0.14
library(stringr)  # version 1.4.0
source("confintLMMFast.R") # function to produce confidence intervals
source("TestFixefFAST.R") # function to test Fixed Effect with bootstrap
source("BCaML.R") # function to produce BCa confidence intervals with lmer() and rlmer()
source("BCaVarCompRob.R") # function to produce BCa confidence intervals with varComprob()
source("BCaHeavy.R") # function to produce BCa confidence intervals with heavyLme()
source("BAOBABEtude2.R") # function for baobab
source("CIfunction_paramNheavyLme.R")



# 3) Import balanced dataset
#baobab
options(echo=TRUE) 
i <- commandArgs(trailingOnly = TRUE)[10]
 i=str_split(i,"=")
rseed=as.numeric(i[[1]][2])
print(rseed)

res0=Simul2(NTot = NTot,
                   extp1 = extp1,
                   extRANI1 = extRANI1,
                   extRANP1 = extRANP1,
                   bGxT = bGxT,
                   typeExt = typeExt,
                   estimator = estimator,
                   bootstrap = bootstrap,
                   test = test,
                   rseed = rseed)


nameparam = paste0("estimator",estimator,"bootstrap",bootstrap,"NTot",NTot, "extr.resid",extp1,"extr.ranS",extRANP1,"interaction",bGxT,"typeExt",typeExt,"test",test,"rseed",rseed)

save(res0,file=paste("output/",nameparam,".RData",sep=""))
