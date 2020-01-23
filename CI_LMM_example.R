# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the sleepstudy example and can be adapted to your own balanced dataset
# From Mason, F., Cantoni, R., & Ghisletta, P. (submitted). Robust estimation and confidence intervals in linear mixed models.
# ---------------------------------------------------------------------

###############
# preliminary setup

# 1) Set working directory (select the folder "CIfunctions")
setwd("..../CIfunctions")

# 2) Load the packages
library(MASS) # version 7.3-51.1
library(robustvarComp) # version 0.1-2
library(robustlmm) # version 2.3
library(heavy) # version 0.38.19
library(lme4) # version 1.1-20
library(doParallel) # version 1.0.14
source("confintLMM.R") # function to produce confidence intervals

# 3) Import balanced dataset
sleepstudy

# Identify dataset, time and participant variables
Dataset = sleepstudy
time = sleepstudy$Days
participant = sleepstudy$Subject

###############
# LMM estimations

# 4a) Estimation with varComprob() (see corresponding helpfile for more details)

# Build the argument "groups" of the varComprob() function
n = length(unique(participant)) # the number of participants
J = length(unique(time)) # the number of repeated observations per participant
groups = cbind(rep(1:J, each=n),rep((1:n), J)) # a numeric matrix with two columns used to group the observations according to participant.

# Build the argument "varcov" of the varComprob() function
z1 = rep(1, J) #Value for intercept (=1) for the J observations by clusters
z2 = unique(time) # Value for the time variable

K = list() # the "varcov" object
K[[1]] = tcrossprod(z1,z1) # Matrix for intercept
K[[2]] = tcrossprod(z2,z2) # Matrix for time variable
K[[3]] = tcrossprod(z1,z2) + tcrossprod(z2,z1) # Matrix of interaction Intercept by time variable
names(K) = c("sigma2_Intercept", "sigma2_Time", "Covariance")

# Define the formula of the model
model.formula = Reaction ~ 1 + Days

# Estimation with S-estimator
model.S = varComprob(model.formula, groups = groups, data = sleepstudy, varcov = K, control = varComprob.control(lower = c(0,0,-Inf), method = "S", psi = "rocke")) 

# Estimation with composite-TAU estimator
model.cTAU = varComprob(model.formula, groups = groups,data = sleepstudy, varcov = K, control = varComprob.control(lower = c(0, 0, -Inf))) 

# 4b) Estimation with rlmer() (see corresponding helpfile for more details)

# Estimation with SMDM
model.SMDM = rlmer(Reaction ~ 1 + Days + (Days|Subject), data = sleepstudy, rho.sigma.e = psi2propII(smoothPsi, k = 2.28), rho.sigma.b = chgDefaults(smoothPsi, k = 5.11, s = 10))                   

# Estimation with SMDMvar
model.SMDMvar = rlmer(Reaction ~ 1 + Days + (Days|Subject), data = sleepstudy, method = "DASvar", rho.sigma.e = psi2propII(smoothPsi, k = 2.28), rho.sigma.b = chgDefaults(smoothPsi, k = 5.11, s = 10))                   

# 4c) Estimation with heavyLme() (see corresponding helpfile for more details)

# Estimation with the multivarite t-ML
model.tML = heavyLme(Reaction ~ 1 + Days, random = ~ 1 + Days, groups = ~ Subject, data = sleepstudy) 

# 4d) Estimation with lmer() (see corresponding helpfile for more details)

# Estimation with ML
model.ML = lmer(Reaction ~ 1 + Days + (Days|Subject), data = sleepstudy, REML = F)

# Estimation with REML 
model.REML = lmer(Reaction ~ 1 + Days + (Days|Subject), data = sleepstudy)

###############
# Confidence interval estimations

# 5) Confidence Intervals
confint.LMM(model = model.ML, Data = Dataset, id = participant, Time = time, method = "parametric", B = 999, level = .95)

# ARGUMENTS
# model: an object of class varComprob (or varComprob.fit or varComprob.S), lmerMod, rlmerMod or heavyLme
# Data: The data.frame object containing the data
# id: The grouping variable vector
# Time: the time variable vector
# method: "parametric" for classical Normal-parametric bootstrap percentile CI, "t-parametric" for Student-parametric bootstrap percentile CI (available with object of class heavyLme),
# "wild" for the wild  bootstrap percentile CI and "Wald" for z-Wald CI.
# B: number of bootstrap samples, positive integer (with "parametric", "t-parametric" or "wild") > 1
# level: confidence level < 1

# VALUE
# A matrix holding columns for point estimates and, lower and upper boudaries confidence intervals for each parameter.
