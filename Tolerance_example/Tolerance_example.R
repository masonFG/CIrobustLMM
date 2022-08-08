# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: August, 2021
# R version: 3.6.0
#  R code for analyses presented in

#Mason, F., Cantoni, E., & Ghisletta, P. (submitted). Robust estimation and confidence intervals in linear mixed models.

#This script contains the code to reproduce the sleepstudy example and can be adapted to your own balanced dataset.
# ---------------------------------------------------------------------

###############
# preliminary setup

# 1) Set working directory (select the folder "CIfunctions")
#setwd("..../CIfunctions")

# 2) Load the packages
library(MASS) # version 7.3-51.1
library(robustvarComp) # version 0.1-2
library(robustlmm) # version 2.3
library(heavy) # version 0.38.19
library(lme4) # version 1.1-20
library(lmerTest) # version 3.1-2
library(doParallel) # version 1.0.14
library(stringr)  # version 1.4.0
source("confintLMMFast.R") # function to produce confidence intervals
source("TestFixefFAST.R") # function get p_value based on bootstrap
source("BCaML.R") # function get p_value based on bootstrap
source("BCaVarCompRob.R") # function get p_value based on bootstrap

bdd <- read.csv("H:/These/CIrobustLMM_suite/ToleranceGIT/Dataset.txt")


# Identify dataset, time and participant variables
Dataset = bdd
participant = bdd$id

###############
# LMM estimations

# 4a) Estimation with varComprob() (see corresponding helpfile for more details)

# Define the within-subject variable (by example the time variable)
WITHIN = bdd$time

# Build the argument "groups" of the varComprob() function
n = length(unique(participant)) # the number of participants
J = length(unique(WITHIN)) # the number of repeated observations per participant
groups = cbind(rep(1:J, each=n),rep((1:n), J)) # a numeric matrix with two columns used to group the observations according to participant.

# Build the argument "varcov" of the varComprob() function
z1 = rep(1, J) #Value for intercept (=1) for the J observations by clusters
z2 = unique(WITHIN) # Value for the time variable

K = list() # the "varcov" object
K[[1]] = tcrossprod(z1,z1) # Matrix for intercept
K[[2]] = tcrossprod(z2,z2) # Matrix for time variable
K[[3]] = tcrossprod(z1,z2) + tcrossprod(z2,z1) # Matrix of interaction Intercept by time variable
names(K) = c("sigma2_u0", "sigma2_u1", "Covariance")

# Define the formula of the two nested models
model.formula = y ~ 1 + group*time

# Estimation with S-estimator
model.S  = varComprob(model.formula, groups = groups, data = Dataset, varcov = K, control = varComprob.control(lower = c(0,0,-Inf), method = "S", psi = "rocke")) 

wildS <- BCabootvarCompRob(model = model.S, data = Dataset, clusterid = participant, methodCI = "wild", B = 100, confint.level = .95, BCa = T)
wildS$Percentile
wildS[[1]]$BCa.interval[[1]]

paramS <- BCabootvarCompRob(model = model.S, data = Dataset, clusterid = participant, methodCI = "parametric", B = 100, confint.level = .95, BCa = T)
paramS$Percentile
paramS[[1]]$BCa.interval[[1]]

# Estimation with composite-TAU estimator
model.cTAU  = varComprob(model.formula, groups = groups, data = Dataset, varcov = K, control = varComprob.control(lower = c(0, 0, -Inf))) 

wildcTAU <- BCabootvarCompRob(model = model.cTAU, data = Dataset, clusterid = participant, methodCI = "wild", B = 100, confint.level = .95, BCa = T)
wildcTAU$Percentile
wildcTAU[[1]]$BCa.interval[[1]]

paramcTAU <- BCabootvarCompRob(model = model.cTAU, data = Dataset, clusterid = participant, methodCI = "parametric", B = 100, confint.level = .95, BCa = T)
paramcTAU$Percentile
paramcTAU[[1]]$BCa.interval[[1]]

# 4b) Estimation with rlmer() (see corresponding helpfile for more details)

# Estimation with DAStau
model.DAStau = rlmer(y ~ 1 + group*time + (time|id), data = Dataset, rho.sigma.e = psi2propII(smoothPsi, k = 2.28), rho.sigma.b = chgDefaults(smoothPsi, k = 5.11, s = 10))                   

wildDAStau <- BCaboot(model = model.DAStau, data = Dataset, clusterid = participant, methodCI = "wild", B = 10, confint.level = .95, BCa = T)
wildDAStau$Percentile
wildDAStau[[1]]$BCa.interval[[1]]

paramDAStau <- BCaboot(model = model.DAStau, data = Dataset, clusterid = participant, methodCI = "parametric", B = 10, confint.level = .95, BCa = T)
paramDAStau$Percentile
paramDAStau[[1]]$BCa.interval[[1]]



# 4c) Estimation with lmer() (see corresponding helpfile for more details)

# Estimation with ML
model.ML = lmer(y ~ 1 + group*time + (time|id), data = Dataset, REML = F)



wildML <- BCaboot(model = model.ML, data = Dataset, clusterid = participant, methodCI = "wild", B = 10, confint.level = .95, BCa = T)
wildML$Percentile
wildML[[1]]$BCa.interval[[1]]

paramML <- BCaboot(model = model.ML, data = Dataset, clusterid = participant, methodCI = "parametric", B = 10, confint.level = .95, BCa = T)
paramML$Percentile
paramML[[1]]$BCa.interval[[1]]

