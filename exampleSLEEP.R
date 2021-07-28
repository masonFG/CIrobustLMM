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
library(Rcpp)

# 3) Import balanced dataset
sleepstudy

# Identify dataset and participant variables
Dataset = sleepstudy
participant = sleepstudy$Subject

###############
# LMM estimations

# 4a) Estimation with varComprob() (see corresponding helpfile for more details)

# Define the within-subject variable (by example the time variable)
WITHIN = sleepstudy$Days

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

# Define the formula of the model
model.formula = Reaction ~ 1 + Days

# Estimation with S-estimator
model.S = varComprob(model.formula, groups = groups, data = sleepstudy, varcov = K, control = varComprob.control(lower = c(0,0,-Inf), method = "S", psi = "rocke"))

# Estimation with composite-TAU estimator
model.cTAU = varComprob(model.formula, groups = groups, data = sleepstudy, varcov = K, control = varComprob.control(lower = c(0, 0, -Inf)))

# 4b) Estimation with rlmer() (see corresponding helpfile for more details)

# Estimation with DAStau
model.DAStau = rlmer(Reaction ~ 1 + Days + (Days|Subject), data = sleepstudy, rho.sigma.e = psi2propII(smoothPsi, k = 2.28), rho.sigma.b = chgDefaults(smoothPsi, k = 5.11, s = 10))

# Estimation with DASvar
model.DASvar = rlmer(Reaction ~ 1 + Days + (Days|Subject), data = sleepstudy, method = "DASvar", rho.sigma.e = psi2propII(smoothPsi, k = 2.28), rho.sigma.b = chgDefaults(smoothPsi, k = 5.11, s = 10))

# 4c) Estimation with lmer() (see corresponding helpfile for more details)

# Estimation with ML
model.ML = lmer(Reaction ~ 1 + Days + (Days|Subject), data = sleepstudy, REML = F)

# Estimation with REML
model.REML = lmer(Reaction ~ 1 + Days + (Days|Subject), data = sleepstudy)

# 5) Confidence Intervals
model = model.S; data = Dataset; clusterid = participant; methodCI = "wild"; B = 10; confint.level = .95; BCa = T
source("confintLMMfast.R")
source("BCaML.R")
CI <- BCaboot(model = model.ML, data = Dataset, clusterid = participant, methodCI = "wild", B = 10, confint.level = .95, BCa = T)
CI$Percentile
CI[[1]]$BCa.interval[[1]]

source("BCaVarCompRob.R")
CIv <- BCabootvarCompRob(model = model.cTAU, data = Dataset, clusterid = participant, methodCI = "wild", B = 100, confint.level = .95, BCa = T)
CIv$Percentile
CIv[[1]]$BCa.interval[[1]]
