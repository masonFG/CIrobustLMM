# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the sleepstudy example with the parametric bootstrap on composite TAU-Estimator
# ---------------------------------------------------------------------

# 1) Load the packages and functions
library(MASS) # version 7.3-51.1
library(robustvarComp) # version 0.1-2
library(lme4) # version 1.1-20
source("CIfunction_paramcTAU.R")

# 2) Import balanced dataset and define time and participant variable
sleepstudy

# Identify time and participant variables
time = sleepstudy$Days
participant = sleepstudy$Subject

# 3) Estimation with varComprob() (data must be sorted by cluster, see corresponding helpfile for more details)

# Build the argument "groups" of the varComprob() function
n = length(unique(participant)) # the number of participants
J = length(unique(time)) # the number of observations per participant
groups = cbind(rep(1:J, each=n),rep((1:n), J)) # a numeric matrix with two columns used to group the observations according to participant.

# Build the argument "varcov" of the varComprob() function
z1 = rep(1, J) # VALUE for intercept (=1) for the J observations by clusters
z2 = unique(time) # VALUE for the time variable

K = list() # the "varcov" object
K[[1]] = tcrossprod(z1,z1) # Matrix for intercept
K[[2]] = tcrossprod(z2,z2) # Matrix for time variable
K[[3]] = tcrossprod(z1,z2) + tcrossprod(z2,z1) # Matrix of interaction Intercept by time variable
names(K) = c("sigma2_Intercept", "sigma2_Time", "Covariance")

# Define the formula of the model
model.formula = Reaction ~ 1 + Days
model.cTAU = varComprob(model.formula, groups = groups,data = sleepstudy, varcov = K, control = varComprob.control(lower = c(0, 0, -Inf))) # Estimation with the composite TAU-estimator

# 4) Wald-z 95% Confidence Intervals
summ = summary(model.cTAU)
Wald_CI.cTAU = matrix(c(model.cTAU$beta[1] - summ$zTable[1,2]*qnorm(.975), model.cTAU$beta[1] + summ$zTable[1,2]*qnorm(.975),
                        model.cTAU$beta[2] - summ$zTable[2,2]*qnorm(.975), model.cTAU$beta[2] + summ$zTable[2,2]*qnorm(.975),
                        model.cTAU$eta[1] - sqrt(diag(model.cTAU$vcov.eta))[1]*qnorm(.975), model.cTAU$eta[1] + sqrt(diag(model.cTAU$vcov.eta))[1]*qnorm(.975),
                        model.cTAU$eta[2] - sqrt(diag(model.cTAU$vcov.eta))[2]*qnorm(.975), model.cTAU$eta[2] + sqrt(diag(model.cTAU$vcov.eta))[2]*qnorm(.975),
                        model.cTAU$eta[3] - sqrt(diag(model.cTAU$vcov.eta))[3]*qnorm(.975), model.cTAU$eta[3] + sqrt(diag(model.cTAU$vcov.eta))[3]*qnorm(.975)), 5, 2,
                        dimnames = list(c("Intercept", "Time", "Sigma2_intercept", "Sigma2_time", "Covariance"), c("lower bound", "upper bound")))

# 5) Percentile Confidence Intervals wih the parametric bootstrap
parametric_cTAU(model = model.cTAU, Time = time, B = 999, level = .95)

# ARGUMENTS
# model: an object of class varComprob or varComprob.fit or varComprob.S
# Time: the time variable vector
# B: number of bootstrap samples, positive integer
# level: confidence level < 1

# VALUE
# A matrix with columns giving lower and upper confidence limits for each parameter.