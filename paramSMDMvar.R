# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the sleepstudy example with the parametric bootstrap on SMDMvar
# ---------------------------------------------------------------------

# 1) Load the packages and functions
library(lme4) # version 1.1-20
library(robustlmm) # version 2.3
library(MASS) # version 7.3-51.1
source("CIfunction_paramREML.R")

# 2) Import balanced dataset
sleepstudy

# 3) Estimation with rlmer() (see corresponding helpfile for more details)
model.SMDMvar = rlmer(Reaction ~ 1 + Days + (Days|Subject), data = sleepstudy, method = "DASvar", rho.sigma.e = psi2propII(smoothPsi, k = 2.28), rho.sigma.b = chgDefaults(smoothPsi, k = 5.11, s = 10))                   

# 4) Percentile Confidence Intervals wih the parametric bootstrap
param_lmer(model = model.SMDMvar, B = 999, level = .95)

# ARGUMENTS
# model: an object of class rlmerMod
# B: number of bootstrap samples, positive integer
# level: confidence level < 1

# VALUE
# A matrix with columns giving lower and upper confidence limits for each parameter.