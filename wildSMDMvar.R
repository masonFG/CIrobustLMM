# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the sleepstudy example with the wild bootstrap on SMDMvar
# ---------------------------------------------------------------------

# 1) Load the packages and functions
library(lme4) # version 1.1-20
library(robustlmm) # version 2.3
library(MASS) # version 7.3-51.1
library(doParallel) # version 1.0.14
source("CIfunction_wildREML.R")

# 2) Import balanced dataset
sleepstudy

# 3) Estimation with rlmer() (see corresponding helpfile for more details)
model.SMDMvar = rlmer(Reaction ~ 1 + Days + (Days|Subject), data = sleepstudy, method = "DASvar", rho.sigma.e = psi2propII(smoothPsi, k = 2.28), rho.sigma.b = chgDefaults(smoothPsi, k = 5.11, s = 10))                   

# 4) Wald-z 95% Confidence Intervals
summ=summary(model.SMDMvar)
Wald_CI.SMDMvar = t(matrix(c(fixef(model.SMDMvar)[1] - summ$coefficients[1,2]*qnorm(.975), fixef(model.SMDMvar)[1] + summ$coefficients[1,2]*qnorm(.975),
                             fixef(model.SMDMvar)[2] - summ$coefficients[2,2]*qnorm(.975), fixef(model.SMDMvar)[2] + summ$coefficients[2,2]*qnorm(.975)), 2, 2,
                             dimnames = list(c("lower bound", "upper bound"), c("Intercept", "Time"))))

# 5) Percentile Confidence Intervals wih the parametric bootstrap
wild_lmer(model = model.SMDMvar, B = 999, level = .95)

# ARGUMENTS
# model: an object of class rlmerMod 
# B: number of bootstrap samples, positive integer
# level: confidence level < 1

# VALUE
# A matrix with columns giving lower and upper confidence limits for each parameter.