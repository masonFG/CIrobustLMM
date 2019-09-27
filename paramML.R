# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the sleepstudy example with the parametric bootstrap on ML
# ---------------------------------------------------------------------

# 1) Load the packages and functions
library(lme4) # version 1.1-20
library(MASS) # version 7.3-51.1
library(doParallel) # version 1.0.14
source("CIfunction_paramML.R")

# 2) Import balanced dataset
sleepstudy

# 3) Estimation with lmer() (see corresponding helpfile for more details)
model.ML = lmer(Reaction ~ 1 + Days + (Days|Subject), data = sleepstudy, REML = F)

# 4) Wald-z 95% Confidence Intervals
summ=summary(model.ML)
Wald_CI.ML = t(matrix(c(fixef(model.ML)[1] - summ$coefficients[1,2]*qnorm(.975), fixef(model.ML)[1] + summ$coefficients[1,2]*qnorm(.975),
                        fixef(model.ML)[2] - summ$coefficients[2,2]*qnorm(.975), fixef(model.ML)[2] + summ$coefficients[2,2]*qnorm(.975)), 2, 2,
                        dimnames = list(c("lower bound", "upper bound"), c("Intercept", "Time"))))

# 5) Percentile Confidence Intervals wih the parametric bootstrap
param_lmerML(model = model.ML, B = 999, level = .95)

# ARGUMENTS
# model: an object of class lmerMod
# B: number of bootstrap samples, positive integer
# level: confidence level < 1

# VALUE
# A matrix with columns giving lower and upper confidence limits for each parameter.