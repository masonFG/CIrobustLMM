# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the sleepstudy example with the Normal-parametric bootstrap on multivariate t-ML
# ---------------------------------------------------------------------

# 1) Load the packages and functions
library(lme4) # version 1.1-20
library(heavy) # version 0.38.19
library(MASS) # version 7.3-51.1
library(stringr) # version 1.3.1
source("CIfunction_paramNheavyLme.R")

# 2) Import balanced dataset
sleepstudy

# 3) Estimation with heavyLme() (see corresponding helpfile for more details)
model.tML = heavyLme(Reaction ~ 1 + Days, random = ~ 1 + Days, groups = ~ Subject, data = sleepstudy)                   

# 4) Percentile Confidence Intervals wih the normal-parametric bootstrap
paramN_heavyLme(model = model.tML, Data = sleepstudy , B = 999, level = .95)

# ARGUMENTS
# model : An object of class heavyLme representing the linear mixed-effects model fit
# B: number of bootstrap samples, positive integer
# level: confidence level < 1

# VALUE
# A matrix with columns giving lower and upper confidence limits for each parameter.