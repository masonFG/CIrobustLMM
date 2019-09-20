# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier - this script contains the sleepstudy example with the wild bootstrap on REML
# ---------------------------------------------------------------------

# 1) Load the packages and functions
library(lme4) # version 1.1-20
library(MASS) # version 7.3-51.1
source("CIfunction_wildREML.R")

# 2) Import balanced dataset
sleepstudy

# 3) Estimation with lmer() (see corresponding helpfile for more details)
model.REML = lmer(Reaction ~ 1 + Days + (Days|Subject), data = sleepstudy)                   

# 4) Percentile Confidence Intervals wih the wild bootstrap
wild_lmer(model = model.REML, B = 999, level = .95)

# ARGUMENTS
# model: an object of class lmerMod
# B: number of bootstrap samples, positive integer
# level: confidence level < 1

# VALUE
# A matrix with columns giving lower and upper confidence limits for each parameter.