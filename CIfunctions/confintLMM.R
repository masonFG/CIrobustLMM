# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier
# ---------------------------------------------------------------------

confint.LMM <- function(model, Data, id, Time, method, B, level){
  
  
  #Object of class varComprob.S (varComprob with S)
  if(class(model)=="varComprob.S"){
    if(method == "parametric"){
      source("CIfunction_paramS.R")
      result = parametric_S(model = model, Time = Time, B = B, level = level)
      return(result)
    }
    if(method == "wild"){
      source("CIfunction_wildS.R")
      result = wild_S(model = model, id = id, Time = Time, B = B, level = level)
      return(result)
    }
    if(method == "Wald"){
      summ = summary(model)
      alpha = 1 - level
      result = t(matrix(c(model.S$beta[1] - summ$zTable[1,2]*qnorm(level+alpha/2), model.S$beta[1] + summ$zTable[1,2]*qnorm(level+alpha/2),
                            model.S$beta[2] - summ$zTable[2,2]*qnorm(level+alpha/2), model.S$beta[2] + summ$zTable[2,2]*qnorm(level+alpha/2),
                            model.S$eta[1] - sqrt(diag(model$vcov.eta))[1]*qnorm(level+alpha/2), model$eta[1] + sqrt(diag(model$vcov.eta))[1]*qnorm(level+alpha/2),
                            model.S$eta[2] - sqrt(diag(model$vcov.eta))[2]*qnorm(level+alpha/2), model$eta[2] + sqrt(diag(model$vcov.eta))[2]*qnorm(level+alpha/2),
                            model.S$eta[3] - sqrt(diag(model$vcov.eta))[3]*qnorm(level+alpha/2), model$eta[3] + sqrt(diag(model$vcov.eta))[3]*qnorm(level+alpha/2)), 2, 5,
                          dimnames = list(c("lower bound", "upper bound"), c("Intercept", "Time", "Sigma2_intercept", "Sigma2_time", "Covariance"))))
      
      return(result)
    }else{
      print("Error! Probably an incorrect value for method argument")
    }
  }
  
  #Object of class varComprob.compositeTau (varComprob with cTAU)
  if(class(model)=="varComprob.compositeTau"){
    if(method == "parametric"){
      source("CIfunction_paramcTAU.R")
      result = parametric_cTAU(model = model, Time = Time, B = B, level = level)
      return(result)
    }
    if(method == "wild"){
      source("CIfunction_wildcTAU.R")
      result = wild_cTAU(model = model, id = id, Time = Time, B = B, level = level)
      return(result)
    }
    if(method == "Wald"){
      summ = summary(model)
      alpha = 1 - level
      result = t(matrix(c(model.S$beta[1] - summ$zTable[1,2]*qnorm(level+alpha/2), model.S$beta[1] + summ$zTable[1,2]*qnorm(level+alpha/2),
                          model.S$beta[2] - summ$zTable[2,2]*qnorm(level+alpha/2), model.S$beta[2] + summ$zTable[2,2]*qnorm(level+alpha/2),
                          model.S$eta[1] - sqrt(diag(model$vcov.eta))[1]*qnorm(level+alpha/2), model$eta[1] + sqrt(diag(model$vcov.eta))[1]*qnorm(level+alpha/2),
                          model.S$eta[2] - sqrt(diag(model$vcov.eta))[2]*qnorm(level+alpha/2), model$eta[2] + sqrt(diag(model$vcov.eta))[2]*qnorm(level+alpha/2),
                          model.S$eta[3] - sqrt(diag(model$vcov.eta))[3]*qnorm(level+alpha/2), model$eta[3] + sqrt(diag(model$vcov.eta))[3]*qnorm(level+alpha/2)), 2, 5,
                        dimnames = list(c("lower bound", "upper bound"), c("Intercept", "Time", "Sigma2_intercept", "Sigma2_time", "Covariance"))))
      
      return(result)
    }else{
      print("Error! Probably an incorrect value for method argument")
    }
  }
  
  
  #Object of class lmerMod (lmer)
  if(class(model)=="lmerMod"){
    if(method == "parametric"){
    if(model@resp$REML == 0){
      source("CIfunction_paramML.R")
      result = param_lmerML(model = model, B = B, level = level)
      return(result)
    }else{
      source("CIfunction_paramREML.R")
      result = param_lmer(model = model, B = B, level = level)
      return(result)
    }
    }
    if(method == "wild"){
      if(model@resp$REML == 0){
        source("CIfunction_wildML.R")
        result = wild_lmerML(model = model, B = B, level = level)
        return(result)
      }else{
        source("CIfunction_wildREML.R")
        result = wild_lmer(model = model, B = B, level = level)
        return(result)
      }
    }
    if(method == "Wald"){
        summ=summary(model)
        alpha = 1 - level
        result = t(matrix(c(fixef(model)[1] - summ$coefficients[1,2]*qnorm(level+alpha/2), fixef(model)[1] + summ$coefficients[1,2]*qnorm(level+alpha/2),
                                fixef(model)[2] - summ$coefficients[2,2]*qnorm(level+alpha/2), fixef(model)[2] + summ$coefficients[2,2]*qnorm(level+alpha/2)), 2, 2,
                              dimnames = list(c("lower bound", "upper bound"), c("Intercept", "Time"))))
        return(result)
    }else{
      print("Error! Probably an incorrect value for method argument")
    }
  }
  
  
  #Object of class rlmerMod (rlmer)
  if(class(model)=="rlmerMod"){
    if(method == "parametric"){
      source("CIfunction_paramREML.R")
      result = param_lmer(model = model, B = B, level = level)
      return(result)
    }
    if(method == "wild"){
      source("CIfunction_wildREML.R")
      result = wild_lmer(model = model, B = B, level = level)
      return(result)
    }
    if(method == "Wald"){
      summ = summary(model)
      alpha = 1 - level
      result = t(matrix(c(fixef(model)[1] - summ$coefficients[1,2]*qnorm(level+alpha/2), fixef(model)[1] + summ$coefficients[1,2]*qnorm(level+alpha/2),
                                fixef(model)[2] - summ$coefficients[2,2]*qnorm(level+alpha/2), fixef(model)[2] + summ$coefficients[2,2]*qnorm(level+alpha/2)), 2, 2,
                              dimnames = list(c("lower bound", "upper bound"), c("Intercept", "Time"))))
      return(result)
    }else{
        print("Error! Probably an incorrect value for method argument")
      }
  }
  
  
  #Object of class heavyLme
  if(class(model)=="heavyLme"){
    if(method == "parametric"){
      source("CIfunction_paramNheavyLme.R")
      result = paramN_heavyLme(model = model, Data = Data, B = B, level = level)
      return(result)
    }
    if(method == "t-parametric"){
      source("CIfunction_paramStheavyLme.R")
      result = paramSt_heavyLme(model = model, Data = Data, B = B, level = level)
      return(result)
    }
    if(method == "wild"){
      source("CIfunction_wildheavyLme.R")
      result = wild_heavyLme(model = model, Data = Data, B = B, level = level)
      return(result)
    }
    if(method == "Wald"){
      summ = summary(model)
      alpha = 1 - level
      result = t(matrix(c(coefficients(model)[1] - summ$coefficients[1,2]*qnorm(level+alpha/2), coefficients(model)[1] + summ$coefficients[1,2]*qnorm(level+alpha/2),
                          coefficients(model)[2] - summ$coefficients[2,2]*qnorm(level+alpha/2), coefficients(model)[2] + summ$coefficients[2,2]*qnorm(level+alpha/2)), 2, 2,
                        dimnames = list(c("lower bound", "upper bound"), c("Intercept", "Time"))))
      
      return(result)
    }else{
      print("Error! Probably an incorrect value for method argument")
    }
  }
}
