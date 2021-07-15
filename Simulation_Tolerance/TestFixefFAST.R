# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier
# ---------------------------------------------------------------------

TestFixef <- function(model, model0, Data, id, Time, method, B, level){
  
  
  #Object of class varComprob.S (varComprob with S)
  if(class(model)[1]=="varComprob.S"){
    if(method == "parametric"){
      source("TestFixef_paramMLS.R")
      result = parametric_S(model = model, model0 = model0, id = id, Time = Time, B = B, level = level)
      result = list(Estimates=result$`estimation`,p_value=result[[2]])
      return(result)
    }
    if(method == "wild"){
      source("TestFixef_wildMLS.R")
      result = wild_REMLvarComp(model = model, model0 = model0, id = id, Time = Time, B = B, level = level)
      result = list(Estimates=result$`estimation`,p_value=result[[2]])
      return(result) 
    }
    if(method == "Wald"){
      summ = summary(model)
      alpha = 1 - level
      inf = c(model$beta - summ$zTable[,2]*qnorm(level+alpha/2),model$eta - sqrt(diag(model$vcov.eta))*qnorm(level+alpha/2))
      sup = c(model$beta + summ$zTable[,2]*qnorm(level+alpha/2),model$eta + sqrt(diag(model$vcov.eta))*qnorm(level+alpha/2))
      result = cbind(inf,sup)
      colnames(result) = c("lower bound", "upper bound")
      return(result)
    }else{
      print("Error! Probably an incorrect value for method argument")
    }
  }
  
  #Object of class varComprob.compositeTau (varComprob with cTAU)
  if(class(model)[1]=="varComprob.compositeTau"){
    if(method == "parametric"){
      source("TestFixef_paramMLcTAU.R")
      result = parametric_cTAU(model = model, model0 = model0, id = id, Time = Time, B = B, level = level)
      result = list(Estimates=result$`estimation`,p_value=result[[2]])
      return(result)
    }
    if(method == "wild"){
      source("TestFixef_wildMLcTAU.R")
      result = wild_REMLvarComp(model = model, model0 = model0, id = id, Time = Time, B = B, level = level)
      result = list(Estimates=result$`estimation`,p_value=result[[2]])
      return(result)
    }
    if(method == "Wald"){
      summ = summary(model)
      alpha = 1 - level
      inf = c(model$beta - summ$zTable[,2]*qnorm(level+alpha/2),model$eta - sqrt(diag(model$vcov.eta))*qnorm(level+alpha/2))
      sup = c(model$beta + summ$zTable[,2]*qnorm(level+alpha/2),model$eta + sqrt(diag(model$vcov.eta))*qnorm(level+alpha/2))
      result = cbind(inf,sup)
      colnames(result) = c("lower bound", "upper bound")
      return(result)
    }else{
      print("Error! Probably an incorrect value for method argument")
    }
  }
  
  
  #Object of class lmerMod (lmer)
  if(class(model)[1]=="lmerModLmerTest"){
    if(method == "parametric"){
      if(model@resp$REML == 0){
        source("TestFixef_paramML.R")
        result = param_lmer(model = model, model0 = model0, B = B, level = level)
        result = list(Estimates=result$`estimation`,p_value=result[[2]])
        return(result)
      }else{
        source("TestFixef_paramML.R")
        result = param_lmer(model = model, model0 = model0, B = B, level = level)
        result = list(Estimates=result$`estimation`,p_value=result[[2]])
        return(result)
      }
    }
    if(method == "wild"){
      if(model@resp$REML == 0){
        source("TestFixef_wildML.R")
        result = wild_lmer(model = model, model0 = model0, B = B, level = level)
        result = list(Estimates=result$`estimation`,p_value=result[[2]])
        return(result)
      }else{
        source("TestFixef_wildML.R")
        result = wild_lmer(model = model, model0 = model0, B = B, level = level)
        result = list(Estimates=result$`estimation`,p_value=result[[2]])
        return(result)
      }
    }
    if(method == "Wald"){
      summ=summary(model)
      alpha = 1 - level
      inf = c(fixef(model) - summ$coefficients[,2]*qnorm(level+alpha/2))
      sup = c(fixef(model) + summ$coefficients[,2]*qnorm(level+alpha/2))
      result = cbind(inf,sup)
      colnames(result) = c("lower bound", "upper bound")
      return(result)
    }else{
      print("Error! Probably an incorrect value for method argument")
    }
  }
  
  
  #Object of class rlmerMod (rlmer)
  if(class(model)[1]=="rlmerMod"){
    if(method == "parametric"){
      source("TestFixef_paramDAStau.R")
      result = param_lmer(model = model, model0 = model0, B = B, level = level)
      result = list(Estimates=result$`estimation`,p_value=result[[2]])
      return(result)
    }
    if(method == "wild"){
      source("TestFixef_wildML.R")
      result = wild_lmer(model = model, model0 = model0, B = B, level = level)
      result = list(Estimates=result$`estimation`,p_value=result[[2]])
      return(result)
    }
    if(method == "Wald"){
      summ = summary(model)
      alpha = 1 - level
      inf = c(fixef(model) - summ$coefficients[,2]*qnorm(level+alpha/2))
      sup = c(fixef(model) + summ$coefficients[,2]*qnorm(level+alpha/2))
      result = cbind(inf,sup)
      colnames(result) = c("lower bound", "upper bound")
      return(result)
    }else{
      print("Error! Probably an incorrect value for method argument")
    }
  }
  
  
  #Object of class heavyLme
  if(class(model)[1]=="heavyLme"){
    if(method == "parametric"){
      source("TestFixef_paramNheavyLme.R")
      result = paramN_heavyLme(model = model, model0 = model0, Data = Data, B = B, level = level)
      result = list(Estimates=result$`estimation`,p_value=result[[2]])
      return(result)
    }
    if(method == "wild"){
      source("TestFixef_wildMLheavyLme.R")
      result = wild_heavyLme(model = model, model0 = model0, Data = Data, B = B, level = level)
      result = list(Estimates=result$`estimation`,p_value=result[[2]])
      return(result)
    }
    if(method == "Wald"){
      summ = summary(model)
      alpha = 1 - level
      inf = c(coefficients(model) - summ$coefficients[,2]*qnorm(level+alpha/2))
      sup = c(coefficients(model) + summ$coefficients[,2]*qnorm(level+alpha/2))
      result = cbind(inf,sup)
      colnames(result) = c("lower bound", "upper bound")
      return(result)
    }else{
      print("Error! Probably an incorrect value for method argument")
    }
  }
}

