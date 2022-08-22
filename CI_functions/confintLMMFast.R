# ---------------------------------------------------------------------
# Author: Fabio Mason
# Date: September, 2019
# R version: 3.5.1
# Comment: R code for analyses presented in titres-du-papier
# ---------------------------------------------------------------------

#' Title
#'
#' @param model 
#' @param Data 
#' @param id 
#' @param method 
#' @param B 
#' @param level 
#'
#' @return
#' @export
#'
#' @examples
confint.LMM <- function(model, Data, id, method, B, level){
  
  
  #Object of class varComprob.S (varComprob with S)
  if(class(model)[1]=="varComprob.S"){
    if(method == "parametric"){
      source("CIfunction_paramMLvarCompRob.R")
      result = parametric_MLvarCompRob(model = model, id = id, B = B, level = level)
      return(result)
    }
    if(method == "wild"){
      source("CIfunction_wildvarCompRob.R")
      result = wild_MLvarCompRob(model = model, id = id, B = B, level = level)
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
      source("CIfunction_paramMLvarCompRob.R")
      result = parametric_MLvarCompRob(model = model, id = id, B = B, level = level)
      return(result)
    }
    if(method == "wild"){
      source("CIfunction_wildvarCompRob.R")
      result = wild_MLvarCompRob(model = model, id = id, B = B, level = level)
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
  if(class(model)[1]=="lmerModLmerTest" | class(model)[1]=="lmerMod"){
    if(method == "parametric"){
      if(model@resp$REML == 0){
        source("CIfunction_paramREML.R")
        result = param_lmer(model = model, B = B, level = level)
        return(result)
      }else{
        source("CIfunction_paramREML.R")
        result = param_lmer(model = model, B = B, level = level)
        return(result)
      }
    }
    if(method == "wild"){
      if(model@resp$REML == 0){
        source("CIfunction_wildREML.R")
        result = wild_lmer(model = model, B = B, level = level)
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
      source("CIfunction_paramREML.R")
      result = param_lmer(model = model, B = B, level = level)
      return(result)
    }
    if(method == "wild"){
      source("CIfunction_wildDAStau.R")
      result = wild_rlmer(model = model, B = B, level = level)
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
  
}

