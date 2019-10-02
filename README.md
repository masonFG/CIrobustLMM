CIrobustLMM
 
 Contains R scripts to calculate percentile confidence intervals based on bootstrap schemes or with Wald-z for robust and classical linear mixed models. Two bootstrap schemes are implemented : the wild bootstrap (Giannerini & Modugno, 2013) and the parametric bootstrap. Both of them are adapted for objects of class merMod (from package lme4), varComprob (or varComprob.fit or varComprob.S from package robustvarComp), heavyLme (from package heavy) and rlmerMod (from package robustlmm). An example of longitudinal linear mixed model on the sleepstudy data is provided.

Contents
 
 The current repository contains a folder with the functions files (named with prefix "CIfunction_") which must be downloaded to execute 
  
  "Functions":
  - CIfunction_paramML.R: the function of the parametric bootstrap with ML estimation (via lmer from lme4) on the bootstrap samples. Available for objects of class merMod and rlmerMod. 
  - CIfunction_paramREML.R: the function of the parametric bootstrap with REML estimation (via lmer from lme4) on the bootstrap samples. Available for objects of class merMod and rlmerMod. 
  - CIfunction_paramcTAU.R: the function of the parametric bootstrap with composite-TAU robust estimation of Agostinelli and Yohai (2016) (via varComprob from robustvarComp) on the bootstrap samples. Available for objects of class varComprob (or varComprob.fit or varComprob.S).
  - CIfunction_paramS.R: the function of the parametric bootstrap with the robust S estimation of Copt and Victoria-Feser (2006) (via varComprob from robustvarComp) on the bootstrap samples. Available for objects of class varComprob (or varComprob.fit or varComprob.S).
  - CIfunction_paramNheavyLme.R: the function of the Normal-parametric bootstrap with the robust multivariate t-ML estimation of Pinheiro et al. (2001) (via heavyLme from heavy) on the bootstrap samples. Available for objects of class heavyLme.
  - CIfunction_paramStheavyLme.R: the function of a Student-parametric bootstrap (residuals and random effect are generated from t distributions) with the robust multivariate t-ML estimation of Pinheiro et al. (2001) (via heavyLme from heavy) on the bootstrap samples. Available for objects of class heavyLme.
  - CIfunction_wildML.R: the function of the wild bootstrap with ML estimation (via lmer from lme4) on the bootstrap samples. Available for objects of class merMod and rlmerMod. 
  - CIfunction_wildREML.R: the function of the wild bootstrap with REML estimation (via lmer from lme4) on the bootstrap samples. Available for objects of class merMod and rlmerMod. 
  - CIfunction_wildcTAU.R: the function of the wild bootstrap with composite-TAU robust estimation of Agostinelli and Yohai (2016) (via varComprob from robustvarComp) on the bootstrap samples. Available for objects of class varComprob (or varComprob.fit or varComprob.S).
  - CIfunction_wildS.R: the function of the wild bootstrap with the robust S estimation of Copt and Victoria-Feser (2006) (via varComprob from robustvarComp) on the bootstrap samples. Available for objects of class varComprob (or varComprob.fit or varComprob.S).
  - CIfunction_wildheavyLme.R: the function of the wild bootstrap with the robust multivariate t-ML estimation of Pinheiro et al. (2001) (via heavyLme from heavy) on the bootstrap samples. Available for objects of class heavyLme.
  
  "Running":
  - paramML.R: Example of parametric bootstrap CI (from CIfunction_paramML.R) applied on the sleepstudy data fitted with ML method.
  - paramREML.R: Example of parametric bootstrap CI (from CIfunction_paramREML.R) applied on the sleepstudy data fitted with REML method.
  - paramcTAU.R: Example of parametric bootstrap CI (from CIfunction_paramcTAU.R) applied on the sleepstudy data fitted with robust composite-TAU method.
  - paramS.R: Example of parametric bootstrap CI (from CIfunction_paramcTAU.R) applied on the sleepstudy data fitted with robust  S method.
  - paramSMDM.R: Example of parametric bootstrap CI (from CIfunction_paramREML.R) applied on the sleepstudy data fitted with robust SMDM method (Koller, 2013).
  - paramSMDMvar.R: Example of parametric bootstrap CI (from CIfunction_paramREML.R) applied on the sleepstudy data fitted with the robust approximated version of SMDM method (Koller, 2013).
  - paramN_tML.R: Example of Normal-parametric bootstrap CI (from CIfunction_paramNheavyLme.R) applied on the sleepstudy data fitted with the robust multivariate t-ML method (Pinheiro et al., 2001).
  - paramSt_tML.R: Example of Student-parametric bootstrap CI (from CIfunction_paramStheavyLme.R) applied on the sleepstudy data fitted with the robust multivariate t-ML method (Pinheiro et al., 2001).
  - wildML.R: Example of wild bootstrap CI (from CIfunction_wildML.R) applied on the sleepstudy data fitted with ML method
  - wildREML.R: Example of wild bootstrap CI (from CIfunction_wildREML.R) applied on the sleepstudy data fitted with REML method.
  - wildcTAU.R: Example of wild bootstrap CI (from CIfunction_wildcTAU.R) applied on the sleepstudy data fitted with robust composite-TAU method.
  - wildS.R: Example of wild bootstrap CI (from CIfunction_wildcTAU.R) applied on the sleepstudy data fitted with robust  S method.
  - wildSMDM.R: Example of wild bootstrap CI (from CIfunction_wildREML.R) applied on the sleepstudy data fitted with robust SMDM method (Koller, 2013).
  - wildSMDMvar.R: Example of wild bootstrap CI (from CIfunction_wildREML.R) applied on the sleepstudy data fitted with the robust approximated version of SMDM method (Koller, 2013).
  - wild_tML.R: Example of wild bootstrap CI (from CIfunction_wildheavyLme.R) applied on the sleepstudy data fitted with the robust multivariate t-ML method (Pinheiro et al., 2001).
 
References

Agostinelli, C., & Yohai, V. J. (2016). Composite robust estimators for linear mixed
models. Journal of the American Statistical Association, 111 (516), 1764–1774.

Bates, D. e. a. (2014). lme4: Linear mixed-effects models using eigen and s4. R package
version, 1 (7), 1–23.

Copt, S., & Victoria-Feser, M.-P. (2006). High-breakdown inference for mixed linear
models. Journal of the American Statistical Association, 101 (473), 292–300.

Koller, M. (2013). Robust estimation of linear mixed models (Doctoral dissertation, ETH
Zurich). doi: 10.3929/ethz-a-007632241

Koller, M. (2016). robustlmm: An r package for robust estimation of linear mixed-effects
models. Journal of Statistical Software, 75 (6). Retrieved from
http://www.jstatsoft.org/v75/i06/ doi: 10.18637/jss.v075.i06

Modugno, L., & Giannerini, S. (2013, Nov). The wild bootstrap for multilevel models.
Communications in Statistics - Theory and Methods, 44 (22), 4812–4825. doi:
10.1080/03610926.2013.802807

Osorio, F. (2016). heavy: Robust estimation using heavy-tailed distributions. R package
version 0.38. ed.

Pinheiro, J. C., Liu, C., & Wu, Y. N. (2001, Jun). Efficient algorithms for robust
estimation in linear mixed-effects models using the multivariate t distribution.
Journal of Computational and Graphical Statistics, 10 (2), 249–276. doi:
10.1198/10618600152628059
  

