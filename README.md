CIrobustLMM

Contains a set of R functions and scripts to compute percentile confidence intervals based on bootstrap schemes or with Wald-z, for linear mixed models estimated with both classical and robust methods. The scripts allow for balanced longitudinal data designs. Two bootstrap schemes are implemented: the wild bootstrap and the parametric bootstrap. Both are adapted for objects of class merMod (from package lme4), varComprob (or varComprob.fit or varComprob.S from package robustvarComp), heavyLme (from package heavy) and rlmerMod (from package robustlmm). A fully worked out example of longitudinal linear mixed model on the sleepstudy data is provided.

Contents

The main function is confint.lmm(). 

The file "CI_LMM_example.R" contains an example of its application on the sleepstudy data. It includes estimation with lmer() for methods ML and REML (e.g. Bates, 2014), rlmer() with options "DAStau" and "DASvar" (e.g. Koller, 2016), varComprob() for methods S and composite-tau (see Agostinelli & Yohai, 2016) and heavyLme() with the multivariate t-ML method (Pinheiro et al., 2001). Users can adapt this document to their own dataset.
The folder CIfunctions contains the auxiliary functions files, which must be downloaded and sourced before using the function confint.lmm().  More specifically, the file "confintLMM.R" contains the function confint.lmm(), which calls others functions contained in the folder (the files called by the prefix "CIfunctions_") depending on which estimation method is selected and which type of confidence interval is chosen.



  
  
 
References

Agostinelli, C., & Yohai, V. J. (2016). Composite robust estimators for linear mixed
models. Journal of the American Statistical Association, 111 (516), 1764–1774.

Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting linear mixed-effects models
using lme4. Journal of Statistical Software, 67 (1), 1–48. doi: 10.18637/jss.v067.i01

Copt, S., & Victoria-Feser, M.-P. (2006). High-breakdown inference for mixed linear
models. Journal of the American Statistical Association, 101 (473), 292–300.

Koller, M. (2013). Robust estimation of linear mixed models (Doctoral dissertation, Eidgenössische Technische Hochschule
Zurich). doi: 10.3929/ethz-a-007632241

Koller, M. (2016). robustlmm: An r package for robust estimation of linear mixed-effects
models. Journal of Statistical Software, 75 (6). doi: 10.18637/jss.v075.i06

Modugno, L., & Giannerini, S. (2013). The wild bootstrap for multilevel models.
Communications in Statistics - Theory and Methods, 44 (22), 4812–4825. doi:
10.1080/03610926.2013.802807

Osorio, & F. (2019). heavy: Robust estimation using heavy-tailed distributions [Computer
software manual]. Retrieved from https://CRAN.R-project.org/package=heavy
(R package version 0.38.196)

Pinheiro, J. C., Liu, C., & Wu, Y. N. (2001). Efficient algorithms for robust
estimation in linear mixed-effects models using the multivariate t distribution.
Journal of Computational and Graphical Statistics, 10 (2), 249–276. doi:
10.1198/10618600152628059
  

