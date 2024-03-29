CIrobustLMM

Contains a set of R functions and scripts to compute percentile confidence intervals based on bootstrap schemes, for linear mixed models estimated with both classical and robust methods. The scripts allow for balanced longitudinal data designs (and unbalanced usign rlmer (e.g. Koller, 2016) associated to the wild bootstrap). Two bootstrap schemes are implemented: the wild bootstrap and the parametric bootstrap. Both are adapted for objects of class merMod (from package lme4), varComprob (or varComprob.fit or varComprob.S from package robustvarComp) and rlmerMod (from package robustlmm). A fully worked out example of longitudinal linear mixed model on the sleepstudy data is provided.

Contents

The main function is confint.lmm(). 

The file "exampleSLEEP.R" contains an example of its application on the sleepstudy data. It includes estimation with lmer() for methods ML and REML (e.g. Bates, 2014), rlmer() with options "DAStau" and "DASvar" and varComprob() for methods S and composite-tau (see Agostinelli & Yohai, 2016). Users can adapt this document to their own dataset.
The folder CIfunctions contains the auxiliary functions files, which must be downloaded and sourced before using the function confint.lmm().  More specifically, the file "confintLMM.R" contains the function confint.lmm(), which calls others functions contained in the folder (the files called by the prefix "CIfunctions_") depending on which estimation method is selected and which type of confidence interval is chosen.

The folder "sleepstudySimulationGithub" contains all documents to replicate the simulation in Parametric and Semi-Parametric Bootstrap Based Confidence Intervals for Robust Linear Mixed Models (Mason et al., 2021).


  
  
 
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

Mason, F., Cantoni, E., & Ghisletta, P. (2021). Parametric and Semi-Parametric Bootstrap-Based Confidence Intervals for Robust Linear Mixed Models. 
Methodology, 17(4), 271-295. doi: https://doi.org/10.5964/meth.6607

Modugno, L., & Giannerini, S. (2013). The wild bootstrap for multilevel models.
Communications in Statistics - Theory and Methods, 44 (22), 4812–4825. doi:
10.1080/03610926.2013.802807


  

