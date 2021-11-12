# ExtremalDependence_Unreplicated
This project is the work done to model extremal dependence in unreplicated soil data.

Script descriptions

Exponential Factor Copula Scripts
1. CensoredLikelihood_nonspatial.R - centralized likelihood function. Calls other likelihood functions. 
2. logCopula_nonspatial.R - called by (1.). Computes likelihood of fully-censored data. 
3. logLikelihood_nonspatial.R - called by (1.). Computes likelihood of non-censored data. 
4. logPartial_nonspatial.R - called by (1.). Computes likelihood of partially-censored data. 
5. Tools.R - Common functions for likelihood computation. 
6. Fit_NonSpatialExpCopulaModel.R - Fits Exponential Factor Copula model to unreplicated data using (1.). Data is found in pairs_uniform.RData.

Multivariate Generalized Pareto Extremal Dependence Model Scripts
7. MGPD_Dep.R - centralized likelihood function for MGPD dependence model. 
8. RevExp_U_Functions.R - Fitting MGPD dependence model to unreplicated data. 
9. CommonFunctions.R - commonfunctions for (8.)


