# RSFSA_R/Glmnet: R code for Random Subset Feature Selection Algorithm with Glmnet Niche Model from Tracy et al. (2018, submitted)

A Fire Vignette R program (RSFSA_GlmnetBinomialQuad_R_LoBrnFireData.R in RCode/FireVignette folder) and associated functions (functions beginning in "Glmnet" in RCode/Functions folder; other functions in RSFSA_R/MaxEnt/RCode/Functions folder) and csv data (in RSFSA_R/MaxEnt/FireVignetteData folder) are provided for running the first part of the program without needing the 90 North American environmental rasters at 1 km resolution. The first part of the program serves to demonstrate the random subset feature selection algorithm with Glmnet. The last part of the program invovles projecting Glmnet models and would require environmental rasters not provided. Some instructions and comments are provided at the beginning and throughout the vignette program. The Fire Vignette R program is set up to run provided Low Burn Severity csv files (11), but it can be modified to run the provided 11 csv files each for Moderate and High Burn Severities.

The program RSFSA_GlmnetBinomialQuad_R_LoBrnFireData_Compare.R is designed to compare with MaxEnt output by using training data and random subsets generated with a MaxEnt version of the program (can be modified to generate its own data).

Tracy JL, Trabucco A, Lawing AM, Giermakowski T, Tchakerian M, Drus GM, Coulson RN
 (2018) Random subset feature selection of ecological niche models for wildfire activity in western North America. Ecological Modelling (submitted)
