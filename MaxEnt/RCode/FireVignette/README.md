# RSFSA_R/MaxEnt: R code for Random Subset Feature Selection Algorithm with MaxEnt Niche Model from Tracy et al. (2018, submitted)

Order of running programs:

1) RSFSA_MaxEnt_R_LoBrnFireData.R - models low burn severity wildfire activity data, ranking models by AUC
2) RSFSA_MaxEnt_R_LoBrnFireDataAICc.R - uses output of above program to re-rank model by AICc instead of AUC
3) RSFSA_MaxEnt_R_ModBrnFireData.R -models moderate burn severity wildfire activity data
4) RSFSA_MaxEnt_R_ModBrnFireDataAICc.R - uses output of above program to re-rank model by AICc instead of AUC
5) RSFSA_MaxEnt_R_HiBrnFireData.R -models high burn severity wildfire activity data
6) RSFSA_MaxEnt_R_HiBrnFireDataAICc.R - uses output of above program to re-rank model by AICc instead of AUC

