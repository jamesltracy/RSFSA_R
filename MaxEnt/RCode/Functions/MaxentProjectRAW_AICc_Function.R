MaxentProjectRAW.AICc <-
function(Species, VariableNamesIn, PredictorsIn, OutGridID, SubsetVariableNumber, VariableSubset, TrainSWD, TrainPresAbsID, PresenceDat.df, BetaMult, SetNameIn, OutDirectIn..., Output, CatVarsPrefix, MaxentArgsIn) {
  ##############################################################################
  # This R function outputs the Maxent model built upon all 
  # presence data with no evaluation statistics
  ##############################################################################
  #
  ##############################################################################
  # The model associated presence environmental data are processed using the R program
  # "SpeciesPresenceAbsencePointProcessingSTB18k.R" together with the ArcPython program
  # "SpeciesPresenceAbsencePointProcessingSTB18k.py"
  ##############################################################################
  #
  ##############################################################################
  # This section loads libraries, sets working directory and reads environmental rasters
  ##############################################################################
  #MaxentArgsIn=MaxentBaseArgsIn
  # Record start time of program
  starttime <- Sys.time()
  # load needed packages of raster, rgdal, dismo, rjava, and maptools (printouts not shown)
  # This sections reads in the lat/long data and formats it
  library(dismo)
  library(maptools)
  library(raster)
  library(rgdal)
  library(sp)
  library(ENMeval)
  #
  # Set default values for LongOutput, if optional values left out of function call
  MaxentBaseArgs <- c(paste0("betamultiplier=", BetaMult, "writebackgroundpredictions=true"))
  if(missing(MaxentArgsIn)) {
    MaxentArgs1 <- MaxentBaseArgs
  } else {
    MaxentArgs1 <- MaxentArgsIn
  }
  if(missing(CatVarsPrefix)) {
    AnyCategoricalVars=FALSE
    } else {
    AnyCategoricalVars=TRUE
  }
  #
  if(missing(OutGridID)) { OutGridID="" }
  # Create definition for a geographical projection
  crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
  #
  VariableNamesSel <- c(unlist(VariableSubset))
  # Split VarNames of VariableSubsets into separate variables
  VarNames <- unlist(strsplit(VariableNamesSel, "-"))
  SubsetVarNum <- length(VarNames)
  ## Keep only data for selected variables
  head(TrainSWD)
  MaxentTrainDataK <- as.data.frame(TrainSWD[,VarNames])
  colnames(MaxentTrainDataK) <- VarNames
  head(MaxentTrainDataK)
  nrow(MaxentTrainDataK)
  # Assemble environmental rasters
  names(PredictorsIn) <- toupper(names(PredictorsIn))
  names(PredictorsIn) <- gsub("_NS", "", names(PredictorsIn))
  #plot(PredictorsIn[[1]])
  # Subset raster stack of predictors by VarNames
  predictors <- subset(PredictorsIn, VarNames)
  ##
  #plot(predictors[[41]])
  # Develop maxent model from specified presence and background data
  #
  ## Account for categorical variables specified by CarVarsPrefix, if any
  if(AnyCategoricalVars==TRUE) {
    MaxentCatArg <- paste0("togglelayertype=", CatVarsPrefix)
    MaxentArgs <- c(MaxentArgs1, MaxentCatArg)
  } else {
    MaxentArgs <- MaxentArgs1
  }
  #
  MaxentOut <- maxent(MaxentTrainDataK, TrainPresAbsID, args=MaxentArgs, path=OutDirectIn)
  ######
  maxent.score <- predict(MaxentOut, predictors, args=c("outputformat=raw"))
  # Save for Arc as a geoTIFF grid
  setwd(OutDirectIn)
  writeRaster(maxent.score, paste0(Species, "Maxent", OutGridID, SubsetVariableNumber, "Vars_Beta", BetaMult), format = "GTiff", overwrite=TRUE)
  #
  ### Calculate AICc with ENMeval package
  # Obtain number of parameters in Maxent model
  nparam <- get.params(MaxentOut)
  # Derive occurence data
  MaxentPresData.df <- PresenceDat.df[,2:3]
  MaxentPresData.csv <- write.csv(MaxentPresData.df, file="MaxentPresData.csv", row.names=FALSE)
  occ <- read.table("MaxentPresData.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
  #
  MaxentAICcResults.df <- calc.aicc(nparam, occ, maxent.score)
 #####################################
  ## Calculate AICc_bg with point values
  # Obtain number of parameters in Maxent model
  nparam <- get.params(MaxentOut)
  #
  ## From ENMeval Package documentation: AICc is the Akaike Information Criterion corrected for small 
  ## sample sizes calculated as: (2 * K - 2 * logLikelihood) + (2 * K) * (K + 1)=(n - K - 1)
  ## where K is the number of parameters in the model (i.e., number of non-zero parameters in Maxent
  ## lambda file) and n is the number of occurrence localities.
  ## The logLikelihood is sum(log(vals/total))
  ## vals is vector of Maxent raw values at occurence localities
  ## total is the sum of Maxent raw values across the entire study area
  ##
  setwd(OutDirectIn)
  # Read in values for presence data output by maxent run
  presencevals.df <- data.frame(read.csv(paste0(MaxentOut@path, "/species_samplePredictions.csv")), header = TRUE, sep=',')
  # Keep only RAW values as vector
  vals <- presencevals.df[,4]
  head(vals)
  n <- length(vals) # number of occurence localities
  # total is sum of maxent raw values across entire study area, includes background and occurrence localities
  # Read in values for backround data output by maxent run
  backgroundvals.df <- data.frame(read.csv(paste0(MaxentOut@path, "/species_backgroundPredictions.csv")), header = TRUE, sep=',')
  head(backgroundvals.df)
  # Keep only RAW values as vector
  backgroundvals <- backgroundvals.df[,3]
  # Calculate sum of all values
  totalocc <- sum(vals) # sum from occurrence localities
  totalbg <- sum(backgroundvals)  # sum from background localities
  total <- totalocc + totalbg  # grand total sum
  #
  logLikelihood <- sum(log(vals/total))
  K <- nparam
  AICc_bg <- (2*K - 2*logLikelihood) + (2*K)*(K+1)/(n-K-1)
  MaxentAICcResults.df$AICc_bg <- AICc_bg
  MaxentAICcResults.df$NumDVars <- K
  #############################################################################
  # Join above output with other data
  MaxentID.df <- as.data.frame(matrix(c(unlist(VariableSubset), unlist(SubsetVariableNumber), unlist(SetNameIn)), ncol=3, nrow=1), stringsAsFactors=FALSE)
  #str(MaxentID.df)
  colnames(MaxentID.df) <- c("VariableSubset", "SubsetVariableNumber", "SetName")
  MaxentAICcStats.df <- as.data.frame(cbind(MaxentID.df, MaxentAICcResults.df), stringsAsFactors=FALSE)
  #str(MaxentAICcStats.df)
  return(MaxentAICcStats.df)
}



