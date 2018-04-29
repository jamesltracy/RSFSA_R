MaxentProject <-
function(Species, VariableNamesIn, PredictorsIn, OutGridID, SubsetVariableNumber, VariableSubset, TrainSWD, TrainPresAbsID, Threshold, BetaMult, OutDirectIn..., Output, CatVarsPrefix, MaxentArgsIn) {
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
  #BetaMult=1
  # Record start time of program
  starttime <- Sys.time()
  # load needed packages of raster, rgdal, dismo, rjava, and maptools (printouts not shown)
  # This sections reads in the lat/long data and formats it
  library(dismo)
  library(maptools)
  library(raster)
  library(rgdal)
  library(sp)
  #
  # Set default values for LongOutput, if optional values left out of function call
  MaxentBaseArgs <- c(paste0("betamultiplier=", BetaMult), "writebackgroundpredictions=false")
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
  MaxentOut <- maxent(MaxentTrainDataK, TrainPresAbsID, args=MaxentArgs, path=OutDirectSub)
  ######
  # Project maxent model on native range
  #maxent.scoreraw1 <- maxent.scoreraw
  maxent.score1 <- predict(MaxentOut, predictors)
  # Multiply raw grid by 1000 and convert to integer
  maxent.score <- calc(maxent.score1, function(x) as.integer(x * 1000) )
  # Project maxent model on introduced range (if any)
  # Save for Arc as a geoTIFF grid
  setwd(OutDirectIn)
  writeRaster(maxent.score, paste0(Species, "Maxent", OutGridID, SubsetVariableNumber, "Vars_Beta", BetaMult), format = "GTiff", overwrite=TRUE)
  #
  ################################################################################
  ### Calibrate model using using mean threshold of mean model value (e.g., envelope.score)
  ### at maximum TSS from evaluation runs of k-fold data
  ##  Multiply threshold by 1000
  ThresholdK <- Threshold * 1000
  maxent.scorecal <- calc(maxent.score, function(x) ifelse(x < ThresholdK, 0, 1) )
  #plot(maxent.scorecal)
  # SubsetVariableNumber <- 19
  #OutGridID <- "Full"
  writeRaster(maxent.scorecal, paste0(Species, "Maxent", OutGridID, SubsetVariableNumber, "Vars_Beta", BetaMult, "Cal"), format = "GTiff", overwrite=TRUE)
  #
  return(maxent.scorecal)
}



