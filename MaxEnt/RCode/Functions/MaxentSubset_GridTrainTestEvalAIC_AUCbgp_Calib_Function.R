MaxentSubset.GridTrainTestEvalAIC_AUCbgp_Calib <-
function(Species, VariableNamesIn, PredictorsIn, SubsetVariableNumber, TotPres, Run, SetName, VariableSubset, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrppin, kfoldgrpain, CRS.In, OutDirectIn..., CatVarsPrefix, MaxentArgsIn, Output) {
  #CatVarsPrefix <- "ROADS_CAT"
  #MaxentArgsIn <- MaxentBaseArgsIn
  #kfoldgrpp <- kfoldpres
  #PseudoabsenceDat.df <- PseudoabsenceDat.dfin
  #NOTE: Function parameters after the "...", such as TestDataType, have to be set with an equal sign in the function call, such as ScoreType="mean"
  setwd(OutDirectIn)
  library(dismo)
  library(raster)
  library(ENMeval)
  library(rgdal)
  #
  # Make sure VariableNames is a data frame
  VariableNamesIn <- data.frame(VariableNamesIn, stringsAsFactors=FALSE)
  #
  RowNames.df <- data.frame(rownames(VariableSubset), stringsAsFactors=FALSE)  # Save row names to use in output
  SubsetSize.mat <- matrix(c("Singlets", "Doublets", "Triplets", "Quartets", "Quintets", "Sextets", "Septets", "Octets", "Nonets",
  "Dectets", "Undectets", "Duodectets","Tredectets", "Quattuordectets", "Quindectets", "Sexdectets", "Septendectets", "Octodectets", "Novemdectets",
  "Vigetets", "Unvigetets", "Duovigetets", "Trevigetets", "Quattuorvigetets", "Quinvigetets", "Sexvigetets", "Septenvigetets", "Octovigetet",
  "Novemvigetets", "Trigetets", "Untrigetets", "Duotrigetets", "Tretrigetets", "Quottuortrigetets", "Quintrigetets",
  "Sextrigetets", "Septentrigetets", "Octotrigetets", "Novemtrigetets", "Quadragetets", "Unquadragetets", "Duoquadragetets", "Trequadragetets",
  "Quattuorquadragetets", "Quinquadragetets", "Sexquadragetets", "Octoquadragetets", "Octoquadragetets", "Novemquadragetets", "Quinquagetets",
  "Unquinquagetets", "Duoquinquagetets", "Trequinguagetets", "Quattuorquinquagetets", "Quinquinquagetets",
  "Sexquinquagetets", "Septenquinquagetets", "Octoquinquagetets", "Novemquinquagetets", "Sexagetets"), ncol=1, nrow=60, byrow=TRUE, dimnames=list(c
   (seq(1:60)), c("Subset")))
  SubsetSize.df <- as.data.frame(SubsetSize.mat, stringsAsFactors=FALSE)
  Subset <- SubsetSize.df[SubsetVariableNumber,]
  if(is.na(Subset)) {
    Subset <- ""
  }
  TotVars <- nrow(VariableNamesIn)
  #
  # Set default values for LongOutput, if optional values left out of function call
  if(missing(CatVarsPrefix)) { CatVarsPrefix="" }
  MaxentBaseArgs <- c("betamultiplier=1.0", "writebackgroundpredictions=true")
  if(missing(MaxentArgsIn)) {
    MaxentArgs1 <- MaxentBaseArgs
  } else {
    MaxentArgs1 <- MaxentArgsIn
  }
  if(CatVarsPrefix=="") {
    AnyCategoricalVars=FALSE
    } else {
    AnyCategoricalVars=TRUE
  }
  #Output=FALSE
  if(missing(Output)) { Output=TRUE }
  #
  crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
  if(missing(CRS.In)) {CRS.In=crs.geo}
  #
  ModelType <- paste("Maxent", Subset, sep="")
  Model <- paste("Maxent for", SubsetVariableNumber, "Feature Subset")
  tail(VariableSubset)
  #VariableSubsets[2998,]
  ################
  #i=1
  ## Use function to round up from .5 from http://stackoverflow.com/questions/12688717/round-up-from-5-in-r
  Round2 <- function(x, n) {
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5
    z = trunc(z)
    z = z/10^n
    z*posneg
  }
  #
  # Split VarNames of VariableSubsets into separate variables
  #str(VariableSubsets)
  VariableNamesSel <- c(unlist(VariableSubset[1,1]))
  # Check if dash used to separate variables
  if(grepl("-", VariableNamesSel)==TRUE) {
    VarNames <- unlist(strsplit(VariableNamesSel, "-"))
  } else {
    VarNames <- unlist(VariableNamesSel)
  }
  SubsetVarNum <- length(VarNames)
  # Specify model training data
  MaxentPresTrainData <- PresenceDat.df[kfoldgrppin != Run, ]
  head(MaxentPresTrainData)
  ncol(MaxentPresTrainData)
  nrow(MaxentPresTrainData)
  MaxentAbsTrainData <- BackgroundDat.df
  head(MaxentAbsTrainData)
  tail(MaxentAbsTrainData)
  ncol(MaxentAbsTrainData)
  nrow(MaxentAbsTrainData)
  if(SubsetVarNum == 1) {
    MaxentPresTrainData1 <- data.frame(MaxentPresTrainData[,VarNames])
    colnames(MaxentPresTrainData1) <- VarNames
    head(MaxentPresTrainData1)
    MaxentAbsTrainData1 <- data.frame(MaxentAbsTrainData[,VarNames])
    colnames(MaxentAbsTrainData1) <- VarNames
    head(MaxentAbsTrainData1)
    TrainSWD <- rbind(MaxentPresTrainData1, MaxentAbsTrainData1)
    head(TrainSWD)
  } else {
    TrainSWD <- rbind(MaxentPresTrainData[,4:ncol(MaxentPresTrainData)], MaxentAbsTrainData[,4:ncol(MaxentAbsTrainData)])
    head(TrainSWD)
    nrow(TrainSWD)
  }
  TrainPresID <- data.frame(rep(1,nrow(MaxentPresTrainData)))
  colnames(TrainPresID) <- "ID"
  nrow(TrainPresID)
  TrainAbsID <- data.frame(rep(0,nrow(MaxentAbsTrainData)))
  colnames(TrainAbsID) <- "ID"
  TrainPresAbsID <- rbind(TrainPresID, TrainAbsID)
  head(TrainPresAbsID)
  tail(TrainPresAbsID)
  nrow(TrainPresAbsID)
  #
  ## Keep only data for selected variables
  head(TrainSWD)
  MaxentTrainDataK <- as.data.frame(TrainSWD[,VarNames])
  colnames(MaxentTrainDataK) <- VarNames
  head(MaxentTrainDataK)
  nrow(MaxentTrainDataK)
  ##
  ############ Run Maxent using SWD
  # Create subdirectory for Maxent output
  #SubsetVarNum=8
  #system.file("java", package="dismo")
  ## Account for categorical variables specified by CarVarsPrefix, if any
  if(AnyCategoricalVars==TRUE) {
    MaxentCatArg <- paste0("togglelayertype=", CatVarsPrefix)
    MaxentArgs <- c(MaxentArgs1, MaxentCatArg)
  } else {
    MaxentArgs <- MaxentArgs1
  }
  MaxentOut <- maxent(MaxentTrainDataK, TrainPresAbsID, args=MaxentArgs, path=OutDirectIn)
  # Save Maxent model
  setwd(OutDirectIn)
  saveRDS(MaxentOut, "MaxentModel.rds")
  #####################################
  ## Calculate AICc_bg with point values
  # Obtain number of parameters in Maxent model
  params <- MaxentOut@lambdas[1:(length(MaxentOut@lambdas)-4)]
  paramslist <- list()
  count = 0
  for(m in 1:length(params)) {
    #m=1
    params_m <- as.data.frame(unlist(strsplit(params[m], ",")), stringsAsFactors=FALSE)
    colnames(params_m) <- paste0("Param",m)
    if(as.numeric(params_m[2,1])!=0) {
      count = count + 1
      paramslist[[count]] <- params_m
    }
  }
  params.df <- as.data.frame(do.call(cbind, paramslist), stringsAsFactors=FALSE)
  nparams <- ncol(params.df)
  #
  ## From ENMeval Package documentation: AICc is the Akaike Information Criterion corrected for small
  ## sample sizes calculated as: (2 * K - 2 * logLikelihood) + (2 * K) * (K + 1)=(n - K - 1)
  ## where K is the number of parameters in the model (i.e., number of non-zero parameters in Maxent
  ## lambda file) and n is the number of occurrence localities.
  ## The logLikelihood is sum(log(vals/total))
  ## vals is vector of Maxent raw values at occurence localities
  ## total is the sum of Maxent raw values across the entire study area
  ##
  # Read in values for presence data output by maxent run
  presencevals.df <- data.frame(read.csv(paste0(MaxentOut@path, "/species_samplePredictions.csv")), header = TRUE, sep=',')
  # Change any RAW values of zero to 0.001
  presencevals.df[,4][presencevals.df[,4]==0] <- 0.001
  nrow(presencevals.df)
  # Keep only RAW values as vector
  vals <- presencevals.df[,4]
  head(vals)
  n <- length(vals) # number of occurence localities
  # total is sum of maxent raw values across entire study area, includes background and occurrence localities
  # Read in values for backround data output by maxent run
  backgroundvals.df <- data.frame(read.csv(paste0(MaxentOut@path, "/species_backgroundPredictions.csv")), header = TRUE, sep=',')
  head(backgroundvals.df)
  # Change any RAW values of zero to 0.001
  backgroundvals.df[,3][backgroundvals.df[,3]==0] <- 0.001
  # Keep only RAW values as vector
  backgroundvals <- backgroundvals.df[,3]
  # Calculate sum of all values
  totalocc <- sum(vals) # sum from occurrence localities
  totalbg <- sum(backgroundvals)  # sum from background localities
  total <- totalocc + totalbg  # grand total sum
  #
  logLikelihood <- sum(log(vals/total))
  K <- nparams
  AICc_bg <- (2*K - 2*logLikelihood) + (2*K)*(K+1)/(n-K-1)
  NumDVars <- K
  ###
  #############################################################################
  #str(MaxentOut)
  # Copy Maxent output to desired directory
  #flist <- list.files(MaxentOut@path, full.names = TRUE)
  #file.copy(flist, OutDirectSub, overwrite=TRUE)
  ##
  #####################################
  ## Run projection of model
  # Assemble environmental rasters
  names(PredictorsIn) <- toupper(names(PredictorsIn))
  names(PredictorsIn) <- gsub("_NS", "", names(PredictorsIn))
  #plot(PredictorsIn[[1]])
  system.time(maxent.score <- predict(MaxentOut, PredictorsIn))
  #dev.new()
  #plot(maxent.score)
  # Multiply continuous model grid by 1000 and convert to integer
  maxent.score1 <- calc(maxent.score, function(x) as.integer(x * 1000) )
  # Save grid
  setwd(OutDirectIn)
  writeRaster(maxent.score1, paste0(Species, "Maxent", OutGridID, SubsetVariableNumber, "Vars_Beta", BetaMult), format = "GTiff", overwrite=TRUE)
  #
  #################################################
  ### Calculate AICc with ENMeval package
  # Obtain number of parameters in Maxent model
  setwd(OutDirectIn)
  nparam <- get.params(MaxentOut)
  # Derive training occurence data
  PresTrainData.df1 <- PresenceDat.df[kfoldgrppin != Run, ]
  PresTrainData.df <- PresTrainData.df1[,2:3]
  head(PresTrainData.df)
  nrow(PresTrainData.df)
  PresTrainData.csv <- write.csv(PresTrainData.df, file="PresTrainData.csv", row.names=FALSE)
  occ <- read.table("PresTrainData.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
  #
  MaxentAICcResults.df <- calc.aicc(nparam, occ, maxent.score)
  AICc <- MaxentAICcResults.df$AICc
  ####################################################
  # Obtain number of parameters in Maxent model
  params <- MaxentOut@lambdas[1:(length(MaxentOut@lambdas)-4)]
  paramslist <- list()
  count = 0
  for(m in 1:length(params)) {
    #m=7
    params_m <- as.data.frame(unlist(strsplit(params[m], ",")), stringsAsFactors=FALSE)
    colnames(params_m) <- paste0("Param",m)
    if(as.numeric(params_m[2,1])!=0) {
      count = count + 1
      paramslist[[count]] <- params_m
    }
  }
  params.df <- as.data.frame(do.call(cbind, paramslist), stringsAsFactors=FALSE)
  nparams <- ncol(params.df)
  ##
  ## Obtain number of environmental variables used in parameters
  envarslist <- list()
  count = 0
  for(n in 1:ncol(MaxentTrainDataK)) {
    if(length(grep(colnames(MaxentTrainDataK)[n], params.df[1,]))>0) {
      count = count + 1
      envarslist[[count]] <- colnames(MaxentTrainDataK)[n]
    }
  }
  envars.df <- as.data.frame(do.call(rbind, envarslist), stringsAsFactors=FALSE)
  nenvars <- nrow(envars.df)
  if(nenvars > 1) {
    EnvVarsUsed <- paste0(unlist(envarslist), collapse="-")
  } else if(nenvars ==1) {
    EnvVarsUsed <- unlist(envarslist)
  } else {
    EnvVarsUsed <- " "
  }
  ##
  ############ Calculation check varying K relationship with AIC  ###################################
#  AICcalcs <- list()
#  CountK <- 0
#  for(q in 29:100) {
#    #q=29
#    K <- q
#    nparam <- K
#    AICc_Kres.df <- calc.aicc(nparam, occ, maxent.score)
#    AICc_K <- as.vector(AICc_Kres.df$AICc)
#    AICc_bg_K <- (2*K - 2*logLikelihood) + (2*K)*(K+1)/(n-K-1)
#    AICcalcs.df1 <- data.frame(t(as.matrix(c(nparam, AICc_K, AICc_bg_K))))
#    colnames(AICcalcs.df1) <- c("nparam", "AICc_K", "AICc_bg_K")
#    CountK <- CountK + 1
#    AICcalcs[[CountK]] <- AICcalcs.df1
#  }
#  AICcalcs.df <- do.call(rbind, AICcalcs)
#  write.table(AICcalcs.df, sep = ",", col.names=NA, file=paste0(Species, "_AIC_Vs_K_calc_check.csv", sep=""))
  #######################################################
  ##
  ########### Evaluate model using training and test data
  MaxentSubsetEvalL <- list()
  TestDataTypes <- c("Train", "Test", "Test_bgp")
  for(k in 1:3) {
    #k=1
    TestDataType <- TestDataTypes[k]
    if(k==1) {
    # Specify model testing data
      MaxentPresTestData <- PresenceDat.df[kfoldgrppin != Run, ]
      head(MaxentPresTestData)
      nrow(MaxentPresTestData)
      MaxentAbsTestData <- PseudoabsenceDat.df[kfoldgrpain != Run, ]
      head(MaxentAbsTestData)
      nrow(MaxentAbsTestData)
      tail(MaxentAbsTestData)
    } else if (k==2) {
      MaxentPresTestData <- PresenceDat.df[kfoldgrppin == Run, ]
      head(MaxentPresTestData)
      nrow(MaxentPresTestData)
      tail(MaxentPresTestData)
      MaxentAbsTestData <- PseudoabsenceDat.df[kfoldgrpain == Run, ]
      head(MaxentAbsTestData)
      tail(MaxentAbsTestData)
      nrow(MaxentAbsTestData)
    } else {
      MaxentPresTestData <- PresenceDat.df[kfoldgrppin == Run, ]
      head(MaxentPresTestData)
      nrow(MaxentPresTestData)
      MaxentAbsTestData <- rbind(PresenceDat.df, BackgroundDat.df)
      head(MaxentAbsTestData)
      nrow(MaxentAbsTestData)
      tail(MaxentAbsTestData)
    }
    ########################
    ### Query projection grid to get values for test presence and absence points
    ## First presence points
    setwd(OutDirectIn)
    head(MaxentPresTestData)
    # Convert point data.frame to SpatialPointsDataFrame
    # First, specify xy coordinates
    xy <- MaxentPresTestData[,c("coords.x1", "coords.x2")]
    # Create spatial points data frame
    PresTestPnt.spdf <- SpatialPointsDataFrame(coords=xy, data=MaxentPresTestData, proj4string=CRS.In)
    ## Extract values of presence test points from maxent model grid
    system.time(prespred.df <- data.frame(extract(maxent.score, PresTestPnt.spdf)))
    colnames(prespred.df) <- "MaxentScore"
    # Rejoin coordinates to prespred.df and save as shapefile for checking
    prespred.spdf <- SpatialPointsDataFrame(coords=xy, data=prespred.df, proj4string=CRS.In)
    # Write shapefile including the correct projection
    writeOGR(prespred.spdf, ".", paste0(Species, "PresenceTestMaxentVals"), driver="ESRI Shapefile", overwrite=TRUE)
    ## Then absence points
    head(MaxentAbsTestData)
    # Convert point data.frame to SpatialPointsDataFrame
    # First, specify xy coordinates
    xy <- MaxentAbsTestData[,c("coords.x1", "coords.x2")]
    # Create spatial points data frame
    AbsTestPnt.spdf <- SpatialPointsDataFrame(coords=xy, data=MaxentAbsTestData, proj4string=CRS.In)
    ## Extract values of presence test points from maxent model grid
    system.time(abspred.df <- data.frame(extract(maxent.score, AbsTestPnt.spdf)))
    colnames(abspred.df) <- "MaxentScore"
    # Rejoin coordinates to abspred.df and save as shapefile for checking
    abspred.spdf <- SpatialPointsDataFrame(coords=xy, data=abspred.df, proj4string=CRS.In)
    # Write shapefile including the correct projection
    writeOGR(abspred.spdf, ".", paste0(Species, "AbsenceTestMaxentVals"), driver="ESRI Shapefile", overwrite=TRUE)
    #############################################################################################
    # This section evaluates the Maxent model using the PresenceAbsence package
    #############################################################################################
    #### Create a dataset with model predictions for presence and absence points
    # Use extracted prediction values for presence and absence points for each of three
    # models previously calculated in loop
    #
    library(gtools)
    ## Create directory of output for class pair run
    #setwd(paste("C:/Users/JLTracy/Documents/R/win-library/3.0/10minClimIntEvalSTB1000", "/", output2, sep=""))
    # For presence data, assign a column of "1" under the name "OBSERVED" to indicate presence data
    # and assign the model results a name "MaxentScoreN" where N is the name of the rep
    # Also assign a column "id" for the row numbers to use in merging the data frames later
    presa.df <- data.frame(c(rep(1, nrow(prespred.df))))
    names(prespred.df) <- c("Maxent")
    names(presa.df) <- c("OBSERVED")
    pres.df <- data.frame(cbind(id=1:nrow(presa.df), presa.df, prespred.df))
    nrow(pres.df)
    # Repeat above process with absence data, but assign "OBSERVED" a value of 0
    absa.df <- data.frame(c(rep(0, nrow(abspred.df))))
    names(abspred.df) <- c("Maxent")
    names(absa.df) <- c("OBSERVED")
    abs.df <- data.frame(cbind(id=1:nrow(absa.df), absa.df, abspred.df))
    # For each model output, merge presence and absence data using "id' column as guide when all=TRUE
    # NOTE: PresenceAbsence package cannot handle several models at one time if the sample sizes differ
    # so have to analyze each model output separately
    presabspred <- rbind(pres.df, abs.df)
    tail(presabspred)
    head(presabspred)
    # Drop the id column used in merging for each dataset
    presabspred$id <- NULL
    # Make a column of data with the species name with same number of rows as data from each model
    SPECIES <- data.frame(c(rep(Species, nrow(presabspred))))
    names(SPECIES) <- c("SPECIES")
    # Make final dataset SPDATA by putting together SPECIES with extracted environmental data.
    SPDATA <- data.frame(SPECIES, presabspred)
    head(SPDATA)
    #SPDATA[100:160,]
    ################################################################################
    ### Run this block of code to evaluate model results with PresenceAbsence package
    ################################################################################
    library(PresenceAbsence)
    #setwd(paste("C:/Users/JLTracy/Documents/R/win-library/3.0/10minClimIntEvalSTB1000", "/", output2, sep=""))
    #starttime <- Sys.time()
    #### FOR OLD WORLD DATA EVALUATION STATISTICS
    ### Define variables for later use.
    accurun <- list()
    accusum <- matrix(data=NA, ncol=9, nrow=1, byrow=TRUE, dimnames=list(NULL, c("MaxTSS", "Specificity_maxTSS", "Sensitivity_maxTSS", "AUC", "MaxKappa", "ThresholdMaxTSS", "AICc_bg", "NumDVars", "AICc")))
    accusum <- matrix(data=NA, ncol=11, nrow=1, byrow=TRUE, dimnames=list(NULL, c("MaxTSS", "Specificity_maxTSS", "Sensitivity_maxTSS", "AUC", "MaxKappa", "ThresholdMaxTSS", "AICc_bg", "NumDVars", "NumEnVars", "EnvVarsUsed", "AICc")))
    species <- as.character(unique(SPDATA$SPECIES))
    model.names <- as.character(names(SPDATA)[-c(1, 2)])
    N.models <- ncol(SPDATA) - 2
    N.sp <- length(species)
    N.obs <- length(SPDATA$SPECIES[SPDATA$SPECIES == species[1]])
    Obs.prev <- table(SPDATA$SPECIES, SPDATA$OBSERVED)[, 2]/N.obs
    Obs.prev <- Round2(Obs.prev, 2)
    ### Mainly just run this code
    graphics.off()
    sp <- 1
    # Read in dataset for loop
    DATA <- SPDATA[SPDATA$SPECIES == species[sp], ]
    head(DATA)
    #
    # To assess accuracy per threshold, use limited threshold available for
    #  model based upon number of environmental layers in model
    # ("NumGrids")
    #NumGrids <- max(40, SubsetVarNum)
    #PossThresholds <- seq(1/NumGrids,1,length=NumGrids)
    PossThresholds <- 100
    #accu <- data.frame(presence.absence.accuracy(SPDATA, which.model = 1, threshold = PossThresholds, st.dev=FALSE))
    # accu <- presence.absence.accuracy(DATA, which.model = 1, threshold = c(0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.975, 0.98, 0.99, 0.999999))
    #accu <- presence.absence.accuracy(DATA, which.model = 1, threshold = 100, st.dev=FALSE)
    # print(paste("Species:", species[sp], "Model:", model.names))  not used
    accu <- data.frame(presence.absence.accuracy(DATA, which.model = 1, threshold = 100, st.dev=FALSE))
    # print(paste("Species:", species[sp], "Model:", model.names))  not used
    head(accu)
    maxSSS <-  data.frame(accu$sensitivity + accu$specificity)
    names(maxSSS) <- c("maxSSS")
    head(maxSSS)
    TSS <-  data.frame(accu$sensitivity + accu$specificity - 1)
    names(TSS) <- c("TSS")
    accurun <- data.frame(accu, maxSSS, TSS)
    head(accurun)
    accurun$Conditions <- paste("Maxent", Subset, sep="")
    maxKappa <- max(accurun$Kappa)
    maxTSS <- max(accurun$TSS)
    AUC <- max(accurun$AUC)
    # Find and average thresholds at TSS = maxTSS. In the case of tied optimal
    # thresholds, we select the mean threshold producing maximum TSS following
    # (Freeman and Moisen 2008). But, in the case of discrete thresholds as found
    # in envelope models, if the mean optimal threshold does not represent an
    # actual discrete threshold, we select the nearest discrete threshold to the
    # mean among 3 or more thresholds, or the smaller of two adjacent discrete thresholds.
    ThresholdsMaxTSS <- accurun$threshold[which(accurun$TSS == maxTSS)]
    ThresholdMaxTSS <- mean(ThresholdsMaxTSS)
    #ThresholdMaxTSSM <- mean(ThresholdsMaxTSS)
    ## Following commented code for envelope score
    #if (length(ThresholdsMaxTSS) < 3) {
    #ThresholdMaxTSS <- min(ThresholdsMaxTSS)
    #} else {  ThresholdMaxTSS <- PossThresholds[which(abs(PossThresholds - ThresholdMaxTSSM)== min(abs(PossThresholds - ThresholdMaxTSSM)))]
    #}
    # Calculate specificity and sensitivity at maxTSS
    Specificity_maxTSS <- accurun$specificity[which(accurun$TSS == maxTSS)]
    Specificity_maxTSS <- mean(Specificity_maxTSS)
    Sensitivity_maxTSS <- accurun$sensitivity[which(accurun$TSS == maxTSS)]
    Sensitivity_maxTSS <- mean(Sensitivity_maxTSS)
    #CheckTSS <- Specificity_maxTSS + Sensitivity_maxTSS - 1 # should equal maxTSS
    accusum[1,1] <- max(maxTSS, 0.0001)
    accusum[1,2] <- max(Specificity_maxTSS, 0.0001)
    accusum[1,3] <- max(Sensitivity_maxTSS, 0.0001)
    accusum[1,4] <- max(AUC, 0.0001)
    accusum[1,5] <- max(maxKappa, 0.0001)
    accusum[1,6] <- max(mean(ThresholdMaxTSS), 0.0001)
    accusum[1,7] <- AICc_bg
    accusum[1,8] <- NumDVars
    accusum[1,9] <- nenvars
    accusum[1,10] <- EnvVarsUsed
    accusum[1,11] <- AICc
    # Save Threshold value and multiply by 1000 for grid calibration
    ThresholdK <-  (max(ThresholdMaxTSS, 0.0001))*1000
    ###############################
    #endtime <- Sys.time()
    #durtime <- endtime - starttime
    # Save evaluation statistics to .csv file
    #setwd(paste("C:/Users/JLTracy/Documents/R/win-library/3.0/10minClimIntEvalSTB1000", "/", output2, sep=""))
    accusum.df <- data.frame(accusum, stringsAsFactors=FALSE)
    accusum.df$Model <- ModelType
    #
    # Save variable matrix coordinates and Old World evaluation statistics to a vector
    MaxentEvalStats.df <- data.frame(t(as.matrix(c(VariableNamesSel, accusum.df$EnvVarsUsed, SubsetVariableNumber, TestDataType, TotPres, SetName, Run, accusum.df$MaxTSS, accusum.df$MaxKappa, accusum.df$AUC, accusum.df$Specificity_maxTSS, accusum.df$Sensitivity_maxTSS, accusum.df$ThresholdMaxTSS, accusum.df$AICc_bg, accusum.df$NumDVars, accusum.df$NumEnVars, accusum.df$AICc))))
    #
    setwd(OutDirectIn)
    # Save threshold verus evaluation statitistic data for run
    write.table(SPDATA, file=paste(Species, TestDataType, "SPDATA_Run", Run, ".csv", sep=""), sep=",")
    write.table(accurun, file=paste(Species, TestDataType, "StatsVsThresholds_Run", Run, ".csv", sep=""), sep=",", col.names=NA)
    #print(paste("Species:", species[sp], "Model:", model.names))
    #print(accurun)
    bmp(paste(Species, TestDataType,"TSSvsThreshold_Run", Run, ".bmp", sep=""))
    #dev.new()
    plot(accurun$threshold, accurun$TSS, type="l")
    dev.off()
    # Plot a ROC plot
    # save following plot of grid as .bmp file
    bmp(paste(Species, TestDataType, "MaxentROC_Run", Run, ".bmp", sep=""))
    #dev.new()
    auc.roc.plot(DATA, color = TRUE, legend.cex = 1.2, main = "")
    dev.off()
    # Calculate optimal thresholds by various methods and save output to text file
    #outthresholds <- capture.output(optimal.thresholds(DATA, opt.methods = 1:12, threshold=10001))
    #out1 <- paste(outthresholds)
    #cat(out1, file=paste(Species, "MaxentResults", TestDataType, "Thresholds.txt", sep=""), sep="\n", append=TRUE)
    # Can calculate the confusion matrix for a given threshold (not necessary)
    confmatrix <- cmx(DATA, threshold = ThresholdMaxTSS, na.rm = FALSE)
    # Save confusion matrix at threshold of maxTSS
    write.table(confmatrix, file=paste(Species, "MaxentResults", TestDataType, "ConfusionMatrixatMaxTSSThreshold_Run", Run, ".csv", sep=""), sep=",", col.names=c("Obs Pres","Obs Abs"), row.names=c("Pred Pres", "Pred Abs"))
    #Save evaluation stastics to file
    write.table(accusum.df, file=paste(Species, "MaxentResults", TestDataType, "Stats_Run", Run, ".csv", sep=""), sep=",", col.names=NA)
    graphics.off()
    # Save output evaluation statistics and summary statistics to text file
    # Output codes of variables used
    outA <- paste("\n\nSpecies:", Species, "\nVariable Subset:", VariableNamesSel, "\nModel:", Model, "\nScoring Alogorithm: Maxent")
    cat(outA, file=paste(Species, "MaxentResults", TestDataType, "Summary_Run", Run, ".txt", sep=""), sep="\n", append=TRUE)
    #NotUsedVars <- capture.output(print(DeleteBands))
    #outB <- paste("\nVariables Not Used: Code-", NotUsedVars, "; Variable-", GridNamesDrop)
    #cat(outB, file=paste(Species, InfClass, "MaxentResults", output2, "Summary.txt", sep=""), sep="\n", append=TRUE)
    #UsedVars <- capture.output(print(SRandNumbers))
    #outB <- paste("\n  Code Nos. for Variables Used:", UsedVars)
    #cat(outB, file=paste("MaxentResults", output2, "Summary.txt", sep=""), sep="\n", append=TRUE)
    #out1 <- paste("\nVariable Names")
    #cat(out1, file=paste(Species, "MaxentResults", GridNamesKeepSelT, "Summary.txt", sep=""), sep="\n", append=TRUE)
    #EnvVariables <- as.matrix(GridNamesKeepSel)
    #out2 <- capture.output(print(EnvVariables))
    #cat(out2, file=paste(Species, "MaxentResults", GridNamesKeepSelT, "Summary.txt", sep=""), sep="\n", append=TRUE)
    out3 <- paste("\n  Evaluation Statistics per Run")
    cat(out3, file=paste(Species, "MaxentResults", TestDataType, "Summary_Run", Run, ".txt", sep=""), sep="\n", append=TRUE)
    out4 <- capture.output(print(accusum.df))
    cat(out4, file=paste(Species, "MaxentResults", TestDataType, "Summary_Run", Run, ".txt", sep=""), sep="\n", append=TRUE)
    #
    #############
    MaxentSubsetEvalStats.df1 <- MaxentEvalStats.df
    head(MaxentSubsetEvalStats.df1)
    tail(MaxentSubsetEvalStats.df1)
    nrow(MaxentSubsetEvalStats.df1)
    #str(ESWrapperOutput.df)
    colnames(MaxentSubsetEvalStats.df1) <- c("VarNames", "EnvVarsUsed", "SubsetVariableNumber", "DataType", "TotalPresPnts", "SetName", "Run", "TSS", "Kappa", "AUC", "Spec", "Sens", "ThreshMxTSS", "AICc_bg", "NumDVars", "NumEnVars", "AICc")
    rownames(MaxentSubsetEvalStats.df1) <- c(seq(1:nrow(MaxentEvalStats.df)))
    MaxentSubsetEvalStats.df1 <- as.data.frame(MaxentSubsetEvalStats.df1, stingsAsFactors=FALSE)
    head(MaxentSubsetEvalStats.df1)
    ncol(MaxentSubsetEvalStats.df1)
    # Save output matrix
    # Convert sixth through 15th columns from character to numeric
    MaxentSubsetEvalStats.df1[,7:17] <- sapply(MaxentSubsetEvalStats.df1[,7:17], function(x) as.numeric(as.character(x)))
    if(k < 3) {
      MaxentSubsetEvalL[[k]] <- MaxentSubsetEvalStats.df1
    } else {
      MaxentSubsetEvalL[[4]] <- MaxentSubsetEvalStats.df1
    }
    #Calculate difference between training and test statistics (overfitting)
    if(k==2) {
      MaxentSubsetEvalL[[3]] <- MaxentSubsetEvalL[[2]]
      MaxentSubsetEvalL[[3]][,8:17] <- MaxentSubsetEvalL[[1]][,8:17] - MaxentSubsetEvalL[[2]][,8:17]
      MaxentSubsetEvalL[[3]][,4] <- "Diff"
    }
  }
  #####################################################################################
  MaxentSubsetEvalStats.df <- do.call(rbind, MaxentSubsetEvalL)
  #
  ################################################################################
  ### Calibrate model using using above calculated ThresholdK
  ### at maximum TSS from above evaluation 
  maxent.scorecal <- calc(maxent.score1, function(x) ifelse(x < ThresholdK, 0, 1) )
  #plot(maxent.scorecal)
  # SubsetVariableNumber <- 19
  #OutGridID <- "Full"
  writeRaster(maxent.scorecal, paste0(Species, "Maxent", OutGridID, SubsetVariableNumber, "Vars_Beta", BetaMult, "Cal"), format = "GTiff", overwrite=TRUE)
  ##############
  #
  if(Output==TRUE) {
    setwd(OutDirectIn)
    Sets <- nrow(MaxentSubsetEvalStats.df)
    write.table(MaxentSubsetEvalStats.df, file=paste(Species, "Maxent_TrainTest_", TotVars, "_TotVars_", SubsetVariableNumber, "Var_", Sets, Subset, "_", SetName, "_", Run, ".csv", sep=""), sep=",", col.names=NA)
  }
  #
  gc()
  return(MaxentSubsetEvalStats.df)
}
