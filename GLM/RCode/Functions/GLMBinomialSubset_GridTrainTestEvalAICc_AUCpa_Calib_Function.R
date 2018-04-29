GLMBinomialSubset.GridTrainTestEvalAICc_AUCpa_Calib <-
function(Species, VariableNamesIn, PredictorsIn, SubsetVariableNumber, TotPres, Run, SetName, VariableSubset, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrppin, kfoldgrpain, CRS.In, OutDirectIn..., Output) {
  #kfoldgrppain <- NObskfoldgrppa
  #NOTE: Function parameters after the "...", such as TestDataType, have to be set with an equal sign in the function call, such as ScoreType="mean"
  #FunctDirectIn <- "C:/Users/James/Documents/R/win-library/"
  #
  setwd(OutDirectIn)
  library(AICcmodavg)
  library(dismo)
  library(raster)
  #
  crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
  if(missing(CRS.In)) {CRS.In=crs.geo}
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
  #
  ModelType1 <- paste("GLMBinomial", Subset, sep="")
  Model <- paste("GLM Binomial for", SubsetVariableNumber, "Feature Subset")
  tail(VariableSubset)
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
  ##
  ## Designate final training data and add count column as factor
  GLMPresTrainData <- PresenceDat.df[kfoldgrppin != Run,]
  GLMPresTrainData$Count <- as.factor(1)
  head(GLMPresTrainData)
  #str(GLMPresTrainData)
  GLMAbsTrainData <- PseudoabsenceDat.df[kfoldgrpain != Run,]
  GLMAbsTrainData$Count <- as.factor(0)
  ################
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
  #####################################
  ##
  VariableNamesSel <- c(unlist(VariableSubsetsIn[i,]))
  # Split VarNames of VariableSubsetsIn into separate variables
  if(grepl("-",VariableNamesSel)) {
    VarNames <- unlist(strsplit(VariableNamesSel, "-"))
  } else {
    VarNames <- unlist(VariableNamesSel)
  }
  SubsetVarNum <- length(VarNames)
  ###
  ######################
  # Join presence and absence data
  PresAbsTrainDat.df <- rbind(GLMPresTrainData, GLMAbsTrainData)
  nrow(PresAbsTrainDat.df)
  #
  x <- as.matrix(PresAbsTrainDat.df[,c(1:ncol(PresAbsTrainDat.df)-1)])
  nrow(x)
  head(x)
  # Subset by environmental variables to be used in model from VariableSubset
  x <- x[, c(VarNames), drop=FALSE]
  x <- apply(x, 2, as.numeric)
  head(x)
  #str(x1)
  ################
  # Specify model training data
  xTrainData.df <- data.frame(x)
  head(xTrainData.df)
  nrow(xTrainData.df)
  #####################
  # Identify response variable for training data
  yTrainData <- PresAbsTrainDat.df$Count
  length(yTrainData)
  #
  # Prepare formula for GLM
  fla <- paste("yTrainData ~", paste(VarNames, collapse=" + "))
  # Run GLM
  GLMBinomialModel <- glm(as.formula(fla), data=xTrainData.df, family="binomial")
  # Save GLM model
  setwd(OutDirectIn)
  saveRDS(GLMBinomialModel, "GLMBinomialModel.rds")
  #
  # Calculate AICc and extract coefficients
  AIC <- AIC(GLMBinomialModel)
  AICc <- AICc(GLMBinomialModel)
  #str(GLMBinomialModel)
  # Extract data frame of non-zero coefficients from model
  Coefficients.df1 <- data.frame(as.matrix(GLMBinomialModel$coefficients), stringsAsFactors=FALSE)
  colnames(Coefficients.df1) <- "Coefficients"
  Coefficients.df <- Coefficients.df1[ which(Coefficients.df1$Coefficients!=0), , drop=FALSE]
  # Caculate AICc using formula from Johnny Heineken at  https://stats.stackexchange.com/questions/25817/is-it-possible-to-calculate-aic-and-bic-for-lasso-regression-models
  deviance.fit <- GLMBinomialModel$deviance
  fit.nulldev <- GLMBinomialModel$null.deviance
  k <- nrow(Coefficients.df) - 1
  n <- length(yTrainData)
  tLL <- fit.nulldev - deviance.fit
  AICc2 <- -tLL+2*k+2*k*(k+1)/(n-k-1)
  #AICcor <- (2*(k^2 + k))/(n-k-1)
  #AICc4 <- AIC + AICcor
  ####################################################
  #
  VarNamesUsed <- as.vector(rownames(Coefficients.df)[2:nrow(Coefficients.df)])
  VarsUsed <- gsub("\\^2", "", VarNamesUsed)
  VarsUsed <- unique(VarsUsed)
  nenvars <- length(VarsUsed)
  if(nenvars > 1) {
    EnvVarsUsed <- paste0(VarsUsed, collapse="-")
  } else if(nenvars ==1) {
    EnvVarsUsed <- VarsUsed
  } else {
    EnvVarsUsed <- " "
  }
  ##
  #####################################
  ## Calculate AICc_bg with point values
  # Obtain number of parameters in GLM model
  nparams <- (nrow(Coefficients.df) - 1)
  #
  ## From ENMeval Package documentation: AICc is the Akaike Information Criterion corrected for small
  ## sample sizes calculated as: -1*((2 * K - 2 * logLikelihood) + (2 * K) * (K + 1)=(n - K - 1))
  ## (-1 factor is correction) where K is the number of parameters in the model (i.e., number of non-zero parameters in GLM
  ## lambda file) and n is the number of occurrence localities.
  ## The logLikelihood is sum(log(vals/total))
  ## vals is vector of GLM raw presence values at occurence localities
  ## total is the sum of GLM raw values across the entire study area
  ##
  head(PresAbsTrainDat.df)
  nrow(PresAbsTrainDat.df)
  ######
  xPresTrainData.df1 <- PresAbsTrainDat.df[1:nrow(GLMPresTrainData),]
  nrow(xPresTrainData.df1)
  # Select only terms used in ensemble model
  xPresData.dfIn <- xPresTrainData.df1[, c(VarNamesUsed), drop=FALSE]
  head(xPresData.dfIn)
  # Obtain values for training presence data output by GLM run
  presvals.df <- setNames(data.frame(predict(GLMBinomialModel, xPresData.dfIn, type="response")), "GLMScore")
  head(presvals.df)
  #presvals.df <- GLMBinomialResp_PredMatrix(CoefficientsIn, xPresData.dfIn)
  #str(xTrainData)
  # Obtain values for background data output by GLM run
  # Select only terms used in ensemble model
  xData.dfIn <- BackgroundDat.df[, c(VarNamesUsed), drop=FALSE]
  head(xData.dfIn)
  backgroundvals.df <- setNames(data.frame(predict(GLMBinomialModel, xData.dfIn, type="response")), "GLMScore")
  #backgroundvals.df <- cvGLMEnsemble_PredMatrix_BinomialResp(CoefficientsIn, xData.dfIn)
  head(backgroundvals.df)
  ## Join together presence and background GLMscores
  PresBckgrndData.df <- rbind(presvals.df, backgroundvals.df)
  # Change any values of zero to 0.001
  PresBckgrndData.df[,1][PresBckgrndData.df[,1]==0] <- 0.001
  ## Sum all values
  SumVal <- apply(PresBckgrndData.df,2,sum)
  ## Divide all values by SumVal
  SumValDivFunc <- function(x) {x/SumVal}
  PresBckgrndDataRAW.df <- data.frame(apply(PresBckgrndData.df,2,SumValDivFunc))
  head(PresBckgrndDataRAW.df)
  nrow(PresBckgrndDataRAW.df)
  # Keep values for presence data as vector
  vals <- PresBckgrndDataRAW.df[1:nrow(GLMPresTrainData),1]
  head(vals)
  n <- length(vals) # number of occurence localities
  # Keep values for background  data as vector
  backgroundvals <- PresBckgrndDataRAW.df[(nrow(GLMPresTrainData)+1):nrow(PresBckgrndDataRAW.df),1]
  length(backgroundvals)
  # total is sum of GLM raw values across entire study area, includes background and occurrence localities
  # Calculate sum of all values
  totalocc <- sum(vals) # sum from occurrence localities
  totalbg <- sum(backgroundvals)  # sum from background localities
  total <- totalocc + totalbg  # grand total sum
  #
  logLikelihood <- sum(log(vals/total))
  K <- nparams
  AICc_bg <- -1*((2*K - 2*logLikelihood) + (2*K)*(K+1)/(n-K-1))
  NumDVars <- K
  ###
  #####################################
  ## Run projection of model
  # Assemble environmental rasters
  PredictorsIn <- Predictors[[VarNames]]
  # Project model onto raster stack of environmental variables
  ## Took 1.6 hrs for four of 10 variables over funnel
  system.time(GLMBinomial.score <- predict(PredictorsIn, GLMBinomialModel))
  # Multiply continuous model grid by 1000 and convert to integer
  system.time(GLMBinomial.score1 <- calc(GLMBinomial.score, function(x) as.integer(x * 1000) ))
  # Save grid
  setwd(OutDirectIn)
  writeRaster(GLMBinomial.score1, paste0(Species, "GLMBinomial", OutGridID, SubsetVariableNumber), format = "GTiff", overwrite=TRUE)
  # GLMBinomial.score1 <- raster(paste0(Species, "GLMBinomial", OutGridID, SubsetVariableNumber, ".tif"))
  # str(GLMBinomial.score)
  # GLMBinomial.score <- calc(GLMBinomial.score1, function(x) x/1000)
  ##
  ########### Evaluate model using training and test data
  GLMSubsetEvalL <- list()
  TestDataTypes <- c("Train", "Test", "Test_bgp")
  for(k in 1:3) {
    #k=3
    TestDataType <- TestDataTypes[k]
    if(k==1) {
      # Specify model testing data
      PresTestData <- PresenceDat.df[kfoldgrppin != Run, ]
      AbsTestData <- PseudoabsenceDat.df[kfoldgrpain != Run, ]
    } else if (k==2) {
      PresTestData <- PresenceDat.df[kfoldgrppin == Run, ]
      AbsTestData <- PseudoabsenceDat.df[kfoldgrpain == Run, ]
    } else {
      PresTestData <- PresenceDat.df[kfoldgrppin == Run, ]
      AbsTestData <- PseudoabsenceDat.df[kfoldgrpain == Run, ]
      # Combine test presence and absence data and  background data
      AbsTestData <- rbind(PresTestData, AbsTestData, BackgroundDat.df)
    }
    ########################
    ### Query projection grid to get values for test presence and absence points
    ## First presence points
    setwd(OutDirectIn)
    head(PresTestData)
    # Convert point data.frame to SpatialPointsDataFrame
    # First, specify xy coordinates
    xy <- PresTestData[,c("coords.x1", "coords.x2")]
    # Create spatial points data frame
    PresTestPnt.spdf <- SpatialPointsDataFrame(coords=xy, data=PresTestData, proj4string=CRS.In)
    ## Extract values of presence test points from GLM model grid
    system.time(prespred.df <- data.frame(extract(GLMBinomial.score, PresTestPnt.spdf)))
    colnames(prespred.df) <- "GLMScore"
    # Rejoin coordinates to prespred.df and save as shapefile for checking
    prespred.spdf <- SpatialPointsDataFrame(coords=xy, data=prespred.df, proj4string=CRS.In)
    # Write shapefile including the correct projection
    writeOGR(prespred.spdf, ".", paste0(Species, "PresenceTestGLMVals"), driver="ESRI Shapefile", overwrite=TRUE)
    ## Then absence points
    head(AbsTestData)
    # Convert point data.frame to SpatialPointsDataFrame
    # First, specify xy coordinates
    xy <- AbsTestData[,c("coords.x1", "coords.x2")]
    # Create spatial points data frame
    AbsTestPnt.spdf <- SpatialPointsDataFrame(coords=xy, data=AbsTestData, proj4string=CRS.In)
    ## Extract values of presence test points from GLM model grid
    system.time(abspred.df <- data.frame(extract(GLMBinomial.score, AbsTestPnt.spdf)))
    colnames(abspred.df) <- "GLMScore"
    # Rejoin coordinates to abspred.df and save as shapefile for checking
    abspred.spdf <- SpatialPointsDataFrame(coords=xy, data=abspred.df, proj4string=CRS.In)
    # Write shapefile including the correct projection
    writeOGR(abspred.spdf, ".", paste0(Species, "AbsenceTestGLMVals"), driver="ESRI Shapefile", overwrite=TRUE)
    #############################################################################################
    # This section evaluates the GLM model using the PresenceAbsence package
    #############################################################################################
    #### Create a dataset with model predictions for presence and absence points
    # Use extracted prediction values for presence and absence points for each of three
    # models previously calculated in loop
    #
    library(gtools)
    ## Create directory of output for class pair run
    #setwd(paste("C:/Users/JLTracy/Documents/R/win-library/3.0/10minClimIntEvalSTB1000", "/", output2, sep=""))
    # For presence data, assign a column of "1" under the name "OBSERVED" to indicate presence data
    # and assign the model results a name "GLMScoreN" where N is the name of the rep
    # Also assign a column "id" for the row numbers to use in merging the data frames later
    presa.df <- data.frame(c(rep(1, nrow(prespred.df))))
    names(prespred.df) <- c("GLM")
    names(presa.df) <- c("OBSERVED")
    pres.df <- data.frame(cbind(id=1:nrow(presa.df), presa.df, prespred.df))
    nrow(pres.df)
    # Repeat above process with absence data, but assign "OBSERVED" a value of 0
    absa.df <- data.frame(c(rep(0, nrow(abspred.df))))
    names(abspred.df) <- c("GLM")
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
    accusum <- matrix(data=NA, ncol=12, nrow=1, byrow=TRUE, dimnames=list(NULL, c("MaxTSS", "Specificity_maxTSS", "Sensitivity_maxTSS", "AUC", "MaxKappa", "ThresholdMaxTSS", "AICc_bg", "NumDVars", "NumEnVars", "EnvVarsUsed", "AICc", "AICc2")))
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
    accurun$Conditions <- paste("GLM", Subset, sep="")
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
    accusum[1,12] <- AICc2
    # Save Threshold value and multiply by 1000 for grid calibration
    ThresholdK <-  (max(ThresholdMaxTSS, 0.0001))*1000
    ###############################
    #endtime <- Sys.time()
    #durtime <- endtime - starttime
    # Save evaluation statistics to .csv file
    #setwd(paste("C:/Users/JLTracy/Documents/R/win-library/3.0/10minClimIntEvalSTB1000", "/", output2, sep=""))
    accusum.df <- data.frame(accusum, stringsAsFactors=FALSE)
    accusum.df$Model <- ModelType1
    #
    # Save variable matrix coordinates and Old World evaluation statistics to a vector
    GLMEvalStats.df <- data.frame(t(as.matrix(c(VariableNamesSel, accusum.df$EnvVarsUsed, SubsetVariableNumber, TestDataType, TotPres, SetName, Run, accusum.df$MaxTSS, accusum.df$MaxKappa, accusum.df$AUC, accusum.df$Specificity_maxTSS, accusum.df$Sensitivity_maxTSS, accusum.df$ThresholdMaxTSS, accusum.df$AICc_bg, accusum.df$NumDVars, accusum.df$NumEnVars, accusum.df$AICc, accusum.df$AICc2))))
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
    bmp(paste(Species, TestDataType, "GLMROC_Run", Run, ".bmp", sep=""))
    #dev.new()
    auc.roc.plot(DATA, color = TRUE, legend.cex = 1.2, main = "")
    dev.off()
    # Calculate optimal thresholds by various methods and save output to text file
    #outthresholds <- capture.output(optimal.thresholds(DATA, opt.methods = 1:12, threshold=10001))
    #out1 <- paste(outthresholds)
    #cat(out1, file=paste(Species, "GLMResults", TestDataType, "Thresholds.txt", sep=""), sep="\n", append=TRUE)
    # Can calculate the confusion matrix for a given threshold (not necessary)
    confmatrix <- cmx(DATA, threshold = ThresholdMaxTSS, na.rm = FALSE)
    # Save confusion matrix at threshold of maxTSS
    write.table(confmatrix, file=paste(Species, "GLMResults", TestDataType, "ConfusionMatrixatMaxTSSThreshold_Run", Run, ".csv", sep=""), sep=",", col.names=c("Obs Pres","Obs Abs"), row.names=c("Pred Pres", "Pred Abs"))
    #Save evaluation stastics to file
    write.table(accusum.df, file=paste(Species, "GLMResults", TestDataType, "Stats_Run", Run, ".csv", sep=""), sep=",", col.names=NA)
    graphics.off()
    # Save output evaluation statistics and summary statistics to text file
    # Output codes of variables used
    outA <- paste("\n\nSpecies:", Species, "\nVariable Subset:", VariableNamesSel, "\nModel:", Model, "\nScoring Alogorithm: GLM")
    cat(outA, file=paste(Species, "GLMResults", TestDataType, "Summary_Run", Run, ".txt", sep=""), sep="\n", append=TRUE)
    #NotUsedVars <- capture.output(print(DeleteBands))
    #outB <- paste("\nVariables Not Used: Code-", NotUsedVars, "; Variable-", GridNamesDrop)
    #cat(outB, file=paste(Species, InfClass, "GLMResults", output2, "Summary.txt", sep=""), sep="\n", append=TRUE)
    #UsedVars <- capture.output(print(SRandNumbers))
    #outB <- paste("\n  Code Nos. for Variables Used:", UsedVars)
    #cat(outB, file=paste("GLMResults", output2, "Summary.txt", sep=""), sep="\n", append=TRUE)
    #out1 <- paste("\nVariable Names")
    #cat(out1, file=paste(Species, "GLMResults", GridNamesKeepSelT, "Summary.txt", sep=""), sep="\n", append=TRUE)
    #EnvVariables <- as.matrix(GridNamesKeepSel)
    #out2 <- capture.output(print(EnvVariables))
    #cat(out2, file=paste(Species, "GLMResults", GridNamesKeepSelT, "Summary.txt", sep=""), sep="\n", append=TRUE)
    out3 <- paste("\n  Evaluation Statistics per Run")
    cat(out3, file=paste(Species, "GLMResults", TestDataType, "Summary_Run", Run, ".txt", sep=""), sep="\n", append=TRUE)
    out4 <- capture.output(print(accusum.df))
    cat(out4, file=paste(Species, "GLMResults", TestDataType, "Summary_Run", Run, ".txt", sep=""), sep="\n", append=TRUE)
    #
    #############
    GLMSubsetEvalStats.df1 <- GLMEvalStats.df
    head(GLMSubsetEvalStats.df1)
    tail(GLMSubsetEvalStats.df1)
    nrow(GLMSubsetEvalStats.df1)
    #str(ESWrapperOutput.df)
    colnames(GLMSubsetEvalStats.df1) <- c("VarNames", "EnvVarsUsed", "SubsetVariableNumber", "DataType", "TotalPresPnts", "SetName", "Run", "TSS", "Kappa", "AUC", "Spec", "Sens", "ThreshMxTSS", "AICc_bg", "NumDVars", "NumEnVars", "AICc", "AICc2")
    rownames(GLMSubsetEvalStats.df1) <- c(seq(1:nrow(GLMEvalStats.df)))
    GLMSubsetEvalStats.df1 <- as.data.frame(GLMSubsetEvalStats.df1, stingsAsFactors=FALSE)
    head(GLMSubsetEvalStats.df1)
    ncol(GLMSubsetEvalStats.df1)
    # Save output matrix
    # Convert sixth through 15th columns from character to numeric
    GLMSubsetEvalStats.df1[,7:18] <- sapply(GLMSubsetEvalStats.df1[,7:18], function(x) as.numeric(as.character(x)))
    if(k < 3) {
      GLMSubsetEvalL[[k]] <- GLMSubsetEvalStats.df1
    } else {
      GLMSubsetEvalL[[4]] <- GLMSubsetEvalStats.df1
    }
    #Calculate difference between training and test statistics (overfitting)
    if(k==2) {
      GLMSubsetEvalL[[3]] <- GLMSubsetEvalL[[2]]
      GLMSubsetEvalL[[3]][,8:18] <- GLMSubsetEvalL[[1]][,8:18] - GLMSubsetEvalL[[2]][,8:18]
      GLMSubsetEvalL[[3]][,4] <- "Diff"
    }
  }
  #####################################################################################
  GLMSubsetEvalStats.df <- do.call(rbind, GLMSubsetEvalL)
  #
  ################################################################################
  ### Calibrate model using using above calculated ThresholdK
  ### at maximum TSS from above evaluation 
  GLMBinomial.scorecal <- calc(GLMBinomial.score1, function(x) ifelse(x < ThresholdK, 0, 1) )
  #plot(GLM.scorecal)
  # SubsetVariableNumber <- 19
  #OutGridID <- "Full"
  setwd(OutDirectIn)
  writeRaster(GLMBinomial.scorecal, paste0(Species, "GLMBinomial", OutGridID, SubsetVariableNumber, "Cal"), format = "GTiff", overwrite=TRUE)
  ##############
  #
  if(Output==TRUE) {
    setwd(OutDirectIn)
    Sets <- nrow(GLMSubsetEvalStats.df)
    write.table(GLMSubsetEvalStats.df, file=paste(Species, "GLM_TrainTest_", TotVars, "_TotVars_", SubsetVariableNumber, "Var_", Sets, Subset, "_", SetName, "_", Run, ".csv", sep=""), sep=",", col.names=NA)
  }
  #
  gc()
  return(GLMSubsetEvalStats.df)
}
