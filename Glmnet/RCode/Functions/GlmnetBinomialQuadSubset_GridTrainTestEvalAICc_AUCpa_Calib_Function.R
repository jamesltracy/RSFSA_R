GlmnetBinomialQuadSubset.GridTrainTestEvalAICc_AUCpa_Calib <-
function(Species, VariableNamesIn, PredictorsIn, SubsetVariableNumber, TotPres, Run, SetName, VariableSubset, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrppin, kfoldgrpain, cvglmnetfoldsin, CRS.In, OutDirectIn, FunctDirectIn..., Output, fitAICModel) {
  #kfoldgrppain <- NObskfoldgrppa
  #NOTE: Function parameters after the "...", such as TestDataType, have to be set with an equal sign in the function call, such as ScoreType="mean"
  #FunctDirectIn <- "C:/Users/James/Documents/R/win-library/"
  setwd(FunctDirectIn)
  # Read in User Defined Functions
  source("cvglmnet_PredRaster_BinomialQuadResp_Function.R")
  #
  setwd(OutDirectIn)
  library(glmnet)
  CVGlmnetRuns <- 1
  if(missing(Output)) {Output=TRUE}
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
  ModelType1 <- paste("glmnetBinomialQuad", Subset, sep="")
  Model <- paste("glmnet Binomial-Quadratic for", SubsetVariableNumber, "Feature Subset")
  tail(VariableSubset)
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
  #######################################
  ## Designate final training data and add count column as factor
  GlmnetPresTrainData <- PresenceDat.df[kfoldgrppin != Run,]
  GlmnetPresTrainData$Count <- 1
  head(GlmnetPresTrainData)
  #str(GlmnetPresTrainData)
  GlmnetAbsTrainData <- PseudoabsenceDat.df[kfoldgrpain != Run,]
  GlmnetAbsTrainData$Count <- 0
  ############################################################################
  ## Loop through different kfold groups with equal prevalence (% presence data)
  ## and save coefficients and lambda.min
  CoefficientsList <- list()
  ModelStatsList <- list()
  ModelStats.df <- data.frame(t(as.matrix(c(0,0,0,0))))
  colnames(ModelStats.df) <- c("ModelNo", "deviance.min", "lambda.min", "AICc2")
  ######################
  ta <- Sys.time()
  for(q in 1:CVGlmnetRuns) {
    #q=1
    # Create kfold groups with equal prevalence (percent presence) for
    # presence and absence portions of PresAbsDatTrain.df
    NTrainkfoldgrpp <- kfold(GlmnetPresTrainData, cvglmnetfoldsin)
    NTrainkfoldgrpa <- kfold(GlmnetAbsTrainData, cvglmnetfoldsin)
    # Join presence and absence data and join kfold groups
    PresAbsTrainDat.df <- rbind(GlmnetPresTrainData, GlmnetAbsTrainData)
    head(PresAbsTrainDat.df)
    tail(PresAbsTrainDat.df)
    nrow(PresAbsTrainDat.df)
    ###
    NTrainkfoldgrppa <- c(NTrainkfoldgrpp, NTrainkfoldgrpa)
    length(NTrainkfoldgrppa)
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
    xTrainData.df <- x
    head(xTrainData.df)
    nrow(xTrainData.df)
    ######
    ## Add raw quadratic versions for main effect variables (no interaction)
    xTrainDataQ.df <- cbind(xTrainData.df, xTrainData.df)
    colcount <- 1
    for(i in 1:ncol(xTrainData.df)) {
      #i=3
        # Use poly function to obtained orthonalized main effect and squared quadratic value
        xTrainDataQ.df[,colcount:(colcount+1)] <- poly(xTrainData.df[,i], 2, raw=TRUE)
        head(xTrainDataQ.df)
        colnames(xTrainDataQ.df)[colcount:(colcount+1)] <- c(colnames(xTrainData.df)[i], paste0(colnames(xTrainData.df)[i], "^2"))
        colcount <- colcount + 2
    }
    head(xTrainDataQ.df)
    xTrainDataQ <- as.matrix(xTrainDataQ.df)
    head(xTrainDataQ)
    nrow(xTrainDataQ)
    #####################
    # Identify response variable for training data
    yTrainData <- PresAbsTrainDat.df$Count
    length(yTrainData)
    #
    # Run cv glmnet binomial model for response probability of presence using cross validation
    # to select the model coefficients with the minimum deviance (default)
    #
    cv_glmnetBinomialModel = cv.glmnet(xTrainDataQ, yTrainData, foldid=NTrainkfoldgrppa, family = "binomial", nfolds = cvglmnetfoldsin, alpha=alphain)
    # If CVGlmnetRuns is one, save model
    if(CVGlmnetRuns==1) {
      # Save Glmnet model
      setwd(OutDirectIn)
      saveRDS(cv_glmnetBinomialModel, "cv_glmnetBinomialQuadModel.rds")
    }
    #
    # Retrieve coefficients from selected Binomial model
    s_in = "lambda.min"
    Coefficients.df1 <- data.frame(as.matrix(coef(cv_glmnetBinomialModel, s=s_in)), stringsAsFactors=FALSE)
    colnames(Coefficients.df1) <- "Coefficients"
    #str(cv_glmnetBinomialModel)
    lambda.min <- cv_glmnetBinomialModel$lambda.min
    deviance.min <- min(cv_glmnetBinomialModel$cvm)
    ####################################################
    ###
    # Caculate AICc using formula from Johnny Heineken at  https://stats.stackexchange.com/questions/25817/is-it-possible-to-calculate-aic-and-bic-for-lasso-regression-models
    n <- cv_glmnetBinomialModel$glmnet.fit$nobs
    ## cv.glmnet deviance is divided by number of observation (see https://stackoverflow.com/questions/43468665/poisson-deviance-glmnet), so need to multiply by n
    deviance.fit <- n*(cv_glmnetBinomialModel$cvm[match(cv_glmnetBinomialModel$lambda.min, cv_glmnetBinomialModel$lambda)])
    fit.nulldev <- cv_glmnetBinomialModel$glmnet.fit$nulldev
    k <- cv_glmnetBinomialModel$glmnet.fit$df[match(cv_glmnetBinomialModel$lambda.min, cv_glmnetBinomialModel$lambda)]
    tLL <- fit.nulldev - deviance.fit
    AICc2 <- -tLL+2*k+2*k*(k+1)/(n-k-1)
    ###
    ModelStats.df[1,1] <- q
    ModelStats.df[1,2] <- deviance.min
    ModelStats.df[1,3] <- lambda.min
    ModelStats.df[1,4] <- AICc2
    CoefficientsList[[q]] <- Coefficients.df1
    ModelStatsList[[q]] <- ModelStats.df
  }
  tb <- Sys.time()
  difftime(tb, ta, units = "mins")
  #############################
  #
  ModelStatsResults.df <- do.call(rbind,ModelStatsList)
  ModelStatsResults.df <- ModelStatsResults.df[order(ModelStatsResults.df$deviance.min),]
  nrow(ModelStatsResults.df)
  CoefficientSets.df <- do.call(cbind, CoefficientsList)
  AvgCoefficients.df1 <- data.frame(apply(CoefficientSets.df,1,mean))
  colnames(AvgCoefficients.df1) <- "Coefficients"
  # Retrieve variable names
  VarNamesGLM <- data.frame(rownames(AvgCoefficients.df1), stringsAsFactors=FALSE)
  # Join variable names with coefficients
  AvgCoefficients.df <- cbind(VarNamesGLM, AvgCoefficients.df1)
  colnames(AvgCoefficients.df) <- c("VarNames", "Coefficients")
  # Remove VarNames column
  AvgCoefficients.df <- AvgCoefficients.df[, 2, drop=FALSE]
  #
  AICc2 <- mean(ModelStatsResults.df$AICc2)
  #cv_glmnetBinomialModel <- readRDS("cv_glmnetBinomialQuadModel.rds")
  #####################################################
  #Save coefficients
  setwd(OutDirectIn)
  write.table(AvgCoefficients.df, file=paste("Coefficients", ModelType1, "_", SubsetVariableNumber, "Var_", Loop, ".csv", sep=""), sep=",", col.names=NA)
  #
  VarNamesUsed <- as.vector(rownames(AvgCoefficients.df)[2:nrow(AvgCoefficients.df)])
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
  #checkModel <- readRDS("cv_glmnetBinomialQuadModel.rds")
  #coef(checkModel, s= "lambda.min")
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
  #####################################
  ## Calculate AICc_bg with point values
  # Obtain number of parameters in Glmnet model
  CoefficientsObj <- coef(cv_glmnetBinomialModel, s = s_in)
  nparams <- (length(CoefficientsObj@x) - 1)
  #
  ## From ENMeval Package documentation: AICc is the Akaike Information Criterion corrected for small
  ## sample sizes calculated as: (2 * K - 2 * logLikelihood) + (2 * K) * (K + 1)=(n - K - 1)
  ## where K is the number of parameters in the model (i.e., number of non-zero parameters in Glmnet
  ## lambda file) and n is the number of occurrence localities.
  ## The logLikelihood is sum(log(vals/total))
  ## vals is vector of Glmnet raw values at occurence localities
  ## total is the sum of Glmnet raw values across the entire study area
  ##
  # Obtain values for training presence data output by glmnet run
  xPresTrainDataQ <- xTrainDataQ[1:nrow(GlmnetPresTrainData),]
  # Select only terms used in ensemble model
  xPresDataQ <- xPresTrainDataQ[, c(VarNamesUsed), drop=FALSE]
  #str(xPresDataQ)
  head(xPresDataQ)
  presvals.df <- data.frame(predict(cv_glmnetBinomialModel, newx = xPresDataQ, type = "response", s = s_in))
  #str(xPresTrainDataQ)
  # total is sum of Glmnet raw values across entire study area, includes background and occurrence localities
  # Obtain values for training presence/absence data output by glmnet run
  BackgroundDat.df2 <- subset(BackgroundDat.df, select=c(VarNames))
  head(BackgroundDat.df2)
  ######
  ## Add raw quadratic versions for main effect variables of background (no interaction)
  BackgroundDatQ.df2 <- cbind(BackgroundDat.df2, BackgroundDat.df2)
  head(BackgroundDatQ.df2)
  colcount <- 1
  for(i in 1:ncol(BackgroundDat.df2)) {
    #i=1
      # Use poly function to obtained orthonalized main effect and squared quadratic value
      BackgroundDatQ.df2[,colcount:(colcount+1)] <- poly(BackgroundDat.df2[,i], 2, raw=TRUE)
      head(BackgroundDatQ.df2)
      colnames(BackgroundDatQ.df2)[colcount:(colcount+1)] <- c(colnames(BackgroundDat.df2)[i], paste0(colnames(BackgroundDat.df2)[i], "^2"))
      colcount <- colcount + 2
  }
  head(BackgroundDatQ.df2)
  BackgroundDatQ <- as.matrix(BackgroundDatQ.df2)
  head(BackgroundDatQ)
  nrow(BackgroundDatQ)
  #####################
  # Select only terms used in ensemble model
  xDataQ <- BackgroundDatQ[, c(VarNamesUsed), drop=FALSE]
  head(xDataQ)
  #
  backgroundvals.df <- data.frame(predict(cv_glmnetBinomialModel, newx = xDataQ, type = "response", s = s_in))
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
  vals <- PresBckgrndDataRAW.df[1:nrow(GlmnetPresTrainData),1]
  head(vals)
  n <- length(vals) # number of occurence localities
  # Keep values for background  data as vector
  backgroundvals <- PresBckgrndDataRAW.df[(nrow(GlmnetPresTrainData)+1):nrow(PresBckgrndDataRAW.df),1]
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
  #############################################################################
  #str(GlmnetOut)
  # Copy Glmnet output to desired directory
  #flist <- list.files(GlmnetOut@path, full.names = TRUE)
  #file.copy(flist, OutDirectSub, overwrite=TRUE)
  ##
  #####################################
  ## Run projection of model
  # Assemble environmental rasters
  PredictorsIn <- Predictors[[VarNames]]
  # Project model onto raster stack of environmental variables
  ## Took 1.6 hrs for four of 10 variables over funnel
  system.time(glmnetBinomial.score <- cvglmnet_PredRaster_BinomialQuadResp(cv_glmnetBinomialModel, PredictorsIn, Loop, OutDirectIn, s = s_in))
  # Multiply continuous model grid by 1000 and convert to integer
  system.time(glmnetBinomial.score1 <- calc(glmnetBinomial.score, function(x) as.integer(x * 1000) ))
  # Save grid
  setwd(OutDirectIn)
  writeRaster(glmnetBinomial.score1, paste0(Species, "glmnetBinomial", OutGridID, SubsetVariableNumber), format = "GTiff", overwrite=TRUE)
  # glmnetBinomial.score1 <- raster(paste0(Species, "glmnetBinomial", OutGridID, SubsetVariableNumber, ".tif"))
  # str(glmnetBinomial.score)
  # glmnetBinomial.score <- calc(glmnetBinomial.score1, function(x) x/1000)
  ########################################################################
  GlmnetSubsetEvalL <- list()
  TestDataTypes <- c("Train", "Test", "Test_bgp")
  for(k in 1:3) {
    #k=1
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
    ## Extract values of presence test points from Glmnet model grid
    system.time(prespred.df <- data.frame(extract(glmnetBinomial.score, PresTestPnt.spdf)))
    colnames(prespred.df) <- "GlmnetScore"
    head(prespred.df)
    #mean(prespred.df[,1])
    # Rejoin coordinates to prespred.df and save as shapefile for checking
    prespred.spdf <- SpatialPointsDataFrame(coords=xy, data=prespred.df, proj4string=CRS.In)
    # Write shapefile including the correct projection
    writeOGR(prespred.spdf, ".", paste0(Species, "PresenceTestGlmnetVals"), driver="ESRI Shapefile", overwrite=TRUE)
    ## Then absence points
    head(AbsTestData)
    # Convert point data.frame to SpatialPointsDataFrame
    # First, specify xy coordinates
    xy <- AbsTestData[,c("coords.x1", "coords.x2")]
    # Create spatial points data frame
    AbsTestPnt.spdf <- SpatialPointsDataFrame(coords=xy, data=AbsTestData, proj4string=CRS.In)
    ## Extract values of presence test points from Glmnet model grid
    system.time(abspred.df <- data.frame(extract(glmnetBinomial.score, AbsTestPnt.spdf)))
    colnames(abspred.df) <- "GlmnetScore"
    #mean(abspred.df[,1])
    # Rejoin coordinates to abspred.df and save as shapefile for checking
    abspred.spdf <- SpatialPointsDataFrame(coords=xy, data=abspred.df, proj4string=CRS.In)
    # Write shapefile including the correct projection
    writeOGR(abspred.spdf, ".", paste0(Species, "PresenceTestGlmnetVals"), driver="ESRI Shapefile", overwrite=TRUE)
    #############################################################################################
    # This section evaluates the Glmnet model using the PresenceAbsence package
    #############################################################################################
    #### Create a dataset with model predictions for presence and absence points
    # Use extracted prediction values for presence and absence points for each of three
    # models previously calculated in loop
    #
    library(gtools)
    ## Create directory of output for class pair run
    #setwd(paste("C:/Users/JLTracy/Documents/R/win-library/3.0/10minClimIntEvalSTB1000", "/", output2, sep=""))
    # For presence data, assign a column of "1" under the name "OBSERVED" to indicate presence data
    # and assign the model results a name "GlmnetScoreN" where N is the name of the rep
    # Also assign a column "id" for the row numbers to use in merging the data frames later
    presa.df <- data.frame(c(rep(1, nrow(prespred.df))))
    names(prespred.df) <- c("Glmnet")
    names(presa.df) <- c("OBSERVED")
    pres.df <- data.frame(cbind(id=1:nrow(presa.df), presa.df, prespred.df))
    nrow(pres.df)
    # Repeat above process with absence data, but assign "OBSERVED" a value of 0
    absa.df <- data.frame(c(rep(0, nrow(abspred.df))))
    names(abspred.df) <- c("Glmnet")
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
    accusum <- matrix(data=NA, ncol=11, nrow=1, byrow=TRUE, dimnames=list(NULL, c("MaxTSS", "Specificity_maxTSS", "Sensitivity_maxTSS", "AUC", "MaxKappa", "ThresholdMaxTSS", "AICc_bg", "NumDVars", "NumEnVars", "EnvVarsUsed", "AICc2")))
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
    accurun$Conditions <- paste("Glmnet", Subset, sep="")
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
    accusum[1,11] <- AICc2
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
    GlmnetEvalStats.df <- data.frame(t(as.matrix(c(VariableNamesSel, accusum.df$EnvVarsUsed, SubsetVariableNumber, TestDataType, TotPres, SetName, Run, accusum.df$MaxTSS, accusum.df$MaxKappa, accusum.df$AUC, accusum.df$Specificity_maxTSS, accusum.df$Sensitivity_maxTSS, accusum.df$ThresholdMaxTSS, accusum.df$AICc_bg, accusum.df$NumDVars, accusum.df$NumEnVars, accusum.df$AICc2))))
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
    bmp(paste(Species, TestDataType, "GlmnetROC_Run", Run, ".bmp", sep=""))
    #dev.new()
    auc.roc.plot(DATA, color = TRUE, legend.cex = 1.2, main = "")
    dev.off()
    # Calculate optimal thresholds by various methods and save output to text file
    #outthresholds <- capture.output(optimal.thresholds(DATA, opt.methods = 1:12, threshold=10001))
    #out1 <- paste(outthresholds)
    #cat(out1, file=paste(Species, "GlmnetResults", TestDataType, "Thresholds.txt", sep=""), sep="\n", append=TRUE)
    # Can calculate the confusion matrix for a given threshold (not necessary)
    confmatrix <- cmx(DATA, threshold = ThresholdMaxTSS, na.rm = FALSE)
    # Save confusion matrix at threshold of maxTSS
    write.table(confmatrix, file=paste(Species, "GlmnetResults", TestDataType, "ConfusionMatrixatMaxTSSThreshold_Run", Run, ".csv", sep=""), sep=",", col.names=c("Obs Pres","Obs Abs"), row.names=c("Pred Pres", "Pred Abs"))
    #Save evaluation stastics to file
    write.table(accusum.df, file=paste(Species, "GlmnetResults", TestDataType, "Stats_Run", Run, ".csv", sep=""), sep=",", col.names=NA)
    graphics.off()
    # Save output evaluation statistics and summary statistics to text file
    # Output codes of variables used
    outA <- paste("\n\nSpecies:", Species, "\nVariable Subset:", VariableNamesSel, "\nModel:", Model, "\nScoring Alogorithm: Glmnet")
    cat(outA, file=paste(Species, "GlmnetResults", TestDataType, "Summary_Run", Run, ".txt", sep=""), sep="\n", append=TRUE)
    #NotUsedVars <- capture.output(print(DeleteBands))
    #outB <- paste("\nVariables Not Used: Code-", NotUsedVars, "; Variable-", GridNamesDrop)
    #cat(outB, file=paste(Species, InfClass, "GlmnetResults", output2, "Summary.txt", sep=""), sep="\n", append=TRUE)
    #UsedVars <- capture.output(print(SRandNumbers))
    #outB <- paste("\n  Code Nos. for Variables Used:", UsedVars)
    #cat(outB, file=paste("GlmnetResults", output2, "Summary.txt", sep=""), sep="\n", append=TRUE)
    #out1 <- paste("\nVariable Names")
    #cat(out1, file=paste(Species, "GlmnetResults", GridNamesKeepSelT, "Summary.txt", sep=""), sep="\n", append=TRUE)
    #EnvVariables <- as.matrix(GridNamesKeepSel)
    #out2 <- capture.output(print(EnvVariables))
    #cat(out2, file=paste(Species, "GlmnetResults", GridNamesKeepSelT, "Summary.txt", sep=""), sep="\n", append=TRUE)
    out3 <- paste("\n  Evaluation Statistics per Run")
    cat(out3, file=paste(Species, "GlmnetResults", TestDataType, "Summary_Run", Run, ".txt", sep=""), sep="\n", append=TRUE)
    out4 <- capture.output(print(accusum.df))
    cat(out4, file=paste(Species, "GlmnetResults", TestDataType, "Summary_Run", Run, ".txt", sep=""), sep="\n", append=TRUE)
    #
    #############
    GlmnetSubsetEvalStats.df1 <- GlmnetEvalStats.df
    head(GlmnetSubsetEvalStats.df1)
    tail(GlmnetSubsetEvalStats.df1)
    nrow(GlmnetSubsetEvalStats.df1)
    #str(ESWrapperOutput.df)
    colnames(GlmnetSubsetEvalStats.df1) <- c("VarNames", "EnvVarsUsed", "SubsetVariableNumber", "DataType", "TotalPresPnts", "SetName", "Run", "TSS", "Kappa", "AUC", "Spec", "Sens", "ThreshMxTSS", "AICc_bg", "NumDVars", "NumEnVars", "AICc2")
    rownames(GlmnetSubsetEvalStats.df1) <- c(seq(1:nrow(GlmnetEvalStats.df)))
    GlmnetSubsetEvalStats.df1 <- as.data.frame(GlmnetSubsetEvalStats.df1, stingsAsFactors=FALSE)
    head(GlmnetSubsetEvalStats.df1)
    ncol(GlmnetSubsetEvalStats.df1)
    # Save output matrix
    # Convert sixth through 15th columns from character to numeric
    GlmnetSubsetEvalStats.df1[,7:17] <- sapply(GlmnetSubsetEvalStats.df1[,7:17], function(x) as.numeric(as.character(x)))
    if(k < 3) {
      GlmnetSubsetEvalL[[k]] <- GlmnetSubsetEvalStats.df1
    } else {
      GlmnetSubsetEvalL[[4]] <- GlmnetSubsetEvalStats.df1
    }
    #Calculate difference between training and test statistics (overfitting)
    if(k==2) {
      GlmnetSubsetEvalL[[3]] <- GlmnetSubsetEvalL[[2]]
      GlmnetSubsetEvalL[[3]][,8:17] <- GlmnetSubsetEvalL[[1]][,8:17] - GlmnetSubsetEvalL[[2]][,8:17]
      GlmnetSubsetEvalL[[3]][,4] <- "Diff"
    }
  }
  #####################################################################################
  GlmnetSubsetEvalStats.df <- do.call(rbind, GlmnetSubsetEvalL)
  #
  ### Calculate AICc_bg_reg save to data set
  ## Read in the regression relationship between AICc_bg and AICc from an AUC selection run
  if(!missing(fitAICModel)) {
    ## Calculate AICc_bg_reg using regression equation
    xDatIn <- data.frame(AICc2 = GlmnetSubsetEvalStats.df$AICc2)
    head(xDatIn)
    GlmnetSubsetEvalStats.df$AICc_bg_reg <- predict(fitAIC, newdata=xDatIn)
  }
  #
  ################################################################################
  ### Calibrate model using using above calculated ThresholdK
  ### at maximum TSS from above evaluation 
  glmnetBinomial.scorecal <- calc(glmnetBinomial.score1, function(x) ifelse(x < ThresholdK, 0, 1) )
  #plot(Glmnet.scorecal)
  # SubsetVariableNumber <- 19
  #OutGridID <- "Full"
  setwd(OutDirectIn)
  writeRaster(glmnetBinomial.scorecal, paste0(Species, "GlmnetBinomial", OutGridID, SubsetVariableNumber, "Cal"), format = "GTiff", overwrite=TRUE)
  ##############
  #
  if(Output==TRUE) {
    setwd(OutDirectIn)
    Sets <- nrow(GlmnetSubsetEvalStats.df)
    write.table(GlmnetSubsetEvalStats.df, file=paste(Species, "Glmnet_TrainTest_", TotVars, "_TotVars_", SubsetVariableNumber, "Var_", Run, ".csv", sep=""), sep=",", col.names=NA)
  }
  #
  gc()
  return(GlmnetSubsetEvalStats.df)
}
