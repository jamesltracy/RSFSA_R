GLMBinomialMultiSubset.WrapperTrainTestEvalAICc <-
function(Species, VariableNamesIn, SubsetVariableNumber, TotPres, VariableSubsetsIn, SetRunID, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrpp, kfoldgrpa, OutDirectIn, FunctDirectIn..., Output, DataSetType) {
  setwd(OutDirectIn)
  library(AICcmodavg)
  library(dismo)
  library(foreach)
  library(doParallel)
  # Make sure VariableNames is a data frame
  VariableNamesIn <- data.frame(VariableNamesIn, stringsAsFactors=FALSE)
  #
  RowNames.df <- data.frame(rownames(VariableSubsetsIn), stringsAsFactors=FALSE)  # Save row names to use in output
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
  if(missing(DataSetType)) { DataSetType="" }
  if(missing(Output)) { Output=TRUE }
  #
  ModelType1 <- paste("GLMBinomial", Subset, sep="")
  Model <- paste("GLMBinomial for", SubsetVariableNumber, "Feature Subset")
  tail(VariableSubsetsIn)
  #VariableSubsetsIn[2998,]
  # Register cluster for running parallel with doParallel pkg
  # Use a minimum of one core or the number of available cores minus one
  cores <- detectCores() # detect cores running
  workcores <- max(1, cores-1)
  cl <- makeCluster(workcores)  # make cluster of workcores
  #cl <- makeCluster(1)
  registerDoParallel(cl) # register cluster
  getDoParWorkers() # check clusters registered
  #registerDoSEQ()
  #
  if(nrow(VariableSubsetsIn)<2) {
    ivec <- 1
  } else {
    ivec <- c(seq(1, nrow(VariableSubsetsIn), 1))
  }
  #ivec <- c(1,2,3)
  #ivec <- c(1,2,3,4)
  #ivec <- c(1)
  ###
  # Set up foreach to output to ten different list elements corresponding to
  # the two Old World and New World sets of five evaluation statistics
  ## Designate wrapper training data and add count column as factor
  GLMPresTrainData <- PresenceDat.df[kfoldgrpp==1,]
  GLMPresTrainData$Count <- 1
  head(GLMPresTrainData)
  #str(GLMPresTrainData)
  GLMAbsTrainData <- PseudoabsenceDat.df[kfoldgrpp==1,]
  GLMAbsTrainData$Count <- 0
  ##############################################################
  t1 <- Sys.time()
  #VariablesSubsets <- VariableKeepSubsets
  GLMSubsetEvalStats <- foreach(i=ivec, .combine='rbind') %dopar%  {
    library(dismo)
    library(raster)
    library(AICcmodavg)
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
    ## Define function to predict matrix with GLM
    GLMBinomialResp_PredictMatrix <-
    function(CoefficientsIn, xData.dfIn) {
      #CoefficientsIn <- AvgCoefficients.df
      #xDataIn <- PresAbsDat.df
      # Define function makeglm.R: Creates a "fake" glm object with specific coefficients that you can use for predicting without fitting a model first
      # https://gist.github.com/MrFlick/ae299d8f3760f02de6bf
      ##
      makeglm <- function(formula, ..., family, data=NULL) {
          dots <- list(...)
          out<-list()
          tt <- terms(formula, data=data)
          if(!is.null(data)) {
              mf <- model.frame(tt, data)
              vn <- sapply(attr(tt, "variables")[-1], deparse)

              if((yvar <- attr(tt, "response"))>0)
                  vn <- vn[-yvar]
                  xlvl <- lapply(data[vn], function(x) if (is.factor(x))
                 levels(x)
              else if (is.character(x))
                 levels(as.factor(x))
              else
                  NULL)
              attr(out, "xlevels") <- xlvl[!vapply(xlvl,is.null,NA)]
              attr(tt, "dataClasses") <- sapply(data[vn], stats:::.MFclass)
          }
          out$terms <- tt
          coef <- numeric(0)
          stopifnot(length(dots)>1 & !is.null(names(dots)))
          for(i in seq_along(dots)) {
              if((n<-names(dots)[i]) != "") {
                  v <- dots[[i]]
                  if(!is.null(names(v))) {
                      coef[paste0(n, names(v))] <- v
                  } else {
                      stopifnot(length(v)==1)
                      coef[n] <- v
                  }
              } else {
                  coef["(Intercept)"] <- dots[[i]]
              }
          }
          out$coefficients <- coef
          out$rank <- length(coef)
          if (!missing(family)) {
              out$family <- if (class(family) == "family") {
                  family
              } else if (class(family) == "function") {
                  family()
              } else if (class(family) == "character") {
                  get(family)()
              } else {
                  stop(paste("invalid family class:", class(family)))
              }
              out$qr <- list(pivot=seq_len(out$rank))
              out$deviance <- 1
              out$null.deviance <- 1
              out$aic <- 1
              class(out) <- c("glm","lm")
          } else {
              class(out) <- "lm"
              out$fitted.values <- predict(out, newdata=dd)
              out$residuals <- out$mf[attr(tt, "response")] - out$fitted.values
              out$df.residual <- nrow(data) - out$rank
              out$model <- data
              #QR doesn't work
          }
          out
      }
      ###################################################
      setwd(OutDirectIn)
      ##
      # Retrieve variables in regression
      VarNamesUsed <- as.vector(rownames(CoefficientsIn)[2:nrow(CoefficientsIn)])
      # Make sure input data matches VarNamesUsed
      head(xData.dfIn)
      nrow(xData.dfIn)
      xData.dfIn <- subset(xData.dfIn, select=VarNamesUsed)
      # Replace "^" with "P"
      colnames(xData.dfIn) <- gsub("\\^2", "P2", colnames(xData.dfIn))
      ## Prepare to construct regression formula from CoefficientsIn
      Coefficients.mat1 <- t(as.matrix(CoefficientsIn))
      Coefficients.df <- data.frame(Coefficients.mat1, stringsAsFactors=FALSE)
      colnames(Coefficients.df) <- colnames(Coefficients.mat1)
      colnames(Coefficients.df) <- gsub("\\^2", "P2", colnames(Coefficients.df))
      #str(Coefficients.df)
      ## Create dummy fake training data for binomial GLM to use in makeglm function
      ## for assigning coefficients to GLM model
      FakeYTrainData <- sample(0:1, nrow(xData.dfIn), replace=T)
      # Prepare GLM formula for makeglm
      fla <- paste("FakeYTrainData ~", paste(colnames(Coefficients.df[2:ncol(Coefficients.df)]), collapse=" + "))
      ###########
      ## Create character representation of coefficients for makeglm function including intercept and variables with coefficients
      CoefficientsVector <- c(seq(1:ncol(Coefficients.df)))
      for(i in CoefficientsVector) {
        #i=2
        if(i==1) {
          CoefficientsVector[i] <- paste(Coefficients.df[1,1])
        } else {
          CoefficientsVector[i] <- paste0(colnames(Coefficients.df)[i], "=", Coefficients.df[1,i])
        }
      }
      ###########
      CoefForm <- paste(CoefficientsVector, collapse=", ")
      ## Use makeglm function to calculate dummy GLM and assign coefficients
      GLMBinomialModel <- eval(parse(text = paste("makeglm(as.formula(fla), family=binomial, data=xData.dfIn,", CoefForm, ")")))
      ## Use GLM model with assigned coefficients to predict input matrix data values
      GLMscore.df <- data.frame("GLMScore" = predict(GLMBinomialModel, newdata=xData.dfIn, type="response"))
      head(GLMscore.df)
      ##
      return(GLMscore.df)
    }
    ###############
    # Split VarNames of VariableSubsetsIn into separate variables
    #str(VariableSubsets)
    VariableNamesSel <- c(unlist(VariableSubsetsIn[i,]))
    # Split VarNames of VariableSubsetsIn into separate variables
    if(grepl("-",VariableNamesSel)) {
      VarNames <- unlist(strsplit(VariableNamesSel, "-"))
    } else {
      VarNames <- unlist(VariableNamesSel)
    }
    SubsetVarNum <- length(VarNames)
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
    xPresTrainData.df <- xPresTrainData.df1[, c(VarsUsed), drop=FALSE]
    head(xPresTrainData.df)
    #
    # Select only terms used in ensemble model
    xPresData.dfIn <- xPresTrainData.df[, c(VarNamesUsed), drop=FALSE]
    head(xPresData.dfIn)
    # Obtain values for training presence data output by GLM run
    CoefficientsIn <- Coefficients.df
    presvals.df <- setNames(data.frame(predict(GLMBinomialModel, xPresData.dfIn, type="response")), "GLMScore")
    head(presvals.df)
    #presvals.df <- GLMBinomialResp_PredictMatrix(CoefficientsIn, xPresData.dfIn)
    #str(xTrainData)
    # Obtain values for background data output by GLM run
    # Select only terms used in ensemble model
    xData.dfIn <- BackgroundDat.df[, c(VarNamesUsed), drop=FALSE]
    head(xData.dfIn)
    backgroundvals.df <- setNames(data.frame(predict(GLMBinomialModel, xData.dfIn, type="response")), "GLMScore")
    #backgroundvals.df <- cvGLMEnsemble_PredictMatrix_BinomialResp(CoefficientsIn, xData.dfIn)
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
    #############################################################################
    ########### Evaluate model using training and test data
    GLMSubsetEvalL <- list()
    TestDataTypes <- c("WrapperTrain", "WrapperTest")
    for(TestDataType in TestDataTypes) {
      #TestDataType="WrapperTrain"
      if(TestDataType=="WrapperTrain") {
        k <- 1
      } else {
        k <- 2
      }
      # Specify model testing data
      GLMPresTestData <- PresenceDat.df[kfoldgrpp == k, ]
      head(GLMPresTestData)
      nrow(GLMPresTestData)
      ###
      GLMAbsTestData <- PseudoabsenceDat.df[kfoldgrpa == k, ]
      ########################
      ### Predict model values for test matrix to get values for test presence and absence points
      ## First presence points
      head(GLMPresTestData)
      # Subset PresTestData by VarNames
      PresTestData.df <- GLMPresTestData[, c(VarsUsed), drop=FALSE]
      head(PresTestData.df)
      #
      # Select only terms used in ensemble model
      xData.dfIn <- PresTestData.df[, c(VarNamesUsed), drop=FALSE]
      head(xData.dfIn)
      #CoefficientsIn <- Coefficients.df
      prespred.df <- setNames(data.frame(predict(GLMBinomialModel, xData.dfIn, type="response")), "GLMScore")
      #prespred.df <- GLMBinomialResp_PredictMatrix(CoefficientsIn, xData.dfIn)
      #colnames(prespred.df) <- "GLMScore"
      ## Then absence points
      head(GLMAbsTestData)
      # Subset PresTestData by VarNames
      AbsTestData.df <- GLMAbsTestData[, c(VarsUsed), drop=FALSE]
      head(AbsTestData.df)
      #
      # Select only terms used in ensemble model
      xData.dfIn <- AbsTestData.df[, c(VarNamesUsed), drop=FALSE]
      head(xData.dfIn)
      #
      abspred.df <- setNames(data.frame(predict(GLMBinomialModel, xData.dfIn, type="response")), "GLMScore")
      #abspred.df <- GLMBinomialResp_PredictMatrix(CoefficientsIn, xData.dfIn)
      #colnames(abspred.df) <- "GLMScore"
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
      accusum <- matrix(data=NA, ncol=12, nrow=1, byrow=TRUE, dimnames=list(NULL, c("MaxTSS", "Specificity_maxTSS", "Sensitivity_maxTSS", "AUC", "MaxKappa", "ThresholdMaxTSS", "AICc_bg", "AICc", "AICc2",  "NumDVars", "NumEnVars", "EnvVarsUsed")))
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
      accusum[1,8] <- AICc
      accusum[1,9] <- AICc2
      accusum[1,10] <- NumDVars
      accusum[1,11] <- nenvars
      accusum[1,12] <- EnvVarsUsed
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
      GLMSubsetEvalL[[k]] <- c(VariableNamesSel, accusum.df$EnvVarsUsed, SubsetVariableNumber, TestDataType, accusum.df$MaxTSS, accusum.df$MaxKappa, accusum.df$AUC, accusum.df$Specificity_maxTSS, accusum.df$Sensitivity_maxTSS, accusum.df$ThresholdMaxTSS, accusum.df$AICc_bg, accusum.df$AICc, accusum.df$AICc2, accusum.df$NumDVars, accusum.df$NumEnVars)
      #
    }
    return(GLMSubsetEvalL)
  }
  ###############################################################################
  t2 <- Sys.time()
  difftime(t2, t1, units = "mins")
  ###
  #GLMSubsetEvalStats <- GLMSubsetEvalL
  GLMSubsetEvalStats.df <- data.frame(do.call(rbind, GLMSubsetEvalStats), stringsAsFactors=FALSE)
  head(GLMSubsetEvalStats.df)
  tail(GLMSubsetEvalStats.df)
  nrow(GLMSubsetEvalStats.df)
  ncol(GLMSubsetEvalStats.df)
  colnames(GLMSubsetEvalStats.df) <- c("VarNames", "EnvVarsUsed", "SubsetVariableNumber", "DataType", "TSS", "Kappa", "AUC", "Spec", "Sens", "ThreshMxTSS", "AICc_bg", "AICc", "AICc2", "NumDVars", "NumEnVars")
  rownames(GLMSubsetEvalStats.df) <- c(seq(1:nrow(GLMSubsetEvalStats.df)))
  #str(GLMSubsetEvalStats.df1)
  # Save output
  # Convert second column and fourth through 11th columns from character to numeric
  GLMSubsetEvalStats.df[,c(3,5:15)] <- sapply(GLMSubsetEvalStats.df[,c(3,5:15)], function(x) as.numeric(as.character(x)))
  if(Output==TRUE) {
    setwd(OutDirectIn)
    Sets <- (nrow(GLMSubsetEvalStats.df)/2)
    if(DataSetType!="") {
      write.table(GLMSubsetEvalStats.df, file=paste0(Species, "GLMResults_WrapperTrainTest_", TotVars, "TotVars_", SubsetVariableNumber, "Vars_", DataSetType, "_", Sets, "_", SetRunID, ".csv"), sep=",", col.names=NA)
    } else {
      write.table(GLMSubsetEvalStats.df, file=paste0(Species, "GLMResults_WrapperTrainTest_", TotVars, "TotVars_", SubsetVariableNumber, "Vars_", Sets, "_", SetRunID, ".csv"), sep=",", col.names=NA)
    }
    #
    out5 <- paste("\n  GLM Response")
    cat(out5, file=paste(Species, "GLMWrapperResultsSummary_", TotVars, "TotVars_", ".txt", sep=""), sep="\n", append=TRUE)
  }
  #
  stopCluster(cl)
  registerDoSEQ()
  closeAllConnections()
  gc()
  # Delete temp directories created for individual GLM runs
  if(TempDir!="") { unlink(TempDir, recursive=TRUE) }
  #
  return(GLMSubsetEvalStats.df)
}
