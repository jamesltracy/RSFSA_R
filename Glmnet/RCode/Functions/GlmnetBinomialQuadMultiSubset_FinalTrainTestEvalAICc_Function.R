GlmnetBinomialQuadMultiSubset.FinalTrainTestEvalAICc <-
function(Species, VariableNamesIn, SubsetVariableNumber, TotPres, VariableSubsetsIn, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrpp, kfoldgrpa, cvglmnetfoldsin, alphain, OutDirectIn, FunctDirectIn..., SetRunIDIn, CVGlmnetRuns, Output, DataSetType) {
  setwd(OutDirectIn)
  library(foreach)
  library(doParallel)
  #FunctDirectIn <- "C:/Users/James/Documents/R/win-library/"
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
  # Set default value for CVGlmnetRuns if not entered (number runs of cv.glmnet)
  if(missing(DataSetType)) { DataSetType="" }
  if(missing(CVGlmnetRuns)) { CVGlmnetRuns=1 }
  #Output=FALSE
  if(missing(Output)) { Output=TRUE }
  if(missing(SetRunIDIn)) { SetRunIDIn="" }
  #
  ModelType1 <- paste("GlmnetBinomialQuad", Subset, sep="")
  Model <- paste("GlmnetBinomialQuad for", SubsetVariableNumber, "Feature Subset")
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
  ## Designate final training data and add count column as factor
  GlmnetPresTrainData <- PresenceDat.df[kfoldgrpp==3,]
  GlmnetPresTrainData$Count <- 1
  head(GlmnetPresTrainData)
  #str(GlmnetPresTrainData)
  GlmnetAbsTrainData <- PseudoabsenceDat.df[kfoldgrpp==3,]
  GlmnetAbsTrainData$Count <- 0
  ##############################################################
  t1 <- Sys.time()
  #VariablesSubsets <- VariableKeepSubsets
  GlmnetSubsetEvalStats <- foreach(i=ivec, .combine='rbind') %dopar%  {
    library(dismo)
    library(raster)
    library(glmnet)
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
    ## Define function to predict matrix with Glmnet
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
    ###
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
    # Remove Variables with coefficient of zero
    AvgCoefficients.df1 <- AvgCoefficients.df[!AvgCoefficients.df$Coefficients==0,]
    # Remove VarNames column
    AvgCoefficients.df <- AvgCoefficients.df1[, 2, drop=FALSE]
    #
    AICc2 <- mean(ModelStatsResults.df$AICc2)
    ## If there are no coefficients return NA for this iteration
#    GlmnetSubsetEvalL <- list()
#    if(nrow(AvgCoefficients.df)==1) {
#      k <- 2
#      GlmnetSubsetEvalL[[k]] <- c(rep(NA,14))
#      return(GlmnetSubsetEvalL)
#    }
    #####################################################
    #Save coefficients
    #setwd(OutDirectIn)
    #write.table(AvgCoefficients.df, file=paste("Coefficients", ModelType1, "_", "TrainSetEnsAvg", KeepModels, "of", CVGlmnetRuns, "_", SubsetVariableNumber, "Var_", Loop, ".csv", sep=""), sep=",", col.names=NA)
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
    #####################################
    ## Calculate AICc_bg with point values
    # Obtain number of parameters in Glmnet model
    nparams <- (nrow(AvgCoefficients.df) - 1)
    #
    ## From ENMeval Package documentation: AICc is the Akaike Information Criterion corrected for small
    ## sample sizes calculated as: (2 * K - 2 * logLikelihood) + (2 * K) * (K + 1)=(n - K - 1)
    ## where K is the number of parameters in the model (i.e., number of non-zero parameters in Glmnet
    ## lambda file) and n is the number of occurrence localities.
    ## The logLikelihood is sum(log(vals/total))
    ## vals is vector of Glmnet raw values at occurence localities
    ## total is the sum of Glmnet raw values across the entire study area
    ##
    head(PresAbsTrainDat.df)
    nrow(PresAbsTrainDat.df)
    ######
    xPresTrainData.df1 <- PresAbsTrainDat.df[1:nrow(GlmnetPresTrainData),]
    nrow(xPresTrainData.df1)
    xPresTrainData.df <- xPresTrainData.df1[, c(VarsUsed), drop=FALSE]
    head(xPresTrainData.df)
    #
    ## Add raw quadratic versions for main effect variables (no interaction)
    xPresTrainDataQ.df <- cbind(xPresTrainData.df, xPresTrainData.df)
    colcount <- 1
    for(i in 1:ncol(xPresTrainData.df)) {
      #i=3
        # Use poly function to obtained orthonalized main effect and squared quadratic value
        xPresTrainDataQ.df[,colcount:(colcount+1)] <- poly(xPresTrainData.df[,i], 2, raw=TRUE)
        head(xPresTrainDataQ.df)
        colnames(xPresTrainDataQ.df)[colcount:(colcount+1)] <- c(colnames(xPresTrainData.df)[i], paste0(colnames(xPresTrainData.df)[i], "^2"))
        colcount <- colcount + 2
    }
    head(xPresTrainDataQ.df)
    # Select only terms used in ensemble model
    xPresData.dfIn <- xPresTrainDataQ.df[, c(VarNamesUsed), drop=FALSE]
    head(xPresData.dfIn)
    # Obtain values for training presence data output by glmnet run
    CoefficientsIn <- AvgCoefficients.df
    presvals.df <- GLMBinomialResp_PredictMatrix(CoefficientsIn, xPresData.dfIn)
    #str(xTrainData)
    # Obtain values for training presence/absence data output by glmnet run
    BackgroundDat.df2 <- BackgroundDat.df[,c(VarsUsed), drop=FALSE]
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
    #####################
    # Select only terms used in ensemble model
    xData.dfIn <- BackgroundDatQ.df2[, c(VarNamesUsed), drop=FALSE]
    head(xData.dfIn)
    # Obtain values for training presence/absence data output by glmnet run
    CoefficientsIn <- AvgCoefficients.df
    backgroundvals.df <- GLMBinomialResp_PredictMatrix(CoefficientsIn, xData.dfIn)
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
    ########### Evaluate model using training and test data
    GlmnetSubsetEvalL <- list()
    TestDataTypes <- c("FinalTrain", "FinalTest")
    for(TestDataType in TestDataTypes) {
      #TestDataType="FinalTrain"
      if(TestDataType=="FinalTrain") {
        k <- 3
      } else {
        k <- 4
      }
      # Specify model testing data
      GlmnetPresTestData <- PresenceDat.df[kfoldgrpp == k, ]
      head(GlmnetPresTestData)
      nrow(GlmnetPresTestData)
      ###
      GlmnetAbsTestData <- PseudoabsenceDat.df[kfoldgrpa == k, ]
      ########################
      ### Predict model values for test matrix to get values for test presence and absence points
      ## First presence points
      setwd(OutDirectIn)
      head(GlmnetPresTestData)
      # Subset PresTestData by VarNames
      PresTestData.df <- GlmnetPresTestData[, c(VarsUsed), drop=FALSE]
      head(PresTestData.df)
      #
      ## Add raw quadratic versions for main effect variables (no interaction)
      PresTestDataQ.df <- cbind(PresTestData.df, PresTestData.df)
      colcount <- 1
      for(i in 1:ncol(PresTestData.df)) {
        #i=3
          # Use poly function to obtained orthonalized main effect and squared quadratic value
          PresTestDataQ.df[,colcount:(colcount+1)] <- poly(PresTestData.df[,i], 2, raw=TRUE)
          head(PresTestDataQ.df)
          colnames(PresTestDataQ.df)[colcount:(colcount+1)] <- c(colnames(PresTestData.df)[i], paste0(colnames(PresTestData.df)[i], "^2"))
          colcount <- colcount + 2
      }
      head(PresTestDataQ.df)
      # Select only terms used in ensemble model
      xData.dfIn <- PresTestDataQ.df[, c(VarNamesUsed), drop=FALSE]
      head(xData.dfIn)
      CoefficientsIn <- AvgCoefficients.df
      prespred.df <- GLMBinomialResp_PredictMatrix(CoefficientsIn, xData.dfIn)
      colnames(prespred.df) <- "GlmnetScore"
      ## Then absence points
      head(GlmnetAbsTestData)
      # Subset PresTestData by VarNames
      AbsTestData.df <- GlmnetAbsTestData[, c(VarsUsed), drop=FALSE]
      head(AbsTestData.df)
      #
      ## Add raw quadratic versions for main effect variables (no interaction)
      AbsTestDataQ.df <- cbind(AbsTestData.df, AbsTestData.df)
      colcount <- 1
      for(i in 1:ncol(AbsTestData.df)) {
        #i=3
          # Use poly function to obtained orthonalized main effect and squared quadratic value
          AbsTestDataQ.df[,colcount:(colcount+1)] <- poly(AbsTestData.df[,i], 2, raw=TRUE)
          head(AbsTestDataQ.df)
          colnames(AbsTestDataQ.df)[colcount:(colcount+1)] <- c(colnames(AbsTestData.df)[i], paste0(colnames(AbsTestData.df)[i], "^2"))
          colcount <- colcount + 2
      }
      head(AbsTestDataQ.df)
      # Select only terms used in ensemble model
      xData.dfIn <- AbsTestDataQ.df[, c(VarNamesUsed), drop=FALSE]
      head(xData.dfIn)
      #
      abspred.df <- GLMBinomialResp_PredictMatrix(CoefficientsIn, xData.dfIn)
      colnames(abspred.df) <- "GlmnetScore"
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
      accusum <- matrix(data=NA, ncol=11, nrow=1, byrow=TRUE, dimnames=list(NULL, c("MaxTSS", "Specificity_maxTSS", "Sensitivity_maxTSS", "AUC", "MaxKappa", "ThresholdMaxTSS", "AICc_bg", "AICc2", "NumDVars", "NumEnVars", "EnvVarsUsed")))
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
      accusum[1,8] <- AICc2
      accusum[1,9] <- NumDVars
      accusum[1,10] <- nenvars
      accusum[1,11] <- EnvVarsUsed
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
      GlmnetSubsetEvalL[[k]] <- c(VariableNamesSel, accusum.df$EnvVarsUsed, SubsetVariableNumber, TestDataType, accusum.df$MaxTSS, accusum.df$MaxKappa, accusum.df$AUC, accusum.df$Specificity_maxTSS, accusum.df$Sensitivity_maxTSS, accusum.df$ThresholdMaxTSS, accusum.df$AICc_bg, accusum.df$AICc2, accusum.df$NumDVars, accusum.df$NumEnVars)
      #
    }
    return(GlmnetSubsetEvalL)
  }
  ###############################################################################
  t2 <- Sys.time()
  difftime(t2, t1, units = "mins")
  ###
  #GlmnetSubsetEvalStats <- GlmnetSubsetEvalL
  GlmnetSubsetEvalStats.df <- data.frame(do.call(rbind, GlmnetSubsetEvalStats), stringsAsFactors=FALSE)
  head(GlmnetSubsetEvalStats.df)
  tail(GlmnetSubsetEvalStats.df)
  nrow(GlmnetSubsetEvalStats.df)
  ncol(GlmnetSubsetEvalStats.df)
  colnames(GlmnetSubsetEvalStats.df) <- c("VarNames", "EnvVarsUsed", "SubsetVariableNumber", "DataType", "TSS", "Kappa", "AUC", "Spec", "Sens", "ThreshMxTSS", "AICc_bg", "AICc2", "NumDVars", "NumEnVars")
#  # Omit any rows with NA values
#  GlmnetSubsetEvalStats.df <- na.omit(GlmnetSubsetEvalStats.df)
#  any(is.na(GlmnetSubsetEvalStats.df))
#  # Keep only NumberModSets rows
#  GlmnetSubsetEvalStats.df <- GlmnetSubsetEvalStats.df[1:NumberModSets,]
  rownames(GlmnetSubsetEvalStats.df) <- c(seq(1:nrow(GlmnetSubsetEvalStats.df)))
  #str(GlmnetSubsetEvalStats.df1)
  # Save output
  # Convert second column and fourth through 11th columns from character to numeric
  GlmnetSubsetEvalStats.df[,c(3,5:14)] <- sapply(GlmnetSubsetEvalStats.df[,c(3,5:14)], function(x) as.numeric(as.character(x)))
  if(Output==TRUE) {
    setwd(OutDirectIn)
    Sets <- nrow(GlmnetSubsetEvalStats.df)/2
    if(DataSetType!="") {
      if(SetRunIDIn!="") {
        write.table(GlmnetSubsetEvalStats.df, file=paste0(Species, "GlmnetResults_TrainTest_", TotVars, "TotVars_", SubsetVariableNumber, "Vars_", DataSetType, "_", Sets, "_", SetRunIDIn, ".csv"), sep=",", col.names=NA)
      } else {
        write.table(GlmnetSubsetEvalStats.df, file=paste0(Species, "GlmnetResults_TrainTest_", TotVars, "TotVars_", SubsetVariableNumber, "Vars_", DataSetType, "_", Sets, ".csv"), sep=",", col.names=NA)
      }
    } else {
      if(SetRunIDIn!="") {
        write.table(GlmnetSubsetEvalStats.df, file=paste0(Species, "GlmnetResults_FinalTrainTest_", TotVars, "TotVars_", SubsetVariableNumber, "Vars_", Sets, "_", SetRunIDIn, ".csv"), sep=",", col.names=NA)
      } else {
        write.table(GlmnetSubsetEvalStats.df, file=paste0(Species, "GlmnetResults_FinalTrainTest_", TotVars, "TotVars_", SubsetVariableNumber, "Vars_", Sets, ".csv"), sep=",", col.names=NA)
      }
    }
    out5 <- paste("\n  cvglmnetfolds: ", cvglmnetfoldsin, "    alpha: ", alphain, "    CVGlmnetRuns: ", CVGlmnetRuns)
    cat(out5, file=paste(Species, "GlmnetWrapperResultsSummary_", TotVars, "TotVars_", ".txt", sep=""), sep="\n", append=TRUE)
  }
  #
  stopCluster(cl)
  registerDoSEQ()
  closeAllConnections()
  gc()
  # Delete temp directories created for individual Glmnet runs
  if(TempDir!="") { unlink(TempDir, recursive=TRUE) }
  #
  return(GlmnetSubsetEvalStats.df)
}
