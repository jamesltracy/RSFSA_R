MaxentMultiSubset.FinalTrainTestEval <-
function(Species, VariableNamesIn, SubsetVariableNumber, TotPres, VariableSubsets, FinalTrainSWD, FinalTrainPresAbsID, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrpp, kfoldgrpa, OutDirectIn..., CatVarsPrefix, MaxentArgsIn, LongOutDirect, DataSetType, Output) {
  #(Species, VariableNamesIn, SubsetVariableNumber, TotPres, VariableSubsets, TrainSWD, TrainPresAbsID, TestSWD, TestPresAbsID, OutDirectIn, TestDataType)
  #CatVarsPrefix <- "ROADS_CAT"
  #NOTE: Function parameters after the "...", such as TestDataType, have to be set with an equal sign in the function call, such as ScoreType="mean"
  setwd(OutDirectIn)
  library(foreach)
  library(doParallel)
  # Make sure VariableNames is a data frame
  VariableNames <- data.frame(VariableNames, stringsAsFactors=FALSE)
  #
  RowNames.df <- data.frame(rownames(VariableSubsets), stringsAsFactors=FALSE)  # Save row names to use in output
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
  TotVars <- nrow(VariableNames)
  # Set default values for LongOutput, if optional values left out of function call
  if(missing(CatVarsPrefix)) { CatVarsPrefix="" }
  MaxentBaseArgsPar <- c("betamultiplier=1.0", "writebackgroundpredictions=false")
  if(missing(MaxentArgsIn)) {
    MaxentArgs1 <- MaxentBaseArgsPar
  } else {
    MaxentArgs1 <- MaxentArgsIn
  }
  if(CatVarsPrefix=="") {
    AnyCategoricalVars=FALSE
    } else {
    AnyCategoricalVars=TRUE
  }
  if(missing(LongOutDirect)) { LongOutDirect="None" }
  if(missing(DataSetType)) { DataSetType="" }
  #Output=FALSE
  if(missing(Output)) { Output=TRUE }
  #
  ModelType <- paste("Maxent", Subset, sep="")
  Model <- paste("Maxent for", SubsetVariableNumber, "Feature Subset")
  tail(VariableSubsets)
  #VariableSubsets[2998,]
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
  if(nrow(VariableSubsets)<2) {
    ivec <- 1
  } else {
    ivec <- c(seq(1, nrow(VariableSubsets), 1))
  }
  #ivec <- c(1,2,3)
  #tvec <- c(1,2,3,4)
  # Set up foreach to output to ten different list elements corresponding to
  # the two Old World and New World sets of five evaluation statistics
  ##############################################################
  t1 <- Sys.time()
  #VariablesSubsets <- VariableKeepSubsets
  MaxentSubsetEvalStats <- foreach(i=ivec, .combine='rbind') %dopar%  {
    library(dismo)
    library(raster)
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
    VariableNamesSel <- c(unlist(VariableSubsets[i,]))
    # Split VarNames of VariableSubsets into separate variables
    if(grepl("-",VariableNamesSel)) {
      VarNames <- unlist(strsplit(VariableNamesSel, "-"))
    } else {
      VarNames <- unlist(VariableNamesSel)
    }
    SubsetVarNum <- length(VarNames)
    ## Keep only data for selected variables
    head(FinalTrainSWD)
    MaxentTrainDataK <- FinalTrainSWD[,VarNames]
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
    MaxentOut <- maxent(MaxentTrainDataK, FinalTrainPresAbsID, args=MaxentArgs)
    #str(MaxentOut)
    # Copy Maxent output to desired directory
    #flist <- list.files(MaxentOut@path, full.names = TRUE)
    #file.copy(flist, OutDirectSub, overwrite=TRUE)
    ###
    ########### Evaluate model using training and test data
    MaxentSubsetEvalL <- list()
    TestDataTypes <- c("FinalTrain", "FinalTest")
    for(TestDataType in TestDataTypes) {
      #TestDataType="Train"
      if(TestDataType=="FinalTrain") {
        k <- 3
      } else {
        k <- 4
      }
      MaxentPresTestData <- PresenceDat.df[kfoldgrpp == k, ]
      head(MaxentPresTestData)
      nrow(MaxentPresTestData)
      MaxentAbsTestData <- PseudoabsenceDat.df[kfoldgrpa == k, ]
      head(MaxentAbsTestData)
      if(SubsetVarNum == 1) {
        MaxentPresTestData1 <- data.frame(MaxentPresTestData[,4])
        colnames(MaxentPresTestData1) <- colnames(MaxentPresTestData)[4]
        head(MaxentPresTestData1)
        MaxentAbsTestData1 <- data.frame(MaxentAbsTestData[,4])
        colnames(MaxentAbsTestData1) <- colnames(MaxentAbsTestData)[4]
        head(MaxentAbsTestData1)
        TestSWD <- rbind(MaxentPresTestData1, MaxentAbsTestData1)
        head(TestSWD)
      } else {
        TestSWD <- rbind(MaxentPresTestData[,4:ncol(MaxentPresTestData)], MaxentAbsTestData[,4:ncol(MaxentAbsTestData)])
        head(TestSWD)
      }
      TestPresID <- data.frame(rep(1,nrow(MaxentPresTestData)))
      colnames(TestPresID) <- "ID"
      TestAbsID <- data.frame(rep(0,nrow(MaxentAbsTestData)))
      colnames(TestAbsID) <- "ID"
      TestPresAbsID <- rbind(TestPresID, TestAbsID)
      if(SubsetVarNum == 1) {
        MaxentTestDataK <- TestSWD
      } else {
        MaxentTestDataK <- TestSWD[,VarNames]
        head(MaxentTestDataK)
        nrow(MaxentTestDataK)
      }
      ## Assign predictor values from above trained model to test data
      MaxentPred <- data.frame(predict(MaxentOut, MaxentTestDataK))
      colnames(MaxentPred) <- "MaxentScore"
      tail(MaxentPred)
      # Divide prediction data into presence and absence
      PresRows <- length(TestPresAbsID[which(TestPresAbsID$ID==1),])
      prespred.df <- data.frame(MaxentPred[1:PresRows,1])
      colnames(prespred.df) <- "MaxentScore"
      #
      AbsRows <- length(TestPresAbsID[which(TestPresAbsID$ID==0),])
      abspred.df <- data.frame(MaxentPred[(PresRows+1):(PresRows+AbsRows),1])
      colnames(abspred.df) <- "MaxentScore"
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
      # and assign the model results a name "EnvScoreN" where N is the name of the rep
      # Also assign a column "id" for the row numbers to use in merging the data frames later
      presa.df <- data.frame(c(rep(1, nrow(prespred.df))))
      names(prespred.df) <- c("EnvScore")
      names(presa.df) <- c("OBSERVED")
      pres.df <- data.frame(cbind(id=1:nrow(presa.df), presa.df, prespred.df))
      # Repeat above process with absence data, but assign "OBSERVED" a value of 0
      absa.df <- data.frame(c(rep(0, nrow(abspred.df))))
      names(abspred.df) <- c("EnvScore")
      names(absa.df) <- c("OBSERVED")
      abs.df <- data.frame(cbind(id=1:nrow(absa.df), absa.df, abspred.df))
      # For each model output, merge presence and absence data using "id' column as guide when all=TRUE
      # NOTE: PresenceAbsence package cannot handle several models at one time if the sample sizes differ
      # so have to analyze each model output separately
      presabspred <- merge(pres.df, abs.df, all=TRUE)
      # Drop the id column used in merging for each dataset
      presabspred$id <- NULL
      # Make a column of data with the species name with same number of rows as data from each model
      SPECIES <- data.frame(c(rep(Species, nrow(presabspred))))
      names(SPECIES) <- c("SPECIES")
      # Make final dataset SPDATA by putting together SPECIES with extracted environmental data.
      SPDATA <- data.frame(SPECIES, presabspred)
      ################################################################################
      ### Run this block of code to evaluate model results with PresenceAbsence package
      ################################################################################
      library(PresenceAbsence)
      #setwd(paste("C:/Users/JLTracy/Documents/R/win-library/3.0/10minClimIntEvalSTB1000", "/", output2, sep=""))
      #starttime <- Sys.time()
      #### FOR OLD WORLD DATA EVALUATION STATISTICS
      ### Define variables for later use.
      accurun <- list()
      accusum <- matrix(data=NA, ncol=6, nrow=1, byrow=TRUE, dimnames=list(NULL, c("MaxTSS", "Specificity_maxTSS", "Sensitivity_maxTSS", "AUC", "MaxKappa", "ThresholdMaxTSS")))
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
      #
      # To assess accuracy per threshold, use limited threshold available for
      # Envelope Score model based upon number of environmental layers in model
      # ("NumGrids")
      NumGrids <- max(40, nrow(VariableNames))
      PossThresholds <- seq(1/NumGrids,1,length=NumGrids)
      #accu <- presence.absence.accuracy(SPDATA, which.model = 1, threshold = PossThresholds, st.dev=FALSE)
      # accu <- presence.absence.accuracy(DATA, which.model = 1, threshold = c(0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.975, 0.98, 0.99, 0.999999))
      #accu <- presence.absence.accuracy(DATA, which.model = 1, threshold = 100, st.dev=FALSE)
      # print(paste("Species:", species[sp], "Model:", model.names))  not used
      accu <- data.frame(presence.absence.accuracy(DATA, which.model = 1, threshold = 100, st.dev=FALSE))
      # print(paste("Species:", species[sp], "Model:", model.names))  not used
      #accu
      maxSSS <-  data.frame(accu$sensitivity + accu$specificity)
      names(maxSSS) <- c("maxSSS")
      #maxSSS
      TSS <-  data.frame(accu$sensitivity + accu$specificity - 1)
      names(TSS) <- c("TSS")
      accurun <- data.frame(accu, maxSSS, TSS)
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
      accusum[1,6] <- max(ThresholdMaxTSS, 0.0001)
      ###############################
      #endtime <- Sys.time()
      #durtime <- endtime - starttime
      # Save evaluation statistics to .csv file
      #setwd(paste("C:/Users/JLTracy/Documents/R/win-library/3.0/10minClimIntEvalSTB1000", "/", output2, sep=""))
      accusum.df <- data.frame(accusum)
      accusum.df$Model <- ModelType
      #
      # Save variable matrix coordinates and Old World evaluation statistics to a vector
      MaxentSubsetEvalL[[k]] <- c(VariableNamesSel, SubsetVariableNumber, TestDataType, accusum.df$MaxTSS, accusum.df$MaxKappa, accusum.df$AUC, accusum.df$Specificity_maxTSS, accusum.df$Sensitivity_maxTSS, accusum.df$ThresholdMaxTSS)
      #
      if(LongOutDirect!="None") {
        setwd(paste(LongOutDirect))
        # Save threshold verus evaluation statitistic data for run
        write.table(SPDATA, file=paste(Species, VariableNamesSel, "SPDATA.csv", sep=""), sep=",")
        write.table(accurun, file=paste(Species, "EnvScrResults", VariableNamesSel, "StatsVsThresholds_Run", ".csv", sep=""), sep=",", col.names=NA)
        print(paste("Species:", species[sp], "Model:", model.names))
        print(accurun)
        bmp(paste(Species, VariableNamesSel,"TSSvsThreshold", ".bmp", sep=""))
        plot(accurun$threshold, accurun$TSS, type="l")
        dev.off()
        # Plot a ROC plot
        # save following plot of grid as .bmp file
        bmp(paste(Species, VariableNamesSel, "EnvScoreROC", ".bmp", sep=""))
        auc.roc.plot(DATA, color = TRUE, legend.cex = 1.2, main = "")
        dev.off()
        # Calculate optimal thresholds by various methods and save output to text file
        outthresholds <- capture.output(optimal.thresholds(DATA, opt.methods = 1:12, threshold=10001))
        out1 <- paste(outthresholds)
        cat(out1, file=paste(Species, "EnvScrResults", VariableNamesSel, "Thresholds.txt", sep=""), sep="\n", append=TRUE)
        # Can calculate the confusion matrix for a given threshold (not necessary)
        confmatrix <- cmx(DATA, threshold = ThresholdMaxTSS, na.rm = FALSE)
        # Save confusion matrix at threshold of maxTSS
        write.table(confmatrix, file=paste(Species, "EnvScrResults", VariableNamesSel, "ConfusionMatrixatMaxTSSThreshold_Run", ".csv", sep=""), sep=",", col.names=c("Obs Pres","Obs Abs"), row.names=c("Pred Pres", "Pred Abs"))
        #Save evaluation stastics to file
        write.table(accusum.df, file=paste(Species, "EnvScrResults", VariableNamesSel, "Stats", ".csv", sep=""), sep=",", col.names=NA)
        graphics.off()
        # Save output evaluation statistics and summary statistics to text file
        # Output codes of variables used
        outA <- paste("\nDate/Time:", Sys.time(), "\nProgram Start Time:", starttime, "\nProgram End Time:", endtime, "\nProgram Execution Time:", durtime, "\n\nSpecies:", Species, "\nVariable Subset:", VariableNamesSel, "\nModel:", Model, "\nScoring Alogorithm: Maxent")
        cat(outA, file=paste(Species, "EnvScrResults", VariableNamesSel, "Summary.txt", sep=""), sep="\n", append=TRUE)
        #NotUsedVars <- capture.output(print(DeleteBands))
        #outB <- paste("\nVariables Not Used: Code-", NotUsedVars, "; Variable-", GridNamesDrop)
        #cat(outB, file=paste(Species, InfClass, "EnvScrResults", output2, "Summary.txt", sep=""), sep="\n", append=TRUE)
        #UsedVars <- capture.output(print(SRandNumbers))
        #outB <- paste("\n  Code Nos. for Variables Used:", UsedVars)
        #cat(outB, file=paste("STBEnvScrResults", output2, "Summary.txt", sep=""), sep="\n", append=TRUE)
        #out1 <- paste("\nVariable Names")
        #cat(out1, file=paste(Species, "EnvScrResults", GridNamesKeepSelT, "Summary.txt", sep=""), sep="\n", append=TRUE)
        #EnvVariables <- as.matrix(GridNamesKeepSel)
        #out2 <- capture.output(print(EnvVariables))
        #cat(out2, file=paste(Species, "EnvScrResults", GridNamesKeepSelT, "Summary.txt", sep=""), sep="\n", append=TRUE)
        out3 <- paste("\n  Evaluation Statistics per Run")
        cat(out3, file=paste(Species, "EnvScrResults", VariableNamesSel, "Summary.txt", sep=""), sep="\n", append=TRUE)
        out4 <- capture.output(print(accusum.df))
        cat(out4, file=paste(Species, "EnvScrResults", VariableNamesSel, "Summary.txt", sep=""), sep="\n", append=TRUE)
        }
    }
    return(MaxentSubsetEvalL)
  }
  ###############################################################################
  t2 <- Sys.time()
  durtime <- difftime(t1, t2, units = "mins")
  ###
  MaxentSubsetEvalStats.df <- data.frame(do.call(rbind, MaxentSubsetEvalStats), stringsAsFactors=FALSE)
  head(MaxentSubsetEvalStats.df)
  tail(MaxentSubsetEvalStats.df)
  nrow(MaxentSubsetEvalStats.df)
colnames(MaxentSubsetEvalStats.df) <- c("VarNames", "SubsetVariableNumber", "DataType", "TSS", "Kappa", "AUC", "Spec", "Sens", "ThreshMxTSS")
  rownames(MaxentSubsetEvalStats.df) <- c(seq(1:nrow(MaxentSubsetEvalStats.df)))
  #str(MaxentSubsetEvalStats.df1)
  # Save output
  # Convert second column and fourth through ninth columns from character to numeric
  MaxentSubsetEvalStats.df[,c(2,4:9)] <- sapply(MaxentSubsetEvalStats.df[,c(2,4:9)], function(x) as.numeric(as.character(x)))
  if(Output==TRUE) {
    setwd(OutDirectIn)
    Sets <- (nrow(MaxentSubsetEvalStats.df)/2)
    if(DataSetType!="") {
      write.table(MaxentSubsetEvalStats.df, file=paste0(Species, "MaxentResults_TrainTest_", TotVars, "TotVars_", SubsetVariableNumber, "Vars_", DataSetType, "_", Sets, ".csv"), sep=",", col.names=NA)
    } else {
      write.table(MaxentSubsetEvalStats.df, file=paste0(Species, "MaxentResults_TrainTest_", TotVars, "TotVars_", SubsetVariableNumber, "Vars_", Sets, ".csv"), sep=",", col.names=NA)
    }
    out5 <- capture.output(print(MaxentArgsIn))
    cat(out5, file=paste(Species, "MaxentTrainTestResultsSummary_", TotVars, "TotVars_", ".txt", sep=""), sep="\n", append=TRUE)
  }
  #
  stopCluster(cl)
  registerDoSEQ()
  closeAllConnections()
  gc()
  # Delete temp directories created for individual maxent runs
  if(TempDir!="") { unlink(TempDir, recursive=TRUE) }
  #
  return(MaxentSubsetEvalStats.df)
}
