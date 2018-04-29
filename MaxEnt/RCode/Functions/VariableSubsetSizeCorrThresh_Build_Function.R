VariableSubsetSizeCorrThresh.Build <-
function(CorrelationsIn, SearchLoopsIn, SubsetVariableNumber, MaxModSetsIn, CorrThreshIn, SetNameIn, OutDirectIn..., SpeciesIn, GridNamesIn) {
  ###############################################################################
  ## From a given set of variables, find sets of variables that have correlations
  ## below the maximum correlation threshold of CorrThresIn (usually 0.7) that
  ## are of a subset size of variables yielding the most sets close to MaxModSetsIn
  ###############################################################################
  #CorrelationsIn <- AbsValsCor.mat
  #CorrThreshIn <- 0.7
  #SetNameIn <- "Check"
  library(foreach)
  library(doParallel)
  #
  setwd(OutDirectIn)
  ##
  # Register cluster for running parallel with doParallel pkg
  # Use a minimum of one core or the number of available cores minus one
  cores <- detectCores() # detect cores running
  workcores <- max(1, cores-1)
  cl <- makeCluster(workcores)  # make cluster of workcores
  #cl <- makeCluster(1)
  registerDoParallel(cl) # register cluster
  getDoParWorkers() # check clusters registered
  #registerDoSEQ()
  ##
  ## Use parallel processing to find SearchLoopIn sets of SubsetvariableNumber variables below 0.7 correlation
  SetSize <- SubsetVariableNumber
  ## Find variables to be used in sets
  if(missing(SpeciesIn)) {SpeciesIn <- ""}
  #
  if(missing(GridNamesIn)) {
    VariableNames <- data.frame(colnames(CorrelationsIn), stringsAsFactors=FALSE)
  } else {
    VariableNames <- data.frame(GridNamesIn[,1, drop=FALSE])
    # Convert VariableNames from factor to character
    VariableNames <- apply(VariableNames,2, as.character)
  }
  ##
  ivec <- seq(1, SearchLoopsIn,1)
  SubsetBuildOutput <- foreach(p=ivec, .combine='cbind') %dopar% {
    # Randomize VariableNames
    VariableNamesR <- data.frame(VariableNames[sample(nrow(VariableNames)),], stringsAsFactors=FALSE)
    VarSet <- as.matrix(VariableNamesR[1,1])
    for(i in 2:nrow(VariableNamesR)) {
      #i=2
      CheckVar <- VariableNamesR[i,1]
      CorrCheck <- CorrelationsIn[VarSet, CheckVar]
      if(max(CorrCheck) < CorrThreshIn) {
        VarSet <- rbind(VarSet,CheckVar)
      }
      if(i==nrow(VariableNamesR) && nrow(VarSet) < SetSize) {
        for(j in (nrow(VarSet)+1):SetSize) {
        VarSet <- rbind(VarSet,NA)
        }
      }
      if(nrow(VarSet)==SetSize) {break}
    }
    rownames(VarSet) <- seq(1:nrow(VarSet))
    return(VarSet)
  }
  #################################################
  head(SubsetBuildOutput)
  #SetList[1:10]
  SetList.mat1 <- as.matrix(SubsetBuildOutput)
  #SetList.mat1[,1:10]
  SetList.mat2 <- SetList.mat1[ , apply(SetList.mat1, 2, function(x) !any(is.na(x)))]
  #
  if(ncol(SetList.mat2) < 1) {
    VarNamesSubsetsFinal <- 0
    } else {
    # Sort each column of matrix
    SetList.mat3 <- apply(SetList.mat2,2,sort,decreasing=F)
    #SetList.mat3[,1:10]
    # Only keep unique columns
    # Include matrix transpositions for proper working of unique function on rows
    SetList.mat <- t(unique(t(SetList.mat3)))
    ncol(SetList.mat)
    #
    #SetList.mat[,c(903,1208)]
    #SetList.mat2[,c(1,2)]
    # Format variable lists for subsets
    MaxSets <- min(ncol(SetList.mat), MaxModSetsIn)
    #
    if(ncol(SetList.mat)<2) {
      ivec <- 1
    } else {
      ivec <- c(seq(1, MaxSets, 1))
    }
    #
    ##################################
    t1 <- Sys.time()
    #VariablesSubsets <- VariableKeepSubsets
    VarNamesSubsetsFinal <- foreach(i=ivec, .combine='rbind') %dopar%  {
      VarNamesSubsets <- data.frame(paste0(unlist(SetList.mat[,i]), collapse="-"), stringsAsFactors=FALSE)
      return(VarNamesSubsets)
    }
    t2 <- Sys.time()
    ##################################
    #VarNamesSubsetsFinal[1:50,]
    colnames(VarNamesSubsetsFinal) <- "VarNames"
    head(VarNamesSubsetsFinal)
  }
  #
  stopCluster(cl)
  registerDoSEQ()
  closeAllConnections()
  gc()
  # Print summary output
  out1 <- paste("Found ", ncol(SetList.mat2), " subsets with ", SubsetVariableNumber, " variables from a total of")
  cat(out1,file=paste0(SetNameIn, "_CorrThresh", SetSize, "VarSets.txt"), sep="\n", append=TRUE)
  out2 <- paste(nrow(VariableNames), " variables with maximum correlation of ", CorrThreshIn, " in ", SearchLoopsIn, " loops")
  cat(out2, file=paste0(SetNameIn, "_CorrThresh", SetSize, "VarSets.txt"), sep="\n", append=TRUE)
  # Save ranked variable sets
  setwd(OutDirectIn)
  write.table(VarNamesSubsetsFinal, file=paste0(SpeciesIn, "Random_", MaxModSetsIn, "_", SetSize, "VariableSets.csv"), sep=",", col.names=NA)
  # Return ranked variable sets
  return(VarNamesSubsetsFinal)
}
