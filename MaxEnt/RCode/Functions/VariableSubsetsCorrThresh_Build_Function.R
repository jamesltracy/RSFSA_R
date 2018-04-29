VariableSubsetsCorrThresh.Build <- 
function(CorrelationsIn, SearchLoopsIn, MaxModSetsIn, CorrThreshIn, SetNameIn, OutDirectIn) {
  ###############################################################################
  ## From a given set of variables, find sets of variables that have correlations
  ## below the maximum correlation threshold of CorrThresIn (usually 0.7) that 
  ## are of a subset size of variables yielding the most sets close to MaxModSetsIn
  ###############################################################################
  #CorrelationsIn <- AbsValsCor.mat
  #StartNumVarsIn <- 10
  #MinModSetsIn <- 12
  #MaxModSetsIn <- 13
  #CorrThreshIn <- 0.7
  #SetNameIn <- "Check"
  library(caret)
  setwd(OutDirectIn)
  ##################################################################################
  ### Run a loop checking CorrelationsIn between variables starting at either end of
  ### randomly ordered variables. Loop will iterate to find the lowest possible
  ### CorrelationsIn for the given number of variables, NumVars
  ################################################################################
  ## Set number of loop (loops)
  TotVars <- nrow(CorrelationsIn)
  #
  ## Use function to round up from .5 from http://stackoverflow.com/questions/12688717/round-up-from-5-in-r
  Round2 <- function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
  }
  ## Use findCorrelation function of caret package to estimate maximum number of variables
  ## in a subset that can meet criterion of having correlation less than CorrThreshIn
  ## Add to the number of variables in case of an underestimate
  VariableNamesCut <- findCorrelation(AbsValsCor.mat, cutoff=CorrThreshIn, names=TRUE, exact=TRUE)
  VariableNamesKeep <- data.frame(VariableNamesSelect[!(VariableNamesSelect[,1]) %in% VariableNamesCut,], stringsAsFactors=FALSE)
  colnames(VariableNamesKeep) <- "VarNames"
  StartNumVars <- Round2((nrow(VariableNamesKeep) + 0.4*(nrow(VariableNamesKeep))),0)
  #
  ## Find total possible combinations of NumVars variables in a set of TotVars variables
  #format(choose(TotVars,NumVars),scientific=FALSE)
  #
  ### Run loop to output sets of given number of variables with lowest correlation among them
  ### The search starts with a single random variable, and the correlation of each succeeding
  ### randomly chosen variable is compared with the previous to make sure it is under the
  ### threshold "CorrThreshIn" until the desired subset number of variables (NumvVars) is reached
  #
  ### To optimize the speed of finding variable sets with the minimum maximum correlation,
  ### Will need to adjust down the correlation threshold "CorrThreshIn" until it is only about
  ### 0.01 above the minimum of the maximum thresholds observed in about 10000 runs (SearchLoops)
  #
  NumVarsInSet <- seq(StartNumVars, 2, -1)
  ###
  t1 <- Sys.time()
  Count <- 0
  PrevSet <- 0
  SelectvarsL <- list()
  for(h in 1:length(NumVarsInSet)) {
    #h=4
    NumVarsIn <- NumVarsInSet[h]
    ## Define and initialize variables for loop finding sets of variables at lowest correlation.
    RowNum <- NumVarsIn + 2
    # Matrix mselectvar will store values of set of n="NumVars" variables from all loops
    # and store them with the maximum correlation among them
    mselectvar <- matrix(data=NA, ncol=SearchLoopsIn, nrow=RowNum, byrow=TRUE)
    # Matrix mcorrvals is used for checking accuracy of stored CorrelationsIn between individual
    # sets of variables for each i loop below
    mcorrvals <- matrix(data=NA, ncol=SearchLoopsIn, nrow=RowNum, byrow=TRUE)
    # Define progres text bar for loops
    #pb <- txtProgressBar(1, SearchLoopsIn, style=3)
    #
    #########################
    for (n in 1:SearchLoopsIn) {
    #PID <- Sys.getpid()
    #setTxtProgressBar(pb, n)
    # Matrix selectvars stores a set of variables and their maximum correlation from
    # each large loop (n)
    selectvars <- matrix(data="none", ncol=1, nrow=NumVarsIn, byrow=TRUE, dimnames=list(NULL,
      c(paste(NumVarsIn, "selected variables with CorrelationsIn <", CorrThreshIn, "using randomly ordered variables at replicate", n))))
    # Matrix "corrvaluestore" stores the maximum of CorrelationsIn among variables
    # found in the loop below
    corrvaluestore <- matrix(data=-1, ncol=1, nrow=NumVarsIn, byrow=TRUE)
    # "variables" selects variables V1 through the total (TotVars)
    variables <- CorrelationsIn[,1]
    # "variablesr" randomizes the order of variables
    variablesr <- variables[order(runif(length(variables)))]
    # The matrix "corrvalues" stores correlation values between one variable and the
    # preceding randomly chose variables in the k loop below
    corrvalues <- matrix(data=-1, ncol=1, nrow=NumVarsIn-1, byrow=TRUE, dimnames=list(NULL,
     c("Correlation Values")))
    # "IniLook is the first randomly chosen variable of a set
    IniLook <- names(variablesr[1])
    # "IniLook" is stored as the first variable of a set in the matrix "selectvars"
    selectvars[1,] <- IniLook
    corrvaluestore[1,] <- -1
    t=1
    # The following loops process through and choose the first set of variables
    # trying to find a set with maximum correlation less than the threshold (CorrThreshIn)
    # and if the effort fails, a value of "failed" is where the last variable would be
    # in the matrix "selectvars"
    for (i in 2:NumVarsIn) {
      if(selectvars[NumVarsIn,]=="failed") break
        for (j in 2:TotVars) {
         VarLook <- names(variablesr[j])
         if(i==2) { corrvalues[1] <- CorrelationsIn[VarLook, IniLook] }
         else {
           for (k in 1:t) {
             corrvalues[k] <- CorrelationsIn[VarLook, selectvars[k,]]
           }
          }
          maxcorr <- max(corrvalues)
          if (maxcorr < CorrThreshIn) {
             selectvars[i,] <- c(VarLook)
             corrvaluestore[i,] <- maxcorr
             t=t+1
          break }
          if(j==TotVars & maxcorr >= CorrThreshIn) {
            selectvars[NumVarsIn,]="failed"
            mselectvar[1,n] <- n
            mselectvar[2,n] <- maxcorr
            mselectvar[3:RowNum,n] <- selectvars
            mcorrvals[1,n] <- n
            mcorrvals[2,n] <- VarLook
            mcorrvals[3:RowNum,n] <- corrvaluestore
           break }
        }
    }
      if(selectvars[NumVarsIn,]!="none") {
         mselectvar[1,n] <- n
         mcorrvals[1,n] <- n
         mcorrvals[2,n] <- VarLook
         mselectvar[3:RowNum,n] <- selectvars
         mcorrvals[3:RowNum,n] <- corrvaluestore
         if(selectvars[NumVarsIn,]=="failed") {
           mselectvar[2,n] <- maxcorr
           # For checking computations at certain conditions
    #         if(maxcorr==1) {
    #           checkvars <- variablesr
    #           loopcheck <- n
    #           varlookcheck <- VarLook
    #           maxcorrcheck <- maxcorr
    #         }
         }
         else { mselectvar[2,n] <- max(corrvaluestore) }
      }
    }
    #
    # For checking computations at certain loop iterations
    #mselectvar[,1:10]
    #mcorrvals[,831]
    ncol(mselectvar)
    #
    # Delete unfilled matrix columns if program interrupted
    #mselectvar[,1:50]
    #mselectvar <- mselectvar[,-c(17489:500000)]
    # If program interrupted, reset loops to actual columns in matrix mselectvar
    #SearchLoopsIn <- 17488
    #
    #mselectvar <- as.matrix(read.csv("mselectvarSep3aPID4168loop10000.csv", header = TRUE, sep=',', row.names=1))
    #
    # Delete columns of variable sets for which correlation threshold not met
    delcol <- which(mselectvar[RowNum,]=="failed")
    if(length(delcol)>0) {
      mselectvar2 <- as.matrix(mselectvar[,-delcol])
    } else {
      mselectvar2 <- mselectvar
    }
    #mselectvar2[,1:50]
    origcols <- ncol(mselectvar2)
    #
    # Find minimum correlation among selected sets of variables
    if(origcols==0) {
      mincorr <- 0
      maxcorr <- 0
    } else {
      mincorr <- min(as.numeric(mselectvar2[2,]))
      maxcorr <- max(as.numeric(mselectvar2[2,]))
    }  
    #
    #
    mselectvarmaxb <- as.matrix(mselectvar2[-1:-2,])
    #mselectvarmaxb[,1:50]
    # Sort columns
    mselectvarmaxbs <- apply(mselectvarmaxb, 2, sort)
    #mselectvarmaxbs[,1:50]
    # Calculate percentage of duplicate variable sets in matrices
    # Include matrix transpositions for proper working of unique function on rows
    mselectvarmaxsu <- t(unique(t(mselectvarmaxbs)))
    fincols <- ncol(mselectvarmaxsu)
    percduplicates <- (1 - fincols/origcols) * 100
    #mselectvarmaxsu[,1:50]
    ncol(mselectvarmaxsu)
    #
    # Print summary output
    out1 <- paste(fincols, "subsets found with ", NumVarsIn, "variables from a total of")
    cat(out1,file=paste0(SetNameIn, "CorrThreshVarSets.txt"), sep="\n", append=TRUE)
    out2 <- paste(TotVars, "variables with maximum correlation of", maxcorr, "in", SearchLoopsIn, " loops")
    cat(out2, file=paste0(SetNameIn, "CorrThreshVarSets.txt"), sep="\n", append=TRUE)
    #print(" ")
    #print(paste0("For ", SetNameIn, ":"))
    #print(out1)
    #print(out2)
    #
    Count <- Count + 1
    SelectvarsL[[Count]] <- mselectvarmaxsu
    if(ncol(mselectvarmaxsu) >= MaxModSetsIn) { break }
    # if set is smaller than the previous one then break
    if(Count > 1) {
      PrevSet <- 1
      if(length(SelectvarsL[[Count-1]])> length(SelectvarsL[[Count]])) { break }
    }
  }
  t2 <- Sys.time()
  ###############################################################################
  difftime(t2,t1, units = "mins")
  #
  # Use data from previous iteration of loop, if has larger number of sets
  if(PrevSet==1) {
    mselectvarmaxsu = SelectvarsL[[h-1]]
  }
  # Transform each set of variables into a VariableSubset
  # Keep only first up to MaxModSetsIn sets
  VarNamesSubsetsL <- list()
  if(ncol(mselectvarmaxsu)>0) {
    for(j in 1:ncol(mselectvarmaxsu)) 
    #j=1
    VarNamesSubsetsL[j] <- as.data.frame(paste0(unlist(mselectvarmaxsu[,j]), collapse="-"), stringsAsFactors=FALSE)
    VarNamesSubsets <- data.frame(do.call(rbind, VarNamesSubsetsL), stringsAsFactors=FALSE)
    head(VarNamesSubsets)
    tail(VarNamesSubsets)
    nrow(VarNamesSubsets)
    colnames(VarNamesSubsets) <- "VarNames"
    if(fincols>MaxModSetsIn) {
      Max <- MaxModSetsIn
    } else {
      Max <- fincols
    }
    VarNamesSubsetsF <- as.data.frame(VarNamesSubsets[1:Max,], stringsAsFactors=FALSE)
    tail(VarNamesSubsetsF)
    colnames(VarNamesSubsetsF) <- "VarNames"
  } else {
    VarNamesSubsetsF <- 0
  }    
  # Save output of variable sets at given minimum correlation (duplicates deleted)
  #write.csv(mselectvarmaxsu, file=paste("mselectvarmax", output, "loop", n, ".csv", sep=""))
  return(VarNamesSubsetsF)
}  
  
