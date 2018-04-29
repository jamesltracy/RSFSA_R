##############################################################################
# This R program implements a Randoms Subset Feature Selection Algorithm for
# MaxEnt niche modelling described in Tracy et al. (in prep.)
#
# Tracy JL, Trabucco A, Lawing AM, Giermakowski T, Tchakerian M, Drus GM, Coulson RN
#   (in prep.) Random subset feature selection of ecological niche models for wildfire
#    activity in western North America
#
# The program specifically processes occurrence and absence data or large
# wildfires of low burn severity for the western United States and projects
# burned areas for all of western North America current, 2050 and 2070 climates
##############################################################################
#
##################################################################################
# NOTE: Some sections are commented since pre-processed presence, pseudoabsence,
# and background data is provided in the FireVignetteData folder without having
# to query raster environmental layers
# In addition, several k-fold partition schemes are commented since they are provided
# in FireVignetteData folder
# Modify the three input directories to match where you put input files on your system
# 1) For R Functions: C:/Users/James/Documents/R/win-library/
# 2) For names of variables: E:/30SecClimIntAnthTopoNorm
# 3) For csv files with extracted environmental data and k-fold partition schemes
#    E:/FireFrequency
# Also, various output files are created under E:/FireFrequency
# Code after line 1546 requires 90 North American environmental layers at 1 km resolution
# that are not provided for current and future climates
#################################################################################
#
##############################################################################
# This section loads libraries, sets working directory, and establishes source
# file locations and output file naming conventions
##############################################################################
# Specify  directory for saved functions
FunctDirect <- "C:/Users/James/Documents/R/win-library/"
setwd(FunctDirect)
# Read in User Defined Functions
source("SpatialFilter_Function.R")
source("Round2_Function.R")
source("MaxentSubset_TrainTestEval_Function.R")
source("MaxentSubset_GridTrainTestEvalAIC_AUCbgp_Calib_Function.R")
source("MaxentSubset_TrainTestEvalK_Function.R")
source("MaxentMultiSubset_WrapperTrainEvalAIC_Function.R")
source("MaxentMultiSubset_WrapperTestEvalAIC_Function.R")
source("MaxentMultiSubset_WrapperTrainTestEvalAIC_Function.R")
source("MaxentMultiSubset_FinalTrainTestEval_Function.R")
source("MaxentMultiSubset_FinalTrainTestEvalAIC_Function.R")
source("MaxentProject_Function.R")
source("MaxentProjectRAW_AICc_Function.R")
source("EvalStatVars_Summary_Function.R")
source("HeatMapDataFrame_Plot_Function.R")
source("VariableSubsetsCorrThresh_Build_Function.R")
source("VariableSubsetSizeCorrThresh_Build_Function.R")
source("GamesHowell.Test.Padj_Function.R")
#startSocketServer(port=8889)
#
# Load needed packages of raster, rgdal, dismo, rjava, and maptools (printouts not shown)
library(dismo)
library(maptools)
library(rgdal)
library(sp)
library(PresenceAbsence)
library(raster)
library(car)
library(caret)
library(ENMeval)
library(devtools)
library(MCDM)
library(dplyr)
#install_github('johnbaums/rmaxent')  # for installing rmaxent
#library(rmaxent)
#library(rasterVis)
#library(viridis)
# Packages used in functions:
library(plyr)
# Create definition for a geographical projection
#crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
#crs.geo <- CRS("+init=epsg:26912")  # geographical, datum NAD83 12N
############################################################################################
## Enter name of image area and type of classifier
## Specify prefix for cateogrical variables (or give complete name if only single
## categorical variable)
#CatVarsPrefix <- "ROADS_CAT"
# If no categorical variables being considered for maxent, need to remove CatVarsPrefix variable
# Specify beta regularization factor MaxEnt base arguments and for tagging output grids
BetaMult <- 2
## Specify base arguments for Maxent settings for evaluation run
#
MaxentBaseArgsIn <- c(paste0("betamultiplier=", BetaMult), "linear=false", "product=false", "threshold=false", "writebackgroundpredictions=true")
#
## Specify maxent arguments for jackknife cross validation run
MaxentEvalArgs1 <- c(paste0("betamultiplier=", BetaMult), "linear=false", "product=false", "threshold=false", "writebackgroundpredictions=true", "outputgrids=FALSE", "replicates=3", "replicatetype=crossvalidate", "responsecurves=true", "jackknife=true")
##
## Check if have categorical variable to add to maxent arguments
if(exists("CatVarsPrefix")==TRUE) {
  MaxentCatArg <- paste0("togglelayertype=", CatVarsPrefix)
  MaxentEvalArgs <- c(MaxentEvalArgs1, MaxentCatArg)
  } else {
  MaxentEvalArgs <- MaxentEvalArgs1
}
# if it exists in the active R session using the remove statement below
#remove(CatVarsPrefix)
# Specify  directory for general data and shapefiles
InDirect <- "E:/FireFrequency"
setwd(InDirect)
# Specify which climatic data set is being used
DataSet <- "All"
# Specify maximum correlation allowed among variables to reduce multicollinearity
CorrThresh <- 0.7
SearchLoops <- 20000 # maximum number of loops to use creating potential variable sets meeting correlation threshold
MaxModSets <- 3000 # maximum number of variable sets needed to keep below correlation threshold
NumberModSets <- 250 # maximum number of variable sets on which to calculate evaluation statistics per subset size
#
# Specify number of observations (NObs): = number of kfold divisions for full variable subset evaluation; and
# = number of top subsets and number of random subsets compared in RSFSA evaluations
NObs <- 10  ## Number of replicates for first two exploration loops
FullNObs <- 250 ## Number of replicates for intensive search looop
NumGrids <- 90  # Total number of variable grids in study
## Specify location of temporary directory where R stores Maxent output if want
## deleted to avoid accumulation of many maxent output folders when running features selection algorithm
#TempDir <- "C:/Users/james/AppData/Local/Temp/raster/maxent"
TempDir <- ""
## Specify extent and various descriptors for maxent run
Extent <- "SWNA"
PresLimit <- 20000 # Limited number of presence points through random sample after spatial thinning
SpatFiltBuff <- 10 # units in km
SpecFileName <- "losevfreq1"
Species <- paste0("LoBrn")
SpecShort <- paste0("LoBr")
SpeciesGen <- "Fire"
#Species <- paste0("LoBrnLTE16yrInt", PresLimit, "ThinR", SpatFiltBuff,"k")
#Species2 <- "LoBrnLTE16yrInt"
Species2 <- "LoBrn"
MaxBackgrndPseudoabsLimit <- 20000 # specify limit for background and pseudoabsence data
# Specify resolution of rasters in square km
ResolutionEnv <- 1.0  ## in units of square kilometers
Climate <- "Current"
TagIn <- "R3" # Using data from a previous run set
Tag <- "R5"
LongSpecies <- "LoBurnSeverityLessorEqualto16yrInterval"
ClassRaster <- "mtlfalllosumr"  # name of classification raster of burn frequencies
ExtentRaster <- "fireevalareaz" # name of background evaluation extent raster
# Specify ranks for AICc and AUC in Multi-Object Optimization ranking
# Specify ranks for AICc and AUC in Multi-Object Optimization ranking
AICcrank <- 0
AUCrank <- 1
if(AICcrank!=0) {
  AICTxt <- "AIC"
} else {
  AICTxt <- ""
}
#
if(AUCrank!=0) {
  AUCTxt <- "AUC"
} else {
  AUCTxt <- ""
}
RankStatistic <- paste0(AUCTxt, AICTxt)

EvalType <- "psa_corrfilt0.7"
AnalysisType <- paste0("RSFSA", RankStatistic)
# Specify Feature Selection Algorithm"
FSAType <- ""
# Specify units of buffers
Units <- "km"
# Specify buffers for spatial thinning and pseudoabsenct points from presence points
SpatFiltBuffs <- "10km"
# Specify if spatial filtering
SpatFilt <- "Filtered" # there is spatial filtering
PsAbsBuff <- 20  # in km
PsAbsBuffs <- "20km"
Conditions <- paste0(PsAbsBuff, "mPsAbsBuff_", SpatFiltBuff, "mSpatFiltBuff")
Conditions2 <- paste("Only Native Data for Model Input")
# Set projection
CRS.WGS84 <- CRS("+init=epsg:4326")
## Specify directory where to be made environmental grids masked to background area extent
## are to be stored
MaskGrids <- paste0(InDirect, "/", "PredictorsMsk", "/", SpeciesGen)
## Read in list of variable names
GridDirect <- "E:/30SecClimIntAnthTopoNorm" # Set environmental grid
setwd(GridDirect)
GridNamesFile <- "Set90ClimaticTopoAnthIndicesNamesR.csv"
GridNames <- as.matrix(read.csv(GridNamesFile))
#str(BandNames)
GridNamesL <- as.list(GridNames)
GridNamesL2 <- toupper(GridNamesL)
VariableNames <- data.frame(GridNamesL2, stringsAsFactors=FALSE)
colnames(VariableNames) <- "VarNames"
VariableNames$VarNames <- gsub("_NS", "", VariableNames$VarNames)
TotVars <- nrow(VariableNames)
ModelName <- paste0("Maxent", DataSet, "_", AnalysisType)
Model <- paste(ModelName, TotVars, "Feature Subset")
##############
## Read in csv file of sets of variables for loop below
VariableSetTypes <- read.csv("VariableSets90.csv", stringsAsFactors=FALSE, na.strings = c("NA", ""))
###############
# Specify names of non-climate grids located in MaskGrids
NonClimateGridNamesFile <- "Set33TopoAnthIndicesNamesR.csv"
NonClimateGridNames <- as.matrix(read.csv(NonClimateGridNamesFile))
NonClimateGridNamesL <- as.list(NonClimateGridNames)
# Specify list of original grid names for all climate variables
FutureGridNamesFile <- "Set57ClimaticIndicesNamesR.csv"
FutureGridNames <- as.matrix(read.csv(paste0(GridDirect, "/", FutureGridNamesFile)))
FutureGridNamesL <- as.list(FutureGridNames)
### Specify files and directories related to future climate
FutClimDir <- paste0("F:/GIS/FutureClimate/", SpeciesGen)
## Establish subdirectory extensions for variables of future climate scenarios
FutClimSubDirList <- c("2050he26", "2050he85", "2070he26", "2070he85")
#
#### IMPORTANT- Specify output names for files and directories
# Identifying part of output rasters
## Specify run id consisting of date
Run <- 1
Set <- Sys.Date()
output1 <- paste0(Species, ModelName, TotVars, "Varsof", NumGrids, "_", Extent, Tag)
# Identifying part of output directory
# Create directory for storing envelope score related grids
dir.create(paste0(InDirect, "/", output1))
OutDirect <- paste0(InDirect, "/", output1)
#
dir.create(paste0(OutDirect, "/", EvalType))
OutDirectpsa <- paste0(OutDirect, "/", EvalType)
#
#########################################################################################
#
##### Create point shapefile of presence points from classification raster
## Read in raster of burn frequency classes
#setwd(InDirect)
#mtlfallsumr <- raster(ClassRaster)
### Convert burn frequency classes greater than zero in raster to spatial point data frames, restrict to 20,000 max
### and spatially thin to 10 km minimum distance
#system.time(burnpres.spdf <- rasterToPoints(mtlfallsumr, fun=function(x){x>0}, spatial=TRUE))
#nrow(burnpres.spdf)
#head(burnpres.spdf)
#tail(burnpres.spdf)
## Write to shapefile for generation of pseudoabsence buffer
#writeOGR(burnpres.spdf, InDirect, paste0(SpecFileName, "unthin"), driver = "ESRI Shapefile", overwrite=TRUE)
#########
### Randomly thin points to 20,000
#burnpresPresLimit.spdf <- burnpres.spdf[sample(1:length(burnpres.spdf), size=min(PresLimit, nrow(burnpres.spdf))), ]
#head(burnpresPresLimit.spdf)
#nrow(burnpresPresLimit.spdf)
#### Spatially thin presence points to minimum distance of SpatFiltBuff (generally 10 km) (takes 9 minutes with 10,000 points)
#system.time(burnpresthin.spdf <- SpatialFilter(burnpresPresLimit.spdf , dist=SpatFiltBuff, mapUnits=F))
#nrow(burnpresthin.spdf)
#head(burnpresthin.spdf)
##plot(burnpresthin.spdf, add=TRUE, col='blue')
## Convert to data frame
#burnpresthin.df <- data.frame(burnpresthin.spdf)
#head(burnpresthin.df)
## Reformat
#burnpresthin.df <- burnpresthin.df[,c(2,3,1)]
#nrow(burnpresthin.df)
#write.table(burnpresthin.df, file=paste0(LongSpecies, "_Thinned", SpatFiltBuff, "km_", PresLimit, "LPresence.csv"))
### Convert to shapefile for display in ArcGIS
#coordinates(burnpresthin.df)=~x+y
#proj4string(burnpresthin.df)=CRS.WGS84
## Convert first to spatial points data frame
#burnpresthin.spdf <- SpatialPointsDataFrame(burnpresthin.df, data.frame(id=1:length(burnpresthin.df)))
#head(burnpresthin.spdf)
## Write to shapefile
#writeOGR(burnpresthin.spdf , InDirect, paste0(SpecFileName, "_", PresLimit, "Thin", SpatFiltBuff,"km"), driver = "ESRI Shapefile")
#####################################################
### Plot presence data with country background
## Read ESRI shapefile into CountryShp below.
#AdminShp <- readShapePoly("admin00.shp", proj4string=CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
##
## xlim sets the LONGITUDE limits and ylim the LATITUDE limits
#dev.new()
#plot(AdminShp, xlim=c(((burnpresPresLimit.spdf@bbox[1,1])-2),((burnpresPresLimit.spdf@bbox[1,2])+2)), ylim=c(((burnpresPresLimit.spdf@bbox[2,1])-2),((burnpresPresLimit.spdf@bbox[2,2])+2)), axes=TRUE, col='light yellow')
#plot(burnpres.spdf, add=TRUE, col='red')
#plot(burnpresPresLimit.spdf, add=TRUE, col='blue')
#plot(burnpresthin.spdf, add=TRUE, col='green')
##
### Optional, use zm() to zoom into points on plot
##library(zoom)
##zm()
####################################################################
#### In ArcGIS, use Pythons script PseudoabsenceBackgroundPointGeneration_PresenceNoThin.py
#### to generate a raster for deriving Pseudoabsence points
#########################################################################################
####
################################################################################
#### Generate and spatially filter pseudoabsence points
####
### Read in above created raster in ArcGIS for pseudoabsence points
#PSARaster <- raster(paste0(SpecFileName, "psr"))
### Convert burn frequency classes greater than zero in raster to spatial point data frames, restrict to 20,000 max
### and spatially thin to 10 km minimum distance
#system.time(psaraw.spdf <- rasterToPoints(PSARaster, spatial=TRUE))
#nrow(psaraw.spdf)
#head(psaraw.spdf)
#tail(psaraw.spdf)
#########
### Randomly thin points to 20,000
#psarawPresLimit.spdf <- psaraw.spdf[sample(1:length(psaraw.spdf), size=max(PresLimit, nrow(burnpres.spdf))), ]
#head(psarawPresLimit.spdf)
#nrow(psarawPresLimit.spdf)
#### Spatially thin presence points to minimum distance of SpatFiltBuff (generally 10 km) (takes 9 minutes with 10,000 points)
#system.time(psathin.spdf <- SpatialFilter(psarawPresLimit.spdf , dist=SpatFiltBuff, mapUnits=F))
#nrow(psathin.spdf)
#head(psathin.spdf)
##plot(burnpresthin.spdf, add=TRUE, col='blue')
## Convert to data frame
#psathin.df <- data.frame(psathin.spdf)
#head(psathin.df)
## Reformat
#psathin.df <- psathin.df[,c(2,3,1)]
#nrow(psathin.df)
#write.table(psathin.df, file=paste0(SpecFileName, "psa.csv"))
### Convert to shapefile for display in ArcGIS
#coordinates(psathin.df)=~x+y
#proj4string(psathin.df)=CRS.WGS84
## Convert first to spatial points data frame
#psathin.spdf <- SpatialPointsDataFrame(psathin.df, data.frame(id=1:length(psathin.df)))
#head(psathin.spdf)
## Write to shapefile
#writeOGR(psathin.spdf, InDirect, paste0(SpecFileName, "psa"), driver = "ESRI Shapefile", overwrite=TRUE)
###############################
### Read in evaluation area raster created in ArcGIS for generating background points
#bkgrndRaster <- raster(ExtentRaster)
### Convert burn frequency classes greater than zero in raster to spatial point data frames, restrict to 20,000 max
### and spatially thin to 10 km minimum distance
#system.time(bkgrndraw.spdf <- rasterToPoints(bkgrndRaster, spatial=TRUE))
#nrow(bkgrndraw.spdf)
#head(bkgrndraw.spdf)
#tail(bkgrndraw.spdf)
#########
### Randomly thin points to 20,000
#bkgrndrawPresLimit.spdf <- bkgrndraw.spdf[sample(1:length(bkgrndraw.spdf), size=max(PresLimit, nrow(burnpres.spdf))), ]
#head(bkgrndrawPresLimit.spdf)
#nrow(bkgrndrawPresLimit.spdf)
#### Spatially thin presence points to minimum distance of SpatFiltBuff (generally 10 km) (takes 9 minutes with 10,000 points)
#system.time(bkgrndthin.spdf <- SpatialFilter(bkgrndrawPresLimit.spdf , dist=SpatFiltBuff, mapUnits=F))
#nrow(bkgrndthin.spdf)
#head(bkgrndthin.spdf)
##plot(burnpresthin.spdf, add=TRUE, col='blue')
## Convert to data frame
#bkgrndthin.df <- data.frame(bkgrndthin.spdf)
#head(bkgrndthin.df)
## Reformat
#bkgrndthin.df <- bkgrndthin.df[,c(2,3,1)]
#nrow(bkgrndthin.df)
#write.table(bkgrndthin.df, file=paste0(SpecFileName, "bkgrnd.csv"))
### Convert to shapefile for display in ArcGIS
#coordinates(bkgrndthin.df)=~x+y
#proj4string(bkgrndthin.df)=CRS.WGS84
## Convert first to spatial points data frame
#bkgrndthin.spdf <- SpatialPointsDataFrame(bkgrndthin.df, data.frame(id=1:length(bkgrndthin.df)))
#head(bkgrndthin.spdf)
## Write to shapefile
#writeOGR(bkgrndthin.spdf, InDirect, paste0(SpecFileName, "bkgrnd"), driver = "ESRI Shapefile")
##
####################################################################
#### Extract variable values for presence, pseudoabsence and background
#### poins for running Samples With Data (SWD) evaluations
#############################################################################################
############
#### Use Arc Python script such as BatchGridClippingSWFL.py to clip
#### environmental grids to background evaluation extent created
#### shapefile created above and put them in directory MaskGrids to be read in here
#### NOTE: Maxent prediction does not work with R grids, must use Arc Grid or Arc Ascii Grid
######################################################################
#setwd(MaskGrids)
#Predictors <- stack(GridNamesL)
#names(Predictors) <- toupper(names(Predictors))
#names(Predictors) <- gsub("_NS", "", names(Predictors))
##plot(Predictors[[89]])
##
#### Read in Pseudoabsence and Background points generated in ArcGIS for each species
### NOTE: If shapefiles are open in ArcGIS, a lock may give errors when reading them back in and
### it may be necessary to close ArcGIS.
#setwd(InDirect)
#PseudoabsencePntShp <- readOGR(".", paste0(SpecFileName, "psa"))  # NOTE: readShapePoly was giving error, so used readOGR instead
##str(PseudoabsencePntShp)
##
#BackgroundPntShp <- readOGR(".", paste0(SpecFileName, "bkgrnd"))
### Read back in spatially thinned presence points
#PresThinShp <- readOGR(".", paste0(SpecFileName, "_", PresLimit, "Thin", SpatFiltBuff,"km"))  # NOTE: readShapePoly was giving so used readOGR instead
#AdminShp <- readShapePoly("admin00.shp", proj4string=CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
###
##str(PresThinShp)
#PresThinShp@data$Latitude
#dev.new()
#plot(AdminShp, xlim=c(((PseudoabsencePntShp@bbox[1,1])-2),((PseudoabsencePntShp@bbox[1,2])+2)), ylim=c(((PseudoabsencePntShp@bbox[2,1])-2),((PseudoabsencePntShp@bbox[2,2])+2)), axes=TRUE, col='light yellow')
#plot(BackgroundPntShp, add=TRUE, col='yellow')
#plot(PresThinShp, add=TRUE, col='red')
#plot(PseudoabsencePntShp, add=TRUE, col='blue')
##
###############################################################################################
##
#### Extract environmental data for presence points
## Extract Predictors for spatial points data frame of thinned roadkill points
#setwd(InDirect)
### Read back in spatially thinned presence points
#PresThinShp <- readOGR(".", paste0(SpecFileName, "_", PresLimit, "Thin", SpatFiltBuff,"km"))  # NOTE: readShapePoly so used readOGR instead
##
#PresenceDat.mat1 <- data.frame(extract(Predictors, PresThinShp))
#PresenceDat.mat1[1:10,]
## Re-associate data with species and x y coordinates
#PresenceDat.mat <- merge(coordinates(PresThinShp), PresenceDat.mat1, by="row.names")
#head(PresenceDat.mat)
#tail(PresenceDat.mat)
## Replace first column of row.names with LongSpecies name
#colnames(PresenceDat.mat)[1] <- c("Species")
#PresenceDat.mat$Species <- LongSpecies
## Convert to data frame
#PresenceDat.df <- data.frame(PresenceDat.mat, stringsAsFactors=FALSE)
#tail(PresenceDat.df)
#nrow(PresenceDat.df)
## Omit any rows with NA values
#PresenceDat.df <- na.omit(PresenceDat.df)
#any(is.na(PresenceDat.df))
## Save extracted data
#setwd(InDirect)
#write.table(PresenceDat.df, sep = ",", col.names=NA, file=paste0(Species, "_presevaldat", DataSet, Tag, "_", TotVars, "Vars.csv"))
######################
## Extract Predictors from pseudoabsence points (may take about 11 minutes)
#setwd(InDirect)
#system.time(PseudoabsenceDat.mat1 <- data.frame(extract(Predictors, PseudoabsencePntShp)))
#PseudoabsenceDat.mat1[1:10,]
## Re-associate data with species and x y coordinates
#PseudoabsenceDat.mat <- merge(coordinates(PseudoabsencePntShp), PseudoabsenceDat.mat1, by="row.names")
#tail(PseudoabsenceDat.mat)
## Replace first column of row.names with LongSpecies name
#colnames(PseudoabsenceDat.mat)[1] <- c("Species")
#PseudoabsenceDat.mat$Species <- "Pseudoabsence"
## Convert to data frame
#PseudoabsenceDat.df <- data.frame(PseudoabsenceDat.mat, stringsAsFactors=FALSE)
#tail(PseudoabsenceDat.df)
#nrow(PseudoabsenceDat.df)
## Omit any rows with NA values
#PseudoabsenceDat.df <- na.omit(PseudoabsenceDat.df)
## Save extracted data
#colnames(PseudoabsenceDat.df) <- gsub("_NS", "", colnames(PseudoabsenceDat.df))
#write.table(PseudoabsenceDat.df, sep = ",", col.names=NA, file=paste0(Species2, "_", PsAbsBuff, Units, "PsAbsBuff_pseudoabsevaldat", DataSet, Tag, "_", TotVars, "Vars.csv"))
###
###########################
## Extract Predictors from background points (may take about 11 minutes)
#system.time(BackgroundDat.mat1 <- data.frame(extract(Predictors, BackgroundPntShp)))
#BackgroundDat.mat1[1:10,]
## Re-associate data with species and x y coordinates
#BackgroundDat.mat <- merge(coordinates(BackgroundPntShp), BackgroundDat.mat1, by="row.names")
#tail(BackgroundDat.mat)
## Replace first column of row.names with LongSpecies name
#colnames(BackgroundDat.mat)[1] <- c("Species")
#BackgroundDat.mat$Species <- "Background"
## Convert to data frame
#BackgroundDat.df <- data.frame(BackgroundDat.mat, stringsAsFactors=FALSE)
#tail(BackgroundDat.df)
#nrow(BackgroundDat.df)
## Omit any rows with NA values
#BackgroundDat.df <- na.omit(BackgroundDat.df)
## Save extracted data
#setwd(InDirect)
#write.table(BackgroundDat.df, sep = ",", col.names=NA, file=paste0(Species, "_backgrounddat", DataSet, Tag, "_", TotVars, "Vars.csv", sep=""))
###
#########################################################################################
# Read back in Samples With Data (SWD) formatted data if needed
setwd(InDirect)
BackgroundDat.df <- as.data.frame(read.csv(paste0(Species, "_backgrounddat", DataSet, TagIn, "_", TotVars, "Vars.csv", sep=""), row.names=1))
nrow(BackgroundDat.df)
head(BackgroundDat.df)
tail(BackgroundDat.df)
if(nrow(BackgroundDat.df)>MaxBackgrndPseudoabsLimit) {
  BackgroundDat.df <- BackgroundDat.df[1:MaxBackgrndPseudoabsLimit ,]
}
##
PseudoabsenceDat.df <- data.frame(read.csv(paste0(Species2, "_", PsAbsBuff, Units, "PsAbsBuff_pseudoabsevaldat", DataSet, TagIn, "_", TotVars, "Vars.csv"), sep=",", row.names=1, stringsAsFactors=FALSE))
head(PseudoabsenceDat.df)
tail(PseudoabsenceDat.df)
if(nrow(PseudoabsenceDat.df)>MaxBackgrndPseudoabsLimit) {
  PseudoabsenceDat.df <- PseudoabsenceDat.df[1:MaxBackgrndPseudoabsLimit,]
}
nrow(PseudoabsenceDat.df)
##
PresenceDat.df <- data.frame(read.csv(paste0(Species, "_presevaldat", DataSet, TagIn, "_", TotVars, "Vars.csv"), sep=",", row.names=1, stringsAsFactors=FALSE))
head(PresenceDat.df)
nrow(PresenceDat.df)
TotPres <- nrow(PresenceDat.df)
#####
#####################################################################################################
### Run Random Subset Features Selection Algorithm (RFSA) evaluations
####################################################################################
### Create and save uneven four kfold partition scheme for multiple rounds of test with presence
### NOTE: Only run partition scheme once for set of Presence and Pseudoabsence data
#Create a kfold partition with 1/3 data for Wrapper training, 1/6 for Wrapper testing and 1/3 for final training and 1/6 for
# final testing both presence and pseudoabsence data
# Find lengths of thirds and sixths for PresenceDat.df
#PresWrapperTrainNum <- Round2(0.3333*nrow(PresenceDat.df),0)
#PresFinalTrainNum <- PresWrapperTrainNum
#PresWrapperTestNum <- Round2(0.1667*nrow(PresenceDat.df),0)
#PresFinalTestNum <- PresWrapperTestNum
#PresFinalTestNum <- nrow(PresenceDat.df) - PresWrapperTrainNum - PresFinalTrainNum - PresWrapperTestNum
#kfoldgrppl <- split(1:nrow(PresenceDat.df), sample(rep(1:4, c(PresWrapperTrainNum, PresWrapperTestNum, PresFinalTrainNum, PresFinalTestNum))))
#kfoldgrpp <- rep(1,nrow(PresenceDat.df))
#kfoldgrpp[kfoldgrppl[[1]]]=1
#kfoldgrpp[kfoldgrppl[[2]]]=2
#kfoldgrpp[kfoldgrppl[[3]]]=3
#kfoldgrpp[kfoldgrppl[[4]]]=4
#occurrences <- table(unlist(kfoldgrpp))
## Repeat for Pseudoabsence data
#PsaWrapperTrainNum <- Round2(0.3333*nrow(PseudoabsenceDat.df),0)
#PsaFinalTrainNum <- PsaWrapperTrainNum
#PsaWrapperTestNum <- Round2(0.1667*nrow(PseudoabsenceDat.df),0)
#PsaFinalTestNum <- PsaWrapperTestNum
#PsaFinalTestNum <- nrow(PseudoabsenceDat.df) - PsaWrapperTrainNum - PsaFinalTrainNum - PsaWrapperTestNum
#kfoldgrpal <- split(1:nrow(PseudoabsenceDat.df), sample(rep(1:4, c(PsaWrapperTrainNum, PsaWrapperTestNum, PsaFinalTrainNum, PsaFinalTestNum))))
#kfoldgrpa <- rep(1,nrow(PseudoabsenceDat.df))
#kfoldgrpa[kfoldgrpal[[1]]]=1
#kfoldgrpa[kfoldgrpal[[2]]]=2
#kfoldgrpa[kfoldgrpal[[3]]]=3
#kfoldgrpa[kfoldgrpal[[4]]]=4
#occurrences <- table(unlist(kfoldgrpa))
##### Save kfold partition scheme for later testing
#setwd(InDirect)
#write.table(data.frame(kfoldgrpp), file=paste0(Species, "PresenceDatRSFSAKfold", ".csv"), sep=",")
#write.table(data.frame(kfoldgrpa), file=paste0(Species, PsAbsBuffs, "PseudoabsenceDatRSFSAKfold", ".csv"), sep=",")
##
####################################################################################
# Read in k-fold partition scheme
setwd(InDirect)
kfoldgrpp <- as.vector(unlist(read.csv(paste0(Species, "PresenceDatRSFSAKfold", ".csv"))))
length(kfoldgrpp)
kfoldgrpa <- as.vector(unlist(read.csv(paste0(Species, PsAbsBuffs, "PseudoabsenceDatRSFSAKfold", ".csv"))))
if(nrow(PseudoabsenceDat.df)>MaxBackgrndPseudoabsLimit) {
  NObskfoldgrpa <- NObskfoldgrpa[1:MaxBackgrndPseudoabsLimit ]
}
length(kfoldgrpa)
#
SubsetSizes <- c(3, 6, 8, 10, 12, 15, 20, 25)
#SubsetSizes <- c(25, 30, 35, 40)
SetType <- 10 # use full set
#SubsetSizes <- c(8,12,15,20,25)
#SubsetSizes <- c(3)
length(SubsetSizes)
OutCount <- 0
## Loop takes 80 to 87 min on new laptop
####################################################################################
t1 <- Sys.time()
MaxentKeepEvalAllStatsL <- list()
MaxentKeepEvalAllStatSummL <- list()
for(SubsetLoop in SubsetSizes) {
  #SubsetLoop <- 3
  # Generate random subsets of 2 of 19 variables
  OutCount <- OutCount + 1
  VariableNamesIn <- VariableNames
  TotVars <- nrow(VariableNamesIn)
  SubsetVariableNumber <- SubsetLoop
  NumberModSets <- 250
  #NumberModSets <- 20
  TestDataType <- "Wrapper"
  OutDirectIn <- OutDirectpsa
  RunType <- "RndFSAFilt"
  Run <- 1
  #
  VariableNamesSelect <- data.frame(na.omit(toupper(VariableSetTypes[,SetType])), stringsAsFactors=FALSE)
  if(SetType < 10) {
    SetNumber <- paste0("0", SetType)
  } else {
    SetNumber <- as.character(SetType)
  }
  SetName <- paste0(SetNumber, "_", colnames(VariableSetTypes[SetType]))
  SetNameF <- colnames(VariableSetTypes[SetType])
  colnames(VariableNamesSelect) <- SetName
  InitVars <- nrow(VariableNamesSelect)
  #
  ## Using background values data, find variables with correlation below 0.7
  ## First calculate correlation matrix
  if(OutCount == 1) {
    # Delete first three columns to just keep environmental data
    BackgroundDatCorr.df <- BackgroundDat.df[,-1:-3]
    head(BackgroundDatCorr.df)
    # Subset data by chosen variables
    BackgroundDatCorr.df <- BackgroundDatCorr.df[, unlist(VariableNamesSelect)]
    head(BackgroundDatCorr.df)
    ### Run Spearman Rank Correlation and produce a correlation matrix on all variables
    AbsValsCor.mat <- abs(cor(BackgroundDatCorr.df, method="spearman"))
    # Find mean value of correlations in set
    LoAbsValsCor.mat <- (lower.tri(AbsValsCor.mat)*AbsValsCor.mat)
    MeanCor <- mean(LoAbsValsCor.mat)
    # Save correlation matrix
    write.table(AbsValsCor.mat, sep = ",", col.names=NA, file=paste0(Species, "_DataSet", SetName, "_CorrelationBackgroundEnvironVariables.csv"))
    ## Plot correlation matrix
    DataFrameName <- data.frame(AbsValsCor.mat)
    PlotName <- paste0(LongSpecies, "_", SetNameF, "_CorrelationBackgroundEnvironVariables.tiff")
    PlotTitle <- paste0(LongSpecies, ": VariableCorrelations ", SetNameF, " from Background\n                                                    (Mean Correlation of ", Round2(MeanCor,2), ")\n")
    ## Call Heat Map Function
    HeatMapDataFrame.Plot(DataFrameName, PlotName, PlotTitle, OutDirectIn)
  }
  ## Find subsets meeting correlation criteria
  CorrelationsIn <- AbsValsCor.mat
  SearchLoopsIn <- SearchLoops
  MaxModSetsIn <- NumberModSets # maximum number of sets needed to keep
  CorrThreshIn <- CorrThresh
  SetNameIn <- SetName
  OutDirectIn <- OutDirectpsa
  #
  system.time(VariableSubsets <- VariableSubsetSizeCorrThresh.Build(CorrelationsIn, SearchLoopsIn, SubsetVariableNumber, MaxModSetsIn, CorrThreshIn, SetNameIn, OutDirectIn))
  tail(VariableSubsets)
  Sets <- nrow(VariableSubsets)
  ### Make Heat Map of above matrix for correlation values
  # Convert correlation matrix to data frame
  ## Use function to round up from .5 from http://stackoverflow.com/questions/12688717/round-up-from-5-in-r
  Round2 <- function(x, n) {
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5
    z = trunc(z)
    z = z/10^n
    z*posneg
  }
  ###
  #
  #FSAType <- paste0("Wrap")
  SetName2 <- paste0(SetName, InitVars, "CorrFilt", SubsetVariableNumber)
  TestDataType <- "Wrapper"
  OutDirectIn <- OutDirectpsa
  RunType <- paste0("RndFSAFilt_CorrThresh", CorrThresh)
  Run <- 1
  #########################
  # Save row names to use in output
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
  #TotVars <- 260
  #
  #NumberModSets <- 500
  #FSAType <- "NoWrapper"
  output2 <- paste0(NumberModSets,"Setsof", Subset, SubsetVariableNumber, "of", TotVars, "Vars")
  OutDirect2 <- paste0(OutDirectpsa, "/", output2)
  dir.create(OutDirect2)
  setwd(OutDirect2)
  OutDirectIn <- OutDirect2
  #
  ## Evaluate random subsets for NumberModSets number of subsets
  # Takes 4.48 minutes for 1000 runs
  #ncol(combn(19,5))
  #  t1a <- Sys.time()
  # Designate Wrapper Training data for Wrapper test
  MaxentPresTrainData <- PresenceDat.df[kfoldgrpp==1,]
  head(MaxentPresTrainData)
  nrow(MaxentPresTrainData)
  MaxentAbsTrainData <- BackgroundDat.df
  nrow(MaxentAbsTrainData)
  head(MaxentAbsTrainData)
  WrapperTrainSWD <- rbind(MaxentPresTrainData[,4:ncol(MaxentPresTrainData)], MaxentAbsTrainData[,4:ncol(MaxentAbsTrainData)])
  head(WrapperTrainSWD)
  tail(WrapperTrainSWD)
  WrapperTrainPresID <- data.frame(rep(1,nrow(MaxentPresTrainData)))
  colnames(WrapperTrainPresID) <- "ID"
  WrapperTrainAbsID <- data.frame(rep(0,nrow(MaxentAbsTrainData)))
  colnames(WrapperTrainAbsID) <- "ID"
  WrapperTrainPresAbsID <- rbind(WrapperTrainPresID, WrapperTrainAbsID)
  head(WrapperTrainPresAbsID)
  tail(WrapperTrainPresAbsID)
  VariableSubsetsIn <- VariableSubsets
  SetRunID <- ""
  #
#  NOTE: The next line can be commented if just re-ranking the subsets by a different evaluation statistic
  system.time(MaxentFiltSubsetEvalStats.df <- MaxentMultiSubset.WrapperTestEvalAIC(Species, VariableNamesIn, SubsetVariableNumber, TotPres, VariableSubsetsIn, SetRunID, WrapperTrainSWD, WrapperTrainPresAbsID, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrpp, kfoldgrpa, OutDirectIn, MaxentArgsIn=MaxentBaseArgsIn))
  ######################
  ## Rank top NObs subsets by Wrapper TSS and obtain training and final testing TSS for evaluation of each set
  ## and for ensemble set
  ##
  # If necessary, read back in evaluation stats
  ## NOTE: The next two lines are for running the loop with reranking and testing only
  #OutDirect2Orig <- gsub("AIC", "AUC", OutDirect2)
  #setwd(OutDirect2Orig)
  #setwd(OutDirect2)
  Run <- 1
  Sets <- NumberModSets
  TestDataType <- "Wrapper"
  MaxentFiltSubsetEvalStats.df <- data.frame(read.csv(paste0(Species, "MaxentResults_Wrapper_", TotVars, "TotVars_", SubsetVariableNumber, "Vars_", Sets, "_", SetRunID, ".csv"), sep=",", row.names=1, stringsAsFactors=FALSE))
  head(MaxentFiltSubsetEvalStats.df)
  tail(MaxentFiltSubsetEvalStats.df)
  # Delete any subsets with "Inf" for AICc_bg (means too many parameters in model relative to observations)
  MaxentFiltSubsetEvalStats.df2 <- subset(MaxentFiltSubsetEvalStats.df, AICc_bg!="Inf")
  NewSets1 <- nrow(MaxentFiltSubsetEvalStats.df2)
  ## Use Multi-Objective Optimization from MCDM package to sort by two criteria of AUC and AICc_bg
  decision.mat <- as.matrix(MaxentFiltSubsetEvalStats.df2[,c(11,7)])
  head(decision.mat)
  # Normalize AICc_bg and AUC from zero to 1
  normalizef <- function(x) c((x+(max(x)-min(x))-max(x))/(max(x)-min(x)))
  decision.matn <- as.matrix(apply(decision.mat,2,normalizef))
  head(decision.matn)
  # Assign weights to criteria
  weightscrit <- c(AICcrank, AUCrank)
  # Assign whether cost "min", or benefit "max"
  cb <- c("min", "max")
  # Rank criteria
  MMOORA.Rank <- MMOORA(decision.matn, weightscrit, cb)
  head(MMOORA.Rank)
  MaxentFiltSubsetEvalStats.df4 <- MaxentFiltSubsetEvalStats.df2
  # Join MMOORA rank to original data
  MaxentFiltSubsetEvalStats.df4$MMOORA_Rank <- MMOORA.Rank[,8]
  head(MaxentFiltSubsetEvalStats.df4)
  # Sort data by MMOORA rank
  MaxentFiltSubsetEvalStats.df4s <- MaxentFiltSubsetEvalStats.df4[order(MaxentFiltSubsetEvalStats.df4$MMOORA_Rank),]
  head(MaxentFiltSubsetEvalStats.df4s)
  ncol(MaxentFiltSubsetEvalStats.df4s)
  # Delete rank column
  MaxentFiltSubsetEvalStats.df5 <- MaxentFiltSubsetEvalStats.df4s[,-14]
  head(MaxentFiltSubsetEvalStats.df5)
  # Keep top NObs selected subsets
  MaxentTopNObsFiltSubsetEvalStats.df <- data.frame(MaxentFiltSubsetEvalStats.df5[1:NObs,], stringsAsFactors=FALSE)
  MaxentTopNObsFiltSubsetEvalStats.df
  TopNObsSubsets <- data.frame(MaxentTopNObsFiltSubsetEvalStats.df[,1], stringsAsFactors=FALSE)
  ##############################
  # Keep random NObs random subsets
  MaxentRandNObsFiltSubsetEvalStats.df <- data.frame(MaxentFiltSubsetEvalStats.df2[sample(nrow(MaxentFiltSubsetEvalStats.df2), NObs), ], stringsAsFactors=FALSE)
  RandNObsSubsets <- data.frame(MaxentRandNObsFiltSubsetEvalStats.df[,1], stringsAsFactors=FALSE)
  colnames(RandNObsSubsets) <- "VarNames"
  # Use difference between evaluation statistics between train and test data to evaluate overfitting
  ###############################################################################
  ## Test Top NObs Ranking Subsets and Random NObs Subsets
  ###############################################################################
  #####################################
  setwd(OutDirectpsa)
  SubsetRuns <- list(TopNObsSubsets, RandNObsSubsets)
  DataSetTypes <- c(paste0("Top",NObs), paste0("Random", NObs))
  MaxentKeepEvalStatsL <- list()
  #
  for (j in 1:2) {
    #N3Subsets <- SubsetRuns[[1]]
    #j=1
    VariableSubsets <- data.frame(SubsetRuns[[j]], stringsAsFactors=FALSE)
    DataSetType <- DataSetTypes[j]
    # Identify Wrapper data evaluation statistics from above in loop
    if(DataSetType==paste0("Top", NObs)) {
      MaxentTestSubsetEvalStats.df <- MaxentTopNObsFiltSubsetEvalStats.df
    } else {
      MaxentTestSubsetEvalStats.df <- MaxentRandNObsFiltSubsetEvalStats.df
    }
    MaxentTestSubsetEvalStats.df$DataType <- paste0(DataSetType, "WrapperTest")
    ##
    ## Evaluate random subsets for NumberModSets number of subsets
    system.time(MaxentTrainSubsetEvalStats.df <- MaxentMultiSubset.WrapperTrainEvalAIC(Species, VariableNamesIn, SubsetVariableNumber, TotPres, VariableSubsets, WrapperTrainSWD, WrapperTrainPresAbsID, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrpp, kfoldgrpa, OutDirectIn, MaxentArgsIn=MaxentBaseArgsIn, DataSetType=DataSetType))
    ## Separate out training and testing data
    MaxentTrainSubsetEvalStats.df$DataType <- paste0(DataSetType, "WrapperTrain")
    # Calculate difference between test and train statistics for overfitting
    MaxentDiffTestTrainSubsetEvalStats.df <- MaxentTestSubsetEvalStats.df
    MaxentDiffTestTrainSubsetEvalStats.df[,5:13] <- MaxentTrainSubsetEvalStats.df[,5:13] - MaxentTestSubsetEvalStats.df[,5:13]
    MaxentDiffTestTrainSubsetEvalStats.df$DataType <- paste0(DataSetType, "WrapperDiffTestTrain")
    #
    MaxentSubsetEvalStats.df <- do.call(rbind, list(MaxentTrainSubsetEvalStats.df, MaxentTestSubsetEvalStats.df, MaxentDiffTestTrainSubsetEvalStats.df))
    MaxentSubsetEvalStats.df$SetName <- SetName
    MaxentSubsetEvalStats.df$Set <- DataSetType
    MaxentSubsetEvalStats.df$TotSubsets <- Sets
    MaxentSubsetEvalStats.df$PresTrainPoints <- nrow(PresenceDat.df[kfoldgrpp==1,])
    MaxentSubsetEvalStats.df$Run <- Run
    MaxentSubsetEvalStats.df <- MaxentSubsetEvalStats.df[, colnames(MaxentSubsetEvalStats.df)[c(1:4,14:18,5:13)]]
    #
    MaxentKeepEvalStatsL[[j]] <- MaxentSubsetEvalStats.df
  }
  ###############
  ## Designate Final Training data for Final Test
  MaxentPresTrainData <- PresenceDat.df[kfoldgrpp==3,]
  head(MaxentPresTrainData)
  nrow(MaxentPresTrainData)
  MaxentAbsTrainData <- BackgroundDat.df
  nrow(MaxentAbsTrainData)
  head(MaxentAbsTrainData)
  FinalTrainSWD <- rbind(MaxentPresTrainData[,4:ncol(MaxentPresTrainData)], MaxentAbsTrainData[,4:ncol(MaxentAbsTrainData)])
  head(FinalTrainSWD)
  tail(FinalTrainSWD)
  FinalTrainPresID <- data.frame(rep(1,nrow(MaxentPresTrainData)))
  colnames(FinalTrainPresID) <- "ID"
  FinalTrainAbsID <- data.frame(rep(0,nrow(MaxentAbsTrainData)))
  colnames(FinalTrainAbsID) <- "ID"
  FinalTrainPresAbsID <- rbind(FinalTrainPresID, FinalTrainAbsID)
  head(FinalTrainPresAbsID)
  tail(FinalTrainPresAbsID)
  #
  setwd(OutDirectpsa)
  SubsetRuns <- list(TopNObsSubsets, RandNObsSubsets)
  DataSetTypes <- c(paste0("Top",NObs), paste0("Random", NObs))
  MaxentKeepEvalStatsL2 <- list()
  #
  for (j in 1:2) {
    #N3Subsets <- SubsetRuns[[1]]
    #j=1
    VariableSubsets <- data.frame(SubsetRuns[[j]], stringsAsFactors=FALSE)
    DataSetType <- DataSetTypes[j]
    ## Evaluate random subsets for NumberModSets number of subsets
    system.time(MaxentTrainTestSubsetEvalStats.df <- MaxentMultiSubset.FinalTrainTestEvalAIC(Species, VariableNamesIn, SubsetVariableNumber, TotPres, VariableSubsets, FinalTrainSWD, FinalTrainPresAbsID, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrpp, kfoldgrpa, OutDirectIn, MaxentArgsIn=MaxentBaseArgsIn, DataSetType=DataSetType))
    ## Separate out training and testing data
    MaxentTrainSubsetEvalStats.df <- MaxentTrainTestSubsetEvalStats.df[which(MaxentTrainTestSubsetEvalStats.df$DataType=="FinalTrain"),]
    MaxentTrainSubsetEvalStats.df$DataType <- paste0(DataSetType, "FinalTrain")
    MaxentTestSubsetEvalStats.df <- MaxentTrainTestSubsetEvalStats.df[which(MaxentTrainTestSubsetEvalStats.df$DataType=="FinalTest"),]
    MaxentTestSubsetEvalStats.df$DataType <- paste0(DataSetType, "FinalTest")
    # Calculate difference between test and train statistics for overfitting
    MaxentDiffTestTrainSubsetEvalStats.df <- MaxentTestSubsetEvalStats.df
    MaxentDiffTestTrainSubsetEvalStats.df[,5:13] <- MaxentTrainSubsetEvalStats.df[,5:13] - MaxentTestSubsetEvalStats.df[,5:13]
    MaxentDiffTestTrainSubsetEvalStats.df$DataType <- paste0(DataSetType, "FinalDiffTestTrain")
    #
    MaxentSubsetEvalStats.df <- do.call(rbind, list(MaxentTrainSubsetEvalStats.df, MaxentTestSubsetEvalStats.df, MaxentDiffTestTrainSubsetEvalStats.df))
    MaxentSubsetEvalStats.df$SetName <- SetName
    MaxentSubsetEvalStats.df$Set <- DataSetType
    MaxentSubsetEvalStats.df$TotSubsets <- Sets
    MaxentSubsetEvalStats.df$PresTrainPoints <- nrow(PresenceDat.df[kfoldgrpp==3,])
    MaxentSubsetEvalStats.df$Run <- Run
    MaxentSubsetEvalStats.df <- MaxentSubsetEvalStats.df[, colnames(MaxentSubsetEvalStats.df)[c(1:4,14:18,5:13)]]
    #
    MaxentKeepEvalStatsL2[[j]] <- MaxentSubsetEvalStats.df
  }
  ###############
  MaxentKeepEvalStats.df1 <- do.call(rbind, MaxentKeepEvalStatsL)
  MaxentKeepEvalStats.df2 <- do.call(rbind, MaxentKeepEvalStatsL2)
  MaxentKeepEvalStats.df3 <- rbind(MaxentKeepEvalStats.df1, MaxentKeepEvalStats.df2)
  MaxentKeepEvalAllStatsL[[OutCount]] <- MaxentKeepEvalStats.df3
  ### Summarize Output
  ###Calculate Mean and Standard Deviation Values
  RunType <- "Final"
  EvaluationStats <- MaxentKeepEvalStats.df3
  ncol(EvaluationStats)
  EvaluationStats$Model <- ModelName
  EvaluationStats <- EvaluationStats[, colnames(EvaluationStats)[c(19,1:18)]]
  ##
  #FSAType <- RankStatistic
  OutName <- paste0(Species, ModelName, "FullSetSummaryStats_", SubsetVariableNumber, "of", TotVars, "Vars", FSAType,"_", RunType, ".csv")
  SortGroups <- c("Model", "DataType", "Set", "SetName", "SubsetVariableNumber", "TotSubsets")
  ## All statistics, and only statistics, should be at and after the column specified below
  StatVarFirstColumn <- 9
  OutDirectIn <- OutDirectpsa
  EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
  head(EvaluationStatsSummary.df)
  MaxentKeepEvalAllStatSummL[[OutCount]] <- EvaluationStatsSummary.df
  ######
  MaxentKeepEvalStatsAllSets.df <- do.call(rbind, MaxentKeepEvalAllStatSummL)
  ncol(MaxentKeepEvalStatsAllSets.df)
  ## Save data
  write.table(MaxentKeepEvalStatsAllSets.df, file=paste0(Species, ModelName, "_FullSet", "SummaryStats_of", min(SubsetSizes), "to", max(SubsetSizes), "Vars", FSAType,"_", RunType, ".csv"), sep=",", col.names=NA)
  ### Summarize and save raw data as loop progresses
  MaxentKeepEvalStats.df <- do.call(rbind, MaxentKeepEvalAllStatsL)
  # Save data
  setwd(OutDirectpsa)
  write.table(MaxentKeepEvalStats.df , file=paste0(Species, ModelName, NObs, "_FullSetTrainVsTestStats_", min(SubsetSizes), "to", max(SubsetSizes), "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff.csv"), sep=",", col.names=NA)
}
#
####################################
t2 <- Sys.time()
##################
difftime(t2,t1, units = "mins")
setwd(OutDirectpsa)
RunType <- "Final"
# Read data back in and sort
MaxentKeepEvalSumStats.df1 <- data.frame(read.csv(paste0(Species, ModelName, "_FullSet", "SummaryStats_of", min(SubsetSizes), "to", max(SubsetSizes), "Vars", FSAType,"_", RunType, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
## Sort data by Statistic and then by DataType for ease of plotting in excel
MaxentKeepEvalSumStats.dfs1 <- arrange(MaxentKeepEvalSumStats.df1, Statistic, DataType)
MaxentKeepEvalSumStats.dfs <- arrange(MaxentKeepEvalSumStats.dfs1 , DataType)
head(MaxentKeepEvalSumStats.dfs)
#
write.table(MaxentKeepEvalSumStats.dfs, file=paste0(Species, ModelName, "_FullSet", "SummaryStats_of", min(SubsetSizes), "to", max(SubsetSizes), "Vars", FSAType,"_", RunType, ".csv"), sep=",", col.names=NA)



######################################################################################
## Assemble data for bar charts and ANOVA
##########################################################################################
setwd(OutDirectpsa)
# Read in results for first FSA run to get output on the leave one out variable run with all variables
#FSAType <- RankStatistic
#SubsetSizes <- c(3,6,8,10,12,15,20,25)
#
MaxentKeepEvalStats.df1 <- data.frame(read.csv(paste0(Species, ModelName, NObs, "_FullSetTrainVsTestStats_", min(SubsetSizes), "to", max(SubsetSizes), "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
nrow(MaxentKeepEvalStats.df1)
ncol(MaxentKeepEvalStats.df1)
head(MaxentKeepEvalStats.df1)
tail(MaxentKeepEvalStats.df1)
# Keep only DataTypes including string "FinalTest"
MaxentKeepEvalStats.df2 <- MaxentKeepEvalStats.df1[grep("FinalTest", MaxentKeepEvalStats.df1$DataType), ]
MaxentKeepEvalStats.df2[1:50,]
tail(MaxentKeepEvalStats.df2)
# Create DataTypeN combining SubsetVariableNumber and DataType
MaxentKeepEvalStats.df2$DataTypeN <- paste0(MaxentKeepEvalStats.df2$DataType, "_",MaxentKeepEvalStats.df2$SubsetVariableNumber)
head(MaxentKeepEvalStats.df2)
ncol(MaxentKeepEvalStats.df2)
# Merge SubsetVariableNumber with SetName
MaxentKeepEvalStats.df2$SetName <- paste0(MaxentKeepEvalStats.df2$SetName, "_", as.character(MaxentKeepEvalStats.df2$SubsetVariableNumber))
# Reorder Columns
MaxentKeepEvalStats.df2 <- MaxentKeepEvalStats.df2[,c(19,1:18)]
###Calculate Mean and Standard Deviation Values
RunType <- "Final"
EvaluationStats <- MaxentKeepEvalStats.df2
# Delete any subsets with "Inf" for AICc_bg (means too many parameters in model relative to observations)
EvaluationStats <- subset(EvaluationStats, AICc_bg!="Inf")
# Replace NAN values with zero
EvaluationStats[is.na(EvaluationStats)] <- 0
#
OutName <- paste0(Species, RankStatistic, "Maxent_VariousSetsof_", SubsetVariableNumber, "_", NObs, "N_FinalTestStatSummary.csv")
SortGroups <- c("DataTypeN", "DataType", "SetName", "SubsetVariableNumber")
## All statistics, and only statistics, should be at and after the column specified below
StatVarFirstColumn <- 7
OutDirectIn <- OutDirectpsa
EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
head(EvaluationStatsSummary.df)
##
EvalStatsOut <- EvaluationStatsSummary.df
## Sort data by DataType for ease of plotting in excel
EvalStatsOuts <- arrange(EvalStatsOut, DataType, Statistic, SubsetVariableNumber)
head(EvalStatsOuts)
write.table(EvalStatsOuts, file=paste0(Species, RankStatistic, "MaxentFullSet", NObs, "StatSetRunsbySubsetSizes_FinalTestSummaryTable.csv"), sep=",", col.names=NA)
##
#



#####################################################################################################
### Run Random Subset Features Selection Algorithm (RFSA) evaluations for 10,000 subsets of chosen subset size
####################################################################################
## Loop takes 2.6 hours for 10,200 NumberModsets of 8 variables for flycatcher on new laptop
## Loop takes 2.8 hours for 10,200 NumberModsets of 10 variables for flycatcher on new laptop
## Loop takes 3.75 hours for 10,200 NumberModsets of 15 variables for flycatcher on new laptop
## Loop takes 4.2 hours for 10,150 NumberModsets of 15 variables for flycatcher on new laptop
## Loop takes 3.6 hours for 9,000 NumberModsets of 12 variables for flycatcher on new laptop
## Loop takes 3.07 hours for 9,000 NumberModSets of 6 variables for lo brn fire on new laptop
## Loop takes 10.6 hours for 9,000 NumberModSets of 25 variables for lo brn fire on new laptop
###############################################################################
OutCount <- 0
##
MaxentKeepEvalAllStatsL <- list()
MaxentKeepEvalAllStatSummL <- list()
NumberModSetsIn <- 3150 # Number sets per rep: extra sets if needed due to AICc values deleted- recommend no more than 5050 sets at once in case of crash
NumberModSets <- 3000  # Number sets per rep
#NumberModSetsIn <- 33 # extra sets if needed due to AICc values deleted
#NumberModSets <- 30
SetType.Vec <- c(10)
#SetType.Vec <- c(1)
#SetRunsList <- c(1)
#NumberModSetsIn <- 5200 # extra sets if needed due to AICc values deleted
#NumberModSets <- 5000
#
FinalModelVariables <- 15
#choose(90,6)
#choose(8,4)
SubsetVariableNumber <- FinalModelVariables
Reps <- 3 # Number of randomizations for Wrapper/final data in feature selection
RepList <- seq(1:Reps)
FullNObs <- 250
PlusNum <- 100 # extra random subsets
#FullNObs <- NumberModSets/3
#PlusNum <- 1
#########################################################
t1 <- Sys.time()
for(j in SetType.Vec) {
#for(j in 1:2) {
  #j=10
  VariableNamesIn <- VariableNames
  TotVars <- nrow(VariableNamesIn)
  #
  VariableNamesSelect <- data.frame(na.omit(toupper(VariableSetTypes[,j])), stringsAsFactors=FALSE)
  if(j < 10) {
    SetNumber <- paste0("0", j)
  } else {
    SetNumber <- as.character(j)
  }
  SetName <- paste0(SetNumber, "_", colnames(VariableSetTypes[j]))
  SetNameF <- colnames(VariableSetTypes[j])
  colnames(VariableNamesSelect) <- SetName
  InitVars <- nrow(VariableNamesSelect)
  #
  ## Using background values data, find variables with correlation below 0.7
  # Delete first three columns to just keep environmental data
  BackgroundDatCorr.df <- BackgroundDat.df[,-1:-3]
  head(BackgroundDatCorr.df)
  # Subset data by chosen variables
  BackgroundDatCorr.df <- BackgroundDatCorr.df[, unlist(VariableNamesSelect)]
  head(BackgroundDatCorr.df)
  ### Run Spearman Rank Correlation and produce a correlation matrix on all variables
  AbsValsCor.mat <- abs(cor(BackgroundDatCorr.df, method="spearman"))
  # Find mean value of correlations in set
  LoAbsValsCor.mat <- (lower.tri(AbsValsCor.mat)*AbsValsCor.mat)
  MeanCor <- mean(LoAbsValsCor.mat)
  # Save correlation matrix
  write.table(AbsValsCor.mat, sep = ",", col.names=NA, file=paste0(Species, "_DataSet", SetName, "_CorrelationBackgroundEnvironVariables.csv"))
  ##
  CorrelationsIn <- AbsValsCor.mat
  SearchLoopsIn <- SearchLoops
  MaxModSetsIn <- NumberModSetsIn * Reps # maximum number of sets needed to keep
  CorrThreshIn <- 0.7
  SetNameIn <- SetName
  OutDirectIn <- OutDirectpsa
  #
  system.time(VariableSubsets <- VariableSubsetSizeCorrThresh.Build(CorrelationsIn, SearchLoopsIn, SubsetVariableNumber, MaxModSetsIn, CorrThreshIn, SetNameIn, OutDirectIn))
  tail(VariableSubsets)
  # Find number of variables in subsets
  VariableNamesEx <- c(unlist(VariableSubsets[1,1]))
  # Check if dash used to separate variables
  if(grepl("-", VariableNamesEx)==TRUE) {
    VarNamesEx <- unlist(strsplit(VariableNamesEx, "-"))
  } else {
    VarNamesEx <- unlist(VariableNamesEx)
  }
#  ### Make Heat Map of above matrix for correlation values
#  # Convert correlation matrix to data frame
#  ## Use function to round up from .5 from http://stackoverflow.com/questions/12688717/round-up-from-5-in-r
#  Round2 <- function(x, n) {
#    posneg = sign(x)
#    z = abs(x)*10^n
#    z = z + 0.5
#    z = trunc(z)
#    z = z/10^n
#    z*posneg
#  }
#  DataFrameName <- data.frame(AbsValsCor.mat)
#  PlotName <- paste0(LongSpecies, "_", SetNameF, "_CorrelationBackgroundEnvironVariables.tiff")
#  PlotTitle <- paste0(LongSpecies, ": VariableCorrelations ", SetNameF, " from Background\n                                                    (Mean Correlation of ", Round2(MeanCor,2), ")\n")
#  ## Call Heat Map Function
#  HeatMapDataFrame.Plot(DataFrameName, PlotName, PlotTitle, OutDirectIn)
  ###
  #
  #FSAType <- RankStatistic
  SetName2 <- paste0(SetName, InitVars, "CorrFilt", SubsetVariableNumber)
  TestDataType <- "Wrapper"
  #
  RunType <- paste0("RndFSAFilt_CorrThresh", CorrThresh)
  Run <- 1
  #
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
  ##########################
  for(Rep in RepList) {
    #Rep=2
    ### Create and save a different set for each of three reps of uneven four kfold partition scheme for multiple rounds of test with presence
    ### NOTE: Only run partition scheme once for set of Presence and Pseudoabsence data
    #Create a kfold partition with 1/3 data for Wrapper training, 1/6 for Wrapper testing and 1/34 for final training and 1/6 for
    # final testing both presence and pseudoabsence data
    # Find lengths of thirds and sixths for PresenceDat.df
#    PresWrapperTrainNum <- Round2(0.3333*nrow(PresenceDat.df),0)
#    PresFinalTrainNum <- PresWrapperTrainNum
#    PresWrapperTestNum <- Round2(0.1667*nrow(PresenceDat.df),0)
#    PresFinalTestNum <- PresWrapperTestNum
#    PresFinalTestNum <- nrow(PresenceDat.df) - PresWrapperTrainNum - PresFinalTrainNum - PresWrapperTestNum
#    kfoldgrppl <- split(1:nrow(PresenceDat.df), sample(rep(1:4, c(PresWrapperTrainNum, PresWrapperTestNum, PresFinalTrainNum, PresFinalTestNum))))
#    kfoldgrpp <- rep(1,nrow(PresenceDat.df))
#    kfoldgrpp[kfoldgrppl[[1]]]=1
#    kfoldgrpp[kfoldgrppl[[2]]]=2
#    kfoldgrpp[kfoldgrppl[[3]]]=3
#    kfoldgrpp[kfoldgrppl[[4]]]=4
#    occurrences <- table(unlist(kfoldgrpp))
#    # Repeat for Pseudoabsence data
#    PsaWrapperTrainNum <- Round2(0.3333*nrow(PseudoabsenceDat.df),0)
#    PsaFinalTrainNum <- PsaWrapperTrainNum
#    PsaWrapperTestNum <- Round2(0.1667*nrow(PseudoabsenceDat.df),0)
#    PsaFinalTestNum <- PsaWrapperTestNum
#    PsaFinalTestNum <- nrow(PseudoabsenceDat.df) - PsaWrapperTrainNum - PsaFinalTrainNum - PsaWrapperTestNum
#    kfoldgrpal <- split(1:nrow(PseudoabsenceDat.df), sample(rep(1:4, c(PsaWrapperTrainNum, PsaWrapperTestNum, PsaFinalTrainNum, PsaFinalTestNum))))
#    kfoldgrpa <- rep(1,nrow(PseudoabsenceDat.df))
#    kfoldgrpa[kfoldgrpal[[1]]]=1
#    kfoldgrpa[kfoldgrpal[[2]]]=2
#    kfoldgrpa[kfoldgrpal[[3]]]=3
#    kfoldgrpa[kfoldgrpal[[4]]]=4
#    occurrences <- table(unlist(kfoldgrpa))
#    ### Save kfold partition scheme for later testing
#    setwd(InDirect)
#    write.table(data.frame(kfoldgrpp), file=paste0(Species, "PresenceDatRSFSAKfold_Rep", Rep, ".csv"), sep=",")
#    write.table(data.frame(kfoldgrpa), file=paste0(Species, PsAbsBuffs, "PseudoabsenceDatRSFSAKfold_Rep", Rep, ".csv"), sep=",")
#    ##
    # Read in k-fold partition scheme
    setwd(InDirect)
    kfoldgrpp <- as.vector(unlist(read.csv(paste0(Species, "PresenceDatRSFSAKfold_Rep", Rep, ".csv"))))
    length(kfoldgrpp)
    kfoldgrpa <- as.vector(unlist(read.csv(paste0(Species, PsAbsBuffs, "PseudoabsenceDatRSFSAKfold_Rep", Rep, ".csv"))))
    if(nrow(PseudoabsenceDat.df)>MaxBackgrndPseudoabsLimit) {
      NObskfoldgrpa <- NObskfoldgrpa[1:MaxBackgrndPseudoabsLimit ]
    }
    length(kfoldgrpa)
    ################################################
    ## Evaluate random subsets for Sets number of subsets
    # Takes 4.48 minutes for 1000 runs
    #ncol(combn(19,5))
    #  t1a <- Sys.time()
    #
    MaxentPresTrainData <- PresenceDat.df[kfoldgrpp==1,]
    head(MaxentPresTrainData)
    MaxentAbsTrainData <- BackgroundDat.df
    head(MaxentAbsTrainData)
    WrapperTrainSWD <- rbind(MaxentPresTrainData[,4:ncol(MaxentPresTrainData)], MaxentAbsTrainData[,4:ncol(MaxentAbsTrainData)])
    head(WrapperTrainSWD)
    tail(WrapperTrainSWD)
    WrapperTrainPresID <- data.frame(rep(1,nrow(MaxentPresTrainData)))
    colnames(WrapperTrainPresID) <- "ID"
    WrapperTrainAbsID <- data.frame(rep(0,nrow(MaxentAbsTrainData)))
    colnames(WrapperTrainAbsID) <- "ID"
    WrapperTrainPresAbsID <- rbind(WrapperTrainPresID, WrapperTrainAbsID)
    head(WrapperTrainPresAbsID)
    tail(WrapperTrainPresAbsID)
    ##############
    ## Partition VariableSubsets according to Rep
    RepSubsets <- nrow(VariableSubsets)
    RepSubsetMinList <- c(1, (RepSubsets/3)+1, ((RepSubsets/3)*2)+1)
    RepSubsetMaxList <- c((RepSubsets/3), ((RepSubsets/3)*2), RepSubsets)
    RepSubsetMin <- RepSubsetMinList[Rep]
    RepSubsetMax <- RepSubsetMaxList[Rep]
    VariableSubsetsRep <- as.data.frame(VariableSubsets[RepSubsetMin:RepSubsetMax,], stringsAsFactors=FALSE)
    colnames(VariableSubsetsRep) <- "VarNames"
    MaxentWrapperTrainTestSubsetEvalStatsL <- list()
    SetCount <- 0
    Sets <- nrow(VariableSubsetsRep)
    #
    output2 <- paste0(NumberModSets,"Setsof", Subset, SubsetVariableNumber, "of", TotVars, "Vars", "_Rep", Rep)
    OutDirect2 <- paste0(OutDirectpsa, "/", output2)
    dir.create(OutDirect2)
    setwd(OutDirect2)
    OutDirectIn <- OutDirect2
    #
    ### Specify three equal RunSize groups of VariableSubsetsRep of no more than around 5,050
    RunSize <- Round2(Sets/3,0)
    NumberModSetsIn <- Sets
    SetRunsList <- list(seq(1,RunSize,1), seq(RunSize+1,RunSize*2,1), seq(RunSize*2+1,Sets,1)) # Groups of VariableSubsets
    ### NOTE: Can comment this loop if just re-ranking wrapper models by different statistic such as AICc
#    ####################
    for(SetRun in SetRunsList) {
      #SetRun <- SetRunsList[[1]]
      VariableSubsetsIn <- as.data.frame(VariableSubsetsRep[SetRun,], stringsAsFactors=FALSE)
      colnames(VariableSubsetsIn) <- "VarNames"
      SetCount <- SetCount + 1
      SetRunID <- SetCount
      system.time(MaxentWrapperTrainTestSubsetEvalStats.df1 <- MaxentMultiSubset.WrapperTrainTestEvalAIC(Species, VariableNamesIn, SubsetVariableNumber, TotPres, VariableSubsetsIn, SetRunID, WrapperTrainSWD, WrapperTrainPresAbsID, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrpp, kfoldgrpa, OutDirectIn, MaxentArgsIn=MaxentBaseArgsIn))
    }
    ######################
  #} # For ending loop here
  #t2 <- Sys.time()
  ######################################
  #difftime(t2,t1, units = "mins")
    # Read back in data
    ### NOTE: The next two lines are for running the loop with reranking and testing only
#    OutDirect2Orig <- gsub("AIC", "AUC", OutDirect2)
#    setwd(OutDirect2Orig)
    setwd(OutDirect2)
    SetRunIDList <- c(1:3)
    #SetRunIDList <- c(1)
    #RowMultList <- c(2)
    RowCount <- 0
    MaxentWrapperTrainTestSubsetEvalStatsL <- list()
    for(SetRunID in SetRunIDList) {
      #SetRunID<-1
      #NumberSets <- length(SetRunsList[[SetRunID]])
      # For re-reading in data
      NumberSets <- length(SetRunsList[[SetRunID]])
      MaxentWrapperTrainTestSubsetEvalStats.df1 <-  as.data.frame(read.csv(paste0(Species, "MaxentResults_WrapperTrainTest_", TotVars, "TotVars_", SubsetVariableNumber, "Vars_", NumberSets, "_", SetRunID, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
      tail(MaxentWrapperTrainTestSubsetEvalStats.df1)
      # Convert first, second and fourth columns from factor to character
      MaxentWrapperTrainTestSubsetEvalStats.df1[,c(1,2,4)] <- sapply(MaxentWrapperTrainTestSubsetEvalStats.df1[,c(1,2,4)], function(x) as.character(x))
      Rows <- nrow(MaxentWrapperTrainTestSubsetEvalStats.df1)
      #str(MaxentWrapperTrainTestSubsetEvalStats.df1)
      #Renumber rows
      if(RowCount==0) {
        RowCount <- Rows
      } else {
        RowCount <- RowCount + Rows
        row.names(MaxentWrapperTrainTestSubsetEvalStats.df1) <- seq(RowCount-Rows+1, RowCount, 1)
      }
      MaxentWrapperTrainTestSubsetEvalStatsL[[SetRunID]] <- MaxentWrapperTrainTestSubsetEvalStats.df1
    }
    MaxentWrapperTrainTestSubsetEvalStats.df <- do.call(rbind, MaxentWrapperTrainTestSubsetEvalStatsL)
    head(MaxentWrapperTrainTestSubsetEvalStats.df)
    ## Separate out training and testing data
    MaxentTrainSubsetEvalStats.df <- MaxentWrapperTrainTestSubsetEvalStats.df[which(MaxentWrapperTrainTestSubsetEvalStats.df$DataType=="WrapperTrain"),]
    MaxentTestSubsetEvalStats.df <- MaxentWrapperTrainTestSubsetEvalStats.df[which(MaxentWrapperTrainTestSubsetEvalStats.df$DataType=="WrapperTest"),]
    # Calculate difference between test and train statistics for overfitting
    MaxentDiffTestTrainSubsetEvalStats.df <- MaxentTrainSubsetEvalStats.df
    MaxentDiffTestTrainSubsetEvalStats.df[,5:13] <- MaxentTrainSubsetEvalStats.df[,5:13] - MaxentTestSubsetEvalStats.df[,5:13]
    MaxentDiffTestTrainSubsetEvalStats.df$DataType <- "WrapperDiffTestTrain"
    #
    # Join AUCDiff with each data set
    MaxentDiffTestTrainSubsetEvalStats.df$AUCdiff <- MaxentDiffTestTrainSubsetEvalStats.df$AUC
    MaxentTrainSubsetEvalStats.df$AUCdiff <- MaxentDiffTestTrainSubsetEvalStats.df$AUC
    MaxentTestSubsetEvalStats.df$AUCdiff <- MaxentDiffTestTrainSubsetEvalStats.df$AUC
    #
    MaxentWrapperTrainTestSubsetEvalStats.df <- do.call(rbind, list(MaxentTrainSubsetEvalStats.df, MaxentTestSubsetEvalStats.df, MaxentDiffTestTrainSubsetEvalStats.df))
    ######################
    ## Loop through different numbers of total subsets from 250 to 500 to 1,000 for ranked and random subsets
    ChosenModSetSizes <- c(Sets)
    #ChosenModSetSizes <- c(1000, 5000, 10000)
    for(ChosenModSets in ChosenModSetSizes) {
      #ChosenModSets=Sets
      ## Rank top NObs subsets by Wrapper TSS and obtain training and final testing TSS for evaluation of each set
      ## and for ensemble set
      ##
      ## NOTE: The next two lines are for running the loop with reranking and testing only
      #OutDirect2Orig <- gsub("AIC", "AUC", OutDirect2)
      #setwd(OutDirect2Orig)
      #setwd(OutDirect2)
      OutCount <- OutCount + 1
      Run <- 1
      #length(unique(MaxentFiltSubsetEvalStats.df$VarNames)) # check for duplicate VarNames
      ## Separate out WrapperTest from data
      MaxentFiltSubsetEvalStats.df <- MaxentWrapperTrainTestSubsetEvalStats.df[which(MaxentWrapperTrainTestSubsetEvalStats.df$DataType=="WrapperTest"),]
      # Delete any subsets with "Inf" for AICc_bg (means too many parameters in model relative to observations)
      MaxentFiltSubsetEvalStats.df2 <- subset(MaxentFiltSubsetEvalStats.df, AICc_bg!="Inf")
      NewSets1 <- nrow(MaxentFiltSubsetEvalStats.df2)
      tail(MaxentFiltSubsetEvalStats.df2)
      ## Keep only first ChosenModSets number of sets
      MaxentFiltSubsetEvalStats.df3 <- MaxentFiltSubsetEvalStats.df2[1:ChosenModSets,]
      head(MaxentFiltSubsetEvalStats.df3)
      tail(MaxentFiltSubsetEvalStats.df3)
      nrow(MaxentFiltSubsetEvalStats.df3)
      ## Use Multi-Objective Optimization from MCDM package to sort by two criteria of AUC and AICc_bg
      decision.mat <- as.matrix(MaxentFiltSubsetEvalStats.df3[,c(11,7,14)])
      head(decision.mat)
      # Normalize AICc_bg and AUC from zero to 1
      normalizef <- function(x) c((x+(max(x)-min(x))-max(x))/(max(x)-min(x)))
      decision.matn <- as.matrix(apply(decision.mat,2,normalizef))
      head(decision.matn)
      ######## Make rank schemes
      RankStatistic <- "AUC"
      AnalysisType <- paste0("RSFSA", RankStatistic)
      ModelName <- paste0("Maxent", DataSet, "_", AnalysisType)
      Model <- paste0(ModelName, TotVars, "Feature Subset")
      ##
      AICcrankList <- c(0)
      AUCrankList <- c(1)
      AUCdiffrankList <- c(0)
      MaxentTopNObsFiltSubsetEvalStatsL <- list()
      TopNObsSubsetsL <- list()
      for(m in 1:length(AUCrankList)) {
        # Specify ranks for AICc and AUC in Multi-Object Optimization ranking
        AICcrank <- AICcrankList[m]
        AUCrank <- AUCrankList[m]
        AUCdiffrank <- AUCdiffrankList[m]
        ##########
        # Assign weights to criteria
        weightscrit <- c(AICcrank, AUCrank, AUCdiffrank)
        # Assign whether cost "min", or benefit "max"
        cb <- c("min", "max", "min")
        # Rank criteria
        MMOORA.Rank <- MMOORA(decision.matn, weightscrit, cb)
        head(MMOORA.Rank)
        MaxentFiltSubsetEvalStats.df4 <- MaxentFiltSubsetEvalStats.df3
        # Join MMOORA rank to original data
        MaxentFiltSubsetEvalStats.df4$MMOORA_Rank <- MMOORA.Rank[,8]
        head(MaxentFiltSubsetEvalStats.df4)
        # Sort data by MMOORA rank
        MaxentFiltSubsetEvalStats.df4s <- MaxentFiltSubsetEvalStats.df4[order(MaxentFiltSubsetEvalStats.df4$MMOORA_Rank),]
        head(MaxentFiltSubsetEvalStats.df4s)
        tail(MaxentFiltSubsetEvalStats.df4s)
        nrow(MaxentFiltSubsetEvalStats.df4s)
        ncol(MaxentFiltSubsetEvalStats.df4s)
        # Delete rank column
        MaxentFiltSubsetEvalStats.df5 <- MaxentFiltSubsetEvalStats.df4s[,-15]
        head(MaxentFiltSubsetEvalStats.df5)
        # Keep top FullNObs selected subsets
        MaxentTopNObsFiltSubsetEvalStatsL[[m]] <- data.frame(MaxentFiltSubsetEvalStats.df5[1:FullNObs,], stringsAsFactors=FALSE)
        TopNObsSubsetsL[[m]] <- data.frame(MaxentFiltSubsetEvalStats.df5[1:FullNObs,1], stringsAsFactors=FALSE)
      }
      #######################
      MaxentTopNObsFiltSubsetEvalStats.df <- do.call(rbind, MaxentTopNObsFiltSubsetEvalStatsL)
      nrow(MaxentTopNObsFiltSubsetEvalStats.df)
      TopNObsSubsets <- do.call(rbind, TopNObsSubsetsL)
      ##############################
      # Keep random FullNObs plus PlusNum random subsets
      MaxentRandNObsFiltSubsetEvalStats.df1 <- data.frame(MaxentFiltSubsetEvalStats.df3[sample(nrow(MaxentFiltSubsetEvalStats.df3), FullNObs+PlusNum), ], stringsAsFactors=FALSE)
      ## Clean random data from aberrant AICc_bg values
      # Delete any subsets with "Inf" for AICc_bg (means too many parameters in model relative to observations)
      MaxentRandNObsFiltSubsetEvalStats.df2 <- subset(MaxentRandNObsFiltSubsetEvalStats.df1, AICc_bg!="Inf")
      # Delete any subsets with negative AICc_bg which appear to be spurious values
      MaxentRandNObsFiltSubsetEvalStats.df3 <- subset(MaxentRandNObsFiltSubsetEvalStats.df2, AICc_bg>-1)
      # Keep only FullNObs plus PlusNum/2, or around 50 sets
      MaxentRandNObsFiltSubsetEvalStats.df <- MaxentRandNObsFiltSubsetEvalStats.df3[1:(FullNObs+Round2(PlusNum/2,0)),]
      tail(MaxentRandNObsFiltSubsetEvalStats.df)
      nrow(MaxentRandNObsFiltSubsetEvalStats.df)
      ncol(MaxentRandNObsFiltSubsetEvalStats.df)
      #
      RandNObsSubsets <- data.frame(MaxentRandNObsFiltSubsetEvalStats.df[,1], stringsAsFactors=FALSE)
      colnames(RandNObsSubsets) <- "VarNames"
      # Use difference between evaluation statistics between train and test data to evaluate overfitting
      ###############################################################################
      ## Test Top FullNObs Ranking Subsets and Random FullNObs Subsets
      ###############################################################################
      #####################################
      setwd(OutDirect2)
      #SubsetRuns <- list(TopNObsSubsets[1:NObs,], TopNObsSubsets[(NObs+1):(NObs*2),], TopNObsSubsets[(NObs*2+1):(NObs*3),], TopNObsSubsets[(NObs*3+1):(NObs*4),],RandNObsSubsets[1:NObs,], RandNObsSubsets[(NObs+1):(NObs*2),], RandNObsSubsets[(NObs*2+1):(NObs*3),])
      #DataSetTypes <- c(paste0("TopAUC",NObs), paste0("TopAICc",NObs),  paste0("TopAUCdiff",NObs), paste0("TopAUC_AICc_AUCdiff",NObs), paste0("Random", NObs, "A"), paste0("Random", NObs, "B"), paste0("Random", NObs, "C"))
      SubsetRuns <- list(TopNObsSubsets, RandNObsSubsets)
      DataSetTypes <- c(paste0("TopAUC",FullNObs), paste0("Random", FullNObs))
      MaxentKeepEvalStatsL1 <- list()
      #
      for (j in 1:length(DataSetTypes)) {
        #N3Subsets <- SubsetRuns[[1]]
        #j=1
        VariableSubsetsIn <- data.frame(SubsetRuns[[j]], stringsAsFactors=FALSE)
        DataSetType <- DataSetTypes[j]
        # Identify Wrapper data evaluation statistics from above in loop
        if(DataSetType==paste0("TopAUC", FullNObs)) {
          MaxentTestSubsetEvalStats.df <- MaxentTopNObsFiltSubsetEvalStats.df
        } else {
          MaxentTestSubsetEvalStats.df <- MaxentRandNObsFiltSubsetEvalStats.df
        }
        MaxentTestSubsetEvalStats.df$DataType <- paste0(DataSetType, "WrapperTest")
        ##
        ## Join Wrapper Test data back with Wrapper Train and WrapperDiffTestTrain
        ## Subset original data using VarNames
        selectedRows <- (MaxentWrapperTrainTestSubsetEvalStats.df$VarNames %in% MaxentTestSubsetEvalStats.df$VarNames)
        MaxentTestSubsetEvalStats.df2 <- MaxentWrapperTrainTestSubsetEvalStats.df[selectedRows,]
        ## Take apart by DataType, reorder each by rank and reassemble
        MaxentTrainSubsetEvalStats.df3 <- MaxentTestSubsetEvalStats.df2[which(MaxentTestSubsetEvalStats.df2$DataType=="WrapperTrain"),]
        MaxentTrainSubsetEvalStats.df4 <- MaxentTrainSubsetEvalStats.df3[match(unlist(MaxentTestSubsetEvalStats.df$VarNames), MaxentTrainSubsetEvalStats.df3$VarNames), ]
        #
        MaxentTestSubsetEvalStats.df3 <- MaxentTestSubsetEvalStats.df2[which(MaxentTestSubsetEvalStats.df2$DataType=="WrapperTest"),]
        MaxentTestSubsetEvalStats.df4 <- MaxentTestSubsetEvalStats.df3[match(unlist(MaxentTestSubsetEvalStats.df$VarNames), MaxentTestSubsetEvalStats.df3$VarNames), ]
        #
        MaxentDiffTestTrainSubsetEvalStats.df3 <- MaxentTestSubsetEvalStats.df2[which(MaxentTestSubsetEvalStats.df2$DataType=="WrapperDiffTestTrain"),]
        MaxentDiffTestTrainSubsetEvalStats.df4 <- MaxentDiffTestTrainSubsetEvalStats.df3[match(unlist(MaxentTestSubsetEvalStats.df$VarNames), MaxentDiffTestTrainSubsetEvalStats.df3$VarNames), ]
        # Reassemble
        MaxentSubsetEvalStats.df <- rbind(MaxentTrainSubsetEvalStats.df4, MaxentTestSubsetEvalStats.df4, MaxentDiffTestTrainSubsetEvalStats.df4)
        MaxentSubsetEvalStats.df$DataType <- paste0(DataSetType, MaxentSubsetEvalStats.df$DataType)
        MaxentSubsetEvalStats.df$SetName <- SetName
        MaxentSubsetEvalStats.df$Set <- DataSetType
        MaxentSubsetEvalStats.df$TotSubsets <- ChosenModSets
        MaxentSubsetEvalStats.df$PresTrainPoints <- nrow(PresenceDat.df[kfoldgrpp==1,])
        MaxentSubsetEvalStats.df$Rep <- Rep
        ncol(MaxentSubsetEvalStats.df)
        MaxentSubsetEvalStats.df <- MaxentSubsetEvalStats.df[, colnames(MaxentSubsetEvalStats.df)[c(1:4,15:19,5:14)]]
        #
        MaxentKeepEvalStatsL1[[j]] <- MaxentSubsetEvalStats.df
      }
      ###############
      ## Designate Final Training data for Final Test
      MaxentPresTrainData <- PresenceDat.df[kfoldgrpp==3,]
      head(MaxentPresTrainData)
      nrow(MaxentPresTrainData)
      MaxentAbsTrainData <- BackgroundDat.df
      nrow(MaxentAbsTrainData)
      head(MaxentAbsTrainData)
      FinalTrainSWD <- rbind(MaxentPresTrainData[,4:ncol(MaxentPresTrainData)], MaxentAbsTrainData[,4:ncol(MaxentAbsTrainData)])
      head(FinalTrainSWD)
      tail(FinalTrainSWD)
      FinalTrainPresID <- data.frame(rep(1,nrow(MaxentPresTrainData)))
      colnames(FinalTrainPresID) <- "ID"
      FinalTrainAbsID <- data.frame(rep(0,nrow(MaxentAbsTrainData)))
      colnames(FinalTrainAbsID) <- "ID"
      FinalTrainPresAbsID <- rbind(FinalTrainPresID, FinalTrainAbsID)
      head(FinalTrainPresAbsID)
      tail(FinalTrainPresAbsID)
      #
      setwd(OutDirect2)
      #SubsetRuns <- list(TopNObsSubsets[1:NObs,], TopNObsSubsets[(NObs+1):(NObs*2),], TopNObsSubsets[(NObs*2+1):(NObs*3),], TopNObsSubsets[(NObs*3+1):(NObs*4),],RandNObsSubsets[1:NObs,], RandNObsSubsets[(NObs+1):(NObs*2),], RandNObsSubsets[(NObs*2+1):(NObs*3),])
      #DataSetTypes <- c(paste0("TopAUC",NObs), paste0("TopAICc",NObs),  paste0("TopAUCdiff",NObs), paste0("TopAUC_AICc_AUCdiff",NObs), paste0("Random", NObs, "A"), paste0("Random", NObs, "B"), paste0("Random", NObs, "C"))
      SubsetRuns <- list(TopNObsSubsets, RandNObsSubsets)
      DataSetTypes <- c(paste0("TopAUC",FullNObs), paste0("Random", FullNObs))
      MaxentKeepEvalStatsL2 <- list()
      #
      for (j in 1:length(DataSetTypes)) {
        #j=1
        VariableSubsetsIn <- data.frame(SubsetRuns[[j]], stringsAsFactors=FALSE)
        DataSetType <- DataSetTypes[j]
        ## Evaluate random subsets for NumberModSets number of subsets
        system.time(MaxentTrainTestSubsetEvalStats.df1 <- MaxentMultiSubset.FinalTrainTestEvalAIC(Species, VariableNamesIn, SubsetVariableNumber, TotPres, VariableSubsetsIn, FinalTrainSWD, FinalTrainPresAbsID, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrpp, kfoldgrpa, OutDirectIn, MaxentArgsIn=MaxentBaseArgsIn, DataSetType=DataSetType))
        # Delete any subsets with "Inf" for AICc_bg (means too many parameters in model relative to observations)
        MaxentTrainTestSubsetEvalStats.df2 <- subset(MaxentTrainTestSubsetEvalStats.df1, AICc_bg!="Inf")
        nrow(MaxentTrainTestSubsetEvalStats.df2)
        # Delete any subsets with negative values of AICc_bg which appear to be spurious values (equal to or less than -1)
        MaxentTrainTestSubsetEvalStats.df <- subset(MaxentTrainTestSubsetEvalStats.df2, AICc_bg>-1)
        nrow(MaxentTrainTestSubsetEvalStats.df)
        ## Separate out training and testing data
        MaxentTrainSubsetEvalStats.df <- MaxentTrainTestSubsetEvalStats.df[which(MaxentTrainTestSubsetEvalStats.df$DataType=="FinalTrain"),]
        MaxentTrainSubsetEvalStats.df$DataType <- paste0(DataSetType, "FinalTrain")
        MaxentTestSubsetEvalStats.df <- MaxentTrainTestSubsetEvalStats.df[which(MaxentTrainTestSubsetEvalStats.df$DataType=="FinalTest"),]
        MaxentTestSubsetEvalStats.df$DataType <- paste0(DataSetType, "FinalTest")
        # Calculate difference between test and train statistics for overfitting
        MaxentDiffTestTrainSubsetEvalStats.df <- MaxentTestSubsetEvalStats.df
        MaxentDiffTestTrainSubsetEvalStats.df[,5:13] <- MaxentTrainSubsetEvalStats.df[,5:13] - MaxentTestSubsetEvalStats.df[,5:13]
        MaxentDiffTestTrainSubsetEvalStats.df$DataType <- paste0(DataSetType, "FinalDiffTestTrain")
        #
        # Join AUCDiff with each data set
        MaxentDiffTestTrainSubsetEvalStats.df$AUCdiff <- MaxentDiffTestTrainSubsetEvalStats.df$AUC
        MaxentTrainSubsetEvalStats.df$AUCdiff <- MaxentDiffTestTrainSubsetEvalStats.df$AUC
        MaxentTestSubsetEvalStats.df$AUCdiff <- MaxentDiffTestTrainSubsetEvalStats.df$AUC
        #
        MaxentSubsetEvalStats.df <- do.call(rbind, list(MaxentTrainSubsetEvalStats.df, MaxentTestSubsetEvalStats.df, MaxentDiffTestTrainSubsetEvalStats.df))
        MaxentSubsetEvalStats.df$SetName <- SetName
        MaxentSubsetEvalStats.df$Set <- DataSetType
        MaxentSubsetEvalStats.df$TotSubsets <- ChosenModSets
        MaxentSubsetEvalStats.df$PresTrainPoints <- nrow(PresenceDat.df[kfoldgrpp==3,])
        MaxentSubsetEvalStats.df$Rep <- Rep
        MaxentSubsetEvalStats.df <- MaxentSubsetEvalStats.df[, colnames(MaxentSubsetEvalStats.df)[c(1:4,15:19,5:14)]]
        #
        MaxentKeepEvalStatsL2[[j]] <- MaxentSubsetEvalStats.df
      }
      ###############
      MaxentKeepEvalStats.df1 <- do.call(rbind, MaxentKeepEvalStatsL1)
      nrow(MaxentKeepEvalStats.df1)
      MaxentKeepEvalStats.df2 <- do.call(rbind, MaxentKeepEvalStatsL2)
      nrow(MaxentKeepEvalStats.df2)
      MaxentKeepEvalStats.df3 <- rbind(MaxentKeepEvalStats.df1, MaxentKeepEvalStats.df2)
      nrow(MaxentKeepEvalStats.df3)
      MaxentKeepEvalAllStatsL[[OutCount]] <- MaxentKeepEvalStats.df3
      ### Summarize Output
      ###Calculate Mean and Standard Deviation Values
      setwd(OutDirectpsa)
      RunType <- "Final"
      EvaluationStats <- MaxentKeepEvalStats.df3
      ncol(EvaluationStats)
      EvaluationStats$Model <- ModelName
      EvaluationStats <- EvaluationStats[, colnames(EvaluationStats)[c(20,1:19)]]
      #str(EvaluationStats)
      #FSAType <- RankStatistic
      OutName <- paste0(Species, ModelName, FullNObs, "_", NumberModSets,"SetNumber_SubsetSummaryStats_", SubsetVariableNumber, "of", TotVars, "Vars", FSAType,"_", RunType, ".csv")
      SortGroups <- c("Model", "Rep", "DataType", "Set", "SetName", "SubsetVariableNumber", "TotSubsets")
      ## All statistics, and only statistics, should be at and after the column specified below
      StatVarFirstColumn <- 9
      OutDirectIn <- OutDirectpsa
      EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
      head(EvaluationStatsSummary.df)
      MaxentKeepEvalAllStatSummL[[OutCount]] <- EvaluationStatsSummary.df
      ######
      MaxentKeepEvalStatsAllSets.df <- do.call(rbind, MaxentKeepEvalAllStatSummL)
      ncol(MaxentKeepEvalStatsAllSets.df)
      ## Save data
      write.table(MaxentKeepEvalStatsAllSets.df, file=paste0(Species, ModelName, FullNObs, "_", NumberModSets,"SetNumber_", "SummaryStats_of", FinalModelVariables, "of", TotVars, "Vars", FSAType,"_", RunType, ".csv"), sep=",", col.names=NA)
      ### Summarize and save raw data as loop progresses
      MaxentKeepEvalStats.df <- do.call(rbind, MaxentKeepEvalAllStatsL)
      # Save data
      write.table(MaxentKeepEvalStats.df , file=paste0(Species, ModelName, FullNObs, "_", NumberModSets,"SetNumber_TrainVsTestStats_of", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff", ".csv"), sep=",", col.names=NA)
      #
     }
  }
}
#
####################################
t2 <- Sys.time()
####################################
difftime(t2,t1, units = "mins")
setwd(OutDirectpsa)
# Read data back in and sort
MaxentKeepEvalStats.df1 <- data.frame(read.csv(paste0(Species, ModelName, FullNObs, "_", NumberModSets,"SetNumber_TrainVsTestStats_of", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
nrow(MaxentKeepEvalStats.df1)
ncol(MaxentKeepEvalStats.df1)
head(MaxentKeepEvalStats.df1)
tail(MaxentKeepEvalStats.df1)
#
#EvaluationStats <- MaxentKeepEvalStats.df1
#   write.table(EvaluationStatsSummary.df, file=paste0(Species, ModelName, "_", NumberModSets,"SetNumber_", "SummaryStats_of", FinalModelVariables, "of", TotVars, "Vars", FSAType,"_", RunType, ".csv"), sep=",", col.names=NA)

######################################################################################
## Assemble data for bar charts and Welch t test
##########################################################################################
setwd(OutDirectpsa)
# Read in results for first FSA run to get output on the leave one out variable run with all variables
#FSAType <- RankStatistic
#NumberModSets <- 10000
#
# Keep only DataTypes including string "FinalTest"
MaxentKeepEvalStats.df2 <- MaxentKeepEvalStats.df1[grep("FinalTest", MaxentKeepEvalStats.df1$DataType), ]
MaxentKeepEvalStats.df2[1:50,]
head(MaxentKeepEvalStats.df2)
tail(MaxentKeepEvalStats.df2)
###Calculate Mean and Standard Deviation Values
RunType <- "Final"
EvaluationStats <- MaxentKeepEvalStats.df2
OutName <- paste0(Species, RankStatistic, "Maxent_", NumberModSets, "Setsof_", SubsetVariableNumber, "_", FullNObs, "N_StatSummary.csv")
SortGroups <- c("Rep", "DataType", "TotSubsets", "SubsetVariableNumber")
## All statistics, and only statistics, should be at and after the column specified below
StatVarFirstColumn <- 6
OutDirectIn <- OutDirectpsa
EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
head(EvaluationStatsSummary.df)
##
EvalStatsOut <- EvaluationStatsSummary.df
## Sort data by DataType for ease of plotting in excel
EvalStatsOuts <- arrange(EvalStatsOut, Rep, DataType, Statistic, DataType)
EvalStatsOuts <- arrange(EvalStatsOut, Statistic)
head(EvalStatsOuts)
write.table(EvalStatsOuts, file=paste0(Species, RankStatistic, "Maxent_", SetName, "_", FullNObs, "Stat", SubsetVariableNumber, "VariableRunsby", NumberModSets, "SetsSummaryTable.csv"), sep=",", col.names=NA)
##
#
#############################################
### Conduct statistical tests for Rep
##################################################################################
# Loop through various Welch t test comparisons for each ranking scheme and the three random sets
StatTypeNames <- c("AUC", "AICc_bg", "AUCdiff")
RankTypeNames <- c(paste0("TopAUC", FullNObs))
RandomTypeNames <- c(paste0("Random", FullNObs))
#DataType <- "Test"
##
for(Rep in RepList) {
  # Keep only data for Rep
  #Rep = 1
  RepNum <- Rep
  MaxentKeepEvalStats.df3 <- subset(MaxentKeepEvalStats.df2, Rep==RepNum)
  head(MaxentKeepEvalStats.df3)
  for(StatType in StatTypeNames) {
    #StatType <- "AUC"
    for(RankType in RankTypeNames) {
      #RankType <- paste0("TopAUC", FullNObs)
      MaxentEvalStatsL <- list()
      Count <- 0
      for(RandomType in RandomTypeNames) {
        #RandomType <- paste0("Random", FullNObs)
        ## Keep only DataType with RankType  or RandomType
        MaxentTestStats <- subset(MaxentKeepEvalStats.df3, Set %in% c(RankType, RandomType))
        Count <- Count + 1
        if(StatType==StatTypeNames[1]) {
          # Identify set as Top or Random
          MaxentTestStats$SetType <- "Top"
          ncol(MaxentTestStats)
          head(MaxentTestStats)
          MaxentTestStats <- MaxentTestStats[,c(1:6, 20, 7:19)]
          MaxentTestStats <- within(MaxentTestStats, SetType[grepl("Rand", DataType)] <- 'Random')
          # Perform summary stats for excel plotting
          ###Calculate Mean and Standard Deviation Values
          RunType <- "Final"
          EvaluationStats <- MaxentTestStats
          OutName <- paste0(Species, NumberModSets, SetNameF, "Setsof", SubsetVariableNumber, "Variables", RankType, "Vs", RandomType, "_", StatType, FullNObs, "N_StatSummary", "_Rep", Rep, ".csv")
          SortGroups <- c("DataType", "SetName", "SetType", "TotSubsets", "SubsetVariableNumber")
          ## All statistics, and only statistics, should be at and after the column specified below
          StatVarFirstColumn <- 8
          OutDirectIn <- OutDirectpsa
          EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
          head(EvaluationStatsSummary.df)
          ##
          EvalStatsOut <- EvaluationStatsSummary.df
          ## Sort data by DataType for ease of plotting in excel
          EvalStatsOuts <- arrange(EvalStatsOut, DataType, Statistic, DataType)
          EvalStatsOuts <- arrange(EvalStatsOut, Statistic)
          head(EvalStatsOuts)
          #
          MaxentEvalStatsL[[Count]] <- EvalStatsOuts
        }
        #
        setwd(OutDirectpsa)
        ## Use Welch correction t test
        RankWelchttest <- t.test(as.formula(paste0(StatType, " ~ DataType")), data=MaxentTestStats)
        #str(RankWelchttest)
        ## Apply Bonferroni Correction by mulitplying the p-value by the number of t-tests per StatType
        PValueBonfCor <- length(RandomTypeNames)*RankWelchttest$p.value
        # Save above Welch t test output
        out1 <- paste0(RankType, " vs ", RandomType, ": Welch t test for ", StatType, " from ", NumberModSets, SetNameF, " Sets of ", SubsetVariableNumber, " Variables")
        cat(out1, file=paste0(Species, "WelchttestStatsby", FullNObs, "from", NumberModSets, SetNameF, "Setsof", SubsetVariableNumber, "Variables", RankType, "Vs", RandomType, "_", StatType, "_Rep", Rep, ".txt"), sep="\n", append=TRUE)
        out2 <- capture.output(print(RankWelchttest))
        cat(out2, file=paste0(Species, "WelchttestStatsby", FullNObs, "from", NumberModSets, SetNameF, "Setsof", SubsetVariableNumber, "Variables", RankType, "Vs", RandomType, "_", StatType, "_Rep", Rep, ".txt"), sep="\n", append=TRUE)
        out3 <- paste0("p-value with Bonferroni Correction for ", length(RandomTypeNames), " Welch t tests: ", PValueBonfCor)
        cat(out3, file=paste0(Species, "WelchttestStatsby", FullNObs, "from", NumberModSets, SetNameF, "Setsof", SubsetVariableNumber, "Variables", RankType, "Vs", RandomType, "_", StatType, "_Rep", Rep, ".txt"), sep="\n", append=TRUE)
      }
      if(StatType==StatTypeNames[1]) {
        MaxentEvalStatsAll.df <- do.call(rbind, MaxentEvalStatsL)
        ## Sort data
        MaxentEvalStatsAll.dfs <- arrange(MaxentEvalStatsAll.df, SetType, Statistic, SetType)
        write.table(MaxentEvalStatsAll.dfs, file=paste0(Species, "WelchttestStatsby", FullNObs, "from", NumberModSets, SetNameF, "Setsof", SubsetVariableNumber, "Vars", RankType, "vsRandom_DataTable", "_Rep", Rep, ".csv"), sep=",", col.names=NA)
    }
    }
  }
}
##########
##############################################################################



#####################################################################################################
### Loop through all top 30 selected subsets of final selected SubsetVariableNumber
### to generate maxent projection maps and rank by regional variables to choose four top selected
### maps for further analysis
#####################################################################################################
## This loop takes about 80 minutes for 10 variable subsets of selected 10 and random 10 models
## This loop takes 2.2 hours for 23 models of 12 variable subsets
## This loop takes 4.1 hours for 45 models of 12 variable subsets
## This loop takes about 9.3 hours for 30 models of 15 variables subsets for fire
## This loop takes about 5.31 hours for 4 models of 6 variables subsets for fire
#
# Specify final model variable
FinalModelVariables <- 15
#FinalModelVariables <- 3
NObsJ <- 3
NumberModSets <- 3000
#NumProjections <- 9  # Specify number of models to project
#NumberModSets <- 30
NumProjections <- 15
## Read in clipped current climate and other grids
setwd(MaskGrids)
Predictors <- stack(GridNamesL)
names(Predictors) <- toupper(names(Predictors))
names(Predictors) <- gsub("_NS", "", names(Predictors))
#
#############################
j = 10
#j=1
if(j < 10) {
  SetNumber <- paste0("0", j)
} else {
  SetNumber <- as.character(j)
}
#
SetName <- paste0(SetNumber, "_", colnames(VariableSetTypes[j]))
SetNameF <- colnames(VariableSetTypes[j])
setwd(OutDirectpsa)
# Read in results for first FSA run to get output on the leave one out variable run with all variables
#FSAType <- RankStatistic
#NumberModSets <- 33
#SubsetSizes <- c(3,8,10,12,15,20,25)
#
### Assemble top 15 ranked variable subsets by AUCWrappertest from each RSFSA rep
MaxentKeepEvalStats.df1 <- data.frame(read.csv(paste0(Species, ModelName, FullNObs, "_", NumberModSets,"SetNumber_TrainVsTestStats_of", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
nrow(MaxentKeepEvalStats.df1)
ncol(MaxentKeepEvalStats.df1)
head(MaxentKeepEvalStats.df1)
tail(MaxentKeepEvalStats.df1)
### Keep only first NumProjections/3 variable subsets from each Rep
## Keep only FinalTest data
MaxentKeepEvalStats.df2 <- MaxentKeepEvalStats.df1[grep(paste0("FinalTest"), MaxentKeepEvalStats.df1$DataType), ]
nrow(MaxentKeepEvalStats.df2)
head(MaxentKeepEvalStats.df2)
tail(MaxentKeepEvalStats.df2)
ModelTypeNames <- c(paste0("TopAUC",FullNObs), paste0("Random", FullNObs))
MaxentKeepEvalStatsL <- list()
Count <- 0
RepList <- seq(1:3)
##################
for(Rep in RepList) {
  #Rep=1
  for(ModelType in ModelTypeNames) {
    #ModelType <- paste0("TopAUC",FullNObs) 
    # Keep only data for Rep and ModelType
    MaxentTestKeepEvalStats.df3 <- MaxentKeepEvalStats.df2[which(MaxentKeepEvalStats.df2$Rep==Rep & MaxentKeepEvalStats.df2$DataType==paste0(ModelType, "FinalTest")),]       
    # Keep only first NumProjections/3 number of rows
    Count <- Count + 1
    if(NumProjections > 6) {
      MaxentKeepEvalStatsL[[Count]] <- MaxentTestKeepEvalStats.df3[1:(NumProjections/3),]
    } else {
      MaxentKeepEvalStatsL[[Count]] <- MaxentTestKeepEvalStats.df3[1:2, ]
    }
  }
}
#############
MaxentKeepEvalStats.df4 <- do.call(rbind, MaxentKeepEvalStatsL)
## Sort data by DataType
MaxentKeepEvalStats.df4s <- arrange(MaxentKeepEvalStats.df4, DataType)
nrow(MaxentKeepEvalStats.df4s)
## Save data
setwd(OutDirectpsa)
write.table(MaxentKeepEvalStats.df4s, file=paste0(Species, ModelName, FullNObs, "_", NumberModSets, "TopSetsFinalTestStats_of", FinalModelVariables, "of", TotVars, "NumVars_.csv"), sep=",", col.names=NA)
###
## Loop through sets 1 through 10
#ModelTypeNames <- c(paste0("TopAUC",FullNObs), paste0("Random", FullNObs))
#ModelTypeNames <- c(paste0("TopAUC",FullNObs))
ModelTypeNames <- c(paste0("Random", FullNObs))
MaxentKeepEvalStatsL <- list()
Count <- 0
#Count <- 21
t1 <- Sys.time()
#############################################################################
for(ModelType in ModelTypeNames) {
  #ModelType <- paste0("TopAUC",FullNObs)
  # Keep only DataTypes including string of ModelType for Final Test
  MaxentKeepEvalStats.df5 <- MaxentKeepEvalStats.df4s[grep(paste0(ModelType, "FinalTest"), MaxentKeepEvalStats.df4s$DataType), ]
  #
  if(NumProjections < 6) {
    # Sort by RSFSA Wrapper
    MaxentKeepEvalStats.df5 <- arrange(MaxentKeepEvalStats.df5, -MaxentKeepEvalStats.df5[[RankStatistic]])
  }
  #
  for(k in 1:NumProjections) {
  #for(k in 7:NumProjections) {
  #for(k in 1) {
    #k=15
    Loop <- k
    ######################################################################################
    ## Run Threshold Calibration Projection for Top Selected Variable Model for Given Variable Set (j)
    ######################################################################################
    ##
    #### Create and save NObsJ usually three kfold partition scheme for multiple rounds of test with presence
    #### pseudoabsence and background points using all data
    #### NOTE: Only create and save partition scheme once in order to reuse
#    setwd(InDirect)
#    library(dismo)
#    #
#    NObskfoldgrpp <- kfold(PresenceDat.df, NObsJ)
#    NObskfoldgrpa <- kfold(PseudoabsenceDat.df, NObsJ)
#    ## Save kfold partition scheme for later testing
#    setwd(InDirect)
#    write.table(data.frame(NObskfoldgrpp), file=paste0(Species, "PresenceDat", NObsJ, "Kfold", ".csv"), sep=",")
#    write.table(data.frame(NObskfoldgrpa), file=paste0(Species, PsAbsBuffs, "PseudoabsenceDat", NObsJ, "Kfold", ".csv"), sep=",")
    ##
    #############################################################################################
    ##
    #################
    ## NOTE: DO NOT RESORT BY THIS AUCfinaltest- THIS WILL INFLATE AUC- AUCWrappertest already used to sort
    #################
    VariableSubset <- data.frame(as.character(MaxentKeepEvalStats.df5[k,1]), stringsAsFactors=FALSE)
    colnames(VariableSubset) <- "VarNames"
    VarNames <- unlist(strsplit(as.character(VariableSubset[,1]),"-"))
    SubsetVariableNumber <- length(VarNames)
    Rep <- MaxentKeepEvalStats.df5[k,8]
    #
    #VarNames <- VarNames[-14]
    #SubsetVariableNumber <- length(VarNames)
    #VariableSubset <- as.data.frame(paste0(unlist(VarNames), collapse="-"), stringsAsFactors=FALSE)
    #colnames(VariableSubset) <- "VarNames"
    #########
    ##
    ## Use 3-Fold run and use only one fold for presence data to derive threshold (save time): 2/3 train, 1/3 test
    ## But use ALL of the pseudoabsence data for testing (not the fold)
    ######################
    #### FIRST: Find the threshold values for the Maxent model using a training grid
    #### model
    ## Read in three-fold partition scheme
    setwd(InDirect)
    NObskfoldgrpp <- as.vector(unlist(read.csv(paste0(Species, "PresenceDat", NObsJ, "Kfold", ".csv"))))
    length(NObskfoldgrpp)
    NObskfoldgrpa <- as.vector(unlist(read.csv(paste0(Species, PsAbsBuffs, "PseudoabsenceDat", NObsJ, "Kfold", ".csv"))))
    if(nrow(PseudoabsenceDat.df)>MaxBackgrndPseudoabsLimit) {
      NObskfoldgrpa <- NObskfoldgrpa[1:MaxBackgrndPseudoabsLimit ]
    }
    length(NObskfoldgrpa)
    #
    Run=1
    ##
    VariableNamesIn <- VariableNames
    #
    ###########################################
    ##
    ##
    setwd(OutDirectpsa)
    output3 <- paste0("MaxentProj_", ModelType, "_Subsetof", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_ThresholdCalib", Loop)
    #output3 <- paste0("MaxentProj_", ModelType, "_Subsetof", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_ThresholdCalib", Loop, "Minus", 14)
    dir.create(paste0(OutDirectpsa, "/", output3))
    OutDirectSub <- paste0(OutDirectpsa, "/", output3)
    OutDirectIn <- OutDirectSub
    ##
    kfoldgrppin <- NObskfoldgrpp
    kfoldgrpain <- NObskfoldgrpa
    #
  #  # If necessary, read back in r grids
    #setwd(MaskGrids)
  #  Predictors <- stack(GridNamesL)
  #  names(Predictors) <- toupper(names(Predictors))
  #  names(Predictors) <- gsub("_NS", "", names(Predictors))
    PredictorsIn <- Predictors
    #PredictorsIn[[75:90]]
    #plot(PredictorsIn[[84]])
    #writeRaster(PredictorsIn[[86]], paste0("Test86"), format = "GTiff", overwrite=TRUE)
    #
    OutGridID <- paste0(ModelType, "_", Loop, "ThreshCal_", SetNameF)
    # Takes 19 minutes for 25 variable grid
    #
    SetNameIn <- k
    CRS.In <- CRS.WGS84
    #
    Time <- system.time(MaxentKeepEvalStats.df <- MaxentSubset.GridTrainTestEvalAIC_AUCbgp_Calib(Species, VariableNamesIn, PredictorsIn, SubsetVariableNumber, TotPres, Run, SetName, VariableSubset, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrppin, kfoldgrpain, CRS.In, OutDirectIn, MaxentArgsIn=MaxentBaseArgsIn,  Output=TRUE))
    #
    head(MaxentKeepEvalStats.df)
    # Save results
    setwd(OutDirectIn)
    write.table(MaxentKeepEvalStats.df , file=paste0(Species, "Maxent", DataSet, "_TrainVsTestStatsThreshCal_of", NumberModSets, "_", SetName, SubsetVariableNumber, "_", Loop, ".csv"), sep=",", col.names=NA)
    #
    #MaxentMod <- readRDS("MaxentModel.rds") # Reads in stored Maxent model from above function
    Count <- Count + 1
    ## Add ModelType to output
    MaxentKeepEvalStats.df$ModelType <- ModelType
    ncol(MaxentKeepEvalStats.df)
    # Re-order columns
    MaxentKeepEvalStats.df <- MaxentKeepEvalStats.df[,c(1:3, 17, 3:16)]
    # Replace Run with Rep
    MaxentKeepEvalStats.df$Rep <- Rep
    ncol(MaxentKeepEvalStats.df)
    MaxentKeepEvalStats.df <- MaxentKeepEvalStats.df[,c(1:7, 18, 8:17)]
    MaxentKeepEvalStatsL[[Count]] <-  MaxentKeepEvalStats.df
  }
}
t2 <- Sys.time()
#############################################################################
difftime(t2,t1, units = "mins")
###
MaxentKeepEvalStatsAll.df <- do.call(rbind, MaxentKeepEvalStatsL)
## Sort by DataType
MaxentKeepEvalStatsAll.dfs  <- arrange(MaxentKeepEvalStatsAll.df, DataType)
## Save data
setwd(OutDirectpsa)
write.table(MaxentKeepEvalStatsAll.dfs, file=paste0(Species, "Maxent", DataSet, "_TrainVsTestStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", EvalType, "_TopVsRandom", NumProjections, "Sets.csv"), sep=",", col.names=NA)

##################################################################################
### Analyze above data  for linear regression between AICc and AICc_bg
## Read above data back in
# Read data back in and sort
#SubsetVariableNumber <- 8
#NumberModSets <- 10000
setwd(OutDirectpsa)
MaxentKeepEvalGridStats.df <- data.frame(read.csv(paste0(Species, "Maxent", DataSet, "_TrainVsTestStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", EvalType, "_TopVsRandom", NumProjections, "Sets.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
##
## Keep only DataType of Test
MaxentKeepEvalGridStats.df1 <- subset(MaxentKeepEvalGridStats.df, DataType=="Test")
# Test for significant linear regression between AICc and AICcbg
AICc_AICcbg.lm <- lm(AICc ~ AICc_bg, data=MaxentKeepEvalGridStats.df1)
AICc_AICcbg.lmRes <- summary(AICc_AICcbg.lm)
# Save above output
out1 <- paste0("Linear regression F test for AICc Versus AICc_bg from", NumberModSets, SetNameF, "Setsof", SubsetVariableNumber, "Variables")
cat(out1, file=paste0(Species, RankStatistic, "AICcVsAICc_bg_LinearReg_by", NumProjections, "from", NumberModSets, SetNameF, "Setsof", SubsetVariableNumber, "VariablesAICcVsAICc_bg.txt"), sep="\n", append=TRUE)
out3 <- capture.output(print(AICc_AICcbg.lmRes))
cat(out3, file=paste0(Species, RankStatistic, "AICcVsAICc_bg_LinearReg_by", NumProjections, "from", NumberModSets, SetNameF, "Setsof", SubsetVariableNumber, "VariablesAICcVsAICc_bg.txt"), sep="\n", append=TRUE)



####################################################################################
### Use Arc Python script ENMBreedingAreaZonalStatisticsLooopSWFL.py to calculate percent 
### overlap of presence projection inside breeding range with total breeding range area
### and calculate percent overlap of presence projection outside breeding range with total breeding range area
### within NAD_1983_Albers projection
####################################################################################
## Read in above produced csv files from python and calculate percent overlaps of maxent models
## with breeding range for portions of maxent model both inside and outside breeding range
##
FinalModelVariables <- 15
NumProjections <- 30
NObsJ <- 3
NumberModSets <- 3000
#############################
j = 10
#j=1
if(j < 10) {
  SetNumber <- paste0("0", j)
} else {
  SetNumber <- as.character(j)
}
#
SetName <- paste0(SetNumber, "_", colnames(VariableSetTypes[j]))
SetNameF <- colnames(VariableSetTypes[j])
## Read in variable set data
SubsetVariableNumber <- FinalModelVariables
setwd(OutDirectpsa)
MaxentKeepEvalStats.df1 <- data.frame(read.csv(paste0(Species, ModelName, FullNObs, "_", NumberModSets, "TopSetsFinalTestStats_of", FinalModelVariables, "of", TotVars, "NumVars_.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
tail(MaxentKeepEvalStats.df1)
nrow(MaxentKeepEvalStats.df1)
#
####################################
## First, read in csv file for area of breeding range
setwd(InDirect)
BreedRangeArea.df <- data.frame(read.csv(paste0(Species2, "presareatable.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
## Calculate area in Megahectares of breeding range
BreedingRangeArea <- (BreedRangeArea.df$AREA * 0.0001)/ 1000000 # convert from sq. meters to hectares then megehectares
## Read in Maxent model csv cell counts from directory where grids stored
#ModelTypeNames <- c(paste0("TopAUC", FullNObs), paste0("Random", FullNObs))
ModelTypeNames <- c(paste0("TopAUC",FullNObs))
MaxentKeepEvalStatsL <- list()
##############################################################################
Count <- 0
t1 <- Sys.time()
for(ModelType in ModelTypeNames) {
  #ModelType <- paste0("TopAUC", FullNObs)
  # Keep only DataTypes including string of ModelType for Final Test
  MaxentKeepEvalStats.df2 <- MaxentKeepEvalStats.df1[grep(paste0(ModelType, "FinalTest"), MaxentKeepEvalStats.df1$DataType), ]
  nrow(MaxentKeepEvalStats.df2)
  head(MaxentKeepEvalStats.df2)
  ##
  for(k in 1:NumProjections) {
  #for(k in 1) {
    #k=30
    Loop <- k
    # Keep only data for variable subset k of loop
    MaxentKeepEvalStats.df3 <- MaxentKeepEvalStats.df2[Loop,]
    ##########
    setwd(OutDirectpsa)
    output3 <- paste0("MaxentProj_", ModelType, "_Subsetof", NumberModSets, "_", SetNameF, FinalModelVariables, "_ThresholdCalib", Loop)
    OutDirectSub <- paste0(OutDirectpsa, "/", output3)
    ##
    setwd(OutDirectSub)
    ## Read in csv files of maxent grid cell counts inside breeding area
    MaxentProjBreedRangeArea.df <- data.frame(read.csv(paste0(Species, Loop, "mxt", FinalModelVariables, "apresareatable.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
    ## Calculate Megahectares area overlapping breeding range
    BreedOverlapArea <- (MaxentProjBreedRangeArea.df$AREA * 0.0001)/ 1000000 # convert from sq. meters to hectares then megehectares
    ## Calculate percent of breeding range occupied by maxent projection, Regional Fraction Index (RFI)
    RegionalFractionIndex <- BreedOverlapArea/BreedingRangeArea
    ## Read in csv files of maxent grid cell counts outside breeding area
    MaxentProjNonBreedRangeArea.df <- data.frame(read.csv(paste0(Species, Loop, "mxt", FinalModelVariables, "anpresareatable.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
    ## Calculate Megahectares area overlapping breeding range
    NonBreedArea <- (MaxentProjNonBreedRangeArea.df$AREA * 0.0001)/ 1000000 # convert from sq. meters to hectares then megehectares
    ## Calculate percent of breeding range represented by maxent projection outside of breeding range in Background Evaluation Extent (BEE), Regional Excess Index (REI)
    RegionalExcessIndex <- NonBreedArea/BreedingRangeArea
    ## Calculate Regional Confinement Index (RCI)
    RegionalBoundingIndex <- 1 - (RegionalExcessIndex/RegionalFractionIndex)
    ## Add data to evaluation statistics
    ## Add Breeding Range Area in Megahectares
    MaxentKeepEvalStats.df3$BRA <- BreedingRangeArea
    ## Then add RFI, REI, and RCI
    MaxentKeepEvalStats.df3$RFI <- RegionalFractionIndex
    MaxentKeepEvalStats.df3$REI <- RegionalExcessIndex
    MaxentKeepEvalStats.df3$RBI <- RegionalBoundingIndex
    # Record model number/rank from original Wrapper ranking
    MaxentKeepEvalStats.df3$ModelNumber <- Loop
    ncol(MaxentKeepEvalStats.df3)
    MaxentKeepEvalStats.df3 <- MaxentKeepEvalStats.df3[,c(1:5, 22, 6:21)]
    ## Save output to list
    Count <- Count + 1
    MaxentKeepEvalStatsL[[Count]] <- MaxentKeepEvalStats.df3
  }
}
t2 <- Sys.time()
#############################################################################
difftime(t2,t1, units = "mins")
###
MaxentKeepEvalStatsAll.df <- do.call(rbind, MaxentKeepEvalStatsL)
ncol(MaxentKeepEvalStatsAll.df)
##################
## Use Multi-Objective Optimization from MCDM package to sort by two criteria of AUC and AICc_bg
## Assign rank weightings
RFIrank <- 1/3
REIrank <- 1/3
RCIrank <- 1/3
decision.mat <- as.matrix(MaxentKeepEvalStatsAll.df[,c(20:22)])
head(decision.mat)
# Normalize AICc_bg and AUC from zero to 1
normalizef <- function(x) c((x+(max(x)-min(x))-max(x))/(max(x)-min(x)))
decision.matn <- as.matrix(apply(decision.mat,2,normalizef))
head(decision.matn)
# Assign weights to criteria
weightscrit <- c(RFIrank, REIrank, RCIrank)
# Assign whether cost "min", or benefit "max"
# NOTE: for presence/absence grid use "max" for RFIrank; for breeding range area grid of pseudoabsence data, use "min" for RFIrank
cb <- c("max", "min", "max") 
# Rank criteria
MMOORA.Rank <- MMOORA(decision.matn, weightscrit, cb)
head(MMOORA.Rank)
# Join MMOORA rank to original data
MaxentKeepEvalStatsAll.df$MMOORA_Rank <- MMOORA.Rank[,8]
head(MaxentKeepEvalStatsAll.df)
# Sort data by MMOORA rank
MaxentKeepEvalStatsAll.dfs <- arrange(MaxentKeepEvalStatsAll.df, MMOORA_Rank)
## Save data
setwd(OutDirectpsa)
write.table(MaxentKeepEvalStatsAll.dfs, file=paste0(Species, "Maxent", DataSet, "_TrainVsTestandAreaStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", EvalType, "_TopVsRandom", NumProjections, "Sets.csv"), sep=",", col.names=NA)


#####################################################################################################
### For FinalNumberModels selected models by regional indices, determine the number of times each 
### variable is used 
#####################################################################################################
## Retrieve top selected variable from selected subset of 90 variables
# Specify variable number and row number of selected set
FinalNumberModels <- 4 # Select number of top models to keep by Regional Indices
FinalModelVariables <- 15  # Chosen variable set size based trends in AUC and AICc
SubsetVariableNumber <- FinalModelVariables
NumProjections <- 30
NumberModSets <- 3000
#
j = 10
#j=1
if(j < 10) {
  SetNumber <- paste0("0", j)
} else {
  SetNumber <- as.character(j)
}
#
SetName <- paste0(SetNumber, "_", colnames(VariableSetTypes[j]))
SetNameF <- colnames(VariableSetTypes[j])
# Read in data with rankings
setwd(OutDirectpsa)
MaxentKeepEvalAreaStatsAll.dfs <- data.frame(read.csv(paste0(Species, "Maxent", DataSet, "_TrainVsTestandAreaStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", EvalType, "_TopVsRandom", NumProjections, "Sets.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
# Keep only rankings for model type
MaxentKeepEvalAreaStatsAll.df2 <- MaxentKeepEvalAreaStatsAll.dfs[which(MaxentKeepEvalAreaStatsAll.dfs$DataType==paste0("TopAUC", FullNObs, "FinalTest")),]
FinalModelRowNumberList.df <- as.data.frame(MaxentKeepEvalAreaStatsAll.df2$ModelNumber)
FinalModelRowNumberList <- FinalModelRowNumberList.df[1:FinalNumberModels,]
#
# Keep only DataTypes including string for Top FinalTest
TopSelSubsets <- as.data.frame(MaxentKeepEvalAreaStatsAll.df2[1:FinalNumberModels,1], stringsAsFactors=FALSE)
colnames(TopSelSubsets) <- "VarNames"
#
###############################
### Deconstruct VarNames for top subsets and talley counts of each variable among top NumProjections models
VarNamesL <- list()
for(m in 1:FinalNumberModels) {
  #m=2
  VarNamesL[[m]] <- unlist(strsplit(as.character(TopSelSubsets[m,1]), "-"))
}
VarNamesCount.df <- as.data.frame(table(unlist(VarNamesL)), stringsAsFactors=FALSE)
## Save results
write.table(VarNamesCount.df, file=paste0(Species, ModelName, FullNObs, "_", NumberModSets, "Sets_VariableCountsofTop", FinalNumberModels, "Models_", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff.csv"), sep=",", col.names=NA)
#

#######################################################################################
### Obtain mean permutation with n for all variables in top models
### and perform Multi-Objective Optimization Ranking for mean permutation and n (count)
#######################################################################################
##
#ModelTypeNames <- c(paste0("TopAUC",FullNObs))
#ModelTypeNames <- c(paste0("Random",FullNObs))
#ModelType <- ModelTypeNames[1]
#NumProjections <- 6
#SubsetVariableNumber <- 10
PermImp.List <- list()
for(p in 1:NumProjections) {
  #p=1
  setwd(OutDirectpsa)
  output3 <- paste0("MaxentProj_", ModelType, "_Subsetof", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_ThresholdCalib", p)
  #output3 <- paste0("MaxentProj_", ModelType, "_Subsetof", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_ThresholdCalib", Loop, "Minus", 14)
  OutDirectSub <- paste0(OutDirectpsa, "/", output3)
  setwd(OutDirectSub)
  maxentResults.df <- data.frame(read.csv("maxentResults.csv", header = TRUE, sep=','), stringsAsFactors=FALSE)
  PermImp.df <- maxentResults.df[,(7+SubsetVariableNumber+1):(6+(2*SubsetVariableNumber)+1)]
  colnames(PermImp.df) <- gsub(".permutation.importance", "", colnames(PermImp.df))
  PermImp.df2 <- data.frame(t(as.matrix(PermImp.df)), stringsAsFactors=FALSE)
  PermImp.df3 <- data.frame(cbind(rownames(PermImp.df2), PermImp.df2[,1]), stringsAsFactors=FALSE)
  PermImp.df3[,2] <- as.numeric(PermImp.df3[,2])
  #str(PermImp.df3)
  colnames(PermImp.df3) <- c("VarNames", "Permutation_Importance")
  PermImp.List[[p]] <- PermImp.df3
}
##
PermutationImp.df <- do.call(rbind, PermImp.List)
# Drop all permutation importance values of zero since variable not used
PermutationImp.df <- PermutationImp.df[apply(PermutationImp.df!=0,1,all),]
## Sort by VarName
PermutationImp.dfs <- PermutationImp.df[order(PermutationImp.df$VarNames),]
#str(PermutationImp.dfs)
## Calculate mean and standard deviation per VarName
PermutationImp.means <- aggregate(Permutation_Importance ~ VarNames, PermutationImp.dfs, mean)
PermutationImp.sd <- aggregate(Permutation_Importance ~ VarNames, PermutationImp.dfs, sd)
# Replace NAs for sd with zero
PermutationImp.sd[is.na(PermutationImp.sd)] <- 0
# Calculate count for each variable
PermutationImp.count <- aggregate(Permutation_Importance ~ VarNames, PermutationImp.dfs, FUN = function(x){NROW(x)})
PermutationImp <- cbind(PermutationImp.means, PermutationImp.sd[,2], PermutationImp.count[,2])
colnames(PermutationImp) <- c("VarNames", "MeanPermImp", "SDPermImp", "Count")
## Sort by MeanPermImp
PermutationImp <- PermutationImp[order(-PermutationImp$MeanPermImp),]
## Create column with joint ranking by MeanPermImp and Count using Multi-Objective Optimization ranking
## Use Multi-Objective Optimization from MCDM package to sort by two criteria of MeanPermImp and Count
decision.mat <- as.matrix(PermutationImp[,c(2,4)])
head(decision.mat)
# Normalize variable mean permutation importance and variable frequency in models from zero to 1
normalizef <- function(x) c((x+(max(x)-min(x))-max(x))/(max(x)-min(x)))
decision.matn <- as.matrix(apply(decision.mat,2,normalizef))
head(decision.matn)
# Assign weights to criteria
PermImprank <- 0.6
Countrank <- 0.4
weightscrit <- c(PermImprank, Countrank)
# Assign whether cost "min", or benefit "max"
cb <- c("max", "max")
# Rank criteria
MMOORA.Rank <- MMOORA(decision.matn, weightscrit, cb)
head(MMOORA.Rank)
# Join MMOORA rank to original data
PermutationImp$MMOORA_Rank <- MMOORA.Rank[,8]
# Sort data by MMOORA rank
PermutationImp2 <- PermutationImp[order(PermutationImp$MMOORA_Rank),]
head(PermutationImp2)
# Obtain lower case versions of variable names
PermutationImp2[,1] <- tolower(PermutationImp2[,1])
# Save results
setwd(OutDirectpsa)
write.table(PermutationImp2, file=paste0(Species, ModelName, "_Top", NumProjections, "Models_", SubsetVariableNumber, "Vars_JointVariableRankingsbyMeanPermImp", PermImprank, "andCount", Countrank, ".csv"), sep=",", col.names=NA)
#

#######################################################################################

#################
## NOTE: DO NOT RESORT BY THIS AUCfinaltest- THIS WILL INFLATE AUC- AUCWrappertest already used to sort
#################
#####################################################################################################
### Run maxent 3-fold cross-valication evaluation mode with jacknife option to assess variable importance with all variables
### for top model
#####################################################################################################

## Takes about 17 minutes with 15 variables and three models
## Takes 8 minutes with 12 variables and three models
#############################################
NumberTopModels <- seq(1:1) # Specify number of top models on which to run jacknife analysis
TopModelRowNumberList <- FinalModelRowNumberList.df[1:NumberTopModels,]
##############
t1 <- Sys.time()
for(NumberModel in 1:NumberTopModels) {
  #NumberModel <- 1
  TopModel <- FinalModelRowNumberList.df[NumberModel,]
  VariableSubset <- data.frame(as.character(TopSelSubsets[NumberModel,1]), stringsAsFactors=FALSE)
  colnames(VariableSubset) <- "VarNames"
  VarNames <- strsplit(as.character(VariableSubset[,1]),"-")
  SubsetVariableNumber <- length(unlist(VarNames))
  #########
  #
  VariableNamesIn <- VariableNames
  # Specify SubsetVariableNumber
  # Specify model training data
  MaxentPresTrainData <- PresenceDat.df
  head(MaxentPresTrainData)
  MaxentAbsTrainData <- BackgroundDat.df
  head(MaxentAbsTrainData)
  TrainSWD <- rbind(MaxentPresTrainData[,4:ncol(MaxentPresTrainData)], MaxentAbsTrainData[,4:ncol(MaxentAbsTrainData)])
  head(TrainSWD)
  TrainPresID <- data.frame(rep(1,nrow(MaxentPresTrainData)))
  colnames(TrainPresID) <- "ID"
  TrainAbsID <- data.frame(rep(0,nrow(MaxentAbsTrainData)))
  colnames(TrainAbsID) <- "ID"
  TrainPresAbsID <- rbind(TrainPresID, TrainAbsID)
  head(TrainPresAbsID)
  tail(TrainPresAbsID)
  #
  # Subset training and testing data by selected variables
  TrainSWD <- TrainSWD[,unlist(VarNames)]
  head(TrainSWD)
  ######
  # Create subdirectory for Maxent output
  output3 <- paste0("MaxentEvaluationMode_FinalSubset_", SetNameF, SubsetVariableNumber, "_Set", TopModel)
  dir.create(paste0(OutDirectpsa, "/", output3))
  OutDirectSub <- paste0(OutDirectpsa, "/", output3)
  setwd(OutDirectSub)
  #
  #system.file("java", package="dismo")
  ## Check if have categorical variable to add to maxent arguments
  if(exists("CatVarsPrefix")==TRUE) {
    MaxentCatArg <- paste0("togglelayertype=", CatVarsPrefix)
    MaxentEvalArgs <- c(MaxentEvalArgs1, MaxentCatArg)
    } else {
    MaxentEvalArgs <- MaxentEvalArgs1
  }
  # NOTE: the below jackknife run takes about 7.5 minutes with 20 variables
  # Takes 85 seconds for 8 variables
  ## Takes 42 minutes for 86 variables
  system.time(MaxentOut <- maxent(TrainSWD, TrainPresAbsID, args=MaxentEvalArgs, path=OutDirectSub))
}
t2 <- Sys.time()
###################################################################################
difftime(t2,t1, units = "mins")


######################################################################################################
#### Run Future Climate Models using top 3 of 20 selected models and top 3 of 20 random models
#### ranked by regional indices
#####################################################################################
## Takes 38 minutes to loop through 3 models with four climate scenarios for flycatcher
## Loop takes 2.7 hours for 3 selected and 3 random models with four climate scenarios each for flycatcher
## Loop takes about 55 minutes for 3 selected models with 12 variables each and four climate scenarios
## Loop takes 2.8 hours for 7 models with 12 variables and 4 scenarios
# Specify variable number and row number of selected set
FinalModelVariables <- 15  # Chosen variable set size based trends in AUC and AICc
SubsetVariableNumber <- FinalModelVariables
#FinalModelRowNumberList <- seq(1,10,1) # Rank of model by AUCpsa_Wrappertest; chosen based on visual assessment and ranks of RFI, REI, and RBI among top models
#TopFinalModelRowNumberList <- c(15,14,5) # Rank of model by RFI, REI, and RBI among top models
#RandomFinalModelRowNumberList <- c(19,5,10) # Rank of model by RFI, REI, and RBI among random models
#FinalModelRowNumberList <- c(1,2,5,6,7,8,10) # Rank of model by AUCpsa_Wrappertest; chosen based on visual assessment and ranks of RFI, REI, and RBI among top models
NumberModSets <- 3000
j = 10
#j=1
if(j < 10) {
  SetNumber <- paste0("0", j)
} else {
  SetNumber <- as.character(j)
}
#
SetName <- paste0(SetNumber, "_", colnames(VariableSetTypes[j]))
SetNameF <- colnames(VariableSetTypes[j])
###################
## Read in stack of 33 clipped non-climate grids
setwd(MaskGrids)
NonClimPredictors <- stack(NonClimateGridNamesL, proj4string=CRS.WGS84)
names(NonClimPredictors) <- toupper(names(NonClimPredictors))
names(NonClimPredictors) <- gsub("_NS", "", names(NonClimPredictors))
#
## Go to root directory for clipped future climate grids
setwd(FutClimDir)
NumProjections <- 9
FinalNumberModels <- 9 # Select number of top models to keep by Regional Indices
#ModelTypeNames <- c(paste0("TopAUC",FullNObs), paste0("Random", FullNObs))
ModelTypeNames <- c(paste0("TopAUC",FullNObs))
#ModelTypeNames <- c(paste0("Random",FullNObs))
#FutClimSubDirList <- c("2070he26")
FutClimSubDirList <- c("2050he26", "2050he85", "2070he26", "2070he85")
####################
t1 <- Sys.time()
for(FutClimSubDir in FutClimSubDirList) {
  #FutClimSubDir <- FutClimSubDirList[2]
  ## Read in future climate grids
  FutClimGridDir <- paste0(FutClimDir, "/", FutClimSubDir)
  setwd(FutClimGridDir)
  ClimPredictors <- stack(FutureGridNamesL, proj4string=CRS.WGS84)
  names(ClimPredictors) <- toupper(names(ClimPredictors))
  names(ClimPredictors) <- gsub("_NS", "", names(ClimPredictors))
  ## Merge NonClimate and Future Climate grids in one stack
  PredictorsAll <- stack(ClimPredictors, NonClimPredictors)
  names(PredictorsAll)
  #plot(PredictorsAll[[1]])
  for(ModelType in ModelTypeNames) {
    #ModelType <- paste0("TopAUC", FullNObs)
    # Find variable subsets for top ranked three models by regional indices of projections for ModelType
    # Read in data with rankings
    setwd(OutDirectpsa)
    MaxentKeepEvalAreaStatsAll.df <- data.frame(read.csv(paste0(Species, "Maxent", DataSet, "_TrainVsTestandAreaStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", EvalType, "_TopVsRandom", NumProjections, "Sets.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
    # Keep only rankings for model type
    MaxentKeepEvalAreaStatsAll.df2 <- MaxentKeepEvalAreaStatsAll.df[which(MaxentKeepEvalAreaStatsAll.df$DataType==paste0(ModelType, "FinalTest")),]
    FinalModelRowNumberList.df <- as.data.frame(MaxentKeepEvalAreaStatsAll.df2$ModelNumber)
    #FinalModelRowNumberList <- FinalModelRowNumberList.df[4,]
    FinalModelRowNumberList <- FinalModelRowNumberList.df[1:FinalNumberModels,]
    #FinalModelRowNumberList <- seq(6,9,1)
    #FinalModelRowNumberList <- c(3,5,6,7,8,9)
    for(FinalModelRowNumber in FinalModelRowNumberList) {
      #FinalModelRowNumber <- FinalModelRowNumberList[1]
      #FinalModelRowNumber <- 20
      ## Read in saved calibration threshold and variables of final chosen maxent model
      setwd(OutDirectpsa)
      output3 <- paste0("MaxentProj_", ModelType, "_Subsetof", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_ThresholdCalib", FinalModelRowNumber)
      #SubsetVariableNumber <- 14
      #output3 <- paste0("MaxentProj_", ModelType, "_Subsetof", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_ThresholdCalib", Loop, "Minus", 14)
      OutDirectSub <- paste0(OutDirectpsa, "/", output3)
      setwd(OutDirectSub)
      ## Read in calibration threshold
      MaxentKeepEvalStats.df <- data.frame(read.csv(paste0(Species, "Maxent", DataSet, "_TrainVsTestStatsThreshCal_of", NumberModSets, "_", SetName, SubsetVariableNumber, "_", FinalModelRowNumber, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
      #
      SubsetTestMeans <- MaxentKeepEvalStats.df[which(MaxentKeepEvalStats.df$DataType=="Test"),]
      ## Find mean threshold at maximum TSS (from above output)
      Threshold <- SubsetTestMeans$ThreshMxTSS
      ThresholdK <- Threshold*1000
      ## Retrieve variable names from chosen subset
      SubsetVariableNumber <- SubsetTestMeans[,2]
      VariableSubset <- as.character(SubsetTestMeans[,1])
      VarNames <- strsplit(as.character(VariableSubset),"-")
      ## Read in saved Maxent model
      MaxentMod <- readRDS("MaxentModel.rds")
      #### Subset Predictors by Variables
      PredictorsFut <- subset(PredictorsAll, unlist(VarNames))
      ################
      ## Specify output directories
      output3 <- paste0("MaxentProj_", ModelType, "_Subsetof", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_ThresholdCalib", FinalModelRowNumber, "_", FutClimSubDir)
      #output3 <- paste0("MaxentProj_", ModelType, "_Subsetof", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_ThresholdCalib", FinalModelRowNumber, "_Minus14_", FutClimSubDir)
      dir.create(paste0(OutDirectpsa, "/", output3))
      OutDirectSub <- paste0(OutDirectpsa, "/", output3)
      # Generate grid model projection for future climate grids
      system.time(maxent.score1 <- predict(MaxentMod, PredictorsFut))
      #plot(maxent.score1)
      # Multiply raw grid by 1000 and convert to integer
      maxent.score <- calc(maxent.score1, function(x) as.integer(x * 1000) )
      #plot(maxent.score)
      # Save for Arc as a geoTIFF grid
      OutGridID <- paste0(ModelType, "_", FinalModelRowNumber, "ThreshCal_", SetNameF, "_", FutClimSubDir)
      setwd(OutDirectSub)
      writeRaster(maxent.score, paste0(Species, "Maxent", OutGridID, "_", SubsetVariableNumber, "Vars_Beta", BetaMult), format = "GTiff", overwrite=TRUE)
      #
      ################################################################################
      ### Calibrate model using using above calculated ThresholdK at maximum TSS
      ### for current climate model
      maxent.scorecal <- calc(maxent.score, function(x) ifelse(x < ThresholdK, 0, 1) )
      #plot(maxent.scorecal)
      writeRaster(maxent.scorecal, paste0(Species, "Maxent", OutGridID, "_", SubsetVariableNumber, "Vars_Beta", BetaMult, "Cal"), format = "GTiff", overwrite=TRUE)
      }
    }
}
t2 <- Sys.time()
#############################################################################
difftime(t2,t1, units = "mins")

#################################################################################
### Analyze regional indices for top four models for current and future climates
#################################################################################
##
####################################################################################
### Use Arc Python script ENMBreedingAreaZonalStatisticsLooopSWFLFuture3ModsTopVsRandom.py to calculate percent
### overlap of presence projection inside breeding range with total breeding range area
### and calculate percent overlap of presence projection outside breeding range with total breeding range area
### within NAD_1983_Albers projection
####################################################################################
## Read in above produced csv files from python and calculate percent overlaps of maxent models
## with breeding range for portions of maxent model both inside and outside breeding range
##
FinalModelVariables <- 15
#TopFinalModelRowNumberList <- c(3,6,18) # Rank of model by RFI, REI, and RBI among top models
#RandomFinalModelRowNumberList <- c(1,8,13) # Rank of model by RFI, REI, and RBI among random models
NObsJ <- 3
FinalNumberModels <- 4 # Select number of top models to keep by Regional Indices
NumberModSets <- 3000
NumProjections <- 30
ModelTypeList <- c(paste0("TopAUC", FullNObs), paste0("RandAUC", FullNObs))
#############################
j = 10
#j=1
if(j < 10) {
  SetNumber <- paste0("0", j)
} else {
  SetNumber <- as.character(j)
}
#
SetName <- paste0(SetNumber, "_", colnames(VariableSetTypes[j]))
SetNameF <- colnames(VariableSetTypes[j])
## Read in variable set data
SubsetVariableNumber <- FinalModelVariables
##########
# Read in variable information for models
setwd(OutDirectpsa)
MaxentKeepEvalStats.df1 <- data.frame(read.csv(paste0(Species, ModelName, FullNObs, "_", NumberModSets,"SetNumber_TrainVsTestStats_of", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
# Read in data with VarNames and Model Number Rankings
# Find variable subsets for top ranked three models by regional indices of projections for ModelType
MaxentKeepEvalAreaStatsAll.df <- data.frame(read.csv(paste0(Species, "Maxent", DataSet, "_TrainVsTestandAreaStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", EvalType, "_TopVsRandom", NumProjections, "Sets.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
#
####################################
## First, read in csv file for area of breeding range
setwd(InDirect)
BreedRangeArea.df <- data.frame(read.csv(paste0(Species2, "presareatable.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
## Calculate area in Megahectares of breeding range
BreedingRangeArea <- (BreedRangeArea.df$AREA * 0.0001)/ 1000000 # convert from sq. meters to hectares then megehectares
## Read in Maxent model csv cell counts from directory where grids stored
MaxentKeepEvalStatsL <- list()
##############################################################################
Count <- 0
#ModelTypeNames <- c(paste0("TopAUC",FullNObs), paste0("Random", FullNObs))
ModelTypeNames <- c(paste0("TopAUC",FullNObs))
t1 <- Sys.time()
###########
for(ModelType in ModelTypeNames) {
  #ModelType <- paste0("TopAUC", FullNObs)
  # Find variable subsets for top ranked three models by regional indices of projections for ModelType
  # Read in data with rankings
  setwd(OutDirectpsa)
  MaxentKeepEvalAreaStatsAll.dfs <- data.frame(read.csv(paste0(Species, "Maxent", DataSet, "_TrainVsTestandAreaStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", EvalType, "_TopVsRandom", NumProjections, "Sets.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
  # Keep only rankings for model type
  MaxentKeepEvalAreaStatsAll.df2 <- MaxentKeepEvalAreaStatsAll.dfs[which(MaxentKeepEvalAreaStatsAll.dfs$DataType==paste0("TopAUC", FullNObs, "FinalTest")),]
  FinalModelRowNumberList.df <- as.data.frame(MaxentKeepEvalAreaStatsAll.df2$ModelNumber)
  FinalModelRowNumberList <- FinalModelRowNumberList.df[1:FinalNumberModels,]
  # Find subset evaluation stats for ModelType
  MaxentKeepEvalStats.df2 <- MaxentKeepEvalStats.df1[which(MaxentKeepEvalStats.df1$SubsetVariableNumber==FinalModelVariables & MaxentKeepEvalStats.df1$DataType==paste0(ModelType, "FinalTest")),]
  ####
  for(FinalModelRowNumber in FinalModelRowNumberList) {
    #FinalModelRowNumber=17
    for(FutClimSubDir in FutClimSubDirList) {
      #FutClimSubDir <- FutClimSubDirList[[1]]
      SelectedSubsetData <- MaxentKeepEvalAreaStatsAll.df2[which(MaxentKeepEvalAreaStatsAll.df2$ModelNumber==FinalModelRowNumber),]
      ## Subset original evaluation data using VarNames
      selectedRows <- (MaxentKeepEvalStats.df2$VarNames %in% SelectedSubsetData$VarNames)
      MaxentTestSubsetEvalStats.df <- MaxentKeepEvalStats.df2[selectedRows,]
      #
      TopModelSet.df <- MaxentTestSubsetEvalStats.df[,1:8]
      ## Specify input directories
      output3 <- paste0("MaxentProj_", ModelType, "_Subsetof", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_ThresholdCalib", FinalModelRowNumber, "_", FutClimSubDir)
      OutDirectSub <- paste0(OutDirectpsa, "/", output3)
      ##
      setwd(OutDirectSub)
      ## Read in csv files of maxent grid cell counts inside breeding area
      MaxentProjBreedRangeArea.df <- data.frame(read.csv(paste0(SpecShort, FutClimSubDir, "presareatable.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
      ## Calculate Megahectares area overlapping breeding range
      BreedOverlapArea <- (MaxentProjBreedRangeArea.df$AREA * 0.0001)/ 1000000 # convert from sq. meters to hectares then megehectares
      ## Calculate percent of breeding range occupied by maxent projection, Regional Fraction Index (RFI)
      RegionalFractionIndex <- BreedOverlapArea/BreedingRangeArea
      ## Read in csv files of maxent grid cell counts outside breeding area
      MaxentProjNonBreedRangeArea.df <- data.frame(read.csv(paste0(SpecShort, FutClimSubDir, "npresareatable.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
      ## Calculate Megahectares area overlapping breeding range
      NonBreedArea <- (MaxentProjNonBreedRangeArea.df$AREA * 0.0001)/ 1000000 # convert from sq. meters to hectares then megehectares
      ## Calculate percent of breeding range represented by maxent projection outside of breeding range in Background Evaluation Extent (BEE), Regional Excess Index (REI)
      RegionalExcessIndex <- NonBreedArea/BreedingRangeArea
      ## Calculate Regional Bounding Index (RBI)
      RegionalBoundingIndex <- 1 - (RegionalExcessIndex/RegionalFractionIndex)
      ## Add data to model information
      ## Add Breeding Range Area in Megahectares
      MaxentKeepEvalStats.df3 <- TopModelSet.df
      MaxentKeepEvalStats.df3$DataType <- ModelType
      MaxentKeepEvalStats.df3$BRA <- BreedingRangeArea
      ## Then add RFI, REI, and RBI
      MaxentKeepEvalStats.df3$RFI <- RegionalFractionIndex
      MaxentKeepEvalStats.df3$REI <- RegionalExcessIndex
      MaxentKeepEvalStats.df3$RBI <- RegionalBoundingIndex
      # Record model number/rank from original Wrapper ranking
      MaxentKeepEvalStats.df3$FutClimate <- FutClimSubDir
      MaxentKeepEvalStats.df3$ModelNumber <- FinalModelRowNumber
      ncol(MaxentKeepEvalStats.df3)
      MaxentKeepEvalStats.df3 <- MaxentKeepEvalStats.df3[,c(1, 13, 2:8, 14, 9:12)]
      ## Save output to list
      Count <- Count + 1
      MaxentKeepEvalStatsL[[Count]] <- MaxentKeepEvalStats.df3
    }
  }
}
t2 <- Sys.time()
#############################################################################
difftime(t2,t1, units = "mins")
## Save data
MaxentKeepEvalStatsL[[1]]
MaxentKeepEvalStatsFuture.df <- do.call(rbind, MaxentKeepEvalStatsL)
ncol(MaxentKeepEvalStatsFuture.df)
nrow(MaxentKeepEvalStatsFuture.df)
head(MaxentKeepEvalStatsFuture.df)
#
setwd(OutDirectpsa)
write.table(MaxentKeepEvalStatsFuture.df, file=paste0(Species, "MaxentFutureClimate", "_AreaStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", "_TopVsRandom", NumProjections, "Sets.csv"), sep=",", col.names=NA)
###################################################################################

######################################################################################
## Assemble data for bar charts and ANOVA for current and future climates
##########################################################################################
setwd(OutDirectpsa)
##
NumProjections <- 30
## Read in regional index area data and combine for analysis
#FinalModelRowNumberList <- c(9,4,3)
# Read in regional indices for current models
MaxentAreaStatsCurrent.df <- data.frame(read.csv(paste0(Species, "Maxent", DataSet, "_TrainVsTestandAreaStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", EvalType, "_TopVsRandom", NumProjections, "Sets.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
# Reformat to match future climate data
CurrentRegionalIndices.df <- MaxentAreaStatsCurrent.df[,c(1:5,7:9,6,19:22)]
CurrentRegionalIndices.df$Climate <- "Current"
ncol(CurrentRegionalIndices.df)
CurrentRegionalIndices.df1 <- CurrentRegionalIndices.df[,c(1, 14, 2:13)]
# Keep only data for 3 top selected and random models
#ModelTypeNames <- c(paste0("TopAUC",FullNObs), paste0("Random", FullNObs))
ModelTypeNames <- c(paste0("TopAUC",FullNObs))
CurrentRegionalIndicesL <- list()
Count <- 0
##################
for(ModelType in ModelTypeNames) {
    # Keep only rankings for model type
    CurrentRegionalIndices.df2 <- CurrentRegionalIndices.df1[which(CurrentRegionalIndices.df1$DataType==paste0(ModelType, "FinalTest")),]
    Count <- Count + 1
    CurrentRegionalIndicesL[[Count]] <- CurrentRegionalIndices.df2[1:FinalNumberModels,]
}
CurrentRegionalIndices.df <- do.call(rbind, CurrentRegionalIndicesL)
ncol(CurrentRegionalIndices.df)
# Read in regional indices for future models
MaxentAreaStatsFuture.df <- data.frame(read.csv(paste0(Species, "MaxentFutureClimate", "_AreaStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", "_TopVsRandom", NumProjections, "Sets.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
colnames(MaxentAreaStatsFuture.df)[2] <- "Climate"
head(MaxentAreaStatsFuture.df)
ncol(MaxentAreaStatsFuture.df)
# Combine current and future regional indices together
MaxentAreaStatsClimate.df <- rbind(CurrentRegionalIndices.df, MaxentAreaStatsFuture.df)
MaxentAreaStatsClimate.dfs <- arrange(MaxentAreaStatsClimate.df, Climate, DataType, ModelNumber, Climate)
# Create ClimateDataType combining Climate and Set
MaxentAreaStatsClimate.dfs$ClimDataType <- paste0(MaxentAreaStatsClimate.dfs$Climate, "_", MaxentAreaStatsClimate.dfs$Set)
# Rearrange columns
MaxentAreaStatsClimate.dfs <- MaxentAreaStatsClimate.dfs[,c(1:3,15,5:14)]
head(MaxentAreaStatsClimate.dfs)
ncol(MaxentAreaStatsClimate.dfs)
## Save combined data
write.table(MaxentAreaStatsClimate.dfs, file=paste0(Species, "MaxentCurrentFutureClimate", "_AreaStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", "_TopVsRandom", NumProjections, "Sets.csv"), sep=",", col.names=NA)
###
### Calculate Mean and Standard Deviation Values
EvaluationStats <- MaxentAreaStatsClimate.dfs
OutName <- paste0(Species, "MaxentProjection_Top3CurrentFutureModelAreaStats_", SubsetVariableNumber, "Vars_", "_Summary.csv")
SortGroups <- c("ClimDataType", "Set", "SetName", "SubsetVariableNumber")
## All statistics, and only statistics, should be at and after the column specified below
StatVarFirstColumn <- 10
OutDirectIn <- OutDirectpsa
EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
head(EvaluationStatsSummary.df)
##
EvalStatsOut <- EvaluationStatsSummary.df
## Sort data by ClimDataType for ease of plotting in excel
EvalStatsOuts <- arrange(EvalStatsOut, ClimDataType, Statistic, SubsetVariableNumber)
EvalStatsOuts <- arrange(EvalStatsOut, Statistic)
head(EvalStatsOuts)
write.table(EvalStatsOuts, file=paste0(Species, "MaxentProjection_Top3CurrentFutureModelAreaStats_", SubsetVariableNumber, "Vars_", "_SummaryTable.csv"), sep=",", col.names=NA)

##################################################################################
### Loop through various Welch t test comparisons for comparing statistics for
### Current Climate with each of the four future scenarios for both selected and random
### models. Use Holm correction to penalize p-value for four comparisons of each type
##
StatTypeNames <- c("RFI", "REI", "RBI")
#ModelTypeNames <- c(paste0("TopAUC", FullNObs), paste0("Random", FullNObs))
ModelTypeNames <- c(paste0("TopAUC",FullNObs))
ClimateTypeNames <- c("2050he26", "2050he85", "2070he26", "2070he85")
NumberComparisons <- length(ClimateTypeNames)
##
for(StatType in StatTypeNames) {
  #StatType <- "RFI"
  for(ModelType in ModelTypeNames) {
    #ModelType <- paste0("TopAUC", FullNObs)
    MaxentEvalStatsL <- list()
    Count <- 0
    ## Keep only ModelType data
    MaxentTestStats1 <- MaxentAreaStatsClimate.dfs[which(MaxentAreaStatsClimate.dfs$Set==ModelType),]
    ## Use Welch correction one way ANOVA for unequal variances for primary test among all possible climate comparisons
    WelchANOVA <- oneway.test(as.formula(paste0(StatType, " ~ ClimDataType")), data=MaxentTestStats1, na.action=na.omit, var.equal=FALSE)
    # Save above Welch's ANOVA output
    out1 <- paste0("Climate Projection Analysis for ", ModelType, " Models: Welch ANOVA for ", StatType)
    cat(out1, file=paste0(Species, ModelType, "_Climate_WelchANOVA_forCurrentandFutureModelsof_", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
    out2 <- capture.output(print(WelchANOVA))
    cat(out2, file=paste0(Species, ModelType, "_Climate_WelchANOVA_forCurrentandFutureModelsof_", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
    # Conduct specific pre-planned climate comparisons between pairs of current climate with each future climate
    for(ClimateType in ClimateTypeNames) {
      #ClimateType <- "2050he26"
      ## Keep only Current Climate data and ClimateType
      MaxentTestStats <- MaxentTestStats1[which(MaxentTestStats1$Climate=="Current" | MaxentTestStats1$Climate==ClimateType),]
      N <- nrow(MaxentTestStats)/2
      ## Use Welch correction t test
      RankWelchttest <- t.test(as.formula(paste0(StatType, " ~ Climate")), data=MaxentTestStats)
      #str(RankWelchttest)
      ## Apply Holm Correction to the p-value based on the total number of comparions for group (NumberComparisons)
      PValueCor <- p.adjust(RankWelchttest$p.value, method= "holm", n = NumberComparisons)
      # Save above Welch t test output
      out <- paste0("   ")
      cat(out, file=paste0(Species, "WelchttestStatsforCurrentVsFutureClimates", "_",ModelType, "_", N, "Setsof", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
      out1 <- paste0("Current vs ", ClimateType, ": Welch t test of ", StatType, " for ", N, " ", ModelType, " Models of ", SubsetVariableNumber, " Variables")
      cat(out1, file=paste0(Species, "WelchttestStatsforCurrentVsFutureClimates", "_",ModelType, "_", N, "Setsof", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
      out2 <- capture.output(print(RankWelchttest))
      cat(out2, file=paste0(Species, "WelchttestStatsforCurrentVsFutureClimates", "_",ModelType, "_", N, "Setsof", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
      out3 <- paste0("p-value with Holm Correction for ", NumberComparisons, " Welch t tests: ", PValueCor)
      cat(out3, file=paste0(Species, "WelchttestStatsforCurrentVsFutureClimates", "_",ModelType, "_", N, "Setsof", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
    }
  }
}



##################################################################################

#################################################################################
### Analyze regional indices for top four models for current and future climates
### for 46 Ecoregions
#################################################################################
##
#####################################################################################
### First, run through the series of python scripts below to reproject maxent model
### raster to albers projection for both the BEE and Projection areas
### ENMAlbersProjectionBEELoopFire.py
### ENMAlbersProjectionPROJLoopFireCurrent.py
### ENMAlbersProjectionPROJLoopFireFuture.py
####################################################################################
##
####################################################################################
### Use Arc Python script ENMBEEAreaZonalStatisticsLoopFireEcoregion.py to calculate percent
### overlap of historic presence inside ecoregion within NAD_1983_Albers projection
####################################################################################
## Read in above produced csv files from python and calculate percent overlaps of maxent models
## with breeding range for portions of maxent model both inside and outside breeding range
##
FinalModelVariables <- 15
NumProjections <- 30
NObsJ <- 3
NumberModSets <- 3000
#############################
j = 10
#j=1
if(j < 10) {
  SetNumber <- paste0("0", j)
} else {
  SetNumber <- as.character(j)
}
#
SetName <- paste0(SetNumber, "_", colnames(VariableSetTypes[j]))
SetNameF <- colnames(VariableSetTypes[j])

####################################
MaxentKeepEvalStatsL <- list()
##############################################################################
t1 <- Sys.time()
setwd(InDirect)
EcoregionHistoricData.df <- as.data.frame(t(as.matrix(c(1,0,0,0,0))))
colnames(EcoregionHistoricData.df) <- c("EcoRegNum", "BEEEcoRegArea", "ProjEcoRegArea", "BEEFirePresArea", "BEEPercFireEcoreg")
for(EcoRegNum in 1:53) {
  #EcoRegNum <- 47
  if(EcoRegNum < 47) {
    ## Read in csv files of historic grid cell counts inside ecoregion area
    HistoricPresenceRangeArea.df <- data.frame(read.csv(paste0(Species, "presareatable", EcoRegNum, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
    ## Calculate Megahectares presence area within Ecoregion
    if(nrow(HistoricPresenceRangeArea.df)==0) {
      BEEFirePresArea <- 0
    } else {
      BEEFirePresArea <- (HistoricPresenceRangeArea.df$AREA * 0.0001)/ 1000000 # convert from sq. meters to hectares then megehectares
    }
    ## Read in csv files of ecoregion area within Background Evaluation Extent
    BEEEcoregionArea.df <- data.frame(read.csv(paste0("BEEareatable_ecoregion", EcoRegNum, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
    ## Calculate Megahectares nonpresence area within Ecoregion
    BEEEcoregionArea <- (BEEEcoregionArea.df$AREA * 0.0001)/ 1000000 # convert from sq. meters to hectares then megehectares
    ## Calculate Percent Ecoregion with Historic fire
    BEEPercFireEcoreg <- BEEFirePresArea/BEEEcoregionArea
  } else {
    BEEFirePresArea <- NA
    BEEEcoregionArea <- NA
    BEEPercFireEcoreg <- NA
  }
  ###
  ## Read in csv files of ecoregion area within Background Evaluation Extent
  ProjEcoregionArea.df <- data.frame(read.csv(paste0("Projareatable_ecoregion", EcoRegNum, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
  ## Calculate Megahectares nonpresence area within Ecoregion
  ProjEcoregionArea <- (ProjEcoregionArea.df$AREA * 0.0001)/ 1000000 # convert from sq. meters to hectares then megehectares
  ### Save data to dataframe
  if(EcoRegNum>1) {EcoregionHistoricData.df[EcoRegNum,] <- EcoregionHistoricData.df[1,]}
  EcoregionHistoricData.df[EcoRegNum,1] <- EcoRegNum
  EcoregionHistoricData.df[EcoRegNum,2] <- BEEEcoregionArea
  EcoregionHistoricData.df[EcoRegNum,3] <- ProjEcoregionArea
  EcoregionHistoricData.df[EcoRegNum,4] <- BEEFirePresArea
  EcoregionHistoricData.df[EcoRegNum,5] <- BEEPercFireEcoreg
}
t2 <- Sys.time()
#############################################################################
difftime(t2,t1, units = "mins")
###
## Save data
setwd(OutDirectpsa)
write.table(EcoregionHistoricData.df, file=paste0(Species, "EcoregionHistoricFireAreaData.csv"), sep=",", col.names=NA)
####################################################################################

####################################################################################
### Use Arc Python script ENMProjectionAreaZonalStatisticsLoopCurrentFireEcoregion.py to calculate percent 
### overlap of current presence projection inside Ecoregion within NAD_1983_Albers projection
####################################################################################
## Read in above produced csv files from python and calculate percent overlaps of maxent models
## with breeding range for portions of maxent model both inside and outside breeding range
##
FinalNumberModels <- 4
ModelType <- paste0("TopAUC",FullNObs)
FinalModelVariables <- 15
NumProjections <- 30
NObsJ <- 3
NumberModSets <- 3000
#############################
j = 10
#j=1
if(j < 10) {
  SetNumber <- paste0("0", j)
} else {
  SetNumber <- as.character(j)
}
#
SetName <- paste0(SetNumber, "_", colnames(VariableSetTypes[j]))
SetNameF <- colnames(VariableSetTypes[j])
## Read in variable set data
SubsetVariableNumber <- FinalModelVariables
setwd(OutDirectpsa)
MaxentKeepEvalAreaStatsAll.dfs <- data.frame(read.csv(paste0(Species, "Maxent", DataSet, "_TrainVsTestandAreaStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", EvalType, "_TopVsRandom", NumProjections, "Sets.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
# Keep only rankings for model type
MaxentKeepEvalAreaStatsAll.df2 <- MaxentKeepEvalAreaStatsAll.dfs[which(MaxentKeepEvalAreaStatsAll.dfs$DataType==paste0("TopAUC", FullNObs, "FinalTest")),]
FinalModelRowNumberList.df <- as.data.frame(MaxentKeepEvalAreaStatsAll.df2$ModelNumber)
FinalModelRowNumberList <- FinalModelRowNumberList.df[1:FinalNumberModels,]
# Keep only final models
MaxentKeepEvalAreaStatsAll.df3 <- MaxentKeepEvalAreaStatsAll.df2[1:4,]
ncol(MaxentKeepEvalAreaStatsAll.df3)
# Read in Ecoregion Areas and Historical Percent Burn
EcoregionHistoricData.df <-  data.frame(read.csv(paste0(Species, "EcoregionHistoricFireAreaData.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
#
##############################################################################
Count <- 0
MaxentKeepEvalStatsL <- list()
t1 <- Sys.time()
##
for(FinalModelRowNumber in FinalModelRowNumberList) {
  #FinalModelRowNumber=17
  ## Keep only model data for current model
  for(EcoRegNum in 1:53) {
    #EcoRegNum<-7
    MaxentKeepEvalAreaStatsAll.df4 <- MaxentKeepEvalAreaStatsAll.df3[which(MaxentKeepEvalAreaStatsAll.df3$ModelNumber==FinalModelRowNumber),]
    ## First, read in csv file for projected current area of ecoregion
    ## Specify input directories
    output3 <- paste0("MaxentProj_", ModelType, "_Subsetof", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_ThresholdCalib", FinalModelRowNumber)
    OutDirectSub <- paste0(OutDirectpsa, "/", output3)
    ##
    setwd(OutDirectSub)
    #########################
    ## Read in csv files of maxent grid cell counts inside projection area
    ProjEcoregionCurrRangeArea.df <- data.frame(read.csv(paste0(Species2, FinalModelRowNumber, "mxt", SubsetVariableNumber,"bProjpresareatable", EcoRegNum, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE) 
    #EcoRegNum <- 1
    ProjEcoregionArea <- EcoregionHistoricData.df[EcoRegNum, 3]
    ####################################
    ## Calculate area in Megahectares of Ecoregion
    if(nrow(ProjEcoregionCurrRangeArea.df)==0) {
      ProjEcoregionProjPresArea <- 0
    } else {
      ProjEcoregionProjPresArea <- (ProjEcoregionCurrRangeArea.df$AREA * 0.0001)/ 1000000 # convert from sq. meters to hectares then megehectares
    }
    #
    ## Calculate percent of ecoregion occupied by maxent current projection
    ProjPercEcoregion <- ProjEcoregionProjPresArea/ProjEcoregionArea
    ################################
    if(EcoRegNum < 47) {
      ## Read in csv files of maxent grid cell counts inside breeding area
      BEEEcoregionCurrRangeArea.df <- data.frame(read.csv(paste0(Species2, FinalModelRowNumber, "mxt", SubsetVariableNumber,"aBEEpresareatable", EcoRegNum, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE) 
      #EcoRegNum <- 1
      BEEEcoregionArea <- EcoregionHistoricData.df[EcoRegNum, 2]
      ####################################
      ## Calculate area in Megahectares of Ecoregion
      if(nrow(BEEEcoregionCurrRangeArea.df)==0) {
        BEEEcoregionProjPresArea <- 0
      } else {
        BEEEcoregionProjPresArea <- (BEEEcoregionCurrRangeArea.df$AREA * 0.0001)/ 1000000 # convert from sq. meters to hectares then megehectares
      }
      #
      ## Calculate percent of ecoregion occupied by maxent current projection
      BEEPercEcoregion <- BEEEcoregionProjPresArea/BEEEcoregionArea
      BEEHistPercEcoregion <- EcoregionHistoricData.df[EcoRegNum,5]
      ## Calculate percent of current projection relative to historic percent of ecoregion
      DiffBEEPercProjHistEcoregion <- BEEHistPercEcoregion - BEEPercEcoregion
    } else {
      BEEEcoregionProjPresArea <- 0
      BEEPercEcoregion <- 0
      DiffBEEPercProjHistEcoregion <- 0
    }
    ##############################
    ## Add data to evaluation statistics
    ## Add Breeding Range Area in Megahectares
    MaxentKeepEvalAreaStatsAll.df4$BEEENMArea <- BEEEcoregionProjPresArea
    MaxentKeepEvalAreaStatsAll.df4$BEEPercENMArea <- BEEPercEcoregion
    MaxentKeepEvalAreaStatsAll.df4$BEEDiffPercENMCurrArea <-  DiffBEEPercProjHistEcoregion
    MaxentKeepEvalAreaStatsAll.df4$ProjENMArea <- ProjEcoregionProjPresArea
    MaxentKeepEvalAreaStatsAll.df4$ProjPercENMArea <- ProjPercEcoregion
    MaxentKeepEvalAreaStatsAll.df4$EcoRegNum <- EcoRegNum
    MaxentKeepEvalAreaStatsAll.df4$Climate <- "Current"
    MaxentKeepEvalAreaStatsAll.df4 <-MaxentKeepEvalAreaStatsAll.df4[,c(29:30,1:28)]
    ## Save output to list
    Count <- Count + 1
    MaxentKeepEvalStatsL[[Count]] <- MaxentKeepEvalAreaStatsAll.df4
  }
}
t2 <- Sys.time()
#############################################################################
difftime(t2,t1, units = "mins")
###
MaxentKeepEvalStatsAll.df <- do.call(rbind, MaxentKeepEvalStatsL)
ncol(MaxentKeepEvalStatsAll.df)
## Save data
setwd(OutDirectpsa)
write.table(MaxentKeepEvalStatsAll.df, file=paste0(Species, "EcoregionCurrentFireAreaData.csv"), sep=",", col.names=NA)
#

####################################################################################
### For each ecoregion, use Arc Python script ENMProjectionAreaZonalStatisticsLooopFireFuture4ModsEcoregion.py to calculate percent
### overlap of presence projection inside ecoregion within NAD_1983_Albers projection
####################################################################################
## Read in above produced csv files from python and calculate percent overlaps of maxent models
## with breeding range for portions of maxent model both inside and outside breeding range
##
FinalModelVariables <- 15
#TopFinalModelRowNumberList <- c(3,6,18) # Rank of model by RFI, REI, and RBI among top models
#RandomFinalModelRowNumberList <- c(1,8,13) # Rank of model by RFI, REI, and RBI among random models
NObsJ <- 3
FinalNumberModels <- 4 # Select number of top models to keep by Regional Indices
NumberModSets <- 3000
NumProjections <- 30
ModelTypeList <- c(paste0("TopAUC", FullNObs), paste0("RandAUC", FullNObs))
#############################
j = 10
#j=1
if(j < 10) {
  SetNumber <- paste0("0", j)
} else {
  SetNumber <- as.character(j)
}
#
SetName <- paste0(SetNumber, "_", colnames(VariableSetTypes[j]))
SetNameF <- colnames(VariableSetTypes[j])
## Read in variable set data
SubsetVariableNumber <- FinalModelVariables
##########
setwd(OutDirectpsa)
# Read in variable information for models
MaxentKeepEvalAreaStatsAll.dfs  <- data.frame(read.csv(paste0(Species, "Maxent", DataSet, "_TrainVsTestandAreaStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", EvalType, "_TopVsRandom", NumProjections, "Sets.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
#
MaxentKeepEvalStats.df1 <- data.frame(read.csv(paste0(Species, ModelName, FullNObs, "_", NumberModSets,"SetNumber_TrainVsTestStats_of", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
# Read in data with VarNames and Model Number Rankings
# Find variable subsets for top ranked three models by regional indices of projections for ModelType
MaxentKeepEvalAreaStatsAll.df <- data.frame(read.csv(paste0(Species, "Maxent", DataSet, "_TrainVsTestandAreaStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", EvalType, "_TopVsRandom", NumProjections, "Sets.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
head(MaxentKeepEvalAreaStatsAll.df)
# Keep only rankings for model type
MaxentKeepEvalAreaStatsAll.df2 <- MaxentKeepEvalAreaStatsAll.dfs[which(MaxentKeepEvalAreaStatsAll.dfs$DataType==paste0("TopAUC", FullNObs, "FinalTest")),]
FinalModelRowNumberList.df <- as.data.frame(MaxentKeepEvalAreaStatsAll.df2$ModelNumber)
FinalModelRowNumberList <- FinalModelRowNumberList.df[1:FinalNumberModels,]
# Keep only final models
MaxentKeepEvalAreaStatsAll.df3 <- MaxentKeepEvalAreaStatsAll.df2[1:4,]
ncol(MaxentKeepEvalAreaStatsAll.df3)
### Read in Maxent model ENM area projections for current climate
MaxentCurrentClimateProjStats.df <- data.frame(read.csv(paste0(Species, "EcoregionCurrentFireAreaData.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
## Convert NA values to zero
MaxentCurrentClimateProjStats.df[is.na(MaxentCurrentClimateProjStats.df)] <- 0
#
MaxentKeepEvalStatsL <- list()
##############################################################################
Count <- 0
#ModelTypeNames <- c(paste0("TopAUC",FullNObs), paste0("Random", FullNObs))
ModelTypeNames <- c(paste0("TopAUC",FullNObs))
t1 <- Sys.time()
###########
for(ModelType in ModelTypeNames) {
  #ModelType <- paste0("TopAUC", FullNObs)
  ################
  for(EcoRegNum in 1:53) {
    #EcoRegNum <- 46
    # Read in table of projected ecoregion area
    setwd(InDirect)
    ProjAreaTableEcoregion.df <- data.frame(read.csv(paste0("Projareatable_ecoregion", EcoRegNum, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
    ProjEcoregionArea <- (ProjAreaTableEcoregion.df$AREA * 0.0001)/ 1000000 # convert from sq. meters to hectares then megehectares
    #############
    for(FinalModelRowNumber in FinalModelRowNumberList) {
      #FinalModelRowNumber=17
       ## Keep only model data for current model
       for(FutClimSubDir in FutClimSubDirList) {
        #FutClimSubDir <- FutClimSubDirList[[1]]
        MaxentKeepEvalAreaStatsAll.df4 <- MaxentKeepEvalAreaStatsAll.df3[which(MaxentKeepEvalAreaStatsAll.df3$ModelNumber==FinalModelRowNumber),]
        ## Specify input directories
        output3 <- paste0("MaxentProj_", ModelType, "_Subsetof", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_ThresholdCalib", FinalModelRowNumber, "_", FutClimSubDir)
        OutDirectSub <- paste0(OutDirectpsa, "/", output3)
        ##
        setwd(OutDirectSub)
        ## Read in csv files of maxent grid cell counts inside Ecoregion
        EcoregionFutRangeArea.df <- data.frame(read.csv(paste0(SpecShort, FutClimSubDir, "Projpresareatable", EcoRegNum, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
        ## Calculate Megahectares area overlapping breeding range
        if(nrow(EcoregionFutRangeArea.df)==0){
          EcoregionFutProjPresArea <- 0
        } else {
          EcoregionFutProjPresArea <- (EcoregionFutRangeArea.df$AREA * 0.0001)/ 1000000 # convert from sq. meters to hectares then megehectares
        }
        ## Add data to model information
        ## Calculate percent of ecoregion occupied by maxent current projection
        FutProjPercEcoregion <- EcoregionFutProjPresArea/ProjEcoregionArea
        HistPercEcoregion <- EcoregionHistoricData.df[EcoRegNum,5]
        ## Calculate percent of current projection relative to historic percent of ecoregion
        CurrProjPercEcoregion <- MaxentCurrentClimateProjStats.df[(MaxentCurrentClimateProjStats.df$ModelNumber==FinalModelRowNumber & MaxentCurrentClimateProjStats.df$EcoRegNum==EcoRegNum),30]
        DiffPercProjFutCurrEcoregion <- FutProjPercEcoregion - CurrProjPercEcoregion
        ## Add data to evaluation statistics
        ## Add Breeding Range Area in Megahectares
        MaxentKeepEvalAreaStatsAll.df4$ProjENMArea <- EcoregionFutProjPresArea
        MaxentKeepEvalAreaStatsAll.df4$ProjPercENMArea <- FutProjPercEcoregion 
        MaxentKeepEvalAreaStatsAll.df4$ProjDiffPercFutCurrENMCurrArea <-  DiffPercProjFutCurrEcoregion
        MaxentKeepEvalAreaStatsAll.df4$EcoRegNum <- EcoRegNum
        MaxentKeepEvalAreaStatsAll.df4$Climate <- FutClimSubDir
        MaxentKeepEvalAreaStatsAll.df4 <-MaxentKeepEvalAreaStatsAll.df4[,c(27:28,1:26)]
        ## Save output to list
        Count <- Count + 1
        MaxentKeepEvalStatsL[[Count]] <- MaxentKeepEvalAreaStatsAll.df4
      }
    }
  }
}
t2 <- Sys.time()
#############################################################################
difftime(t2,t1, units = "mins")
## Save data
MaxentKeepEvalStatsL[[1]]
MaxentKeepEvalStatsFuture.df <- do.call(rbind, MaxentKeepEvalStatsL)
ncol(MaxentKeepEvalStatsFuture.df)
nrow(MaxentKeepEvalStatsFuture.df)
head(MaxentKeepEvalStatsFuture.df)
#
setwd(OutDirectpsa)
write.table(MaxentKeepEvalStatsFuture.df, file=paste0(Species, "EcoregionFutureFireAreaData.csv"), sep=",", col.names=NA)
###################################################################################

######################################################################################
## Assemble data for bar charts and ANOVA for current and future climates
##########################################################################################
setwd(OutDirectpsa)
##
NumProjections <- 30
## Read in regional index area data and combine for analysis
#FinalModelRowNumberList <- c(9,4,3)
# Read in regional indices for current models
MaxentAreaStatsCurrent.df <- data.frame(read.csv(paste0(Species, "EcoregionCurrentFireAreaData.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
head(MaxentAreaStatsCurrent.df)
# Add columns to match future data
MaxentAreaStatsCurrent.df$ProjDiffPercFutCurrENMCurrArea <- 0
# Read in regional indices for future models
MaxentAreaStatsFuture.df <- data.frame(read.csv(paste0(Species, "EcoregionFutureFireAreaData.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
head(MaxentAreaStatsFuture.df)
# Add columns to match current data
MaxentAreaStatsFuture.df$BEEENMArea <- 0
MaxentAreaStatsFuture.df$BEEPercENMArea <- 0
MaxentAreaStatsFuture.df$BEEDiffPercENMCurrArea <- 0
ncol(MaxentAreaStatsFuture.df)
# Rearrange columns to match current data
MaxentAreaStatsFuture.df <- MaxentAreaStatsFuture.df[,c(1:25, 29:31, 26:28)]
# Combine current and future regional indices together
MaxentAreaStatsClimate.df <- rbind(MaxentAreaStatsCurrent.df, MaxentAreaStatsFuture.df)
MaxentAreaStatsClimate.dfs <- arrange(MaxentAreaStatsClimate.df, Climate, DataType, ModelNumber, EcoRegNum, Climate)
head(MaxentAreaStatsClimate.dfs)
# Create ClimateEcoregionType combining Climate and EcoRegNum
MaxentAreaStatsClimate.dfs$ClimEcoRegType <- paste0(MaxentAreaStatsClimate.dfs$Climate, "_", MaxentAreaStatsClimate.dfs$EcoRegNum)
ncol(MaxentAreaStatsClimate.dfs)
# Rearrange columns
MaxentAreaStatsClimate.dfs <- MaxentAreaStatsClimate.dfs[,c(1:4,32,5:31)]
head(MaxentAreaStatsClimate.dfs)
ncol(MaxentAreaStatsClimate.dfs)
#str(MaxentAreaStatsClimate.dfs)
## Save combined data
write.table(MaxentAreaStatsClimate.dfs, file=paste0(Species, "MaxentCurrentFutureClimate", "_EcoregionAreaStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", "_TopVsRandom", NumProjections, "Sets.csv"), sep=",", col.names=NA)
###
### Calculate Mean and Standard Deviation Values
EvaluationStats <- MaxentAreaStatsClimate.dfs
nrow(EvaluationStats)
OutName <- paste0(Species, "MaxentProjection_Top4CurrFutModelEcoregionAreaStats_", SubsetVariableNumber, "Vars_", "_Summary.csv")
SortGroups <- c("ClimEcoRegType", "Climate", "EcoRegNum", "SubsetVariableNumber")
## All statistics, and only statistics, should be at and after the column specified below
StatVarFirstColumn <- 23
OutDirectIn <- OutDirectpsa
EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
head(EvaluationStatsSummary.df)
##
EvalStatsOut <- EvaluationStatsSummary.df
## Sort data by ClimDataType for ease of plotting in excel
EvalStatsOuts <- arrange(EvalStatsOut, Climate, EcoRegNum, SubsetVariableNumber, Statistic)
EvalStatsOuts <- arrange(EvalStatsOut, EcoRegNum)
head(EvalStatsOuts)
write.table(EvalStatsOuts, file=paste0(Species, "MaxentProjection_Top4CurFutModelEcoregionAreaStats_", SubsetVariableNumber, "Vars_", "_SummaryTable.csv"), sep=",", col.names=NA)

##################################################################################
### Loop through various Welch t test comparisons for comparing statistics for
### Current Climate with each of the four future scenarios for both selected and random
### models. Use Holm correction to penalize p-value for four comparisons of each type
##
StatTypeNames <- c("ProjENMArea")
#ModelTypeNames <- c(paste0("TopAUC", FullNObs), paste0("Random", FullNObs))
ModelTypeNames <- c(paste0("TopAUC",FullNObs))
ClimateTypeNames <- c("2050he26", "2050he85", "2070he26", "2070he85")
NumberComparisons <- length(ClimateTypeNames)
##
for(StatType in StatTypeNames) {
  #StatType <- "ProjArea"
  for(ModelType in ModelTypeNames) {
    #ModelType <- paste0("TopAUC", FullNObs)
    MaxentEvalStatsL <- list()
    Count <- 0
    ## Keep only ModelType data
    MaxentTestStats2 <- MaxentAreaStatsClimate.dfs[which(MaxentAreaStatsClimate.dfs$Set==ModelType),]
    head(MaxentTestStats2)
    tail(MaxentTestStats2)
    for(EcoRegNum in 1:53) {
      #EcoRegNum <- 1
      ## Keep only EcoRegNum
      MaxentTestStats1 <- MaxentTestStats2[which(MaxentTestStats2$EcoRegNum==EcoRegNum),]
      head(MaxentTestStats1)
      ## Use Welch correction one way ANOVA for unequal variances for primary test among all possible climate comparisons
      WelchANOVA <- oneway.test(as.formula(paste0(StatType, " ~ ClimEcoRegType")), data=MaxentTestStats1, na.action=na.omit, var.equal=FALSE)
      # Save above Welch's ANOVA output
      out1 <- paste0("Climate EcoregionProjection Analysis for ", ModelType, " Models of Ecoregion", EcoRegNum, ": Welch ANOVA for ", StatType)
      cat(out1, file=paste0(Species, ModelType, "_Climate_WelchANOVA_forCurrandFutEcoregion_", EcoRegNum, "_Modelsof_", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
      out2 <- capture.output(print(WelchANOVA))
      cat(out2, file=paste0(Species, ModelType, "_Climate_WelchANOVA_forCurrandFutEcoregion_", EcoRegNum, "_Modelsof_", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
      # Conduct specific pre-planned climate comparisons between pairs of current climate with each future climate
      for(ClimateType in ClimateTypeNames) {
        #ClimateType <- "2050he26"
        ## Keep only Current Climate data and ClimateType
        MaxentTestStats <- MaxentTestStats1[which(MaxentTestStats1$Climate=="Current" | MaxentTestStats1$Climate==ClimateType),]
        N <- nrow(MaxentTestStats)/2
        ## Use Welch correction t test
        RankWelchttest <- t.test(as.formula(paste0(StatType, " ~ Climate")), data=MaxentTestStats)
        #str(RankWelchttest)
        ## Apply Holm Correction to the p-value based on the total number of comparions for group (NumberComparisons)
        PValueCor <- p.adjust(RankWelchttest$p.value, method= "holm", n = NumberComparisons)
        # Save above Welch t test output
        out <- paste0("   ")
        cat(out, file=paste0(Species, "WelchttestStatsforCurrVsFutEcoregion_", EcoRegNum, "_Climates", "_",ModelType, "_", N, "Setsof", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
        out1 <- paste0("Current vs ", ClimateType, ": Welch t test of ", StatType, " for ", N, " ", ModelType, " Models of ", SubsetVariableNumber, " Variables")
         cat(out1, file=paste0(Species, "WelchttestStatsforCurrVsFutEcoregion_", EcoRegNum, "_Climates", "_",ModelType, "_", N, "Setsof", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
        out2 <- capture.output(print(RankWelchttest))
        cat(out2, file=paste0(Species, "WelchttestStatsforCurrVsFutEcoregion_", EcoRegNum, "_Climates", "_",ModelType, "_", N, "Setsof", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
        out3 <- paste0("p-value with Holm Correction for ", NumberComparisons, " Welch t tests: ", PValueCor)
        cat(out3, file=paste0(Species, "WelchttestStatsforCurrVsFutEcoregion_", EcoRegNum, "_Climates", "_",ModelType, "_", N, "Setsof", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
      }
    }
  }
}

################################################################################
### Extract variables for 2070 RCP 8.5 scenario from background points
#################################################################################
## Read in stack of 33 clipped non-climate grids
setwd(MaskGrids)
NonClimPredictors <- stack(NonClimateGridNamesL, proj4string=CRS.WGS84)
names(NonClimPredictors) <- toupper(names(NonClimPredictors))
names(NonClimPredictors) <- gsub("_NS", "", names(NonClimPredictors))
##
###########################################################
## Read in Predictors of 2070 8.5 scenario
### Specify files and directories related to future climate
FutClimDir <- paste0("F:/GIS/FutureClimate/", SpeciesGen)
## Establish subdirectory extensions for variables of future climate scenarios
FutClimSubDirList <- c("2050he26", "2050he85", "2070he26", "2070he85")
#
FutClimGridDir <- paste0(FutClimDir, "/", FutClimSubDirList[[4]])
setwd(FutClimGridDir)
ClimPredictors <- stack(FutureGridNamesL, proj4string=CRS.WGS84)
names(ClimPredictors) <- toupper(names(ClimPredictors))
names(ClimPredictors) <- gsub("_NS", "", names(ClimPredictors))
## Merge NonClimate and Future Climate grids in one stack
PredictorsAll <- stack(ClimPredictors, NonClimPredictors)
names(PredictorsAll)
## Extract data
setwd(InDirect)
# Read in background points
BackgroundPntShp <- readOGR(".", paste0(SpecFileName, "bkgrnd"))
#
system.time(BackgroundDat.mat1 <- data.frame(extract(PredictorsAll, BackgroundPntShp)))
BackgroundDat.mat1[1:10,]
# Re-associate data with species and x y coordinates
BackgroundDat.mat <- merge(coordinates(BackgroundPntShp), BackgroundDat.mat1, by="row.names")
tail(BackgroundDat.mat)
# Replace first column of row.names with LongSpecies name
colnames(BackgroundDat.mat)[1] <- c("Species")
BackgroundDat.mat$Species <- "Background"
# Convert to data frame
BackgroundDat.df <- data.frame(BackgroundDat.mat, stringsAsFactors=FALSE)
tail(BackgroundDat.df)
nrow(BackgroundDat.df)
# Omit any rows with NA values
BackgroundDat.df <- na.omit(BackgroundDat.df)
# Save extracted data
setwd(InDirect)
write.table(BackgroundDat.df, sep = ",", col.names=NA, file=paste0(Species, "_backgrounddat", FutClimSubDirList[[4]], "_", DataSet, Tag, "_", TotVars, "Vars.csv", sep=""))
##
########################## For 2050 HE GCM RCP8.5
## Read in Predictors of 2050 8.5 scenario
### Specify files and directories related to future climate
FutClimDir <- paste0("F:/GIS/FutureClimate/", SpeciesGen)
## Establish subdirectory extensions for variables of future climate scenarios
FutClimSubDirList <- c("2050he26", "2050he85", "2070he26", "2070he85")
#
FutClimGridDir <- paste0(FutClimDir, "/", FutClimSubDirList[[2]])
setwd(FutClimGridDir)
ClimPredictors <- stack(FutureGridNamesL, proj4string=CRS.WGS84)
names(ClimPredictors) <- toupper(names(ClimPredictors))
names(ClimPredictors) <- gsub("_NS", "", names(ClimPredictors))
## Merge NonClimate and Future Climate grids in one stack
PredictorsAll <- stack(ClimPredictors, NonClimPredictors)
names(PredictorsAll)
## Extract data
setwd(InDirect)
# Read in background points
BackgroundPntShp <- readOGR(".", paste0(SpecFileName, "bkgrnd"))
#
system.time(BackgroundDat.mat1 <- data.frame(extract(PredictorsAll, BackgroundPntShp)))
BackgroundDat.mat1[1:10,]
# Re-associate data with species and x y coordinates
BackgroundDat.mat <- merge(coordinates(BackgroundPntShp), BackgroundDat.mat1, by="row.names")
tail(BackgroundDat.mat)
# Replace first column of row.names with LongSpecies name
colnames(BackgroundDat.mat)[1] <- c("Species")
BackgroundDat.mat$Species <- "Background"
# Convert to data frame
BackgroundDat.df <- data.frame(BackgroundDat.mat, stringsAsFactors=FALSE)
tail(BackgroundDat.df)
nrow(BackgroundDat.df)
# Omit any rows with NA values
BackgroundDat.df <- na.omit(BackgroundDat.df)
# Save extracted data
setwd(InDirect)
write.table(BackgroundDat.df, sep = ",", col.names=NA, file=paste0(Species, "_backgrounddat", FutClimSubDirList[[2]], "_", DataSet, Tag, "_", TotVars, "Vars.csv", sep=""))
##
########################## For 2070 HE GCM RCP2.6
## Read in Predictors of 2070 2.6 scenario
### Specify files and directories related to future climate
FutClimDir <- paste0("F:/GIS/FutureClimate/", SpeciesGen)
## Establish subdirectory extensions for variables of future climate scenarios
FutClimSubDirList <- c("2050he26", "2050he85", "2070he26", "2070he85")
#
FutClimGridDir <- paste0(FutClimDir, "/", FutClimSubDirList[[3]])
setwd(FutClimGridDir)
ClimPredictors <- stack(FutureGridNamesL, proj4string=CRS.WGS84)
names(ClimPredictors) <- toupper(names(ClimPredictors))
names(ClimPredictors) <- gsub("_NS", "", names(ClimPredictors))
## Merge NonClimate and Future Climate grids in one stack
PredictorsAll <- stack(ClimPredictors, NonClimPredictors)
names(PredictorsAll)
## Extract data
setwd(InDirect)
# Read in background points
BackgroundPntShp <- readOGR(".", paste0(SpecFileName, "bkgrnd"))
#
system.time(BackgroundDat.mat1 <- data.frame(extract(PredictorsAll, BackgroundPntShp)))
BackgroundDat.mat1[1:10,]
# Re-associate data with species and x y coordinates
BackgroundDat.mat <- merge(coordinates(BackgroundPntShp), BackgroundDat.mat1, by="row.names")
tail(BackgroundDat.mat)
# Replace first column of row.names with LongSpecies name
colnames(BackgroundDat.mat)[1] <- c("Species")
BackgroundDat.mat$Species <- "Background"
# Convert to data frame
BackgroundDat.df <- data.frame(BackgroundDat.mat, stringsAsFactors=FALSE)
tail(BackgroundDat.df)
nrow(BackgroundDat.df)
# Omit any rows with NA values
BackgroundDat.df <- na.omit(BackgroundDat.df)
# Save extracted data
setwd(InDirect)
write.table(BackgroundDat.df, sep = ",", col.names=NA, file=paste0(Species, "_backgrounddat", FutClimSubDirList[[3]], "_", DataSet, Tag, "_", TotVars, "Vars.csv", sep=""))
##
########################################################################################################



####################################################################################
### Calculate the MaxEnt model area for the entire Projection Area for
### comparison among variable set sizes
### Use Arc Python script ENMProjectionAreaZonalStatisticsLoopFfireCurrentProjectionArea.py to
### calculate areas within NAD_1983_Albers projection
####################################################################################
## Read in above produced csv files from python for areas within projection area
##
#############################################################################
setwd(InDirect)
MaxEntProjectionAreaData.df <- as.data.frame(t(as.matrix(c(Species,0,0,0,0))), stringsAsFactors=FALSE)
colnames(MaxEntProjectionAreaData.df) <- c("Species", "AreaType", "SubsetVariableNumber", "ModelNum", "MxtProjArea")
MaxEntProjectionAreaData.df[,3:5] <- as.numeric(as.character(MaxEntProjectionAreaData.df[,3:5]))
Count <- 0
ModelTypeList <- c("TopAUC250")
#NumVarList <- c(6,10,15)
NumVarList <- c(15)
#ModelNumberList <- seq(1,9,1)
ModelNumberList <- seq(1,9,1)
AreaTypeList <- c("BEE", "ProjNoBEE")
#
t1 <- Sys.time()
for(ModelType in ModelTypeList){
  for(AreaType in AreaTypeList) {
    for(NumVar in NumVarList)  {
      for(Model in ModelNumberList){
        ## Read in csv file with area of MaxEnt model in projection area
        MaxEntProjectionArea.df <- data.frame(read.csv(paste0(Species, Model, "Mxt", NumVar, "bProjPresAreaTable", AreaType, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
        ## Calculate Megahectares presence area within Ecoregion
        if(nrow(MaxEntProjectionArea.df)==0) {
          MxtProjresArea <- 0
        } else {
          MxtProjPresArea <- (MaxEntProjectionArea.df$AREA * 0.0001)/ 1000000 # convert from sq. meters to hectares then megehectares
        }
        ### Save data to dataframe
        Count <- Count + 1
        MaxEntProjectionAreaData.df[Count,1] <- Species
        MaxEntProjectionAreaData.df[Count,2] <- AreaType
        MaxEntProjectionAreaData.df[Count,3] <- NumVar
        MaxEntProjectionAreaData.df[Count,4] <- Model
        MaxEntProjectionAreaData.df[Count,5] <- as.numeric(MxtProjPresArea)
      }
    }
  }
}
t2 <- Sys.time()
#############################################################################
difftime(t2,t1, units = "mins")
###
## Save data
setwd(OutDirectpsa)
#write.table(MaxEntProjectionAreaData.df, file=paste0(Species, "MaxEntProjectionAreaData.csv"), sep=",", col.names=NA)
write.table(MaxEntProjectionAreaData.df, file=paste0(Species, "MaxEnt15VarProjectionAreaData.csv"), sep=",", col.names=NA)
##########
###Calculate Mean and Standard Deviation Values
#
EvaluationStats <- MaxEntProjectionAreaData.df
ncol(EvaluationStats)
EvaluationStats$Model <- "MaxEnt"
EvaluationStats <- EvaluationStats[, colnames(EvaluationStats)[c(6,1:5)]]
##
#FSAType <- RankStatistic
#OutName <- paste0(Species, "Mxt", "ProjectionAreaSummaryStats.csv")
OutName <- paste0(Species, "Mxt15Var", "ProjectionAreaSummaryStats.csv")
SortGroups <- c("AreaType", "Model", "Species", "SubsetVariableNumber")
## All statistics, and only statistics, should be at and after the column specified below
StatVarFirstColumn <- 6
OutDirectIn <- OutDirectpsa
EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
head(EvaluationStatsSummary.df)

####################################################################################

##################################################################################
### Loop through various Welch t test comparisons for comparing statistics for
### Current Climate with each of the four future scenarios for both selected and random
### models. Use Holm correction to penalize p-value for four comparisons of each type
##
StatTypeNames <- c("RFI", "REI", "RBI")
#ModelTypeNames <- c(paste0("TopAUC", FullNObs), paste0("Random", FullNObs))
ModelTypeNames <- c(paste0("TopAUC",FullNObs))
ClimateTypeNames <- c("2050he26", "2050he85", "2070he26", "2070he85")
NumberComparisons <- length(ClimateTypeNames)
##
for(StatType in StatTypeNames) {
  #StatType <- "RFI"
  for(ModelType in ModelTypeNames) {
    #ModelType <- paste0("TopAUC", FullNObs)
    MaxentEvalStatsL <- list()
    Count <- 0
    ## Keep only ModelType data
    MaxentTestStats1 <- MaxentAreaStatsClimate.dfs[which(MaxentAreaStatsClimate.dfs$Set==ModelType),]
    ## Use Welch correction one way ANOVA for unequal variances for primary test among all possible climate comparisons
    WelchANOVA <- oneway.test(as.formula(paste0(StatType, " ~ ClimDataType")), data=MaxentTestStats1, na.action=na.omit, var.equal=FALSE)
    # Save above Welch's ANOVA output
    out1 <- paste0("Climate Projection Analysis for ", ModelType, " Models: Welch ANOVA for ", StatType)
    cat(out1, file=paste0(Species, ModelType, "_Climate_WelchANOVA_forCurrentandFutureModelsof_", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
    out2 <- capture.output(print(WelchANOVA))
    cat(out2, file=paste0(Species, ModelType, "_Climate_WelchANOVA_forCurrentandFutureModelsof_", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
    # Conduct specific pre-planned climate comparisons between pairs of current climate with each future climate
    for(ClimateType in ClimateTypeNames) {
      #ClimateType <- "2050he26"
      ## Keep only Current Climate data and ClimateType
      MaxentTestStats <- MaxentTestStats1[which(MaxentTestStats1$Climate=="Current" | MaxentTestStats1$Climate==ClimateType),]
      N <- nrow(MaxentTestStats)/2
      ## Use Welch correction t test
      RankWelchttest <- t.test(as.formula(paste0(StatType, " ~ Climate")), data=MaxentTestStats)
      #str(RankWelchttest)
      ## Apply Holm Correction to the p-value based on the total number of comparions for group (NumberComparisons)
      PValueCor <- p.adjust(RankWelchttest$p.value, method= "holm", n = NumberComparisons)
      # Save above Welch t test output
      out <- paste0("   ")
      cat(out, file=paste0(Species, "WelchttestStatsforCurrentVsFutureClimates", "_",ModelType, "_", N, "Setsof", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
      out1 <- paste0("Current vs ", ClimateType, ": Welch t test of ", StatType, " for ", N, " ", ModelType, " Models of ", SubsetVariableNumber, " Variables")
      cat(out1, file=paste0(Species, "WelchttestStatsforCurrentVsFutureClimates", "_",ModelType, "_", N, "Setsof", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
      out2 <- capture.output(print(RankWelchttest))
      cat(out2, file=paste0(Species, "WelchttestStatsforCurrentVsFutureClimates", "_",ModelType, "_", N, "Setsof", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
      out3 <- paste0("p-value with Holm Correction for ", NumberComparisons, " Welch t tests: ", PValueCor)
      cat(out3, file=paste0(Species, "WelchttestStatsforCurrentVsFutureClimates", "_",ModelType, "_", N, "Setsof", SubsetVariableNumber, "Variables_", StatType, ".txt"), sep="\n", append=TRUE)
    }
  }
}



