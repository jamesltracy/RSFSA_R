HeatMapDataFrame.Plot <-
function(DataFrameName, PlotName, PlotTitle, OutDirectIn..., PlotWidth, PlotHeight, PointSize, XAxisPntProp, YAxisPntProp, Angle)  {
  #NOTE: Function parameters after the "...", such as PlotWidth, have to be set with an equal sign in the function call, such as "PlotWidth=10"
  # Set default value for plot options if left out of function call
  if(missing(PlotWidth)) { PlotWidth=10 }
  if(missing(PlotHeight)) { PlotHeight=7 }
  if(missing(PointSize)) { PointSize=8 }
  if(missing(XAxisPntProp)) { XAxisPntProp=0.7 }
  if(missing(YAxisPntProp)) { YAxisPntProp=0.7 }
  if(missing(Angle)) { Angle=300 }
  ## Generate heat map
  library(ggplot2)
  library(plyr)
  library(reshape2)
  library(scales)
  # Reformat data
  ## Convert to matrix and sort rows by diagonal values in matrix
  DiagMatrix.df <- data.frame(diag(as.matrix(DataFrameName)))
  # Sort the diagonal
  DiagMatrix.dfs <- as.data.frame(DiagMatrix.df[order(-1*DiagMatrix.df[,1]),,drop=FALSE])
  # Name original data frame
  DataMatrix.df <- DataFrameName
  head(DataMatrix.df)
  # Sort entire data frame rows then columns by rownames of DiagMatrix.dfs
  # First sort rows
  DataMatrix.dfs1 <- DataMatrix.df[order(match(rownames(DataMatrix.df), rownames(DiagMatrix.dfs))),]
  head(DataMatrix.dfs1)
  # Then sort columns
  DataMatrix.dfs <- DataMatrix.dfs1[, order(match(colnames(DataMatrix.dfs1), rownames(DiagMatrix.dfs)))]
  head(DataMatrix.dfs)
  tail(DataMatrix.dfs)
  #
  DataMatrixt.df <- data.frame(t(as.matrix(DataMatrix.dfs)))
  DataMatrix <- data.frame(rownames(DataMatrixt.df))
  rownames(DataMatrix) <- rownames(DataMatrixt.df)
  colnames(DataMatrix) <- "Classes"
  DataMatrix <- cbind(DataMatrix, DataMatrixt.df)
  ## Generate heat map
  # Melt data for plotting
  DataMatrixm <- melt(DataMatrix)
  #str(DataMatrixm)
  # Make order of factor levels of $Classes in melt product above the same as in original matrix
  DataMatrixm$Classes = with(DataMatrixm, factor(Classes, levels = rownames(DataMatrix)))
  # Calculate midpoint of values for color palette
  middle <- (max(as.matrix(DataMatrix.df))+ min(as.matrix(DataMatrix.df)))/2
  ##
  library(RGtk2)
  library(Cairo)
  library(cairoDevice)
  ### NOTE: May not see anything on plot output, but check directory where tiff saved by ggsave below
  # Open up Cairo device for higher quality plot, may need to adjust width and height based on plot
  Cairo(width = PlotWidth, height = PlotHeight, pointsize = PointSize, surface = c("screen", "png", "pdf", "ps", "svg"), filename = NULL)
  #
  DataMatrixHeatPlot <- ggplot(DataMatrixm, aes(variable, Classes)) + geom_tile(aes(fill = value),
       colour = "white") + scale_fill_gradient2(low="#3333FF", mid= "#FFFF00", high="#FF3333", midpoint=middle)
  #DataMatrixHeatPlot <- ggplot(DataMatrixm, aes(variable, Classes)) + geom_tile(aes(fill = value),
  #     colour = "white") + scale_fill_gradientn(colours = rainbow(7))
  # Clean up heat map
  base_size <- PointSize
  ## For black label colors plotting
  DataMatrixHeatPlot <- DataMatrixHeatPlot + theme_grey(base_size = base_size) + labs(x = "",
       y = "") + scale_x_discrete(expand = c(0, 0)) + labs(title = paste(PlotTitle)) +
       scale_y_discrete(expand = c(0, 0)) + theme(axis.ticks = element_blank(),
       axis.text.x = element_text(size = base_size * XAxisPntProp, angle = Angle, hjust = 0, colour = "black", face = "bold"),
       axis.text.y = element_text(size = base_size * YAxisPntProp, colour = "black", face = "bold"), title = element_text(face = "bold")) + coord_fixed()
  #
  # Need to check how plot looks and possibly adults width and height above on Cairo statement
  setwd(OutDirectIn)
  ggsave(paste(PlotName))
  dev.off()
}
