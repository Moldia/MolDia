######################################################################
##                        Select Grid of interest                   ##
######################################################################
"RCA_GridSelect"
#' Select grid of interest from a tissue.
#' @description Select grid of interest from a tissue
#' @param data Input data in class RCA_class. Output of \link[MolDia]{readRCA}.
#' @param grid_id Grid to select
#' 
#' @examples 
#' ex_data <- readRCA(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),cellid = "CellId", centX = "centroid_x", centY = "centroid_y", rpc = 3)
#' ex_data <- readRCA(file = "D:\\Project Sweden\\SU Project\\Working code\\Data\\All tissue\\CellBlobs2_1.csv",
#'                    cellid = "CellID", centX = "centroidX", centY = "centroidY", rpc =1)
#' 
#' 
#' @export
RCA_GridSelect <- function(data,grid_id)
{
  ## Main data
  mydata <- data@location
  
  ## Find the convex point
  hpts <- grDevices::chull(mydata)
  hpts <- mydata[hpts, ]
  
  ## Sort points by anto clockwise
  anti_hpts <- contoureR::orderPoints(x = hpts$centroid_x, y = hpts$centroid_y, clockwise = FALSE)
  hpts <- hpts[anti_hpts,]
  
  ## Spatial point area
  myspa <- ppp(x = mydata$centroid_x, y = mydata$centroid_y, poly=list(x=hpts$centroid_x, y=hpts$centroid_y))
  plot(myspa)
  
}