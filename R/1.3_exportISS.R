######################################################################
##                       Export RCA ISS  data                       ##
######################################################################
"exportISS"
#' Export ISS data in different formated file.
#' @description Export ISS data to different formated file which are compitible for other software or R package.
#' 
#' @param file Input data in class MolDiaISS. Output of \link[MolDia]{readISS}. See details.
#' @param toexport To which formate data will export. See details. 
#' 
#' @details At present only non-segmentated file (file class \link[MolDia]{MolDiaISS_nonsegment})is supported for 
#'          exporting file formate compitible with R package wholebrainsoftware.
#'          
#'          toexport only support "WB" argument at this moment, which is compitible with  wholebrainsoftware.
#'          
#' @examples 
#' ###### Reading non-segmentated file
#' mdata   <- readISS(file = system.file("extdata", "nonSeg_QT_0.35_0_details.csv", package="MolDia"), segment = FALSE,
#'                  centX = "PosX", centY = "PosY", nogene = "NNNN", gene= c("Gdf7","WNT1","Pak3","Tfap2a"))
#' 
#' mdata_1 <- readISS(file = mdata, gene= c("Gdf7","WNT1","Pak3","Tfap2a"))
#' mdata_2 <- readISS(file = mdata, nogene = c("Gdf7","WNT1"))
#' 
#' mywb_1  <- exportISS(file = mdata_1, toexport= "WB")
#' mywb_2  <- exportISS(file = mdata_2, toexport= "WB")
#' 
#' @export
 exportISS <- function(file, toexport= "WB")
 {
   if(toexport == "WB")
   {
     ## Check if data is in correct class
     if(class(file)%in%c("MolDiaISS_nonsegment") == FALSE ) stop("Please check input data is in class MolDiaISS_nonsegment", call. = FALSE)
     
     ## Define x and y
     x <- file@location$centroid_x
     y <- file@location$centroid_y
   
     ## Define intensity, area
     intensity <- file@quality$MinQuality
     area      <- file@quality$MinAnchor
   
     ## Define counter.x, counter.y and counter.ID
     contour.x  <- numeric(0)
     contour.y  <- numeric(0)
     contour.ID <- NULL 
     
     ## Return result
     
     res <- list(x=x,y=y,intensity=intensity, area=area, contour.x=contour.x, contour.y=contour.y, contour.ID=contour.ID)
     res <- list(res)
     names(res) <- "soma"
   }
   return(res)
 }
 