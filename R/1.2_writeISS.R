######################################################################
##                       Write RCA ISS  data                         ##
######################################################################
"writeISS"
#' Write ISS data in CSV formate
#' @description Write ISS data in CSV formate.
#' 
#' @param file File name in "MolDiaISS" class (Output of \link[MolDia]{readISS}).
#' @param output A character string naming a file.
#' @param expr Write expression data. TRUE or FALSE needed. Default is TRUE.
#' @param location Write location data. TRUE or FALSE needed. Default is TRUE.
#'        
#' @author Mohammad Tanvir Ahamed
#' @return Output will be a CSV formated file.
#' 
#' @examples 
#' 
#' ##### Reading ISS data from CSV formate
#' data_1 <- readISS(file = system.file("extdata", "CellBlobs_QT_0.35.csv", package="MolDia"),
#'                   cellid = "CellID", centX = "centroidX", centY = "centroidY")
#' writeISS(file = data_1, output  = "test.csv", expr = F, location = T)
#' @export
 writeISS <- function(file, output = "", expr = TRUE, location = TRUE)
 {
   ## Check data class
   if(class(file) !="MolDiaISS" ) stop("Please input data in class MolDiaISS", call. = FALSE)
   
   ## Check if all data is available
   if(expr == FALSE && location == FALSE) stop("Please selecet data to write in csv file", call. = FALSE)
   
   ## Read expression data
   if(expr == TRUE) expr_data <- file@data
   else(expr_data <- NA)
   
   ## Read location data
   if (location ==TRUE) loca_data <- file@location
   else (loca_data <- NA)
   
   ## Combine expression data and location data
   final_data <- cbind(expr_data,loca_data)
   final_data<- final_data[ , ! apply( final_data , 2 , function(x) all(is.na(x)) ) ]

   ## Write CSV data
   write.csv(x = final_data, file = output)
      
 } 