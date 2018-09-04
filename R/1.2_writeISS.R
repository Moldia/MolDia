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
#' @param cellid String to naming cell. Default is "CellID".
#' @param centX Name of X co-ordinate in file. Default is "centroidX".
#' @param centY Name of Y co-ordinate in file  Default is "centroidY".
#' @param scientific Fomate of location data. Either a logical specifying 
#'        whether elements of a real or complex vector should be encoded in scientific format, or an integer. Default is FALSE.
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
 writeISS <- function(file, output = "", expr = TRUE, location = TRUE, cellid = "CellID", 
                      centX = "centroidX", centY = "centroidY", scientific = FALSE)
 {
   ## Check data class
   if(class(file) !="MolDiaISS" ) stop("Please input data in class MolDiaISS", call. = FALSE)
   
   ## Check if all data is available
   if(expr == FALSE && location == FALSE) stop("Please selecet data to write in csv file", call. = FALSE)
   
   ## Read expression data
   if(expr == TRUE) expr_data <- file@data
   else(expr_data <- NA)
   
   ## Read location data
   if (location ==TRUE) 
    {
     loca_data <- file@location
     names(loca_data)<- c(centX, centY)
     loca_data <- as.matrix(apply(loca_data, 2, function(i) format(i, scientific = scientific)))
    }
   else (loca_data <- NA)
   
   ## Combine expression data and location data
   final_data <- cbind(expr_data,loca_data)
   final_data <- final_data[ , ! apply( final_data , 2 , function(x) all(is.na(x)) ) ]
   
   ## Renaming rowname
   rn <- rownames(final_data)
   rn <- as.numeric(sapply(strsplit(rn,"_"),"[[",2))
   rn <- data.frame(rn)
   names(rn) <- cellid
   final_data <- cbind(rn, final_data)
   rownames(final_data) <- NULL

   ## Write CSV data
   write.csv(x = final_data, file = output, row.names = FALSE)
      
 } 