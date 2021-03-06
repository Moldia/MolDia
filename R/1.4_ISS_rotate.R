######################################################################
##                       Rotate ISS data                            ##
######################################################################
"ISS_rotate"
#' Rotate and flip ISS data
#'
#' @description Rotate by angle and flip by axix of ISS data
#' @param data Input data in class MolDiaISS. Output of \link[MolDia]{readISS}.
#' @param rot Rotation angle (Degree). Default is zero. Negative angle define the anti-clock wise rotation.
#' @param flipby Flip the data points by x or y axis. Default is NULL. Possible value is NULL, "x", "y".  
#'
#' @examples
#' # Read ISS data
#' mydata  <- readISS(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),
#'                    cellid = "CellId", centX = "centroid_x", centY = "centroid_y")
#' 
#' # Flip/rotate ISS data
#' mydata_x  <- ISS_rotate (data = mydata, rot = 0, flipby = "x")
#' mydata_y  <- ISS_rotate (data = mydata, rot = 0, flipby = "y")
#' mydata_90 <- ISS_rotate (data = mydata, rot = 90)
#' 
#' # Plot ISS data
#' res     <- ISS_map(data = mydata, what = "cell", main = "Left hippocampus")
#' res1    <- ISS_map(data = mydata_x, what = "cell", main = "Left hippocampus")
#' res1    <- ISS_map(data = mydata_y, what = "cell", main = "Left hippocampus")
#' res1    <- ISS_map(data = mydata_90, what = "cell", main = "Left hippocampus")
#'
#' @export
ISS_rotate <- function(data, rot = 0, flipby = NULL)
{
  ## Location data
  loc_data <- data@location 
  loc_names<- colnames(loc_data)
  
  ## Rotation angle
  theta <- - rot*pi/180 
  
  ## Rotate data by specific angle
  loc_data <- spdep::Rotation(xy = loc_data,angle = theta)
  
  ## Flip data
  if(length(flipby)== 0 ) loc_data <- loc_data
  if(length(flipby) > 0)
  {
    if(flipby == "x") 
    {
      loc_data_y <- abs(loc_data[,2]-max(loc_data[,2]))
      loc_data <- cbind(loc_data[,1],loc_data_y)
    }
    if(flipby == "y") 
    {
      loc_data_x <- abs(loc_data[,1]-max(loc_data[,1]))
      loc_data <- cbind(loc_data_x, loc_data[,2])
    }
  }
  
  ## Re-name location
  colnames(loc_data)<- loc_names
  
  ## Restort location information
  data@location <- data.frame(loc_data)
  
  ## Return data 
  return(data)
  
}