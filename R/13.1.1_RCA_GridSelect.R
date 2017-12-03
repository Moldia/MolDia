######################################################################
##                        Select Grid of interest                   ##
######################################################################
"RCA_GridSelect"
#' Select grid of interest from a tissue.
#' @description Select grid of interest from a tissue
#' @param data Input data in class RCA_class. Output of \link[MolDia]{readRCA}.
#' @param grid_id Grid to select. Default id NULL.
#' @param nx,ny Numbers of rectangular quadrats in the x and y directions
#' @param gridtype type of grid to plot. Default is "hexa"
#' 
#' @examples
#' ex_data <- readRCA(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),cellid = "CellId", centX = "centroid_x", centY = "centroid_y", rpc = 3)
#' 
#' ex_data <- readRCA(file = "D:\\Project Sweden\\SU Project\\Working code\\Data\\All tissue\\CellBlobs2_1.csv",
#'                    cellid = "CellID", centX = "centroidX", centY = "centroidY", rpc =1)
#'  
#' mygrid  <- RCA_GridSelect(data = ex_data, nx = 5,gridtype = "hexa")
#' mygrid  <- RCA_GridSelect(data = ex_data, nx = 8,gridtype = "rect", grid_id = c(6,16,20,8,17,21))
#' 
#' 
#' @export
RCA_GridSelect <- function(data, nx = 6, ny = nx, gridtype = "hexa", grid_id = NULL)
{
  ## Save main data
  main_data <- data
  
  ## Main data
  mydata <- data@location
  
  ## Find the convex point
  hpts <- grDevices::chull(mydata)
  hpts <- mydata[hpts, ]
  
  ## Sort points by anto clockwise
  anti_hpts <- contoureR::orderPoints(x = hpts$centroid_x, y = hpts$centroid_y, clockwise = FALSE)
  hpts <- hpts[anti_hpts,]
  
  ## Spatial point area
  myspa <- spatstat::ppp(x = mydata$centroid_x, y = mydata$centroid_y, 
                         poly=list(x=hpts$centroid_x, y=hpts$centroid_y), check = FALSE)
  plot(myspa, main = "Grid on Tissue")
  
  ## Divides window into quadrats and counts the numbers
  if(gridtype == "rect")
    {
    myspa1 <- spatstat::quadrats(myspa, nx = nx, ny = ny)
    #names(myspa1$tiles) <- paste("Tile",c(1:length(myspa1$tiles)))
    names(myspa1$tiles) <- c(1:length(myspa1$tiles))
    }
  
  if(gridtype == "hexa")
  {
    max_range <- max(apply(apply(mydata,2,range),2,diff))
    myspa1    <- spatstat::hextess(myspa, s = abs(max_range)/nx)
    names(myspa1$tiles) <- c(1:length(myspa1$tiles))
  }
  plot(myspa1, add= TRUE,col= "red",do.labels=TRUE, labelargs = list(col = "red"))
  
  ## Select grid of interest
  if(length(grid_id) > 0 ){ 
  
  pp <- myspa1
  tt <- pp$tiles
  tile_id <- paste0("Grid_id_", grid_id)
  for(i in 1:length(seq_along(grid_id)))
  {
    nam <- tile_id[i]
    assign(nam, tt[[grid_id[i]]])
  }
  eval(parse(text = paste("kk1 <- spatstat::union.owin(",paste0(tile_id,collapse=", "),")")))
  isin <- spatstat::inside.owin(x = mydata$centroid_x, y = mydata$centroid_y,w=kk1)
  
  ## Plot grid of interest
  point_in <- ex_data@location[isin,]
  #point_in1 <- spatstat::ppp(x = point_in$centroid_x, y = point_in$centroid_y,window = kk1, check = FALSE)
  
  ## Return result
  final_data <-  main_data@data[rownames(point_in), , drop = FALSE]
  final_data <- final_data[,colSums(final_data)>0, drop = FALSE]
  main_data@data <- final_data
  if(length(main_data@norm.data)  > 0 ) main_data@norm.data  <- main_data@norm.data[rownames(point_in), , drop = FALSE]
  if(length(main_data@scale.data) > 0 ) main_data@scale.data <- main_data@scale.data[rownames(point_in), drop = FALSE] 
  main_data@location <- main_data@location[rownames(point_in),]
  main_data@gene <- colnames(final_data)
  
  RCA_map(main_data, main = "Selected grid from Tissue")
  }
  
  return(main_data)
}

