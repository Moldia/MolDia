######################################################################
##                        Select Grid of interest                   ##
######################################################################
"ISS_roi"
#' Select grid/region of interest (ROI) from a tissue.
#' @description Select grid of interest from a tissue
#' @param data Input data in class MolDiaISS. Output of \link[MolDia]{readISS}.
#' @param gridtype type of grid to plot. Default is "rect". See details.
#' @param grid_id Grid to select. Default id NULL.
#' @param nx,ny Numbers of rectangular quadrats in the x and y directions
#' @param roifile Name of the file that contain ROI. csv formate. Most important is, data points should be in clockwise or
#'                anti-clockwise order.
#' @param roi.id Column name in roifile file that contain ROI id
#' @param roi.x,roi.y X and Y axis name in roifile file.
#  @param clockwise Order of data points in supplied roifile. 
#' @param main Main title.
#' 
#' @details gridtype parameter can have 4 values :  "hexa", "rect", "roifile" and "roi"
#'          
#'          1. "hexa" will create hexagonal grid and number of grid is based on nx parameter and grid_id will select grid of interest.
#'          
#'          
#'          2. "rect" will create rectangular grid and number of grid is based on nx parameter and grid_id will select grid of interest.
#'          
#'          3. "roifile" will require a file input in csv formate in the parameter "roifile". The input file must have at least 3 column of
#'          ROI id, x-axis, y-axis. The parameter "roi.id", "roi.x" and "roi.y" will be the name input of corresponding column.
#'          
#'          4. "roi" will select ROI more interactively. One cal select ROI on the image by pointer. 
#' 
#' @examples
#' # Read ISS data
#' ex_data <- readISS(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),cellid = "CellId", centX = "centroid_x", centY = "centroid_y", rpc = 3)
#' 
#' # Select hexagonal and rectangular ROI
#' mygrid  <- ISS_roi(data = ex_data, nx = 8,gridtype = "hexa", main = "Hexagonal ROI on tissue")
#' mygrid  <- ISS_roi(data = ex_data, nx = 8,gridtype = "rect", main = "Rectangular ROI on tissue")
#' mygrid  <- ISS_roi(data = ex_data, nx = 8,gridtype = "rect", grid_id = c(6,16,20,8,17,21))
#' 
#' # Selected ROI from a pre selected ROI file
#' ex_data <- readISS(file = system.file("extdata", "CellBlobs_ROI.csv", package="MolDia"),
#'                    cellid = "CellID", centX = "centroidX", centY = "centroidY")
#' mygrid  <- ISS_roi(data = ex_data, nx = 6, gridtype = "roifile",
#'                             roifile = system.file("extdata", "polygon_coordinates.csv", package="MolDia"),
#'                             roi.id = "Polygon_id", roi.x ="x_coordiates" , roi.y = "y_coordinates", grid_id = c(1,2,5,6),
#'                             main = "Selected ROI on tissue")
#' # Select ROI interactively by user
#' mygrid  <- ISS_roi(data = ex_data, gridtype = "roi")
#' 
#' @export
ISS_roi <- function(data, gridtype = "rect", nx = 6, ny = nx, grid_id = NULL, 
                    roifile = NULL, roi.id = NULL, roi.x = NULL, roi.y = NULL, # clockwise=TRUE,  
                    main = "ROI on tissue")
{
  ## Save main data
  main_data <- data
  
  ## Main data
  mydata <- data@location
  
  ## Find the convex point
  hpts <- grDevices::chull(mydata)
  hpts <- mydata[hpts, ]
  
  ## Sort points by anti clockwise
  anti_hpts <- contoureR::orderPoints(x = hpts$centroid_x, y = hpts$centroid_y, clockwise = FALSE)
  hpts <- hpts[anti_hpts,]
  
  ## Spatial point area
  myspa <- spatstat::ppp(x = mydata$centroid_x, y = mydata$centroid_y, 
                         poly=list(x=hpts$centroid_x, y=hpts$centroid_y), check = FALSE)
  plot(myspa, main = main)
  
  ## Divides window into quadrats and counts the numbers
  if(gridtype == "rect")
  {
    myspa1 <- spatstat::quadrats(myspa, nx = nx, ny = ny)
    #names(myspa1$tiles) <- paste("Tile",c(1:length(myspa1$tiles)))
    names(myspa1$tiles) <- c(1:length(myspa1$tiles))
    plot(myspa1, add= TRUE,col= "red",do.labels=TRUE, labelargs = list(col = "red"))
  }
  
  if(gridtype == "hexa")
  {
    max_range <- max(apply(apply(mydata,2,range),2,diff))
    myspa1    <- spatstat::hextess(myspa, s = abs(max_range)/nx)
    names(myspa1$tiles) <- c(1:length(myspa1$tiles))
    plot(myspa1, add= TRUE,col= "red",do.labels=TRUE, labelargs = list(col = "red"))
  }
  
  if(gridtype == "roifile")
  {
    if(length(roifile)==0) stop("Please select location of ROI file in csv formate",call. = TRUE)
    ## Read ROI file in CSV formate
    roi <- read.csv(file = roifile)
    
    ## Split ROI file by region id 
    roi <- split(roi, roi[,roi.id])
    
    ## Define clockwise and anto clockwise rotation
    anticlock <- function(x,y) 
      {
        N<- length(x)
        ts <- 0
        for(i in 1: (N-1))
          { ts<- ts + ((x[i]-x[i+1])*(y[i]+y[i+1])) }
        return(ts)
      }
    
    roi <- lapply(roi, function(i)
    {
      mm  <- anticlock(i[,roi.x],i[,roi.y])
      if(sign(mm)==-1) i <- i[seq(dim(i)[1],1),]
      i
    })
    
    
    roi <- lapply(roi, function(i)
    { 
      hpts <- i #subset(i, select=-c(Polygon.id))
      hpts[,roi.id] <- NULL
      #anti_hpts <- contoureR::orderPoints(x = hpts[,roi.x], y = hpts[,roi.y], clockwise = FALSE)
      #hpts <- hpts[anti_hpts,]
      
      #if(clockwise == TRUE){ hpts <- hpts[rev(1:nrow(hpts)),] }
      #if(clockwise == FALSE){ hpts <- hpts }
      
      hpts <- as.matrix(hpts)
      hpts <- apply(hpts,2,list)
      hpts <- lapply(hpts, unlist)
      hpts <- spatstat::owin(poly =  list(x = hpts[[roi.x]],y = hpts[[roi.y]] ))
      #hpts <- spatstat::owin (poly = hpts)
      hpts
    }
    )
    
    ## Ploting selected region
    # myspa1 <-lapply(roi, plot, add=T, border = "red", lwd = 2)
    mycol  <- palette(rainbow(length(roi)))
    myspa1 <- lapply(seq_along(roi), function (i){plot(roi[[i]], add=T, border = mycol[i], lwd = 3)} )
    
    ## label each polygon
    labs<-names(roi)
    for (i in 1: length(roi)) {xy<-spatstat::centroid.owin((poly = roi[[i]]));
    text(xy$x,xy$y, labels = labs[i], col = mycol[i], cex = 1.5)}
    
    ## Add legend
    legend("topleft", legend= paste0("ROI ", 1:length(roi)), #text.col = mycol,
           col= mycol, lty=1, lwd = 4, cex=1)
    
    ## Convert window object to Tessellation
    myspa1 <-spatstat::as.tess(roi)
  }
  
  if(gridtype == "roi")
  {
    myspa1 <- spatstat::clickpoly(add = TRUE, col = 2, lwd = 2)
    myspa1 <- spatstat::as.tess(myspa1)
  }
  
  
  ## Select grid of interest
  if(length(grid_id) > 0 | gridtype == "roi")
    { 
    
    if(gridtype != "roi")
    {
      pp <- myspa1
      tt <- pp$tiles
      tile_id <- paste0("Grid_id_", grid_id)
      for(i in 1:length(seq_along(grid_id)))
      { 
        nam <- tile_id[i]
        assign(nam, tt[[grid_id[i]]])
      }
      eval(parse(text = paste("kk1 <- spatstat::union.owin(",paste0(tile_id,collapse=", "),")")))
    }
    
    if(gridtype == "roi") kk1 <- myspa1$tiles[[1]]
    isin <- spatstat::inside.owin(x = mydata$centroid_x, y = mydata$centroid_y,w=kk1)
    
    ## Plot grid of interest
    point_in <- main_data@location[isin,]
    #point_in1 <- spatstat::ppp(x = point_in$centroid_x, y = point_in$centroid_y,window = kk1, check = FALSE)
    
    ## Return result
    final_data <-  main_data@data[rownames(point_in), , drop = FALSE]
    final_data <- final_data[,colSums(final_data)>0, drop = FALSE]
    main_data@data <- final_data
    if(length(main_data@norm.data)  > 0 ) main_data@norm.data  <- main_data@norm.data[rownames(point_in), , drop = FALSE]
    if(length(main_data@scale.data) > 0 ) main_data@scale.data <- main_data@scale.data[rownames(point_in), drop = FALSE] 
    main_data@location <- main_data@location[rownames(point_in),]
    main_data@gene <- colnames(final_data)
    
    ISS_map(main_data, main = main)
  }
  return(main_data)
}