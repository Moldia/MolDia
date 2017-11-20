######################################################################
##                         RCA spatial cluster                      ##
######################################################################
"RCA_spatial"
#' Spatial clustering for ISS.
#' 
#' @description ISS spatial clustering
#' 
#' @examples 
#' ## Reading data
#' ex_data <- readRCA(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),
#'                    cellid = "CellId", centX = "centroid_x", centY = "centroid_y")
#' myres <- RCA_spatial(ex_data, neighbor.points = 4)
#'
#'@export
RCA_spatial <- function(data, neighbor.points = 9, spatial.info = "phisical")
{
  ######### Nearest neighbour cell
  ## Physical distance
  if(spatial.info == "phisical")
    {
    res1 <- fields::rdist(data@location)
    colnames(res1) <- rownames(data@location)
    rownames(res1) <- rownames(data@location)
    res1
    }
  
  ## Genetic distance
  if(spatial.info == "genetic")
    {
    res1 <- as.matrix(stats::dist(data@data))
    }
  
  res1           <- apply(res1,1,function(i)
    {
      kk<- sort(i, index.return = TRUE)$ix [2:(neighbor.points+1)]
    kk
    })
  res1 <- t(res1)
  rownames(res1) <- match(rownames(data@location),rownames(res1))
  rn <- apply(res1, 1, length)
  rownames(res1) <- rn
  res1 <- cbind(rn,res1)
  res1
  
  ## Binary gene expression
  bin_gene <- 1*(data@data > 0)
  rownames(bin_gene) <- match(rownames(data@location),rownames(bin_gene))
  rn <- rownames(bin_gene)
  bin_gene <- cbind(rn, bin_gene)
  bin_gene
  
  ## Return result
  res <- list(bin_gene,res1)
  return(res)
  
}

