"ISS_lineage"
#'
#'
#'
#'
#' @note Need to develop writeISS function to write back csv file from the S4 object
#'
#' @examples 
#' adult <- readISS(file = system.file("extdata", "lineage_adult.csv", package="MolDia"),  
#'                  cellid = "CellID",centX = "centroidX", centY = "centroidY", 
#'                  gene = c("Actb.L", "Adar1", "Adar2","Aldh1l1", "Blcap.edQR","Blcap.edYC" ))
#' p7    <- readISS(file = system.file("extdata", "linrage_p07.csv", package="MolDia"),  
#'                  cellid = "CellID",centX = "centroidX", centY = "centroidY")
#' p0    <- readISS(file = system.file("extdata", "lineage_p00.csv", package="MolDia"),  
#'                  cellid = "CellID",centX = "centroidX", centY = "centroidY")
#' e15   <- readISS(file = system.file("extdata", "lineage_e15.csv", package="MolDia"),  
#'                  cellid = "CellID",centX = "centroidX", centY = "centroidY")
#'                  
#'  ISS_map (data = adult, what = "cell")
#'  ISS_map (data = p7, what = "cell")
#'  ISS_map (data = p0, what = "cell")
#'  ISS_map (data = e15, what = "cell")



