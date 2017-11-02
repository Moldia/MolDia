######################################################################
##                       RCA Data preprocess                        ##
######################################################################
"RCA_scale"
#' Pre-process RCA data.
#' @description Pre-process RCA data interms of normalization, scalling and centering
#' @param data Input data in class RCA_class. Output of \link[MolDia]{readRCA}
#' @param normalization.method Method for normalization. Default is log-normalization (LogNormalize).
#'        Use NULL if one dont apply normalization. More methods to be added very shortly.
#' @param scale.factor Sets the scale factor for cell-level normalization.
#' @param do.scale Whether to scale the data. See details.
#' @param do.center Whether to center the data. See details.
#' @param display.progress Displays a progress bar
#'
#' @details Setting center to TRUE will center the expression for each gene by subtracting the average expression for that gene.
#'          Setting scale to TRUE will scale the expression level for each gene by dividing the centered gene expression levels
#'          by their standard deviations if center is TRUE and by their root mean square otherwise.
#'
#' @return Return a object in RCA_class with value in slot norm.data or scale.data
#'
#' @references
#' Van den Berg, R. A., Hoefsloot, H. C., Westerhuis, J. A., Smilde, A. K., & van der Werf, M. J. (2006). Centering, scaling, and
#' transformations: improving the biological information content of metabolomics data. BMC Genomics, 7, 142.
#' http://doi.org/10.1186/1471-2164-7-142
#'
#' @examples
#' data_1 <- readRCA(file = system.file("extdata", "CellBlobs_QT_0.35.csv", package="MolDia"),
#'                   cellid = "CellID")
#' scale_data_1 <- RCA_preprocess(data_1)
#'
#' @export
RCA_preprocess <- function(data, normalization.method = "LogNormalize",
                           scale.factor = 10000, do.scale = TRUE,
                           do.center = TRUE,display.progress = FALSE)
{
  ## Save main data
  main_data <- data

  ## Create SEURAT object
  res       <- Seurat::CreateSeuratObject(raw.data = t(data@data), normalization.method = normalization.method,
                                                     project = "RCA_data_amalysis", min.cells = 0, min.genes = 0, is.expr = 0,
                                                     scale.factor = scale.factor, do.scale = do.scale,
                                                     do.center = do.center,display.progress = display.progress)

  ## Repalce required values
  if(length(normalization.method) == 0)   {data@norm.data  <- data@data}
  else {data@norm.data <- t(as.matrix(res@data))}
  if(any(c(do.scale,do.center) == TRUE))  data@scale.data <- as.matrix(t(res@scale.data))

  ## Return result
  return(data)
}
