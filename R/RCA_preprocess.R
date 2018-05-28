######################################################################
##                       RCA Data preprocess                        ##
######################################################################
"RCA_preprocess"
#' Pre-process RCA data.
#' @description Pre-process RCA data interms of normalization, scalling and centering
#' @param data Input data in class MolDiaISS. Output of \link[MolDia]{readISS}
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
#' @return Return a object in MolDiaISS with value in slot norm.data or scale.data
#'
#' @references
#' Van den Berg, R. A., Hoefsloot, H. C., Westerhuis, J. A., Smilde, A. K., & van der Werf, M. J. (2006). Centering, scaling, and
#' transformations: improving the biological information content of metabolomics data. BMC Genomics, 7, 142.
#' http://doi.org/10.1186/1471-2164-7-142
#'
#' @examples
#' mydata  <- readISS(file = system.file("extdata", "CellBlobs_QT_0.35.csv", package="MolDia"),
#'                   cellid = "CellID", centX = "centroidX", centY = "centroidY")
#' mydata_process <- RCA_preprocess(data_1, normalization.method = "QuantileNormalize", do.scale = FALSE,do.center = FALSE)
#'
#'
#'
#'
#' @export
RCA_preprocess <- function(data, normalization.method = "LogNormalize",
                           scale.factor = 10000, do.scale = FALSE,
                           do.center = FALSE,display.progress = FALSE)
{
  ## Save main data
  main_data <- data
  
  ## No normalization
  if(length(normalization.method) == 0 ) data@norm.data <- as.matrix(data@data)
  
  ## Different type normalization
  if(length(normalization.method) == 1 )
  {
    ## Check normalization method is in category
    method_type <- c("LogNormalize","RankNormalize", "QuantileNormalize")
    if(normalization.method %in% method_type == FALSE) stop("Please check available methods for normalization", call. = FALSE)
    
    ## Log normalized
    if(normalization.method == "LogNormalize")
    {
      res       <- Seurat::CreateSeuratObject(raw.data = t(data@data), normalization.method = normalization.method,
                                              project = "RCA_data_amalysis", min.cells = 0, min.genes = 0, is.expr = 0,
                                              scale.factor = scale.factor, do.scale = do.scale,
                                              do.center = do.center,display.progress = display.progress)
      
      data@norm.data <- t(as.matrix(res@data))
    }
    
    ## Ranked normalized
    if(normalization.method == "RankNormalize")
    {
      p  <-  main_data@data
      p1 <- lapply(apply(p,1,list),unlist)
      p1 <- lapply(p1,function(y) data.table::frank(x= y,ties.method = "dense")) #/ length(y))
      p1 <- do.call(rbind,p1)
      colnames(p1) <- colnames(p)
      
      data@norm.data <- as.matrix(p1)
    }
    
    ## Quantile normalized
    if(normalization.method == "QuantileNormalize")
    {
      p  <- as.matrix(main_data@data)
      p1<- t(preprocessCore::normalize.quantiles(t(p),copy=TRUE))
      colnames(p1) <- colnames(p)
      rownames(p1) <- rownames(p)
      
      data@norm.data <- as.matrix(p1)
    }
  }
  
  ## Scale or center data 
  if(any(c(do.scale,do.center) == TRUE))
    {
     #res   <- Seurat::CreateSeuratObject(raw.data = t(data@norm.data))
     res   <- Seurat::CreateSeuratObject(raw.data = t(data@data))
     res   <- Seurat::ScaleData(object = res,do.scale = do.scale,do.center = do.center, check.for.norm = FALSE)
     data@scale.data <- as.matrix(t(res@scale.data))
    }  
  
  ## Return result
  return(data)
}


