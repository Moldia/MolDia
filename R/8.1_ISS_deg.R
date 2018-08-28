######################################################################
##                  Find Differntially express gene                 ##
######################################################################
"ISS_deg"
#' Find differentially express (DE) gene by different methods
#' 
#' @param data Input data matrix. Column is gene and row is cells.
#' @param DEGmethod Methods to find DE genes.
#' 
#' @note This function is in-complete and needs further development by adding 
#'       different methods of finding DE optimized for ISS data.
#' 
#' @examples 
#' ## Reading data
#' ex_data <- readISS(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),
#'                    cellid = "CellId", centX = "centroid_x", centY = "centroid_y")
#' 
#' res <- ISS_deg (data = ex_data)
#' 
#' @export
ISS_deg <- function(data, DEGmethod = "seurat")
{
  ## Check Method
  metho <- c("seurat", "Test")
  if(DEGmethod %in% metho == FALSE ) stop("Please check available method", call. = FALSE)
  
  ## Apply SEURAT
  if(DEGmethod=="seurat")
  {
    myres<- ISS_deg_seurat (data = data )
    return(myres)
  }
  
  ## Apply other methods : coming soon
  
}

################################ SEURAT
ISS_deg_seurat <- function(data)
{
  mdata     <- data
  DE_SEURAT <- Seurat::CreateSeuratObject(t(mdata))
  DE_SEURAT <- Seurat::FindVariableGenes(object = DE_SEURAT, do.plot = FALSE,display.progress = FALSE)
  DE_SEURAT <- Seurat::FindVariableGenes(object = DE_SEURAT, do.plot = FALSE,display.progress = FALSE,
                                     x.low.cutoff = 0, x.high.cutoff = max(DE_SEURAT@hvg.info$gene.dispersion), 
                                     y.cutoff = 0, do.contour=F)
  
  DE_gene   <- DE_SEURAT@var.genes
  if(length(DE_gene) == 0) stop("There is no differentially express gene", call. = FALSE)
  return(DE_gene)
}




