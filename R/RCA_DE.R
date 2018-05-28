######################################################################
##                  Find Differntially express gene                 ##
######################################################################
"RCA_DE"
#' Find differentially express gene by different methods
#' 
#' @param data Input data matrix. Column is gene and row is cells.
#' @param DEGmethod Methods to find DE genes.
#' 
#' @export
RCA_DE <- function(data, DEGmethod = "seurat")
{
  ## Check Method
  metho <- c("seurat", "Test")
  if(DEGmethod %in% metho == FALSE ) stop("Please check available method", call. = FALSE)
  
  ## Apply SEURAT
  if(DEGmethod=="seurat")
  {
    myres<- RCA_DE_SEURAT (data = data )
    return(myres)
  }
  
  ## Apply other methods : coming soon
  
}

################################ SEURAT
RCA_DE_SEURAT <- function(data)
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




