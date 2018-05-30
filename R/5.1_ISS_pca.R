######################################################################
##                            ISS PCA                               ##
######################################################################
"ISS_pca"
#'
#'Principle component analyis on ISS data
#'
#' @description Principle component analysis on any data in class MolDiaISS.
#' @param data Input data in class MolDiaISS. Output of \link[MolDia]{readISS}.
#' @param pc Desired percent of variance to be explained by PCA. Default is 1 which means 100 percent variation explained.
#' @param DEGmethod Methods for finding differentially expressed (DE) genes.
#' 
#' @examples 
#' ## Reading data
#' left_hypo <- readISS(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),
#'                      cellid = "CellId", centX = "centroid_x", centY = "centroid_y")
#'                   
#' ## Arrange marker gene
#' data(marker_gene)
#' mark_gene <- list(genr = marker_gene$genr, neuron = c(marker_gene$genr_neuro,
#'                                                       marker_gene$genr_neuro_pyra1,
#'                                                       marker_gene$genr_neuro_pyra2,
#'                                                       marker_gene$genr_neuro_inter1,
#'                                                       marker_gene$genr_neuro_inter2,
#'                                                       marker_gene$genr_neuro_inter3,
#'                                                       marker_gene$genr_neuro_inter4,
#'                                                       marker_gene$genr_neuro_inter5,
#'                                                       marker_gene$genr_neuro_inter6),
#'                                                       nonneuron = marker_gene$genr_nonneuro)
#' 
#' ## Barplot of Neuronal marker gene and extract those cells only
#' neuron_group <- ISS_barplot(data = left_hypo, gene = mark_gene, gene.target = 2,
#'                             at.least.gene = 2, gene.show = 2)
#' 
#' ## Preprocess data
#' neuron_group <- ISS_preprocess(data = neuron_group, normalization.method = "LogNormalize",
#'                                do.scale = TRUE, do.center = TRUE)
#'                                
#' ## Apply principle component analysis
#' res <- ISS_pca(data = neuron_group, pc = 0.9)
#' 
#' @export
ISS_pca <- function(data, pc = 1, DEGmethod=NULL)
{
  ## Save main data
  main_data <- data
  
  ## Check class of data
  if(class(data)%in%c("MolDiaISS") ==FALSE) stop("Check input data in class 'RCA_class'",call. = FALSE)
  
  # Create Seurat object
  SEURAT_clus   <- Seurat::CreateSeuratObject(t(main_data@data))
  #SEURAT_clus  <- Seurat::ScaleData(object = SEURAT_clus, do.scale = TRUE, do.center = TRUE, check.for.norm = FALSE)
  if(length(main_data@scale.data)== 0 ) stop("Please scale data first with ISS_preprocess function ", call. = FALSE)
  SEURAT_clus@scale.data <- t(main_data@scale.data) # t(main_data@scale.data[rownames(main_data),])
  
  ## Check pc
  if(pc <= 0   | pc > 1) stop("Desired percent of variance explained by PCA should be between 0 and 1 (0 percent to 100 percent)", call. = FALSE)
  
  ## Find Differential express gene
  if (length(DEGmethod) == 0) {
    gene_de <- colnames(main_data@data)
  } else {
    gene_de <- ISS_deg(data, DEGmethod)
  }
  
  if(length(gene_de) == 0 )
  {
    message("These is no differentially expressed gene. So all gene will be used for PCA")
    gene_de <- NULL
  }
  
  ## Find Optimal PC
    npc   <- ncol(data@data) - 1
    npc   <- withCallingHandlers(suppressWarnings(irlba::prcomp_irlba(x = t(data@scale.data), n=npc, 
                                                                      fastpath = TRUE, verbose = FALSE)))
    kk<- npc
    npc   <- summary(npc)$importance[3,]
    pcuse <- max(which(npc<=pc))
    cat("Principle component: ", round(npc[pcuse]*100,2), "% of variation has explained by" ,pcuse, "principle component", "\n")
    pc    <- 1:pcuse

  ## Stopping criteria for principle component
    if(pcuse <= 3) stop("Number of principle component should be greater than 3. Please increase the value of pc", call. = FALSE)  
    
  ## Run PCA
    SEURAT_clus   <- withCallingHandlers(suppressWarnings(Seurat::RunPCA(object  = SEURAT_clus, pc.genes = gene_de, pcs.compute = length(pc),
                                                                         do.print = FALSE, fastpath=TRUE, verbose = FALSE)))
    
  ## save data
    main_data@pca.data <- SEURAT_clus@dr
    
  ## Return data
    return(main_data)
}
