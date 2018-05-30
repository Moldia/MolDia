######################################################################
##                      ISS cluster validation                      ##
######################################################################
"ISS_tsne"
#'
#'Dimentionality reduction to 2D by tSNE
#'
#' @description Any data in class MolDiaISS clusteded or not clustered used to reduce dimention to 2D by RCA-tsne.
#' @param data Input data in class MolDiaISS. Output of \link[MolDia]{readISS}.
#' @param pc Desired percent of variance to be explained by PCA. Default is 1 which means 100 percent variation explained.
#' @param perplexity Numeric; Perplexity parameter. See  \link[Rtsne]{Rtsne}
#' 
#' 
#' @return 2D dataframe of points in slot @tsne.data.
#'
#' @examples
#' ## Reading data
#' left_hypo <- readISS(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),
#'                   cellid = "CellId", centX = "centroid_x", centY = "centroid_y")
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
#'                                             nonneuron = marker_gene$genr_nonneuro)
#'
#' ## Barplot of Neuronal marker gene and extract those cells only
#' neuron_group <- ISS_barplot(data = left_hypo, gene = mark_gene, gene.target = 2,
#'                             at.least.gene = 2, gene.show = 2)
#' 
#' ## Data preprocessing
#' neuron_group <- ISS_preprocess(data = neuron_group, normalization.method = "LogNormalize",
#'                                do.scale = TRUE, do.center = TRUE)
#'                                
#' #### Dimention reduction by tSNE on non-clustered data
#' tsne_noclust <- ISS_tsne(data = neuron_group, pc = 0.7)
#' # Plot tSNE
#' result <- ISS_map(data = tsne_noclust, what = "tsneAll")
#' # Plot selected gene on tSNE plot
#' result <- ISS_map(data = tsne_noclust, what = "tsne", gene =tsne_noclust@gene[1:2] )
#' 
#' 
#' #### Dimention reduction by tSNE on clustered data
#' # Cluster data based on SEURAT pipeline
#' neuron_group_clust  <- ISS_cluster(data = neuron_group, pc = 0.7, resolution = 0.3, method = "seurat")
#' # Dimention reduction by tSNE
#' tsne_clust   <- ISS_tsne(data = neuron_group_clust, pc= 0.9, perplexity= 100)
#' Plot cluster on tSNE plot
#' result <- ISS_map(data = tsne_clust, what = "tsneAll")
#' 
#' @export
ISS_tsne <- function(data, pc = 1, perplexity = 30)
{
  # Check class of data
  if(class(data)%in%c("MolDiaISS") ==FALSE) stop("Check input data in class 'MolDiaISS'", call. = FALSE)

  # Save main data 
  mdata <- data
  
  # Run tSNE on normalize / non normalized data
  data  <- mdata@data
  
  # Run PCA and get optimum pc
  res <- ISS_pca(data = mdata, pc = pc)
  pc <- 1:ncol(res@pca.data$pca@cell.embeddings)

  # Create SEURAT object
  RCAtsne   <- Seurat::CreateSeuratObject(t(data))
  #RCAtsne   <- Seurat::ScaleData(object = RCAtsne, do.scale = TRUE, do.center = TRUE, check.for.norm = FALSE)
  if(length(mdata@scale.data)== 0 ) stop("Please scale data first with ISS_preprocess function ", call. = FALSE)
  RCAtsne@scale.data <- t(mdata@scale.data[rownames(data),])
  RCAtsne@dr         <- res@pca.data # Seurat::RunPCA(object  = RCAtsne, pc.genes = colnames(data),do.print = FALSE)
  
  # Running tSNE
  set.seed(12345)
  cat("tSNE: Running dimentionality reduction by tSNE..\n")
  RCAtsne   <- Seurat::RunTSNE(object = RCAtsne , dims.use = pc, do.fast = TRUE, check_duplicates = FALSE, perplexity = perplexity )
  cat("      Finished running tSNE.\n")

  # Assign cluster information on tSNE plot
  #if(length(mdata@cluster)> 0 ) RCAtsne@ident <- mdata@cluster
  #RCAtsne_1 <- Seurat::TSNEPlot(RCAtsne, do.label = do.label, label.size = 8)
  RCAtsne   <- RCAtsne@dr$tsne@cell.embeddings

  # Assign tSNE data to slot
  mdata@tsne.data <- data.frame(RCAtsne)
  return(mdata)
}
