######################################################################
##                       RCA define cluster marker                  ##
######################################################################
"RCA_marker"
#' Find Cluster marker of in-situ RCA data.
#' @description Find Cluster marker of in-situ RCA data and plot significant gene cluster wise by barplot and heatmap.
#'
#' @param data Input data in class RCA_class. Output of \link[MolDia]{readISS}.
#' @param topgene Desired number of top gene ineach cluster to show in summary result.
#' @param test.use Denotes which test to use. Seurat currently implements "bimod"
#'        (likelihood-ratio test for single cell gene expression, McDavid et al., Bioinformatics, 2013, default), "roc"
#'        (standard AUC classifier), "t" (Students t-test), and "tobit" (Tobit-test for differential gene expression,
#'        as in Trapnell et al., Nature Biotech, 2014). For details see \link[Seurat]{FindAllMarkers}
#' @param only.pos Only return positive markers (TRUE by default)
#' @param marker.sig Lower value will identify less significant marker. Default is 0.005
#' @param main Title of the plot.
#'
#' @return A list of cluster with putative ranked markers and associated statistics in slot cluster.marker of main data.
#'         Also clusterwise figure in barplot and heatmap of desired top genes.
#'
#' @examples
#' ## Reading data
#' data_3 <- readISS(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),
#'                   cellid = "CellId", centX = "centroid_x", centY = "centroid_y")
#'
#' ## Arrange marker gene
#' data(marker_gene)
#' marker_gene <- marker_gene
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
#' neuron_group <- RCA_barplot(data = data_3, gene = mark_gene, gene.target = 2,
#'                             at.least.gene = 2, gene.show = 2)
#' ## Data preprocessing
#' neuron_group <- RCA_preprocess(data = neuron_group, normalization.method = "LogNormalize",
#'                                do.scale = TRUE, do.center = TRUE)
#'
#' ## Cluster data based on SEURAT pipeline
#' neuron_group_clust  <- RCA_cluster (data = neuron_group, method = "seurat",
#'                                     pc = 0.9, resolution = 0.1)
#' ## Get cluster marker
#' neuron_group_clust_marker <- RCA_marker(data = neuron_group_clust, topgene =15,
#'                                         test.use="bimod", main = "")
#'                                         
#'  RCA_map (data=neuron_group_clust_marker, what = "cluster")                                       
#'
#' @export
RCA_marker <- function(data, topgene= 5, test.use = "bimod", marker.sig = 0.005,
                       only.pos = TRUE, main = "")
{
  # Main data
  main_data <- data
  # Check if data is clustered or not
  if (length(data@cluster) == 0) stop("Please cluster your data befor finding marker gene of cluster", call. = FALSE)

  # Create SEURAT object
  RCA_Seurat_obj            <- Seurat::CreateSeuratObject(t(main_data@data))
  RCA_Seurat_obj@scale.data <- t(main_data@scale.data)
  RCA_Seurat_obj@ident      <- main_data@cluster

  ## Find cluster marker
  set.seed(12345)
  RCA_markers   <- Seurat::FindAllMarkers(object = RCA_Seurat_obj, only.pos = only.pos,
                                          test.use = test.use ,thresh.use = marker.sig,
                                          min.pct = marker.sig, min.diff.pct = marker.sig,  min.cells.gene = 3, min.cells.group = 3)
  RCA_markers_1 <- lapply(split(RCA_markers,RCA_markers$cluster), function(i){i[order( abs(i$pct.1), decreasing = TRUE),]})
  res           <- RCA_markers_1
  RCA_markers_1 <- lapply(RCA_markers_1, function(i) utils::head(i,topgene))
  RCA_markers_1 <- unique(unlist(lapply(RCA_markers_1, "[[", "gene")))

  ## List of cells name by each cluster
  data_identity <- data.frame(main_data@cluster)
  colnames(data_identity)<- "cluster"
  mmdata <- split(data_identity, data_identity)
  mmdata_2 <- mmdata

  ## Get reads/cell on each cluster only for found marker gene
  mmdata <- lapply(mmdata, function(object)
  {
    mmdata <- merge(object, main_data@data, by = "row.names")
    rownames(mmdata) <-mmdata$Row.names
    mmdata$Row.names <- NULL
    mmdata$cluster   <- NULL
    mmdata <- colSums(mmdata)[RCA_markers_1] / nrow(mmdata)
    mmdata
  })
  mmdata_1 <- do.call(cbind, mmdata)
  mmdata   <- log2(as.matrix(mmdata_1+1))
  dfp1     <- reshape::melt(mmdata)
  dfp1$X1  <- factor(dfp1$X1, levels = RCA_markers_1)

  p <- ggplot2::ggplot(dfp1, ggplot2::aes_string(x = "X2", y= "value", fill = "X1")) +
    ggplot2::scale_x_discrete(limits= (dfp1$X2)) +
    ggplot2::geom_bar(stat="identity",position = "dodge") +
    ggplot2::labs(x = "Cell clusters", y = "log(Reads per cell+1)", title= paste0(main, "-All common marker gene (Barplot)") ) +
    ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position="right") +
    ggplot2::geom_vline(xintercept = 0.5:(length(levels(data@cluster)) -1.5), colour="grey") +
    ggplot2::coord_cartesian(xlim = c(0, length(unique(dfp1$X2))-1))

  #### Heatmap of cluster based on marker gene
  p1 <- Seurat::DoHeatmap(object = RCA_Seurat_obj, genes.use = RCA_markers_1, group.by = "ident",
                  slim.col.label = TRUE, remove.key = TRUE, rotate.key = FALSE, use.scaled = TRUE, title = paste0(main, "-All common marker gene (Heatmap)"))

  ## Multiple plot
  print(p)
  print(p1)

  # Return result
  res <- lapply(res, function(object)
    {
     myres <- object[,c("cluster","gene")]
     myres
    })
  main_data@cluster.marker <- res
  return(main_data)
  #mmdata_1
}
