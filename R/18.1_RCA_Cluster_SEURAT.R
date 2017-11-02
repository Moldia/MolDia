######################################################################
##                  RCA define cluster: SEURAT                      ##
######################################################################
#' "RCA_seruat_cluster"
#' Cluster in-situ RCA data by SEURAT.
#'
#' @param data Input data in class RCA_class. Output of \link[MolDia]{readRCA}.
#' @param pc Desired percent of variance to be explained by PCA. Default in 0.9.
#' @param cluster_id Re-cluster clustreded data. Numeric input. Default is NULL.
#' @param resolution Value of the resolution parameter, use a value above (below)
#'                   1.0 if you want to obtain a larger (smaller) number of communities.
#'                   Default is 0.3.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm;
#'        2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm). Default is 1.
#' @keywords internal
#'
#' @examples
#' ## Reading data
#' data_3 <- readRCA(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),
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
#' neuron_group <- RCA_barplot(data = data_3, gene = mark_gene, gene.target = 2,
#'                             at.least.gene = 8, gene.show = 2)
#' ## Data preprocessing
#' neuron_group <- RCA_preprocess(data = neuron_group, normalization.method = "LogNormalize",
#'                                do.scale = TRUE, do.center = TRUE)
#'
#' ## Cluster data based on SEURAT pipeline
#' # neuron_group_clust  <- RCA_seruat_cluster (data = neuron_group, pc = 0.9, resolution = 0.4)
#'
#' ## Re-cluster specific cluster
#' # re_clust  <- RCA_seruat_cluster (data = neuron_group_clust, pc = 0.9,
#' #                                  cluster_id = 0, resolution = 0.5)
#'
#'
RCA_seruat_cluster <- function(data, pc = 0.9, cluster_id = NULL,
                            resolution = 0.3, algorithm = 1,
                            do.norm = TRUE, do.scale = TRUE)
{
  ## Save main data
  main_data <- data
  # Check class of data
  if(class(data)%in%c("RCA_class") ==FALSE) stop("Check input data in class 'RCA_class'",call. = FALSE)

  if (length(cluster_id) == 0 )
    {
      ## Select data
      data_2 <- data@data
      data   <- data_2
    }

    if (length(cluster_id) > 0 )
    {
      if((class(cluster_id) == "numeric") == FALSE) stop("cluster_id value should be numeric", call. = FALSE)
      if((class(cluster_id) == "numeric") == TRUE && length(data@cluster) == 0 ) stop("Please cluster data first and then select desired cluste to re-clusterr", call. = FALSE)
      if((class(cluster_id) == "numeric") == TRUE && length(data@cluster) > 0)
      {
        ## Select cell in desired cluster
        data_1 <- data@cluster
        data_1 <- names(data_1[which(data_1 == as.character(cluster_id))])
        data_2 <- data@data
        data   <- data_2[data_1,]
      }
    }

    # Create Seurat object
    SEURAT_clus   <- Seurat::CreateSeuratObject(t(data))
    #SEURAT_clus   <- Seurat::ScaleData(object = SEURAT_clus, do.scale = TRUE, do.center = TRUE, check.for.norm = FALSE)
    SEURAT_clus@scale.data <- t(main_data@scale.data[rownames(data),])
    SEURAT_clus   <- Seurat::RunPCA(object  = SEURAT_clus, pc.genes = colnames(data),do.print = FALSE)

    # Find optimal PCA component
    #myPCA  <- SEURAT_clus@dr$pca@cell.embeddings
    sdev   <- SEURAT_clus@dr$pca@sdev
    pcuse  <- cumsum(log2(sdev)^2 / sum(log2(sdev)^2))
    print(pcuse)
    pcuse  <- max(which(pcuse<=pc))
    if(pcuse <= 3) stop("PC is too low", call. = FALSE)
    cat("Number of principle component to be used :", pcuse, "\n")

    ## Find cluster
    cat ("Running data clustering..\n")
    SEURAT_clus   <- Seurat::FindClusters(object = SEURAT_clus, reduction.type = "pca", dims.use = 1:pcuse,
                                          plot.SNN = FALSE, print.output = 0, save.SNN = T,n.iter = 10,
                                          algorithm = algorithm, resolution = resolution, modularity.fxn = 1)
    ## return RCA object
    newdata <- as.data.frame(t(SEURAT_clus@raw.data))
    #res <- methods::new("RCA_class",
    #                      data  = data,
    #                      location = main_data@location[rownames(newdata),],
    #                      gene = colnames(newdata),
    #                      cluster = SEURAT_clus@ident)

    main_data@data <- data
    if(length(main_data@norm.data)  > 0 ) main_data@norm.data  <- main_data@norm.data[rownames(newdata),colnames(newdata)]
    if(length(main_data@scale.data) > 0 ) main_data@scale.data <- main_data@scale.data[rownames(newdata),colnames(newdata)]
    main_data@location <- main_data@location[rownames(newdata),]
    main_data@gene <- colnames(newdata)
    main_data@cluster <- SEURAT_clus@ident

    return(main_data)
}



