######################################################################
##                  RCA define cluster: SEURAT                      ##
######################################################################
#' "ISS_cluster_seruat"
#' Cluster ISS data by SEURAT.
#'
#' @param data Input data in class MolDiaISS. Output of \link[MolDia]{readISS}.
#' @param pc Desired percent of variance (0 to 1) to be explained by PCA. Default in NULL (All PC will use).
#' @param cluster_id Re-cluster clustreded data. Numeric input. Default is NULL.
#' @param resolution Value of the resolution parameter, use a value above (below)
#'                   1.0 if you want to obtain a larger (smaller) number of communities.
#'                   Default is 0.3.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm;
#'        2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm). Default is 1.
#' @param DEGmethod Methods to find DE genes.
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @keywords internal
#'
#' @examples
#' ## Reading data
#' data_3 <- readISS(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),
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
#' neuron_group <- ISS_barplot(data = data_3, gene = mark_gene, gene.target = 2,
#'                             at.least.gene = 8, gene.show = 2)
#' ## Data preprocessing
#' neuron_group <- ISS_preprocess(data = neuron_group, normalization.method = "LogNormalize",
#'                                do.scale = TRUE, do.center = TRUE)
#'
#' ## Cluster data based on SEURAT pipeline
#' neuron_group_clust  <- ISS_cluster_seruat (data = neuron_group, pc = 0.9, resolution = 0.4)
#'
#' ## Re-cluster specific cluster
#' # re_clust  <- ISS_cluster_seruat (data = neuron_group_clust, pc = 0.9,
#' #                                  cluster_id = 0, resolution = 0.5)
#'
ISS_cluster_seruat <- function(data, pc = NULL, cluster_id = NULL,
                            resolution = 0.3, algorithm = 1,
                            # do.norm = TRUE, do.scale = TRUE,
                            DEGmethod = "seurat", k.param = 30)
{
  ## Save main data
  main_data <- data
  # Check class of data
  if(class(data)%in%c("MolDiaISS") ==FALSE) stop("Check input data in class 'MolDiaISS'",call. = FALSE)

  if (length(cluster_id) == 0 )
    {
      ## Select data
      data_2 <- data@norm.data
      data   <- data_2
      
      # Run PCA and get optimum pc
      res <- ISS_pca(data = main_data, pc = pc)
      pc <- 1:ncol(res@pca.data$pca@cell.embeddings)
      pcs.compute <- length(pc)
      
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
      
      # Select data in cluster
      selected_data <- ISS_clustSelect(data = main_data , cluster_id = cluster_id)
      
      # Run PCA and get optimum pc
      res <- ISS_pca(data = selected_data, pc = pc)
      pc <- 1:ncol(res@pca.data$pca@cell.embeddings)
      pcs.compute <- length(pc)
    }
  
    # Create Seurat object
    SEURAT_clus   <- Seurat::CreateSeuratObject(t(data))
    #SEURAT_clus   <- Seurat::ScaleData(object = SEURAT_clus, do.scale = TRUE, do.center = TRUE, check.for.norm = FALSE)
    if(length(main_data@scale.data)== 0 ) stop("Please scale data first with ISS_preprocess function ", call. = FALSE)
    SEURAT_clus@scale.data <- t(main_data@scale.data[rownames(data),])
    SEURAT_clus@dr         <- res@pca.data
    
    # Find Differential express gene
    if(length(DEGmethod) == 0) 
      {
      gene_de <- colnames(data)
      }
    else
      {
      gene_de <- ISS_deg(data, DEGmethod)
      }
    
    
    ## Find cluster
    cat ("Data clustering method: SEURAT\n")
    cat ("Running data clustering.....")
    SEURAT_clus   <- suppressMessages(Seurat::FindClusters(object = SEURAT_clus, genes.use = gene_de, dims.use = 1:pcs.compute, 
                                                           algorithm = algorithm, resolution = resolution, k.param = k.param,
                                                           reduction.type = "pca", plot.SNN = FALSE, print.output = FALSE, 
                                                           save.SNN = FALSE, n.iter = 10, modularity.fxn = 1,
                                                           temp.file.location = tempfile() ))
    cat ("FINISHED \n")
    
    ## return RCA object
    newdata <- as.data.frame(t(SEURAT_clus@raw.data))
    main_data@data <- main_data@data[rownames(newdata),colnames(newdata)] #data.frame(newdata)
    if(length(main_data@norm.data)  > 0 ) main_data@norm.data  <- main_data@norm.data[rownames(newdata),colnames(newdata)]
    if(length(main_data@scale.data) > 0 ) main_data@scale.data <- main_data@scale.data[rownames(newdata),colnames(newdata)]
    main_data@location <- main_data@location[rownames(newdata),]
    main_data@gene <- colnames(newdata)
    main_data@cluster <- SEURAT_clus@ident

    return(main_data)
}



