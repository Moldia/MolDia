######################################################################
##                      RCA cluster validation                      ##
######################################################################
"RCA_tsne"
#'
#'Dimentionality reduction to 2D by tSNE
#'
#' @description Any data in class MolDiaISS clusteded or not clustered used to reduce dimention to 2D by RCA-tsne.
#' @param data Input data in class MolDiaISS. Output of \link[MolDia]{readISS}.
#' @param clus As factor. Each cell names has id. Cell order of data and cell order of clus should be same. (work for class data.frame)
#' @param pc Desired percent of variance to be explained by PCA. Default in NULL , equavalant with 0.9.
#' @param perplexity Numeric; Perplexity parameter. See  \link[Rtsne]{Rtsne}
#' 
#' ## param do.label Label tsne plot or not
#'
#' @return 2D dataframe of points.
#'
#' @examples
#' ## Reading data
#' data_3 <- readISS(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),
#'                   cellid = "CellId", centX = "centroid_x", centY = "centroid_y")
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
#' neuron_group_clust  <- RCA_cluster(data = neuron_group, pc = 0.9, resolution = 0.3)
#'
#' ## Dimention reduction by tSNE on non-clustered data
#' tsne_noclust <- RCA_tsne(data = neuron_group)
#'
#' ## Dimention reduction by tSNE on clustered data
#' tsne_clust   <- RCA_tsne(data = neuron_group_clust, pc= 0.9, perplexity= 100)
#'
#' ## Get cluster marker
#' tsne_clust <- RCA_marker(data = tsne_clust, topgene =15, marker.sig = 0.005,
#'                                         test.use="bimod", main = "")
#' ## Plot tSNE data
#' result <- RCA_map(data = tsne_clust, what = "tsne")
#'
#' @export
RCA_tsne <- function(data, clus = NULL, pc = NULL, perplexity = 100)
{
  # Check class of data
  if(class(data)%in%c("MolDiaISS") ==FALSE) stop("Check input data in class 'MolDiaISS'",call. = FALSE)

  # Save main data 
  mdata <- data
  
  # Run tSNE on normalize / non normalized data
  data  <- mdata@data

  # Create SEURAT object
  RCAtsne   <- Seurat::CreateSeuratObject(t(data))
  #RCAtsne   <- Seurat::ScaleData(object = RCAtsne, do.scale = TRUE, do.center = TRUE, check.for.norm = FALSE)

  if(length(mdata@scale.data)== 0 ) stop("Please scale data first with RCA_preprocess function ", call. = FALSE)
  #if(length(mdata@scale.data)> 0 ) stop("Please scale data first", call. = FALSE)
  RCAtsne@scale.data <- t(mdata@scale.data[rownames(data),])
  RCAtsne   <- Seurat::RunPCA(object  = RCAtsne, pc.genes = colnames(data),do.print = FALSE)
  
  # Find optimal PCA component
  if(length(pc) > 0 )
  {
    # Find optimal PCA component
    #myPCA  <- SEURAT_clus@dr$pca@cell.embeddings
    sdev   <- RCAtsne@dr$pca@sdev
    pcuse  <- cumsum(log2(sdev)^2 / sum(log2(sdev)^2))
    pcuse  <- max(which(pcuse<=pc))
    #if(pcuse <= 3) stop("PC is too low", call. = FALSE)
    #cat("Number of principle component to be used :", pcuse, "\n")
    cat("Principle component: ", round(npc[pcuse]*100,2), "% of variation has explained by" ,pcuse, "principle component", "\n")
    pc <- 1:pcuse
  }
  if(length(pc) == 0 )
  {
    ## Find number of optimal principle component that explain 90 percent of variaiance
    if(length(colnames(data)) < 20){ ## Number of PC should be less than number of source of variation i.e. genes.
      npc <- ncol(data) -1
      npc   <- withCallingHandlers(suppressWarnings(irlba::prcomp_irlba(RCAtsne@data, n=npc, 
                                                                        fastpath = TRUE, verbose = FALSE)))}
    else{
      npc   <- withCallingHandlers(suppressWarnings(irlba::prcomp_irlba(RCAtsne@data, n=20, 
                                                                        fastpath = TRUE, verbose = FALSE)))}

    # How many pc to be used: Cut off is 0.90
    npc   <- summary(npc)$importance[3,]
    if(max(npc) < 0.9) {pcuse <- length(npc)
    }else {pcuse <- which(npc > 0.90)[1]}
    #cat("Number of principle component to be used :", pcuse, "\n")
    cat("Principle component: ", round(npc[pcuse]*100,2), "% of variation has explained by" ,pcuse, "principle component", "\n")
    pc <- 1:pcuse
  }
  
  # Stopping criteria for principle component
  if(pcuse <= 3) stop("Number of principle component should be greater than 3. Please increase the value of pc", call. = FALSE)
  
  # Find optimal PCA component
  #sdev   <- RCAtsne@dr$pca@sdev
  #pcuse  <- cumsum(log2(sdev)^2 / sum(log2(sdev)^2))
  #pcuse  <- max(which(pcuse<=pc))
  #if(pcuse <= 3) stop("PC is too low")
  #cat("Number of principle component to be used :", pcuse, "\n")

  # Running tSNE
  set.seed(12345)
  cat("Running dimentionality reduction by tSNE..")
  RCAtsne   <- Seurat::RunTSNE(object = RCAtsne , dims.use = pc, do.fast = TRUE, check_duplicates = FALSE,perplexity = perplexity )

  # Assign cluster information on tSNE plot
  #if(length(mdata@cluster)> 0 ) RCAtsne@ident <- mdata@cluster
  #RCAtsne_1 <- Seurat::TSNEPlot(RCAtsne, do.label = do.label, label.size = 8)
  RCAtsne   <- RCAtsne@dr$tsne@cell.embeddings

  # Assign tSNE data to slot
  mdata@tsne.data <- data.frame(RCAtsne)
  return(mdata)
}
