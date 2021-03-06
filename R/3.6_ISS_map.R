######################################################################
##                         RCA mapping                              ##
######################################################################
### Different methods for ploting MolDiaISS
"ISS_map"
#' Plot function for MolDiaISS data.
#'
#' @description Map RCA data based on cell, cluster or tSNE
#' @param data Input data in class MolDiaISS. Output of \link[MolDia]{readISS}.
#' @param what What to plot. Values can be "cell", "cluster", "tsne", "tsneAll" and "vlnplot". Default is "cell". See details.
#' @param xlab Label of x-axis
#' @param ylab Label of y-axis
#' @param ptsize Point size
#' @param pchuse Pch for plotting
#' @param main Main title
#' @param image Plot image. Default is TRUE.
#' @param live Plot interactive image. Default is FALSE.
#' @param label.topgene Active only when "what = cluster or tsne". Number of genes to label each cluster.
#'        Only work when data is clustered and clusted marker has identified in cluster.marker slot of input data.
#' @param gene Gene of interest. When live= TRUE, limited to max 20 genes.
#' @param cluster_id Which cluster to plot. Only work when what = "cluster".
#' @param same.y.lims Set all the y-axis limits to the same values.
#' @param adjust.use A multiplicate bandwidth adjustment. This makes it possible to adjust the 
#'        bandwidth while still using the a bandwidth estimator. For exampe, adjust = 1/2 means use half of the default bandwidth.
#' @param log a character string which contains "x" if the x axis is to be logarithmic, "y" if the y axis is to be logarithmic and
#'        "xy" or "yx" if both axes are to be logarithmic.
#' 
#' @details what parameter can have the value "cell", "gene", "cluster", "tsne", "tsneAll" and "vlnplot".
#'          
#'          "cell" will plot all cells with selected/all genes togather.
#'          
#'          "gene" will plot all/selected gene separately one  by one on the tissue.   
#'          
#'          "cluster" will plot all cluster information.
#'          
#'          "tsne" and "tsneAll" will plot the Dimention reduction by tSNE
#'          
#'          "vlnplot" will plot the violin plot on cluster data.
#'     
#'
#' @examples 
#' ## Reading ISS data
#' left_hypo <- readISS(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"), 
#'                         cellid = "CellId", centX = "centroid_x", centY = "centroid_y")
#'
#' ############################# Plot cell                                                  
#' ## Plot all cell with selected/all genes togather
#' res       <- ISS_map(data = left_hypo, what = "cell", gene = left_hypo@gene[1:4], main = "Plot all cells gene")
#'
#' ############################# Plot gene                                                 
#' ## Plot each single gene
#' res       <- ISS_map(data = left_hypo, what = "gene", gene = left_hypo@gene[1:4], main = "Plot selected genes")
#' 
#' ############################# Data Pre-processing
#' left_hypo <- ISS_preprocess(data = left_hypo, normalization.method = "LogNormalize",
#'                             do.scale = TRUE, do.center = FALSE)
#'                             
#' ############################# Plot tsne and tsneAll
#' ## Dimention reduction by tSNE on non-clustered data
#' left_hypo <- ISS_tsne(data = left_hypo, perplexity= 30, pc = 0.7)
#' # Plot tSNE
#' result <- ISS_map(data = left_hypo, what = "tsneAll")
#' # Plot selected gene on tSNE plot
#' result <- ISS_map(data = left_hypo, what = "tsne", gene = left_hypo@gene[1:4] )
#' 
#' ## Dimention reduction by tSNE on clustered data
#' # Cluster data based on SEURAT pipeline
#' left_hypo_clust  <- ISS_cluster(data = left_hypo, pc = 0.7, resolution = 0.3, method = "seurat")
#' # Dimention reduction by tSNE
#' left_hypo_clust   <- ISS_tsne(data = left_hypo_clust, pc= 0.9, perplexity= 100)
#' Plot cluster on tSNE plot
#' result <- ISS_map(data = left_hypo_clust, what = "tsneAll")
#' 
#' ############################# Plot cluster 
#' # Cluster data based on SEURAT pipeline
#' left_hypo <- ISS_cluster(data = left_hypo, pc = 0.7, resolution = 0.08, method = "seurat")
#' res       <- ISS_map(data = left_hypo, what = "cluster")
#' res       <- ISS_map(data = left_hypo, what = "cluster", cluster_id = 1:2)
#' 
#' ############################# Plot violin plot
#' res       <- ISS_map(data = left_hypo_clust, what = "vlnplot", gene = left_hypo@gene[4:7], same.y.lims = F, adjust.use = 1)
#' 
#' 
#' ###### Reading non-segmentated file
#' data_nonsegment  <- readISS(file = system.file("extdata", "nonSeg_QT_0.35_0_details.csv", package="MolDia"), segment = FALSE,
#'                             centX = "PosX", centY = "PosY", nogene = "NNNN")
#' res   <- ISS_map(data = data_nonsegment, what = "gene")
#' res   <- ISS_map(data = data_nonsegment, what = "gene", gene= c("Gdf7","WNT1","Pak3","Tfap2a"))
#' 
#' @export
ISS_map <- function(data, what = "cell", xlab = "centroid_x", ylab = "centroid_y", main = "Plot title", ptsize = 1,pchuse = 16,
                    image = TRUE, live = FALSE, label.topgene = NULL, gene = NULL, cluster_id = NULL, same.y.lims = FALSE,
                    adjust.use = 0.5, log = "")
{
  ## Check if data is in correct class
  if(class(data)%in%c("MolDiaISS","MolDiaISS_nonsegment") == FALSE ) stop("Please check input data is in class MolDiaISS or MolDiaISS_nonsegment", call. = FALSE)
  
  ## If data is in class MolDiaISS_nonsegment
  if(class(data) == "MolDiaISS_nonsegment")
   {
    
    ## Check what in category
    what_type <- c("gene")
    if(what %in% what_type == FALSE ) stop("Please check available plot type in `what` argument", call. = FALSE)
    
    if(what == "gene")
    {
      ## Data with all gene
      if(length(gene) == 0 )
      {
       genes <- data@gene
       data  <- data@location
      }
      
      ## Data with selected gene
      if(length(gene) > 0 )
      {
        if(all(gene%in%data@gene)==FALSE) stop("Selected gene not present. Check gene name in 'gene'.", call. = FALSE)
        genes <- data@gene[which(data@gene%in%gene)]
        data  <- data@location[names(genes),] # my_file[which(my_file$Gene%in%gene),]
      }
      
      ## Check image and live has same value
      if (image == live ) stop("Argument value of `live` or `image should be different", call. = FALSE)
      
      ## Create ggplot object
      p <- ggplot2::ggplot(data, ggplot2::aes_string(x= "centroid_x", y= "centroid_y")) +
           ggplot2::geom_point(ggplot2::aes(colour= genes)) +
           ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)), fill=ggplot2::guide_legend(nrow = 10)) +
           ggplot2::labs(x = "", y = "", title= main) +
           { if(log == "y")  ggplot2::scale_y_continuous(limits = c(min(data$centroid_y),max(data$centroid_y)), trans = "log10") } +
           { if(log == "x")  ggplot2::scale_x_continuous(limits = c(min(data$centroid_x),max(data$centroid_x)), trans = "log10") } +
           { if(log == "xy") ggplot2::scale_y_continuous(limits = c(min(data$centroid_y),max(data$centroid_y)), trans = "log10") } +
           { if(log == "xy") ggplot2::scale_x_continuous(limits = c(min(data$centroid_x),max(data$centroid_x)), trans = "log10") } +
           { if(log == "" )  ggplot2::scale_x_continuous(limits = c(min(data$centroid_x),max(data$centroid_x))) } + 
           { if(log == "" )  ggplot2::scale_y_continuous(limits = c(min(data$centroid_y),max(data$centroid_y))) } +
           #ggplot2::scale_x_continuous(limits = c(min(data$centroid_x),max(data$centroid_x)))+ #,expand=c(0,0)) +
           #ggplot2::scale_y_continuous(limits = c(min(data$centroid_y),max(data$centroid_y)))+ #,expand=c(0,0)) +
           ggplot2::theme_void() + 
           ggplot2::theme(legend.title=ggplot2::element_blank(), legend.position="right")
       
      ## Select which plot to show
      if(live)  
      {
        if(length(unique(genes)) > 20) stop("Gene number more than 20 is too slow for interactive vizualization. So inactive at this moment", call. = FALSE)
        q <- plotly::ggplotly(p)
      }
      
      if(image) 
      {
        q <- p
      }
      #print(p)
    }
    
    ## Print Image
    print(q)
   }
  
  ## If data is in class MolDiaISS
  if(class(data) == "MolDiaISS")
  {
  
  ## Save main data
  mainData <- data
  
  ## Check gene
  if(length(gene) > 20) stop("Gene number more than 20 is too slow for interactive vizualization. So inactive at this moment", call. = FALSE)
  
  ## If gene of interest is present
  if(length(gene)!=0)
  {
    ## Gene name
    gene <- gene
    ## Check if all Gene is present in data
    if(all(gene %in% colnames(data@data)) == FALSE ) stop("Gene of interest not present in data", call. = FALSE)
    data@data <- data@data[,gene,drop = FALSE]
    if(length(data@norm.data) !=0) data@norm.data  <- data@norm.data[,gene,drop = FALSE]
    if(length(data@scale.data)!=0) data@scale.data <- data@scale.data[,gene,drop = FALSE]
    data@gene <- gene
  }
  
  ## Gene name
  if(length(gene)==0) gene <- data@gene
  
  ## Check what in category
  what_type <- c("cell","gene","cluster","tsne", "tsneAll","vlnplot")
  if(what %in% what_type == FALSE ) stop("Please check available plot type in `what` argument", call. = FALSE)
  
  ## Check image and live has same value
  if (image == live ) stop("Argument value of `live` or `image should be different", call. = FALSE)
  
  ######## What = "cell"
  if(what == "cell")
  {
    ## If data is clustered or not
    if(length(data@cluster) > 0)
    {
      ## Select cluster to plot
      if(length(cluster_id) > 0 )
      { data@cluster <-  droplevels(data@cluster[which(data@cluster %in% cluster_id )]) } 
      
      ## Select cluster information data
      data1 <- data.frame(data@cluster)
      colnames(data1) <- "cluster"
    }
    
    if(length(data@cluster) == 0)
    {
      data1 <- data.frame(rep("",nrow(data@data)))
      colnames(data1) <- "cluster"
      rownames(data1) <- rownames(data@data)
    } 
    
    ## Merge expression and location information
    data<-merge(data@data, data@location, by = "row.names")
    rownames(data) <- data$Row.names
    data$Row.names <- NULL
    data <- reshape::melt (data, id = c("centroid_x","centroid_y"))
    data$value1 <- data$value>0
    data       <- data[data$value1,]
    data$value1 <- NULL
    colnames(data) <- c("centroid_x", "centroid_y", "Genes", "Expression")
    
    ## Create gggplot object
    p <- ggplot2::ggplot(data, ggplot2::aes_string(x= "centroid_x", y= "centroid_y")) +
         ggplot2::geom_point(ggplot2::aes(colour= Genes, size = Expression), alpha=0.75) +
         ggplot2::theme(legend.title=ggplot2::element_blank(), legend.position="right") +
         ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)),
                         fill=ggplot2::guide_legend(nrow = 10)) +
         ggplot2::labs(x = "", y = "", title= main) +
         ggplot2::scale_x_continuous(limits = c(min(data$centroid_x),max(data$centroid_x)),expand=c(0,0)) +
         ggplot2::scale_y_continuous(limits = c(min(data$centroid_y),max(data$centroid_y)),expand=c(0,0)) +
         ggplot2::theme_void()
    
    ## Select which plot to show
    if(live)  q <- plotly::ggplotly(p)
    if(image) q <- p
  }
  
  
  ######## What = "gene"
  if(what == "gene")
  {
    if(length(data@location) == 0 ) stop("Single cell no spatial information", call. = FALSE)
    
    ## Create SEURAT object 
    RCAtsne   <- Seurat::CreateSeuratObject(raw.data = t(mainData@data))
    #if(length(data@scale.data)== 0 ) stop("Please scale data first with ISS_preprocess function ", call. = FALSE)
    
    cell.embeddings<- as.matrix(data@location) 
    reduction.key  <- "tSNE_"
    pca.obj <- new(
      Class = "dim.reduction",
      #gene.loadings = gene.loadings,
      cell.embeddings = cell.embeddings,
      #sdev = sdev,
      key = reduction.key
    )
    RCAtsne@dr$tsne<- pca.obj
    
    # Select color
    colplt<- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")
    
    # Selct feature for tsne plot for all gene
    if(length(gene) == 0 ) {fplot <- data@gene
    }else { fplot <- gene }
    
    # Plot feature
    mm <- Seurat::FeaturePlot(object = RCAtsne, features.plot = c(fplot), cols.use = colplt, pt.size = ptsize, pch.use = pchuse, 
                              do.return = F, no.axes = TRUE)
    return(mm)
  }
  
  ######## What = "cluster"
  if(what == "cluster")
  {
    ## If data is already clustered and label cluster by gene
    if(length(label.topgene)>0 )
    {
      if(length(data@cluster.marker) == 0 ) stop("Please define cluster marker on data first", call. = FALSE)
      new.cluster.id <- data@cluster.marker
      new.cluster.id <- lapply(seq_along(new.cluster.id), function(i)
      {
        myres <- stats::na.omit(as.vector(new.cluster.id[[i]][,"gene"][1:label.topgene]))
        myres <- paste0(c(names(new.cluster.id[i]),myres), collapse = "_")
        names(myres)<- names(new.cluster.id)[i]
        myres
      })
      new.cluster.id<- unlist(new.cluster.id)
      
      ## Assign new cluster name
      new.cluster <- data@cluster
      levels(new.cluster) <- new.cluster.id
      
      ## Replace new cluster in main data
      data@cluster <- new.cluster
    }
    
    ## If data is clustered or not
    if(length(data@cluster) > 0)
    {
      ## Select cluster to plot
      if(length(cluster_id) > 0 )
      { data@cluster <-  droplevels(data@cluster[which(data@cluster %in% cluster_id )]) } 
      
      ## Select cluster information data
      data1 <- data.frame(data@cluster)
      colnames(data1) <- "cluster"
    }
    
    if(length(data@cluster) == 0)
    {
      data1 <- data.frame(rep("",nrow(data@data)))
      colnames(data1) <- "cluster"
      rownames(data1) <- rownames(data@data)
    }
    
    ## Merge gene expression data with cluster data
    data1 <- merge(data@data, data1, by = "row.names")
    rownames(data1) <- data1$Row.names
    data1$Row.names <- NULL
    
    ## Merge location data with cluster data
    data <- merge(data@location,data1, by = "row.names")
    rownames(data) <- data$Row.names
    data$Row.names <- NULL
    
    ## Melt data
    data <- reshape::melt(data, c("centroid_x","centroid_y","cluster"))
    
    ## Create gggplot object
    if(length(gene)<= 6 )
    {
      p<- ggplot2::ggplot(data,ggplot2::aes_string(x= "centroid_x",y= "centroid_y")) +
        ggplot2::geom_point(ggplot2::aes(colour=factor(cluster), shape = variable, size = ptsize, shape = pchuse,
                                         alpha=log2(value+1))) +
        ggplot2::scale_alpha(range = c(0, 1)) +
        ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position="right") +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)),
                        fill=ggplot2::guide_legend(nrow = 10)) +
        ggplot2::labs(x = "", y = "", colour = "Cluster", shape = "Gene",alpha = "log(Reads)",title= main)+
        ggplot2::scale_x_continuous(limits = c(min(data$centroid_x),max(data$centroid_x)),expand=c(0,0)) +
        ggplot2::scale_y_continuous(limits = c(min(data$centroid_y),max(data$centroid_y)),expand=c(0,0)) +
        ggplot2::theme_void()
      
      ## Select which plot to show
      if(live)  q <- plotly::ggplotly(p)
      if(image) q <- p
      
      ## Select which plot to show
      if(live)
      {
        if(length(gene) > 6) stop("Gene number more than 6 is too slow for interactive vizualization. So inactive at this moment", call. = FALSE)
        q <- plotly::ggplotly(p)
      }
      if(image) q <- p
    }
    
    if(length(gene) > 6 )
    {
      p<- ggplot2::ggplot(data,ggplot2::aes_string(x= "centroid_x",y= "centroid_y")) +
        ggplot2::geom_point(ggplot2::aes(colour=factor(cluster)),size = ptsize, shape = pchuse) +
        ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position="right") +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)),
                        fill=ggplot2::guide_legend(nrow = 10)) +
        ggplot2::labs(x = "", y = "", colour = "Cluster", title= main) +
        ggplot2::scale_x_continuous(limits = c(min(data$centroid_x),max(data$centroid_x)),expand=c(0,0)) +
        ggplot2::scale_y_continuous(limits = c(min(data$centroid_y),max(data$centroid_y)),expand=c(0,0)) +
        ggplot2::theme_void()
      
      ## Select which plot to show
      if(live)  q <- plotly::ggplotly(p)
      if(image) q <- p
    }
  }
  
  ######## What = "tsne"
  if(what == "tsne")
  {
    # Check if tsne applied
    if(length(data@tsne.data) == 0 ) stop("Please perform tSNE on data first", call. = FALSE)
    
    ## Create SEURAT object 
    RCAtsne   <- Seurat::CreateSeuratObject(raw.data = t(mainData@data))
    if(length(data@scale.data)== 0 ) stop("Please scale data first with ISS_preprocess function ", call. = FALSE)

    cell.embeddings<- as.matrix(data@tsne.data) 
    #cell.embeddings<- as.matrix(left_hypo@location) 
    reduction.key  <- "tSNE_"
    pca.obj <- new(
      Class = "dim.reduction",
      #gene.loadings = gene.loadings,
      cell.embeddings = cell.embeddings,
      #sdev = sdev,
      key = reduction.key
    )
    RCAtsne@dr$tsne<- pca.obj
    
    # Select color 
    colplt<- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")
    
    # Selct feature for tsne plot for all gene
    if(length(gene) == 0 ) {fplot <- data@gene
    }else { fplot <- gene }
    
    # Plot feature
    mm<- Seurat::FeaturePlot(object = RCAtsne, features.plot = fplot, cols.use = colplt, pt.size = ptsize, 
                             pch.use = pchuse, no.axes = TRUE)
    return(NULL)
  }
  
  ######## What = "tsneAll"
  if(what == "tsneAll")
  {
    if(length(data@tsne.data) == 0 ) stop("Please perform tSNE on data first", call. = FALSE)
    
    ## If data is already clustered and label cluster by gene
    if(length(label.topgene)>0 )
    {
      if(length(data@cluster.marker) == 0 ) stop("Please define cluster marker on data first", call. = FALSE)
      new.cluster.id <- data@cluster.marker
      new.cluster.id <- lapply(seq_along(new.cluster.id), function(i)
      {
        myres <- stats::na.omit(as.vector(new.cluster.id[[i]][,"gene"][1:label.topgene]))
        myres <- paste0(c(names(new.cluster.id[i]),myres), collapse = "_")
        names(myres)<- names(new.cluster.id)[i]
        myres
      })
      new.cluster.id<- unlist(new.cluster.id)
      
      ## Assign new cluster name
      new.cluster <- data@cluster
      levels(new.cluster) <- new.cluster.id
      
      # Replace new cluster in main data
      data@cluster <- new.cluster
    }
    
    ## If data is clustered or not
    if(length(data@cluster) > 0)
    {
      
      ## Select cluster to plot
      if(length(cluster_id) > 0 )
      { data@cluster <-  droplevels(data@cluster[which(data@cluster %in% cluster_id )]) } 
      
      ## Select cluster information data
      data1 <- data.frame(data@cluster)
      colnames(data1) <- "cluster"
    }
    if(length(data@cluster) == 0)
    {
      data1 <- data.frame(rep(0,nrow(data@data)))
      colnames(data1) <- "cluster"
      rownames(data1) <- rownames(data@data)
    }
    
    ## Merge gene expression data with cluster data
    data1 <- merge(data@data, data1, by = "row.names")
    rownames(data1) <- data1$Row.names
    data1$Row.names <- NULL
    
    ## Merge tsne data with cluster data
    data <- merge(data@tsne.data,data1, by = "row.names")
    rownames(data) <- data$Row.names
    data$Row.names <- NULL
    
    ## Melt data
    data <- reshape::melt(data, c("tSNE_1","tSNE_2","cluster"))
    
    ## Create gggplot object
    if(length(gene)<= 6 )
    {
      p<- ggplot2::ggplot(data,ggplot2::aes_string(x= "tSNE_1",y= "tSNE_2")) +
        ggplot2::geom_point(ggplot2::aes(colour=factor(cluster), size = ptsize,
                                         shape = variable, alpha=log2(value+1))) +
        ggplot2::scale_alpha(range = c(0, 1)) +
        ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position="right") +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)),
                        fill=ggplot2::guide_legend(nrow = 10)) +
        ggplot2::labs(x = "", y = "", colour = "Cluster", shape = "Gene",alpha = "Reads", title= main)+
        ggplot2::scale_x_continuous(limits = c(min(data$tSNE_1),max(data$tSNE_1)),expand=c(0,0)) +
        ggplot2::scale_y_continuous(limits = c(min(data$tSNE_2),max(data$tSNE_2)),expand=c(0,0)) +
        ggplot2::theme_void()
      
      ## Select which plot to show
      if(live)  q <- plotly::ggplotly(p)
      if(image) q <- p
      
      ## Select which plot to show
      if(live)
      { if(length(gene) > 6) stop("Gene number more than 6 is too slow for interactive vizualization. So inactive at this moment", call. = FALSE)
        q <- plotly::ggplotly(p)
      }
      if(image) q <- p
    }
    
    if(length(gene) > 6 )
    {
      p<- ggplot2::ggplot(data,ggplot2::aes_string(x= "tSNE_1",y= "tSNE_2")) +
        ggplot2::geom_point(ggplot2::aes(colour=factor(cluster)), size = ptsize, shape = pchuse) +
        ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position="right") +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)),
                        fill=ggplot2::guide_legend(nrow = 10)) +
        ggplot2::labs(x = "", y = "", colour = "Cluster", title= main) +
        ggplot2::scale_x_continuous(limits = c(min(data$tSNE_1),max(data$tSNE_1)),expand=c(0,0)) +
        ggplot2::scale_y_continuous(limits = c(min(data$tSNE_2),max(data$tSNE_2)),expand=c(0,0)) +
        ggplot2::theme_void()
      
      ## Select which plot to show
      if(live)  q <- plotly::ggplotly(p)
      if(image) q <- p
    }
  }
  ######## What = "vlnplot"
  if(what == "vlnplot")
  {
    # Main data
    data  <- mainData
    
    # Check if data is clustered or not
    if (length(data@cluster) == 0) stop("Please cluster your data befor finding marker gene of cluster", call. = FALSE) 
    
    # Create SEURAT object
    RCA_Seurat_obj            <- Seurat::CreateSeuratObject(t(data@data))
    RCA_Seurat_obj@scale.data <- t(data@scale.data)
    RCA_Seurat_obj@ident      <- data@cluster
    
    # Selct feature for violin plot
    if(length(gene) == 0 ) {fplot <- data@gene
    }else { fplot <- gene }
    
    # same y limit
    if(same.y.lims == TRUE)  ymax <- max(data@data[,fplot])
    if(same.y.lims == FALSE) ymax <- NULL
    
    # Violin plot
    for(i in 1: length(fplot))
    {
      v <- Seurat::VlnPlot(object = RCA_Seurat_obj,features.plot = fplot[i], same.y.lims = TRUE, do.return = TRUE, y.max = ymax,adjust.use = adjust.use) + 
           ggplot2::labs(x = "Cluster", y = "Reads count")
      assign(x = paste0("p",i), value = v)
    }
    
    
    tt<- paste0("cowplot::plot_grid(",paste0("p",1:length(fplot),collapse = " ,"),")")
    q <-  suppressWarnings(eval(parse(text = tt)))
    
  }
  
  ## Print Image
  print(q)
  #return(data)
  }
}