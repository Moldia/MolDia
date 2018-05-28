######################################################################
##                        Select cluster of interest                ##
######################################################################
"RCA_ClusterSelect"
#' Select cluster of interest after clustering.
#' 
#' @description Select cluster of interest after clustering
#' @param data Input data in class RCA_class. Output of \link[MolDia]{readRCA}.
#' @param cluster_id Cluster to select
#' 
#' @export
RCA_ClusterSelect <- function(data, cluster_id = NULL)
{
  ## Check id data is cluster or not 
  if(length(data@cluster) == 0) stop("Please define cluster on data first", call. = FALSE)
  
  ## Select cluster to select
  if(length(cluster_id) == 0) stop("Please define cluster to select",call. = FALSE)
  
  ## Selected cluster and cellid
  data@cluster <-  droplevels(data@cluster[which(data@cluster %in% cluster_id )]) 
  cellid <- names(data@cluster)
  
  ## Select all other dada based on selected cell 
  data@data       <- data@data[cellid,]
  if(length(data@norm.data)  > 0) data@norm.data  <- data@norm.data[cellid,]
  if(length(data@scale.data) > 0) data@scale.data <- data@scale.data[cellid,]
  data@gene       <- data@gene
  data@location   <- data@location[cellid,]
  if(length(data@cluster.marker) > 0) 
   {
    data@cluster.marker <- data@cluster.marker[match(as.character(cluster_id),names(data@cluster.marker))]
   }

  if(length(data@cluster.marker) > 0) 
  { 
    data@cluster.marker <- data@cluster.marker[match(as.character(cluster_id),names(data@cluster.marker))]
  }
  if(length(data@tsne.data) > 0) data@tsne.data  <- data@tsne.data[cellid,]
  
  ## Return data 
  return(data)
}