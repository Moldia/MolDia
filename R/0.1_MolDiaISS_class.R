#' The ISS data class
#' The MolDiaISS class is a main object for all ISS based analysis
#'
#' A MolDiaISS class object has the following slot
#'
#' @slot data Main data for analysis
#' @slot norm.data Normalized data (log-scale).
#' @slot scale.data Scaled (default is z-scoring each gene) expression matrix; used for dimmensional reduction and heatmap visualization.
#' @slot location Insitu location of RCA cell location on X and Y co-ordinate.
#' @slot gene Name of genes in the data.
#' @slot cluster Cluster identification of all cells
#' @slot cluster.marker Marker gene of each cluster.
#' @slot tsne.data tSNE co-ordinate for data.
#' @slot pca.data List of stored dimmensional reductions
#'
#' @section Dataclass: Data class section
#'
#' @name MolDiaISS
MolDiaISS <- methods::setClass(Class = "MolDiaISS",
                               slots = c(data           = "data.frame",
                                         norm.data      = "matrix",
                                         scale.data     = "matrix",
                                         gene           = "vector",
                                         cluster        = "factor",
                                         location       = "data.frame",
                                         cluster.marker = "list",
                                         tsne.data      = "data.frame",
                                         pca.data       = "list"))


methods::setMethod( f = "show", signature = "MolDiaISS",
                    definition = function(object)
                      {
                      cat( "An object of class", class(object),
                           "\nNumber of genes:", ncol(object@data),#sum(colSums(object@data)>=0),
                           "\nNumber of cells:", nrow(object@data),#sum(rowSums(object@data)>=0),
                           "\nNumber of reads:", sum(object@data),
                           if(length(object@cluster) > 0) "\nCluster id:", levels(object@cluster),
                           if(length(object@cluster) > 0) "\nCluster size:",table(object@cluster))
                      invisible(x = NULL)
                      }
                    )
MolDiaISS <- methods::new("MolDiaISS")
