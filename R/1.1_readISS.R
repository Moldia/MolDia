######################################################################
##                           Read RCA data                          ##
######################################################################
"readISS"
#' Read RCA data
#' @description Read RCA data
#'
#' @param file File name in CSV formate, data in "data.frame" class and data in "MolDiaISS" class (Output of \link[MolDia]{readISS}). 
#' @param cellid String to naming cell. Default is "CellID".
#' @param centX Name of X co-ordinate in file. Default is "centroidX"
#' @param centY Name of Y co-ordinate in file  Default is "centroidY"
#' @param genepos Name of genes to consider for gene positive cells. Default is NULL
#' @param geneposOPT Only work when 'genepos' has a value. "AND", "OR" and "NONE" condition for genepos. Default is "OR".
#' @param rpc Total reads per cell to be consider. Default is 1.
#' @param rpg Total reads per gene to be consider. Default is 1.
#' @param gene Gene name to include in data. Default is NULL.
#' @param nogene Gene name to exclude from data. Default is NULL.
#'
#' @author Mohammad Tanvir Ahamed
#'
#' @return Output will be a object in class MolDiaISS. See detail \link[MolDia]{MolDiaISS}
#'
#' @examples
#' ##### Reading ISS data in CSV formate
#' data_1 <- readISS(file = system.file("extdata", "CellBlobs_QT_0.35.csv", package="MolDia"),
#'                   cellid = "CellID", centX = "centroidX", centY = "centroidY")
#' data_2 <- readISS(file = system.file("extdata", "CellBlobs_QT_0.40.csv", package="MolDia"),
#'                   cellid = "CellID",centX = "centroidX", centY = "centroidY")
#' data_3 <- readISS(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),
#'                   cellid = "CellId",centX = "centroid_x", centY = "centroid_y")
#' data_4 <- readISS(file = system.file("extdata", "Hypocampus_right.csv", package="MolDia"),
#'                   cellid = "CellId",centX = "centroid_x", centY = "centroid_y")
#'
#' ## Not RUN
#' ## Read RCA data in dataframe formate
#' #gene <- data_3@gene[1:3]
#' #data(single_cell)
#' #single_cell$CellID <- rownames(single_cell)
#' #Note: Big data. Take long time to load without selected gene
#' #data_sc <- readISS(file = single_cell, cellid = "CellID",gene= gene) 
#'
#'
#' @export
readISS <- function(file, cellid = "CellID", centX = NULL, centY = NULL, genepos= NULL, geneposOPT = "OR",
                    rpc = 1, rpg = 1, gene = NULL, nogene = NULL)
{
  ## Data type : REading data from a specific location
  if(class(file)=="character")
  {
    ## Reading data
    my_file <- data.frame(data.table::fread(input = file , showProgress = TRUE))
    
    ## Delete cell with NA values
    indexNA <- apply(my_file,2,function(i)
    {
      kk<- which(is.na(i))
      kk
    })
    indexNA <- unique(unlist(indexNA))
    
    if(length(indexNA) == 0 ) my_file <- my_file
    if(length(indexNA) > 0 ) 
    {
      cat("Removed" , length(indexNA), "cells having NA\n")
      my_file <- my_file[-indexNA,]
    }
    
    
    ## Check cellid, centx, centy all in data or not
    if(any(c(cellid,centX,centY)%in%colnames(my_file))==FALSE) stop("Pleace check `cellid`, `centX` and `centY`", call. = TRUE)
    
    
    ## Define cell id
    cell_id <- paste0("cellid_",as.character(my_file[,cellid] )); my_file[,cellid] <- NULL;
    row.names(my_file) <- cell_id
    
    ## Data of Reads count
    data_reads <- my_file[,!(names(my_file) %in% c(centX,centY))]
    
    ## Data of location information
    data_loca <- my_file[, (names(my_file) %in% c(centX,centY))]
    data_loca <- data.frame(data_loca)
    colnames(data_loca) <- c("centroid_x", "centroid_y")
    
    ## Selected / un-select genes
    if(length(gene) > 0 )   data_reads <- data_reads[,gene, drop=FALSE]
    if(length(nogene) > 0 ) data_reads <- data_reads[ , -which(colnames(data_reads) %in% nogene), drop=FALSE]
    
    ## Filter by reads per cell (rpc) and reads per gene (rpg)
    data_reads <- data_reads[rowSums(data_reads) >= rpc,, drop=FALSE]
    data_reads <- data_reads[,colSums(data_reads) >= rpg, drop=FALSE]
    
    ## Filter by genes+ cells
    if(length(genepos) > 0 )
    { 
      if(any(genepos %in% colnames(my_file))==FALSE) stop("Pleace check `genepos` ", call. = TRUE)
      if (geneposOPT == "AND")  genepos <- paste(genepos, " > 0" , collapse = " & ") 
      if (geneposOPT == "OR")   genepos <- paste(genepos, " > 0" , collapse = " | ")
      if (geneposOPT == "NONE") genepos <- paste(genepos, " == 0", collapse = " & ")
      data_reads  <- subset(data_reads, eval(parse(text = genepos)))
    } 
    
    ## Deleter emplt genes and emplt cells 
    data_reads <- data_reads[,colSums(data_reads) > 0]
    data_reads <- data_reads[rowSums(data_reads) > 0,]
    
    ## Filtered cell location
    data_loca<- data_loca[rownames(data_reads),]
    
    ## return RCA object
    res <- methods::new("MolDiaISS",
                        data     = data_reads,
                        norm.data = matrix(, nrow = 0, ncol = 0),
                        scale.data = matrix(, nrow = 0, ncol = 0),
                        gene = colnames(data_reads),
                        cluster = factor(matrix(, nrow = 0, ncol = 0)),
                        location = data_loca,
                        cluster.marker = list(),
                        tsne.data = data.frame(matrix(, nrow = 0, ncol = 0)))
  }
  
  ## Data type : Dataframe
  if(class(file)=="data.frame" )
  {
    my_file <- file
    
    ## Selected gene
    #if(length(gene) >0 )
    #  {
    #  gene<- c(gene[gene%in%colnames(my_file)],cellid)
    #  my_file <- my_file[,gene,drop = FALSE]
    #}
    
    ## Define cell ID
    cell_id <- paste0("cellid_",as.vector(my_file[,cellid]))
    rownames(my_file) <- cell_id
    my_file[,cellid] <- NULL
    my_file <- apply(my_file,2,as.numeric)
    my_file <- data.frame(as.matrix(my_file))
    
    
    if(length(centX)== 1 )
    { 
      ## Data of Reads count
      data_reads <- my_file[,!(names(my_file) %in% c(centX,centY))]
      
      ## Data of location information
      data_loca <- my_file[, (names(my_file) %in% c(centX,centY))]
      data_loca <- data.frame(data_loca)
      colnames(data_loca) <- c("centroid_x", "centroid_y")
    }
    if(length(centX)== 0 )
    {
      data_reads <- my_file
      data_loca  <- data.frame(matrix(, nrow = 0, ncol = 0))
    }
    
    ## Selected / un-select genes
    if(length(gene) > 0 )   data_reads <- data_reads[ ,  gene , drop=FALSE]
    if(length(nogene) > 0 ) data_reads <- data_reads[ , -which(colnames(data_reads) %in% nogene), drop=FALSE]
    
    ## Filter by reads per cell (rpc) and reads per gene (rpg)
    data_reads <- data_reads[rowSums(data_reads) >= rpc, , drop = FALSE]
    data_reads <- data_reads[,colSums(data_reads) >= rpg,drop = FALSE]
    
    ## Filter by genes+ cells
    if(length(genepos) > 0 )
    { 
      if(any(genepos %in% colnames(my_file))==FALSE) stop("Pleace check `genepos` ", call. = TRUE)
      if (geneposOPT == "AND")  genepos <- paste(genepos, " > 0" , collapse = " & ") 
      if (geneposOPT == "OR")   genepos <- paste(genepos, " > 0" , collapse = " | ")
      if (geneposOPT == "NONE") genepos <- paste(genepos, " == 0", collapse = " & ")
      data_reads  <- subset(data_reads, eval(parse(text = genepos)))
    }
    
    ## Deleter emplt genes and emplt cells 
    data_reads <- data_reads[,colSums(data_reads) > 0]
    data_reads <- data_reads[rowSums(data_reads) > 0,]
    
    ## Filtered cell location
    if(length(centX)== 1 ) data_loca<- data_loca[rownames(data_reads),]
    
    ## return RCA object
    res <- methods::new("MolDiaISS",
                        data     = data_reads,
                        norm.data = matrix(, nrow = 0, ncol = 0),
                        scale.data = matrix(, nrow = 0, ncol = 0),
                        gene = colnames(data_reads),
                        cluster = factor(matrix(, nrow = 0, ncol = 0)),
                        location = data_loca,
                        cluster.marker = list(),
                        tsne.data = data.frame(matrix(, nrow = 0, ncol = 0)))
  }
  
  
  ## Data type : MolDiaISS class 
  if(class(file)=="MolDiaISS" )
  {
    my_file <- file@data
    
    ## Selected gene
    #if(length(gene) >0 )
    #  {
    #  gene<- c(gene[gene%in%colnames(my_file)],cellid)
    #  my_file <- my_file[,gene,drop = FALSE]
    #}
    
    ## Define cell ID
    #cell_id <- paste0("cellid_",as.vector(my_file[,cellid]))
    #rownames(my_file) <- cell_id
    #my_file[,cellid] <- NULL
    #my_file <- apply(my_file,2,as.numeric)
    #my_file <- data.frame(as.matrix(my_file))
    
    
    if(length(centX)== 1 )
    { 
      ## Data of Reads count
      data_reads <- my_file[,!(names(my_file) %in% c(centX,centY))]
      
      ## Data of location information
      data_loca <- my_file[, (names(my_file) %in% c(centX,centY))]
      data_loca <- data.frame(data_loca)
      colnames(data_loca) <- c("centroid_x", "centroid_y")
    }
    if(length(centX)== 0 )
    {
      data_reads <- my_file
      data_loca  <- file@location #data.frame(matrix(, nrow = 0, ncol = 0))
    }
    
    ## Selected / un-select genes
    if(length(gene) > 0 )   data_reads <- data_reads[ ,  gene , drop=FALSE]
    if(length(nogene) > 0 ) data_reads <- data_reads[ , -which(colnames(data_reads) %in% nogene), drop=FALSE]
    
    ## Filter by reads per cell (rpc) and reads per gene (rpg)
    data_reads <- data_reads[rowSums(data_reads) >= rpc, , drop = FALSE]
    data_reads <- data_reads[,colSums(data_reads) >= rpg,drop = FALSE]
    
    ## Filter by genes+ cells
    if(length(genepos) > 0 )
    { 
      if(any(genepos %in% colnames(my_file))==FALSE) stop("Pleace check `genepos` ", call. = TRUE)
      if (geneposOPT == "AND")  genepos <- paste(genepos, " > 0" , collapse = " & ") 
      if (geneposOPT == "OR")   genepos <- paste(genepos, " > 0" , collapse = " | ")
      if (geneposOPT == "NONE") genepos <- paste(genepos, " == 0", collapse = " & ")
      data_reads  <- subset(data_reads, eval(parse(text = genepos)))
    }
    
    ## Delete empty genes and empty cells 
    data_reads <- data_reads[,colSums(data_reads) > 0]
    data_reads <- data_reads[rowSums(data_reads) > 0,]
    
    ## Filtered cell location
    if(length(centX)== 1 ) data_loca<- data_loca[rownames(data_reads),]
    
    ## return RCA object
    res <- methods::new("MolDiaISS",
                        data     = data_reads,
                        norm.data = matrix(, nrow = 0, ncol = 0),
                        scale.data = matrix(, nrow = 0, ncol = 0),
                        gene = colnames(data_reads),
                        cluster = factor(matrix(, nrow = 0, ncol = 0)),
                        location = data_loca,
                        cluster.marker = list(),
                        tsne.data = data.frame(matrix(, nrow = 0, ncol = 0)))
  }
  
  return(res)
}