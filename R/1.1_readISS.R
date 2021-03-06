######################################################################
##                       Read RCA ISS  data                         ##
######################################################################
"readISS"
#' Read ISS data
#' @description Read ISS data
#'
#' @param file File name in CSV format.Also data formate in "data.frame" class and "MolDiaISS"(Segment) & 'MolDiaISS_nonsegment'(Non-segment) 
#'             class (Output of \link[MolDia]{readISS}).
#' @param segment If the file is segmentated or not. TRUE/FALSE needed. Default if TRUE. See details for file structure,
#' @param cellid String to naming cell. Default is "CellID".
#' @param centX Name of X co-ordinate in file. Default is "centroidX"
#' @param centY Name of Y co-ordinate in file  Default is "centroidY"
#' @param genepos Name of genes to consider for gene positive cells. Default is NULL
#' @param geneposOPT Only work when 'genepos' has a value. "AND", "OR" and "NONE" condition for genepos. Default is "OR".
#' @param rpc Total reads per cell to be consider. Default is 1.
#' @param rpg Total reads per gene to be consider. Default is 1.
#' @param gene Gene names to be consider. Object in vector or list class. In list formated input, every list element is a group of
#'             interested genes. Every list element should have a name. Default is NULL.
#' @param nogene Gene name to exclude from data. Default is NULL.
#' 
#' @details File structure of non-segmentated input: 
#'          
#'          For non-segmentated file input, there should be 9 column. 7 column names are "Read", "Gene", "ParentCell",
#'          "Tile",  "MinAnchor",  "MinQuality" and  "MinAlign". Rest 2 column are x and y axix position and defined by 
#'          'centX' and 'centY' parameter.  
#'
#' @author Mohammad Tanvir Ahamed
#'
#' @return Output will be a object in class MolDiaISS or MolDiaISS_nonsegment. See detail \link[MolDia]{MolDiaISS} 
#'         and \link[MolDia]{MolDiaISS_nonsegment}
#'
#' @examples
#' ##### Reading ISS data from CSV formate
#' data_1 <- readISS(file = system.file("extdata", "CellBlobs_QT_0.35.csv", package="MolDia"),
#'                   cellid = "CellID", centX = "centroidX", centY = "centroidY")
#' data_2 <- readISS(file = system.file("extdata", "CellBlobs_QT_0.40.csv", package="MolDia"),
#'                   cellid = "CellID",centX = "centroidX", centY = "centroidY")
#' data_3 <- readISS(file = system.file("extdMinQualityata", "Hypocampus_left.csv", package="MolDia"),
#'                   cellid = "CellId",centX = "centroid_x", centY = "centroid_y")
#' data_4 <- readISS(file = system.file("extdata", "Hypocampus_right.csv", package="MolDia"),
#'                   cellid = "CellId",centX = "centroid_x", centY = "centroid_y")
#'
#' ###### Merge genes of interest into groups and plot
#' data(marker_gene)
#' data_5 <- readISS(file = data_4, gene = marker_gene)
#' result <- ISS_map(data = data_5, what = "cell", gene = data_5@gene[1:8]) 
#'
#'
#' ## Not RUN
#' ## Read RCA data in dataframe formate
#' #gene <- data_3@gene[1:3]
#' #data(single_cell)
#' #single_cell$CellID <- rownames(single_cell)
#' #Note: Big data. Take long time to load without selected gene
#' #data_sc <- readISS(file = single_cell, cellid = "CellID",gene= gene) 
#' 
#' ###### Reading non-segmentated file
#' data_6 <- readISS(file = system.file("extdata", "nonSeg_QT_0.35_0_details.csv", package="MolDia"), segment = FALSE,
#'                   centX = "PosX", centY = "PosY", nogene = "NNNN")
#' data_7 <- readISS(file = data_6, gene= c("Gdf7","WNT1","Pak3","Tfap2a"))
#' data_8 <- readISS(file = data_7, nogene = c("Gdf7","WNT1"))
#'
#'
#' @export
readISS <- function(file, segment = TRUE, cellid = "CellID", centX = NULL, centY = NULL, genepos= NULL, geneposOPT = "OR",
                    rpc = 1, rpg = 1, gene = NULL, nogene = NULL)
{
  
  ## Data type : Reading data from a specific location
  if(class(file)=="character")
  {
    ## Find segmentation
    if(segment == FALSE)
      { 
       ## Reading data
       my_file <- data.frame(data.table::fread(input = file , showProgress = TRUE))
       my_file <- my_file[,c("Read", "Gene", "ParentCell", "Tile", "MinAnchor", "MinQuality", "MinAlign", centX, centY)]

       ## Data without selected gene
       if(length(nogene) > 0 )
        {
        if(all(nogene%in%my_file$Gene)==FALSE) stop("Selected gene not present. Check gene name in 'nogene'.", call. = FALSE)
        my_file <- my_file[-which(my_file$Gene%in%nogene),]
        }
       
       ## Data with selected gene
       if(length(gene) > 0 )
        {
        if(all(gene%in%my_file$Gene)==FALSE) stop("Selected gene not present. Check gene name in 'gene'.", call. = FALSE)
        my_file <- my_file[which(my_file$Gene%in%gene),]
        }
       
       ## Final data
       res <- my_file
       
       ## Gene data
       gene <- res[, c("Gene")]
       names(gene) <- paste0("gene_", 1:length(gene))
       
       ## Reads data
       reads_data <- res[, c( "Read"), drop = FALSE]
       colnames(reads_data) <- c("Reads")
       rownames(reads_data) <- names(gene)
       
       ## Tile data
       tile_data <- res[, c( "ParentCell","Tile"), drop = FALSE]
       colnames(tile_data) <- c( "ParentCell","Tile")
       rownames(tile_data) <- names(gene)
       
       ## Quality data
       quality_data <- res[, c( "MinAnchor","MinQuality","MinAlign"), drop = FALSE]
       colnames(quality_data) <- c( "MinAnchor","MinQuality","MinAlign")
       rownames(quality_data) <- names(gene)
       
       ## Location data
       options(scipen = 999)
       loca_data <- res[, c("PosX","PosY"), drop = FALSE]
       colnames(loca_data) <- c( "centroid_x", "centroid_y")
       rownames(loca_data) <- names(gene)
       
       ## Return ISS object
       res <- methods::new("MolDiaISS_nonsegment",
                           reads    = reads_data,
                           gene     = gene,
                           tile     = tile_data,
                           quality  = quality_data,
                           location = loca_data)
       }
    
    if(segment == TRUE)
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
       #if(length(gene) > 0 )   data_reads <- data_reads[,gene, drop=FALSE]
       if(length(gene) > 0 )   
       {
         if(class(gene) == "character") data_reads <- data_reads[,gene, drop=FALSE]
         if(class(gene) == "list")      
         {
           data_reads <- ISS_sumstat (data = data_reads, gene = gene, stat = "sum")
           data_reads$total_reads <- NULL # Extra column added by ISS_sumstat function 
         }
       }
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
    
       ## Delete empty genes and empty cells 
       data_reads <- data_reads[,colSums(data_reads) > 0]
       data_reads <- data_reads[rowSums(data_reads) > 0,]
    
       ## Filtered cell location
       data_loca<- data_loca[rownames(data_reads),]
    
       ## Return RCA object
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
    #if(length(gene) > 0 )   data_reads <- data_reads[ ,  gene , drop=FALSE]
    if(length(gene) > 0 )   
    {
      if(class(gene) == "character") data_reads <- data_reads[,gene, drop=FALSE]
      if(class(gene) == "list")      
      {
        data_reads <- ISS_sumstat (data = data_reads, gene = gene, stat = "sum")
        data_reads$total_reads <- NULL # Extra column added by ISS_sumstat function 
      }
    }
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
    #if(length(gene) > 0 )   data_reads <- data_reads[ ,  gene , drop=FALSE]
    if(length(gene) > 0 )   
    {
      if(class(gene) == "character") data_reads <- data_reads[,gene, drop=FALSE]
      if(class(gene) == "list")      
      {
        data_reads <- ISS_sumstat (data = data_reads, gene = gene, stat = "sum")
        data_reads$total_reads <- NULL # Extra column added by ISS_sumstat function 
      }
    }
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
  
  ## Data type : MolDiaISS_nonsegment class
  if(class(file)=="MolDiaISS_nonsegment")
  {
    
    ## Data without selected gene
    if(length(nogene) > 0 )
    {
      if(all(nogene%in%file@gene)==FALSE) stop("Selected gene not present. Check gene name in 'nogene'.", call. = FALSE)
      gene <- file@gene[-which(file@gene%in%nogene)]
    }
    
    ## Data with selected gene
    if(length(gene) > 0 )
    {
      if(all(gene%in%file@gene)==FALSE) stop("Selected gene not present. Check gene name in 'gene'.", call. = FALSE)
      gene <- file@gene[which(file@gene%in%gene)]
    }
    
    ## Data with same gene
    if(length(gene) == 0 )
    {
      gene <- file@gene
    }
    
    ## Data
    reads    <- file@reads[names(gene), ,drop = FALSE] 
    gene     <- file@gene[names(gene)]
    tile     <- file@tile[names(gene), ,drop = FALSE]
    quality  <- file@quality[names(gene), ,drop = FALSE]
    location <- file@location[names(gene), ,drop = FALSE]
      
    ## Return ISS object
    res <- methods::new("MolDiaISS_nonsegment",
                        reads    = reads ,
                        gene     = gene ,    # gene,
                        tile     = tile  ,    # tile_data,
                        quality  = quality ,  # quality_data,
                        location = location ) # loca_data)
  }
  
  return(res)
  #return(gene)
}