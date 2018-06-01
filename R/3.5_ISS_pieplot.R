"ISS_pieplot"
#' Venn-pie chart on ISS data based on genes of interest
#'
#' @param data Input data in class MolDiaISS. Output of \link[MolDia]{readISS}.
#' @param gene Gene of interest. It could be a vector of gene or a list of genes groups. Default is NULL means takes all genes
#'             in the given dataset.
#' @param with_gene From gene of interest , select only cells only with these genes.
#' @param without_gene From gene of interest , select only cells only without these genes.
#' @param segmentcol Hex color code of the overlap segment.
#' @param randomseed Ramdom seed to change color schema.
#' @param colorpalette Default is NULL. A vector of value . Can be applicabile with selected element with gene list. 
#' 
#' @note Maximum 256 genes (256 color) can be ploted with this version at this moment.
#'
#' @return The Venn-pie plot.
#'
#'         Number of cells having specific gene profile.
#' @author Mohammad Tanvir Ahamed
#'
#' @examples
#'########## Reading data
#' ex_data      <- readISS(file = system.file("extdata", "CellBlobs_QT_0.35.csv", package="MolDia"),
#'                        cellid = "CellID",centX = "centroidX", centY = "centroidY")
#'
#' ## Define marker gene group
#' marker_gene <- ex_data@gene
#' c_rna <- marker_gene[grepl(".C",marker_gene)] # Circular RNA
#' l_rna <- setdiff(marker_gene,c_rna)           # Linear RNA
#'
#' ## Neuronal marker
#' Pyramidal   <- c("Nrn1.L","Pcp4.L")
#' Interneuron <- c("Sst.L","Pvalb.L","Ndnf.L","Vip.L","Sncg.L","Cck.L")
#' Other       <- c("Map2.L")
#' Neuron      <- c(Pyramidal,Interneuron,Other)
#'
#' ## Non-Neuronal marker
#' Oligodendrocyte <- c("Plp1.L","Enpp2.L")
#' Astrocytes  <- c("Gfap.L","Mfge8.L","Aldh1l1.L","S100b.L")
#' Nonneuron   <- c(Oligodendrocyte,Astrocytes)
#' marker_gene1 <- list(Neuron = Neuron,
#'                      Nonneuron = Nonneuron,
#'                      Cir_rna = c_rna,
#'                      Lin_rna = c(setdiff(setdiff(l_rna,Neuron),Nonneuron)))
#'
#' ####  Select cell with at least 3 neuronal genes
#' mydata <- ISS_barplot(data = ex_data, gene = marker_gene1, gene.target = 1,
#'                       at.least.gene = 3, show.same.gene = FALSE)
#'
#' #####  Pie chart of all Circular mRNA in cell with at least 3 neuronal genes.
#' genes <- marker_gene1[[1]]
#' res   <- ISS_pieplot (data = mydata, gene = genes, without_gene = "QKI.L", segmentcol = "EEEEEE")
#'
#' @export
#' 
ISS_pieplot <- function(data, gene = NULL, with_gene = NULL, without_gene = NULL, segmentcol = "EEEEEE", randomseed = 10, colorpalette = NULL)
{
  ################################## Data preparation 
  ## Check for at least 2 genes
  if(length(gene) == 1) stop("There should be at least 2 genes to be selected", call. = FALSE)
  
  ## Select main data
  main_data <- data
  data <- data@data
  
  ## Check if gene is in list formate input
  if(is.list(gene) == TRUE)
  {
    ## Check if all gene in the list is in data
    if(all(unlist(gene) %in%colnames(data))==FALSE) stop("Check gene name are same in both input data and supplied gene list", call. = FALSE)
    data <- lapply(seq_along(gene), function(i)
    {
      pp<- data.frame(rowSums(data[,gene[[i]],drop = FALSE]))
      colnames(pp) <- names(gene[i])
      pp
    }
    )
    data <- do.call(cbind,data)
    gene <- names(gene)
  }
  
  ## Select with_gene(+) and without_gene(-) enriched cells
  if(length(with_gene) > 0 )
  {
    data <- data[rowSums(data[,with_gene,    drop = FALSE]) >  0,,drop = FALSE]
  }
  if(length(without_gene) > 0 ) 
  {
    data <- data[rowSums(data[,without_gene, drop = FALSE]) == 0,,drop = FALSE]
  }
  
  ## Check gene name in data
  if (length(gene)==0 ) data <- data
  else {
    if(all(gene %in%colnames(data))==FALSE) stop("Check gene name are same and equel length in both input ", call. = FALSE)
    data <- data[,gene, drop= FALSE]
  }
  
  ## Delete row with cell with no genes
  data <- data[rowSums(data) > 0,]
  data <- data[, apply(data, 2, sum)!=0]
  
  ################################# Data convertion into sunbrust format
  ## Convery data into sunbrust formated file
  data <- data2sunbrust(data)
  res_data <- data[order(data$Count, decreasing = TRUE),]
  rownames(res_data) <- NULL
  
  ## Optimised for sunbrust plot
  data$Gene <- paste0(data$Gene,"- ")
  
  ## Gene names
  #gname <- sort(unique(unlist(strsplit(data$Gene,"-"))))
  #gname <- c(gname[-1],gname[1])
  gname <- c(gene," ")
  
  ## Define color
  set.seed(randomseed)
  #mypalette <- as.list(randomcoloR::distinctColorPalette(length(gname)-1))
  mypalette <- randomcoloR::distinctColorPalette(k = 256)
  if(length(colorpalette) == 0)
    { mypalette <- as.list(mypalette[1:(length(gname)-1)])
  } else 
      { mypalette <-as.list(mypalette[colorpalette])}
  
  mypalette[length(gname)+1] <- segmentcol
  mypalette <- unlist(mypalette)
  #names(mypalette) <- NULL
  
  ## Custom Messege
  custom.message = "function (d) {
  root = d;
  while (root.parent) {
  root = root.parent
  }
  p = (100*d.value/root.value).toPrecision(3);
  msg = p+' %<br/>'+d.value+' of '+root.value;
  return msg; }"
  
  ## Sunbrust plot 
  set.seed(randomseed)
  res <- sunburstR::sunburst(data, count = T, legendOrder = gname[-length(gname)], 
                             colors = list(range  = mypalette,
                                           domain = gname),
                             explanation = custom.message)
  rm(.Random.seed, envir=globalenv())
  print(res)
  
  ## Return
  return(res_data)
}

### Function to convert data frame to sunbrust formate data set
data2sunbrust <- function(data)
{
  ## Gene name
  data <- data[,names(sort(colSums(data), decreasing = T))] 
  gname <- colnames(data)
  #data  <- data[,gname]
  
  ## Get combinition of gene name (Gene count > 0) for each cell and count frequency of gene combinition
  data <- apply(data,1,function(i)
  {
    pp<- paste0(colnames(data)[i>0],collapse = "-")
  })
  data <- data.frame(table(unlist(data)))
  colnames(data) <- c("Gene","Count")
  #data <- data[order(data$Count, decreasing = T),]
  
  ## Sort gene name in each gene name according their high to low occurance in gene combinition
  tot_comb  <- as.vector(data$Gene)
  tot_comb1 <- strsplit(tot_comb,"-")
  tot_comb2 <- unique(unlist(tot_comb1))
  
  g  <- rep(seq_along(tot_comb1), sapply(tot_comb1, length))
  m  <- lapply(tot_comb2, function(x) g[which(unlist(tot_comb1) %in% x)])
  names(m) <- tot_comb2
  m1 <- sort(unlist(lapply(m,function(i) sum(data$Count[i]))), decreasing = TRUE)
  m2 <- names(m1)
  
  ll <- lapply (tot_comb1, function(i) { m2[sort(match(i,m2))] })
  ll1 <- unlist (lapply(ll, paste0, collapse= "-"))
  
  data$Gene <- ll1
  
  ## Return data
  return(data)
}