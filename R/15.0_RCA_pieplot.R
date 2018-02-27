"ISS_pieplot"
#' Venn-pie chart on RCA data based on genes of interest
#'
#' @param data Input data in class RCA_class. Output of \link[MolDia]{readRCA}.
#' @param gene Gene of interest. Default is NULL.
#' @param with_gene From gene of interest , select only cells only with these genes.
#' @param without_gene From gene of interest , select only cells only without these genes.
#'
#' @return The Venn-pie plot.
#'
#'         Number of cells having specific gene profile.
#' @author Mohammad Tanvir Ahamed
#'
#' @examples
#'########## Reading data
#' data_1      <- readRCA(file = system.file("extdata", "CellBlobs_QT_0.35.csv", package="MolDia"),
#'                        cellid = "CellID",centX = "centroidX", centY = "centroidY")
#'
#' ## Define marker gene group
#' marker_gene <- data_1@gene
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
#' mydata <- RCA_barplot(data = data_1, gene = marker_gene1, gene.target = 1,
#'                       at.least.gene = 3, show.same.gene = FALSE)
#'
#' #####  Pie chart of all Circular mRNA in cell with at least 3 neuronal genes.
#' genes <- marker_gene1[[1]]
#' res   <- ISS_pieplot (data = mydata, gene = genes)
#'
#' @export
#' 
ISS_pieplot <- function(data, gene = NULL, with_gene = NULL, without_gene = NULL)
{
  ## Main data 
  main_data <- data
  data <- data@data
  
  ## Select gene with or without names
  if(length(with_gene) > 0 )    data <- data[rowSums(data[,with_gene, drop = FALSE]) > 0,]
  if(length(without_gene) > 0 ) data <- data[rowSums(data[,without_gene, drop = FALSE]) == 0,]
  
  ## Select genes
  if (length(gene)==0 ) data <- data
  else {
    data <- data[,gene, drop= FALSE]}
  
  ## Check gene name in data
  if(all(gene %in%colnames(data))==FALSE) stop("Check gene name are same and equel length in both input ", call. = FALSE)
  
  ## Delete cell with no gene
  data <- data[rowSums(data) > 0,]
  data <- data[,colSums(data) > 0]
  
  ## Gene name
  gname <-colnames(data)
  
  ## Get combinition of gene name (Gene count > 0) for each cell and count frequency of gene combinition
  data <- apply(data,1,function(i)
  {
    pp<- paste0(colnames(data)[i>0],collapse = "-")
  })
  data <- data.frame(table(unlist(data)))
  colnames(data) <- c("Gene","Count")
  data <- data[order(data$Count, decreasing = T),]
  
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
  
  ## Apply Venn-pie plot
  if(length(gene) == 0 )
  {
    set.seed(1)  ## This random seed will keep the constant color reproduce
    res <- sunburstR::sunburst(data, count = T, legendOrder = sort(gname))
  } else {
    set.seed(1)
    cols <- as.list(randomcoloR::distinctColorPalette(length(gene)))
    names(cols) <- NULL
    
  ## Custom Messege
    custom.message = "function (d) {
    root = d;
    while (root.parent) {
    root = root.parent
    }
    p = (100*d.value/root.value).toPrecision(3);
    msg = p+' %<br/>'+d.value+' of '+root.value;
    return msg;
  }"
    
    
    
    set.seed(1)  ## This random seed will keep the constant color reproduce
    res <- sunburstR::sunburst(data, count = T, percent = T,explanation = custom.message,
                               legendOrder = sort(gname),
                               colors = list(range  = cols,
                                             domain = gene))
}
  print(res)
  
  ## return RCA object
  myres <- methods::new("RCA_class",
                        data     = data,
                        location = main_data@location[rownames(data),],
                        gene = colnames(data))
  return(data)
  
}

################################################################################
RCA_pieplot <- function(data, gene = NULL, gene_with = NULL, gene_without = NULL)
  {
  ## Data select
  main_data_1 <- data
  main_data_2 <- data@data
  data        <- main_data_2

  ## Select gene names
  if(length(gene_with) > 0 )    data <- data[rowSums(data[,gene_with, drop = FALSE]) > 0,]
  if(length(gene_without) > 0 ) data <- data[rowSums(data[,gene_without, drop = FALSE]) == 0,]

  ## Check gene names in main data
  if(all(gene %in%colnames(data))==FALSE) stop("Check gene name are same and equel length in both input ", call. = FALSE)

  ## Gene selection
  if (length(gene)==0 ) data <- data
  else {
    data <- data[,gene]}

  ## Delete cell with no gene
  data <- data[rowSums(data) > 0,]
  data <- data[,colSums(data) > 0]

  ## Save data to return
  newdata<- data

  ## Gene name
  gname <-colnames(data)

  ## Get combinition of gene name (Gene count > 0) for each cell and count frequency of gene combinition
  data <- apply(data,1,function(i)
    {
     pp<- paste0(colnames(data)[i>0],collapse = "-")
    })
  data <- data.frame(table(unlist(data)))
  colnames(data) <- c("Gene","Count")
  data <- data[order(data$Count, decreasing = T),]

  ## Sort gene name in each gene name accorging their high to low occurance in gene combinition
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

  ## Apply Venn-pie plot
  if(length(gene) == 0 )
  {
  set.seed(100000009)  ## This random seed will keep the constant color reproduce
  res <- sunburstR::sunburst(data, count = T)
  } else {
    cols <- as.list(randomcoloR::distinctColorPalette(length(gene)))
    names(cols) <- NULL

    set.seed(100000009)  ## This random seed will keep the constant color reproduce
    res <- sunburstR::sunburst(data, count = T,
                               colors = list(range  = cols,
                                             domain = gene))
    }
   print(res)

   ## return RCA object
   myres <- methods::new("RCA_class",
                         data     = newdata,
                         location = main_data_1@location[rownames(newdata),],
                         gene = colnames(newdata))
   return(myres)
}
