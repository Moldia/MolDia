######################################################################
##                   RCA data Filter : Barplot                      ##
######################################################################
"ISS_barplot"
#' Plot barplot on ISS data bsed on different condition
#' @description Plot bar plot on RCA data bsed on different condition
#'
#' @param data Input data in class MolDiaISS. Output of \link[MolDia]{readISS}.
#' @param gene Object in vector or list formate. In list formated input every list element is a group of
#'        interested genes.
#' @param total.expr Scale factor to re-scale the data. Default is 1e4.
#' @param gene.target Only applicable when 'gene' parameter is a list. Which gene group in 'gene' parameter
#'        will consider for special operation. Default is NULL. See details.
#' @param gene.show Only applicable when 'gene' parameter is a list. Which group/groups of gene in 'gene' parameter
#'        to show. Default is NULL (Show all groups). See details
#' @param target.min.count.cell Minimum number of reads per cell to consider in targated gene group in
#'        'gene.target' parameter. Default is 1.
#' @param rest.min.count.cell Minimum number of reads per cell to consider in the genes that not consider in
#'        'gene.target' parameter. Default is 1.
#' @param min.count Only applicable when 'gene.target' has a value. Tharshold to consider as minimum count of a gene
#'        in a single cell to consider as gene expression. Default is 0. See details.
#' @param at.least.gene Only applicable when 'gene.target' has a value. Minimum number of genes expressed in a cell
#'        to consider. Default is 0, means all.
#' @param at.most.gene Only applicable when 'gene.target' has a value. Maximum number of genes expressed in a cell
#'        to consoder
#' @param show.same.gene Consider to show same gene or not. Default is FALSE. See details.
#' @param str.same.gene Define string that will make the difference between different group in same gene name.
#'                      Only active if show.same.gene = TRUE. See detail.
#' @param main Main title of the plot.
#'
#' @details 'gene.target' parameter will only work when 'gene' parameter is a list. 'gene.target' will consider
#'          only one gene group from gene list for special operation. Defaule value is NULL, means it will not
#'          consider any specific gene group.
#'
#'          'gene.show' by default (NULL) consider all gene group in 'gene' parameter. But one can choose
#'          which gene group to be consider. Example: c(3,5) means, 3rd and 5th group of genes will be
#'          consider to show or analysis.
#'
#'          'min.count' is the tharshold to consider as minimum count of a gene per single cell to consides as
#'          expression. Default valus is 1, which means, 1 reads count  in a sigle cell for a gene will consider
#'          as expression of that specific gene in that cell.
#'
#'          'show.same.gene' will consider genes with same name. This case may happen if same gene has different form
#'          of detections. In this case the naming formte will be [gene name].[different form]. Gene name and different
#'          form name should be isolated by "." See 'str.same.gene' parameter for details.
#'
#'          'str.same.gene' is a string or a vector of string that will make the different group in same gene name.
#'          Example: gene Adar has two different miRNA which is circular (Adar.C) and linear (Adar.L). So the value
#'          of str.same.gene will be c(".C",".L").
#'
#' @return Data in class \link[MolDia]{readISS}.
#'
#'         The barplot will explain the total number of reads per gene in 'total.expr' (Default 1e4) cells. The yellow line
#'         inside the bar plot indicate the number of expressed cells and percentage count on each bar is the percent of cells
#'         express that gene.
#' @importFrom graphics barplot lines text legend
#' @author Mohammad Tanvir Ahamed
#' @examples
#' ########## Reading data
#' data_1      <- readISS(file = system.file("extdata", "CellBlobs_QT_0.35.csv", package="MolDia"),
#'                        cellid = "CellID", centX = "centroidX", centY =  "centroidY")
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
#'
#' marker_gene1 <- list(Neuron = Neuron,
#'                      Nonneuron = Nonneuron,
#'                      Cir_rna = c_rna,
#'                      Lin_rna = c(setdiff(setdiff(l_rna,Neuron),Nonneuron)))
#'
#' # Barplot
#' all_data     <- ISS_barplot(data = data_1, gene = marker_gene,gene.target = NULL)
#' neuron_group <- ISS_barplot(data = data_1, gene = marker_gene1, gene.target = 1,
#'                             target.min.count.cell = 2, rest.min.count.cell = 2,
#'                             at.least.gene = 1, gene.show = NULL, main = " Neuron group (>=1 genes)")
#'
#'
#'
#' @export
ISS_barplot <- function(data, gene, total.expr = 1e4, gene.target = NULL, gene.show = NULL, target.min.count.cell = 1, rest.min.count.cell = 1,
                        min.count = 1, at.least.gene = 0, at.most.gene = ncol(data), show.same.gene = FALSE,
                        str.same.gene = c(".C",".L"), main = "")
{
  ## Extract data from MolDiaISS object
  main_data <- data
  data      <- main_data@data

  ## Reproduce same result with same color
  set.seed(100000009)

  ## Get Targated gene
  #if(length(gene.target) == 0) {gene.target <- (unlist(gene))
  #} else {
  #  gene.target <- unlist(gene[gene.target])}


  ## Check if length of minimum gene (at.least.gene) is greater than length of target gene (gene.target)
  if(length(gene.target) >0 )
    {
    if(at.least.gene > length(gene[[gene.target]]))
    {stop("Minimum number of gene (at.least.gene) cant exceed the length of targated gene (gene.target)", call. = FALSE)}
    }
  ## Check if gene name in both input are same and equel
  # Check duplicated gene name in gene list
  if(any(duplicated(unlist(gene)))== TRUE) stop ("Check duplicated gene names in gene list.", call. = FALSE )

  # Check gene name exist in main data
  if (all(unlist(gene)%in%colnames(data))== FALSE)
  {
    gnam<- unlist(gene)[which(!unlist(gene)%in%colnames(data))]
    stop (paste0("Gene name: "), paste0(gnam, collapse = ", "), " not present in data. Check provided gene list.", call. = FALSE)
  }

  ## Select cell that have only provided gene
  data <- data[,unlist(gene), drop = FALSE]
  data <- data[rowSums(data) > 0,,drop= FALSE]

  ## Filter Targated gene
  if(length(gene.target) == 0) data <- data
  if(length(gene.target) > 0)
  {

    gene.target <- unlist(gene[gene.target])

    data_1 <- data[,gene.target, drop = FALSE]
    #data_1 <- data[,gene[[gene.target]], drop = FALSE]
    data_1 <- data_1[rowSums(data_1) >= target.min.count.cell,,drop = FALSE]


    ### Define user required read thrashold to define the level of expression
    data_3 <- rowSums(data_1 >= min.count)
    gene_cell <- names(data_3)[which(data_3 %in% c(at.least.gene:at.most.gene))]
    data     <- data[gene_cell,,drop = FALSE]


  ## Gene reads thrashhold per cel
  #data <- data[rowSums(data) >= min.count.cell,]

  #data_gene.target <- data[, which(colnames(data)%in%gene.target[[1]])]
  data_gene.rest   <- data[,-which(colnames(data)%in%gene.target)]
  data_gene.rest   <- data_gene.rest[rowSums(data_gene.rest) >= rest.min.count.cell,]
  data <- data[rownames(data_gene.rest),]
  }

  ## Select gene group
  if (length(gene.show) == 0 ) {gene.show <- c(1:length(gene))
  }else {gene.show <- gene.show}
  gene <- gene[gene.show]

  ## Split data according to gene group order and again merge them accordingly
  res1 <- lapply(gene,function(i){ res2<- data[,i,drop = FALSE]})
  names(res1)<- NULL
  res1 <- do.call(cbind, res1)

  ## Show only same gene in different group
  if (show.same.gene == TRUE)
  {
    gene_1  <- unlist(lapply(strsplit(colnames(res1),"\\."),"[[",1))
    gene    <- sort(unique(gene_1[duplicated(gene_1)]))
    gene    <- sort(colnames(res1)[which(gene_1%in%gene)])
    res1    <- res1[,gene]

    color_group1 <- randomcoloR::distinctColorPalette(length(str.same.gene))
    color_group2 <- NA

    for(i in 1: length(str.same.gene))
    {
      color_group2[grepl(str.same.gene[i],colnames(res1))]<- color_group1[i]
    }
  }

  ## Define color for each group of gene
  if (show.same.gene == FALSE)
  {
    color_group1 <- randomcoloR::distinctColorPalette(length(gene))
    color_group2 <- rep(color_group1, unlist(lapply(gene,length)))
  }

  ## Final Output
  final_output <- res1
  final_output <- final_output[,colSums(final_output)>0, drop = FALSE]
  final_output <- final_output[rowSums(final_output)>0,, drop = FALSE]

  ## Change total number of cell (number of cell adjust)
  data_totalexp           <- (res1/nrow(res1))*total.expr
  data_totalexp_colsums   <- round(colSums(data_totalexp))

  ## Barplot
  bp1<-barplot(data_totalexp_colsums, col = color_group2, las = 2,
               legend.text = if( show.same.gene == TRUE) str.same.gene else names(gene),
               ylab = "Total reads counts", xlab = "", cex.names = 1,
               ylim = c(1, max(data_totalexp_colsums)+(max(data_totalexp_colsums)/6)),
               args.legend = list(x = "top", bty = "n",fill = color_group1),
               main = paste0(main, " \nCells: ", nrow(final_output), " Genes: ",
                             ncol(final_output),"\nR/C:", round(sum(final_output)/nrow(final_output),1),
                             " R/G:", round(sum(final_output)/ncol(final_output),1),
                             " Sparcity:", round(sum(final_output == 0)/(dim(final_output)[1]*dim(final_output)[2]),2)))

  tnc  <- round((colSums(res1 > 0)/nrow(res1))*total.expr)
  lines(cbind(bp1,tnc),col = "orange", lwd = 2)
  text(x = bp1+0.3, y = data_totalexp_colsums, label = paste0("           ",round((tnc/total.expr)*100,2)," %"),
       adj= 0.5,col = "red", cex = 1,srt= 90, pos = 3)
  
  
  #final_output

  ## return RCA object
  #res <- methods::new("MolDiaISS",
  #                    data     = final_output,
  #                    normdata = main_data@normdata[rownames(final_output),],
  #                    scaledata = main_data@scaledata[rownames(final_output),],
  #                    location = main_data@location[rownames(final_output),],
  #                    gene = colnames(final_output))

  ## Return result
  main_data@data <- final_output
  if(length(main_data@norm.data)  > 0 ) main_data@norm.data  <- main_data@norm.data[rownames(final_output),colnames(final_output)]
  if(length(main_data@scale.data) > 0 ) main_data@scale.data <- main_data@scale.data[rownames(final_output),colnames(final_output)]
  main_data@location <- main_data@location[rownames(final_output),]
  main_data@gene <- colnames(final_output)

  #return(main_data)
  res<- list( count =data_totalexp_colsums, gene = gene, tnc = tnc)
  #return(res)
  #return(gene_group)
  return(main_data)
}