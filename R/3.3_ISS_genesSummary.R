######################################################################
##                       RCA gene Summary                          ##
######################################################################
"ISS_genesSummary"
#' Summary of specific gene of interest in ISS data
#'
#' @description Summary of specific gene of interest of RCA data
#' @param data Input data in class MolDiaISS. Output of \link[MolDia]{readISS}
#' @param gene gene of interest
#'
#' @examples
#' ex_data <- readISS(file = system.file("extdata", "CellBlobs_QT_0.35.csv", package="MolDia"),
#'                   cellid = "CellID", centX ="centroidX", centY="centroidY" )
#' res <- ISS_genesSummary(data = ex_data, gene = "Actb.L")
#'
#' @return Barplot
#'
#' @export
ISS_genesSummary <- function(data, gene = NULL)
  
{
  ## Select gene
  if(length(gene) == 0) stop("Please select genes of interest", call. = FALSE)
  if(length(gene) > 1) stop("Maximum length of gene should be 1", call. = FALSE)
  if(any(gene%in%data@gene)==FALSE) stop("Genes not present in data. Check gene name", call. = FALSE)
  
  ## Main data
  main_data <- data@data
  main_data <- main_data[,order(colnames(main_data))]
  mydata  <- main_data
  mydata1 <- mydata > 0
  
  gene_lst <- as.list(colnames(mydata))
  
  ## Count frequency number of other gene in gene+ cell
  res_0      <- lapply(gene_lst, function(object)
  {
    data1 <- mydata1[mydata1[,object],]
    data2 <- colSums(data1)
    #data2 <- round((data2/max(data2))*100,10)
    data2[which(data2 == max(data2))] <- max(data2)
    data2
  })
  names(res_0) <- unlist(gene_lst)
  res_0 <- res_0[sort(names(res_0))]
  res_0 <- do.call(rbind,res_0)
  
  ## Count percentage of other gene in gene+ cell
  res      <- lapply(gene_lst, function(object)
  {
    data1 <- mydata1[mydata1[,object],]
    data2 <- colSums(data1)
    data2 <- round((data2/max(data2))*100,10)
    data2[which(data2 == max(data2))] <- max(data2)
    data2
  })
  names(res) <- unlist(gene_lst)
  res <- res[sort(names(res))]
  res_1 <- do.call(rbind,res)
  
  ### Reshape data
  # For count
  heat_data_0 <- data.frame(reshape::melt(res_0))
  
  # For precentage
  heat_data <- data.frame(reshape::melt(res_1))
  
  ## Merge both count and percentage data
  heat_data<- merge(heat_data, heat_data_0, by = c("X1", "X2"))
  colnames(heat_data) <- c("X1","X2","value","value_count")
  
  ## Subset if gene of interest
  heat_data_1 <- split(heat_data, heat_data$X2)
  heat_data_1 <- heat_data_1[match(gene, names(heat_data_1))]
  heat_data_1 <- heat_data_1[[1]]
  
  ## Get 100% gene count
  gene100 <- heat_data_1$value_count[which(heat_data_1$X1==heat_data_1$X2)]
  
  ## Create group in genes
  heat_data_1 <- heat_data_1[-which(heat_data_1$X1==heat_data_1$X2),]
  heat_data_1$X2 <- paste0("(B) ", gene, " in all cells")
  
  heat_data_2 <- split(heat_data, heat_data$X1)
  heat_data_2 <- heat_data_2[match(gene, names(heat_data_2))]
  heat_data_2 <- heat_data_2[[1]]
  heat_data_3 <- heat_data_2
  
  heat_data_3$X1 <-  heat_data_2$X2
  heat_data_3$X2 <-  heat_data_2$X1
  heat_data_3$value <- heat_data_2$value
  heat_data_3    <- heat_data_3[-which(heat_data_3$X1==heat_data_3$X2),]
  heat_data_3$X2 <- paste0("(A) All genes in ",gene,"+ Cells")
  
  heat_data_3 <- rbind(heat_data_1,heat_data_3)
  heat_data_3$X2 <- as.factor(heat_data_3$X2)
  heat_data_3 <- data.frame(heat_data_3)
  
  ## GGplot
  p <- ggplot2::ggplot(data = heat_data_3, ggplot2::aes_string(x = "X1", y = "value")) +
    ggplot2::geom_bar(ggplot2::aes_string(fill = "value"),stat="identity") +
    ggplot2::geom_text(stat="identity",ggplot2::aes_string(label="value_count"),angle= 90, hjust = 0) +
    ggplot2::facet_grid(X2 ~ ., scales = "free_y") +
    ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 90, hjust = 1, vjust= 0.5, size = 10, face="bold"),
                   axis.title.x = ggplot2::element_blank()) +
    ggplot2::theme(strip.text.y = ggplot2::element_text(size = 15)) +
    ggplot2::scale_fill_gradientn(name = "Cells %", colours = c("darkgreen", "yellow", "red")) +
    ggplot2::labs (y= "% of Cells", title = paste0("Distibution of ", gene, " gene"), subtitle = paste0("Total ", gene, "+ cells is ", gene100))
  
  print(p)
  
  return(NULL)
}