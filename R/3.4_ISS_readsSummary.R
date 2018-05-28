######################################################################
##                       RCA data Summary                          ##
######################################################################
"ISS_readsSummary"
#' Summary of ISS data and expression
#'
#' @param data Input data in class MolDiaISS. Output of \link[MolDia]{readISS}
#' @param readlimit Total number of reads per cell to consider
#' @param text.size Text size of Gene name.
#' @param intervel.dep Vector of values between 0 and 100 , where you want to set interval. If one want to define the range like
#'        0 to 5, 5 to 20, 20 to 40, 40 to 60, 60 to 100 then the input valuse will be c(5,20,40,60)
#' @param coexp Same gene coexpression. Default is FALSE. TRUE means at least 2 reads.
#'
#' @return Summary of cells after reads delete.
#'
#' @importFrom graphics axis par plot
#' @importFrom stats lm var
#'
#' @author Mohammad Tanvir Ahamed
#'
#' @examples
#' data_1 <- readISS(file = system.file("extdata", "CellBlobs_QT_0.35.csv", package="MolDia"),
#'                   cellid = "CellID", centX = "centroidX", cent = "centroidY")
#' res    <- ISS_readsSummary(data = data_1, readlimit = 10, text.size = 6, intervel.dep = NULL, coexp = TRUE)
#'
#' @export
ISS_readsSummary <- function(data, readlimit = 10, text.size = 6, intervel.dep = NULL, coexp = FALSE)
{
  ## Gurbage cleaning
  gc()

  ## Save main data
  main_data <- data@data
  data   <- main_data
  data1  <- main_data
  
  ## Check maximum readsper cell
  
  rl <- max(rowSums(main_data))
  if(readlimit > rl ) 
  {
    cat("Maximum number of reads/cell is ", rl,"\n")
    readlimit <- rl
  }
  
  res  <- matrix(NA, ncol = 5, nrow =readlimit+1 )
  gene <- matrix(NA, ncol = ncol(data), nrow =readlimit+1 )
  colnames(gene) <- colnames(data)
  for(i in 0 : readlimit){
    if (i == 0 )
    { res [i+1,1] <- i
    res [i+1,2] <- nrow(data)
    res [i+1,3] <- round((nrow(data)/ nrow(data1))*100,2)
    res [i+1,4] <- sum(data)
    res [i+1,5] <- round((sum(data) / sum(data1))*100,2)
    gene[i+1,] <- round((colSums(data)/colSums(data1))*100,2)
    }
    if( i > 0 ){
      mm <- which(rowSums(data1) < i)
      if (length(mm) == 0 ) data <- data1
      else {
        data <- data1[-which(rowSums(data1) < i),]
      }
      res [i+1,1] <- i
      res [i+1,2] <- nrow(data)
      res [i+1,3] <- round((nrow(data)/ nrow(data1))*100,2)
      res [i+1,4] <- sum(data)
      res [i+1,5] <- round((sum(data) / sum(data1))*100,2)
      gene[i+1,] <- round((colSums(data)/colSums(data1))*100,2)
    }
  }
  res <- data.frame(res)
  colnames(res) <- c("readDel", "cell","%Cell", "totalreads","%Reads")
  #res <- cbind(res, gene)

  ## Define data
  data <- data1

  ## Multiple plot
  #op<- par()
  #op <- par(mfrow=c(3,2), xpd = TRUE)

  ## Number of reads per gene
  res1<- apply(data,2,sum)
  #barplot(res1, las = 2, ylim = c(0,max(res1)), main = "Number of Reads per gene", cex.names = 0.9, ylab = "Reads count")
  df <- reshape::melt(res1)
  df$gene<- rownames(df)
  p1 <-ggplot2::ggplot(data=df, ggplot2::aes_string(x= "gene", y= "value")) +
       ggplot2::geom_bar(stat="identity") +
       ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 90, hjust = 1, vjust= 0.5, size = text.size,face="bold")) +
       ggplot2::ggtitle("Number of reads per gene") +
       ggplot2::ylab("Read count") +
       ggplot2::xlab("")
  #p1

  ## Plot frequency of genes in cell
  #res11 <- data >0
  #res11 <- apply(res11,1,sum)
  #hist(res11, breaks = unique(res11), xlab = "Number of genes" , ylab = "% of cells share", main = "Cell share number of genes")

  ## Plot cells per genes : % of cell per gene
  res2 <- data > 0
  res2 <- apply(res2,2,sum)
  res2 <- (res2 / nrow(data))*100
  #barplot(res2, las = 2, ylim = c(0,max(res2)), main = "% of Cells per gene", cex.names = 0.9, ylab = "% of cells" )
  df <- reshape::melt(res2)
  df$gene<- rownames(df)
  p2 <-ggplot2::ggplot(data=df, ggplot2::aes_string(x= "gene", y= "value")) +
       ggplot2::geom_bar(stat="identity") +
       ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 90, hjust = 1, vjust= 0.5, size = text.size,face="bold")) +
       ggplot2::ggtitle("% of cells per gene") +
       ggplot2::ylab("% of cells") +
       ggplot2::xlab("")
  #p2

  ## Plot Reads per cell: Average number of reads per cell per gene
  res3 <- data > 0
  res3 <- apply(res3,2,sum)
  res4 <- res1 / res3
  #barplot(res4, las = 2, ylim = c(0,max(res4)), main = "Average number of reads per cell per gene", cex.names = 0.9)

  df <- reshape::melt(res4)
  df$gene<- rownames(df)
  p3 <- ggplot2::ggplot(data=df, ggplot2::aes_string(x= "gene", y= "value")) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 90, hjust = 1, vjust= 0.5, size = text.size,face="bold")) +
        ggplot2::ggtitle("Average number of reads per cell per gene") +
        ggplot2::ylab("") +
        ggplot2::xlab("")
  #p3

  ## plot relation between cell number with reads number
  #res5 <- data > 0
  #res5 <- apply(res5,2,sum)
  #res6 <- cbind(res1 ,res5)
  #res6 <- res6[order(res6[,1]),]

  #fit <-  lm(res6[,1] ~ res6[,2])
  #r <- round(summary(fit)$r.squared,2)
  #plot(res6, xlab = "Reads per gene", ylab = "Cells per gene", main = bquote(.("Relation between cell & reads per gene") ~ R^2==.(r)))

  # Frequency distrubution of number of reads per cell : Cells frequency distri bution by number of cells
  res5 <- table(rowSums(data))
  #barplot(res5, xlab = "Number of reads", ylab = "Number of cells", main = "Cell's frequency distribution by number of reads")
  df <- reshape::melt(res5)
  p4 <- ggplot2::ggplot(data=df, ggplot2::aes_string(x= "Var.1", y= "value")) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 0, hjust = 1, vjust= 0.5)) +
        ggplot2::ggtitle("Cell frequency distribution by number of cells") +
        ggplot2::ylab("Number of cells") +
        ggplot2::xlab("Reads")
  #p4


  gc()
  ## Ploting genes per reads deleted
  #library(lattice)
  #pp<- levelplot(gene, xlab = "Deleted reads", ylab = "Genes percentage")
  #print(pp)
  #plotrix::color2D.matplot(gene,c(1,0),c(0,0),c(0,1),show.legend = TRUE ,
  #                xlab="",ylab="Reads deleted",main="% of reads after reads delete", axes = FALSE)

  #axis(side = 1, at = (1:ncol(data1))-0.5, labels=colnames(data), cex.axis = 0.9, las = 2)
  #axis(side = 2, at = 0: readlimit, labels= readlimit:0, cex.axis = 0.7, las = 2)

  hm.palette  <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')), space='Lab')
  df <- reshape::melt(gene)
  p5 <- ggplot2::ggplot(data =df, ggplot2::aes_string(x = "X2", y= "X1") ) +
        ggplot2::geom_tile(ggplot2::aes_string(fill = "value")) +
        ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 90, hjust = 1, vjust= 0.5, size = text.size, face="bold"),
                       axis.title.x = ggplot2::element_blank()) +
        ggplot2::coord_cartesian(ylim = c(1, max(df$X1)))+
        ggplot2::scale_fill_gradientn(name = "Reads %",colours = hm.palette(100)) +
        ggplot2::scale_y_continuous(breaks=c(1:max(df$X1)), labels=c(0:(max(df$X1)-1))) +
        ggplot2::ylab("Number of deleted reads") +
        ggplot2::ggtitle("% of reads after delete")
  #p5


  ## Variance analysis
  v<- apply(gene, 1, var)
  #x11()
  #plot(v, xlab = "Deleted reads", ylab = "Variance (On % of reads)", main = "Change in variance among gene")
  #x11()
  #plot(diff(v), xlab = "Deleted reads", ylab = "Change in Var (On % of reads)", main = "Change in variance by deleting per reads" )

  ## Plot % of missing reads vs cell
  pp <- t(res[,c(3,5)])
  #x11()
  #par(mar = par()$mar+c(0,0,0,10), xpd = TRUE)(res[,c(3,5)])
  colnames(pp)<- as.character(res[,1])
  #ll<- barplot(pp, col = c("red","green"), beside = TRUE, main = "Cell vs Reads reduction", xpd = TRUE,
  #             xlab = "Deleted reads", ylab = "Remaining (%)" ) #, args.legend = list(x = "topright", legend = c("Cells","Reads"), bty = "n", horiz = FALSE,xjust= 1,yjust=1))
  #legend("bottom", legend = rownames(pp) , lty=c(1,1), col = c("red","green"),lwd = c(5,5),
  #       bty = "n", xpd = TRUE, horiz = TRUE)#  xjust = 1, yjust = 1,

  ## Cell vs reads reduction
  df <- reshape::melt(pp)
  p6 <- ggplot2::ggplot(data=df, ggplot2::aes_string(x= "X2", y= "value")) +
        ggplot2::geom_bar(ggplot2::aes_string(fill = "X1"), position = "dodge", stat="identity") +
        ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 0, hjust = 1, vjust= -1)) +
        ggplot2::scale_x_continuous(breaks=c(0:max(df$X2)), labels=c(0:(max(df$X2)))) +
        ggplot2::ggtitle("% of cells") +
        ggplot2::ylab("Remaining (%)") +
        ggplot2::xlab("Deleted reads") +
        ggplot2::guides(fill= ggplot2::guide_legend(title=""))
  #p6

  mymultiplot(p1, p2, p3, p4,p5,p6, layout = matrix(c(1,2,3,4,5,6), nrow=3, ncol = 2,  byrow=TRUE))

  #on.exit(par(op))
  #return(gene)
  #return(list(res1,res2,res4,res5,gene,pp))

  ############################# Genes dependency (%)
  mydata  <- main_data
  mydata1 <- mydata > 0

  gene_lst <- as.list(sort(colnames(mydata)))
  res      <- lapply(gene_lst, function(object)
  {
    data1 <- mydata1[mydata1[,object],, drop = F]
    data2 <- colSums(data1)
    data2 <- round((data2/max(data2))*100,10)
    data2[which(data2 == max(data2))] <- max(data2) - 0.0000001 ## decrease max values by on for cut function later to plot correcty
    data2
  })
  names(res) <- unlist(gene_lst)
  res <- res[sort(names(res))]
  res_1 <- do.call(rbind,res)
  res_1 <- ((res_1/100))*100
  res_1 <- res_1[colnames(res_1),]
  #rownames(res_1) <- NULL
  
  if(coexp)
  {
    mydata  <- main_data
    pmdata  <- lapply(gene_lst, function(object)
    {
      mdata <- mydata[,object]
      pp1   <- length(mdata[mdata > 0])
      pp2   <- length(mdata[mdata > 1])
      pp    <- c(pp1,pp2)
      pp    <- (pp[2]/pp[1])*100-0.0000001
    }
    )
    names(pmdata) <- unlist(gene_lst)
    diag_res <- unlist(pmdata)
 
    ### Replace diagonal value 
    diag(res_1)<- diag_res[colnames(res_1)]
  }
  

  #library(reshape)
  heat_data    <- data.frame(reshape::melt(res_1))
  heat_data$X1 <- as.vector(heat_data$X1)
  heat_data$X2 <- as.vector(heat_data$X2)
  heat_data$value <- as.numeric(heat_data$value)
  #heat_data    <- data.frame(heat_data)

  if(length(intervel.dep) != 0)
  {
    heat_data$value <- cut(heat_data$value,breaks = c(0,intervel.dep,100),right = FALSE)
    heat_data <- data.frame(heat_data)
    #library(ggplot2)
    p <- ggplot2::ggplot(data =heat_data, ggplot2::aes_string(x = "X2", y= "X1") ) +
      ggplot2::geom_tile(ggplot2::aes_string(fill = "value")) +
      ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 90, hjust = 1, vjust= 0.5,
                                                          size = text.size, face="bold"),
                     axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     #axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y  = ggplot2::element_text(size = text.size, face="bold")) +
      #ggplot2::scale_y_continuous(breaks=1:length(unique(heat_data$X2)), labels=unique(as.vector(heat_data$X2)),
      #                            sec.axis = ggplot2::sec_axis(~ . * 1, breaks = 1:length(unique(heat_data$X2)),
      #                                                         labels = unique(as.vector(heat_data$X2)))) +
      #ggplot2::scale_y_discrete(position = c( "right")) +
      ggplot2::ggtitle("Genes dependency (%)") +
      ggplot2::coord_cartesian(ylim = c(4, length(unique(heat_data$X1))-3)) +
      ggplot2::scale_fill_brewer(palette = "Spectral", direction = -1) +
      ggplot2::guides(fill=ggplot2::guide_legend(title="Genes %"))
    print(p)
  }
  if (length(intervel.dep) == 0)
  {
    #heat_data$value <- factor(floor(heat_data$value))
    #heat_data <- data.frame(heat_data)
    #library(ggplot2)
    p <- ggplot2::ggplot(data =heat_data, ggplot2::aes_string(x = "X2", y= "X1") ) +
      ggplot2::geom_tile(ggplot2::aes_string(fill = "value")) +
      ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 90, hjust = 1, vjust= 0.5,
                                                          size = text.size, face="bold"),
                     #axis.title.x = ggplot2::element_blank(),
                     #axis.title.y = ggplot2::element_blank(),
                     #axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y  = ggplot2::element_text(size = text.size, face="bold")) +
      #ggplot2::scale_y_continuous(breaks=1:length(unique(heat_data$X2)), labels=sort(unique(heat_data$X2)),
      #ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ .)) +
      #ggplot2::scale_y_discrete(position = c( "right")) +
      #ggplot2::scale_y_discrete(position = c( "left")) +
      #ggplot2::ggtitle("Genes dependency (%)") +
      ggplot2::labs (x = "Distribution (also has)", y= "Condition (Cell having)", title = "Genes dependency (%)") +
      #ggplot2::coord_cartesian(ylim = c(0, length(unique(heat_data$X1))+1), expand = FALSE) +
      #ggplot2::scale_fill_brewer(palette = "Spectral")
      ggplot2::scale_fill_gradientn(name = "Genes %", colours = hm.palette(100))
   #print(p)

  }

  #p<- cowplot::ggdraw(cowplot:::switch_axis_position( p, axis="y", keep="y"))
  print(p)

  #return(res)
  return(heat_data)
}