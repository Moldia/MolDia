######################################################################
##              RCA data Filter: Probability based                  ##
######################################################################
"ISS_filter"
#' Filter ISS data based on poisson distibution
#'
#' @param data Input data in class MolDiaISS. Output of \link[MolDia]{readISS}
#' @param data_mean Expected mean of number of reads press cell. Default is NULL. See details.
#'
#' @description This function estimate the probabity of number of reads per cell in a specific range
#'              with desired mean (Rate or number or reads per cell) by Poisson distibution. By default
#'              current fuction calculate the mean from data.
#' @return Number of reads with peobability
#' @examples
#' data_1 <- readISS(file = system.file("extdata", "CellBlobs_QT_0.35.csv", package="MolDia"),
#'                   cellid = "CellID")
#' res    <- ISS_filter(data = data_1, data_mean = 10)
#'
#' @importFrom stats dpois
#'
#' @export
ISS_filter <- function(data, data_mean =NULL )
  
{
  ## Reading main data
  main_data <- data
  data      <- main_data@data
  ## Define frequency of reads
  freq       <- table(rowSums(data))
  reads_freq <- as.numeric(names(freq))
  
  ## Define mean of reads / reads per cell
  if(length(data_mean)==0) data_mean  <- sum(data)/nrow(data)
  
  ## Define true probability
  reads_prob_true <- dpois(x = reads_freq, lambda = sum(data)/nrow(data))
  reads_prob_true_1 <- freq
  reads_prob_true_1[1:length(reads_prob_true_1)] <- reads_prob_true
  reads_prob_true <- reads_prob_true_1
  
  ## Define probability of each reads frequency
  reads_prob <- dpois(x = reads_freq, lambda = data_mean)
  reads_prob_1 <- freq
  reads_prob_1[1:length(reads_prob_1)] <- reads_prob
  reads_prob <- reads_prob_1
  
  ## Round pribability to desired decimal points
  reads_prob_1 <- round(reads_prob,2)
  
  ## Keep probability greater than 0
  reads_prob_1 <- reads_prob_1[reads_prob_1!=0]
  
  ## Get the reads those has probability greater than 0
  reads_freq_1 <- as.numeric(names(reads_prob_1))
  reads_range  <- range(reads_freq_1)
  reads_prob_2 <- reads_prob [names(reads_prob_1)]
  
  reads_range  <- range(as.numeric(names(reads_prob_2)))
  cum_prob     <- sum(reads_prob_2)
  
  res <- list(cum_prob,reads_range,freq,reads_prob)
  names(res)<- c("prob","reads","freq","point.prob")
  
  ## Plot distribution
  data_1 <- res[[4]]
  data_2 <- res[[4]][which(as.numeric(names(res[[4]]))%in%res[[2]][1]:res[[2]][2])]
  x_axis <- as.numeric(names(data_2))
  x_axis <- c(x_axis[1],x_axis,x_axis[length(x_axis)])
  y_axis <- c(0,as.vector(data_2),0)
  
  r1 <- res[[2]]
  r2 <- round(res[[1]],7)
  
  ## Data output
  data_reads <- main_data@data[rowSums(main_data@data) %in% reads_range[1]:reads_range[2],]
  res <- methods::new("MolDiaISS",
                      data     = data_reads,
                      location = main_data@location[rownames(data_reads),],
                      gene = colnames(data_reads))
  ## Plot
  p<- ggplot2::ggplot(data = data.frame(as.vector(data_1))) +
    ggplot2::geom_line(ggplot2::aes(x=as.numeric(names(reads_prob_true)),as.vector(reads_prob_true))) +
    ggplot2::geom_line(ggplot2::aes(x=as.numeric(names(data_1)),as.vector(data_1))) +
    ggplot2::geom_polygon(data = data.frame(cbind(x_axis,y_axis)), ggplot2::aes_string(x = "x_axis",y = "y_axis"), colour = "red", alpha=I(0.1)) +
    ggplot2::geom_text(ggplot2::aes(x=Inf,y=Inf,hjust=1,vjust=1,label= paste("x = ",round(data_mean,2),"(Reads per cell)\n",
                                                                             "p(",r1[1],"< x <",r1[2],") = ", r2,"\n",
                                                                             "Number of cells: ",nrow(data_reads) ))) +
    ggplot2::labs(x = "Number of Reads/cell", y = "Probability", title="Probability distribution of number of reads/cell")
  
  ## Frequency distibution plot before
  mm1<- data.frame(table(rowSums(main_data@data)))
  p1 <- ggplot2::ggplot(data=mm1, ggplot2::aes_string(x= "Var1", y= "Freq")) +
    ggplot2::geom_bar(stat="identity", width=0.5) +
    ggplot2::labs(x = "Number of reads", y = "Numbe of cells", subtitle = "Before exclude outlier")+
    ggplot2::geom_text(ggplot2::aes(x=Inf,y=Inf,hjust=1,vjust=1,label= paste("x = ",round(sum(main_data@data)/nrow(main_data@data),2),"(Reads/Cell)\n",
                                                                             "Cells number: ",nrow(main_data@data))))
  
  mm2<- data.frame(table(rowSums(data_reads)))
  p2 <- ggplot2::ggplot(data=mm2, ggplot2::aes_string(x= "Var1", y= "Freq")) +
    ggplot2::geom_bar(stat="identity", width=0.5) +
    ggplot2::labs(x = "Number of reads", y = "Numbe of cells", subtitle = " After exclude outlier") +
    ggplot2::geom_text(ggplot2::aes(x=Inf,y=Inf,hjust=1,vjust=1,label= paste("x = ",round(sum(data_reads)/nrow(data_reads),2),"(Reads/Cell)\n",
                                                                             "Cells number: ",nrow(data_reads))))
  
  mymultiplot(p, p1, p2, layout = matrix(c(1,1,2,1,1,3), nrow=2, ncol = 3,  byrow=TRUE))
  
  #print(p1)
  return(res)
}

