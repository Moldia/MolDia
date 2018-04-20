######################################################################
##              ISS : Correlation and ration of genes groups        ##
######################################################################
"ISS_ratiocor"
#' Calculate and plot correlation and ratio of total reads between genes
#' 
#' @description Calculate and plot correlation and ratio of total reads between genes
#' 
#' @param data Data in list of groups with group name.   Individual data is in class RCA_class. Output of \link[MolDia]{readRCA} 
#' @param gene Object in vector or list formate. In list formated input every list element is a group of
#'              interested genes.
#' @param select_gene Select gene of interest, Default is NULL i.e. all genes in gele list.
#' @param errorbar Show error bar or not. Default is TRUE
#' @param plty Show plot type. Available balue is "corr", "ratio" and "both". Default is "both".
#' @param logratio Taking log2 of the ratio of gene's total reads count.
#' @param sig.level confidence level for the returned confidence interval. Currently only used for 
#'             the Pearson product moment correlation coefficient if there are at least 4 complete pairs of observations.
#' @param main Title of the plot 
#' @param stat Mode of operation. Possible value is "sum", "gene" and "present". Default is "sum". 
#' 
#'@export

ISS_ratiocor <- function(data, gene = marker_gene, select_gene = NULL, errorbar= TRUE, plty = "both", 
                         logratio = TRUE, sig.level = 0.05, main = "", stat = "sum")
{
  ## Combine data and calculate correlation and ratio
  mydata <- lapply(seq_along (data), function(i)
  {
    kk1<- ISS_ratiocor_2(data[[i]], group = names(data[i]), gene = gene, stat = stat, logratio = logratio, sig.level = sig.level, main = main)
  })
  names(mydata) <- names(data)
  
  ## Extract name form input data 
  grpname <- mydata
  grpname <- unlist(lapply(grpname, function(i) levels(as.factor(i$group))))
  
  ## Combine all data
  mydata <- do.call(rbind,mydata)
  mydata$group <- factor(mydata$group, levels = grpname)
  
  mydata$X1 <- factor(mydata$X1)
  mydata$X2 <- factor(mydata$X2)
  mydata <- mydata[-which(mydata$X1 == mydata$X2),]
  mydata <- split(mydata, mydata$X1)
  
  ## Select group
  if(length(select_gene) > 0) mydata<- mydata[select_gene]
  
  ## Plot
  p  <- list()
  kk <- for(i in 1 : length(mydata))
  {
    kk<- mydata[i]
    kk1 <- reshape2::melt(kk[1], id.vars = c("X1","X2","group"),measure.vars = c("correlation_mean", "ratio_mean"), variable.name = "mean", value.name = "mean_value")
    kk1$L1 <- NULL
    kk2 <- reshape2::melt(kk[1], id.vars = c("X1","X2","group"),measure.vars = c("correlation_sd", "ratio_sd"), variable.name = "sd", value.name = "sd_value")
    kk2$L1 <- NULL
    kk1 <- list(kk1,kk2)
    pp<- do.call(cbind,kk1)
    kk1 <- pp[,c(1:5,9:10)]
    
    if(plty == "both")
    {
      p[[i]] <- ggplot2::ggplot(kk1 , ggplot2::aes(x = X2, y = mean_value, fill = group)) +
                ggplot2::geom_bar(position = ggplot2::position_dodge(), stat="identity") +
                {if(errorbar) ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), width=.2, position = ggplot2::position_dodge(.9))} +
                ggplot2::labs(x=names(kk), y = "") + 
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)) +
                ggplot2::facet_grid(mean ~ ., scales = "free_y") +
                ggplot2::labs(title = main)
    }
    if(plty == "corr")
    {
      kk1 <- subset(x= kk1, mean == "correlation_mean")
      p[[i]]  <- ggplot2::ggplot(kk1 , ggplot2::aes(x = X2, y = mean_value, fill = group)) +
                 ggplot2::geom_bar(position = ggplot2::position_dodge(), stat="identity") +
                 {if(errorbar) ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), width=.2, position = ggplot2::position_dodge(.9))} +
                 ggplot2::labs(x=names(kk), y = "") + 
                 ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))+
                 ggplot2::labs(title = main)
    }
    if(plty == "ratio")
    {
      kk1 <- subset(x= kk1, mean == "ratio_mean")
      p[[i]]  <- ggplot2::ggplot(kk1 , ggplot2::aes(x = X2, y = mean_value, fill = group)) +
                 ggplot2::geom_bar(position = ggplot2::position_dodge(), stat="identity") +
                 {if(errorbar) ggplot2::geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), width=.2, position = ggplot2::position_dodge(.9))} +
                 ggplot2::labs(x=names(kk), y = "") + 
                 ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))+
                 ggplot2::labs(title = main)
    }
  }
  
  ## Plot
  gridExtra::grid.arrange(grobs = p)
  
  ## Return Data
  return(mydata)
}


ISS_ratiocor_2<- function(data = list(...), group = "Sample_group", gene = gene, stat = stat, logratio = logratio, sig.level = sig.level, main = main)
  {

  ll <- data
  sample_number <- length(ll)
  if(length(ll) == 1 ) ll <- list(ll[[1]],ll[[1]])
  ll <- lapply(ll, ISS_ratiocor_1,gene, stat, logratio, sig.level, main)
  
  
  kk<- list()
  for(i in 1:3) ## How many feature we want
  {
    ## Data
    ll <- ll
    
    ## names of function
    nam <- unique(unlist(lapply(ll,names)))[i]

    ## Assign values
    val      <-lapply(ll,"[[",i)
    val_mean <- apply(simplify2array(val), 1:2, mean)
    assign(paste(nam,"mean",sep="_"),val_mean)
    val_sd   <- apply(simplify2array(val), 1:2, sd)
    assign(paste(nam,"sd",sep="_"),val_sd)
    
    ## Get values
    val_mean   <- paste(nam,"mean",sep="_")
    val_mean_1 <- get (val_mean)
    val_mean_1 <- reshape::melt(val_mean_1)
    
    val_sd   <-paste(nam,"sd",sep="_")
    val_sd_1 <- get (val_sd)
    val_sd_1 <- reshape::melt(val_sd_1)
    
    ## Merge values
    mm <- merge(val_mean_1, val_sd_1, by = c('X1', 'X2'), all = T)
    colnames(mm) <- c("X1","X2",val_mean,val_sd)
    kk[[i]]<- mm
  }
  
  kk1 <- Reduce(function(...) merge(..., by=c("X1","X2"), all=TRUE), kk)
  kk1$group <- group
  kk1$n <- sample_number
  return(kk1)
}

ISS_ratiocor_1 <- function(data, gene = marker_gene, stat = stat, logratio = logratio, sig.level = sig.level, main = main )
{
  ## Loading data 
  my_data <- data@data
  
  ## Special data matrix
  res <- ISS_sumstat(data = my_data, gene = gene, stat = stat)
  
  ## Compute ratio
  mm <- colSums(res)
  if (logratio)  {data_ratio <-log2(t(outer(mm,mm,"/")))
  } else { 
    data_ratio <- t(outer(mm,mm,"/"))}
  
  ## Compute p value of correlation
  pv  <- ggcorrplot::cor_pmat(res, conf.level = (1-sig.level))
  data_cor  <- cor(res)
  
  ## Correlation plot
  p<- ggcorrplot::ggcorrplot(corr = data_cor, hc.order = FALSE, type = "lower", lab = TRUE, tl.cex = 11,tl.srt= 90, 
                 p.mat = pv, sig.level = sig.level, colors = c("red","white","green"), 
                 outline.col = "white", title = paste("Correlation: ", main))
  #print(p)
  ## Return result
  res1 <- list(data_cor, pv, data_ratio)
  names(res1) <- c("correlation", "p_correlation", "ratio")
  return(res1)
}



#########################################################################
##              Calculate summary statstics of data                    ##
## Total reads, Total number of gene present and Gene present in group ##
#########################################################################
"ISS_sumstat"
#' Calculate summary statstics of data
#' 
#' @description Calculate summary statstics (Total reads, Total number of gene present and Gene present in group) of data
#'  
#' @param data Data in dataframe class
#' @param gene Object in vector or list formate. In list formated input every list element is a group of
#'             interested genes
#' @param stat Mode of operation. Possible value is "sum", "gene" and "present". Default is "sum"

ISS_sumstat <- function(data, gene, stat = "sum")
  {
  ##Create empty result data frame
  res <- matrix(NA,nrow = nrow(data), ncol = length(gene))
  colnames(res) <- names(gene)
  rownames(res) <- rownames(data)
  
  #### Main function to calculate summary statisitcs (Sum, Total gene present, Gene group present)
  ## Helper function
  myfun <- function(v,stat = stat)
  {
    if(stat == "sum")     res <- sum(v)
    if(stat == "gene")    res <- sum(v>0)
    if(stat == "present") res <- as.numeric(any(v>0))
    return(res)
  }
  
  ## Main function
  kk<- for(i in 1 : length(gene))
  {
    ## Select gene which is in data 
    gene[[i]] <- gene[[i]][gene[[i]]%in%colnames(data)]
    
    data1 <- data[,gene[[i]], drop = FALSE]
    res1  <- apply(data1,1,function(j)
    {
      kk<- myfun(j,stat)
      return(kk)
    })
    res[,i] <- res1
  }
  
  ## Return result
  res <- data.frame(res)
  res$total_reads <- as.numeric(rowSums(data))
  
  ## Delete column (gene) with no reads
  res <- res[,colSums(res) > 0]
  return(res)
}