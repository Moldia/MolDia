######################################################################
##              RCA data Filter: Probability based                  ##
######################################################################
"RCA_filter"
#' Filter RCA data based on poisson distibution
#'
#' @param data Input data in class RCA_class. Output of \link[MolDia]{readRCA}
#' @param data_mean Expected mean of number of reads press cell. Default is NULL. See details.
#'
#' @description This function estimate the probabity of number of reads per cell in a specific range
#'              with desired mean (Rate or number or reads per cell) by Poisson distibution. By default
#'              current fuction calculate the mean from data.
#' @return Number of reads with peobability
#' @examples
#' data_1 <- readRCA(file = system.file("extdata", "CellBlobs_QT_0.35.csv", package="MolDia"),
#'                   cellid = "CellID")
#' res    <- RCA_filter(data = data_1, data_mean = 10)
#'
#' @importFrom stats dpois
#'
#' @export
RCA_filter <- function(data, data_mean =NULL )

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
  res <- methods::new("RCA_class",
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



######################################################################
##                   RCA data Filter : Barplot                      ##
######################################################################
"RCA_barplot"
#' Plot barplot on RCA data bsed on different condition
#' @description Plot bar plot on RCA data bsed on different condition
#'
#' @param data Input data in class RCA_class. Output of \link[MolDia]{readRCA}.
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
#' @return Data in class \link[MolDia]{readRCA}.
#'
#'         The barplot will explain the total number of reads per gene in 'total.expr' (Default 1e4) cells. The yellow line
#'         inside the bar plot indicate the number of expressed cells and percentage count on each bar is the percent of cells
#'         express that gene.
#' @importFrom graphics barplot lines text legend
#' @author Mohammad Tanvir Ahamed
#' @examples
#' ########## Reading data
#' data_1      <- readRCA(file = system.file("extdata", "CellBlobs_QT_0.35.csv", package="MolDia"),
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
#' all_data     <- RCA_barplot(data = data_1, gene = marker_gene,gene.target = NULL)
#' neuron_group <- RCA_barplot(data = data_1, gene = marker_gene1, gene.target = 1,
#'                             target.min.count.cell = 2, rest.min.count.cell = 2,
#'                             at.least.gene = 1, gene.show = NULL, main = " Neuron group (>=1 genes)")
#'
#'
#'
#' @export
RCA_barplot <- function(data, gene, total.expr = 1e4, gene.target = NULL, gene.show = NULL, target.min.count.cell = 1, rest.min.count.cell = 1,
                        min.count = 1, at.least.gene = 0, at.most.gene = ncol(data), show.same.gene = FALSE,
                        str.same.gene = c(".C",".L"), main = "")
{
  ## Extract data from RCA_class object
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
  #res <- methods::new("RCA_class",
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


######################################################################
##                        Select Grid of interest                   ##
######################################################################
"RCA_GridSelect"
#' Select grid of interest from a tissue.
#' @description Select grid of interest from a tissue
#' @param data Input data in class RCA_class. Output of \link[MolDia]{readRCA}.
#' @param gridtype type of grid to plot. Default is "rect". See details.
#' @param grid_id Grid to select. Default id NULL.
#' @param nx,ny Numbers of rectangular quadrats in the x and y directions
#' @param roifile Name of the file that contain ROI. csv formate
#' @param roi.id Column name in roifile file that contain ROI id
#' @param roi.x,roi.y X and Y axis name in roifile file.
#' 
#' @details gridtype parameter can have 4 values :  "hexa", "rect", "roifile" and "roi"
#'          "hexa" will create hexagonal grid and number of grid is based on nx parameter and grid_id will select grid of interest.
#'          
#'          
#'          "rect" will create rectangular grid and number of grid is based on nx parameter and grid_id will select grid of interest.
#'          
#'          "roifile" will require a file input in csv formate in the parameter "roifile". The input file must have at least 3 column of
#'          ROI id, x-axis, y-axis. The parameter "roi.id", "roi.x" and "roi.y" will be the name input of corresponding column.
#'          
#'          "roi" will select ROI more interactively. One cal select ROI on the image by pointer. 
#' 
#' @examples
#' ex_data <- readRCA(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),cellid = "CellId", centX = "centroid_x", centY = "centroid_y", rpc = 3)
#' mygrid  <- RCA_GridSelect(data = ex_data, nx = 8,gridtype = "hexa")
#' mygrid  <- RCA_GridSelect(data = ex_data, nx = 8,gridtype = "rect", grid_id = c(6,16,20,8,17,21))
#' 
#' ## Selected ROI
#' ex_data <- readRCA(file = system.file("extdata", "CellBlobs_ROI.csv", package="MolDia"),
#'                    cellid = "CellID", centX = "centroidX", centY = "centroidY")
#' mygrid  <- RCA_GridSelect(data = ex_data, nx = 6, gridtype = "roifile",
#'                             roifile = system.file("extdata", "polygon_coordinates.csv", package="MolDia"),
#'                             roi.id = "Polygon_id", roi.x ="x_coordiates" , roi.y = "y_coordinates", grid_id = c(1,2,5,6))
#' mygrid  <- RCA_GridSelect(data = ex_data, gridtype = "roi")
#' 
#' @export
RCA_GridSelect <- function(data, gridtype = "rect", nx = 6, ny = nx, grid_id = NULL, 
                           roifile = NULL, roi.id = NULL, roi.x = NULL, roi.y = NULL)
{
  ## Save main data
  main_data <- data
  
  ## Main data
  mydata <- data@location
  
  ## Find the convex point
  hpts <- grDevices::chull(mydata)
  hpts <- mydata[hpts, ]
  
  ## Sort points by anti clockwise
  anti_hpts <- contoureR::orderPoints(x = hpts$centroid_x, y = hpts$centroid_y, clockwise = FALSE)
  hpts <- hpts[anti_hpts,]
  
  ## Spatial point area
  myspa <- spatstat::ppp(x = mydata$centroid_x, y = mydata$centroid_y, 
                         poly=list(x=hpts$centroid_x, y=hpts$centroid_y), check = FALSE)
  plot(myspa, main = "ROI on Tissue")
  
  ## Divides window into quadrats and counts the numbers
  if(gridtype == "rect")
  {
    myspa1 <- spatstat::quadrats(myspa, nx = nx, ny = ny)
    #names(myspa1$tiles) <- paste("Tile",c(1:length(myspa1$tiles)))
    names(myspa1$tiles) <- c(1:length(myspa1$tiles))
    plot(myspa1, add= TRUE,col= "red",do.labels=TRUE, labelargs = list(col = "red"))
  }
  
  if(gridtype == "hexa")
  {
    max_range <- max(apply(apply(mydata,2,range),2,diff))
    myspa1    <- spatstat::hextess(myspa, s = abs(max_range)/nx)
    names(myspa1$tiles) <- c(1:length(myspa1$tiles))
    plot(myspa1, add= TRUE,col= "red",do.labels=TRUE, labelargs = list(col = "red"))
  }
  
  if(gridtype == "roifile")
  {
    if(length(roifile)==0) stop("Please select location of ROI file in csv formate",call. = TRUE)
    ## Read ROI file in CSV formate
    roi <- read.csv(file = roifile)
    
    ## Split ROI file by region id 
    roi <- split(roi, roi[,roi.id])
    roi <- lapply(roi, function(i)
    { 
      hpts <- i #subset(i, select=-c(Polygon.id))
      hpts[,roi.id] <- NULL
      anti_hpts <- contoureR::orderPoints(x = hpts[,roi.x], y = hpts[,roi.y], clockwise = FALSE)
      hpts <- hpts[anti_hpts,]
      hpts <- apply(hpts,2,list)
      hpts <- lapply(hpts, unlist)
      hpts <- spatstat::owin(poly =  list(x = hpts[[roi.x]],y = hpts[[roi.y]] ))
      hpts
    }
    )
    
    ## Ploting selected region
    myspa1 <-lapply(roi, plot, add=T, border = "red", lwd = 2)
    ## label each polygon
    labs<-names(roi)
    for (i in 1: length(roi)) {xy<-spatstat::centroid.owin((poly = roi[[i]]));
    text(xy$x,xy$y, labels = labs[i], col = "red")}
    ## Convert window object to Tessellation
    myspa1 <-spatstat::as.tess(myspa1)
  }
  
  if(gridtype == "roi")
  {
    myspa1 <- spatstat::clickpoly(add = TRUE, col = 2, lwd = 2)
    myspa1 <- spatstat::as.tess(myspa1)
  }

  
  ## Select grid of interest
  if(length(grid_id) > 0 | gridtype == "roi"){ 
    
    if(gridtype != "roi")
      {
      pp <- myspa1
      tt <- pp$tiles
      tile_id <- paste0("Grid_id_", grid_id)
      for(i in 1:length(seq_along(grid_id)))
        { 
        nam <- tile_id[i]
        assign(nam, tt[[grid_id[i]]])
        }
      eval(parse(text = paste("kk1 <- spatstat::union.owin(",paste0(tile_id,collapse=", "),")")))
      }
    
    if(gridtype == "roi") kk1 <- myspa1$tiles[[1]]
    isin <- spatstat::inside.owin(x = mydata$centroid_x, y = mydata$centroid_y,w=kk1)
    
    ## Plot grid of interest
    point_in <- main_data@location[isin,]
    #point_in1 <- spatstat::ppp(x = point_in$centroid_x, y = point_in$centroid_y,window = kk1, check = FALSE)
    
    ## Return result
    final_data <-  main_data@data[rownames(point_in), , drop = FALSE]
    final_data <- final_data[,colSums(final_data)>0, drop = FALSE]
    main_data@data <- final_data
    if(length(main_data@norm.data)  > 0 ) main_data@norm.data  <- main_data@norm.data[rownames(point_in), , drop = FALSE]
    if(length(main_data@scale.data) > 0 ) main_data@scale.data <- main_data@scale.data[rownames(point_in), drop = FALSE] 
    main_data@location <- main_data@location[rownames(point_in),]
    main_data@gene <- colnames(final_data)
    
    RCA_map(main_data, main = "Selected ROI from Tissue")
  }
  
  return(main_data)
}

