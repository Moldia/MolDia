"ISS_ClustCompare"
#' Compare among ISS clusters based on Random Forest algorithm
#'
#' @description Compare two cluster and find their similarity (Supervised).
#' @param preditorData Data to be used to compare. Input data in class MolDiaISS. Output of \link[MolDia]{readISS}
#'        and must be clustered.
#' @param prediction Data that need to compare. Input data in class MolDiaISS. Output of \link[MolDia]{readISS}
#'        and must be clustered.
#' @param ntree Number of trees to grow. \link[randomForest]{randomForest}. See This should not be set to
#'        too small a number, to ensure that every input row gets predicted at least a few times.
#' @param method Methods to return comarison results. Possible values are "vote" and "prob". Defaut is "prob".
#'
#' @details This function will be apply to compare Two cluster by Random forest algorithm.
#' 
#' @note This function need optimization to run properly. Work fine with small data but not with big data .
#'
#' @importFrom foreach %dopar%
#'
#' @examples
#' ## Read data: Left and right HC
#' left_hypo  <- readISS(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),
#'                   cellid = "CellId",centX = "centroid_x", centY = "centroid_y")
#' right_hypo <- readISS(file = system.file("extdata", "Hypocampus_right.csv", package="MolDia"),
#'                   cellid = "CellId",centX = "centroid_x", centY = "centroid_y")
#'
#' ## Data preprocessing
#' left_hypo  <- ISS_preprocess(data = left_hypo, normalization.method = "LogNormalize",
#'                                do.scale = TRUE, do.center = TRUE)
#' right_hypo <- ISS_preprocess(data = right_hypo, normalization.method = "LogNormalize",
#'                                do.scale = TRUE, do.center = TRUE)
#'
#' ## Cluster data based on SEURAT pipeline
#' left_hypo   <- ISS_cluster(data =  left_hypo, pc = 0.9, resolution = 0.05)
#' right_hypo  <- ISS_cluster(data = right_hypo, pc = 0.9, resolution = 0.05)
#'
#' ## Cluster compare
#' comapre_left_right <- ISS_ClustCompare(preditorData = left_hypo, prediction = right_hypo, ntree = 10)
#' comapre_right_left <- ISS_ClustCompare(preditorData = right_hypo, prediction = left_hypo, ntree = 10)
#'
#'
#' @export
ISS_ClustCompare <- function(preditorData, prediction, ntree = 100, method = "prob")
{
  #### Library
  #library("foreach")
  #library("doSNOW")
  #library("parallel")
  #library("randomForest")

  #### Create predictor data
  predictors          <- as.data.frame((as.matrix(preditorData@data)))
  predictors$subtype  <- as.vector(preditorData@cluster)
  predictors$subtype  <- as.factor(predictors$subtype)

  ## Creat prediction data
  predictionData <- as.data.frame((as.matrix(prediction@data)))
  predictionData$subtype <- as.vector(prediction@cluster)
  predictionData         <- split(x = predictionData, f=predictionData[,"subtype"])
  predictionData         <- lapply(predictionData, function(i)
  {i[,"subtype"]<- NULL
  res <- apply(i,2,as.numeric)
  rownames(res)<- rownames(i)
  res})

  #### Random forest settings
  process_core <-  parallel::detectCores()
  ntree        <- floor(ntree/process_core)

  #### Parallal random forest model
  cl    <- parallel::makeCluster(process_core , type="SOCK")
  doSNOW::registerDoSNOW(cl)
  rf_parallal <- foreach::foreach(ntree = ntree, .combine = randomForest::combine,.multicombine=TRUE, .packages = "randomForest") %dopar%  {
    set.seed(10000001)
    randomForest::randomForest(subtype ~ ., data=predictors, importance=TRUE, proximity=TRUE, ntree = ntree)
  }
  parallel::stopCluster(cl)

  ############################################### Prediction
  #### Prediction of one cluster with another cluster. Compare (Supervised) between two two cluster.

  if(method == "prob"){
    pred       <- lapply(1:length(predictionData), function(i)
    {
      res <- stats::predict(rf_parallal, predictionData[[i]], type="prob")
      clust_mean <- data.frame(as.matrix(apply(res,2,mean),nrow=1))
      colnames(clust_mean)<- names(predictionData[i])
      clust_mean
    })
    clust_mean <- do.call(cbind, pred)}

  if(method == "vote"){
    pred       <- lapply(1:length(predictionData), function(i)
    {
      res <- stats::predict(rf_parallal, predictionData[[i]], type="prob")
      #clust_mean <- data.frame(as.matrix(apply(res,2,mean),nrow=1))
      #colnames(clust_mean)<- names(predictionData[i])

      kk1 <- apply (res,1,function(j)
      {
        pp<- which(j == max(j))
        pp<- names(pp)
        pp
      }
      )
      kk1 <- floor((table(unlist(kk1))/ sum(table(unlist(kk1))))*100)
      #kk1 <- table(unlist(kk1))

      kk<- rep(0,length(colnames(res)))
      names(kk)<- colnames(res)
      kk[names(kk1)] <- kk1
      nam <- names(predictionData[i])
      kk <- data.frame(kk)
      colnames(kk) <- nam
      kk
    })

    names(pred) <- names(predictionData)
    clust_mean <- do.call(cbind,pred)
  }


  #### Assign new cluster based on predictor. (Consider highest probability)
  pred1       <- lapply(1:length(predictionData), function(i)
  {
    res <- stats::predict(rf_parallal, predictionData[[i]], type="prob")
    clust <- apply(res,1,function(j)
    {
      pp<- which(j == max(j))
      pp<- colnames(res)[pp][1]
      pp
    })
    clust <- data.frame(clust)
  })
  names(pred1)<- names(predictionData)

  #### Probability of each cell belong to predictor cluster.
  pred2       <- lapply(1:length(predictionData), function(i) { res <- stats::predict(rf_parallal, predictionData[[i]], type="prob")})
  names(pred2)<- names(predictionData)

  #### Get new cluster in Seurat formate
  pp<- pred1
  pp<- lapply(pp,function(i)
  {i$cell<- rownames(i)
  i})
  pp<- do.call(rbind,pp)
  rownames(pp)<- pp$cell
  pp$cell <- NULL
  cell_name <- rownames(pp)
  pp<- as.vector((pp[,1]))
  names(pp) <- cell_name
  pp<- as.factor(pp)

  prediction@cluster <- pp

  ## ggplot object
  clust_mean <- as.matrix(clust_mean)
  my_data <- reshape::melt(clust_mean)
  hm.palette  <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')), space='Lab')
  p<- ggplot2::ggplot(data = my_data, ggplot2::aes_string(x = "X1", y = "X2" )) +
      ggplot2::geom_tile(ggplot2::aes_string(fill = "value")) +
      ggplot2::scale_fill_gradientn(name = "P(match)",colours = hm.palette(100)) +
      ggplot2::scale_y_continuous(breaks=0:nrow(clust_mean), labels=0: nrow(clust_mean)) +
      ggplot2::scale_x_continuous(breaks=0:ncol(clust_mean), labels=0: ncol(clust_mean)) +
      ggplot2::labs(x = "Predictor", y = "Prediction")

  print(p)

  ## Return result
  res <- list(cluster_sim = clust_mean, cluster_new = pred1, clust_prob = pred2, seurat_new = prediction)
  #return(res)
  return(clust_mean)
}

