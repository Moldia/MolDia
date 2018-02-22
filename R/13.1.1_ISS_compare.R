"ISS_compare"
#' Compare ISS data
#' 
#' @description Compare multiple ISS data
#'
#' 
#' @examples 
#' hc_left  <- readRCA(file = system.file("extdata", "Hypocampus_left.csv", package="MolDia"),  cellid = "CellId",centX = "centroid_x", centY = "centroid_y")
#' hc_right <- readRCA(file = system.file("extdata", "Hypocampus_right.csv", package="MolDia"), cellid = "CellId",centX = "centroid_x", centY = "centroid_y")
#' 
#' kk<- ISS_compare(hc_left,hc_right,hc_left,hc_right, label = F, levelCI = 0.99, live = F, logdata = FALSE)
#' 
#' @export
ISS_compare <- function(..., logdata = FALSE, label = TRUE, levelCI = 0.95, live = TRUE)
{
  options(warn=2)
  ## Main data
  maindata <- list(...)
  ## Stopping criteria
  if(length(maindata) < 2) stop("Provide at least 2 dataset", call. = TRUE)
  
  ## Generate counts per probe per data
  maindata <- lapply(seq_along(maindata), function(i)
  {
    data_1  <- data.frame(colSums(maindata[[i]]@data))
    colnames(data_1) <- paste0("Data.",i)
    data_1$probe <- rownames(data_1)
    data_1
  }
  )
  names(maindata)<- paste(paste0("data",1:length(maindata)))
  
  ## Merge dataset by probe name
  merge.all <- function(x, y) {
    merge(x, y, all=TRUE, by="probe")
  }
  maindata <- Reduce(merge.all, maindata)
  probe <- maindata$probe
  maindata$probe <- NULL
  
  ## Take all possible combinition
  mypairs <- t(combn(length(maindata),2))
  mypairs <- lapply(apply(mypairs,1,list),unlist)
  mypairs_name <- lapply(mypairs, function(i)
  {
    paste(paste0("Data.",i), collapse = "_")
    
  })
  mypairs <- lapply(mypairs, function(i) 
    {
    pp<- maindata[i]
    pp$probe <- probe
    #colnames(pp) <- c("d1","d2","probe")
    pp
    })
  names(mypairs) <- mypairs_name
  
  ## Taking log on data 
  if(logdata){
  mypairs <- lapply(mypairs, function(i)
  {
    i[,c(1,2)] <- log2(i[,c(1,2)])
    i
  }
    )}
  
  ## Plotting Function and apply it 
  myggplot <- lapply(seq_along(mypairs), function(i)
  {
    ## Define variable
    d1 <- colnames(mypairs[[i]])[1]
    d2 <- colnames(mypairs[[i]])[2]
    probe <- colnames(mypairs[[i]])[3]
    ## Find CO-efecient 
    f <- paste(d1, " ~ ", d2)
    m <- suppressMessages(lm(f, data= mypairs[[i]]))

        ## Plot
    pp<- ggplot(data = mypairs[[i]], aes_string(x = d1, y = d2,label = probe)) +
      {if(label) geom_text()}  +
      geom_point(color='blue') +
      suppressMessages(geom_smooth(method = "lm", se = T, show.legend = TRUE, level = levelCI)) + 
      suppressWarnings(labs(title = names(mypairs[i]), subtitle = substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2, 
                                                            list(a = format(coef(m)[1], digits = 2), 
                                                                 b = format(coef(m)[2], digits = 2), 
                                                                 r2 = format(summary(m)$r.squared, digits = 3))))) +
      scale_y_log10() +
      scale_x_log10()
    pp
    
  })
  
  ## Return multiple plot or single plot 
  if(live & length(mypairs)==1) print(ggplotly(myggplot[[1]]))
  else grid.arrange(grobs = myggplot)
  return(mypairs)
}
