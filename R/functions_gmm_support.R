#' Function to find last consecutive true
#'
#' The function...
#'
#' @param row vector of data
#'
#' @return A vector of ranked genes.
#'
last_consecutive_true <- function(row) {
  true_indices <- which(row)
  if (length(true_indices) == 0) {
    return(NA)
  }
  runs <- split(true_indices, cumsum(c(1, diff(true_indices) != 1)))
  last_run <- tail(runs, 1)[[1]]
  return(last_run)
}

#' Function to find threshold based on GMM parameters
#'
#' The function...
#'
#' @param params matrix of GMM decomposition parameters. Rows represent clusters while columns represent GMM parameters with following order: mean, sd, alpha.
#' @param thrs vector of thresholds of GMM decomposition.
#'
#' @return a numeric value with threshold.
#'
#' @import factoextra
#'
km_search<-function(params,thrs){
  for(c in 1:ncol(params)){if(length(unique(params[,c])) > 1){params[,c] <- scale(params[,c])}}
  opt_k <- fviz_nbclust(params, FUNcluster = kmeans, method = "silhouette", k.max = nrow(params)-1, iter.max = 20)
  max_cluster <- as.numeric(opt_k$data$clusters[which.max(opt_k$data$y)])
  clut <- kmeans(params, max_cluster)


  if(length(unique(clut$cluster)) == nrow(params)){
    thr<-max(thrs,na.rm=T)
  } else {
    ord<-(clut$cluster==rownames(clut$centers)[which.max(clut$centers[,1])])
    if (all(diff(ord) >= 0)){
      take<-which(ord)
      thr<-min(thrs[take-1],na.rm=T)

    }else{
      if (ord[length(ord)]==F){
        thr=max(thrs,na.rm = T)
      } else{
        take<-last_consecutive_true(ord)
        thr<-min(thrs[take-1],na.rm=T)
      }
    }

  }
  return(thr)
}
