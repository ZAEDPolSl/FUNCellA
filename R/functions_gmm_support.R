#' Function to find last consecutive true
#'
#' The function...
#'
#' @param row vector of data
#'
#' @return A vector of ranked genes.
#'
#' @importFrom utils tail
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
#' @importFrom stats kmeans
#'
km_search<-function(params,thrs){

  # One or two distribution
  if (nrow(params) <= 2) {
    return(max(thrs, na.rm = TRUE))
  }

  # ---- K-means part ----
  # scaling parameters
  params <- apply(params, 2, function(col) if (length(unique(col)) > 1) scale(col) else col)

  # k-means run
  opt_k <- fviz_nbclust(params, FUNcluster = kmeans, method = "silhouette", k.max = nrow(params)-1, iter.max = 20)
  max_cluster <- as.numeric(opt_k$data$clusters[which.max(opt_k$data$y)])
  clusters <- kmeans(params, max_cluster)

  # ---- threshold with k-means ----
  # Edge case: no of cluster equals no of parameters
  if (length(unique(clusters$cluster)) == nrow(params)) {
    return(max(thrs, na.rm = TRUE))
  }

  ord<-(clusters$cluster==rownames(clusters$centers)[which.max(clusters$centers[,1])])
  if (all(diff(ord) >= 0)){
      take<-which(ord)
      threshold<-min(thrs[take-1],na.rm=T)

  }else{
      if (ord[length(ord)]==F){
        threshold=max(thrs,na.rm = T)
      } else{
        take<-last_consecutive_true(ord)
        threshold<-min(thrs[take-1],na.rm=T)
      }
  }


  return(threshold)
}
