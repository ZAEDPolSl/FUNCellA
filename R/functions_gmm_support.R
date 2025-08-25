#' Find the Last Run of Consecutive TRUE Values
#'
#' This function identifies the indices of the last contiguous (consecutive) sequence
#' of \code{TRUE} values within a logical vector.
#'
#' @param row A logical vector (e.g., output of a comparison or logical condition).
#'
#' @return An integer vector of indices corresponding to the last consecutive \code{TRUE} segment.
#' If no \code{TRUE} values are found, returns \code{NA}.
#'
#' @examples
#' \dontrun{
#' x <- c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)
#' last_consecutive_true(x)
#' # Returns: 5 6 7
#'
#' y <- c(FALSE, FALSE, FALSE)
#' last_consecutive_true(y)
#' # Returns: NA
#'}
#' @importFrom utils tail
#' @keywords internal
last_consecutive_true <- function(row) {
  true_indices <- which(row)
  if (length(true_indices) == 0) {
    return(NA)
  }
  runs <- split(true_indices, cumsum(c(1, diff(true_indices) != 1)))
  last_run <- tail(runs, 1)[[1]]
  return(last_run)
}

#' Determine GMM Threshold Using K-Means Clustering
#'
#' This function selects an optimal threshold from a Gaussian Mixture Model (GMM) decomposition,
#' using k-means clustering on the GMM component parameters (mean, standard deviation, and weight).
#' It uses silhouette analysis to estimate the optimal number of clusters and then identifies
#' the most relevant cluster for defining the threshold.
#'
#' If the number of components is one or two, the function simply returns the maximum threshold.
#' Otherwise, it scales the parameters, performs k-means clustering, and selects the minimum
#' threshold corresponding to the cluster with the highest mean.
#'
#' @param params A numeric matrix of GMM decomposition parameters. Each row corresponds to one component
#'   (cluster), with columns in the order: mean, standard deviation, and alpha (mixing proportion).
#' @param thrs A numeric vector of thresholds produced during GMM decomposition (typically one per component).
#'
#' @return A single numeric value representing the selected GMM threshold.
#'
#' @import factoextra
#' @importFrom stats kmeans
#' @keywords internal
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
