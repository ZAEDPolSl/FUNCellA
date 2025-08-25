#' Thresholding pathway activity scores using k-means clustering
#'
#' This function applies k-means clustering to threshold pathway activity scores
#' for each pathway (row) individually. The optimal number of clusters is
#' estimated using the silhouette method. Samples belonging to the cluster
#' with the highest mean activity are marked as active (1); all others are marked inactive (0).
#'
#' @param df_path A numeric matrix or data.frame of pathway activity scores,
#'   where rows correspond to pathways and columns to samples.
#' @param K Maximum number of clusters to consider when estimating the optimal
#'   number of clusters using the silhouette method (default: 10).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{binary} – matrix of 0/1 assignments (cluster with relatvie pathway activation = 1).
#'   \item \code{clusters} – matrix of full cluster assignments (1…k).
#' }
#'
#'
#' @import factoextra
#' @import cli
#' @importFrom stats kmeans
#'
#' @export
thr_KM<-function(df_path,K = 10){

  idprog <- cli_progress_bar("Calculating thresholds", total=nrow(df_path))
  ret <- matrix(0,nrow=nrow(df_path),ncol=ncol(df_path))
  binary_mat  <- matrix(0, nrow = nrow(df_path), ncol = ncol(df_path))
  cluster_mat <- matrix(NA_integer_, nrow = nrow(df_path), ncol = ncol(df_path))

  for (i in seq_len(nrow(df_path))) {
    cli_progress_update(id = idprog)

    sample_values <- as.numeric(df_path[i, , drop = TRUE])
    sample_matrix <- matrix(sample_values, ncol = 1)

    opt_k <- fviz_nbclust(sample_matrix, FUNcluster = kmeans, method = "silhouette",
                          k.max = K, iter.max = 20)

    max_cluster <- as.numeric(opt_k$data$clusters[which.max(opt_k$data$y)])
    clut <- kmeans(sample_matrix, max_cluster)

    label <- which.max(clut$centers)
    binary_mat[i, clut$cluster == label] <- 1
    cluster_mat[i, ] <- clut$cluster


  }

  cli_progress_done(idprog)
  cli_alert_success('Thresholds calculated')

  rownames(binary_mat)<-rownames(cluster_mat)<-rownames(df_path)
  colnames(binary_mat)<-colnames(cluster_mat)<-colnames(df_path)
  return(list(
    binary   = binary_mat,
    clusters = cluster_mat
  ))
}
